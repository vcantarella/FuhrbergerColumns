using CairoMakie
using DataFrames, XLSX, Statistics
using Dates
using DataInterpolations
using QuadGK
include("model_data_structures.jl")
include("prepare_lab_data_model_m2.jl")

# Load the prepared data
function get_data()
    return v_interp, c_ins, all_ds, q_disch, v_da
end
v_interp, c_ins, all_ds, q_disch, v_da = get_data()

# --- Helper Functions (Lagrangian) ---

function make_v_func(v_da::VDataA)
    function v_inst(t)
        for i in eachindex(v_da.end_times)
            if t ≤ v_da.end_times[i]
                return v_da.v[i]
            end
        end
        return v_da.v[end]
    end
    return v_inst 
end

function make_c_in_func(c_in::CinData, component_index::Int)
    function cin(t)
        if isempty(c_in.t_in)
             # Handle constant inflow (e.g. Col 4 if not updated)
             return c_in.c_in[1][component_index]
        end
        for i in eachindex(c_in.t_in)
            if t ≤ c_in.t_in[i]
                return c_in.c_in[i][component_index]
            end
        end
        return c_in.c_in[end][component_index]
    end
    return cin
end

# --- Analysis ---

const L_COLUMN = 0.08 

results = DataFrame(Column=Int[], Component=String[], Conc_Avg=Float64[], Conc_Std=Float64[], Time_Avg_d=Float64[], Tau_d=Float64[], Cin_Eff=Float64[], Rate_mM_d=Float64[])

# Component mapping: Name -> Index in c_ins
# Updated based on user instruction and prepare_bc_model_m2.jl structure.
# NO3: Index 1 (Ensures Col 4 is 0.0)
# DOC: Index 5
# DIC: Index 7
components = [
    ("NO3", 1), 
    ("DOC", 5),
    ("DIC", 7)
]

# Mapping for all_ds fields
get_data_field(ds, name) = 
    if name == "NO3" ds.no3
    elseif name == "DOC" ds.doc
    elseif name == "DIC" ds.dic
    else nothing end

println("Processing data...")

for c in 1:4
    v_inst = make_v_func(v_da[c])
    
    # Calculate X(t) = integral of v(t)
    X(t) = quadgk(v_inst, 0, t)[1]
    
    # Pre-calculate for interpolation
    # Span enough time. 
    max_t = 30 * 24 * 3600 # 30 days
    dense_t = 0:3600:max_t
    dense_x = [X(t) for t in dense_t]
    T_interp = DataInterpolations.LinearInterpolation(dense_t, dense_x)
    
    for (name, cin_idx) in components
        data = get_data_field(all_ds[c], name)
        
        # Filter NaNs
        valid_idx = .!isnan.(data.conc)
        times = data.t[valid_idx]
        concs = data.conc[valid_idx]
        
        if length(times) < 3
            println("Not enough data for Col $c $name")
            continue
        end
        
        # Last 3 points
        # Ensure sorted by time
        sorted_p = sortperm(times)
        times = times[sorted_p]
        concs = concs[sorted_p]
        
        last_t = times[end-2:end]
        last_c = concs[end-2:end]
        
        # Use Median and corresponding time
        # Create pairs to keep time associated with concentration
        pairs = collect(zip(last_c, last_t))
        # Sort by concentration
        sorted_pairs = sort(pairs, by = first)
        
        # Median index (for 3 points, it's 2)
        mid_idx = div(length(sorted_pairs) + 1, 2)
        med_c, med_t = sorted_pairs[mid_idx]
        
        std_c = std(last_c)
        
        # Residence Time using time of median concentration
        dist_at_t = X(med_t)
        target_x = dist_at_t - L_COLUMN
        
        if target_x < 0
            tau = NaN
            cin_val = NaN
            rate = NaN
        else
            t_in = T_interp(target_x)
            tau = med_t - t_in
            
            if cin_idx != -1
                cin_func = make_c_in_func(c_ins[c], cin_idx)
                cin_val = cin_func(t_in) * 1000 # Convert M to mM
            else
                cin_val = 0.0
            end
            
            # Rate Calculation: R = (C_out - C_in) / tau
            # Units: mM / d
            rate = (med_c - cin_val) / (tau / 86400) 
        end
        
        push!(results, (c, name, med_c, std_c, med_t/86400, tau/86400, cin_val, rate))
    end
end

println("\n--- Raw Rates ---")
display(results)

# Corrected Rates (Using Column 4 as Control)
# Assuming Col 4 represents abiotic/background processes.
# Corrected Rate = Rate_Col - Rate_Col4

println("\n--- Corrected Rates (Col X - Col 4) ---")
corrected_results = DataFrame(Column=Int[], Component=String[], Rate_Corrected=Float64[], Rate_Control=Float64[])

# Get Col 4 rates
col4_rates = Dict()
for r in eachrow(results)
    if r.Column == 4
        col4_rates[r.Component] = r.Rate_mM_d
    end
end

for r in eachrow(results)
    if r.Column != 4
        control_rate = get(col4_rates, r.Component, NaN)
        corrected_rate = r.Rate_mM_d - control_rate
        push!(corrected_results, (r.Column, r.Component, corrected_rate, control_rate))
    end
end

display(corrected_results)

println("\n--- Aggregated Corrected Rates by Species (Min, Med, Max) ---")
summary_stats = combine(groupby(corrected_results, :Component)) do df
    (
        Rate_Min = minimum(df.Rate_Corrected),
        Rate_Med = median(df.Rate_Corrected),
        Rate_Max = maximum(df.Rate_Corrected),
        Rate_Control = first(df.Rate_Control) # All rows in group have same control rate
    )
end

display(summary_stats)

# Optional: Plotting with ErrorBars (as requested implicitly)
# The user asked to "Add errorbars calculation", which is Conc_Std.
# We have it in `results`.
CSV.write("data/no3_doc_dic_rates.csv", summary_stats)