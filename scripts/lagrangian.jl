"""
    lagrangian.jl

This script performs a Lagrangian simulation of nitrate transport and reduction in sediment columns.
It calculates the trajectory of fluid parcels based on time-variable flow velocities, determines
the residence time for each parcel, and applies a zero-order reaction model to estimate
nitrate concentrations at the outlet.

The script compares the model results with experimental data.

# Dependencies:
- `model_data_structures.jl`: Defines data structures for flow and concentration data.
- `prepare_lab_data_model_m2.jl`: Loads and processes experimental data (velocity, concentrations).
"""

using DrWatson
using CairoMakie
using DataInterpolations
using QuadGK
using Statistics
using CSV
using DataFrames

# --- Configuration & Data Loading ---

# Include definitions of data structures and data preparation scripts
# `prepare_lab_data_model_m2.jl` includes `prepare_bc_model_m2.jl`, which loads:
# - `v_da`: Dictionary of VDataA (Velocity Data with Start/End times)
# - `c_ins`: Dictionary of CinData (Inflow Concentrations)
# - `all_ds`: Dictionary of experimental datasets (concentrations at outlet)
include("model_data_structures.jl")
include("prepare_lab_data_model_m2.jl")

"""
    load_tracer_parameters(path::String)

Loads optimized tracer parameters and averages them for column 4 (if needed).
"""
function load_tracer_parameters(path::String)
    params = load(path)["tracer_params"]
    # Use average porosity/dispersivity of columns 1-3 for column 4 as a fallback/estimate
    params[4] = [mean([params[k][1] for k in 1:3]),
                 mean([params[k][2] for k in 1:3])]
    return params
end

tracer_params = load_tracer_parameters("data/optimized_tracer_params_m2.jld2")

"""
    get_model_data()

Retrieves necessary data structures from the global scope (populated by included scripts).
Returns velocity data, inflow concentration data, and experimental measurements.
"""
function get_model_data()
    # v_da, c_ins, and all_ds are global variables from `prepare_lab_data_model_m2.jl`
    return v_da, c_ins, all_ds
end

# Load data into local variables for clarity
velocity_data, inflow_data, experimental_data = get_model_data()

# --- Model Functions ---

"""
    make_v_func(v_da::VDataA)

Creates a function `v(t)` that returns the instantaneous flow velocity at time `t`.
The velocity is assumed to be piecewise constant.
"""
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

"""
    make_c_in_func(c_in::CinData)

Creates a function `c_in(t)` that returns the inflow nitrate concentration at time `t`.
The concentration is assumed to be piecewise constant based on the experimental schedule.
"""
function make_c_in_func(c_in::CinData)
    function cin(t)
        for i in eachindex(c_in.t_in)
            if t ≤ c_in.t_in[i]
                return c_in.c_in[i][1] # Index 1 corresponds to NO3-
            end
        end
        return c_in.c_in[end][1]
    end
    return cin
end

"""
    calculate_concentration_out(t, c_in_func, residence_time_func, reaction_rate)

Calculates the concentration at the outlet at time `t`.
Model: C_out(t) = C_in(t - τ(t)) + r * τ(t)
where τ(t) is the residence time and r is the zero-order reaction rate.
"""
function calculate_concentration_out(t, c_in_func, residence_time_func, reaction_rate)
    τ = residence_time_func(t)
    # Ensure concentration doesn't go below zero (though simple linear model might allow it)
    # Here we just apply the formula:
    return c_in_func(t - τ) + reaction_rate * τ
end

# --- Simulation & Plotting ---

# Plot settings
colors = [:blue, :orange, :green, :red]

# Figure 1: Nitrate Output and Reaction Rates
fig_conc = Figure(size = (500, 500))
axn = Axis(fig_conc[1, 1], 
    title = "a. Nitrate Outflows",
    titlealign = :left,
    ylabel = "NO₃⁻ [mmol L⁻¹]",
    yticks = 0:0.5:2.5,
    xticks = 5:5:28
)
axr = Axis(fig_conc[2, 1],
    title = "b. Estimated Reaction Rates",
    titlealign = :left,
    ylabel = "Rate [mmol L⁻¹ day⁻¹]",
    xlabel = "Time [days]",
    xticks = 5:5:28,
    limits = (3,27,1,4.3)
)
linkxaxes!(axn, axr)

# Figure 2: Velocity
fig_vel = Figure(size = (700, 400))
axq = Axis(fig_vel[1,1], 
    title = "Velocity",
    titlealign= :left,
    xlabel = "Time [days]", 
    ylabel = "v [m s⁻¹]",
    xticks = 5:5:28
)

# Column length (m)
const L_COLUMN = 0.08 

# Store calculated rates for comparison table
model_rates_collection = []

# Loop over columns 1 to 3
for c in 1:3
    # 1. Setup Velocity and Inflow Functions
    v_inst = make_v_func(velocity_data[c])
    c_in = make_c_in_func(inflow_data[c])

    # 2. Lagrangian Trajectory Calculation
    # X(t) represents the position of a fluid parcel that entered at t=0? 
    # Actually, X(t) here is defined as integral of v from 0 to t.
    X(t) = quadgk(v_inst, 0, t)[1]

    # Pre-calculate X(t) for interpolation to speed up inverse lookup
    # We span the entire experimental duration
    dense_t = 1:(3*3600):(27*24*60*60) # Every 3 hours
    dense_x = [X(t) for t in dense_t]
    
    T_interp = DataInterpolations.LinearInterpolation(dense_t, dense_x)

    # Calculate minimum time before any fluid could have exited (plug flow)
    mint = T_interp(L_COLUMN) 

    # Residence time function τ(t)
    τ(t) = t - T_interp(X(t) - L_COLUMN)

    # 3. Estimate Reaction Rate from Data
    col_data = experimental_data[c]
    no3_exp = col_data.no3
    
    # Filter data points that occur after breakthrough (t > mint)
    valid_indices = findall(t -> t > mint, no3_exp.t)
    valid_t = no3_exp.t[valid_indices]
    valid_conc_exp = no3_exp.conc[valid_indices] .* 1e-3 # Convert mmol/L to mol/L for calculation

    calculated_rates = Float64[]
    for i in eachindex(valid_t)
        t_val = valid_t[i]
        tau_val = τ(t_val)
        c_in_val = c_in(t_val - tau_val)
        c_out_val = valid_conc_exp[i]
        
        # Rate = (C_out - C_in) / tau
        # Rate is in mol/L/s
        r_val = (c_out_val - c_in_val) / tau_val
        push!(calculated_rates, r_val)
    end
    
    mean_rate = isempty(calculated_rates) ? 0.0 : mean(calculated_rates)
    println("Column $c: Mean Rate = $mean_rate mol/L/s")

    # Calculate rate in mol per kg of sand
    # Formula: r_sand = (r_pw * phi) / ((1 - phi) * rho_grain)
    # phi = tracer_params[c][1]
    # rho_grain = 2.65 kg/L
    phi = tracer_params[c][1]
    rho_grain = 2.65
    mean_rate_sand = (-mean_rate * phi * 86400) / ((1 - phi) * rho_grain)
    println("Column $c: Mean Rate = $mean_rate_sand mol/kg_sand/day")
    
    push!(model_rates_collection, (c, -mean_rate, mean_rate_sand))

    # 4. Calculate Model Output with MEAN rate
    analysis_t = dense_t[dense_t .> mint]
    c_out_values = [calculate_concentration_out(t, c_in, τ, mean_rate) for t in analysis_t]
    
    # 5. Plotting
    
    # Concentrations
    scatter!(axn, no3_exp.t ./ (24*60*60), no3_exp.conc, 
        label = "Column $c", color = colors[c], markersize = 8)

    plot_t = collect(analysis_t) ./ (3600*24) # Convert to days
    lines!(axn, plot_t, c_out_values .* 1e3, # Convert to mmol/L for plot
        label = "Column $c (Model)", color = colors[c]) 
        
    if c == 3
        lines!(axn, collect(analysis_t)./(3600*24), c_in.(analysis_t) .* 1e3, 
            linestyle = :dash, label = "Inflow\n concentration", color = :black)
    end

    # Reaction Rates
    # Convert mol/L/s to mmol/L/day: * 1e3 * 86400
    conv_factor = 1e3 * 86400
    scatter!(axr, valid_t ./ (24*60*60), -calculated_rates .* conv_factor, 
        color = colors[c], markersize = 8, label = "Calc. Rate Col $c")
    lines!(axr, valid_t ./ (24*60*60), -calculated_rates .* conv_factor, 
        color = colors[c], label = "Calc. Rate Col $c")    
    hlines!(axr, [-mean_rate * conv_factor], color = colors[c], linestyle = :dash, label = "Mean Rate Col $c")

    # Velocity Plot
    v_plot = v_inst.(analysis_t)
    lines!(axq, plot_t, v_plot, color = colors[c], label = "Column $c")
end

# Finalize Plots

Legend(fig_conc[3,1], axn, framevisible = false,
    nbanks=2, merge = true, orientation = :horizontal, labelsize = 10)

save(plotsdir("nitrate_and_rates.png"), fig_conc)
save(plotsdir("velocity_separate.png"), fig_vel)

# export the rates data
df_export = DataFrame(
    column = [data[1] for data in model_rates_collection],
    rate_mol_L_s = [data[2] for data in model_rates_collection],
    rate_mol_kg_d = [data[3] for data in model_rates_collection],
)

CSV.write("data/no3_rates.csv", df_export)