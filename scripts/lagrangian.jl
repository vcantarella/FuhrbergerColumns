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
fig_height = 300 * 3
fig = Figure(size = (1000, fig_height))

# Axis for Nitrate Output
axn = Axis(fig[1:2, 1], 
    title = "a. Nitrate Outflows",
    titlealign = :left,
    ylabel = "NO₃⁻ [mmol L⁻¹]",
    yticks = 0:5e-1:2.1,
    xticks = 5:5:28
)

# Axis for Velocity
axq = Axis(fig[3,1], 
    title = "b. Velocity",
    titlealign= :left,
    xlabel = "Time [days]", 
    ylabel = "v [m s⁻¹]",
    xticks = 5:5:28
)

# Column length (m)
const L_COLUMN = 0.08 

# Reaction rates (approximate zero-order rates for each column)
# These values seem to be fitted or estimated previously.
reaction_rates = Dict(
    1 => -2.7e-8,
    2 => -3.3e-8,
    3 => -3.2e-8
)

# Loop over columns 1 to 3
for c in 1:3
    # 1. Setup Velocity and Inflow Functions
    v_inst = make_v_func(velocity_data[c])
    c_in = make_c_in_func(inflow_data[c])

    # 2. Lagrangian Trajectory Calculation
    # X(t) represents the position of a fluid parcel that entered at t=0? 
    # Actually, X(t) here is defined as integral of v from 0 to t.
    # This represents the total distance traveled by a parcel introduced at t=0 by time t.
    # Or more accurately, it's the cumulative displacement field.
    X(t) = quadgk(v_inst, 0, t)[1]

    # Pre-calculate X(t) for interpolation to speed up inverse lookup
    # We span the entire experimental duration
    dense_t = 1:(3*3600):(27*24*60*60) # Every 3 hours
    dense_x = [X(t) for t in dense_t]
    
    # T(x) is the inverse function: given a distance x, when does the "cumulative flow" reach it?
    # This allows us to find when a parcel reaching L at time t must have entered.
    # Wait, the residence time logic below is: τ(t) = t - T(X(t) - L)
    # X(t) is total distance "flow" has moved since t=0.
    # X(t) - L is the "position" in the cumulative flow frame that is L meters behind the current front.
    # T(X(t) - L) gives the time t_in when the cumulative flow was X(t) - L.
    # So a parcel entering at t_in is now at X(t) - X(t_in) = L. Correct.
    T_interp = DataInterpolations.LinearInterpolation(dense_t, dense_x)

    # Calculate minimum time before any fluid could have exited (plug flow)
    mint = T_interp(L_COLUMN) 

    # Residence time function τ(t)
    # t is the current time (observation time at outlet)
    # t_in = T_interp(X(t) - L) is the time the parcel currently at outlet entered the column.
    τ(t) = t - T_interp(X(t) - L_COLUMN)

    # 3. Calculate Model Output
    r0 = reaction_rates[c]
    
    # Define points to evaluate (only after breakthrough)
    analysis_t = dense_t[dense_t .> mint]
    
    # Calculate concentrations
    c_out_values = [calculate_concentration_out(t, c_in, τ, r0) for t in analysis_t]
    
    # 4. Plotting
    # Experimental Data
    col_data = experimental_data[c]
    no3_exp = col_data.no3
    scatter!(axn, no3_exp.t ./ (24*60*60), no3_exp.conc, 
        label = "Column $c", color = colors[c], markersize = 8)

    # Model Output
    plot_t = collect(analysis_t) ./ (3600*24) # Convert to days
    lines!(axn, plot_t, c_out_values .* 1e3, # Convert to mmol/L or relevant scale? Result in to be mol/L.
        label = "Column $c", color = colors[c]) # Code had *1e3, assuming plot wants µM? 
    if c == 3
        lines!(axn, collect(analysis_t)./(3600*24), c_in.(analysis_t) .* 1e3, 
            linestyle = :dash, label = "Inflow concentration", color = :black)
    end

    # Velocity Plot
    v_plot = v_inst.(analysis_t)
    lines!(axq, plot_t, v_plot, color = colors[c])
end

# Finalize Plot
linkxaxes!(axn, axq)
Legend(fig[4,1], axn, framevisible = false, merge = true, orientation = :horizontal)

fig
save(plotsdir("zero_order_fit.png"), fig)