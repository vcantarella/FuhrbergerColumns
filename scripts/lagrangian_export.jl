"""
    lagrangian_export.jl

This script performs a Lagrangian simulation of nitrate transport and reduction in sediment columns.
It is designed to be shared with collaborators and loads all necessary data from CSV files.

### CSV Files Explanation:
1. `data/velocity_data_export.csv`: 
   Contains piecewise constant flow velocities for each column.
   - `column`: Column index (1-3).
   - `end_time`: Time (seconds since experiment start) when the current velocity interval ends.
   - `velocity`: Instantaneous flow velocity (m/s).

2. `data/inflow_data_export.csv`:
   Contains piecewise constant inflow concentrations.
   - `column`: Column index (1-3).
   - `switch_time`: Time (seconds since experiment start) when the concentration switches to the next value. 
     The last entry for each column is `Inf`.
   - `NO3`: Inflow nitrate concentration (mM).

3. `data/experimental_data_export.csv`:
   Contains measured outlet concentrations for model comparison.
   - `column`: Column index (1-3).
   - `time`: Measurement time (seconds since experiment start).
   - `NO3`: Measured nitrate concentration at the outlet (mM).
"""

using CairoMakie
using DataInterpolations
using QuadGK
using DataFrames
using CSV
using Statistics

# --- Data Structures ---

struct VelocityData
    v::Vector{Float64}        # Vector of velocities (m/s)
    end_times::Vector{Float64} # End times for each velocity interval (s)
end

struct InflowData
    c_in::Vector{Float64}      # Nitrate concentrations (mM)
    t_in::Vector{Float64}      # Switch times (s)
end

struct ExperimentalData
    t::Vector{Float64}         # Measurement times (s)
    conc::Vector{Float64}      # Measured NO3- (mM)
end

# --- Data Loading ---

function load_exported_data()
    # 1. Load Velocity Data
    df_v = CSV.read("data/velocity_data_export.csv", DataFrame)
    velocity_data = Dict{Int, VelocityData}()
    for col in unique(df_v.column)
        sub = filter(r -> r.column == col, df_v)
        velocity_data[col] = VelocityData(sub.velocity, sub.end_time)
    end

    # 2. Load Inflow Data
    df_in = CSV.read("data/inflow_data_export.csv", DataFrame)
    inflow_data = Dict{Int, InflowData}()
    for col in unique(df_in.column)
        sub = filter(r -> r.column == col, df_in)
        # Switch times are all entries except the last one (which is Inf)
        t_switch = filter(x -> x != Inf, sub.switch_time)
        inflow_data[col] = InflowData(sub.NO3, t_switch)
    end

    # 3. Load Experimental Data
    df_exp = CSV.read("data/experimental_data_export.csv", DataFrame)
    experimental_data = Dict{Int, ExperimentalData}()
    for col in unique(df_exp.column)
        sub = filter(r -> r.column == col, df_exp)
        experimental_data[col] = ExperimentalData(sub.time, sub.NO3)
    end

    return velocity_data, inflow_data, experimental_data
end

# Reconstruct data from CSV
velocity_data, inflow_data, experimental_data = load_exported_data()

# --- Model Functions ---

function make_v_func(vd::VelocityData)
    function v_inst(t)
        for i in eachindex(vd.end_times)
            if t ≤ vd.end_times[i]
                return vd.v[i]
            end
        end
        return vd.v[end]
    end
    return v_inst 
end

function make_c_in_func(id::InflowData)
    function cin(t)
        for i in eachindex(id.t_in)
            if t ≤ id.t_in[i]
                return id.c_in[i]
            end
        end
        return id.c_in[end]
    end
    return cin
end

function calculate_concentration_out(t, c_in_func, residence_time_func, reaction_rate)
    τ = residence_time_func(t)
    # C_out(t) = C_in(t - τ) + r * τ
    return c_in_func(t - τ) + reaction_rate * τ
end

# --- Simulation & Plotting ---

colors = [:blue, :orange, :green, :red]
fig = Figure(size = (800, 800))

axn = Axis(fig[1:2, 1], 
    title = "Nitrate Outflow Simulation (Lagrangian Model)",
    ylabel = "NO₃⁻ [mmol L⁻¹]",
    yticks = 0:0.5:2.5
)

axq = Axis(fig[3,1], 
    xlabel = "Time [days]", 
    ylabel = "Velocity [m s⁻¹]"
)

const L_COLUMN = 0.08 
# Simplified zero-order reaction rates [mmol L-1 s-1]
reaction_rates = Dict(1 => -2.7e-8, 2 => -3.3e-8, 3 => -3.2e-8)

for c in 1:3
    # Initialize functions
    v_inst = make_v_func(velocity_data[c])
    c_in = make_c_in_func(inflow_data[c])

    # Lagrangian Trajectory
    X(t) = quadgk(v_inst, 0, t)[1]
    
    # Speed up residence time lookup using interpolation
    dense_t = 1:(3*3600):(27*24*60*60) # Every 3 hours
    dense_x = [X(t) for t in dense_t]
    T_interp = DataInterpolations.LinearInterpolation(dense_t, dense_x)
    
    mint = T_interp(L_COLUMN) 
    τ(t) = t - T_interp(X(t) - L_COLUMN)

    # Evaluate model
    r0 = reaction_rates[c]
    analysis_t = dense_t[dense_t .> mint]
    c_out_values = [calculate_concentration_out(t, c_in, τ, r0) for t in analysis_t]
    
    # --- Plotting ---
    plot_t_days = collect(analysis_t) ./ (3600*24)
    
    # Model Line
    lines!(axn, plot_t_days, c_out_values .* 1e3, label = "Column $c Outflow", color = colors[c])
    
    # Experimental points
    exp = experimental_data[c]
    scatter!(axn, exp.t ./ (24*60*60), exp.conc, color = colors[c], markersize = 8)

    # Inflow (use column 3 as reference)
    if c == 3
        lines!(axn, plot_t_days, [c_in(t)*1e3 for t in analysis_t], 
            linestyle = :dash, color = :black, alpha = 0.5, label="Inflow concentration")
    end

    # Velocity
    lines!(axq, plot_t_days, [v_inst(t) for t in analysis_t], color = colors[c])
end

linkxaxes!(axn, axq)
Legend(fig[4,1], axn, framevisible = false, orientation = :horizontal)

save("plots/lagrangian_simulation_simplified.png", fig)
println("Analysis complete. Result saved to plots/lagrangian_simulation_simplified.png")