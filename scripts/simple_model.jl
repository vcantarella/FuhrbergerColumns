using DrWatson
using CairoMakie
using DataFrames, XLSX, Statistics
using Dates
using StaticArrays
import OrdinaryDiffEq as ODE
using Optimization
using SciMLSensitivity
using SparseConnectivityTracer
using Symbolics
using LinearSolve
using DataInterpolations
# Data structures used in the model building and simulation
include("model_data_structures.jl")
include("prepare_lab_data_model_m2.jl")
tracer_params = load("data/optimized_tracer_params_m2.jld2") # Load the optimized tracer parameters
tracer_params = tracer_params["tracer_params"]
tracer_params[4] = [mean([tracer_params[k][1] for k in 1:3]),
                    mean([tracer_params[k][2] for k in 1:3])] # use the average porosity of columns 1 and 2 for column 4
# Load the prepared data
function prepare_data()
    return v_interp, c_ins, all_ds
end

v_interp, c_ins, all_ds = prepare_data()
# Diffusion coefficients for the different (mobile) components
Deff = @SVector [
    1.0e-9, # NO3-
]
"""
    reactive_transport_builder(v_data::VData, cin_data::CinData, Deff, dx, ϕ, ρₛ, αₗ)
Builds the reactive transport model for the given flow data and inflow concentration data.
    (column specific)
    # Arguments:
    - `interp`: Interpolation function for the flow velocity data.
    - `Deff`: Diffusion coefficients for the different mobile species.
    - `dx`: Spatial step size for the model 1D grid.
    - `ϕ`: Porosity of the column.
    - `ρₛ`: Bulk density of the column material. 2.65 g/cm³ (sand)
    - `αₗ`: Longitudinal dispersivity of the column.
    # Returns:
    - `rhs!(du, u, p, t)`: Function that computes the right-hand side of the ODE system.
    According to DifferentialEquations.jl solver requirements.
"""
function reactive_transport_builder(interp, Deff, dx, ϕ, αₗ, ρₛ)
    De = MVector{length(Deff)}(Deff)
        function rhs!(du, u, p, t)
            @inline v = interp(t)
            De .= Deff .+ αₗ * v
            # unpack the state variables
            c_in = @view u[1, :] # first row is the inflow concentration
            #no3_ = @view u[2:end,1]
            

            # unpack the parameters
            r_s, = p

            n_rows = size(u, 1) # number of spatial rows excluding the inflow row
            nmob = size(De, 1)  # Number of mobile components
            # transport
            # Calculate transport terms directly without temporary arrays
            @inbounds for j in 1:nmob
                # First cell (boundary condition)
                du[2,j] = -v * (u[2,j] - c_in[1,j]) / dx
                
                # Calculate dispersion at first cell - only forward gradient
                grad_fwd = (u[3,j] - u[2,j]) / dx
                du[2,j] += De[j] * grad_fwd / dx  # Remove the gradient difference

                # Interior cells
                for i in 3:n_rows-1
                    # Advection
                    du[i,j] = -v * (u[i,j] - u[i-1,j]) / dx
                    
                    # Dispersion
                    grad_fwd = (u[i+1,j] - u[i,j]) / dx
                    grad_bwd = (u[i,j] - u[i-1,j]) / dx
                    du[i,j] += De[j] * (grad_fwd - grad_bwd) / dx
                end
                
                # Last cell
                du[n_rows,j] = -v * (u[n_rows,j] - u[n_rows-1,j]) / dx
                grad_bwd = (u[n_rows,j] - u[n_rows-1,j]) / dx
                du[n_rows,j] += De[j] * (0.0 - grad_bwd) / dx  # Zero-gradient at boundary
            end
            @inbounds for k in 2:n_rows
                r_no3 = ifelse(u[k, 1] > 0, r_s, 0.0) # constant term for NO3-
                       # Update state variables
                du[k,1] -= r_no3
            end
            #make sure du[1, :] = 0
            du[1, :] .= 0.0
        end
    return rhs!
end

dx = 0.0001 # Spatial step size
L = 0.08 #m (8 cm)  # Spatial locations
fig_height = 300*3
fig = Figure(size = (1000, fig_height))
axn = Axis(fig[1:2, 1], title = "Nitrate Outflows",
    xlabel = "Time (days)", ylabel = "(NO₃⁻) [mmol L⁻¹]",
    yticks = 0:5e-1:2.1,
    xticks = 5:5:28
    )
axdoc = Axis(fig[1, 2], title = "DOC Outflows",
    xlabel = "Time (days)", ylabel = "DOC [mmol L⁻¹]",
    #yticks = 0:50:300
    )
axdic = Axis(fig[2, 2], title = "DIC Outflows",
    xlabel = "Time (days)", ylabel = "DIC [mmol L⁻¹]",
    #yticks = 0:50:300
    )
axso4 = Axis(fig[3, 2], title = "Sulfate Outflows",
    xlabel = "Time (days)", ylabel = "(SO₄²⁻) [mmol L⁻¹]",
    #yticks = 0:2e-1:1.1
    )
axno2 = Axis(fig[3, 1], title = "Nitrite Outflows",
    xlabel = "Time (days)", ylabel = "(NO₂⁻) [mmol L⁻¹]",
    #yticks = 0:2e-1:1.1
    )

colors = [:blue, :orange, :green, :red]
for c in 1:4
# Starting the model for column 1
x = range(0+dx/2, stop=L-dx/2, step=dx)  # Spatial locations
rhs! = reactive_transport_builder(v_interp[c], Deff, dx, tracer_params[c][1],
    tracer_params[c][2], 2.65)
r0 = if c == 1
    2.7e-8
elseif c == 2
    3.3e-8
else
    3.2e-8
end

p0 = [
    r0, # r_s (steady_state rate for NO3-)
]
u0 = zeros(length(x)+1, 1) # 5 mobile components + 2 immobile components (active and inactive biomass)
du0 = copy(zeros(size(u0))) # Initialize the derivative array
c_indata = c_ins[c]

u0[1, :] .= c_indata.c_in[1][1] # initial inflow concentration 2 mM NO3-
rhs!(du0, u0, p0, 0.0) # Calculate the initial derivative
c_indata.t_in
## Create a callback set based on c_indata
cbs = ODE.DiscreteCallback[]
for i in 1:length(c_indata.t_in)
    t_switch = c_indata.t_in[i]
    c_new = c_indata.c_in[i+1][1]
    condition(u, t, integrator) = t == t_switch
    affect!(integrator) = integrator.u[1,1] = c_new
    push!(cbs, ODE.DiscreteCallback(condition, affect!))
end
cbset = ODE.CallbackSet(cbs...)

# @time rhs!(du0, u0, p0, 0.0) # Benchmark the RHS function
# using BenchmarkTools
# @benchmark rhs!($du0, $u0, $p0, 0.0) # Benchmark the RHS function
# @code_warntype rhs!(du0, u0, p0, 0.0) # Profile the RHS function

du0
tspan = (0.0, 27*24*60*60) # 27 days in seconds
old_prob = ODE.ODEProblem(rhs!, u0, tspan, p0)

# Use Symbolics for sparsity detection - handles ifelse properly
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs!(du, u, p0, 100000.0),
    du0, u0)
fixed_rhs! = ODE.ODEFunction(rhs!, jac_prototype=jac_sparsity)
fastprob = ODE.ODEProblem(fixed_rhs!, u0, tspan, p0)

# defining points to stop
tstops = c_indata.t_in
sol = ODE.solve(fastprob, ODE.QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
    callback = cbset,
    tstops = tstops,
    )
sol.t
# Check model outflow:
no3_out = [sol.u[i][end, 1] for i in eachindex(sol.t)]
# no2_out = [sol.u[i][end, 2] for i in eachindex(sol.t)]
# so4_out = [sol.u[i][end, 3] for i in eachindex(sol.t)]
# fe_out = [sol.u[i][end, 4] for i in eachindex(sol.t)]
# lac_out = [sol.u[i][end, 5] for i in eachindex(sol.t)]
# tracer_out = [sol.u[i][end, 6] for i in eachindex(sol.t)]


# check the outflow data
col1 = all_ds[c]
no3 = col1.no3
#no3_std = col1.no3_std

# Plot results


# axf = Axis(fig[2, 1],
#     xlabel = "Time (days)", ylabel = "Concentration (M)")
# axs = Axis(fig[3, 1],
#     xlabel = "Time (days)", ylabel = "Concentration (M)",
#     yticks = 1e-3:5e-4:3e-3)
# ylims!(axs, 9e-4, 3e-3)
plot_t = sol.t ./ (24*60*60) # convert seconds to days
lines!(axn, plot_t, no3_out*1e3, label = "Column $c", color = colors[c])
# lines!(axs, plot_t, tracer_out, label = "NO3- tracer outflow", color = :blue, linestyle = :dash)
# lines!(axn, plot_t, no2_out, label = "NO2- outflow", color = :orange)
# lines!(axs, plot_t, so4_out, label = "SO4-2 outflow", color = :green)
# lines!(axf, plot_t, fe_out, label = "Fe+2 outflow", color = :purple)
# lines!(ax, plot_t, lac_out, label = "Lactate outflow", color = :red)
# scatter!(axn, no2.t ./ (24*60*60), no2.conc*1e-6, label = "Measured NO2- outflow", color = :orange, markersize = 8)
scatter!(axn, no3.t ./ (24*60*60), no3.conc, label = "Column $c", color = colors[c], markersize = 8)
#errorbars!(axn, no3.t ./ (24*60*60), no3.conc/14, no3_std.conc ./14; color = colors[c])
# scatter!(axs, so4.t ./ (24*60*60), so4.conc*1e-6, label = "Measured SO4-2 outflow", color = :green, markersize = 8)
# scatter!(axf, fe.t ./ (24*60*60), fe.conc*1e-6, label = "Measured Fe outflow", color = :purple, markersize = 8)
scatter!(axdoc, col1.doc.t ./ (24*60*60), col1.doc.conc, label = "Column $c", color = colors[c], markersize = 8)
lines!(axdoc, col1.doc.t ./ (24*60*60), col1.doc.conc, label = "Column $c", color = colors[c], linestyle = :dash)
scatter!(axdic, col1.dic.t ./ (24*60*60), col1.dic.conc, label = "Column $c", color = colors[c], markersize = 8)
lines!(axdic, col1.dic.t ./ (24*60*60), col1.dic.conc, label = "Column $c", color = colors[c], linestyle = :dash)
scatter!(axso4, col1.so4.t ./ (24*60*60), col1.so4.conc, label = "Column $c", color = colors[c], markersize = 8)
lines!(axso4, col1.so4.t ./ (24*60*60), col1.so4.conc, label = "Column $c", color = colors[c], linestyle = :dash)
scatter!(axno2, col1.no2.t ./ (24*60*60), col1.no2.conc, label = "Column $c", color = colors[c], markersize = 8)
lines!(axno2, col1.no2.t ./ (24*60*60), col1.no2.conc, label = "Column $c", color = colors[c], linestyle = :dash)
end
fig
# Plot the inflow concentration in time
c_in_plot = []
c_indata = c_ins[1] # using column 1 inflow data (all are very similar)
plot_t = 0:0.1:27
for t in plot_t
    # find the last switch time before t
    idx = findlast(c_indata.t_in .<= t*24*60*60)
    if idx === nothing
        c_in_loc = c_no3_bck1
    else
        c_in_loc = c_indata.c_in[idx+1][1]*1e3
    end
    push!(c_in_plot, c_in_loc)
end
lines!(axn, plot_t, c_in_plot, label = "Inflow concentration", color = :black, linestyle = :dash)
# inflow concentration for the remaining plots
lines!(axdoc, plot_t, zeros(length(plot_t)), label = "Inflow concentration", color = :black, linestyle = :dash)
lines!(axdic, plot_t, ones(length(plot_t)).*30/12, label = "Inflow concentration", color = :black, linestyle = :dash)
lines!(axso4, plot_t, ones(length(plot_t)).*0.74/96, label = "Inflow concentration", color = :black, linestyle = :dash)
lines!(axno2, plot_t, zeros(length(plot_t)), label = "Inflow concentration", color = :black, linestyle = :dash)
xlims!(axn, (4.0, 27.0))
xlims!(axdoc, (4.0, 27.0))
xlims!(axdic, (4.0, 27.0))
xlims!(axso4, (4.0, 27.0))
xlims!(axno2, (4.0, 27.0))
# plots_in_fig = AbstractPlot[]
# labels_in_fig = AbstractString[]
# # for ax in [axn, axf, axs]
# #     pl, lb = Makie.get_labeled_plots(ax, merge=false, unique=false)
# #     append!(plots_in_fig, pl)
# #     append!(labels_in_fig, lb)
# # end

# ulabels = Base.unique(labels_in_fig)
# mergedplots = [[lp for (i, lp) in enumerate(plots_in_fig) if labels_in_fig[i] == ul]
#         for ul in ulabels]

Legend(fig[4, :], axn, framevisible=false, merge=true, orientation = :horizontal)
# linkxaxes!(axn, axf, axs)
resize_to_layout!(fig)
fig
save("outflow_concentrations_m2.png", fig, px_per_unit = 2.0)
