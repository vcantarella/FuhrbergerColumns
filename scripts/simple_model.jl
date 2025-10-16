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
tracer_params = load("data/optimized_tracer_params_m2.jld2") # Load the optimized tracer parameters
tracer_params = tracer_params["tracer_params"]
# Load the prepared data
function prepare_data()
    include("prepare_lab_data_model_m2.jl")
    return v_interp, c_ins, all_ds
end

prepare_data()
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
fig_height = 300
fig = Figure(size = (800, fig_height*3))
for c in 1:3
# Starting the model for column 1
x = range(0+dx/2, stop=L-dx/2, step=dx)  # Spatial locations
rhs! = reactive_transport_builder(v_interp[c], Deff, dx, tracer_params[c][1],
    tracer_params[c][2], 2.65)

p0 = [
    2.8e-8, # r_s (steady_state rate for NO3-)
]
u0 = zeros(length(x)+1, 1) # 5 mobile components + 2 immobile components (active and inactive biomass)
du0 = copy(zeros(size(u0))) # Initialize the derivative array
c_indata = c_ins[c]
u0[1, :] .= c_indata.c_in[1][1] # initial inflow concentration 2 mM NO3-
rhs!(du0, u0, p0, 0.0) # Calculate the initial derivative
c_indata.t_in
# Create a callback to update the inflow concentration at specific times
switch_callback_condition(u, t, integrator) = t == c_indata.t_in[1] # t switch to 1.5 mM
affect!(integrator) = integrator.u[1,1] = c_indata.c_in[2][1]
cb = ODE.DiscreteCallback(switch_callback_condition, affect!)
switch_callback_condition2(u, t, integrator) = t == c_indata.t_in[2] # t switch to 1.5 mM
affect2!(integrator) = integrator.u[1,1] = c_indata.c_in[3][1]
cb2 = ODE.DiscreteCallback(switch_callback_condition2, affect2!)
cbset = ODE.CallbackSet(cb, cb2)

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
sol = ODE.solve(fastprob, QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
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
col1 = all_ds[1]
no3 = col1.no3

# Plot results

axn = Axis(fig[c, 1], title = "Outflow concentrations",
    xlabel = "Time (days)", ylabel = "Concentration (M)",
    yticks = 0:2e-4:1e-3)
# axf = Axis(fig[2, 1],
#     xlabel = "Time (days)", ylabel = "Concentration (M)")
# axs = Axis(fig[3, 1],
#     xlabel = "Time (days)", ylabel = "Concentration (M)",
#     yticks = 1e-3:5e-4:3e-3)
# ylims!(axs, 9e-4, 3e-3)
plot_t = sol.t ./ (24*60*60) # convert seconds to days
lines!(axn, plot_t, no3_out, label = "NO3- outflow", color = :blue)
# lines!(axs, plot_t, tracer_out, label = "NO3- tracer outflow", color = :blue, linestyle = :dash)
# lines!(axn, plot_t, no2_out, label = "NO2- outflow", color = :orange)
# lines!(axs, plot_t, so4_out, label = "SO4-2 outflow", color = :green)
# lines!(axf, plot_t, fe_out, label = "Fe+2 outflow", color = :purple)
# lines!(ax, plot_t, lac_out, label = "Lactate outflow", color = :red)
# scatter!(axn, no2.t ./ (24*60*60), no2.conc*1e-6, label = "Measured NO2- outflow", color = :orange, markersize = 8)
scatter!(axn, no3.t ./ (24*60*60), no3.conc/14*1e-3, label = "Measured NO3- outflow", color = :blue, markersize = 8)
# scatter!(axs, so4.t ./ (24*60*60), so4.conc*1e-6, label = "Measured SO4-2 outflow", color = :green, markersize = 8)
# scatter!(axf, fe.t ./ (24*60*60), fe.conc*1e-6, label = "Measured Fe outflow", color = :purple, markersize = 8)
end
plots_in_fig = AbstractPlot[]
labels_in_fig = AbstractString[]
# for ax in [axn, axf, axs]
#     pl, lb = Makie.get_labeled_plots(ax, merge=false, unique=false)
#     append!(plots_in_fig, pl)
#     append!(labels_in_fig, lb)
# end

ulabels = Base.unique(labels_in_fig)
mergedplots = [[lp for (i, lp) in enumerate(plots_in_fig) if labels_in_fig[i] == ul]
        for ul in ulabels]

Legend(fig[:, 2], mergedplots, ulabels, framevisible=false)
# linkxaxes!(axn, axf, axs)
resize_to_layout!(fig)
fig
save("outflow_concentrations_m2.png", fig, px_per_unit = 2.0)
