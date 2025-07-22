using DrWatson
using CairoMakie
using DataFrames, XLSX, Statistics
using Dates
using StaticArrays
using OrdinaryDiffEq
using Optimization
using SciMLSensitivity
using SparseConnectivityTracer
using Symbolics
using LinearSolve
using DataInterpolations
# Data structures used in the model building and simulation
include("model_data_structures.jl")
tracer_params = load("data/optimized_tracer_params.jld2") # Load the optimized tracer parameters

# Load the prepared data
function prepare_data()
    include("prepare_lab_data_model.jl")
    return v_ds, c_ins, all_ds
end

prepare_data()
# Diffusion coefficients for the different (mobile) components
Deff = @SVector [
    1.0e-9, # NO3-
    1.0e-9, # NO2-
    1.0e-9, # SO4-2
    1.0e-9, # Fe+2
    1.0e-9, # Lactate
    1.0e-9, # NO3- (tracer)
]
"""
    reactive_transport_builder(v_data::VData, cin_data::CinData, Deff, dx, ϕ, ρₛ, αₗ)
Builds the reactive transport model for the given flow data and inflow concentration data.
    (column specific)
    # Arguments:
    - `v_data::VData`: Flow velocity data for the column.
    - `cin_data::CinData`: Inflow concentration data for the column.
    - `Deff`: Diffusion coefficients for the different mobile species.
    - `dx`: Spatial step size for the model 1D grid.
    - `ϕ`: Porosity of the column.
    - `ρₛ`: Bulk density of the column material. 2.65 g/cm³ (sand)
    - `αₗ`: Longitudinal dispersivity of the column.
    # Returns:
    - `rhs!(du, u, p, t)`: Function that computes the right-hand side of the ODE system.
    According to DifferentialEquations.jl solver requirements.
"""
function reactive_transport_builder(v_data::VData, cin_data::CinData, Deff, dx, ϕ, αₗ, ρₛ)
    vs = SVector{length(v_data.v)}(v_data.v)
    # v = mean(vs) # average velocity
    Des = MVector{length(Deff)}(Deff)
    end_times = SVector{length(v_data.end_times)}(v_data.end_times)
    c_ins = SVector{length(cin_data.c_in)}(cin_data.c_in)
    t_ins = SVector{length(cin_data.t_in)}(cin_data.t_in)
    function dynamic_transport(t, vs=vs, end_times=end_times)::Float64
        @inbounds for i in eachindex(vs)
            if t <= end_times[i]
                v = vs[i]
                return v   # This should be inside the if block
            end
        end
        # Add a fallback return outside the loop
        return vs[end]
    end
    interp = SmoothedConstantInterpolation(vs, end_times, d_max=300,extrapolation = ExtrapolationType.Constant)

    function dynamic_c_in(t, c_ins=c_ins, t_ins=t_ins)
        @inbounds for i in eachindex(t_ins)
            if t <= t_ins[i]
                return c_ins[i]
            end
        end
        return c_ins[end]  # Return the last concentration if no match found
    end
    function build_rhs(dynamic_transport, dynamic_c_in, dx, Deff, αₗ, De=Des)
        function rhs!(du, u, p, t)
            @inline v = interp(t)
            De .= Deff .+ αₗ * v
            # unpack the state variables
            no3_ = @view u[:,1]
            no2_ = @view u[:,2]
            so4_ = @view u[:,3]
            fe2_ = @view u[:,4]
            lac_ = @view u[:,5]
            no3t_ = @view u[:,6]  # NO3- tracer
            b = @view u[:,7]  # biomass active fraction (immobile))
            b_s = @view u[:,8]  # biomass sulfate reducers (immobile)

            # unpack the parameters
            μ₁, μ₂, αₑ, k_dec, Ya1, Ya2, Yd1, Yd2, Ka1, Ka2, Kb,
            r_so4_max, μₛ, Ks, Ys, Yds, I_no3, I_no2, f = p

            n_rows = size(u, 1)
            @inline c_in = dynamic_c_in(t)
            nmob = size(De, 1)  # Number of mobile components
            # transport
            # Calculate transport terms directly without temporary arrays
            @inbounds for j in 1:nmob
                # First cell (boundary condition)
                du[1,j] = -v * (u[1,j] - c_in[j]) / dx
                
                # Calculate dispersion at first cell - only forward gradient
                grad_fwd = (u[2,j] - u[1,j]) / dx
                du[1,j] += De[j] * grad_fwd / dx  # Remove the gradient difference
                
                # Interior cells
                for i in 2:n_rows-1
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
            @inbounds for k in 1:n_rows
                # Calculate each term once to avoid repeated computation
                γ = αₑ/(b[k]+Kb)*(μ₁/Yd1*no3_[k]/(Ka1 + no3_[k]) + μ₂/Yd2*no2_[k]/(Ka2 + no2_[k]))^(-1)
                γ = ifelse(γ ≥ 0.99, .99, γ)  # Limit γ to 1
                γₛ = αₑ/(b_s[k]+Kb)*(μₛ/Yds*so4_[k]/(Ks + so4_[k])*I_no2/(I_no2 + no2_[k])*I_no3/(I_no3 + no3_[k]))^(-1)
                γₛ = ifelse(γₛ ≥ 0.99, .99, γₛ)  # Limit γₛ to 1
                r_no3b = μ₁/Ya1 * no3_[k] / (Ka1 + no3_[k]) * γ
                r_no2b = μ₂/Ya2 * no2_[k] / (Ka2 + no2_[k]) * γ
                r_b = ((μ₁ * no3_[k] / (Ka1 + no3_[k]) + 
                    μ₂ * no2_[k] / (Ka2 + no2_[k])) * γ - k_dec) * b[k]  # Decay term
                so4_max = 2.6e-3
                r_so4 = r_so4_max * (1 - so4_[k]/so4_max) # Sulfate dissolution term
                r_fe = (1/7 * r_no3b * b[k] + 3/14 * r_no2b * b[k])*f  # Iron dissolution term
                r_so4n = (2/7 * r_no3b * b[k] + 3/7 * r_no2b * b[k])*f  # Sulfate dissolution term
                r_s = μₛ/Ys * so4_[k] / (Ks + so4_[k]) * I_no2/(I_no2 + no2_[k]) * I_no3/(I_no3 + no3_[k]) * γₛ  # Sulfate reduction term
                r_b_s = (r_s-k_dec) * b_s[k]  # Sulfate reduction term for biomass
                # Update state variables
                du[k,1] -= r_no3b*b[k]
                du[k,2] += (r_no3b - r_no2b)*b[k]
                du[k,3] += r_so4 - r_s * b_s[k] + r_so4n  # sulfate concentration
                du[k,4] += r_fe
                # du[k,5] -= mult_term*r_no2b*b
                du[k,7] = r_b  # biomass active fraction
                du[k,8] = r_b_s  # biomass inactive fraction
            end
        end
    return rhs!
    end
    return build_rhs(dynamic_transport, dynamic_c_in, dx, Deff, αₗ)
end

dx = 0.0005 # Spatial step size
L = 0.08 #m (8 cm)  # Spatial locations


# Starting the model for column 1
x = range(0+dx/2, stop=L-dx/2, step=dx)  # Spatial locations
rhs! = reactive_transport_builder(v_ds[1], c_ins[1], Deff, dx, tracer_params[1][1],
    tracer_params[1][2], 2.65)

γa = 1/25 # stoichiometric coefficient of e-acceptor in the anabolic rctieaction
γc1 = 4/2 # stoichiometric coefficient of e-acceptor in the catabolic rea
γc2 = 4/3 # stoichiometric coefficient of e-acceptor in the catabolic rea
γcs = 1/2
Yd = 0.3 # yield of biomass from the catabolic reaction
p0 = [
    3e-5, # μ₁ (growth rate for NO3-)
    3e-5, # μ₂ (growth rate for NO2-)
    4e-7, # αₑ (effective diffusion coefficient)
    1.15e-6, # k_dec (decay rate)
    1/((1/Yd-1)*γc1), # Ya1 (yield of biomass from NO3-)
    1/((1/Yd-1)*γc2), # Ya2 (yield of biomass from NO2-)
    Yd, # Yd1 (yield of biomass from NO3-)
    Yd, # Yd2 (yield of biomass from NO2-)
    5e-4, # Ka1 (half-saturation constant for NO3-)
    5e-4, # Ka2 (half-saturation constant for NO2-)
    1e-3, # Ks (half-saturation constant for sulfate)
    1e-6, # k_s (decay rate for sulfate)
    8e-6, # μs (growth rate for sulfate reducers)
    1e-6, # Ks (half-saturation constant for sulfate)
    1/((1/Yd-1)*γcs), # Ys (yield of biomass from sulfate)
    Yd, # Yd3 (yield of biomass from sulfate)
    5e-6, # I_no2 (inhibition constant for NO2-)
    5e-6, # I_no3 (inhibition constant for NO3-)
    0.2, # f (fraction of fes2 electron donor)
]

u0 = zeros(length(x), 8) # 5 mobile components + 2 immobile components (active and inactive biomass)
u0[:,7] .= 1e-4 # Initial concentration of inactive biomass
u0[:, 8] .= 1e-4 # Initial concentration of inactive biomass
du0 = copy(zeros(size(u0))) # Initialize the derivative array

rhs!(du0, u0, p0, 0.0) # Calculate the initial derivative
# @time rhs!(du0, u0, p0, 0.0) # Benchmark the RHS function
# using BenchmarkTools
# @benchmark rhs!($du0, $u0, $p0, 0.0) # Benchmark the RHS function
# @code_warntype rhs!(du0, u0, p0, 0.0) # Profile the RHS function

du0
tspan = (0.0, 32*24*60*60) # 32 days in seconds
old_prob = ODEProblem(rhs!, u0, tspan, p0)
# Symbolics for sparsity detection
detector = TracerSparsityDetector()
jac_sparsity2 = ADTypes.jacobian_sparsity((du, u) -> rhs!(du, u, p0, 100000),
    du0, u0, detector) # add the sparsity pattern to speed up the solution
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs!(du, u, p0, 100000),
    du0, u0) # add the sparsity pattern to speed up the
jac_sparsity == jac_sparsity2
fixed_rhs! = ODEFunction(rhs!, jac_prototype=jac_sparsity2)
fastprob = ODEProblem(fixed_rhs!, u0, tspan, p0)

# defining points to stop
tstops = vcat(v_ds[1].end_times, c_ins[1].t_in)
sort!(tstops) # sort the times
sol = solve(fastprob, FBDF(), abstol = 1e-14, reltol = 1e-18,
    maxiters = 100000,
    tstops = tstops,
    )
sol.t
# Check model outflow:
no3_out = [sol.u[i][end, 1] for i in eachindex(sol.t)]
no2_out = [sol.u[i][end, 2] for i in eachindex(sol.t)]
so4_out = [sol.u[i][end, 3] for i in eachindex(sol.t)]
fe_out = [sol.u[i][end, 4] for i in eachindex(sol.t)]
lac_out = [sol.u[i][end, 5] for i in eachindex(sol.t)]
tracer_out = [sol.u[i][end, 6] for i in eachindex(sol.t)]


# check the outflow data
col2 = all_ds[2]
no2 = col2.no2
no3 = col2.no3
so4 = col2.so4
fe = col2.fe
# Plot results
fig = Figure()
axn = Axis(fig[1, 1], title = "Outflow concentrations",
    xlabel = "Time (days)", ylabel = "Concentration (M)",
    yticks = 0:2e-4:1e-3)
axf = Axis(fig[2, 1],
    xlabel = "Time (days)", ylabel = "Concentration (M)")
axs = Axis(fig[3, 1],
    xlabel = "Time (days)", ylabel = "Concentration (M)",
    yticks = 1e-3:5e-4:3e-3)
ylims!(axs, 9e-4, 3e-3)
plot_t = sol.t ./ (24*60*60) # convert seconds to days
lines!(axn, plot_t, no3_out, label = "NO3- outflow", color = :blue)
lines!(axs, plot_t, tracer_out, label = "NO3- tracer outflow", color = :blue, linestyle = :dash)
lines!(axn, plot_t, no2_out, label = "NO2- outflow", color = :orange)
lines!(axs, plot_t, so4_out, label = "SO4-2 outflow", color = :green)
lines!(axf, plot_t, fe_out, label = "Fe+2 outflow", color = :purple)
# lines!(ax, plot_t, lac_out, label = "Lactate outflow", color = :red)
scatter!(axn, no2.t ./ (24*60*60), no2.conc*1e-6, label = "Measured NO2- outflow", color = :orange, markersize = 8)
scatter!(axn, no3.t ./ (24*60*60), no3.conc*1e-6, label = "Measured NO3- outflow", color = :blue, markersize = 8)
scatter!(axs, so4.t ./ (24*60*60), so4.conc*1e-6, label = "Measured SO4-2 outflow", color = :green, markersize = 8)
scatter!(axf, fe.t ./ (24*60*60), fe.conc*1e-6, label = "Measured Fe outflow", color = :purple, markersize = 8)
plots_in_fig = AbstractPlot[]
labels_in_fig = AbstractString[]
for ax in [axn, axf, axs]
    pl, lb = Makie.get_labeled_plots(ax, merge=false, unique=false)
    append!(plots_in_fig, pl)
    append!(labels_in_fig, lb)
end

ulabels = Base.unique(labels_in_fig)
mergedplots = [[lp for (i, lp) in enumerate(plots_in_fig) if labels_in_fig[i] == ul]
        for ul in ulabels]

Legend(fig[:, 2], mergedplots, ulabels, framevisible=false)
linkxaxes!(axn, axf, axs)
fig
save("outflow_concentrations.png", fig, px_per_unit = 2.0)

# Plot alive vs inactive biomass
fig2 = Figure()
ax2 = Axis(fig2[1, 1], title = "Biomass concentrations",
    xlabel = "length (m)", ylabel = "Concentration (M)")
plot_x = x ./ 0.01 # convert cm to m
lines!(ax2, plot_x, sol.u[end][:, 7], label = "Active biomass", color = :blue)
lines!(ax2, plot_x, sol.u[end][:, 8], label = "sULFUR biomass", color = :orange)
#lines!(ax2, plot_x, sol.u[end][:, 7], label = "Inactive biomass", color = :orange)
axislegend(ax2, position = :rc)
fig2