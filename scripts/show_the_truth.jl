using CairoMakie
using StaticArrays
import OrdinaryDiffEq as ODE
using SparseConnectivityTracer
using Symbolics
using LinearSolve
using SciMLSensitivity
# using FiniteDiff
# using DiffEqBayes
# using Turing
# using Distributions
include("model_monod.jl")


dx = 0.0001 # Spatial step size
L = 0.08 #m (8 cm)  # Spatial locations
Q = 2 #ml/hr (convert to m/s)
Q = Q / 3600 * 1e-6 # convert from ml/hr to m3/s
D = 3.5*1e-2 #cm to m diameter of the column
A = π * D^2 / 4 # Cross-sectional area
q = Q/A
ϕ = 0.22 # porosity
αₗ = 1e-3 # dispersivity (m)
v = q/ϕ
De = [1e-9, 0] .+ v*αₗ
x = range(0+dx/2, stop=L-dx/2, step=dx)  # Spatial locations
p_monod = [
    1e-5,
    1e-6,
    1e-6,
    7.7e-4,
    0.22,
]
r0 = 3.3e-8
p0 = [
    r0, # r_s (steady_state rate for NO3-)
]
u0_monod = zeros(length(x)+1, 2) # 5 mobile components + 2 immobile components (active and inactive biomass)
u0_monod[:,2] .= 1e-5
u0_0 = zeros(length(x)+1, 1) # 5 mobile components + 2 immobile components (active and inactive biomass)
du0_monod = copy(zeros(size(u0_monod))) # Initialize the derivative array
du0_0 = copy(zeros(size(u0_0))) # Initialize the derivative array
c_in = [2e-3, 0]
u0_monod[1, :] .= c_in # initial inflow concentration 2 mM NO3-
u0_0[1,:] .= 2e-3
monod! = monod_model(v, De, dx, 1)
zero! = zero_order_model(v, De, dx, 1)

ts = [18, 25] .* 86400 # switch times
cs = [1.5, 1.0] .* 1e-3 #switch concentrations
## Create a callback set based on c_indata
cbs = ODE.DiscreteCallback[]
for i in 1:length(ts)
    t_switch = ts[i]
    condition(u, t, integrator) = t == t_switch
    affect!(integrator) = integrator.u[1,1] = convert(eltype(integrator.u),cs[i])
    push!(cbs, ODE.DiscreteCallback(condition, affect!))
end
cbset = ODE.CallbackSet(cbs...)

tspan = (0.0, 27*24*60*60) # 27 days in seconds
old_prob_monod = ODE.ODEProblem(monod!, u0_monod, tspan, p_monod)
old_prob_zero = ODE.ODEProblem(zero!, u0_0, tspan, p0)

# Use Symbolics for sparsity detection - handles ifelse properly
jac_sparsity_monod = Symbolics.jacobian_sparsity((du, u) -> monod!(du, u, p_monod, 100000.0),
    du0_monod, u0_monod)
fixed_rhs_monod! = ODE.ODEFunction(monod!, jac_prototype=jac_sparsity_monod)
fastprob_monod = ODE.ODEProblem(fixed_rhs_monod!, u0_monod, tspan, p_monod)

jac_sparsity_zero = Symbolics.jacobian_sparsity((du, u) -> zero!(du, u, p0, 100000.0),
    du0_0, u0_0)
fixed_rhs_zero! = ODE.ODEFunction(zero!, jac_prototype=jac_sparsity_zero)
fastprob_zero = ODE.ODEProblem(fixed_rhs_zero!, u0_0, tspan, p0)

# defining points to stop
tstops = ts # a copy
sol_monod = ODE.solve(fastprob_monod, ODE.FBDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
    callback = cbset,
    tstops = tstops,
    )
sol_zero = ODE.solve(fastprob_zero, ODE.FBDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
    callback = cbset,
    tstops = tstops,
    )

# Check model outflow:
no3_monod = [sol_monod.u[i][end, 1] for i in eachindex(sol_monod.t)]
no3_zero = [sol_zero.u[i][end, 1] for i in eachindex(sol_zero.t)]


fig_height = 300
fig = Figure(size = (450, fig_height))
axn = Axis(fig[1, 1], title = "Nitrate Outflows",
    xlabel = "Time (days)", ylabel = "(NO₃⁻) [mmol L⁻¹]",
    yticks = 0:5e-1:2.1,
    xticks = 5:5:28
    )
# plot nitrate concentrations from both models to visually compare
lines!(axn, sol_monod.t ./ 86400, no3_monod .* 1e3, color = :red, label = "Monod - based")
lines!(axn, sol_zero.t ./ 86400, no3_zero .* 1e3, color = :green, label = "Zero - order")
Legend(fig[2, :], axn, framevisible=false, merge=true, orientation = :horizontal)
fig

data_points_t = collect(4:27) .* 86400
data_points_no3 = [sol_zero(t)[end, 1] for t in data_points_t]


sol = ODE.solve(fastprob_monod, ODE.QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
        callback = cbset,
        saveat = data_points_t,
        tstops = tstops,
        sensealg = ForwardDiffSensitivity(;convert_tspan = true)
    )
# find the idexes in sol.t that match data_points_t
index = [findfirst(t .== sol.t) for t in data_points_t]


function monod_prediction_model(p_log)
    # p_log contains log-transformed parameters: 5 for monod, 1 for initial biomass
    p_fit = exp.(p_log[1:5])
    b0_fit = exp(p_log[6])

    # Set up initial conditions with the new biomass value
    u0_fit = zeros(eltype(b0_fit), length(x)+1, 2)
    u0_fit[:,2] .= b0_fit
    u0_fit[1, :] .= [2e-3, 0] # initial inflow concentration

    # Remake and solve the problem
    fprob = remake(fastprob_monod, p = p_fit, u0 = u0_fit)
    sol = ODE.solve(fprob, ODE.QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
        callback = cbset,
        saveat = data_points_t,
        tstops = tstops,
        sensealg = ForwardDiffSensitivity(;convert_tspan = true)
    )

    # Handle solver failures
    if sol.retcode != :Success
        return fill(Inf, length(data_points_t)) # Return infinite error if solve fails
    end

    predicted_no3 = [u[end, 1] for u in sol.u[index]]
    return predicted_no3 .- data_points_no3
end

using ForwardDiff
# calculate the jacobian with ForwardDiff

# direct adjoint calculation investigation:
# discrete adjoint gradient:
function dg(out, u, p, t, i)
    out .= -data_points_no3[i]+u[end, 1]
end 
sol = ODE.solve(fastprob_monod, ODE.QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
        callback = cbset,
        tstops = tstops,
    )
res = adjoint_sensitivities(sol, ODE.Vern9(), t = data_points_t,
    dgdu_discrete = dg, abstol = 1e-8,
    reltol = 1e-8,
    sensealg = GaussAdjoint(checkpointing=true),
    )

# Run the optimization
p_initial = vcat(p_monod, [1e-5]) # Add initial guess for biomass
p_log_initial = log.(p_initial)

jac(p) = ForwardDiff.jacobian((p_log) -> monod_prediction_model(p_log), p)
using nonlinearlstr

fit = nonlinearlstr.lm_trust_region(monod_prediction_model, jac, p_log_initial)

p_fit_log = fit[1]
p_fit_final = exp.(p_fit_log)

println("Fitted Monod parameters: ", p_fit_final)