using CairoMakie
using StaticArrays
import OrdinaryDiffEq as ODE
using SparseConnectivityTracer
using Symbolics
using LinearSolve
using SciMLSensitivity
using FiniteDiff
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
    affect!(integrator) = integrator.u[1,1] = cs[i]
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
sol_monod = ODE.solve(fastprob_monod, ODE.QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
    callback = cbset,
    tstops = tstops,
    )
sol_zero = ODE.solve(fastprob_zero, ODE.QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
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
data_points_no3 = [sol_monod(t)[end, 1] for t in data_points_t]

function f(p)
    u0_monod = zeros(eltype(p),length(x)+1, 2) # 5 mobile components + 2 immobile components (active and inactive biomass)
    u0_monod[:,2] .= p[6]
    c_in = [2e-3, 0]
    u0_monod[1, :] .= c_in #initial inflow concentration 2 mM NO3-
    fprob = remake(fastprob_monod, u0 = u0_monod, p = p)
    sol_monod = ODE.solve(fprob, ODE.QNDF(linsolve = KLUFactorization()), abstol = 1e-8, reltol = 1e-8,
    callback = cbset,
    #tstops = tstops,
    saveat = data_points_t,
    )
    no3_results = zeros(eltype(p), length(data_points_t))
    for i in eachindex(data_points_t)
        ind = findfirst(sol_monod.t .== data_points_t[i])
        no3_results[i] = sol_monod.u[ind][end, 1]
    end
    residuals = data_points_no3 .- no3_results
    return residuals
end

# Define the jacobian via finite differences
jac(p) = FiniteDiff.finite_difference_jacobian(f, p)

p = vcat(p_monod,[1e-5])

f(p)
jac(p)

f_log_p(p) = f(exp.(p))
jac_log_p(p) = FiniteDiff.finite_difference_jacobian(f_log_p, p)

f_log_p(log.(p))
jac_log_p(log.(p))