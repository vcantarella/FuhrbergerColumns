using DrWatson
@quickactivate "FuhrbergerColumns"
using OrdinaryDiffEq
using CairoMakie
using Colors
using Symbolics
using DataFrames, XLSX, Statistics, CSV
using PRIMA
using LinearSolve
using ForwardDiff
using BlackBoxOptim
using SciMLSensitivity

include(srcdir("freds_model.jl"))


# Geometry of the domain

L = 0.08 #length of the domain [m]
dx = 5e-4 #cell spacing [m]
x = 0:dx:L #domain [m]
r = (3.7/2)/100 #radius of the column [m]
A = pi*r^2 #cross-sectional area of the column [m²]


# Initial and boundary conditions

Q = 3.5/3600*1e-6 #flow rate [m³/s]
q = Q/A #specific flow rate [m/s]

## parameters from tracer test:
ϕ = 0.2412 #porosity [-]
αˡ = 6.5e-3 #longitudinal dispersivity [m]

# model velocity andn dispersivity:
v = q/ϕ #velocity [m/s]
Di = [1e-9 1e-9] #diffusion coefficient [m²/s]
De = [αˡ*v αˡ*v] + Di #dispersion coefficient [m²/s]

## initial conditions:
c_in = [0.7e-3 0.0] # inflow concentration of no3- and no2- [mmol/L]
u0 = zeros(length(x), 4) .+ 1e-16
u0[:,3] .= 62e-3 # initial concentration of edc [mmol/Lw] # water equivalent
u0[:,4] .= 1e-3 # initial biomass concentration [-]

# Model time
# Duration, injection time & conc. and flow rate:--------------------------
t_end  = 200 # Experiment end time [h]
tspan = (0,t_end*3600) # Experiment time span [s]
teval = 360

## model
rhs! = create_fredsmodel(v, De, dx, c_in, 2)

## parameters
k_no3 = 6e-9 # reaction rate of no3- [1/s]
k_no3c = 3.3e-8
k_no2 = 9.9e-8 # reaction rate of no2- [1/s]
k_no2c = 1e-7
K_no3 = 3e-3 # half-saturation constant of no3- [mmol/L]
K_no2 = 1e-3 # half-saturation constant of no2- [mmol/L]
K_pyr = 2e-3 # half-saturation constant of pyr [mmol/L]
K_c = 0.5e-3 # half-saturation constant of c [mmol/L]
c_t = 0.07e-3 # total concentration of no3- [mmol/L]
st = 0.4 # steepness of the activation function [-]
p0 = [k_no3, k_no2, k_no3c, k_no2c, K_no3, K_no2, K_pyr, K_c, c_t, st]
lb = [1e-7, 1e-7, 1e-4, 1e-4, 1e-4, 0.1, 1e-8]
ub = [1e-2, 1e-2, 10.0, 10.0, 0.2, 1.0, 1e-2]
## optimizing the problem and ODE solver
old_prob = ODEProblem(rhs!, u0, tspan, p0)
du0 = zeros(size(u0))
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs!(du, u, p0, 1),
    du0, u0) # add the sparsity pattern to speed up the solution
fixed_rhs! = ODEFunction(rhs!, jac_prototype=jac_sparsity)
fastprob = ODEProblem(fixed_rhs!, u0, tspan, p0)
sol = solve(fastprob, QNDF(), saveat=teval, abstol = 1e-12, reltol = 1e-12)

fig = Figure()
ax = Axis(fig[1, 1], title = "Model vs. Data",
   xlabel = "Time (h)", ylabel = "Concentration (mmol/L)",
   #yticks = 0:0.2:1.3
   )
axno2 = Axis(fig[1,1],yaxisposition = :right, ygridvisible = false, xgridvisible = false, ylabel = "[NO2⁻]")

lines!(ax, sol.t./3600, [sol.u[i][end,1] for i in eachindex(sol.u)], label = "NO3- model", color = :blue)
lines!(axno2, sol.t./3600, [sol.u[i][end,2] for i in eachindex(sol.u)], label = "NO2- model", color = :red)
fig

const col_sizex = size(x)[1]
function flat_rhs!(du, u, p, t)
    @inbounds du = reshape(du, col_sizex, 4)
    @inbounds u = reshape(u, col_sizex, 4)
    rhs!(du, u, p, t)
    @inbounds du = reshape(du, col_sizex*4)
    @inbounds u = reshape(u, col_sizex*4)
    nothing
end

du0f = zeros(col_sizex*4)
u0f = u0[:]
flat_rhs!(du0f, u0f, p0, 1)
rhs!(du0, u0, p0, 1)
du0f
du0f = reshape(du0f, col_sizex, 4)
du0f == du0
tspan = (0.0, t_end*3600)
g1(u, p, t) = u[end, 1]
g2(u, p, t) = u[end, 2]
function dg1du(out, u, p, t, i)
    out .= 0.0
    out[end, 1] = 1.0
end
function dg2du(out, u, p, t, i)
    out = zeros(out)
    out[end, 2] = 1.0
end
res = adjoint_sensitivities(sol, QNDF(), t = teval, dgdu_discrete = dg1du, abstol = 1e-8,
    reltol = 1e-8)
length(res)

sensprob = ODEForwardSensitivityProblem(flat_rhs!, u0f, tspan, p0)
sensol = solve(sensprob, DP8(), saveat=teval, abstol = 1e-8, reltol = 1e-8)
sol_u, dp = extract_local_sensitivities(sensol)

# getting and plotting the sensitivities with respect to model rates (1-4)
k_no3dp = dp[1]
k_no3dp = k_no3dp[[161,322],:]
k_no2dp = dp[2]
k_no2dp = k_no2dp[[161,322],:]
k_no3cdp = dp[3]
k_no3cdp = k_no3cdp[[161,322],:]
k_no2cdp = dp[4]
k_no2cdp = k_no2cdp[[161,322],:]

figsens = Figure()
ax_kno3 = Axis(figsens[1, 1], title = "Sensitivity of k_no3-",
   xlabel = "Time (h)", ylabel = L"\frac{-\partial [NO_{3/2}^-]}{\partial k_{NO_3^-}}",
   #yticks = -0.5:0.2:0.5
   )
ax_kno2 = Axis(figsens[1, 2], title = "Sensitivity of k_no2-",
   xlabel = "Time (h)", ylabel = L"\frac{-\partial [NO_{3/2}^-]}{\partial k_{NO_2^-}}",
   #yticks = -0.5:0.2:0.5
   )
ax_kno3c = Axis(figsens[2, 1], title = "Sensitivity of k_no3c-",
    xlabel = "Time (h)", ylabel = L"\frac{-\partial [NO_{3/2}^-]}{\partial k_{NO_3c^-}}",
    #yticks = -0.5:0.2:0.5
    )
ax_kno2c = Axis(figsens[2, 2], title = "Sensitivity of k_no2c-",
    xlabel = "Time (h)", ylabel = L"\frac{-\partial [NO_{3/2}^-]}{\partial k_{NO_2c^-}}",
    #yticks = -0.5:0.2:0.5
    )
lines!(ax_kno3, sensol.t./3600, -k_no3dp[1,:], label = L"$NO_3^-$ sensitivity", color = :blue)
lines!(ax_kno3, sensol.t./3600, -k_no3dp[2,:], label = L"$NO_2^-$ sensitivity", color = :red)
lines!(ax_kno2, sensol.t./3600, -k_no2dp[1,:], label = L"$NO_3^-$ sensitivity", color = :blue)
lines!(ax_kno2, sensol.t./3600, -k_no2dp[2,:], label = L"$NO_2^-$ sensitivity", color = :red)
lines!(ax_kno3c, sensol.t./3600, -k_no3cdp[1,:], label = L"$NO_3^-$ sensitivity", color = :blue)
lines!(ax_kno3c, sensol.t./3600, -k_no3cdp[2,:], label = L"$NO_2^-$ sensitivity", color = :red)
lines!(ax_kno2c, sensol.t./3600, -k_no2cdp[1,:], label = L"$NO_3^-$ sensitivity", color = :blue)
lines!(ax_kno2c, sensol.t./3600, -k_no2cdp[2,:], label = L"$NO_2^-$ sensitivity", color = :red)
Legend(figsens[3, 1:2], ax_kno3, position = :topright,
    merge = true,
    title = "Sensitivity of k_no3- and k_no2-",
    titlefontsize = 12, fontsize = 10, colorbar = false,
    #horizontaldirection
    orientation = :horizontal,
    )
resize_to_layout!(figsens)
figsens
save("plots/sensitivity_k.png", figsens)

