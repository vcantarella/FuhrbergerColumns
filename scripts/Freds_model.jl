using DrWatson
@quickactivate "FuhrbergerColumns"
using OrdinaryDiffEq
using CairoMakie
using Colors
using Symbolics
using DataFrames, XLSX, Statistics, CSV
using PRIMA
using LinearSolve
using ModelingToolkit
using Sundials
using ForwardDiff
using BlackBoxOptim

include(srcdir("freds_model.jl"))


# Geometry of the domain

L = 0.121 #length of the domain [m]
dx = 5e-4 #cell spacing [m]
x = 0:dx:L #domain [m]
r = (3.7/2)/100 #radius of the column [m]
A = pi*r^2 #cross-sectional area of the column [m²]


# Initial and boundary conditions

Q = 48.9/60*1e-9 #flow rate [m³/s]
q = Q/A #specific flow rate [m/s]

## parameters from tracer test:
ϕ = 0.2412 #porosity [-]
αˡ = 6.5e-3 #longitudinal dispersivity [m]

# model velocity andn dispersivity:
v = q/ϕ #velocity [m/s]
Di = [1e-9 1e-9] #diffusion coefficient [m²/s]
De = [αˡ*v αˡ*v] + Di #dispersion coefficient [m²/s]

## initial conditions:
c_in = [1.3 0.0] # inflow concentration of no3- and no2- [mmol/L]
u0 = zeros(4, length(x))
u0[:,3] .= 1000.0 # initial concentration of edc [mmol/Lw] # water equivalent
u0[:,4] .= 1e-2 # initial biomass concentration [-]

# Model time
# Duration, injection time & conc. and flow rate:--------------------------
t_end  = 181.5 # Experiment end time [h]
tspan = (0,t_end*3600) # Experiment time span [s]
teval = 3600

## model
rhs! = create_fredsmodel(v, De, dx, c_in, 2)

## parameters
k_no3 = 1.2e-3 # reaction rate of no3- [1/s]
k_no2 = 1e-2 # reaction rate of no2- [1/s]
K_no3 = 1e-2 # half-saturation constant of no3- [mmol/L]
K_no2 = 1e-1 # half-saturation constant of no2- [mmol/L]
c_t = 0.5 # total concentration of no3- [mmol/L]
st = 0.45 # steepness of the activation function [-]
k_act = 1e-5 # activation rate [1/s]
p0 = [k_no3, k_no2, K_no3, K_no2, c_t, st, k_act]
lb = [1e-7, 1e-7, 1e-4, 1e-4, 1e-4, 0.1, 1e-8]
ub = [1e-2, 1e-2, 10.0, 10.0, 0.2, 1.0, 1e-2]
## optimizing the problem and ODE solver
old_prob = ODEProblem(rhs!, u0, tspan, p0)
de = complete(modelingtoolkitize(old_prob))
fastprob = ODEProblem(de, [], tspan, jac=true, sparse=true)
sol = solve(fastprob, CVODE_BDF(linear_solver=:GMRES), saveat=teval, abstol = 1e-8, reltol = 1e-8)


# Load the data
df23_24 = CSV.File(datadir("exp_pro","freds_processed_data_23-24.csv")) |> DataFrame
df27_28 = CSV.File(datadir("exp_pro","freds_processed_data_27-28.csv")) |> DataFrame

# lets start with the first data set
sheet = "23-24"
sheet_data = df23_24
t = sheet_data[!, "h"]
no3 = sheet_data[!,"NO3-"]
no3_idx = findall(!ismissing, no3)
no3 = convert.(Float64,no3[no3_idx])
t_no3 = t[no3_idx]
no2 = sheet_data[!,"NO2-"]*1e-3
no2_idx = findall(!ismissing, no2)
no2 = convert.(Float64,no2[no2_idx])
t_no2 = t[no2_idx]

function cost_2324(p)
    cost_prob = remake(fastprob, p=p)
    sol = solve(cost_prob, CVODE_BDF(linear_solver=:GMRES), saveat=t.*3600, abstol = 1e-10, reltol = 1e-10)
    return sum(abs2, [reshape(sol.u[i],size(u0))[end,1] for i in eachindex(sol.u)][no3_idx] .- no3) +
     sum(abs2, [reshape(sol.u[i],size(u0))[end,2] for i in eachindex(sol.u)][no2_idx] .- no2.*1e-3)*1000
end

cost_2324(p0)

res = PRIMA.bobyqa(cost_2324, p0, xl = lb, xu = ub, rhobeg = 5e-9,
   iprint = PRIMA.MSG_RHO)

cost_2324(res[1])
p = res[1]
sol = solve(remake(fastprob, p=p), CVODE_BDF(linear_solver=:GMRES), saveat=t.*3600, abstol = 1e-10, reltol = 1e-10)
# Plot the results:
fig = Figure()
ax = Axis(fig[1, 1], title = "Model vs. Data",
   xlabel = "Time (h)", ylabel = "Concentration (mmol/L)",
   yticks = 0:0.2:1.3
   )
axno2 = Axis(fig[1,1],yaxisposition = :right, ygridvisible = false, xgridvisible = false, ylabel = "[NO2⁻]")

lines!(ax, sol.t./3600, [reshape(sol.u[i],size(u0))[end,1] for i in eachindex(sol.u)], label = "NO3- model", color = :blue)
lines!(axno2, sol.t./3600, [reshape(sol.u[i],size(u0))[end,2] for i in eachindex(sol.u)], label = "NO2- model", color = :red)
scatter!(ax, t_no3, no3, label = "NO3- data", color = :blue)
scatter!(axno2, t_no2, no2, label = "NO2- data", color = :red)
lines!(ax, sol.t./3600, [reshape(sol.u[i],size(u0))[1,4] for i in eachindex(sol.u)], label = "B inlet", color = :green)
fig