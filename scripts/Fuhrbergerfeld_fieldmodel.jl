using DrWatson
@quickactivate "FuhrbergerColumns"
using CSV
using DataFrames
using CairoMakie
using OrdinaryDiffEq
using NonlinearSolve
using LinearAlgebra
df = CSV.read(datadir("exp_raw","fuhrberger_field/batch_proxy_fuhrberger.csv"), DataFrame)
t_samp = df.t_yrs
ca_samp = df."Ca+2"
hco3_samp = df."HCO3-"
so4_samp = df."SO4-2"
pH_samp = 7.27 # mean of measurements
H_act = 10^(-pH_samp)
K_cal = 10^(-8.48 + 10.3)
K_gyp = 10^-4.6
P_co2 = 10^-2
K_cal = 10^-6

function reactions!(du, u, p, t)
    ca = u[1]
    hco3 = u[2]
    so4 = u[3]

    k_gyp = p[1] # gypsum diss rate
    k_cal = p[2] # calcite diss rate
    k_cal_p = p[3] # calcite precipitation rate
    k_r = p[4] # rate of respiration (zero order hco3 production model)

    # gypsum dissolution model
    iap_gyp = ca * so4
    r_gyp = k_gyp * (1 - iap_gyp / K_gyp)

    # calcite dissolution/precipitation model
    iap_cal = (ca*hco3^2)/P_co2
    z = tanh(1-iap_cal/K_cal)
    if z > 0
        r_cal = k_cal*z
    else
        r_cal = k_cal_p*z
    end

    # dcdts
    du[1] = r_gyp + r_cal
    du[2] = r_cal + k_r
    du[3] = r_gyp
end


u0 = [ca_samp[1], hco3_samp[1], so4_samp[1]]
p0 = [3e-5, 7e-5, 1e-5, 1e-5]
tspan = (0.0, maximum(t_samp))

prob = ODEProblem(reactions!, y0, tspan, p0)
sol = solve(prob, Tsit5(), saveat=t_samp, abstol=1e-12, reltol=1e-12)

fig = Figure()
ax = Axis(fig[1, 1])

lines!(ax, sol.t, sol[1, :], linestyle=:dash, color=:lightblue, label="Ca+2")
lines!(ax, sol.t, sol[2, :], linestyle=:dash, color=:darkblue, label="HCO3-")
lines!(ax, sol.t, sol[3, :], linestyle=:dash, color=:darkred, label="SO4-2")

scatter!(ax, t_samp, ca_samp, color=:lightblue, marker=:diamond)
scatter!(ax, t_samp, hco3_samp, color=:darkblue, marker=:diamond)
scatter!(ax, t_samp, so4_samp, color=:darkred, marker=:diamond)
axislegend(ax, position=:rt)
fig

# Calculating residuals
function residuals(p, params)
    prob0 = params[1]
    prob = remake(prob0, p=p)
    sol = solve(prob, Tsit5(), saveat=t_samp, abstol=1e-12, reltol=1e-12)
    res = [sol[1, :] - params[2], sol[2, :] - params[3], sol[3, :] - params[4]]
    res = vcat(res...)
end

# Fit with NonlinearSolve

params = (prob, ca_samp, hco3_samp, so4_samp)

nlls_prob = NonlinearLeastSquaresProblem(residuals, p0, params)

res = solve(nlls_prob, TrustRegion(); maxiters = 1000, show_trace = Val(true),
    trace_level = TraceWithJacobianConditionNumber(25))

newp = res.u

prob = remake(prob, p=newp)
sol = solve(prob, Tsit5(), saveat=t_samp, abstol=1e-12, reltol=1e-12)
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (yrs)", ylabel = "Concentration (mmol/L)",
title = "Can Sulphate be explained by just gypsum dissolution (interacting with calcite)?")
lines!(ax, sol.t, sol[1, :], linestyle=:dash, color=:lightblue, label="Ca+2")
lines!(ax, sol.t, sol[2, :], linestyle=:dash, color=:darkblue, label="HCO3-")
lines!(ax, sol.t, sol[3, :], linestyle=:dash, color=:darkred, label="SO4-2")

scatter!(ax, t_samp, ca_samp, color=:lightblue, marker=:diamond)
scatter!(ax, t_samp, hco3_samp, color=:darkblue, marker=:diamond)
scatter!(ax, t_samp, so4_samp, color=:darkred, marker=:diamond)
axislegend(ax, position=:lt)

fig
save(plotsdir("fuhrberger_field_model.png"), fig)