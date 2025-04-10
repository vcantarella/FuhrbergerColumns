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


K_o2 = 3e-9
K_glu_n2 = 50e-3
K_glu_o = 25e-3
K_glu_n = 50e-3
I_o2 = 2e-8
K_no3 = 2e-3
K_no2 = 4e-3
kl_alpha = 9/3600

#Solubility of Henry's law:
H_298 = 1.2e-5
dlnH_dt_1_t = 1700
H_265 = H_298*exp(dlnH_dt_1_t*(1/265.65-1/298))
c_o2 = H_298*21.2e3*1e-3
c_n2 = 6.4e-6*78e3*1e-3
#c_o2 = 0.00027
a_nh4 = 1e-4
R = 8.314*1e3
T = 298
c_co2 = 3.3e-4*1e3*1e-3
print(c_co2)

print(c_o2)

function rhs!(du, u, p, t)
    bio, glu, o2, no3, no2, n2, co2 = @view u[:]
    dbiodt, dgludt, do2dt, dno3dt, dno2dt, dn2dt, dco2dt = @view du[:]

    
    qs_ox = (0.28*0.27+0.0023)./3600
    qs_no3 = (0.27*0.34+0.0034)./3600
    qs_no2 = (0.06*0.26+0.0019)./3600
    
    # Q_ox = co2**5*a_nh4/(glu*o2**(9/2))
    # delta_g_ox = -1999.9+R*T*np.log(Q_ox)
    # Ft_ox = 1-np.exp((delta_g_ox+15)/(R*T))
    # Q_no3 = co2**5*no2**9*a_nh4/(glu*no3**9)
    # delta_g_no3 = -1287.5+R*T*np.log(Q_no3)
    # Ft_no3 = 1-np.exp((delta_g_no3+15)/(R*T))
    # Q_no2 = co2**5*n2**3*a_nh4/(glu*no2**6)
    # delta_g_no2 = -2278.4 +R*T*np.log(Q_no2)
    # Ft_no2 = 1-np.exp((delta_g_no2+15)/(R*T))
    r_ox = qs_ox*bio*(glu/(glu+K_glu_o))*(o2/(o2+K_o2))#*Ft_ox
    r_no3 = qs_no3*bio*(glu/(glu+K_glu_n))*(no3/(no3+K_no3))*(I_o2/(o2+I_o2))
    r_no2 = qs_no2*bio*(glu/(glu+K_glu_n2))*(no2/(no2+K_no2))*(I_o2/(o2+I_o2))
    
    
    #mass transfer for O2:
    #c_o2 = c_o2
    transfer_rate = kl_alpha* (c_o2 - o2)
    
    du[1] = r_ox/0.27- 0.0023/3600/0.27*bio + r_no3/0.34- 0.0034/3600/0.34*bio + r_no2/0.26- 0.0019/3600/0.26*bio
    du[2] = -r_ox - r_no3 - r_no2
    du[3] = -r_ox/0.491 + transfer_rate
    du[4] = -r_no3/0.199
    du[5] = -r_no2/0.393+r_no3/0.199
    du[6] = r_no2/0.393*0.5
    du[7] = r_ox*(0.35/0.27) + r_no3*(0.7/0.34) + r_no2*(0.3/0.26)
end

u0 = [62e-3/24.626, 40e-3,7.2e-3/32, 390/14*1e-3, 1e-12, c_n2, c_co2]
p = [1]
tspan = (0, 36 .*3600)
prob = ODEProblem(rhs!, u0, tspan, p)
sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-8)

kornaros = CSV.read(datadir("exp_raw", "kornaros_fig5.csv"), DataFrame)
t_no3 = kornaros[!, "no3.t"]
t_no3 = convert.(Float64, skipmissing(t_no3))
no3 = kornaros[!, "no3.c"]
no3 = convert.(Float64, skipmissing(no3))
no2 = kornaros[!, "no2.c"]
no2 = convert.(Float64, skipmissing(no2))
t_no2 = kornaros[!, "no2.t"]
t_no2 = convert.(Float64, skipmissing(t_no2))
t_o2 = kornaros[!, "o2.t"]
t_o2 = convert.(Float64, skipmissing(t_o2))
o2 = kornaros[!, "o2.c"]
o2 = convert.(Float64, skipmissing(o2))
bio = kornaros[!, "biomass.c"]
bio = convert.(Float64, skipmissing(bio))
t_bio = kornaros[!, "biomass.t"]
t_bio = convert.(Float64, skipmissing(t_bio))

fig = Figure()
ax = Axis(fig[1, 1], title = "Model vs. Data",
   xlabel = "Time (h)", ylabel = "Concentration (mmol/L)",
   #yticks = (0:0.2:1.3
   )
# axno2 = Axis(fig[1,1],yaxisposition = :right, ygridvisible = false, xgridvisible = false, ylabel = "[O_2]")
lines!(ax, sol.t./3600, [sol.u[i][4] for i in eachindex(sol.u)], label = "NO3-", color = :blue)
lines!(ax, sol.t./3600, [sol.u[i][5] for i in eachindex(sol.u)], label = "NO2-", color = :red)
lines!(ax, sol.t./3600, [sol.u[i][2] for i in eachindex(sol.u)], label = "glu", color = :black)
lines!(ax, sol.t./3600, [sol.u[i][3] for i in eachindex(sol.u)], label = "O2", color = :orange)
lines!(ax, sol.t./3600, [sol.u[i][1] for i in eachindex(sol.u)], label = "Bio", color = :green)
scatter!(ax, t_no3, no3.*1e-3./14, label = "NO3-", color = :blue)
scatter!(ax, t_no2, no2.*1e-3./14, label = "NO2-", color = :red)
scatter!(ax, t_o2, o2.*1e-3./32, label = "O2", color = :orange)
scatter!(ax, t_bio, bio.*1e-3./24.626, label = "Bio", color = :green)
#lines!(ax, sol.t./3600, [reshape(sol.u[i],size(u0))[1,4] for i in eachindex(sol.u)], label = "B inlet", color = :green)
fig
#save("plots/freds_model_23-24_Marc.png", fig)
