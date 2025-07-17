using CairoMakie
using DrWatson
using DataFrames
using CSV
using XLSX
using Dates

# load analytical data
file_path = datadir("exp_raw", "no2_analyticall_results.xlsx")
sheet_name = "sampling"
data = XLSX.readtable(file_path, sheet_name)
df = DataFrame(data)


## Events:
# Start time for the t0 of the experiment for columns 1 and 2
t0_1_2 = DateTime(2025, 05, 12, 20, 0, 0)
c_no3_1_2 = 0.6 # concentration of NO3- in the input solution for columns 1 and 2 [mM]
c_lac_1_2 = 0.0 # concentration of lactate in the input solution for columns 1 and 2 [mM]
# Start time for the t0 of the experiment for columns 3 and 4
t0_3_4 = DateTime(2025, 05, 12, 20, 14, 0)
c_no3_3_4 = 0.6 # concentration of NO3- in the input solution for columns 3 and 4 [mM]
c_lac_3_4 = 0.33 # concentration of lactate in the input solution for columns 3 and 4 [mM]
# Moment when we swithed the bag for columns 3 and 4, removing the lactate solution
t_switch_3_4 = DateTime(2025, 05, 16, 16, 35, 0)
c_lac_3_4_switch = 0.0 # concentration of lactate in the input solution for columns 3 and 4 after the switch [mM]
t_s_34_d = Dates.value(t_switch_3_4 - t0_1_2) / (1000 * 60 * 60 * 24)
# Moment when we switch the bag for all columns, increasing the NO3- concentration
t_switch_all_1 = DateTime(2025, 05, 23, 19, 00, 0)
c_no3_all_1 = 1.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
t_1_d = Dates.value(t_switch_all_1 - t0_1_2) / (1000 * 60 * 60 * 24)
# Moment when the increase the NO3- concentration
t_switch_all_2 = DateTime(2025, 05, 28, 19, 15, 0)
c_no3_all_2 = 2.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
t_2_d = Dates.value(t_switch_all_2 - t0_1_2) / (1000 * 60 * 60 * 24)
# Moment when pumps 1 and 2 stopped due to clogging.
t_stop_1 = DateTime(2025,05,28, 19, 50)
t_stop_2 = DateTime(2025,05,28, 20,10)
# moment when the pumps were restarted
t_restart_1_2 = DateTime(2025,05,29, 08, 45)
# Approximate time when all pumps stoppeds due to energie failure in the building
t_stop_all = DateTime(2025,05,29, 12, 30)
t_stop_d = Dates.value(t_stop_all - t0_1_2) / (1000 * 60 * 60 * 24)
# Moment when the pumps were restarted after the energy failure
t_restart_all = DateTime(2025,05,29, 18, 16)
t_restart_d = Dates.value(t_restart_all - t0_1_2) / (1000 * 60 * 60 * 24)
# Moment when we switch to 4 mM inflow concentration
t_switch_4mM = DateTime(2025, 06, 03, 18, 25, 0)
c_no3_4mM = 4.
t_4_d = Dates.value(t_switch_4mM - t0_1_2) / (1000 * 60 * 60 * 24)

# make a df per column
df1 = df[df.column .== 1, :]
df2 = df[df.column .== 2, :]
df3 = df[df.column .== 3, :]
df4 = df[df.column .== 4, :]
dfs = [df1, df2, df3, df4]
# Make a plot for the most importntant data, skip misisng data
colors = [:red, :blue, :green, :orange]
fig = Figure(size = (1200, 800))
ax1 = Axis(fig[1, 1], title = "pH", ylabel = "pH", xlabel = "t (days)")
ax2 = Axis(fig[2, 1], title = "EC", ylabel = "EC [μS/cm]", xlabel = "t (days)")
ax3 = Axis(fig[1, 2], title = "NO2-", ylabel = "NO2- [μM]", xlabel = "t (days)")
ax4 = Axis(fig[2, 2], title = "Q [μL/min]", ylabel = "Q [μL/min]", xlabel = "t (days)")
ax5 = Axis(fig[2, 3], title = "Fe2+", ylabel = "Fe2+ [μM]", xlabel = "t (days)")
ax6 = Axis(fig[1, 3], title = "NO3-", ylabel = "NO3- [μM]", xlabel = "t (days)")
ax7 = Axis(fig[3, 1], title = "SO4-2", ylabel = "SO4-2 [μM]", xlabel = "t (days)")

# Figure from the early time data
fige = Figure()
ax1e = Axis(fige[1, 1], title = "NO₂⁻", ylabel = "conc [μM]", xlabel = "t (days)")
ax2e = Axis(fige[1, 2], title = "NO₃⁻", ylabel = "conc [μM]", xlabel = "t (days)")
ax3e = Axis(fige[2, 1], title = "Fe²⁺", ylabel = "conc [μM]", xlabel = "t (days)")
ax4e = Axis(fige[2, 2], title = "SO₄²⁻", ylabel = "conc [μM]", xlabel = "t (days)")

# Figure from the early time data
figl = Figure()
ax1l = Axis(figl[1, 1], title = "NO₃⁻ and NO₂⁻", ylabel = "conc [μM]", xlabel = "t (days)")
ax2l = Axis(figl[2, 1], ylabel = "conc [μM]", xlabel = "t (days)")

using Statistics
c_so4_in = mean([238.5, 232.2, 234.6, 234.4])/ 96 *1e3

for i in 1:4
    # get the data for the column
    df = dfs[i]
    # get the data for the column
    pH = df[!, "pH"]
    index = findall(!ismissing, pH)
    pH = pH[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax1, t, pH,  color = colors[i])
    scatter!(ax1, t, pH, color = colors[i])
    # get the data for the column
    EC = df[!, "EC"]
    index = findall(!ismissing, EC)
    EC = EC[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax2, t, EC, label = "Column $i", color = colors[i])
    scatter!(ax2, t, EC, label = "Column $i", color = colors[i])
    # get the data for the column
    NO2 = df[!, "NO2-"]
    index = findall(!ismissing, NO2)
    NO2 = NO2[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax3, t, NO2, label = "Column $i", color = colors[i])
    scatter!(ax3, t, NO2, label = "Column $i", color = colors[i])
    # plot in the early time data
    lines!(ax1e, t[t.< 8], NO2[t.< 8], label = "Column $i", color = colors[i])
    scatter!(ax1e, t[t.< 8], NO2[t.< 8], label = "Column $i", color = colors[i])
    lines!(ax1l, t, NO2, label = "Column $i", color = colors[i])
    scatter!(ax1l, t, NO2, label = "Column $i", color = colors[i])
    # plot in the early time data
    # get the data for the column
    NO3 = df[!, "NO3-"]
    index = findall(!ismissing, NO3)
    NO3 = NO3[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax6, t, NO3, label = "Column $i", color = colors[i])
    scatter!(ax6, t, NO3, label = "Column $i", color = colors[i])
    lines!(ax2e, t[t.< 8], NO3[t.< 8], label = "Column $i", color = colors[i])
    scatter!(ax2e, t[t.< 8], NO3[t.< 8], label = "Column $i", color = colors[i])
    lines!(ax1l, t, NO3, label = "Column $i", color = colors[i])
    scatter!(ax1l, t, NO3, label = "Column $i", color = colors[i])
    # get the data for the column
    SO4 = df[!, "SO4-2"]
    index = findall(!ismissing, SO4)
    SO4 = SO4[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax7, t, SO4, label = "Column $i", color = colors[i])
    scatter!(ax7, t, SO4, label = "Column $i", color = colors[i])
    hlines!(ax7, c_so4_in, color = colors[i], linestyle = :dash, label = "Column $i SO4-2 inflow")
    lines!(ax2l, t, SO4, label = "Column $i", color = colors[i])
    scatter!(ax2l, t, SO4, label = "Column $i", color = colors[i])
    hlines!(ax2l, c_so4_in, color = colors[i], linestyle = :dash, label = "Column $i SO4-2 inflow")
    lines!(ax4e, t[t.< 8], SO4[t.< 8], label = "Column $i", color = colors[i])
    scatter!(ax4e, t[t.< 8], SO4[t.< 8], label = "Column $i", color = colors[i])
    hlines!(ax4e, c_so4_in, color = colors[i], linestyle = :dash, label = "Column $i SO4-2 inflow")
    # plot in the early time data
    # get the data for the column
    Q = df[!, "Q"]
    index = findall(!ismissing, Q)
    Q = Q[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax4, t, Q, label = "Column $i", color = colors[i])
    scatter!(ax4, t, Q, label = "Column $i", color = colors[i])
    # Fe data
    Fe = df[!, "Fe2+"]
    index = findall(!ismissing, Fe)
    Fe = Fe[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax5, t, Fe, label = "Column $i", color = colors[i])
    scatter!(ax5, t, Fe, label = "Column $i", color = colors[i])
    # plot in the early time data
    lines!(ax3e, t[t.< 8], Fe[t.< 8], label = "Column $i", color = colors[i])
    scatter!(ax3e, t[t.< 8], Fe[t.< 8], label = "Column $i", color = colors[i])
end
for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax1l, ax2l]
    # add the vertical lines for the events
    vlines!(ax, t_1_d, color = :black, label = "1mM NO3- switch")
    vlines!(ax, t_2_d, color = :darkblue, label = "2mM NO3- switch")
    vlines!(ax, t_s_34_d, color = :crimson, label = "Switch lactate off in columns 3 and 4")
    vlines!(ax, t_stop_d, color = :brown, label = "Pumps stopped due to clogging")
    vlines!(ax, t_restart_d, color = :lightgreen, label = "Pumps restarted after clogging")
    vlines!(ax, t_4_d, color = :purple, label = "4mM NO3- switch")
end
for ax in [ax1e, ax2e, ax3e, ax4e]
    vlines!(ax, t_s_34_d, color = :crimson, label = "Switch lactate off in columns 3 and 4")
end
ylims!(ax7, 2000, 2700)
ylims!(ax2l, 2000, 2700)
ylims!(ax4e, 2400, 2700)
Legend(fig[3, 2:3], ax2, merge = true,
framevisible = false, nbanks = 2)
# adjust the layout
CairoMakie.resize_to_layout!(fig)
fig
save("lab_data.png", fig)
Legend(fige[3, 1:2], ax4e, merge = true,
framevisible = false, nbanks = 2)
fige
save("lab_data_early.png", fige)
linkxaxes!(ax1l, ax2l)
figl
Legend(figl[1:2, 2], ax2l, merge = true,
framevisible = false, nbanks = 1)
save("lab_data_early_lactate.png", figl)