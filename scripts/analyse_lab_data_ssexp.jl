using CairoMakie
using DrWatson
using DataFrames
using CSV
using XLSX
using Dates
using Statistics

# load analytical data
file_path = datadir("exp_raw", "ssexp_data.xlsx")
dfs = DataFrame[]
for col in ["p1", "p2", "p3", "p4"]
    sheet_name = "general_samples_$col"
    data = XLSX.readtable(file_path, sheet_name)
    df_temp = DataFrame(data)
    df_temp[!, "column"] .= parse(Int, last(col))  # Add a column indicating the column number
    push!(dfs, df_temp)
end
t0 = DateTime(2025, 9, 18, 18, 15)  # Start time for columns 2-4. One was already functioning.
# Columns of each df
# sample: the sample ID
# time_in: start of sampling
# time_out: end of sampling
# calculate: t (unadjusted) [days]: time since start of experiment in days
for df in dfs
    df[!, "avg_time"] .= df[!, "time_in"] .+ (df[!, "time_out"] .- df[!, "time_in"]) ./ 2
    df[!, "t (unadjusted) [days]"] .= Dates.value.(df[!, "avg_time"] .- t0) / (1000 * 60 * 60 * 24)
end

# read the NO2- data and merge it with the main dataframes
df_no2 = DataFrame(XLSX.readtable(file_path, "standard_curve_no2"))
# for each sample in df_no2, find the corresponding sample in each df and add the NO2- value
for df in dfs
    df[!, "NO2-"] .= missing
    for (i, sample) in enumerate(df_no2[!, "Sample"])
        for j in 1:nrow(df)
            if df[j, "Sample"] == sample
                df[j, "no2_mmol_L"] = df_no2[i, "no2- [micromol/L]"] / 1000
            end
        end
    end
end
## Events:
# Start time for the t0 of the experiment for columns 1 and 2

# make a df per column

# Make a plot for the most importntant data, skip misisng data
colors = [:red, :blue, :green, :orange]
fig = Figure(size = (1200, 800))
ax1 = Axis(fig[1, 1], title = "pH", ylabel = "pH", xlabel = "t (days)")
ax2 = Axis(fig[2, 1], title = "EC", ylabel = "EC [μS/cm]", xlabel = "t (days)")
ax3 = Axis(fig[1, 2], title = "NO₂⁻", ylabel = "NO₂⁻ [μM]", xlabel = "t (days)")
ax4 = Axis(fig[2, 2], title = "Q [ml/hr]", ylabel = "Q [ml/hr]", xlabel = "t (days)")
ax5 = Axis(fig[2, 3], title = "Fe²⁺", ylabel = "Fe²⁺ [μM]", xlabel = "t (days)")
ax6 = Axis(fig[1, 3], title = "NO₃⁻", ylabel = "NO₃⁻ [μM]", xlabel = "t (days)")
ax7 = Axis(fig[3, 1], title = "SO₄²⁻", ylabel = "SO₄²⁻ [μM]", xlabel = "t (days)")

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
    NO2 = df[!, "no2_mmol_L"]
    index = findall(!ismissing, NO2)
    NO2 = NO2[index]*1000
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # # plot the data
    lines!(ax3, t, NO2, label = "Column $i", color = colors[i])
    scatter!(ax3, t, NO2, label = "Column $i", color = colors[i])
    # # plot in the early time data
    # lines!(ax1e, t[t.< 8], NO2[t.< 8], label = "Column $i", color = colors[i])
    # scatter!(ax1e, t[t.< 8], NO2[t.< 8], label = "Column $i", color = colors[i])
    # lines!(ax1l, t, NO2, label = "Column $i", color = colors[i])
    # scatter!(ax1l, t, NO2, label = "Column $i", color = colors[i])
    # plot in the early time data
    # get the data for the column
    NO3 = df[!, "conc_mgN_L"]
    index = findall(!ismissing, NO3)
    NO3 = NO3[index]
    NO3 .= NO3 ./14*1000  # convert mgN/L to μM
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
    # SO4 = df[!, "SO4-2"]
    # index = findall(!ismissing, SO4)
    # SO4 = SO4[index]
    # t = df[!, "t (unadjusted) [days]"]
    # t = t[index]
    # # plot the data
    # lines!(ax7, t, SO4, label = "Column $i", color = colors[i])
    # scatter!(ax7, t, SO4, label = "Column $i", color = colors[i])
    # hlines!(ax7, c_so4_in, color = colors[i], linestyle = :dash, label = "Column $i SO4-2 inflow")
    # lines!(ax2l, t, SO4, label = "Column $i", color = colors[i])
    # scatter!(ax2l, t, SO4, label = "Column $i", color = colors[i])
    # hlines!(ax2l, c_so4_in, color = colors[i], linestyle = :dash, label = "Column $i SO4-2 inflow")
    # lines!(ax4e, t[t.< 8], SO4[t.< 8], label = "Column $i", color = colors[i])
    # scatter!(ax4e, t[t.< 8], SO4[t.< 8], label = "Column $i", color = colors[i])
    # hlines!(ax4e, c_so4_in, color = colors[i], linestyle = :dash, label = "Column $i SO4-2 inflow")
    # # plot in the early time data
    # get the data for the column
    Q = df[!, "flow_rate"]
    index = findall(!ismissing, Q)
    Q = Q[index]
    t = df[!, "t (unadjusted) [days]"]
    t = t[index]
    # plot the data
    lines!(ax4, t, Q, label = "Column $i", color = colors[i])
    scatter!(ax4, t, Q, label = "Column $i", color = colors[i])
    # Fe data
    # Fe = df[!, "Fe2+"]
    # index = findall(!ismissing, Fe)
    # Fe = Fe[index]
    # t = df[!, "t (unadjusted) [days]"]
    # t = t[index]
    # # plot the data
    # lines!(ax5, t, Fe, label = "Column $i", color = colors[i])
    # scatter!(ax5, t, Fe, label = "Column $i", color = colors[i])
    # # plot in the early time data
    # lines!(ax3e, t[t.< 8], Fe[t.< 8], label = "Column $i", color = colors[i])
    # scatter!(ax3e, t[t.< 8], Fe[t.< 8], label = "Column $i", color = colors[i])
end
# for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax1l, ax2l]
#     # add the vertical lines for the events
#     vlines!(ax, t_1_d, color = :black, label = "1mM NO3- switch")
#     vlines!(ax, t_2_d, color = :darkblue, label = "2mM NO3- switch")
#     vlines!(ax, t_s_34_d, color = :crimson, label = "Switch lactate off in columns 3 and 4")
#     vlines!(ax, t_stop_d, color = :brown, label = "Pumps stopped due to clogging")
#     vlines!(ax, t_restart_d, color = :lightgreen, label = "Pumps restarted after clogging")
#     vlines!(ax, t_4_d, color = :purple, label = "4mM NO3- switch")
# end
# for ax in [ax1e, ax2e, ax3e, ax4e]
#     vlines!(ax, t_s_34_d, color = :crimson, label = "Switch lactate off in columns 3 and 4")
# end
ylims!(ax7, 2000, 2700)
ylims!(ax2l, 2000, 2700)
ylims!(ax4e, 2400, 2700)
Legend(fig[3, 2:3], ax2, merge = true,
framevisible = false, nbanks = 2)
# adjust the layout
CairoMakie.resize_to_layout!(fig)
fig
save("plots/ssexp_overview.png", fig)