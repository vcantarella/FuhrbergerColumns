using CairoMakie
using DrWatson
using DataFrames
using CSV
using XLSX

# load analytical data
file_path = datadir("exp_raw", "no2_analyticall_results.xlsx")
sheet_name = "sampling"
data = XLSX.readtable(file_path, sheet_name)
df = DataFrame(data)


# make a df per column
df1 = df[df.column .== 1, :]
df2 = df[df.column .== 2, :]
df3 = df[df.column .== 3, :]
df4 = df[df.column .== 4, :]
dfs = [df1, df2, df3, df4]
# Make a plot for the most importntant data, skip misisng data
colors = [:red, :blue, :green, :orange]
fig = Figure()
ax1 = Axis(fig[1, 1], title = "pH", ylabel = "pH", xlabel = "t (days)")
ax2 = Axis(fig[2, 1], title = "EC", ylabel = "EC [μS/cm]", xlabel = "t (days)")
ax3 = Axis(fig[1, 2], title = "NO2-", ylabel = "NO2- [μM]", xlabel = "t (days)")
ax4 = Axis(fig[2, 2], title = "Q [μL/min]", ylabel = "Q [μL/min]", xlabel = "t (days)")
ax5 = Axis(fig[2, 3], title = "Fe2+", ylabel = "Fe2+ [μM]", xlabel = "t (days)")
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
end
Legend(fig[3, 1:3], ax2, position = (0.5, 0.5), orientation = :horizontal, merge = true)
fig
save("lab_data.png", fig)