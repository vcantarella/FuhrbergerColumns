using DrWatson
@quickactivate "FuhrbergerColumns"
using CairoMakie
using DataFrames
using CSV
using Statistics
# Here you may include files from the source directory
# include(srcdir("dummy_src_file.jl"))

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Have fun with your new project!

You can help us improve DrWatson by opening
issues on GitHub, submitting feature requests,
or even opening your own Pull Requests!
"""
)
# Load the data
dfavs = CSV.File(datadir("exp_pro","processed_AVS_data.csv")) |> DataFrame


# AVS plots
# Create a figure
fig = Figure(size = (800, 600))
# pH Axis
ax_avs = Axis(fig[1, 1], xlabel = "column length (cm)", ylabel = "AVS(S-2) concentration (nmol/g)",
    xlabelsize = 14, ylabelsize = 14, xticklabelsize = 12, yticklabelsize = 12)
# plot the data from each sheet
## create a color palette
palette = [:blue, :green, :red]
## iterate over the sheets

column = dfavs[!,"column"]
a = column .== "A"
b = column .== "B"
c = column .== "C"
x = dfavs[a,"x"]
avs_a = collect(skipmissing(dfavs[a,"S-2"]))
avs_b = collect(skipmissing(dfavs[b,"S-2"]))
avs_c = collect(skipmissing(dfavs[c,"S-2"]))
lines!(ax_avs, x, avs_a, label = "Column A", color = palette[1], linestyle = :dash)
lines!(ax_avs, x, avs_b, label = "Column B", color = palette[2], linestyle = :dash)
lines!(ax_avs, x, avs_c, label = "Column C", color = palette[3], linestyle = :dash)
scatter!(ax_avs, x, avs_a, color = palette[1], markersize = 11)
scatter!(ax_avs, x, avs_b, color = palette[2], markersize = 11)
scatter!(ax_avs, x, avs_c, color = palette[3], markersize = 11)

# Averag AVS concentration
avs_avg_a = (mean(skipmissing(dfavs[a,"S-2 r1"])) + mean(skipmissing(dfavs[a,"S-2 r2"])))/2
avs_avg_b = mean(skipmissing(dfavs[b,"S-2 r1"]))
avs_avg_c = mean(skipmissing(dfavs[c,"S-2 r1"]))
lines!(ax_avs, [0, 12], [avs_avg_a, avs_avg_a], label = "Column A average", color = palette[1], linestyle = :solid)
lines!(ax_avs, [0, 12], [avs_avg_b, avs_avg_b], label = "Column B average", color = palette[2], linestyle = :solid)
lines!(ax_avs, [0, 12], [avs_avg_c, avs_avg_c], label = "Column C average", color = palette[3], linestyle = :solid)


# Add legend
Legend(fig[2, 1], ax_avs, "columns", orientation = :horizontal)
# save
save(plotsdir("AVS_plot.png"), fig)