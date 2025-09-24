# You can type the `±` symbol in several ways:

# 1. **In Julia REPL or editors with Julia support**: Type `\pm` and press Tab
# 2. **On Mac**: Press `Option + Shift + =`
# 3. **On Windows**: Hold `Alt` and type `0177` on the numeric keypad
# 4. **Copy and paste**: Copy `±` from anywhere and paste it

# The `\pm<Tab>` method is most common in Julia environments since many editors support LaTeX-style symbol completion. 
using Measurements
using CairoMakie
using Statistics
using CSV
using DataFrames

df50 = CSV.read("data/falcon50mlweights.csv", DataFrame;header=false)
df15 = CSV.read("data/falcon15mlweights.csv", DataFrame;header=false)

# Make a histogram plot of the weights
fig_hist = Figure()
ax_hist50 = Axis(fig_hist[1, 1], xlabel="Weight (g)",
 ylabel="Frequency", title="50 mL Falcon Tube Weights")
ax_hist15 = Axis(fig_hist[2, 1], xlabel="Weight (g)",
 ylabel="Frequency", title="15 mL Falcon Tube Weights")
hist!(ax_hist50, df50[:, 1], bins=12, color=:blue, label="50 mL Tubes")
hist!(ax_hist15, df15[:, 1], bins=12, color=:orange, label="15 mL Tubes")
axislegend(ax_hist50)
axislegend(ax_hist15)
fig_hist
mean_vial_weight50 = mean(df50[:, 1])
vial_error50 = std(df50[:, 1])
vial_weight50 = mean_vial_weight50 ± vial_error50
mean_vial_weight15 = mean(df15[:, 1])
vial_error15 = std(df15[:, 1])
vial_weight15 = mean_vial_weight15 ± vial_error15
water_weights50 = collect(0.1:0.1:20)
water_weights15 = collect(0.1:0.1:10)
sampling_times = collect(5:1:240)

volumes50 = (water_weights50 .+ mean_vial_weight50) .- vial_weight50
volumes15 = (water_weights15 .+ mean_vial_weight15) .- vial_weight15
# error per volume collected:
fig = Figure()
Label(fig[1, 1:2], "Error analysis of volumes and flow rate by using the average vial weight", fontsize=16)
ax = Axis(fig[2, 1], xlabel="Collected Volume (mL)", ylabel="Error in Volume %", title="Falcon Tube 50ml")
lines!(ax, Measurements.value.(volumes50), Measurements.uncertainty.(volumes50) ./ Measurements.value.(volumes50) .* 100,
 color=:blue, label="Error %")
hlines!(ax, [5], color=:red, linestyle=:dash, label="5% Error")
hlines!(ax, [1], color=:orange, linestyle=:dash, label="1% Error")
axislegend(ax)
text!(ax, 3, 20; text="mean weight: $(round(Measurements.value(vial_weight50); digits=2)) ± $(round(Measurements.uncertainty(vial_weight50); digits=2))", color=:black)
ax2 = Axis(fig[2, 2], xlabel="Collected Volume (mL)",
 ylabel="Error in Volume %", title="Falcon Tube 15ml",
 xticks = 0:2:10)
lines!(ax2, Measurements.value.(volumes15), Measurements.uncertainty.(volumes15) ./ Measurements.value.(volumes15) .* 100,
 color=:blue, label="Error %")
hlines!(ax2, [5], color=:red, linestyle=:dash, label="5% Error")
hlines!(ax2, [1], color=:orange, linestyle=:dash, label="1% Error")
axislegend(ax2)
text!(ax2, 2, 7; text="mean weight: $(round(Measurements.value(vial_weight15); digits=2)) ± $(round(Measurements.uncertainty(vial_weight15); digits=2))", color=:black)
fig

# Create meshgrid for volumes and sampling_times
volume_grid = repeat(volumes', length(sampling_times), 1) # shape: (n_times, n_volumes)
time_grid = repeat(sampling_times, 1, length(water_weights)) # shape: (n_times, n_volumes)

flowrates_grid = volume_grid ./ time_grid
error_percentage_grid = Measurements.uncertainty.(flowrates_grid) ./ Measurements.value.(flowrates_grid) .* 100
flowrates = vec(flowrates_grid)
error_percentage = vec(error_percentage_grid)

fig_grid = Figure()
ax_grid = Axis(fig_grid[1, 1], xlabel="Sampling Time (s)", ylabel="Collected Volume (mL)", title="Error Percentage (%)")
hm = heatmap!(ax_grid, Measurements.value.(volumes), sampling_times, error_percentage_grid')
Colorbar(fig_grid[1, 2], hm, label="Error %")
fig_grid

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Flow rate (mL/min)", ylabel="Error Percentage (%)", title="Error Percentage vs Flowrate")
lines!(ax, flowrates, error_percentage, color=:blue, label="Error %")
axislegend(ax)
fig

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Sampling Time (s)", ylabel="Flowrate (mL/min)", title="Flowrate vs Sampling Time with Error Bars")
scatter!(ax, sampling_times, flowrates, color=:blue, label="Flowrate")
errorbars!(ax, sampling_times, Measurements.value.(flowrates), [Measurements.uncertainty(f) for f in flowrates], color=:blue)
axislegend(ax)
fig