using CairoMakie
using Statistics
using Random

# Load empirical vial weights
vial_weights_15 = CSV.read("data/falcon15mlweights.csv", DataFrame;header=false)[:, 1]
vial_weights_50 = CSV.read("data/falcon50mlweights.csv", DataFrame;header=false)[:, 1]

# Choose which vial to use for analysis
vial_weights = vial_weights_15 # or vial_weights_50

water_weights = collect(1:6)
sampling_times = collect(30:30:240)
volumes = water_weights

# Create meshgrid for volumes and sampling_times
volume_grid = repeat(volumes', length(sampling_times), 1)
time_grid = repeat(sampling_times, 1, length(water_weights))

n_boot = 1000
error_grid = zeros(size(volume_grid))

for i in 1:size(volume_grid, 1)
    for j in 1:size(volume_grid, 2)
        boot_flowrates = Float64[]
        for b in 1:n_boot
            vial_sample = rand(vial_weights, 1)[1]
            vol = water_weights[j] + vial_sample - vial_sample # just water_weights[j]
            flowrate = vol / time_grid[i, j]
            push!(boot_flowrates, flowrate)
        end
        error_grid[i, j] = std(boot_flowrates) / mean(boot_flowrates) * 100
    end
end

fig_grid = Figure()
ax_grid = Axis(fig_grid[1, 1], xlabel="Sampling Time (s)", ylabel="Volume (mL)", title="Bootstrap Error Percentage (%)")
hm = heatmap!(ax_grid, sampling_times, volumes, error_grid')
Colorbar(fig_grid[1, 2], hm, label="Error %")
fig_grid
