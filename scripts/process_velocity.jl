include("../src/lowpass_filter.jl")
include("model_data_structures.jl")
using JLD2
using CairoMakie
@load "data/flow_velocity_data_model_m2.jld2" v_da

ads = v_da[4]

v = ads.v
t = ads.start_times .+ (ads.end_times .- ads.start_times) ./ 2

α = 0.8
v_filt = lowpass_filter(v, α)

A = 1.0
P = 0.0
Q = 2e-7^2
R = 3.e-7^2
x = v[1]
v_k = Float64[]
v_k_ = Float64[]
Pk = Float64[]
Pk_ = Float64[]
for n in 1:length(v)
    x = A * x
    P = A * P * A + Q
    K = P * 1' / (1 * P * 1' + R)
    x = x + K * (v[n] - 1 * x)
    P = (1 - K * 1) * P
    push!(v_k, x)
end
v_k

# smooth


# Making a Figure for plotting the flow velocity for each column
fig_v = Figure()
ax = Axis(fig_v[1,1], title="Flow velocity in the model",
xlabel ="Time (d)", ylabel="Flow velocity (m/s)")
flowt = 0:0.0001:27
# Create a line plot for the flow velocity
lines!(ax, t ./ 86400, v_k, label="Column 4", color=:red,
linewidth=2)
# Check the data points
scatter!(ax, t ./ 86400, ads.v, color   = :red, label="Column 4", markersize=10)
axislegend(ax, position = :lt, framevisible = false)
resize_to_layout!(fig_v)
fig_v