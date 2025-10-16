using DrWatson
using DataFrames, XLSX, Statistics
using Dates
using StaticArrays
using CairoMakie
import DataInterpolations as DI
import RegularizationTools as RT
include("model_data_structures.jl")
# load analytical data
file_path = datadir("exp_raw", "ssexp_data.xlsx")
sheet_names = ["general_samples_p$i" for i in 1:4]
datas = [XLSX.readtable(file_path, sheet_name) for sheet_name in sheet_names]
dfs = [DataFrame(data) for data in datas]
# start-time and end-time columns are of type Any because some of the values are in DateTiem and some in Time.
df = vcat(dfs...)
start_time = df[!,"start_time"]
end_time = df[!,"end_time"]
tracer_sheet = "bromide_curve_v2"
df_tr = DataFrame(XLSX.readtable(file_path, tracer_sheet))
re_tr = r"^B(\d?)"
df_tr[!, "column"] = match.(re_tr, df_tr[!,"Sample"]) .|> x -> x.captures[1] |> x -> parse(Int, x)
start_time = vcat(start_time, df_tr[!,"start_time"])
end_time = vcat(end_time, df_tr[!,"end_time"])
flow_rate = vcat(df[!,"flow_rate"], df_tr[!,"flow_rate"])
# Start time for the t0 of the experiment for columns 1 and 2
t0 = DateTime(2025, 09, 18, 18, 15)
re = r"^P(\d?)"
df[!, "column"] = match.(re, df[!,"Sample"]) .|> x -> x.captures[1] |> x -> parse(Int, x)
column = vcat(df[!,"column"], df_tr[!,"column"])
# convert start_time and end_time to seconds since t0_dt
Q0 = Dict(
    1 =>  flow_rate[(.!ismissing.(flow_rate) .& (column .== 1))][1]./ 3600 .* 1e-6,
    2 => flow_rate[(.!ismissing.(flow_rate) .& (column .== 2))][1]./ 3600 .* 1e-6,
    3 => flow_rate[(.!ismissing.(flow_rate) .& (column .== 3))][1]./ 3600 .* 1e-6,
    4 => flow_rate[(.!ismissing.(flow_rate) .& (column .== 4))][1]./ 3600 .* 1e-6
)

#Dead volumes per column
# dvs: dead volumes from the column to the sampler (outlet of the column)
dvs = Dict(1 => 29,
           2 => 27,
           3 => 33.5,
           4 => 33.5,
           ) # in cm
# dvs_t0: dead volumes from the tedlar bag to the column (outlet of the column)
base_dv = 30 + 7.5
dvs_t0 = Dict(1 => base_dv + 38 + 48,
              2 => base_dv + 38 + 25,
              3 => base_dv + 38 + 23,
              4 => 13 + 38 + 22,
              ) # in cm
tube_diam = 0.152 # cm
# Dead volume in cm3
dv = Dict(1 => dvs[1] * π * (tube_diam/2)^2,
          2 => dvs[2] * π * (tube_diam/2)^2,
          3 => dvs[3] * π * (tube_diam/2)^2,
          4 => dvs[4] * π * (tube_diam/2)^2
          )
dv_t0    = Dict(1 => dvs_t0[1] * π * (tube_diam/2)^2,
             2 => dvs_t0[2] * π * (tube_diam/2)^2,
             3 => dvs_t0[3] * π * (tube_diam/2)^2,
             4 => dvs_t0[4] * π * (tube_diam/2)^2
             )
t0s = Dict(1 => t0 + Dates.Second(floor(Int64, dv_t0[1]*1e-6 / Q0[1])),
            2 => t0 + Dates.Second(floor(Int64, dv_t0[2]*1e-6 / Q0[2])),
            3 => t0 + Dates.Second(floor(Int64, dv_t0[3]*1e-6 / Q0[3])),
            4 => t0 + Dates.Second(floor(Int64, dv_t0[4]*1e-6 / Q0[4]))
            )
# Now we have the t0 reference for the discharge data.

disch_ds = Dict()
for i in 1:4
    # only the values where there are no missing values in the flow rate
    bool_index = (.!ismissing.(flow_rate) .& (column .== i))
    # convert the flow rate to μL/min
    Q = flow_rate[bool_index] ./ 3600 .* 1e-6 # convert from ml/hr to m3/s
    end_times = Dates.Second.(end_time[bool_index] .- t0s[i]) # end times in seconds
    # create a QData object and push it to the dictionary
    disch_ds[i] = QData(Q, Dates.value.(end_times))
end

"""
    disch_function(t, Qs, end_times)

Calculate the discharge/flow rate at a given time `t` based on piecewise constant flow rates.

# Arguments
- `t`: Time value for which to determine the flow rate
- `Qs`: Vector of flow rate values corresponding to different time periods
- `end_times`: Vector of end times for each flow rate period (must be sorted in ascending order)

# Returns
- Flow rate value at time `t`. Returns the flow rate for the first time period where `t ≤ end_times[i]`, 
  or the last flow rate value if `t` exceeds all end times.
"""
function disch_function(t, Qs, end_times)
    for i in eachindex(Qs)
        if t <= end_times[i]
            return Qs[i]   # This should be inside the if block
        end
    end
    return Qs[end]
end


# import transport parameters
transp_params = load("data/optimized_tracer_params_m2.jld2")
tracer_params = transp_params["tracer_params"]
# Generate the transport dataset for each column
# Q is the flow rate in μL/min for each column
# And we calculate the v and αₗ*v that are dependent on the current Q


v_ds = Dict{Int, VData}()
v_st = Dict{Int, VDataS}()
v_da = Dict{Int, VDataA}()
## Make a data Interpolation of the velocity data
v_interp = Dict{Int, DI.RegularizationSmooth}()
D = 3.5*1e-2 #cm to m diameter of the column
A = π * D^2 / 4 # Cross-sectional area
for i in 1:4
    # only the values where there are no missing values in the flow rate
    if i < 3
        ϕ = tracer_params[i][1]
        αₗ = tracer_params[i][2]
    else
        ϕ = mean([tracer_params[i][1] for i in 1:3])
        αₗ = mean([tracer_params[i][2] for i in 1:3])
    end
    bool_index = (.!ismissing.(flow_rate) .& (column .== i))
    # convert the flow rate to μL/min
    Q = flow_rate[bool_index] ./ 3600 .* 1e-6 # convert from ml/hr to m3/s
    end_times = Dates.value.(Dates.Second.(end_time[bool_index] .- t0s[i])) # end times in seconds
    start_times = Dates.value.(Dates.Second.(start_time[bool_index] .- t0s[i])) # start times in seconds
    avg_times = start_times .+ (end_times .- start_times) ./ 2
    # create a VData object and push it to the dictionary
    v = Q./(ϕ * A) # flow velocity in m/s
    v = convert.(Float64, v) # convert to Float64
    d = 2
    λ = 1e9
    dense_sample_times = range(minimum(avg_times), stop=maximum(avg_times), length=500)
    Am = DI.RegularizationSmooth(v, avg_times, d; λ=λ, alg = :fixed, extrapolation=ExtrapolationType.Constant)
    
    u = Am(dense_sample_times)
    interp = Am

    #De = v .* αₗ # longitudinal dispersion coefficient in m2/s
    # create a QData object and push it to the vector
    v_ds[i] = VData(v, end_times)
    v_st[i] = VDataS(v, start_times)
    v_da[i] = VDataA(v, start_times, end_times)
    v_interp[i] = interp
end


# Making a Figure for plotting the flow velocity for each column
fig_v = Figure()
ax = Axis(fig_v[1,1], title="Flow velocity in the model",
xlabel ="Time (d)", ylabel="Flow velocity (m/s)")
flowt = 0:0.0001:27
colors = [:blue, :orange, :green, :red]
for i in 1:4
    vds = v_da[i]
    # Calculate the flow velocity for the given time points
    flowv = [v_interp[i](t*24*60*60) for t in flowt]
    # Create a line plot for the flow velocity
    lines!(ax, flowt, flowv, label="Column $i", color=colors[i],
    linewidth=2)
    # Check the data points
    avg_times = vds.start_times .+ (vds.end_times .- vds.start_times) ./ 2
    scatter!(ax, avg_times ./ 86400, vds.v, color   = colors[i], label="Data $i", markersize=10)
end
axislegend(ax, position = :lt, framevisible = false)
resize_to_layout!(fig_v)
fig_v
save("plots/flow_velocity_model_m2.png", fig_v)

## Events:
# Moment when we switch to 1.5 mM inflow concentration
t_switch_1_5mM = DateTime(2025, 10, 06, 17, 10)
c_no3_1_5mM = 1.5
t_switch_1mM = DateTime(2025, 10, 12, 20, 35)
c_no3_1mM = 1.0

c_ins = Dict{Int64, CinData}()
for i in 1:4
    c_no3 = 2e-3
    c_doc = 0.0
    c_so4 = 0e-3 # concentration of SO4-2 in the input solution [mM]
    c_fe = 0.0e-3 # concentration of Fe+2 in the input solution [mM]

    cins = [[c_no3, 1e-16, c_so4, c_fe, c_doc, c_no3],]
    t0switch = []
    t_1 = Dates.value(Dates.Second(t_switch_1_5mM - t0s[i])) # convert days to seconds
    #t_1 += dv_t0[i]*1e-6 / disch_function(t_1, disch_ds[i].Q, disch_ds[i].end_times) # convert days to seconds
    push!(cins, [c_no3_1_5mM*1e-3, 1e-16, c_so4, c_fe, c_doc, c_no3_1_5mM*1e-3])
    push!(t0switch, t_1) # convert days to seconds
    cins = vcat(cins, [[c_no3_1mM*1e-3, 1e-16, c_so4, c_fe, c_doc, c_no3_1mM*1e-3],])
    t_2 = Dates.value(Dates.Second(t_switch_1mM - t0s[i])) # convert days to seconds
    t_2 += dv_t0[i]*1e-6 / disch_function(t_2, disch_ds[i].Q, disch_ds[i].end_times) # convert days to seconds
    push!(t0switch, t_2) # convert days to seconds
    t0switch = convert.(Float64, t0switch) # convert to seconds
    # Add the initial concentration at t0
    c_ins[i] = CinData(cins, t0switch)
end
