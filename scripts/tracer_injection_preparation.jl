using SpecialFunctions
using DrWatson
using CairoMakie
using DataFrames, XLSX, Statistics
using Dates


# load analytical data
file_path = datadir("exp_raw", "no2_analyticall_results.xlsx")
sheet_name = "br_standard_curve"
data = XLSX.readtable(file_path, sheet_name)
df = DataFrame(data)
# start-time and end-time columns are of type Any because some of the values are in DateTiem and some in Time.
# Assume the Time happened on 17.06.2025 and convert them to DateTime
Hour(df[5,"start-time"])
start_time = DateTime[]
end_time = DateTime[]
for i in eachindex(df[!,"start-time"])
    if typeof(df[i, "start-time"]) == DateTime
        push!(start_time, df[i, "start-time"])
    else #type = Time
        push!(start_time, DateTime(2025, 06, 17, Dates.value(Hour(df[i, "start-time"])),
            Dates.value(Minute(df[i, "start-time"]))))
    end
    if typeof(df[i, "end-time"]) == DateTime
        push!(end_time, df[i, "end-time"])
    else #type = Time
        push!(end_time, DateTime(2025, 06, 17, Dates.value(Hour(df[i, "end-time"])),
            Dates.value(Minute(df[i, "end-time"]))))
    end
end
t0_dt = DateTime(2025, 06, 16, 22, 48) # start time of the experiment
# convert start_time and end_time to seconds since t0_dt
start_time_s = [Dates.Second(t - t0_dt) for t in start_time]
end_time_s = [Dates.Second(t - t0_dt) for t in end_time]
Q = convert.(Float64, df[!,"Q"]) # flow rate in μL/min
Q = Q .* 1e-3 / 60 # convert to cm3/s


#Dead volumes per column
dvs = Dict(1 => 35,
           2 => 29,
           3 => 34,
           4 => 30) # in cm
dvs_t0 = Dict(1 => 66,
              2 => 40+57,
              3 => 41+38,
              4 => 41+38.5) # in cm
tube_diam = 0.152 # cm
# Dead volume in cm3
dv = Dict(1 => dvs[1] * π * (tube_diam/2)^2,
          2 => dvs[2] * π * (tube_diam/2)^2,
          3 => dvs[3] * π * (tube_diam/2)^2,
          4 => dvs[4] * π * (tube_diam/2)^2)
dv_t0 = Dict(1 => dvs_t0[1] * π * (tube_diam/2)^2,
             2 => dvs_t0[2] * π * (tube_diam/2)^2,
             3 => dvs_t0[3] * π * (tube_diam/2)^2,
             4 => dvs_t0[4] * π * (tube_diam/2)^2)
t0s = Dict(1 => t0_dt + Dates.Second(floor(Int64, dv_t0[1] / Q[1])),
            2 => t0_dt + Dates.Second(floor(Int64, dv_t0[2] / Q[2])),
            3 => t0_dt + Dates.Second(floor(Int64, dv_t0[3] / Q[3])),
            4 => t0_dt + Dates.Second(floor(Int64, dv_t0[4] / Q[4])))
# calculate the start and end
column = convert.(Int64,df[!,"column"])
start_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in start_time[column .== 1]],
                        2 => [Dates.value(Dates.Second(t - t0s[2])) for t in start_time[column .== 2]],
                        3 => [Dates.value(Dates.Second(t - t0s[3])) for t in start_time[column .== 3]],
                        4 => [Dates.value(Dates.Second(t - t0s[4])) for t in start_time[column .== 4]])
end_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in end_time[column .== 1]],
                      2 => [Dates.value(Dates.Second(t - t0s[2])) for t in end_time[column .== 2]],
                      3 => [Dates.value(Dates.Second(t - t0s[3])) for t in end_time[column .== 3]],
                      4 => [Dates.value(Dates.Second(t - t0s[4])) for t in end_time[column .== 4]])
# save the flow rate data to a file
Qs = Dict(1 => Q[column .== 1],
          2 => Q[column .== 2],
          3 => Q[column .== 3],
          4 => Q[column .== 4])
avg_times_dict = Dict(1 => (start_times_dict[1] .+ end_times_dict[1]) ./ 2 .- dv[1] ./ Qs[1],
    2 => (start_times_dict[2] .+ end_times_dict[2]) ./ 2 .- dv[2] ./ Qs[2],
    3 => (start_times_dict[3] .+ end_times_dict[3]) ./ 2 .- dv[3] ./ Qs[3],
    4 => (start_times_dict[4] .+ end_times_dict[4]) ./ 2 .- dv[4] ./ Qs[4]
)

flow_rate_data = Dict("Qs" => Qs, "start_times_dict" => start_times_dict, "end_times_dict" => end_times_dict)
save("data/br_flow_rate.jld2", flow_rate_data)

Br = df[!,"Br-"] # concentration in mM
## calculate the dead times to correct the start and end times.

Br_dict = Dict(1 => Br[column .== 1],
          2 => Br[column .== 2],
          3 => Br[column .== 3],
          4 => Br[column .== 4])
# Now fix the is missing values in Br (check the missing value in each column and fix also the avg time dict)
Br_ = Dict{Int64, Vector{Float64}}(1 => Br_dict[1][.!ismissing.(Br_dict[1])],
          2 => Br_dict[2][.!ismissing.(Br_dict[2])],
          3 => Br_dict[3][.!ismissing.(Br_dict[3])],
          4 => Br_dict[4][.!ismissing.(Br_dict[4])])
avg_time = Dict{Int64, Vector{Float64}}(1 => avg_times_dict[1][.!ismissing.(Br_dict[1])],
          2 => avg_times_dict[2][.!ismissing.(Br_dict[2])],
          3 => avg_times_dict[3][.!ismissing.(Br_dict[3])],
          4 => avg_times_dict[4][.!ismissing.(Br_dict[4])])
# Now we have the Br and avg_time dictionaries with the missing values removed.
# Ok now we are ready to create the datasets per column

# Now we can plot the data
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time [s]",
    ylabel = "Br⁻ [mM]",
    title = "Bromide breakthrough curves")
for col in 1:4
    # get the data for the column
    Br_col = Br_[col]
    avg_time_col = avg_time[col]
    # plot the data
    #lines!(ax, df.time, df.Br, label = "Column $col")
    scatter!(ax, avg_time_col, Br_col, label = "Column $col")
end
# add a legend
axislegend(ax, position = :lt, framevisible = false)
fig
# save the figure
save("plots/br_breakthrough_data.png", fig)


# save the datasets to a file

save("data/br_breakthrough_data.jld2", Dict("Br" => Br_, "avg_time" => avg_time))

