using DrWatson
using DataFrames, XLSX, Statistics
using Dates
using StaticArrays
include("model_data_structures.jl")
# load analytical data
file_path = datadir("exp_raw", "no2_analyticall_results.xlsx")
sheet_name = "sampling"
data = XLSX.readtable(file_path, sheet_name)
df = DataFrame(data)
# start-time and end-time columns are of type Any because some of the values are in DateTiem and some in Time.
# Assume the Time happened on 17.06.2025 and convert them to DateTime
Hour(df[5,"start time"])
start_time = df[!,"start time"]
end_time = df[!,"end time"]
# Start time for the t0 of the experiment for columns 1 and 2
t0_1_2 = DateTime(2025, 05, 12, 20, 0, 0)
# Start time for the t0 of the experiment for columns 3 and 4
t0_3_4 = DateTime(2025, 05, 12, 20, 14, 0)
# convert start_time and end_time to seconds since t0_dt
Q0 = Dict(
    1 =>  df[(.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== 1),"Q"][1]./ 60 .* 1e-9,
    2 => df[(.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== 2),"Q"][1]./ 60 .* 1e-9,
    3 => df[(.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== 3),"Q"][1]./ 60 .* 1e-9,
    4 => df[(.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== 4),"Q"][1]./ 60 .* 1e-9
)
dvs_t0 = Dict(1 => 66,
              2 => 40+57,
              3 => 41+38,
              4 => 41+38.5) # in cm
tube_diam = 1.52e-3 # m
dead_t0 = Dict(1 => dvs_t0[1]*0.01 * π * (tube_diam/2)^2,
             2 => dvs_t0[2]*0.01 * π * (tube_diam/2)^2,
             3 => dvs_t0[3]*0.01 * π * (tube_diam/2)^2,
             4 => dvs_t0[4]*0.01 * π * (tube_diam/2)^2) #dead volume in m3
# First we have to correct the t0, because it is actually after the dvs_t0.
t0sdict = Dict(1 => t0_1_2 + Dates.Second(round(dvs_t0[1] / 100 * tube_diam^2 * π / 4 / Q0[1])),
            2 => t0_1_2 + Dates.Second(round(dvs_t0[2] / 100 * tube_diam^2 * π / 4 / Q0[2])),
            3 => t0_3_4 + Dates.Second(round(dvs_t0[3] / 100 * tube_diam^2 * π / 4 / Q0[3])),
            4 => t0_3_4 + Dates.Second(round(dvs_t0[4] / 100 * tube_diam^2 * π / 4 / Q0[4])))


# Now we have the t0 reference for the discharge data.

disch_ds = Dict()
for i in 1:4
    # only the values where there are no missing values in the flow rate
    bool_index = (.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== i)
    # convert the flow rate to μL/min
    Q = df[bool_index,"Q"] ./ 60 .* 1e-9 # convert from μL/min to m3/s
    end_times = Dates.Second.(end_time[bool_index] .- t0sdict[i]) # end times in seconds
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

## Calculating the velocity dataset for each column
## Here I include not only the flow rate, but I also consider when the pump has stopped (due to problems)
# Moment when pumps 1 and 2 stopped due to clogging.
t_stop_1 = DateTime(2025,05,28, 19, 50)
t_stop_2 = DateTime(2025,05,28, 20,10)
# moment when the pumps were restarted
t_restart_1_2 = DateTime(2025,05,29, 08, 45)
# Approximate time when all pumps stoppeds due to energie failure in the building
t_stop_all = DateTime(2025,05,29, 12, 30)
# Moment when the pumps were restarted after the energy failure
t_restart_all = DateTime(2025,05,29, 18, 16)

# import transport parameters
transp_params = load("data/optimized_tracer_params.jld2")
tracer_params = transp_params["tracer_params"]
# Generate the transport dataset for each column
# Q is the flow rate in μL/min for each column
# And we calculate the v and αₗ*v that are dependent on the current Q


v_ds = Dict{Int, VData}()
D = 3.5*1e-2 #cm to m diameter of the column
A = π * D^2 / 4 # Cross-sectional area
for i in 1:4
    # only the values where there are no missing values in the flow rate
    ϕ = tracer_params[i][1]
    αₗ = tracer_params[i][2]
    bool_index = (.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== i)
    # convert the flow rate to μL/min
    Q = df[bool_index,"Q"] ./ 60 .* 1e-9 # convert from μL/min to m3/s
    end_times = Dates.value.(Dates.Second.(end_time[bool_index] .- t0sdict[i])) # end times in seconds
    # create a VData object and push it to the dictionary
    # check when the flow rate was 0
    t_stop_1_s = Dates.value(Dates.Second(t_stop_1 - t0sdict[i]))
    t_stop_2_s = Dates.value(Dates.Second(t_stop_2 - t0sdict[i]))
    t_restart_1_2_s = Dates.value(Dates.Second(t_restart_1_2 - t0sdict[i]))
    t_stop_s = Dates.value(Dates.Second(t_stop_all - t0sdict[i]))
    t_restart_s = Dates.value(Dates.Second(t_restart_all - t0sdict[i]))

    if i == 1
        ind_i = findlast(end_times .< t_stop_1_s)
        Q = vcat(Q[1:ind_i],[Q[ind_i],0.],Q[ind_i+1:end]) # take only the values before the stop
        end_times = vcat(end_times[1:ind_i], [t_stop_1_s, t_restart_1_2_s], end_times[ind_i+1:end]) # add the stop time
    elseif i == 2
        ind_i = findlast(end_times .< t_stop_2_s)
        Q = vcat(Q[1:ind_i],[Q[ind_i],0.],Q[ind_i+1:end]) # take only the values before the stop
        end_times = vcat(end_times[1:ind_i], [t_stop_2_s, t_restart_1_2_s], end_times[ind_i+1:end]) # add the stop time
    end
    
    ind = findlast(end_times .< t_stop_s)
    Q = vcat(Q[1:ind],[Q[ind], 0.],Q[ind+1:end]) # take only the values before the stop
    end_times = vcat(end_times[1:ind], [t_stop_s, t_restart_s], end_times[ind+1:end]) # add the stop time
    v = Q./(ϕ * A) # flow velocity in m/s
    v = convert.(Float64, v) # convert to Float64
    #De = v .* αₗ # longitudinal dispersion coefficient in m2/s
    # create a QData object and push it to the vector
    v_ds[i] = VData(v, end_times)
end



## Events:
# Start time for the t0 of the experiment for columns 1 and 2
c_no3_1_2 = 0.6 # concentration of NO3- in the input solution for columns 1 and 2 [mM]
c_lac_1_2 = 0.0 # concentration of lactate in the input solution for columns 1 and 2 [mM]
# Start time for the t0 of the experiment for columns 3 and 4
c_no3_3_4 = 0.6 # concentration of NO3- in the input solution for columns 3 and 4 [mM]
c_lac_3_4 = 0.33 # concentration of lactate in the input solution for columns 3 and 4 [mM]
# Moment when we swithed the bag for columns 3 and 4, removing the lactate solution
t_switch_3_4 = DateTime(2025, 05, 16, 16, 35, 0)
c_lac_3_4_switch = 0.0 # concentration of lactate in the input solution for columns 3 and 4 after the switch [mM]
# Moment when we switch the bag for all columns, increasing the NO3- concentration
t_switch_all_1 = DateTime(2025, 05, 23, 19, 00, 0)
c_no3_all_1 = 1.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
# Moment when the increase the NO3- concentration
t_switch_all_2 = DateTime(2025, 05, 28, 19, 15, 0)
c_no3_all_2 = 2.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
# Moment when we switch to 4 mM inflow concentration
t_switch_4mM = DateTime(2025, 06, 03, 18, 25, 0)
c_no3_4mM = 4.

c_ins = Dict{Int64, CinData}()
for i in 1:4
    c_no3 = 0.6e-3
    c_lac = 0.0
    c_so4 = 2.3e-3 # concentration of SO4-2 in the input solution [mM]
    c_fe = 0.0e-3 # concentration of Fe+2 in the input solution [mM]

    if i == 3 || i == 4
        c_lac = 0.33e-3 # concentration of NO3- in the input solution for columns 3 and 4 [mM]
    end
    cins = [[c_no3, 1e-16, c_so4, c_fe, c_lac, c_no3],]
    t0s = []
    if i == 3 || i == 4
        # For columns 3 and 4, we have a switch in the lactate concentration
        push!(cins, [c_no3, 1e-16, c_so4, c_fe, c_lac_3_4_switch, c_no3])
        t_lac = Dates.value(Dates.Second(t_switch_3_4 - t0sdict[i])) # convert days to seconds
        t_lac += dead_t0[i] / disch_function(t_lac, disch_ds[i].Q, disch_ds[i].end_times)
        push!(t0s, t_lac) # convert days to seconds
    end
    t_1 = Dates.value(Dates.Second(t_switch_all_1 - t0sdict[i])) # convert days to seconds
    t_1 += dead_t0[i] / disch_function(t_1, disch_ds[i].Q, disch_ds[i].end_times) # convert days to seconds
    push!(cins, [c_no3_all_1*1e-3, 1e-16, c_so4, c_fe, c_lac_3_4_switch, c_no3_all_1*1e-3])
    push!(t0s, t_1) # convert days to seconds
    push!(cins, [c_no3_all_2*1e-3, 1e-16, c_so4, c_fe, c_lac_3_4_switch, c_no3_all_2*1e-3])
    t_2 = Dates.value(Dates.Second(t_switch_all_2 - t0sdict[i])) # convert days to seconds
    t_2 += dead_t0[i] / disch_function(t_2, disch_ds[i].Q, disch_ds[i].end_times) # convert days to seconds
    push!(t0s, t_2) # convert days to seconds
    push!(cins, [c_no3_4mM*1e-3, 1e-16, c_so4, c_fe, c_lac_3_4_switch, c_no3_4mM*1e-3])
    t_4 = Dates.value(Dates.Second(t_switch_4mM - t0sdict[i])) # convert days to seconds
    t_4 += dead_t0[i] / disch_function(t_4, disch_ds[i].Q, disch_ds[i].end_times) # convert days to seconds
    push!(t0s, t_4) # convert days to seconds
    t0s = convert.(Float64, t0s) # convert to seconds
    # Add the initial concentration at t0
    c_ins[i] = CinData(cins, t0s)
end
