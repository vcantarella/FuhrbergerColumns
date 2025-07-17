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


Hour(df[5,"start time"])
start_time = df[!,"start time"]
end_time = df[!,"end time"]
# Start time for the t0 of the experiment for columns 1 and 2
t0_1_2 = DateTime(2025, 05, 12, 20, 0, 0)
c_no3_1_2 = 0.6 # concentration of NO3- in the input solution for columns 1 and 2 [mM]
c_lac_1_2 = 0.0 # concentration of lactate in the input solution for columns 1 and 2 [mM]
# Start time for the t0 of the experiment for columns 3 and 4
t0_3_4 = DateTime(2025, 05, 12, 20, 14, 0)
t0_arr = [ifelse(col == 1 || col == 2, t0_1_2, t0_3_4) for col in df[!,"column"]]
# convert start_time and end_time to seconds since t0_dt
start_time_s = Dates.Second.(start_time - t0_arr)
end_time_s = Dates.Second.(end_time - t0_arr)

# import transport parameters
transp_params = load("data/optimized_tracer_params.jld2")
tracer_params = transp_params["tracer_params"]
# Generate the transport dataset for each column
# Q is the flow rate in μL/min for each column
# And we calculate the v and αₗ*v that are dependent on the current Q
struct QData{QT, T}
    column::Int64
    Q::QT # flow velocity in m/s
    end_times::T
end
qvec = QData[]
D = 3.5*1e-2 #cm to m diameter of the column
A = π * D^2 / 4 # Cross-sectional area
for i in 1:4
    # only the values where there are no missing values in the flow rate
    bool_index = (.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== i)
    # convert the flow rate to μL/min
    Q = df[bool_index,"Q"] ./ 60 .* 1e-9 # convert from μL/min to m3/s
    # create a QData object and push it to the vector
    push!(qvec, QData(i, Q, Dates.value.(end_time_s[bool_index])))
end

# Moment when pumps 1 and 2 stopped due to clogging.
t_stop_1 = DateTime(2025,05,28, 19, 50)
t_stop_2 = DateTime(2025,05,28, 20,10)
# moment when the pumps were restarted
t_restart_1_2 = DateTime(2025,05,29, 08, 45)
# Approximate time when all pumps stoppeds due to energie failure in the building
t_stop_all = DateTime(2025,05,29, 12, 30)
t_stop_s = Dates.Second(t_stop_all - t0_1_2)
# Moment when the pumps were restarted after the energy failure
t_restart_all = DateTime(2025,05,29, 18, 16)
t_restart_s = Dates.Second(t_restart_all - t0_1_2)


function flow_function(t, Qs, end_times)
    for i in eachindex(Qs)
        if t <= end_times[i]
            return Qs[i]   # This should be inside the if block
        end
    end
    return Qs[end]
end


#Dead volumes per column
dvs = Dict(1 => 35,
           2 => 29,
           3 => 34,
           4 => 30) # in cm
dvs_t0 = Dict(1 => 66,
              2 => 40+57,
              3 => 41+38,
              4 => 41+38.5) # in cm
tube_diam = 1.52e-3 # m
# First we have to correct the t0, because it is actually after the dvs_t0.
t0s = Dict(1 => t0_1_2 + Dates.Second(round(dvs_t0[1] / 100 * tube_diam^2 * π / 4 / flow_function(0, qvec[1].Q, qvec[1].end_times))),
            2 => t0_1_2 + Dates.Second(round(dvs_t0[2] / 100 * tube_diam^2 * π / 4 / flow_function(0, qvec[2].Q, qvec[2].end_times))),
            3 => t0_3_4 + Dates.Second(round(dvs_t0[3] / 100 * tube_diam^2 * π / 4 / flow_function(0, qvec[3].Q, qvec[3].end_times))),
            4 => t0_3_4 + Dates.Second(round(dvs_t0[4] / 100 * tube_diam^2 * π / 4 / flow_function(0, qvec[4].Q, qvec[4].end_times))))

# Now, for each sample, we calculate the average time of the sample and remove the dead volume.
start_time
column = convert.(Int64,df[!,"column"])
start_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in start_time[column .== 1]],
                        2 => [Dates.value(Dates.Second(t - t0s[2])) for t in start_time[column .== 2]],
                        3 => [Dates.value(Dates.Second(t - t0s[3])) for t in start_time[column .== 3]],
                        4 => [Dates.value(Dates.Second(t - t0s[4])) for t in start_time[column .== 4]])
end_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in end_time[column .== 1]],
                      2 => [Dates.value(Dates.Second(t - t0s[2])) for t in end_time[column .== 2]],
                      3 => [Dates.value(Dates.Second(t - t0s[3])) for t in end_time[column .== 3]],
                      4 => [Dates.value(Dates.Second(t - t0s[4])) for t in end_time[column .== 4]])

avg_times_dict = Dict(1 => (start_times_dict[1] .+ end_times_dict[1]) ./ 2 .-
    dvs[1]/100*π*tube_diam^2/4 ./ flow_function.(end_times_dict[1], Ref(qvec[1].Q), Ref(qvec[1].end_times)),
    2 => (start_times_dict[2] .+ end_times_dict[2]) ./ 2 .-
    dvs[2]/100*π*tube_diam^2/4 ./ flow_function.(end_times_dict[2], Ref(qvec[2].Q), Ref(qvec[2].end_times)),
    3 => (start_times_dict[3] .+ end_times_dict[3]) ./ 2 .- 
    dvs[3]/100*π*tube_diam^2/4 ./ flow_function.(end_times_dict[3], Ref(qvec[3].Q), Ref(qvec[3].end_times)),
    4 => (start_times_dict[4] .+ end_times_dict[4]) ./ 2 .- 
    dvs[4]/100*π*tube_diam^2/4 ./ flow_function.(end_times_dict[4], Ref(qvec[4].Q), Ref(qvec[4].end_times)),
)

# Now we have the average times for each column, we can create the datasets per column
no3 = df[!,"NO3-"] # concentration in mM
no2 = df[!,"NO2-"] # concentration in mM
so4 = df[!,"SO4-2"] # concentration in mM
fe = df[!,"Fe2+"] # concentration in mM
pH = df[!,"pH"] # pH values
ec = df[!,"EC"] # electrical conductivity in μS/cm

struct conc_ds{C, T}
    conc::C # concentration in mM
    t::T # time in seconds since t0
    conc_ds(conc::C, t::T) where {C,T} = size(conc) == size(t) ? new{C,T}(conc, t) : error("Size mismatch between concentration and time")
end

struct ds{I, D}
    column::I
    no3::D
    no2::D
    so4::D
    fe::D
    pH::D
    ec::D
end


# For each column, we create a dataset with the concentrations and the average times
all_ds = Dict{String, ds}()
for i in 1:4
    # only the values where there are no missing values in the concentration
    id1 = column .== i
    no3_col = conc_ds(convert.(Float64, no3[id1][.!ismissing.(no3[id1])]), avg_times_dict[i][.!ismissing.(no3[id1])])
    no2_col = conc_ds(convert.(Float64, no2[id1][.!ismissing.(no2[id1])]), avg_times_dict[i][.!ismissing.(no2[id1])])
    so4_col = conc_ds(convert.(Float64, so4[id1][.!ismissing.(so4[id1])]), avg_times_dict[i][.!ismissing.(so4[id1])])
    fe_col = conc_ds(convert.(Float64, fe[id1][.!ismissing.(fe[id1])]), avg_times_dict[i][.!ismissing.(fe[id1])])
    pH_col = conc_ds(convert.(Float64, pH[id1][.!ismissing.(pH[id1])]), avg_times_dict[i][.!ismissing.(pH[id1])])
    ec_col = conc_ds(convert.(Float64, ec[id1][.!ismissing.(ec[id1])]), avg_times_dict[i][.!ismissing.(ec[id1])])

    # create a dataset for the column
    col_ds = ds(i, no3_col, no2_col, so4_col, fe_col, pH_col, ec_col)
    all_ds["$(i)"] = col_ds
    # save the dataset to a file
    
end

save("data/outflow_dataset.jld2", all_ds)
