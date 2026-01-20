using DrWatson
using DataFrames, XLSX, Statistics
using Dates
using StaticArrays
import DataInterpolations as DI
import RegularizationTools as RT
using JLD2
using CSV

include("model_data_structures.jl")
# load analytical data
file_path = datadir("exp_raw", "ssexp_data.xlsx")
sheet_names = ["general_samples_p$i" for i in 1:4]
datas = [XLSX.readtable(file_path, sheet_name) for sheet_name in sheet_names]
dfs = [DataFrame(data) for data in datas]

df_no2 = DataFrame(XLSX.readtable(file_path, "standard_curve_no2"))
for df in dfs
    df[!, "no2_mmol_L"] = Vector{Union{Float64, Missing}}(missing, nrow(df))
    for (i, sample) in enumerate(df_no2[!, "Sample"])
        for j in 1:nrow(df)
            if df[j, "Sample"] == sample
                df[j, "no2_mmol_L"] = df_no2[i, "no2- [micromol/L]"] / 1000
            end
        end
    end
end
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
t0 = DateTime(2025, 09, 18, 18, 15)
re = r"^P(\d?)"
df[!, "column"] = match.(re, df[!,"Sample"]) .|> x -> x.captures[1] |> x -> parse(Int, x)
column = vcat(df[!,"column"], df_tr[!,"column"])

Q0 = Dict(
    1 =>  flow_rate[(.!ismissing.(flow_rate) .& (column .== 1))][1]./ 3600 .* 1e-6,
    2 => flow_rate[(.!ismissing.(flow_rate) .& (column .== 2))][1]./ 3600 .* 1e-6,
    3 => flow_rate[(.!ismissing.(flow_rate) .& (column .== 3))][1]./ 3600 .* 1e-6,
    4 => flow_rate[(.!ismissing.(flow_rate) .& (column .== 4))][1]./ 3600 .* 1e-6
)

dvs_t0 = Dict(1 => 37.5 + 38 + 48,
              2 => 37.5 + 38 + 25,
              3 => 37.5 + 38 + 23,
              4 => 13 + 38 + 22,
              ) # in cm
tube_diam = 0.152 # cm
dv_t0    = Dict(i => dvs_t0[i] * π * (tube_diam/2)^2 for i in 1:4)
t0s = Dict(i => t0 + Dates.Second(floor(Int64, dv_t0[i]*1e-6 / Q0[i])) for i in 1:4)

disch_ds = Dict()
for i in 1:4
    bool_index = (.!ismissing.(flow_rate) .& (column .== i))
    Q = flow_rate[bool_index] ./ 3600 .* 1e-6 # m3/s
    e_times = Dates.Second.(end_time[bool_index] .- t0s[i])
    disch_ds[i] = QData(Q, Dates.value.(e_times))
end

function disch_function(t, Qs, end_times)
    for i in eachindex(Qs)
        if t <= end_times[i]
            return Qs[i]
        end
    end
    return Qs[end]
end

transp_params = load("data/optimized_tracer_params_m2.jld2")
tracer_params = transp_params["tracer_params"]

v_da = Dict{Int, VDataA}()
D = 3.5*1e-2 
A = π * D^2 / 4 
for i in 1:4
    if i < 3
        ϕ = tracer_params[i][1]
    else
        ϕ = mean([tracer_params[k][1] for k in 1:3])
    end
    bool_index = (.!ismissing.(flow_rate) .& (column .== i))
    Q = flow_rate[bool_index] ./ 3600 .* 1e-6 
    e_times = Dates.value.(Dates.Second.(end_time[bool_index] .- t0s[i]))
    s_times = Dates.value.(Dates.Second.(start_time[bool_index] .- t0s[i]))
    v = Q./(ϕ * A) 
    v_da[i] = VDataA(v, s_times, e_times)
end

# Export velocity data (Simplified)
df_velocity = DataFrame(column = Int[], end_time = Float64[], velocity = Float64[])
for (col, data) in v_da
    for i in 1:length(data.v)
        push!(df_velocity, (col, data.end_times[i], data.v[i]))
    end
end
CSV.write("data/velocity_data_export.csv", df_velocity)

c_no3_bck1 = (124.56+124.61)/2/62 
t_switch_bck2 = DateTime(2025, 09, 22, 10, 20)
c_no3_bck2 = (114.12 + 117.25)/2/62 
t_switch_bck3 = DateTime(2025, 10, 02, 18, 15)
c_no3_bck3 = (127.79+127.51)/2/62 
t_switch_1_5mM = DateTime(2025, 10, 06, 17, 10)
c_no3_1_5mM = 94.31/62 
t_switch_1mM = DateTime(2025, 10, 12, 20, 35)
c_no3_1mM = (63.16 + 63.52)/2/62 

c_ins = Dict{Int64, CinData}()
for i in 1:4
    if i < 4
        c_no3 = c_no3_bck1*1e-3
        cins = [[c_no3],]
        t0switch = Float64[]
        
        switches = [t_switch_bck2, t_switch_bck3, t_switch_1_5mM, t_switch_1mM]
        concs = [c_no3_bck2*1e-3, c_no3_bck3*1e-3, c_no3_1_5mM*1e-3, c_no3_1mM*1e-3]
        
        for (ts, val) in zip(switches, concs)
            t_s = Dates.value(Dates.Second(ts - t0s[i]))
            t_s += dv_t0[i]*1e-6 / disch_function(t_s, disch_ds[i].Q, disch_ds[i].end_times)
            push!(t0switch, t_s)
            push!(cins, [val])
        end
        c_ins[i] = CinData(cins, t0switch)
    else
        c_ins[i] = CinData([[0.0]], Float64[])
    end
end

# Export inflow data (Simplified)
df_inflow = DataFrame(column = Int[], switch_time = Float64[], NO3 = Float64[])
for (col, data) in c_ins
    for i in 1:length(data.c_in)
        t_s = i <= length(data.t_in) ? data.t_in[i] : Inf
        push!(df_inflow, (col, t_s, data.c_in[i][1]))
    end
end
CSV.write("data/inflow_data_export.csv", df_inflow)