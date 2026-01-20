using CairoMakie
using DrWatson
using DataFrames
using CSV
using XLSX
using Dates

include("prepare_bc_model_export.jl")

tube_diam = 1.52e-3 # m
dvs = Dict(1 => 29, 2 => 27, 3 => 33.5, 4 => 33.5) # in cm

df_cp = CSV.read("data/DIC_DOC_processed.csv", DataFrame)
df_cp[!, "t"] .= NaN
for i in 1:4
    index = df_cp[!, "column"] .== i
    t_uncorrected = Dates.value.(Dates.Second.(df_cp[index, "avg_time"] - t0s[i]))
    df_cp[index, "t"] = t_uncorrected .-
        dvs[i]/100*π*tube_diam^2/4 ./ disch_function.(t_uncorrected, Ref(disch_ds[i].Q), Ref(disch_ds[i].end_times))
end

start_times_dict = Dict(i => [Dates.value(Dates.Second(t - t0s[i])) for t in start_time[column .== i]] for i in 1:4)
end_times_dict = Dict(i => [Dates.value(Dates.Second(t - t0s[i])) for t in end_time[column .== i]] for i in 1:4)

avg_times_dict = Dict(i => (start_times_dict[i] .+ end_times_dict[i]) ./ 2 .-
    dvs[i]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[i], Ref(disch_ds[i].Q), Ref(disch_ds[i].end_times)) for i in 1:4)

no3 = vcat(df[!,"no3- [mgN/L]"], df_tr[!,"no3- [mgN/L]"]) 
no3_ic = vcat(df[!,"no3_mg_L"], repeat([missing], size(df_tr[!,"no3- [mgN/L]"], 1))) 

df_exp_final = DataFrame(column = Int[], time = Float64[], NO3 = Float64[])

for i in 1:4
    id1 = column .== i
    t = avg_times_dict[i][.!ismissing.(no3_ic[id1])]
    t_cuv = avg_times_dict[i][.!ismissing.(no3[id1]) .& ismissing.(no3_ic[id1])]
    no3_v = vcat(convert.(Float64, no3_ic[id1][.!ismissing.(no3_ic[id1])])./62,
                convert.(Float64, no3[id1][.!ismissing.(no3[id1]) .& ismissing.(no3_ic[id1])])./14)
    
    a_times = vcat(t, t_cuv)
    idx = sortperm(a_times)
    
    for k in idx
        push!(df_exp_final, (i, a_times[k], no3_v[k]))
    end
end

CSV.write("data/experimental_data_export.csv", df_exp_final)
println("Simplified data export complete.")