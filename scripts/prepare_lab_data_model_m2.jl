using CairoMakie
using DrWatson
using DataFrames
using CSV
using XLSX
using Dates

# 
include("prepare_bc_model.jl")


# Dead volumes per column. This refer only to the dead volume of the outflow tube
dvs = Dict(1 => 35,
           2 => 29,
           3 => 34,
           4 => 30) # in cm
tube_diam = 1.52e-3 # m

column = convert.(Int64,df[!,"column"])
start_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0sdict[1])) for t in start_time[column .== 1]],
                        2 => [Dates.value(Dates.Second(t - t0sdict[2])) for t in start_time[column .== 2]],
                        3 => [Dates.value(Dates.Second(t - t0sdict[3])) for t in start_time[column .== 3]],
                        4 => [Dates.value(Dates.Second(t - t0sdict[4])) for t in start_time[column .== 4]])
end_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0sdict[1])) for t in end_time[column .== 1]],
                      2 => [Dates.value(Dates.Second(t - t0sdict[2])) for t in end_time[column .== 2]],
                      3 => [Dates.value(Dates.Second(t - t0sdict[3])) for t in end_time[column .== 3]],
                      4 => [Dates.value(Dates.Second(t - t0sdict[4])) for t in end_time[column .== 4]])

avg_times_dict = Dict(1 => (start_times_dict[1] .+ end_times_dict[1]) ./ 2 .-
    dvs[1]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[1], Ref(disch_ds[1].Q), Ref(disch_ds[1].end_times)),
    2 => (start_times_dict[2] .+ end_times_dict[2]) ./ 2 .-
    dvs[2]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[2], Ref(disch_ds[2].Q), Ref(disch_ds[2].end_times)),
    3 => (start_times_dict[3] .+ end_times_dict[3]) ./ 2 .- 
    dvs[3]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[3], Ref(disch_ds[3].Q), Ref(disch_ds[3].end_times)),
    4 => (start_times_dict[4] .+ end_times_dict[4]) ./ 2 .-
    dvs[4]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[4], Ref(disch_ds[4].Q), Ref(disch_ds[4].end_times)),
)

# Now we have the average times for each column, we can create the datasets per column
no3 = df[!,"NO3-"] # concentration in mM
no2 = df[!,"NO2-"] # concentration in mM
so4 = df[!,"SO4-2"] # concentration in mM
fe = df[!,"Fe2+"] # concentration in mM
pH = df[!,"pH"] # pH values
ec = df[!,"EC"] # electrical conductivity in μS/cm


# For each column, we create a dataset with the concentrations and the average times
all_ds = Dict{Int, ds}()
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
    all_ds[i] = col_ds
    # save the dataset to a file
    
end

# Making a Figure for each column with the main variables
colors = Dict(
    "NO3-" => :blue,
    "NO2-" => :orange,
    "Fe2+" => :red,
    "SO4-2" => :green
)
for i in 1:4
    fig = Figure()
    ax1 = Axis(fig[1,1], title="Nitrogen species in column $i",
    xlabel ="Time (d)", ylabel="Conc. (μM)",)
    ax2 = Axis(fig[2,1], title="Iron species in column $i",
    xlabel ="Time (d)", ylabel="Conc. (μM)",)
    ax3 = Axis(fig[3,1], title = "Sulfate in column $i",
    xlabel ="Time (d)", ylabel="Conc. (mM)",)
    lines!(ax1, all_ds[i].no3.t ./ 86400, all_ds[i].no3.conc, label="NO₃⁻", color=colors["NO3-"])
    lines!(ax1, all_ds[i].no2.t ./ 86400, all_ds[i].no2.conc, label="NO₂⁻", color=colors["NO2-"])
    lines!(ax2, all_ds[i].fe.t ./ 86400, all_ds[i].fe.conc, label="Fe²⁺", color=colors["Fe2+"])
    lines!(ax3, all_ds[i].so4.t ./ 86400, all_ds[i].so4.conc ./ 1000, label="SO₄²⁻", color=colors["SO4-2"])
    scatter!(ax1, all_ds[i].no3.t ./ 86400, all_ds[i].no3.conc, color=colors["NO3-"], label = "NO₃⁻", markersize=8)
    scatter!(ax1, all_ds[i].no2.t ./ 86400, all_ds[i].no2.conc, color=colors["NO2-"], label = "NO₂⁻", markersize=8)
    scatter!(ax2, all_ds[i].fe.t ./ 86400, all_ds[i].fe.conc, color=colors["Fe2+"], label = "Fe²⁺", markersize=8)
    scatter!(ax3, all_ds[i].so4.t ./ 86400, all_ds[i].so4.conc ./ 1000, color=colors["SO4-2"], label = "SO₄²⁻", markersize=8)
    axislegend(ax1, position = :rt, merge = true)
    axislegend(ax2, position = :rt, merge = true) 
    axislegend(ax3, position = :rt, merge = true)
    ylims!(ax3, 2, 2.8)
    resize_to_layout!(fig)
    save("plots/raw_data_species_column_$i.png", fig)

end


# Making a Figure for plotting the flow velocity for each column
fig_v = Figure()
ax = Axis(fig_v[1,1], title="Flow velocity in the model",
xlabel ="Time (d)", ylabel="Flow velocity (m/s)")
flowt = 15:0.0001:20
colors = [:blue, :orange, :green, :red]
function vfunc(t, vds::VData)
    for i in eachindex(vds.v)
        if t <= vds.end_times[i]
            return vds.v[i]
        end
    end
    return vds.v[end]
end
for i in 1:4
    vds = v_ds[i]
    # Calculate the flow velocity for the given time points
    flowv = [vfunc(t*24*60*60, vds) for t in flowt]
    # Create a line plot for the flow velocity
    lines!(ax, flowt, flowv, label="Column $i", color=colors[i],
    linestyle=:dash, linewidth=3)
end
axislegend(ax, position = :rt)
fig_v
save("plots/flow_velocity_model.png", fig_v)

println("Data preparation complete. Datasets for each column are ready.")