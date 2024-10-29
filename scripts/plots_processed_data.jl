using DrWatson
@quickactivate "FuhrbergerColumns"
using CairoMakie
using DataFrames
using CSV
# Here you may include files from the source directory
# include(srcdir("dummy_src_file.jl"))

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Have fun with your new project!

You can help us improve DrWatson by opening
issues on GitHub, submitting feature requests,
or even opening your own Pull Requests!
"""
)
# Load the data
dfA = CSV.File(datadir("exp_pro","processed_data_A.csv")) |> DataFrame
dfB = CSV.File(datadir("exp_pro","processed_data_B.csv")) |> DataFrame
dfC = CSV.File(datadir("exp_pro","processed_data_C.csv")) |> DataFrame

data = Dict(
    "A" => dfA,
    "B" => dfB,
    "C" => dfC
)

# Extract relevant sheets and columns
sheets = keys(data)
threshold = 0.1

for sheet in sheets
    sheet_data = data[sheet]
    
    # Extract columns and skip missing values
    time = sheet_data[!, "h"]
    no3 = collect(skipmissing(sheet_data[!,"NO3-"]))
    so4 = collect(skipmissing(sheet_data[!,"SO4-2"]))
    s2 = collect(skipmissing(sheet_data[!,"S-2"]))
    fe = collect(skipmissing(sheet_data[!,"Fe"]))
    no2 = collect(skipmissing(sheet_data[!,"NO2-"])).*1e-3
    # Filter time values corresponding to non-missing values

    time_no3 = time[.!ismissing.(sheet_data[!,"NO3-"])]
    time_so4 = time[.!ismissing.(sheet_data[!,"SO4-2"])]
    time_s2 = time[.!ismissing.(sheet_data[!,"S-2"])]
    time_fe = time[.!ismissing.(sheet_data[!,"Fe"])]
    time_no2 = time[.!ismissing.(sheet_data[!,"NO2-"])]
    
    # Create a figure
    fig = Figure(size = (800, 600))
    
    # Plot each ion concentration against time
    ax = Axis(fig[1, 1], title = "Ion Concentrations over Time in $sheet", xlabel = "Time (h)", ylabel = "Concentration")
    lines!(ax, time_no3, no3, label = "NO3-", color = "#60F5AC")
    lines!(ax, time_so4, so4, label = "SO4-2", color = "#A3625A")
    lines!(ax, time_s2, s2, label = "S-2", color = "#E35644")
    lines!(ax, time_fe, fe, label = "Fe", color = "#5A8E74")
    lines!(ax, time_no2, no2, label = "NO2-", color = "#635250")
    scatter!(ax, time_no3, no3, color = "#60F5AC", markersize = 8)
    scatter!(ax, time_so4, so4, color = "#A3625A", markersize = 8)
    scatter!(ax, time_s2, s2, color = "#E35644", markersize = 8)
    scatter!(ax, time_fe, fe, color = "#5A8E74", markersize = 8)
    scatter!(ax, time_no2, no2, color = "#635250", markersize = 8)
    
    # Add grid lines
    
    # Define a threshold for low concentrations

    
    # Add a subplot for low concentrations
    ax_zoom = Axis(fig[2, 1], title = "Zoomed Ion Concentrations", xlabel = "Time (h)", ylabel = "Concentration")
    lines!(ax_zoom, time_no3[no3 .< threshold], no3[no3 .< threshold], label = "NO3-", color = "#7FBDF9")
    lines!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], label = "SO4-2", color = "#A3625A")
    lines!(ax_zoom, time_s2[s2 .< threshold], s2[s2 .< threshold], label = "S-2", color = "#E35644")
    lines!(ax_zoom, time_fe[fe .< threshold], fe[fe .< threshold], label = "Fe", color = "#5A8E74")
    lines!(ax_zoom, time_no2[no2 .< threshold], no2[no2 .< threshold], label = "NO2-", color = "#635250")
    scatter!(ax_zoom, time_no3[no3 .< threshold], no3[no3 .< threshold], color = "#7FBDF9", markersize = 11)
    scatter!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], color = "#A3625A", markersize = 11)
    scatter!(ax_zoom, time_s2[s2 .< threshold], s2[s2 .< threshold], color = "#E35644", markersize = 11)
    scatter!(ax_zoom, time_fe[fe .< threshold], fe[fe .< threshold], color = "#5A8E74", markersize = 11)
    scatter!(ax_zoom, time_no2[no2 .< threshold], no2[no2 .< threshold], color = "#635250", markersize = 11)
    
    # Add grid lines to zoomed plot
    
    # Add legend
    Legend(fig[1, 2], ax, "Ion Concentrations")
    
    # Save the figure
    save(plotsdir("plot_$sheet.png"), fig)
end


# Create a figure
fig = Figure(size = (1600, 600))


# Store axes for linking
main_axes = []
zoom_axes = []

# Iterate over sheets and create subplots
for (i, sheet) in enumerate(sheets)
    sheet_data = data[sheet]
    
    # Extract columns and skip missing values
    time = sheet_data[!, "h"]
    no3 = collect(skipmissing(sheet_data[!,"NO3-"]))
    so4 = collect(skipmissing(sheet_data[!,"SO4-2"]))
    s2 = collect(skipmissing(sheet_data[!,"S-2"]))
    fe = collect(skipmissing(sheet_data[!,"Fe"]))
    no2 = collect(skipmissing(sheet_data[!,"NO2-"])).*1e-3
    
    # Filter time values corresponding to non-missing values
    time_no3 = time[.!ismissing.(sheet_data[!,"NO3-"])]
    time_so4 = time[.!ismissing.(sheet_data[!,"SO4-2"])]
    time_s2 = time[.!ismissing.(sheet_data[!,"S-2"])]
    time_fe = time[.!ismissing.(sheet_data[!,"Fe"])]
    time_no2 = time[.!ismissing.(sheet_data[!,"NO2-"])]
    
    # Create subplots
    ax = Axis(fig[1, i], xlabel = "Time (h)", ylabel = "Concentration (mmol/L)",
    xlabelsize = 14, ylabelsize = 14, xticklabelsize = 12, yticklabelsize = 12)
    lines!(ax, time_no3, no3, label = "NO3-", color = "#7FBDF9")
    lines!(ax, time_so4, so4, label = "SO4-2", color = "#A3625A")
    lines!(ax, time_s2, s2, label = "S-2", color = "#E35644")
    lines!(ax, time_fe, fe, label = "Fe", color = "#5A8E74")
    lines!(ax, time_no2, no2, label = "NO2-", color = "#635250")
    scatter!(ax, time_no3, no3, color = "#60F5AC", markersize = 11)
    scatter!(ax, time_so4, so4, color = "#A3625A", markersize = 11)
    scatter!(ax, time_s2, s2, color = "#E35644", markersize = 11)
    scatter!(ax, time_fe, fe, color = "#5A8E74", markersize = 11)
    scatter!(ax, time_no2, no2, color = "#635250", markersize = 11)
    
    # Add a subplot for low concentrations
    ax_zoom = Axis(fig[2, i], xlabel = "Time (h)", ylabel = "Concentration (mmol/L)",
    xlabelsize = 14, ylabelsize = 14, xticklabelsize = 12, yticklabelsize = 12)
    lines!(ax_zoom, time_no3[no3 .< threshold], no3[no3 .< threshold], label = "NO3-", color = "#7FBDF9")
    lines!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], label = "SO4-2", color = "#A3625A")
    lines!(ax_zoom, time_s2[s2 .< threshold], s2[s2 .< threshold], label = "S-2", color = "#E35644")
    lines!(ax_zoom, time_fe[fe .< threshold], fe[fe .< threshold], label = "Fe", color = "#5A8E74")
    lines!(ax_zoom, time_no2[no2 .< threshold], no2[no2 .< threshold], label = "NO2-", color = "#635250")
    scatter!(ax_zoom, time_no3[no3 .< threshold], no3[no3 .< threshold], color = "#60F5AC", markersize = 11)
    scatter!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], color = "#A3625A", markersize = 11)
    scatter!(ax_zoom, time_s2[s2 .< threshold], s2[s2 .< threshold], color = "#E35644", markersize = 11)
    scatter!(ax_zoom, time_fe[fe .< threshold], fe[fe .< threshold], color = "#5A8E74", markersize = 11)
    scatter!(ax_zoom, time_no2[no2 .< threshold], no2[no2 .< threshold], color = "#635250", markersize = 11)

    # Store axes for linking
    push!(main_axes, ax)
    push!(zoom_axes, ax_zoom)
end

# Link x-axes and y-axes of subplots horizontally
linkxaxes!(main_axes...)
linkyaxes!(main_axes...)
linkxaxes!(zoom_axes...)
linkyaxes!(zoom_axes...)

ga = fig[:, 1] = GridLayout()
gb = fig[:, 2] = GridLayout()
gc = fig[:, 3] = GridLayout()

for (label, layout) in zip(["a.", "b.", "c."], [ga, gb, gc])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
end

# Add legend at the bottom of the figure
Legend(fig[3, 2], main_axes[1], "measurement", orientation = :horizontal)

# Save the figure
save(plotsdir("plot_combined.png"), fig)