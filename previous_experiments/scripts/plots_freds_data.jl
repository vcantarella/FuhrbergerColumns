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
df23_24 = CSV.File(datadir("exp_pro","freds_processed_data_23-24.csv")) |> DataFrame
df27_28 = CSV.File(datadir("exp_pro","freds_processed_data_27-28.csv")) |> DataFrame

data = Dict(
    "23-24" => df23_24,
    "27-28" => df27_28,
)

# Extract relevant sheets and columns
sheets = keys(data)
threshold = 0.1

for sheet in sheets
    sheet_data = data[sheet]
    
    # Extract columns and skip missing values
    ttm = sheet_data[!, "h"]
    no3_ulla = collect(skipmissing(sheet_data[sheet_data[!,:analytical_proc].=="Ulla","NO3-"]))
    no3_fred = collect(skipmissing(sheet_data[sheet_data[!,:analytical_proc].=="Fred","NO3-"]))
    so4 = collect(skipmissing(sheet_data[!,"SO4-2"]))
    no2_ulla = collect(skipmissing(sheet_data[sheet_data[!,:analytical_proc].=="Ulla","NO2-"])).*1e-3
    no2_fred = collect(skipmissing(sheet_data[sheet_data[!,:analytical_proc].=="Fred","NO2-"])).*1e-3
    # Filter ttm values corresponding to non-missing values
    time_ulla = ttm[sheet_data[!,:analytical_proc].=="Ulla"]
    time_fred = ttm[sheet_data[!,:analytical_proc].=="Fred"]
    time_no3_ulla = time_ulla[.!ismissing.(sheet_data[sheet_data[!,:analytical_proc].=="Ulla","NO3-"])]
    time_no3_fred = time_fred[.!ismissing.(sheet_data[sheet_data[!,:analytical_proc].=="Fred","NO3-"])]
    time_so4 = ttm[.!ismissing.(sheet_data[!,"SO4-2"])]
    time_no2_ulla = time_ulla[.!ismissing.(sheet_data[sheet_data[!,:analytical_proc].=="Ulla","NO2-"])]
    time_no2_fred = time_fred[.!ismissing.(sheet_data[sheet_data[!,:analytical_proc].=="Fred","NO2-"])]
    # time_no2 = ttm[.!ismissing.(sheet_data[!,"NO2-"])]
    
    # Create a figure
    fig = Figure(size = (800, 600))
    
    # Plot each ion concentration against ttm
    ax = Axis(fig[1, 1], title = "Ion Concentrations over ttm in $sheet", xlabel = "ttm (h)", ylabel = "Concentration")
    lines!(ax, time_no3_ulla, no3_ulla, label = "NO3- Ulla", color = "#60F5AC", linestyle = :dash)
    lines!(ax, time_no3_fred, no3_fred, label = "NO3- Fred", color = "#60F5AC")
    lines!(ax, time_so4, so4, label = "SO4-2", color = "#A3625A")
    lines!(ax, time_no2_ulla, no2_ulla, label = "NO2-", color = "#635250", linestyle = :dash)
    lines!(ax, time_no2_fred, no2_fred, label = "NO2- Fred", color = "#635250")
    scatter!(ax, time_no3_ulla, no3_ulla, color = "#60F5AC", markersize = 8)
    scatter!(ax, time_no3_fred, no3_fred, color = "#60F5AC", markersize = 8)
    scatter!(ax, time_no2_ulla, no2_ulla, color = "#635250", markersize = 8)
    scatter!(ax, time_no2_fred, no2_fred, color = "#635250", markersize = 8)
    #scatter!(ax, time_no3, no3, color = "#7FBDF9", markersize = 8)
    #scatter!(ax, time_so4, so4, color = "#A3625A", markersize = 8)
    #scatter!(ax, time_no2, no2, color = "#635250", markersize = 8)
    
    # Add grid lines
    
    # Define a threshold for low concentrations

    
    # Add a subplot for low concentrations
    ax_zoom = Axis(fig[2, 1], title = "Zoomed Ion Concentrations", xlabel = "ttm (h)", ylabel = "Concentration")
    lines!(ax_zoom, time_no3_fred[no3_fred .< threshold], no3_fred[no3_fred .< threshold], label = "NO3- Fred", color = "#7FBDF9")
    lines!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], label = "SO4-2", color = "#A3625A")
    lines!(ax_zoom, time_no2_ulla[no2_ulla .< threshold], no2_ulla[no2_ulla .< threshold], label = "NO2- Ulla", color = "#635250", linestyle = :dash)
    lines!(ax_zoom, time_no2_fred[no2_fred .< threshold], no2_fred[no2_fred .< threshold], label = "NO2- Fred", color = "#635250")
    scatter!(ax_zoom, time_no3_fred[no3_fred .< threshold], no3_fred[no3_fred .< threshold], color = "#7FBDF9", markersize = 11)
    scatter!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], color = "#A3625A", markersize = 11)
    scatter!(ax_zoom, time_no2_ulla[no2_ulla .< threshold], no2_ulla[no2_ulla .< threshold], color = "#635250", markersize = 11)
    scatter!(ax_zoom, time_no2_fred[no2_fred .< threshold], no2_fred[no2_fred .< threshold], color = "#635250", markersize = 11)
    
    # Add grid lines to zoomed plot
    
    # Add legend
    Legend(fig[1, 2], ax, "Ion Concentrations")
    
    # Save the figure
    save(plotsdir("plot_fres_$sheet.png"), fig)
end


# Create a figure
fig = Figure(size = (1200, 600))


# Store axes for linking
main_axes = []
zoom_axes = []

# Iterate over sheets and create subplots
for (i, sheet) in enumerate(sheets)
    sheet_data = data[sheet]
    
    # Extract columns and skip missing values
    ttm = sheet_data[!, "h"]
    no3 = collect(skipmissing(sheet_data[!,"NO3-"]))
    so4 = collect(skipmissing(sheet_data[!,"SO4-2"]))
    no2 = collect(skipmissing(sheet_data[!,"NO2-"])).*1e-3
    
    # Filter ttm values corresponding to non-missing values
    time_no3 = ttm[.!ismissing.(sheet_data[!,"NO3-"])]
    time_so4 = ttm[.!ismissing.(sheet_data[!,"SO4-2"])]
    time_no2 = ttm[.!ismissing.(sheet_data[!,"NO2-"])]
    
    # Create subplots
    ax = Axis(fig[1, i], xlabel = "ttm (h)", ylabel = "Concentration (mmol/L)",
    xlabelsize = 14, ylabelsize = 14, xticklabelsize = 12, yticklabelsize = 12)
    lines!(ax, time_no3, no3, label = "NO3-", color = "#7FBDF9")
    lines!(ax, time_so4, so4, label = "SO4-2", color = "#A3625A")
    lines!(ax, time_no2, no2, label = "NO2-", color = "#635250")
    scatter!(ax, time_no3, no3, color = "#7FBDF9", markersize = 11)
    scatter!(ax, time_so4, so4, color = "#A3625A", markersize = 11)
    scatter!(ax, time_no2, no2, color = "#635250", markersize = 11)
    
    # Add a subplot for low concentrations
    ax_zoom = Axis(fig[2, i], xlabel = "ttm (h)", ylabel = "Concentration (mmol/L)",
    xlabelsize = 14, ylabelsize = 14, xticklabelsize = 12, yticklabelsize = 12)
    lines!(ax_zoom, time_no3[no3 .< threshold], no3[no3 .< threshold], label = "NO3-", color = "#7FBDF9")
    lines!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], label = "SO4-2", color = "#A3625A")
    lines!(ax_zoom, time_no2[no2 .< threshold], no2[no2 .< threshold], label = "NO2-", color = "#635250")
    scatter!(ax_zoom, time_no3[no3 .< threshold], no3[no3 .< threshold], color = "#7FBDF9", markersize = 11)
    scatter!(ax_zoom, time_so4[so4 .< threshold], so4[so4 .< threshold], color = "#A3625A", markersize = 11)
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

for (label, layout) in zip(["a. 23-24", "b.27-28", "c."], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
end

# Add legend at the bottom of the figure
Legend(fig[3, 2], main_axes[1], "measurement", orientation = :horizontal)

# Save the figure
save(plotsdir("freds_plot_combined.png"), fig)


# pH and EC plots
# Create a figure
fig = Figure(size = (800, 600))
# pH Axis
ax_pH = Axis(fig[1, 1], xlabel = "ttm (h)", ylabel = "pH",
    xlabelsize = 14, ylabelsize = 14, xticklabelsize = 12, yticklabelsize = 12)
# EC Axis
ax_EC = Axis(fig[2, 1], xlabel = "ttm (h)", ylabel = "EC (µS/cm)",
    xlabelsize = 14, ylabelsize = 14, xticklabelsize = 12, yticklabelsize = 12)
# plot the data from each sheet
## create a color palette
palette = [:blue, :green]
## iterate over the sheets
for (i, sheet) in enumerate(sheets)
    sheet_data = data[sheet]
    ttm = sheet_data[!, "h"]
    pH = collect(skipmissing(sheet_data[!,"pH"]))
    EC = collect(skipmissing(sheet_data[!,"EC"]))
    time_pH = ttm[.!ismissing.(sheet_data[!,"pH"])]
    time_EC = ttm[.!ismissing.(sheet_data[!,"EC"])]
    lines!(ax_pH, time_pH, pH, label = sheet, color = palette[i])
    lines!(ax_EC, time_EC, EC, label = sheet, color = palette[i])
    scatter!(ax_pH, time_pH, pH, color = palette[i], markersize = 11)
    scatter!(ax_EC, time_EC, EC, color = palette[i], markersize = 11)
end
# Add legend
Legend(fig[3, 1], ax_pH, "columns", orientation = :horizontal)
# save
save(plotsdir("freds_pH_EC_plot.png"), fig)