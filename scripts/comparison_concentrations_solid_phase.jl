
using CSV
using DataFrames
using Statistics
using CairoMakie

# --- Constants ---
const m_S = 32.065  # g/mol
const m_C = 12.011  # g/mol

# --- Helpers ---
function parse_depth(d_str)
    ismissing(d_str) && return missing
    try
        parts = split(string(d_str), "-")
        return length(parts) == 2 ? (parse(Float64, parts[1]) + parse(Float64, parts[2])) / 2 : parse(Float64, d_str)
    catch
        return missing
    end
end

function normalize_zone(z_raw, ref_dataset)
    z = lowercase(string(z_raw))
    if occursin("weymann", lowercase(ref_dataset))
        return occursin("autotrophic", z) ? "Sulphidic" : occursin("heterotrophic", z) ? "Non-Sulphidic" : titlecase(z)
    end
    return occursin("non-sulphidic", z) || occursin("non sulphidic", z) ? "Non-Sulphidic" : 
           occursin("transition", z) ? "Transition Zone" : 
           occursin("sulphidic", z) ? "Sulphidic" : titlecase(z)
end

function get_col_robust(df, prefix)
    for n in names(df)
        # Check startswith or contains for robustness
        if occursin(lowercase(prefix), lowercase(n))
            return n
        end
    end
    return nothing
end

# --- Data Loading ---
plot_rows = []

# 1. Eschenbach et al. (2015)
df_esch = CSV.read("data/external/eschenbach_2013_2015.csv", DataFrame)
col_corg_esch = get_col_robust(df_esch, "Corg")
col_stot_esch = get_col_robust(df_esch, "Total-S") # Eschenbach uses "Total-S"

for row in eachrow(df_esch)
    # Filter for FFA
    loc = string(row["Sample Location"])
    !occursin("FFA", loc) && continue

    d_mid = parse_depth(row["Depth interval (m)"])
    ismissing(d_mid) && continue

    z = normalize_zone(row["Sediment Group"], "Eschenbach")
    
    # Corg
    val_c = !isnothing(col_corg_esch) ? row[col_corg_esch] : missing
    if !ismissing(val_c)
        push!(plot_rows, (Reference="Eschenbach et al. (2015)", Zone=z, Depth_Mid=d_mid, Type="C_org", Value=(val_c * 1e-3) / m_C))
    end
    
    # S_total
    val_s = !isnothing(col_stot_esch) ? row[col_stot_esch] : missing
    if !ismissing(val_s)
        push!(plot_rows, (Reference="Eschenbach et al. (2015)", Zone=z, Depth_Mid=d_mid, Type="S_total", Value=(val_s * 1e-3) / m_S))
    end
end

# 2. Weymann et al. (2010)
df_wey = CSV.read("data/external/weymann_et_al_2010.csv", DataFrame)
col_corg_wey = get_col_robust(df_wey, "Org C")
col_stot_wey = get_col_robust(df_wey, "Total S")

for row in eachrow(df_wey)
    d_mid = parse_depth(row["Depth (m)"])
    ismissing(d_mid) && continue

    z = normalize_zone(row["Zone"], "Weymann")
    
    # Corg
    val_c = !isnothing(col_corg_wey) ? row[col_corg_wey] : missing
    if !ismissing(val_c)
        push!(plot_rows, (Reference="Weymann et al. (2010)", Zone=z, Depth_Mid=d_mid, Type="C_org", Value=(val_c * 1e-3) / m_C))
    end
    
    # S_total
    val_s = !isnothing(col_stot_wey) ? row[col_stot_wey] : missing
    if !ismissing(val_s)
        push!(plot_rows, (Reference="Weymann et al. (2010)", Zone=z, Depth_Mid=d_mid, Type="S_total", Value=(val_s * 1e-3) / m_S))
    end
end

df_plot = DataFrame(plot_rows)

# --- Plotting ---
f = Figure(size = (1000, 600))

# Shared Y-axis (Depth), Separate X-axes (Concentrations)
# Grid: [Ax1, Ax2, Legend]
ax_c = Axis(f[1, 1], 
    xlabel = "Organic Carbon (mol C kg⁻¹)", 
    ylabel = "Depth (m)", 
    title = "Organic Carbon", 
    yreversed = true,
    xscale = log10
)
ax_s = Axis(f[1, 2], 
    xlabel = "Total Sulfur (mol S kg⁻¹)", 
    # ylabel = "Depth (m)", # Shared, hide label
    title = "Total Sulfur", 
    yreversed = true,
    xscale = log10
)
linkyaxes!(ax_c, ax_s)
hideydecorations!(ax_s, grid = false)

# Styling
zone_colors = Dict("Non-Sulphidic" => Makie.wong_colors()[1], "Sulphidic" => Makie.wong_colors()[2], "Transition Zone" => Makie.wong_colors()[3])
ref_markers = Dict("Eschenbach et al. (2015)" => :rect, "Weymann et al. (2010)" => :diamond)

# Plot Loops
# C_org
sub_c = filter(r -> r.Type == "C_org", df_plot)
for row in eachrow(sub_c)
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    scatter!(ax_c, row.Value, row.Depth_Mid, color = c, marker = m, markersize = 12, strokewidth = 0.5, strokecolor = :black)
end

# S_total
sub_s = filter(r -> r.Type == "S_total", df_plot)
for row in eachrow(sub_s)
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    scatter!(ax_s, row.Value, row.Depth_Mid, color = c, marker = m, markersize = 12, strokewidth = 0.5, strokecolor = :black)
end

# Legends
# We add "Current Project" as a placeholder or omit it if no data yet. 
# User: "I dont have the data yet... but I will add it shortly." -> Omit for now to avoid confusion.
dataset_leg = [MarkerElement(marker = ref_markers[r], color = :black, markersize = 12) => r for r in keys(ref_markers)]
zone_leg = [MarkerElement(marker = :circle, color = zone_colors[z], markersize = 12) => z for z in sort(collect(keys(zone_colors)))]
Legend(f[1, 3], [first.(dataset_leg), first.(zone_leg)], [last.(dataset_leg), last.(zone_leg)], ["Dataset", "Zone"])

mkpath("plots")
save("plots/comparison_concentrations_styled.png", f)
println("Plot created: plots/comparison_concentrations_styled.png")
