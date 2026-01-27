using CSV
using DataFrames
using Statistics
using CairoMakie

# --- Constants ---
const m_N = 14.007  # g/mol
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
        if startswith(lowercase(n), lowercase(prefix))
            return n
        end
    end
    return nothing
end

# --- Data Loading and Processing ---
plot_rows = []

# 1. Project Data
df_project = CSV.read("data/no3_rates.csv", DataFrame)
for row in eachrow(df_project)
    push!(plot_rows, (
        Reference = "Current Project",
        Zone = "Non-Sulphidic",
        Depth_Mid = 17.5,
        Rate_Min = row.rate_mol_kg_d / 2.3, # Conservative Q10
        Rate_Max = row.rate_mol_kg_d / 1.2  # Optimistic Q10
    ))
end

# 2. Eschenbach et al. (2015)
df_esch = CSV.read("data/external/eschenbach_2013_2015.csv", DataFrame)
col_dcum = get_col_robust(df_esch, "Dcum")
col_insitu = get_col_robust(df_esch, "Dr in situ")

for row in eachrow(df_esch)
    # Filter for FFA as requested
    loc = string(row["Sample Location"])
    !occursin("FFA", loc) && continue

    d_mid = parse_depth(row["Depth interval (m)"])
    ismissing(d_mid) && continue
    
    val_insitu = !isnothing(col_insitu) ? row[col_insitu] : missing
    val_dcum = !isnothing(col_dcum) ? row[col_dcum] : missing
    
    if !ismissing(val_insitu) && !ismissing(val_dcum)
        r_insitu = (val_insitu * 1e-6) / m_N
        r_dcum = (val_dcum * 1e-3) / m_N / 365.0
        push!(plot_rows, (
            Reference = "Eschenbach et al. (2015)",
            Zone = normalize_zone(row["Sediment Group"], "Eschenbach"),
            Depth_Mid = d_mid,
            Rate_Min = min(r_insitu, r_dcum),
            Rate_Max = max(r_insitu, r_dcum)
        ))
    end
end

# 3. Weymann et al. (2010)
df_wey = CSV.read("data/external/weymann_et_al_2010.csv", DataFrame)
col_di = get_col_robust(df_wey, "Di")
col_dmax = get_col_robust(df_wey, "Dmax")

for row in eachrow(df_wey)
    d_mid = parse_depth(row["Depth (m)"])
    ismissing(d_mid) && continue
    
    val_di = !isnothing(col_di) ? row[col_di] : missing
    val_dmax = !isnothing(col_dmax) ? row[col_dmax] : missing
    
    if !ismissing(val_di) && !ismissing(val_dmax)
        r_di = (val_di * 1e-6) / m_N
        r_dmax = (val_dmax * 1e-6) / m_N
        push!(plot_rows, (
            Reference = "Weymann et al. (2010)",
            Zone = normalize_zone(row["Zone"], "Weymann"),
            Depth_Mid = d_mid,
            Rate_Min = min(r_di, r_dmax),
            Rate_Max = max(r_di, r_dmax)
        ))
    end
end

df_plot = DataFrame(plot_rows)

# --- Plotting ---
f = Figure(size = (900, 700))
ax = Axis(f[1, 1], xlabel = "Rate (mol N kg⁻¹ d⁻¹)",
 ylabel = "Depth (m)",
 yreversed = true, xscale = log10)

zone_colors = Dict("Non-Sulphidic" => Makie.wong_colors()[1], "Sulphidic" => Makie.wong_colors()[2], "Transition Zone" => Makie.wong_colors()[3])
ref_markers = Dict("Current Project" => :circle, "Eschenbach et al. (2015)" => :rect, "Weymann et al. (2010)" => :diamond)

for row in eachrow(df_plot)
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    rangebars!(ax, [row.Depth_Mid], [row.Rate_Min], [row.Rate_Max], direction = :x, color = (c, 0.6), linewidth = 2)
    scatter!(ax, [row.Rate_Min, row.Rate_Max], [row.Depth_Mid, row.Depth_Mid], color = c, marker = m, markersize = 10, strokewidth = 0.5, strokecolor = :black)
end

dataset_leg = [MarkerElement(marker = ref_markers[r], color = :black, markersize = 12) => r for r in keys(ref_markers)]
zone_leg = [MarkerElement(marker = :circle, color = zone_colors[z], markersize = 12) => z for z in sort(collect(keys(zone_colors)))]
Legend(f[1, 2], [first.(dataset_leg), first.(zone_leg)], [last.(dataset_leg), last.(zone_leg)], ["Dataset", "Zone"])

mkpath("plots")
save("plots/comparison_rates_styled.png", f)
println("Plot updated: plots/comparison_rates_styled.png")
