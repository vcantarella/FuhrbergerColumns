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
        if occursin(lowercase(prefix), lowercase(n))
            return n
        end
    end
    return nothing
end

# --- Data Containers ---
rates_rows = []
conc_rows = []

# --- 1. Load Project Rates ---
df_project = CSV.read("data/no3_rates.csv", DataFrame)
for row in eachrow(df_project)
    push!(rates_rows, (
        Reference = "Current Project",
        Zone = "Non-Sulphidic",
        Depth_Mid = 17.5,
        Rate_Min = row.rate_mol_kg_d / 2.3, # Conservative Q10
        Rate_Max = row.rate_mol_kg_d / 1.2  # Optimistic Q10
    ))
end

# --- 1b. Load Project Concentrations ---
df_soli = CSV.read("data/soliphase_results.csv", DataFrame)
for row in eachrow(df_soli)
    if !ismissing(row["TOC [mol/kg]"])
        push!(conc_rows, (
            Reference = "Current Project",
            Zone = "Non-Sulphidic",
            Depth_Mid = 17.5,
            Type = "C_org",
            Value = row["TOC [mol/kg]"]
        ))
    end
end

# --- 2. Load Eschenbach et al. (2015) ---
df_esch = CSV.read("data/external/eschenbach_2013_2015.csv", DataFrame)
col_dcum = get_col_robust(df_esch, "Dcum")
col_insitu = get_col_robust(df_esch, "Dr in situ")
col_corg_esch = get_col_robust(df_esch, "Corg")
col_stot_esch = get_col_robust(df_esch, "Total-S")

for row in eachrow(df_esch)
    !occursin("FFA", string(row["Sample Location"])) && continue
    d_mid = parse_depth(row["Depth interval (m)"])
    ismissing(d_mid) && continue
    
    z = normalize_zone(row["Sediment Group"], "Eschenbach")
    ref = "Eschenbach et al.\n (2013,2015)"
    
    # Rates: Collect ANY available rate
    available_rates = Float64[]
    
    if !isnothing(col_insitu) && !ismissing(row[col_insitu])
        push!(available_rates, (row[col_insitu] * 1e-6) / m_N)
    end
    if !isnothing(col_dcum) && !ismissing(row[col_dcum])
        push!(available_rates, (row[col_dcum] * 1e-3) / m_N / 365.0)
    end
    
    if !isempty(available_rates)
        push!(rates_rows, (
            Reference=ref, 
            Zone=z, 
            Depth_Mid=d_mid, 
            Rate_Min=minimum(available_rates), 
            Rate_Max=maximum(available_rates)
        ))
    end
    
    # Concentrations
    if !isnothing(col_corg_esch) && !ismissing(row[col_corg_esch])
        push!(conc_rows, (Reference=ref, Zone=z, Depth_Mid=d_mid, Type="C_org", Value=(row[col_corg_esch] * 1e-3) / m_C))
    end
    if !isnothing(col_stot_esch) && !ismissing(row[col_stot_esch])
        push!(conc_rows, (Reference=ref, Zone=z, Depth_Mid=d_mid, Type="S_total", Value=(row[col_stot_esch] * 1e-3) / m_S))
    end
end

# --- 3. Load Weymann et al. (2010) ---
df_wey = CSV.read("data/external/weymann_et_al_2010.csv", DataFrame)
# Strict column finding for Weymann to avoid partial matches
cols_di = [n for n in names(df_wey) if startswith(n, "Di") || startswith(n, "Di\u200b")]
col_di = isempty(cols_di) ? nothing : first(cols_di)

cols_dmax = [n for n in names(df_wey) if startswith(n, "Dmax") || startswith(n, "Dmax\u200b")]
col_dmax = isempty(cols_dmax) ? nothing : first(cols_dmax)

col_corg_wey = get_col_robust(df_wey, "Org C")
col_stot_wey = get_col_robust(df_wey, "Total S")

for row in eachrow(df_wey)
    d_mid = parse_depth(row["Depth (m)"])
    ismissing(d_mid) && continue
    
    z = normalize_zone(row["Zone"], "Weymann")
    ref = "Weymann et al. (2010)"
    
    # Rates
    available_rates = Float64[]
    
    if !isnothing(col_di) && !ismissing(row[col_di])
        push!(available_rates, (row[col_di] * 1e-6) / m_N)
    end
    if !isnothing(col_dmax) && !ismissing(row[col_dmax])
        push!(available_rates, (row[col_dmax] * 1e-6) / m_N)
    end
    
    if !isempty(available_rates)
        push!(rates_rows, (
            Reference=ref, 
            Zone=z, 
            Depth_Mid=d_mid, 
            Rate_Min=minimum(available_rates), 
            Rate_Max=maximum(available_rates)
        ))
    end
    
    # Concentrations
    if !isnothing(col_corg_wey) && !ismissing(row[col_corg_wey])
        push!(conc_rows, (Reference=ref, Zone=z, Depth_Mid=d_mid, Type="C_org", Value=(row[col_corg_wey] * 1e-3) / m_C))
    end
    if !isnothing(col_stot_wey) && !ismissing(row[col_stot_wey])
        push!(conc_rows, (Reference=ref, Zone=z, Depth_Mid=d_mid, Type="S_total", Value=(row[col_stot_wey] * 1e-3) / m_S))
    end
end

df_rates = DataFrame(rates_rows)
df_conc = DataFrame(conc_rows)

# --- Plotting ---
f = Figure(size = (1200, 700))

# Layout: 3 Columns [Rates | Corg | S_total] + Legend
ax_rates = Axis(f[1, 1], xlabel = "Rate (mol N kg⁻¹ d⁻¹)",
 ylabel = "Depth (m)", title = "a. Denitrification Rates",
 yreversed = true, xscale = log10, titlealign= :left)
ax_c = Axis(f[1, 2], xlabel = "Corg (mol C kg⁻¹)",
 title = "b. Organic Carbon",
 yreversed = true, xscale = log10,
 titlealign = :left)
ax_s = Axis(f[1, 3], xlabel = "Stotal (mol S kg⁻¹)",
 title = "c. Total Sulfur",
  yreversed = true, xscale = log10,
  titlealign = :left)

linkyaxes!(ax_rates, ax_c, ax_s)
hideydecorations!(ax_c, grid = false)
hideydecorations!(ax_s, grid = false)

# Styling
zone_colors = Dict("Non-Sulphidic" => :darkgoldenrod4,
                                    "Sulphidic" => :gold1,
                                    "Transition Zone" => :darkkhaki)
ref_markers = Dict("Current Project" => :circle,
                                    "Eschenbach et al.\n (2013,2015)" => :rect,
                                    "Weymann et al. (2010)" => :diamond)
ref_markersisze = Dict("Current Project" => 19,
                                    "Eschenbach et al.\n (2013,2015)" => 16,
                                    "Weymann et al. (2010)" => 16)

linewidth = 5
strokewidth = 1.0
strokecolor = :grey21 
# 1. Rates Plot
for row in eachrow(df_rates)
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    msize = get(ref_markersisze, row.Reference, 16)
    # Range Bar (only if Min != Max)
    if row.Rate_Min != row.Rate_Max
        rangebars!(ax_rates, [row.Depth_Mid], [row.Rate_Min], [row.Rate_Max], direction = :x, color = (c, 0.8), linewidth = linewidth)
    end
    # Markers
    scatter!(ax_rates, [row.Rate_Min, row.Rate_Max], [row.Depth_Mid, row.Depth_Mid], color = c, marker = m, markersize = msize,
     strokewidth = strokewidth, strokecolor = strokecolor)
end

# 2. Corg Plot
for row in eachrow(filter(r -> r.Type == "C_org", df_conc))
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    msize = get(ref_markersisze, row.Reference, 16)
    scatter!(ax_c, row.Value, row.Depth_Mid, color = c, marker = m, markersize = msize,
     strokewidth = strokewidth, strokecolor = strokecolor)
end

# 3. Stotal Plot
for row in eachrow(filter(r -> r.Type == "S_total", df_conc))
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    msize = get(ref_markersisze, row.Reference, 16)
    scatter!(ax_s, row.Value, row.Depth_Mid, color = c, marker = m, markersize = msize,
     strokewidth = strokewidth, strokecolor = strokecolor)
end

# Legends
dataset_leg = [MarkerElement(marker = ref_markers[r], color = :black, markersize = 12) => r for r in keys(ref_markers)]
zone_leg = [MarkerElement(marker = :circle, color = zone_colors[z], markersize = 12) => z for z in sort(collect(keys(zone_colors)))]
Legend(f[1, 4], [first.(dataset_leg), first.(zone_leg)], [last.(dataset_leg), last.(zone_leg)], ["Dataset", "Zone"])

mkpath("plots")
save("plots/comparison_full_profile.png", f)
println("Full profile plot updated: plots/comparison_full_profile.png")

# --- Export Project Data with 10C Rates ---
println("\n--- Exporting Updated Project Data ---")
# Calculate 10C estimates
df_project.rate_min_10C = df_project.rate_mol_kg_d ./ 2.3
df_project.rate_max_10C = df_project.rate_mol_kg_d ./ 1.2

# Write back to CSV
CSV.write("data/no3_rates_final.csv", df_project)
println("Updated 'data/no3_rates.csv' with 10°C rate estimates.")