
# --- Comparison with Eschenbach 2015 and Weymann et al. 2010---

println("\n--- Comparison with Eschenbach et al. (2015) ---")

df_rates = CSV.read("data/no3_rates.csv")
# Read external data
# Assuming data/external/eschenbach_etal_2015_S1.csv exists relative to project root
df_ext = CSV.read("data/external/eschenbach_etal_2015_S1.csv", DataFrame)

# Filter for FFA
df_ffa = filter(row -> occursin("FFA", row["Sample location"]), df_ext)

# Prepare result table
comparison_rows = NamedTuple{(:Location, :Depth, :Zone, :Type, :Parameter, :Value, :Unit), Tuple{String, Union{String, Missing}, Union{String, Missing}, String, String, Float64, String}}[]

# Atomic masses
const m_N = 14.007
const m_S = 32.065

for row in eachrow(df_ffa)
    loc = row["Sample location"]
    depth = row["Depth interval (m)"]
    zone = row["Aquifer zone"]
    
    # 1. Rates
    # Dcum(365) (mg N kg-1 yr-1) -> mol N kg-1 d-1
    val_dcum = row["Dcum(365) (mg N kg-1 yr-1)"]
    if !ismissing(val_dcum)
        rate = (val_dcum * 1e-3) / m_N / 365.0
        push!(comparison_rows, (Location=loc, Depth=depth, Zone=zone, Type="Rate", Parameter="Dcum(365)", Value=rate, Unit="mol N kg-1 d-1"))
    end
    
    # Dr in situ (µg N kg-1 d-1) -> mol N kg-1 d-1
    val_dr = row["Dr in situ (µg N kg-1 d-1)"]
    if !ismissing(val_dr)
        rate = (val_dr * 1e-6) / m_N
        push!(comparison_rows, (Location=loc, Depth=depth, Zone=zone, Type="Rate", Parameter="Dr in situ", Value=rate, Unit="mol N kg-1 d-1"))
    end

    # SFC (mg S kg-1 yr-1) -> mol S kg-1 d-1
    val_sfc = row["SFC (mg S kg-1 yr-1)"]
    if !ismissing(val_sfc)
        rate = (val_sfc * 1e-3) / m_S / 365.0
        push!(comparison_rows, (Location=loc, Depth=depth, Zone=zone, Type="Rate", Parameter="SFC", Value=rate, Unit="mol S kg-1 d-1"))
    end

    # 2. Amounts (Concentrations in solid phase)
    # SRC, SRCc, SRCs (mg N kg-1) -> mol N kg-1
    for (col_name, param_name) in [("SRC (mg N kg-1)", "SRC"), ("SRCc (mg N kg-1)", "SRCc"), ("SRCs (mg N kg-1)", "SRCs")]
        val_src = row[col_name]
        if !ismissing(val_src)
            amount = (val_src * 1e-3) / m_N
            push!(comparison_rows, (Location=loc, Depth=depth, Zone=zone, Type="Amount", Parameter=param_name, Value=amount, Unit="mol N kg-1"))
        end
    end
end

# Append Model Columns (Rates)
for (col_idx, rate) in model_rates_collection
    push!(comparison_rows, (Location="Column $col_idx", Depth=missing, Zone="Model", Type="Rate", Parameter="Model Rate", Value=rate, Unit="mol N kg-1 d-1"))
end

# Create DataFrame
df_comparison = DataFrame(comparison_rows)

println("Comparison Table:")
println(df_comparison)