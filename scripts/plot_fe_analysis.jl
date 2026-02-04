using XLSX
using DataFrames
using CairoMakie
using Statistics

# Helper to parse sample names
function parse_sample_name(name)
    # Check for Fe pattern: "C1 0-4"
    m_fe = match(r"^C(\d)\s+(\d-\d)", name)
    if !isnothing(m_fe)
        return (column="C" * m_fe.captures[1], depth=m_fe.captures[2])
    end
    
    # Check for TOC pattern: "C1-0-4"
    m_toc = match(r"^C(\d)-(\d-\d)", name)
    if !isnothing(m_toc)
        return (column="C" * m_toc.captures[1], depth=m_toc.captures[2])
    end
    
    return nothing
end

println("--- Loading and Unifying Data ---")

# --- Load Fe Data ---
file_path_fe = "data/exp_raw/fe_analysis.xlsx"
sheet_name_fe = "fe_plate_1"
df_fe = DataFrame(XLSX.readtable(file_path_fe, sheet_name_fe))

# Select relevant columns and clean
df_fe_clean = select(df_fe, "name", "fe2+ [mol/kg]", "fe3+ [mol/kg]")
parsed_fe = parse_sample_name.(convert.(String, df_fe.name))
df_fe_clean.Column = [p !== nothing ? p.column : missing for p in parsed_fe]
df_fe_clean.Depth = [p !== nothing ? p.depth : missing for p in parsed_fe]
# Filter for valid C-columns
filter!(row -> !ismissing(row.Column), df_fe_clean)

# --- Load TOC Data ---
file_path_toc = "data/exp_raw/TOC_Analysis_Fuhrberg.xlsx"
sheet_name_toc = "Tabelle1"
df_toc = DataFrame(XLSX.readtable(file_path_toc, sheet_name_toc))

# Select relevant columns and clean
df_toc_clean = select(df_toc, "Name", "TOC%", "TIC%")
parsed_toc = parse_sample_name.(convert.(String, df_toc.Name))
df_toc_clean.Column = [p !== nothing ? p.column : missing for p in parsed_toc]
df_toc_clean.Depth = [p !== nothing ? p.depth : missing for p in parsed_toc]
# Filter for valid C-columns
filter!(row -> !ismissing(row.Column), df_toc_clean)

# --- Unify ---
# Join on Column and Depth
df_joined = outerjoin(df_fe_clean, df_toc_clean, on=[:Column, :Depth], makeunique=true)

# --- Process Data ---
# 1. Convert TOC/TIC/Carbon % to mol/kg
# 1% = 10 g/kg
# Molar mass of C = 12.011 g/mol
const M_C = 12.011

function percent_to_mol_kg(percent)
    if ismissing(percent)
        return missing
    end
    # Ensure input is Float64
    p_float = try 
        Float64(percent) 
    catch 
        missing 
    end
    if ismissing(p_float) return missing end
    
    return (p_float * 10) / M_C
end

df_joined[!, "TOC [mol/kg]"] = percent_to_mol_kg.(df_joined[!, "TOC%"])
df_joined[!, "TIC [mol/kg]"] = percent_to_mol_kg.(df_joined[!, "TIC%"])

# 2. Group by Column and Depth and calculate means
# Select columns to aggregate
cols_to_mean = ["fe2+ [mol/kg]", "fe3+ [mol/kg]", "TOC [mol/kg]", "TIC [mol/kg]"]

# Group and combine
df_final = combine(groupby(df_joined, [:Column, :Depth])) do gdf
    # Create a dictionary of results
    res = Dict{String, Any}()
    for col in cols_to_mean
        if col in names(gdf)
            vals = skipmissing(gdf[!, col])
            if isempty(vals)
                 res[col] = missing
            else
                 res[col] = mean(vals)
            end
        else
            res[col] = missing
        end
    end
    return DataFrame(res)
end

# Sort
sort!(df_final, [:Column, :Depth])

println("Processed Unified Data (Means, mol/kg):")
println(df_final)

using CSV
CSV.write("data/soliphase_results.csv", df_final)