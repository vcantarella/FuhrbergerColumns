using DataFrames, XLSX
using DrWatson
using DataFrames
using XLSX
using Dates
using Statistics

# load analytical data
file_path = datadir("exp_raw", "ssexp_data.xlsx")
dfs = DataFrame[]
for col in ["p1", "p2", "p3", "p4"]
    sheet_name = "general_samples_$col"
    data = XLSX.readtable(file_path, sheet_name)
    df_temp = DataFrame(data)
    df_temp[!, "column"] .= parse(Int, last(col))  # Add a column indicating the column number
    push!(dfs, df_temp)
end
list_of_samples = String[]
for df in dfs
    samples = df[.!ismissing.(df[!, "no3- [mgN/L]"]), "Sample"]
    append!(list_of_samples, samples)
end
