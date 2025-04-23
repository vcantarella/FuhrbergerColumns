using DrWatson
@quickactivate "FuhrbergerColumns"
using XLSX
using DataFrames
using Dates
using CSV

# Define the path to the Excel file
excel_file_path = datadir("exp_raw", "experimental_data_Fred_raw.xlsx")


# Here you may include files from the source directory
# include(srcdir("dummy_src_file.jl"))

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)

# Function to read and process each sheet in the Excel file
function process_excel_sheets(file_path)
    # Dictionary to store DataFrames for each sheet
    sheets_data = Dict{String, DataFrame}()
    output_dir = datadir("exp_pro")
    
    # Open the Excel file
    XLSX.openxlsx(file_path) do xf
        # Iterate over each sheet and read it as a DataFrame
        for sheet_name in XLSX.sheetnames(xf)
            println("Processing sheet: $sheet_name")
            data = XLSX.readtable(file_path, sheet_name)
            df = DataFrame(data)
            sheets_data[sheet_name] = df
            df."h" = df."Time point"
            df."NO3-" = df."NO3 (micro mol/L)".*1e-3
            df."SO4-2" = df."SO4-2(mg/L)"./96.06
            df."NO2-" = df."NO2 (micromol/L)"
            df."EC" = df."EC (µS/cm)"
            df = select(df, "h", "NO3-","NO2-", "SO4-2", "pH", "EC", "analytical_proc")
            # arrange df by "h"
            sort!(df, :h)
            CSV.write(joinpath(output_dir, "freds_processed_data_$sheet_name.csv"), df)
                
            sheets_data[sheet_name] = df
        end
    end
    return sheets_data
end

# Process the Excel file and get the DataFrames
sheets_data = process_excel_sheets(excel_file_path)

# Example: Print the first few rows of each DataFrame
for (sheet_name, df) in sheets_data
    println("Sheet: $sheet_name")
    println(first(df, 5))  # Print the first 5 rows
end

#save the processed data to a file (serialized?)
save(datadir("exp_pro", "processed_data_fred.jld2"), "processed_data", sheets_data)