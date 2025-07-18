using CSV, DataFrames, Dates, XLSX


# Function to read and process each sheet in the Excel file
function process_excel_sheets(file_path)
    # Dictionary to store DataFrames for each sheet
    sheets_data = Dict{String, DataFrame}()
    
    # Open the Excel file
    XLSX.openxlsx(file_path) do xf
        # Iterate over each sheet and read it as a DataFrame
        for sheet_name in XLSX.sheetnames(xf)
            println("Processing sheet: $sheet_name")
            data = XLSX.readtable(file_path, sheet_name)
            df = DataFrame(data)
            sheets_data[sheet_name] = df
            
            df.x = df."column length"
            df."S-2" = df."s-2 [nano mol/g]"
            df."S-2 r1" = df."T0rep1"
            df."S-2 r2" = df."T0rep2"
            df = select(df, "column","x", "S-2", "S-2 r1", "S-2 r2")
            CSV.write(joinpath(output_dir, "processed_AVS_data.csv"), df)
            sheets_data[sheet_name] = df
        end
    end
    return sheets_data
end

sheets_data = process_excel_sheets(datadir("exp_raw","solid_phase_sulfide.xlsx"))