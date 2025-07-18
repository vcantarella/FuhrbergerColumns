using DrWatson
@quickactivate "FuhrbergerColumns"
using XLSX
using DataFrames
using Dates
using CSV

# Define the path to the Excel file
excel_file_path = datadir("exp_raw", "Final_Table.xlsx")


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
    
    # Open the Excel file
    XLSX.openxlsx(file_path) do xf
        # Iterate over each sheet and read it as a DataFrame
        for sheet_name in XLSX.sheetnames(xf)
            println("Processing sheet: $sheet_name")
            data = XLSX.readtable(file_path, sheet_name)
            df = DataFrame(data)
            sheets_data[sheet_name] = df
            
            # Format the "Day" column as a date
            if "Day" in names(df)
                df.Day = Date.(df.Day)  # Adjust the date format as needed
            end
            
            # Format the "Start Time" and "End Time" columns as time
            if "Start time" in names(df)
                df."Start time" = Time.(string.(df."Start time"), "HH:MM:SS")  # Adjust the time format as needed
            end
            if "End Time" in names(df)
                df."End Time" = Time.(string.(df."End Time"), "HH:MM:SS")  # Adjust the time format as needed
            end
            
            # # Calculate the average date of sampling
            if "Day" in names(df) && "Start time" in names(df) && "End Time" in names(df)
                df."Average Date of Sampling" = [calculate_average_date_of_sampling(day, start_time, end_time) for (day, start_time, end_time) in zip(df.Day, df."Start time", df."End Time")]
            end

            if "Average Date of Sampling" in names(df)
                df."t" = Dates.value.(df."Average Date of Sampling" .- DateTime("2024-06-20T00:00:00"))./(24*60*60*1000)
                df."h" = Dates.value.(df."Average Date of Sampling" .- DateTime("2024-06-20T00:00:00"))./(60*60*1000)
            end
            
            sheets_data[sheet_name] = df
        end
    end
    return sheets_data
end

# Function to calculate the average date of sampling
function calculate_average_date_of_sampling(day::Date, start_time::Time, end_time::Time)
    day_start = day
    day_end = day
    if end_time < start_time
        day_end += Day(1)  # Adjust for day turnover
    end
    time_start = DateTime(day_start, start_time)
    time_end = DateTime(day_end, end_time)
    average_time = time_start + (time_end - time_start) / 2
    return average_time
end

function strip_and_process_data(sheets_data::Dict{String, DataFrame})
    # Processing the concentration data and exporting the tables (processed data) to csv files
    # Define the path to the output directory
    output_dir = datadir("exp_pro")

    processed_data = Dict{String, DataFrame}()
    # convert Fe and S-2 from mg/L to mmol/L, simplify the column names and exporting the relevant columns to a csv file
    for (sheet_name, df) in sheets_data
        df."Fe" = @. ifelse(~ismissing(df."Fe (mg/L)") & (df."Fe (mg/L)"< 0), 0, df."Fe (mg/L)" ./ 55.845)
        df."S-2" = @. ifelse(~ismissing(df."S2- (mg/L)") & (df."S2- (mg/L)"<0),0,df."S2- (mg/L)" ./ 32.065)
        df."NO3-" = df."NO3- (mmol/L)"
        df."SO4-2" = df."SO42- (mmol/L)"
        df."NO2-" = df."NO2-(micromol/L)"
        df = select(df, "Sample","t", "h", "Fe", "S-2", "NO3-","NO2-", "SO4-2", "pH", "EC")
        CSV.write(joinpath(output_dir, "processed_data_$sheet_name.csv"), df)
        processed_data[sheet_name] = df
    end
    return processed_data
end

# Process the Excel file and get the DataFrames
sheets_data = process_excel_sheets(excel_file_path)

# Example: Print the first few rows of each DataFrame
for (sheet_name, df) in sheets_data
    println("Sheet: $sheet_name")
    println(first(df, 5))  # Print the first 5 rows
end

df = sheets_data["A"]

# Processing the concentration data and exporting the tables (processed data) to csv files
# Define the path to the output directory
output_dir = datadir("exp_pro")

processed_data = strip_and_process_data(sheets_data)
# convert Fe and S-2 from mg/L to mmol/L, simplify the column names and exporting the relevant columns to a csv file
# Example: Print the first few rows of each DataFrame
for (sheet_name, df) in processed_data
    println("Sheet: $sheet_name")
    println(first(df, 5))  # Print the first 5 rows
end

#save the processed data to a file (serialized?)
save(datadir("exp_pro", "processed_data.jld2"), "processed_data", processed_data)