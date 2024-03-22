#Concatenate all the excels generated from the GFF files
import os
import pandas as pd

# Path to the folder containing Excel files
folder_path = '/Users/simrankushwaha/Documents/PostDoc/Annotations/474/GFF+Excel'

# Get a list of all Excel files in the folder
excel_files = [file for file in os.listdir(folder_path) if file.endswith('.xlsx')]

# Initialize a DataFrame to hold merged data
merged_data = pd.DataFrame()

# Flag to skip headers from the second file onwards
skip_header = False

for file in excel_files:
    # Read each Excel file into a DataFrame
    df = pd.read_excel(os.path.join(folder_path, file))
    
    # Skip header if not the first file
    if skip_header:
        # Concatenate data, but keep the header row
        merged_data = pd.concat([merged_data, df], ignore_index=True)
    else:
        # Concatenate data including the header row for the first file
        merged_data = pd.concat([merged_data, df])
        skip_header = True  # Set the flag to True after the first file

# Sort the merged data by 'Sequence ID' and then by 'Start'
merged_data = merged_data.sort_values(by=['Sequence_Name', 'Sequence_ID', 'Start'])

# Write the merged data to a new Excel file
output_file_path = '/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/2_Merged_Files.xlsx'
merged_data.to_excel(output_file_path, index=False)
