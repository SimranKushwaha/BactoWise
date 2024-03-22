#Check if any overlapping rows (start and stop of one exists between start and stop of other) add overlap for both rows in status column
import pandas as pd

# Read the Excel file
input_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/6_Grouped_Data_v2.xlsx"
output_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/7_Grouped_Data_v3.xlsx"

df = pd.read_excel(input_file)

# Function to check for overlap
def check_overlap(row1, row2):
    if (row1['Start'] > row2['Start'] and row1['End'] < row2['End']) or \
       (row2['Start'] > row1['Start'] and row2['End'] < row1['End']):
        return True
    return False

# Initialize Status column with empty strings
df['Status'] = ''

# Iterate through rows to check for overlaps
for i, row1 in df.iterrows():
    for j, row2 in df.iterrows():
        if i != j and row1['Sequence_ID'] == row2['Sequence_ID'] and row1['Type'] == row2['Type']:
            if check_overlap(row1, row2):
                df.at[i, 'Status'] = 'Overlap'
                df.at[j, 'Status'] = 'Overlap'

# Save the updated DataFrame to a new Excel file
df.to_excel(output_file, index=False)