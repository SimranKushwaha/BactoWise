# If the Best-Hit column is empty, then add the gene annotation first from Bakta, if Bakta not present, add from Prokka. If both Prokka and Bakta not present,then vBen and then vCom. In all these cases, add "Screen(Maybe)" in Status. All the others Good in Status

import pandas as pd
import numpy as np

input_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/9_Result-v2.xlsx"
output_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/10_Result-v3.xlsx"

# Read the Excel file
df = pd.read_excel(input_file)

# Function to update Best-Gene and Status columns based on conditions
def update_columns(row):
    if pd.isna(row['Best-Gene']):
        if 'Bakta' in row['Sequence_Name']:
            genes = [gene.strip() for gene in str(row['Gene']).split('\\') if gene.strip() and gene.strip() != 'nan']
            if genes:
                row['Best-Gene'] = genes[0]
            row['Status'] = 'Screen(Maybe)'
        elif 'Prokka' in row['Sequence_Name']:
            genes = [gene.strip() for gene in str(row['Gene']).split('\\') if gene.strip() and gene.strip() != 'nan']
            if genes:
                row['Best-Gene'] = genes[0]
            row['Status'] = 'Screen(Maybe)'
        elif 'vBen' in row['Sequence_Name']:
            genes = [gene.strip() for gene in str(row['Gene']).split('\\') if gene.strip() and gene.strip() != 'nan']
            if genes:
                row['Best-Gene'] = genes[0]
            row['Status'] = 'Screen(Maybe)'
        elif 'vCom' in row['Sequence_Name']:
            genes = [gene.strip() for gene in str(row['Gene']).split('\\') if gene.strip() and gene.strip() != 'nan']
            if genes:
                row['Best-Gene'] = genes[0]
            row['Status'] = 'Screen(Maybe)'
    
    # If Status cell is empty, set it to "Good"
    if pd.isna(row['Status']):
        row['Status'] = 'Good'
    
    return row

# Apply the function to each row
df = df.apply(update_columns, axis=1)

# Write the updated data to a new Excel file
df.to_excel(output_file, index=False)