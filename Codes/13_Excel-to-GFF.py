#Create a GFF file for the Excel sheet. 

import pandas as pd

# Read the input Excel file
input_file_path = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/12_Result-v5.xlsx"
df = pd.read_excel(input_file_path)

# Create the new columns with specified data
df['New_Column_2'] = 'SKKPipeline'
df['New_Column_6'] = '.'
df['New_Column_8'] = '.'

# Function to create the concatenation
def create_concat(row):
    new_id = row['New_ID']
    new_type = row['Type']
    product_parts = row['Product'].split('\\')
    product = None
    for part in product_parts:
        if part.strip().lower() != 'hypothetical protein':
            product = part.strip()
            break
    if product is None:
        # If all parts are 'hypothetical protein', choose the last one
        product = product_parts[-1].strip()
    if pd.isna(row['Best-Gene']) or row['Best-Gene'] == 'HP':
        concat_str = f"ID={new_id};gene_biotype={new_type};note={product}"
    else:
        best_gene = row['Best-Gene']
        concat_str = f"ID={new_id};gene_biotype={new_type};Name={best_gene};note={product}"
    return concat_str

# Apply the function to create the new column
df['Concatenated_Column'] = df.apply(create_concat, axis=1)

# Select the required columns
output_df = df[['Sequence_ID', 'New_Column_2', 'Type', 'Start', 'End', 'New_Column_6', 'Direction', 'New_Column_8', 'Concatenated_Column']]

# Write the output to a tab-separated text file with ".gff" extension
output_file_path = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/New-Annotated-SKK.gff"
output_df.to_csv(output_file_path, index=False, sep='\t', header=False)