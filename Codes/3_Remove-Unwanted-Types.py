# Remove everything except CDS	CRISPR	ncRNA	oriC	oriT	oriV	regulatory_region	rRNA	tmRNA	tRNA

import pandas as pd

# Define input and output file paths
input_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/2_Merged_Files.xlsx"
output_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/3_Final_Genes_Files.xlsx"

# Read the Excel file into a pandas DataFrame
df = pd.read_excel(input_file)

# Define a list of valid types
valid_types = ["CDS", "CRISPR", "ncRNA", "oriC", "oriT", "oriV", "regulatory_region", "rRNA", "tmRNA", "tRNA"]

# Filter rows based on the 'Type' column
filtered_df = df[df['Type'].isin(valid_types)]

# Write the filtered DataFrame to a new Excel file
filtered_df.to_excel(output_file, index=False)