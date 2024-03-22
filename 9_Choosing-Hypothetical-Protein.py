# If Type CDS- For rows with all nan or empty consider those as hypothetical protein

import pandas as pd

input_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/8_Result-v1.xlsx"
output_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/9_Result-v2.xlsx"

# Read the Excel file
df = pd.read_excel(input_file)

# Iterate over rows
for index, row in df.iterrows():
    # Check if the value in the "Type" column is "CDS"
    if row["Type"] == "CDS":
        # Check if the value in the "Gene" column is not NaN and not empty
        if not pd.isna(row["Gene"]) and row["Gene"]:
            # Check if Gene column contains only "nan"
            if all(gene == "nan" for gene in row["Gene"].split(" \ ")):
                # Update Best-Gene column
                df.at[index, "Best-Gene"] = "HP"
        # If the value in the "Gene" column is NaN or empty
        else:
            # Update Best-Gene column
            df.at[index, "Best-Gene"] = "HP"

# Write back to Excel
df.to_excel(output_file, index=False)