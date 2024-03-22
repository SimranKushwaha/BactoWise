#Grouping 1 - If Sequence ID, Type, Start, End, Length, Direction is same; join Sequence Name, Gene and Product
import pandas as pd

# Read the Excel file
df = pd.read_excel("/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/4_Modified_Genes_Files.xlsx")

# Group by relevant columns and aggregate Sequence Name, Gene, and Product
grouped_df = df.groupby(['Sequence_ID', 'Type', 'Start', 'End', 'Length', 'Direction']).agg({
    'Sequence_Name': lambda x: ' \\ '.join(x.astype(str)),
    'Gene': lambda x: ' \\ '.join(x.astype(str)),
    'Product': lambda x: ' \\ '.join(x.astype(str))
}).reset_index()

# Sort by Start and Sequence_Name
grouped_df = grouped_df.sort_values(by=['Sequence_ID', 'Start'])

# Write the grouped data to a new Excel file
grouped_df.to_excel("/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/5_Grouped_Data_v1.xlsx", index=False)