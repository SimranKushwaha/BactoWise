#Here we add new ID numbers- So the format of the ID number is Sequence_ID-Type-5 Digit Number. The Number will be in the ascending order of the Type in the genome location.
import pandas as pd

# Load the data from the Excel file
data = pd.read_excel('/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/12_Result-v5.xlsx')

# Sort the data first by Sequence_ID and then by Start
data_sorted = data.sort_values(by=['Sequence_ID', 'Start'])

# Function to generate the ID based on Sequence_ID, Type, and the order
def generate_id(row, type_counts):
    sequence_type = (row['Sequence_ID'], row['Type'])
    if sequence_type not in type_counts:
        type_counts[sequence_type] = 1
    else:
        type_counts[sequence_type] += 1
    id_str = f"{row['Sequence_ID']}_{row['Type']}_{type_counts[sequence_type]:05d}"
    return id_str

# Initialize a dictionary to keep track of the counts for each (Sequence ID, Type) combination
type_counts = {}

# Generate IDs for each row
data_sorted['New_ID'] = data_sorted.apply(lambda row: generate_id(row, type_counts), axis=1)

# Save the result to a new Excel file
data_sorted.to_excel('/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/13_Result-v6.xlsx', index=False)