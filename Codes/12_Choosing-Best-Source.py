# Adding the column sequence_sources to tell of which anotation tool was finally considered for the analysis

import pandas as pd

# Read the Excel file
file_path = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/11_Result-v4.xlsx"
df = pd.read_excel(file_path)

# Create a new column 'Source-Chosen' and initialize it with empty strings
df['Source-Chosen'] = ''

# Iterate over each row
for index, row in df.iterrows():
    # Extract the best gene from the 'Best-Gene' column
    best_gene = row['Best-Gene']
    
    # Check if best_gene is NaN, if so, skip further processing
    if pd.isnull(best_gene):
        continue
    
    # Check if the best gene is 'HP', if so, copy all sequence names to 'Source-Chosen'
    if str(best_gene).strip() == 'HP':
        df.at[index, 'Source-Chosen'] = ', '.join(str(row['Sequence_Name']).split(' \\ ')).replace('#', '').replace('##', '')
        continue
    
    # Extract 'Sequence_Name' and 'Gene' columns, handling NaN values
    sequence_names = str(row['Sequence_Name']).split(' \\ ') if not pd.isnull(row['Sequence_Name']) else []
    genes = str(row['Gene']).split(' \\ ') if not pd.isnull(row['Gene']) else []
    
    # Find the sequence name(s) corresponding to the best gene
    for sequence, gene in zip(sequence_names, genes):
        # Check if the gene is valid and matches the best gene
        if str(gene).strip() == str(best_gene).strip():
            # Append the corresponding sequence name to the 'Source-Chosen' column
            df.at[index, 'Source-Chosen'] += sequence + ', '

# Remove the trailing comma and space from the 'Source-Chosen' column and drop '#' or '##'
df['Source-Chosen'] = df['Source-Chosen'].str.replace(r'#\s*|\#{2}\s*', '', regex=True).str.rstrip(', ')

# Save the updated DataFrame to a new Excel file
output_file_path = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/12_Result-v5.xlsx"
df.to_excel(output_file_path, index=False)