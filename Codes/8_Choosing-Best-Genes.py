#Choose the best gene if it is picked by more than one tool, skip the hits that occurred once and nan
import pandas as pd

# Read the data from Excel file
df = pd.read_excel('/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/7_Grouped_Data_v3.xlsx')

# Define a function to determine the best gene
def find_best_gene(row):
    genes_str = str(row['Gene'])  # Convert NaN to string to handle missing values
    if genes_str == 'nan' or genes_str.strip() == '':  # Handle missing values
        return ''

    genes = [gene.strip('"').strip() for gene in genes_str.split(' \\ ') if gene.strip() != 'nan']
    gene_counts = {}
    for gene in genes:
        if gene:
            if gene in gene_counts:
                gene_counts[gene] += 1
            else:
                gene_counts[gene] = 1

    # Filter out genes occurring only once
    filtered_genes = [gene for gene, count in gene_counts.items() if count > 1]

    if not filtered_genes:  # If there are no genes occurring more than once
        return ''

    # Choose the best gene based on its position in the original string
    for gene in genes:
        if gene in filtered_genes:
            return gene

    return ''

# Apply the function to each row to create the new column 'Best-Gene'
df['Best-Gene'] = df.apply(find_best_gene, axis=1)

# Save the updated data to a new Excel file with the specified name
output_file_path = '/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/8_Result-v1.xlsx'
df.to_excel(output_file_path, index=False)