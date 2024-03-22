#Add a gene synonyms

import pandas as pd

input_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/10_Result-v3.xlsx"
output_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/11_Result-v4.xlsx"

# Read the Excel file
df = pd.read_excel(input_file)

def extract_gene_synonyms(row):
    gene_entry = row['Gene']
    if isinstance(gene_entry, str):  # Check if the entry is a string
        # Split the gene entry by backslashes and handle variations in formatting
        genes = [gene.strip().strip('"') for gene in gene_entry.split('\\') if gene.strip()]
        best_gene = row['Best-Gene']
        # Check if best_gene is a string or NaN
        if isinstance(best_gene, str):
            # Handle case when the best gene contains multiple gene names
            best_genes = [gene.strip().strip('"') for gene in best_gene.split('\\') if gene.strip()]
            synonyms = []
            added_genes = set()  # Keep track of genes already added to synonyms
            for gene in genes:
                if gene not in best_genes and gene != 'nan' and gene not in added_genes:
                    synonyms.append(gene)
                    added_genes.add(gene)
            return ' \\ '.join(synonyms)
        else:
            return ''  # Return an empty string for NaN values in the "Best-Gene" column
    else:
        return ''  # Return an empty string for NaN values in the "Gene" column

# Add a new column for Gene-Putative-Synonyms
df['Gene-Putative-Synonyms'] = df.apply(extract_gene_synonyms, axis=1)

# Save the modified DataFrame to a new Excel file
df.to_excel(output_file, index=False)
