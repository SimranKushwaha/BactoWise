# Remove the "_number" after the gene names and Make them in xxxY form
import pandas as pd

# Read the Excel file
file_path = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/3_Final_Genes_Files.xlsx"
df = pd.read_excel(file_path)

# Define a function to remove "_number" pattern from gene names and format gene names
def format_gene_name(gene_name):
    import re
    if pd.isna(gene_name):
        return ''
    gene_name = re.sub(r'(_\d+)?$', '', str(gene_name))
    if len(gene_name) == 4 and gene_name.isalpha():
        gene_name = gene_name[:3].lower() + gene_name[3].upper()
    return gene_name

# Apply the function to the "Gene" column
df['Gene'] = df['Gene'].apply(format_gene_name)

# Save the modified data back to an Excel file
modified_file_path = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/4_Modified_Genes_Files.xlsx"
df.to_excel(modified_file_path, index=False, na_rep='')