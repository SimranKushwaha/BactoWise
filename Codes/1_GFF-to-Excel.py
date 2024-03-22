#Create Excel file from GFF file
import os
import pandas as pd

def extract_sequence_name(gff_file):
    filename = os.path.splitext(os.path.basename(gff_file))[0]
    with open(gff_file) as gff_handle:
        for line in gff_handle:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) >= 9:
                seq_id = fields[0]
                type_ = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                length = end - start + 1
                direction = fields[6]
                attributes = dict(item.split("=") for item in fields[8].split(";"))
                gene = attributes.get('gene', '')  # Extract gene name
                product = attributes.get('product', '')  # Extract product

                yield filename, seq_id, type_, start, end, length, direction, gene, product

# Example usage
gff_file = "/Users/simrankushwaha/Documents/PostDoc/Annotations/474/GFF+Excel/Bakta.gff"
filename = os.path.splitext(os.path.basename(gff_file))[0]  # Extract filename without extension

data = []
for record in extract_sequence_name(gff_file):
    data.append(record)

# Convert the data to a DataFrame
df = pd.DataFrame(data, columns=["Sequence_Name", "Sequence_ID", "Type", "Start", "End", "Length", "Direction", "Gene", "Product"])

# Write DataFrame to Excel file with the same name as the input GFF file
output_excel = os.path.join(os.path.dirname(gff_file), filename + ".xlsx")
df.to_excel(output_excel, index=False)