#Grouping 2 - If either the Start or the End is same and the Sequence ID and Type is same. We choose the Start and End obtained from most tools, if equal then preference is alphabetically to Sequence Name (Bakta, Prokka, vBen and then vCom). Add Sequence Name, Gene and Product. # indicate same Start and ## indicate same End

import pandas as pd

# Reading the input Excel file
df = pd.read_excel("/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/5_Grouped_Data_v1.xlsx")

# Sorting the DataFrame by specific columns
df_sorted = df.sort_values(by=['Sequence_ID', 'Type', 'Start', 'End'])

# Initialize lists to store merged rows and their indices
merged_rows = []
merged_indices = []

# Iterating through each row in the sorted DataFrame
for i, row in df_sorted.iterrows():
    if i in merged_indices:
        continue  # Skip already merged rows
        
    row1 = row
    
    # Extracting data from the current row
    merged_gene = str(row1['Gene']) 
    merged_product = str(row1['Product']) 
    merged_sequence_name = row1['Sequence_Name']
    merged_sequence_indicator = ()
    merged_start = row1['Start']
    merged_end = row1['End']

    # Iterating through subsequent rows to check for merging conditions
    for j, row2 in df_sorted.loc[i + 1:].iterrows():
        if (row1['Start'] == row2['Start'] and row1['End'] != row2['End']) or \
           (row1['End'] == row2['End'] and row1['Start'] != row2['Start']):
            # Handling cases where Start or End matches
            if row1['Start'] == row2['Start']:
                merged_sequence_indicator = tuple(list(merged_sequence_indicator) + ['#'])  # Indicates same Start
            elif row1['End'] == row2['End']:
                merged_sequence_indicator = tuple(list(merged_sequence_indicator) + ['##'])  # Indicates same End

            # Determining which row's data to prioritize for merging
            num_backslashes_row1 = row1['Sequence_Name'].count('\\')
            num_backslashes_row2 = row2['Sequence_Name'].count('\\')
            
            if num_backslashes_row1 > num_backslashes_row2 or \
               (num_backslashes_row1 == num_backslashes_row2 and 'Bakta' in row1['Sequence_Name']):
                merged_start = row1['Start']
                merged_end = row1['End']
                merged_gene = str(row1['Gene']) + ' \\ ' + str(row2['Gene'])
                merged_product = str(row1['Product']) + ' \\ ' + str(row2['Product'])
                merged_sequence_name = row1['Sequence_Name'] + merged_sequence_indicator[0] + ' \\ ' + row2['Sequence_Name']
                merged_indices.append(i)
            else:
                merged_start = row2['Start']
                merged_end = row2['End']
                merged_gene = str(row2['Gene']) + ' \\ ' + str(row1['Gene'])
                merged_product = str(row2['Product']) + ' \\ ' + str(row1['Product'])
                merged_sequence_name = row2['Sequence_Name'] + merged_sequence_indicator[0] + ' \\ ' + row1['Sequence_Name']
                merged_indices.append(j)
                
            # Add the index of the row to be deleted to the merged_indices list
            merged_indices.append(j)  
        else:
            break  # Exit loop if conditions are not met

    # Append merged row data to the list
    merged_rows.append({
        'Sequence_ID': row1['Sequence_ID'],
        'Type': row1['Type'],
        'Start': merged_start,
        'End': merged_end,
        'Length': row1['Length'],
        'Direction': row1['Direction'],
        'Sequence_Name': merged_sequence_name,
        'Gene': merged_gene,
        'Product': merged_product
    })

# Remove duplicates from merged_indices list
merged_indices = list(set(merged_indices))

# Remove rows with indices in merged_indices list from the sorted DataFrame
df_sorted.drop(merged_indices, inplace=True)

# Create DataFrame without deleted rows
remaining_df = df_sorted

# Create DataFrame for merged rows
merged_df = pd.DataFrame(merged_rows)

# Writing the merged DataFrame to an Excel file
merged_df.to_excel("/Users/simrankushwaha/Documents/PostDoc/Annotations/474/SKK-Pipeline/6_Grouped_Data_v2.xlsx", index=False)