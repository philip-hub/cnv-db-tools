import pandas as pd
import numpy as np
import math
import re

#get dataframes
kar_file_path = "inputs/Data_D1_karyotype.tsv"
kar_data = pd.read_csv(kar_file_path, sep="\t")


data_file_path = "inputs/SJALL003310_D3.tsv"
data_i = pd.read_csv(data_file_path, sep="\t")



data_no_hol = data_i[data_i['Houtlier'] != True]
data = data_no_hol[data_no_hol['cv'] > 30]

# lists from dataframe for kars
kars = kar_data['clone']
ref = kar_data['m']

# lists from dataframe for data
arms = data["arm"]
groups = data["group_tr"]
vaf = data["v"]
positions = data["Pos"]

# Combine relevant columns into a new DataFrame
combined_data = pd.DataFrame({
    'arm': arms,
    'group': groups,
    'v': pd.to_numeric(vaf, errors='coerce'),
    'position': pd.to_numeric(positions, errors='coerce')
})

# Drop rows where any of the relevant columns are NaN
combined_data = combined_data.dropna(subset=['v', 'position'])

r = combined_data.loc[combined_data['arm'].isin(ref), 'v'].median()
     

median_data = combined_data.groupby(['arm', 'group']).agg({
    'v': np.median,
    'position': np.mean
}).reset_index()

# Create an ordered mapping for arms
def extract_chromosome(arm):
    if arm.startswith('chr'):
        if arm == 'chrXp':
            return 1000
        elif arm == 'chrXq':
            return 1001
        elif arm == 'chrYp':
            return 1002
        elif arm == 'chrYq':
            return 1003
        else:
            match = re.match(r'chr(\d+)([pq])', arm)
            if match:
                number = int(match.group(1))
                arm_type = match.group(2)
                # Ensure numeric chromosomes are ordered first
                return number * 2 + (1 if arm_type == 'q' else 0)
    return -1  # Return -1 if not matched

# Create a numeric order for arms & numeric order is correct
median_data['arm_order'] = median_data['arm'].apply(extract_chromosome)
median_data = median_data.sort_values(by=['arm_order', 'group']).reset_index(drop=True)

# Create a new column 'x' which is the sequential order of rows
median_data['x'] = range(len(median_data))

# Calculate median v for each group
group_medians = combined_data.groupby('group')['v'].median().reset_index()
group_medians.rename(columns={'v': 'median_v'}, inplace=True)


median_data = median_data.merge(group_medians, on='group', how='left')
median_data['y'] = median_data['v']

ref_y = median_data['y'].tolist()




with open('y_column.txt', 'w') as file:
    for entry in ref_y:
        file.write(str(entry) + '\n')

final_data = median_data[median_data['y'] < .98]

# Export the DataFrame to a CSV file
csv_file_path = "vaf/vaf_coverage_with_x_and_median.csv"

print(final_data.shape)
final_data.to_csv(csv_file_path, index=False)

# Export the DataFrame to a TSV file
tsv_file_path = "vaf/vaf_coverage_with_x_and_median.tsv"
final_data.to_csv(tsv_file_path, sep='\t', index=False)

print(f"Median data with x column and median v for each group exported to {csv_file_path} and {tsv_file_path}")
