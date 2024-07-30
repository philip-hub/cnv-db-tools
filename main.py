import pandas as pd
import numpy as np
import math
import re

# Function to get the value for ref
def getRef(kars, ref):
    dip_list = [ref[i] for i in range(len(kars)) if kars[i] == 'DIP']
    med_index = math.floor(len(dip_list) / 2)
    return dip_list[med_index]

#get dataframes
kar_file_path = "Data_D1_karyotype.tsv"
kar_data = pd.read_csv(kar_file_path, sep="\t")

data_file_path = "SJALL003310_D3.tsv"
data = pd.read_csv(data_file_path, sep="\t")

# lists from dataframe for kars
kars = kar_data['clone']
ref = kar_data['m']

# lists from dataframe for data
arms = data["arm"]
groups = data["group_tr"]
transcriptions = data["transcription"]
positions = data["Pos"]

# Combine relevant columns into a new DataFrame
combined_data = pd.DataFrame({
    'arm': arms,
    'group': groups,
    'transcription': pd.to_numeric(transcriptions, errors='coerce'),
    'position': pd.to_numeric(positions, errors='coerce')
})

# Drop rows where any of the relevant columns are NaN
combined_data = combined_data.dropna(subset=['transcription', 'position'])

# Group by arm and group then calculate the medians for transcription and position
median_data = combined_data.groupby(['arm', 'group']).agg({
    'transcription': 'median',
    'position': 'median'
}).reset_index()

# Calculate variable r for ref
r = getRef(kars, ref)
print(r)

# Calculate y to be formula on board being log base 2 of (transcription / r)
median_data['y'] = np.log2(median_data['transcription'] / r)

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

# Create a numeric order for arms
median_data['arm_order'] = median_data['arm'].apply(extract_chromosome)

# Ensure the numeric order is correct
median_data = median_data.sort_values(by=['arm_order', 'group']).reset_index(drop=True)

# Create a new column 'x' which is the sequential order of rows
median_data['x'] = range(len(median_data))

# Calculate median transcription for each group
group_medians = combined_data.groupby('group')['transcription'].median().reset_index()
group_medians.rename(columns={'transcription': 'median_transcription'}, inplace=True)

# Merge the median transcription back into the median_data DataFrame
median_data = median_data.merge(group_medians, on='group', how='left')

# Export the DataFrame to a CSV file
csv_file_path = "coverage_with_x_and_median.csv"
median_data.to_csv(csv_file_path, index=False)

# Export the DataFrame to a TSV file
tsv_file_path = "coverage_with_x_and_median.tsv"
median_data.to_csv(tsv_file_path, sep='\t', index=False)

print(f"Median data with x column and median transcription for each group exported to {csv_file_path} and {tsv_file_path}")
