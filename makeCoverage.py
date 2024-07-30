import pandas as pd
import numpy as np
import math
import re

kar_file_path = "Data_D1_karyotype.tsv"
kar_data = pd.read_csv(kar_file_path, sep="\t")

data_file_path = "SJALL003310_D3.tsv"
data = pd.read_csv(data_file_path, sep="\t")

data = data[data['Houtlier'] != True]

kars = kar_data['clone']
ref = kar_data['m']

arms = data["arm"]
lcvs = data["lcv"]
positions = data["Pos"]

combined_data = pd.DataFrame({
    'arm': arms,
    'lcv': pd.to_numeric(lcvs, errors='coerce'),
    'position': pd.to_numeric(positions, errors='coerce')
})

#drop na
combined_data = combined_data.dropna(subset=['lcv', 'position'])

r = combined_data.loc[combined_data['arm'].isin(ref), 'lcv'].median()

ref_arm = kar_data.loc[kar_data['clone'] == 'DIP', 'arm'].tolist()
print(ref_arm)

r = combined_data.loc[[a in ref_arm for a in combined_data['arm'].tolist()], 'lcv'].median()
print(r)

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
                return number * 2 + (1 if arm_type == 'q' else 0)
    return -1  # Return -1 if not matched

combined_data['arm_order'] = combined_data['arm'].apply(extract_chromosome)

combined_data = combined_data.sort_values(by=['arm_order']).reset_index(drop=True)

combined_data['x'] = range(len(combined_data))

combined_data['y'] = np.log2(combined_data['lcv'] / r)

ref_y = combined_data['y'].tolist()


with open('y_column.txt', 'w') as file:
    for entry in ref_y:
        file.write(str(entry) + '\n')

csv_file_path = "coverage_with_x_and_median.csv"
combined_data.to_csv(csv_file_path, index=False)

tsv_file_path = "coverage_with_x_and_median.tsv"
combined_data.to_csv(tsv_file_path, sep='\t', index=False)

print(f"Data with x column and log2 transformed lcv exported to {csv_file_path} and {tsv_file_path}")
