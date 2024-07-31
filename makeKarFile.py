import pandas as pd
import numpy as np
import math
import re # might be needed

# # Function to get the value for ref
# def getRef(kars, ref):
#     dip_list = [ref[i] for i in range(len(kars)) if kars[i] == 'DIP']
#     med_index = math.floor(len(dip_list) / 2)
#     return dip_list[med_index]

kar_file_path = "Data_D1_karyotype.tsv"
kar_data = pd.read_csv(kar_file_path, sep="\t")
data_file_path = "SJALL003310_D3.tsv"
data_i = pd.read_csv(data_file_path, sep="\t")

data = data_i[data_i['Houtlier'] != True]

#combined_data = data[(~data['Houtlier']) & (~data['transcription'].isna()), ['position','lcv']].copy()
 
kars = kar_data['clone']
ref = kar_data['m']
dcn = kar_data['cn']  
arms_kar_data = kar_data['arm']

arms = data["arm"]
transcriptions = data["transcription"]


combined_data = pd.DataFrame({
    'arm': arms,
    'transcription': pd.to_numeric(transcriptions, errors='coerce'),
})

# drop NaN
combined_data = combined_data.dropna(subset=['transcription'])

# median transcription for each arm calculations
arm_medians = combined_data.groupby('arm')['transcription'].median().reset_index()
arm_medians.rename(columns={'transcription': 'arm_median_transcription'}, inplace=True)

ref_arm = kar_data.loc[kar_data['clone'] == 'DIP', 'arm'].tolist()
print(ref_arm)

r = data.loc[[a in ref_arm for a in data['arm'].tolist()], 'lcv'].median()
print(r)


#ref 
#r = getRef(kars, ref)

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
                return number * 2 + (1 if arm_type == 'q' else 0)
    return -1

# Create a numeric order for arms
arm_medians['arm_order'] = arm_medians['arm'].apply(extract_chromosome)

# Ensure the numeric order is correct
arm_medians = arm_medians.sort_values(by='arm_order').reset_index(drop=True)

# Create a new column 'x' which is the sequential order of rows
arm_medians['x'] = range(len(arm_medians))

# Add the column 'dcn' to the arm_medians DataFrame
# Create a dictionary for 'dcn' values from karyotype data
dcn_dict = kar_data.set_index('arm')['dcn'].to_dict()
arm_medians['dcn'] = arm_medians['arm'].map(dcn_dict)

# Calculate column X as (median transcription * 2) / r
arm_medians['X'] = (arm_medians['arm_median_transcription'] * 2) / r

# Prepare final DataFrame with only the required columns
final_df = arm_medians[['arm', 'arm_median_transcription', 'X', 'dcn']]

# Export the final DataFrame to a CSV file
csv_file_path = "karotype_graph.csv"
final_df.to_csv(csv_file_path, index=False)

print(f"Karotype graph data exported to {csv_file_path}")
