"""
Make a dataframe with the following

arm, cn, ai, sample

"""
#imports
import pandas as pd; import numpy as np; import re
import os


folder = "karyotypes/"
output_folder = "multi-output/"
output_file = "multi.tsv"

sample_number = 0
sample_num_in_file_name = False

def getAICN(kar_data):
    ai = kar_data['ai']
    dcn = kar_data['cn']  
    arms_kar_data = kar_data['arm']
    clone = kar_data['clone']
    
    ai_cn_data = pd.DataFrame({
        'arm':arms_kar_data,
        'cn':dcn,
        'ai':ai,
        'clone':clone
    })
    
    return ai_cn_data

combined_df = pd.DataFrame()
i = 0

for file in os.listdir(folder):
    if file.endswith(".tsv"):
        filepath = os.path.join(folder, file)
        print(f'Processing {filepath}')
        
        kar_data = pd.read_csv(filepath, sep="\t")
        ai_cn_data = getAICN(kar_data)
        
        if sample_num_in_file_name:
            sample_number = int(re.search(r'\d+', file).group())  # Assuming filenames have numbers
        else: 
            sample_number = i
        ai_cn_data['sample'] = sample_number
        combined_df = pd.concat([combined_df, ai_cn_data], ignore_index=True)
        i=i+1 #one time spent 2 hours debugging my algorithms hw bc I messed this up

print(combined_df)


os.makedirs(output_folder, exist_ok=True)
output_path = os.path.join(output_folder, output_file)
combined_df.to_csv(output_path, sep='\t', index=False)