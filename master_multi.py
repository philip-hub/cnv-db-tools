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
arm_format = 'chr'

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
                pattern = rf'{re.escape(arm_format)}(\d+)([pq])'
                match = re.match(pattern, arm)
                if match:
                    number = int(match.group(1))
                    arm_type = match.group(2)
                    return number * 2 + (1 if arm_type == 'q' else 0)
        return -1
    


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
    
    ai_cn_data['arm_order'] = ai_cn_data['arm'].apply(extract_chromosome)
    ai_cn_data = ai_cn_data.sort_values(by=['arm_order']).reset_index(drop=True)
    
    return ai_cn_data

combined_df = pd.DataFrame()
i = 0

# for j in range(0,8): #to test scaling and stuff
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
                ai_cn_data['filename'] = file 
                combined_df = pd.concat([combined_df, ai_cn_data], ignore_index=True)
                i=i+1 #one time spent 2 hours debugging my algorithms hw bc I messed this up




    
print(combined_df)


os.makedirs(output_folder, exist_ok=True)
output_path = os.path.join(output_folder, output_file)
combined_df.to_csv(output_path, sep='\t', index=False)