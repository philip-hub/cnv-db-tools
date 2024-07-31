import pandas as pd
import numpy as np
import math
import re # might be needed

# # Function to get the value for ref
# def getRef(kars, ref):
#     dip_list = [ref[i] for i in range(len(kars)) if kars[i] == 'DIP']
#     med_index = math.floor(len(dip_list) / 2)
#     return dip_list[med_index]

kar_file_path = 'SJALL003310_D3_karyotype.tsv'
kar_data = pd.read_csv(kar_file_path, sep="\t")


#combined_data = data[(~data['Houtlier']) & (~data['transcription'].isna()), ['position','lcv']].copy()
 
kars = kar_data['clone']
ai = kar_data['ai']
dcn = kar_data['cn']  
arms_kar_data = kar_data['arm']



new_data = pd.DataFrame({
    "arm2":arms_kar_data,
    "X2":dcn,
    "Y2":ai
})

# Export the final DataFrame to a CSV file
csv_file_path = "karotype_graph.csv"
new_data.to_csv(csv_file_path, index=False)

print(f"Karotype graph data exported to {csv_file_path}")
