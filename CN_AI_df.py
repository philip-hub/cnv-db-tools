import pandas as pd
import numpy as np
import math
import re # might be needed

kar_file_path = 'inputs/SJALL003310_D3_karyotype.tsv'
kar_data = pd.read_csv(kar_file_path, sep="\t")

ai = kar_data['ai']
dcn = kar_data['cn']  
arms_kar_data = kar_data['arm']



new_data = pd.DataFrame({
    "arm2":arms_kar_data,
    "X2":dcn,
    "Y2":ai
})

# Export the final DataFrame to a CSV file
csv_file_path = "CN_AI/karotype_graph.csv"
new_data.to_csv(csv_file_path, index=False)

print(f"Karotype graph data exported to {csv_file_path}")
