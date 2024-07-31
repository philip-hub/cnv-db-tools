import pandas as pd; import numpy as np; import re

#gets midpoint
def proc_Pos (x):
    return (np.min(x) + np.max(x))/2

#get dataframes
kar_file_path = "inputs/Data_D1_karyotype.tsv"
kar_data = pd.read_csv(kar_file_path, sep="\t")

data_file_path = "inputs/SJALL003310_D3.tsv"
data_i = pd.read_csv(data_file_path, sep="\t")

#filter
data_filter = data_i[data_i['cv'] < 20]
print("After Filtering CV: "+str(data_filter.shape))
data = data_filter.dropna(subset=['lcv', 'Pos'])

ref_arms = kar_data.loc[kar_data['clone'] == 'DIP', 'arm'].tolist()
tmp = data.loc[[(a in ref_arms) & (not ho) for a, ho in zip(data['arm'].tolist(), data['Houtlier'].tolist())]]
# what I called r
mlcv = tmp['lcv'].mean()
print("mlcv : "+str(mlcv))
#groups the data
grdata = tmp.groupby(by=['arm', 'group_tr']).agg({'lcv': np.nanmean, 'Pos': proc_Pos, 'cv': len}).reset_index()
print(grdata)
#calulates the desired Y
grdata['y'] = np.log(grdata['lcv'].values / mlcv)


chromorder = dict([(c, o) for c, o in zip(['chr' + str(c) for c in range(1,23)] + ['chrX','chrY'], range(24))])

def particular_sort(series):
    return series.apply(lambda x: chromorder.get(x))

cyto = pd.read_csv ('inputs/hg19_cytoBand.tsv', sep = '\t')

armsizes = cyto.groupby (by = 'chrom')['chromEnd'].max().reset_index()
armsizes.sort_values (by = ['chrom'], key = particular_sort, inplace = True)
armsizes['start'] = np.cumsum (np.concatenate([[0], armsizes['chromEnd'].values[:-1]]))
startDic = armsizes.set_index ('chrom')['start'].to_dict ()

grdata['x'] = [p + startDic[a[:-1]] for p,a in zip (grdata['Pos'].tolist(), grdata['arm'].tolist())]



# # Create an ordered mapping for arms
# def extract_chromosome(arm):
#     if arm.startswith('chr'):
#         if arm == 'chrXp':
#             return 1000
#         elif arm == 'chrXq':
#             return 1001
#         elif arm == 'chrYp':
#             return 1002
#         elif arm == 'chrYq':
#             return 1003
#         else:
#             match = re.match(r'chr(\d+)([pq])', arm)
#             if match:
#                 number = int(match.group(1))
#                 arm_type = match.group(2)
#                 # Ensure numeric chromosomes are ordered first
#                 return number * 2 + (1 if arm_type == 'q' else 0)
#     return -1  # Return -1 if not matched

# # Create a numeric order for arms
# grdata['arm_order'] = grdata['arm'].apply(extract_chromosome)
# grdata = grdata.sort_values(by=['arm_order', 'group_tr']).reset_index(drop=True)
# grdata['x'] = range(len(grdata))

# Export the DataFrame to a CSV file
csv_file_path = "coverage/coverage_with_x_and_median.csv"
grdata.to_csv(csv_file_path, index=False)

# Export the DataFrame to a TSV file
tsv_file_path = "coverage/coverage_with_x_and_median.tsv"
grdata.to_csv(tsv_file_path, sep='\t', index=False)
