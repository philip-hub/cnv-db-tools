import pandas as pd; import numpy as np; import re

#gets midpoint
def proc_Pos (x):
    return (np.min(x) + np.max(x))/2

#get dataframes
kar_file_path = "kars\SJALL003310_D3_karyotype0.tsv"
kar_data = pd.read_csv(kar_file_path, sep="\t")

data_file_path = "data\SJALL003310_D0.tsv"
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
grdata = data.groupby(by=['arm', 'group_tr']).agg({'lcv': np.nanmean, 'Pos': proc_Pos, 'cv': len}).reset_index()
print(grdata)
#calulates the desired Y
grdata['Y1'] = np.log(grdata['lcv'].values / mlcv)


chromorder = dict([(c, o) for c, o in zip(['chr' + str(c) for c in range(1,23)] + ['chrX','chrY'], range(24))])

def particular_sort(series):
    return series.apply(lambda x: chromorder.get(x))

cyto = pd.read_csv ('cyto\hg19_cytoBand0.tsv', sep = '\t')

armsizes = cyto.groupby (by = 'chrom')['chromEnd'].max().reset_index()
armsizes.sort_values (by = ['chrom'], key = particular_sort, inplace = True)
armsizes['start'] = np.cumsum (np.concatenate([[0], armsizes['chromEnd'].values[:-1]]))
startDic = armsizes.set_index ('chrom')['start'].to_dict ()

grdata['X1'] = [p + startDic[a[:-1]] for p,a in zip (grdata['Pos'].tolist(), grdata['arm'].tolist())]


grdata.rename(columns={"arm": "arm1"}, inplace=True)
# Export the DataFrame to a CSV file
csv_file_path = "coverage/coverage_with_x_and_median.csv"
grdata.to_csv(csv_file_path, index=False)

# Export the DataFrame to a TSV file
tsv_file_path = "coverage/coverage_with_x_and_median.tsv"
grdata.to_csv(tsv_file_path, sep='\t', index=False)
