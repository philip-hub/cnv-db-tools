"""
This is the file that exports the final df for the interactive web app

### Coverage will export: X1, Y1, ARM1, M, & R

user options: filter and base log

 advanced options: change the column names of relevant columns


### vaf will export: X2, Y2, ARM2

user options: filter

 advanced options: change the column names of relevant columns


### CN_AI will export: X23, Y3, ARM3

user options: filter

advanced options: change the column names of relevant columns


"""
#imports
import pandas as pd; import numpy as np; import re

#helper functions
def pdConfig(file_path):
    variables = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.split('#')[0].strip()
            if not line:
                continue
            if '=' in line:
                var, value = line.split('=', 1)
                var = var.strip()
                value = value.strip().strip("'")
                variables[var] = value

    return variables


def proc_Pos (x):
    return (np.min(x) + np.max(x))/2




# Example usage
pd_config_path = 'config/pd.txt'
pd_cnames = pdConfig(pd_config_path)
#print(pd_cnames)

cv_cname = pd_cnames["cv_cname"]
lcv_cname = pd_cnames["lcv_cname"]
pos_cname = pd_cnames["pos_cname"]
clone_cname = pd_cnames["clone_cname"]
arm_cname = pd_cnames["arm_cname"]
group_cname = pd_cnames["group_cname"]
chrom_cname = pd_cnames["chrom_cname"]
chrom_end_cname = pd_cnames["chrom_end_cname"]
outlier_cname = pd_cnames["outlier_cname"]


#coverage 
deployed_name = pd_cnames["deployed"]
arm_format = pd_cnames["arm_format"]
X_arm_format = pd_cnames["X_arm_format"]
Y_arm_format = pd_cnames["Y_arm_format"]


#danger zone options
x_coverage = pd_cnames["x_coverage"]
y_coverage = pd_cnames["y_coverage"]
chrom_start_name = pd_cnames["chrom_start"]

#file paths
cyto_path = 'inputs/hg19_cytoBand.tsv'
kar_file_path = "inputs/Data_D1_karyotype.tsv"
data_file_path = "inputs/SJALL003310_D3.tsv"

#get dataframes
kar_data = pd.read_csv(kar_file_path, sep="\t")
data_i = pd.read_csv(data_file_path, sep="\t")





def getCoverage(cyto_path, data_i,cv,lcv,pos,clone,arm, group,chrom,chromEnd,outlier, deployed,arm_format,X_arm_format,Y_arm_format, x_coverage, y_coverage, chrom_start_name):
    #filter
    data_filter = data_i[data_i[cv] < 20]
    print("After Filtering CV: "+str(data_filter.shape))
    data = data_filter.dropna(subset=[lcv, pos])

    ref_arms = kar_data.loc[kar_data[clone] == deployed, arm].tolist()
    tmp = data.loc[[(a in ref_arms) & (not ho) for a, ho in zip(data[arm].tolist(), data[outlier].tolist())]]
    # what I called r
    mlcv = tmp[lcv].mean()
    print("mlcv : "+str(mlcv))
    #groups the data
    grdata = tmp.groupby(by=[arm, group]).agg({lcv: np.nanmean, pos: proc_Pos, cv: len}).reset_index()
    print(grdata)
    #calulates the desired Y
    grdata[y_coverage] = np.log(grdata[lcv].values / mlcv)


    chromorder = dict([(c, o) for c, o in zip([arm_format + str(c) for c in range(1,23)] + [X_arm_format,Y_arm_format], range(24))])

    def particular_sort(series):
        return series.apply(lambda x: chromorder.get(x))

    cyto = pd.read_csv (cyto_path, sep = '\t')

    armsizes = cyto.groupby (by = chrom)[chromEnd].max().reset_index()
    armsizes.sort_values (by = [chrom], key = particular_sort, inplace = True)
    armsizes[chrom_start_name] = np.cumsum (np.concatenate([[0], armsizes[chromEnd].values[:-1]]))
    startDic = armsizes.set_index (chrom)[chrom_start_name].to_dict ()

    grdata[x_coverage] = [p + startDic[a[:-1]] for p,a in zip (grdata[pos].tolist(), grdata[arm].tolist())]

    csv_file_path = "coverage/coverage_with_x_and_median.csv"
    grdata.to_csv(csv_file_path, index=False)

    
    return grdata





getCoverage(cyto_path, data_i,cv_cname,lcv_cname,pos_cname,clone_cname,arm_cname, group_cname,chrom_cname,chrom_end_cname,outlier_cname, deployed_name,arm_format,X_arm_format,Y_arm_format, x_coverage, y_coverage, chrom_start_name)