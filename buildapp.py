import pandas as pd
import numpy as np
import re
import os
from tkinter import Tk, Label, Entry, Button, filedialog
from tkinter.ttk import Progressbar
from tqdm import tqdm

# Helper functions
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

def proc_Pos(x):
    return (np.min(x) + np.max(x)) / 2

def getCoverage(cyto_path, kar_data, data_i, cv, lcv, pos, clone, arm, group, chrom, chromEnd, outlier, deployed, arm_format, X_arm_format, Y_arm_format, x_coverage, y_coverage, chrom_start_name, rai_format, arm_coverage):
    data_filter = data_i[data_i[cv] < 20]
    data = data_filter.dropna(subset=[lcv, pos])
    ref_arms = kar_data.loc[kar_data[clone] == deployed, arm].tolist()
    tmp = data.loc[[(a in ref_arms) & (not ho) for a, ho in zip(data[arm].tolist(), data[outlier].tolist())]]
    mlcv = tmp[lcv].mean()
    grdata = data.groupby(by=[arm, group]).agg({lcv: np.nanmean, pos: proc_Pos, cv: len}).reset_index()
    grdata[y_coverage] = np.log(grdata[lcv].values / mlcv)
    chromorder = dict([(c, o) for c, o in zip([arm_format + str(c) for c in range(1, 23)] + [X_arm_format, Y_arm_format], range(24))])

    def particular_sort(series):
        return series.apply(lambda x: chromorder.get(x))

    cyto = pd.read_csv(cyto_path, sep='\t')
    armsizes = cyto.groupby(by=chrom)[chromEnd].max().reset_index()
    armsizes.sort_values(by=[chrom], key=particular_sort, inplace=True)
    armsizes[chrom_start_name] = np.cumsum(np.concatenate([[0], armsizes[chromEnd].values[:-1]]))
    startDic = armsizes.set_index(chrom)[chrom_start_name].to_dict()
    grdata[x_coverage] = [p + startDic[a[:-1]] for p, a in zip(grdata[pos].tolist(), grdata[arm].tolist())]

    if rai_format:
        columns_to_drop = ['group_tr', 'lcv', 'Pos', 'cv']
        grdata = grdata.drop(columns=columns_to_drop)

    grdata.rename(columns={arm_cname: arm_coverage}, inplace=True)
    return grdata

def getVaf(data_i, cv, outlier, arm, group, v, pos, arm_order, x, y, arm_vaf, arm_format):
    data_no_hol = data_i[data_i[outlier] != True]
    data_fil = data_no_hol[data_no_hol[cv] > 30]
    data = data_fil.dropna(subset=[v, pos])
    median_data = data.groupby([arm, group]).agg({
        v: np.median,
        pos: np.mean
    }).reset_index()

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

    median_data[arm_order] = median_data[arm].apply(extract_chromosome)
    median_data = median_data.sort_values(by=[arm_order, group]).reset_index(drop=True)
    median_data[x] = range(len(median_data))
    median_data[y] = median_data[v]

    if rai_format:
        columns_to_drop = ['group_tr', 'v', 'Pos', 'arm_order']
        median_data = median_data.drop(columns=columns_to_drop)

    median_data.rename(columns={arm_cname: arm_vaf}, inplace=True)
    return median_data

def getAICN(kar_data, ai_cname, cn_cname, arm_cname, x_CN_AI, y_CN_AI, arm_CN_AI):
    ai = kar_data[ai_cname]
    dcn = kar_data[cn_cname]
    arms_kar_data = kar_data[arm_cname]
    new_data = pd.DataFrame({
        arm_CN_AI: arms_kar_data,
        x_CN_AI: dcn,
        y_CN_AI: ai
    })
    return new_data

def process_files(cyto_folder, kars_folder, data_folder, output_folder, config_path):
    cyto_files = sorted([f for f in os.listdir(cyto_folder) if f.endswith('.tsv')])
    kars_files = sorted([f for f in os.listdir(kars_folder) if f.endswith('.tsv')])
    data_files = sorted([f for f in os.listdir(data_folder) if f.endswith('.tsv')])

    pd_cnames = pdConfig(config_path)
    rai_format = pd_cnames["rai_format"]
    cv_cname = pd_cnames["cv_cname"]
    lcv_cname = pd_cnames["lcv_cname"]
    pos_cname = pd_cnames["pos_cname"]
    clone_cname = pd_cnames["clone_cname"]
    arm_cname = pd_cnames["arm_cname"]
    group_cname = pd_cnames["group_cname"]
    chrom_cname = pd_cnames["chrom_cname"]
    chrom_end_cname = pd_cnames["chrom_end_cname"]
    outlier_cname = pd_cnames["outlier_cname"]
    ai_cname = pd_cnames["ai_cname"]
    cn_cname = pd_cnames["cn_cname"]
    deployed_name = pd_cnames["deployed"]
    arm_format = pd_cnames["arm_format"]
    X_arm_format = pd_cnames["X_arm_format"]
    Y_arm_format = pd_cnames["Y_arm_format"]
    arm_coverage = pd_cnames["arm_coverage"]
    vaf_cname = pd_cnames["vaf_cname"]
    m_cname = pd_cnames["m_cname"]
    arm_order = pd_cnames['arm_order']
    arm_vaf = pd_cnames['arm_vaf']
    x_coverage = pd_cnames["x_coverage"]
    y_coverage = pd_cnames["y_coverage"]
    x_vaf = pd_cnames["x_vaf"]
    y_vaf = pd_cnames["y_vaf"]
    chrom_start_name = pd_cnames["chrom_start"]
    ai_cname = pd_cnames["ai_cname"]
    cn_cname = pd_cnames["cn_cname"]
    arm_CN_AI = pd_cnames["arm_CN_AI"]
    x_CN_AI = pd_cnames["x_CN_AI"]
    y_CN_AI = pd_cnames["y_CN_AI"]

    for i, (cyto_file, kars_file, data_file) in tqdm(enumerate(zip(cyto_files, kars_files, data_files)), total=len(cyto_files)):
        cyto_path = os.path.join(cyto_folder, cyto_file)
        kar_file_path = os.path.join(kars_folder, kars_file)
        data_file_path = os.path.join(data_folder, data_file)

        kar_data = pd.read_csv(kar_file_path, sep='\t')
        data_i = pd.read_csv(data_file_path, sep='\t')

        coverage_df = getCoverage(cyto_path, kar_data, data_i, cv_cname, lcv_cname, pos_cname, clone_cname, arm_cname, group_cname, chrom_cname, chrom_end_cname, outlier_cname, deployed_name, arm_format, X_arm_format, Y_arm_format, x_coverage, y_coverage, chrom_start_name, rai_format, arm_coverage)
        vaf_df = getVaf(data_i, cv_cname, outlier_cname, arm_cname, group_cname, vaf_cname, pos_cname, arm_order, x_vaf, y_vaf, arm_vaf, arm_format)
        cn_ai_df = getAICN(kar_data, ai_cname, cn_cname, arm_cname, x_CN_AI, y_CN_AI, arm_CN_AI)
        final_df = pd.concat([coverage_df, vaf_df, cn_ai_df], axis=1)

        output_path = os.path.join(output_folder, f'master{i}.tsv')
        final_df.to_csv(output_path, sep='\t', index=False)

def browse_folder(entry):
    folder_selected = filedialog.askdirectory()
    entry.delete(0, 'end')
    entry.insert(0, folder_selected)

def browse_file(entry):
    file_selected = filedialog.askopenfilename(filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    entry.delete(0, 'end')
    entry.insert(0, file_selected)

def start_processing():
    cyto_folder = cyto_folder_entry.get()
    kars_folder = kars_folder_entry.get()
    data_folder = data_folder_entry.get()
    output_folder = output_folder_entry.get()
    config_path = config_file_entry.get()

    process_files(cyto_folder, kars_folder, data_folder, output_folder, config_path)

# GUI
root = Tk()
root.title("TSV Processor")

Label(root, text="Cyto Folder:").grid(row=0, column=0, sticky='e')
cyto_folder_entry = Entry(root, width=50)
cyto_folder_entry.grid(row=0, column=1)
Button(root, text="Browse", command=lambda: browse_folder(cyto_folder_entry)).grid(row=0, column=2)

Label(root, text="Kars Folder:").grid(row=1, column=0, sticky='e')
kars_folder_entry = Entry(root, width=50)
kars_folder_entry.grid(row=1, column=1)
Button(root, text="Browse", command=lambda: browse_folder(kars_folder_entry)).grid(row=1, column=2)

Label(root, text="Data Folder:").grid(row=2, column=0, sticky='e')
data_folder_entry = Entry(root, width=50)
data_folder_entry.grid(row=2, column=1)
Button(root, text="Browse", command=lambda: browse_folder(data_folder_entry)).grid(row=2, column=2)

Label(root, text="Output Folder:").grid(row=3, column=0, sticky='e')
output_folder_entry = Entry(root, width=50)
output_folder_entry.grid(row=3, column=1)
Button(root, text="Browse", command=lambda: browse_folder(output_folder_entry)).grid(row=3, column=2)

Label(root, text="Config File:").grid(row=4, column=0, sticky='e')
config_file_entry = Entry(root, width=50)
config_file_entry.grid(row=4, column=1)
Button(root, text="Browse", command=lambda: browse_file(config_file_entry)).grid(row=4, column=2)

Button(root, text="Start Processing", command=start_processing).grid(row=5, column=1, pady=10)

root.mainloop()
