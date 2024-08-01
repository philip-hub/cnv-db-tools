import pandas as pd
import numpy as np
import re
import os
import webbrowser
from tkinter import Tk, Label, Entry, Button, filedialog, StringVar, font
from tkinter.ttk import Progressbar
from tqdm import tqdm

# Default configuration settings
default_config = {
    "cv_cname": 'cv',
    "lcv_cname": 'lcv',
    "pos_cname": 'Pos',
    "clone_cname": 'clone',
    "arm_cname": 'arm',
    "group_cname": 'group_tr',
    "outlier_cname": 'Houtlier',
    "chrom_cname": 'chrom',
    "chrom_end_cname": 'chromEnd',
    "ai_cname": 'ai',
    "cn_cname": 'cn',
    "m_cname": 'm',
    "vaf_cname": 'v',
    "arm_order": 'arm_order',
    "deployed": 'DIP',
    "arm_format": 'chr',
    "X_arm_format": 'chrX',
    "Y_arm_format": 'chrY',
    "chrom_start": 'Start',
    "x_coverage": 'X1',
    "y_coverage": 'Y1',
    "arm_coverage": 'arm1',
    "x_vaf": 'X2',
    "y_vaf": 'Y2',
    "arm_vaf": 'arm2',
    "x_CN_AI": 'X3',
    "y_CN_AI": 'Y3',
    "arm_CN_AI": 'arm3',
    "rai_format": True
}

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

def getCoverage(cyto_path, kar_data, data_i, config):
    data_filter = data_i[data_i[config['cv_cname']] < 20]
    data = data_filter.dropna(subset=[config['lcv_cname'], config['pos_cname']])
    ref_arms = kar_data.loc[kar_data[config['clone_cname']] == config['deployed'], config['arm_cname']].tolist()
    tmp = data.loc[[(a in ref_arms) & (not ho) for a, ho in zip(data[config['arm_cname']].tolist(), data[config['outlier_cname']].tolist())]]
    mlcv = tmp[config['lcv_cname']].mean()
    grdata = data.groupby(by=[config['arm_cname'], config['group_cname']]).agg({config['lcv_cname']: np.nanmean, config['pos_cname']: proc_Pos, config['cv_cname']: len}).reset_index()
    grdata[config['y_coverage']] = np.log(grdata[config['lcv_cname']].values / mlcv)
    chromorder = dict([(c, o) for c, o in zip([config['arm_format'] + str(c) for c in range(1, 23)] + [config['X_arm_format'], config['Y_arm_format']], range(24))])

    def particular_sort(series):
        return series.apply(lambda x: chromorder.get(x))

    cyto = pd.read_csv(cyto_path, sep='\t')
    armsizes = cyto.groupby(by=config['chrom_cname'])[config['chrom_end_cname']].max().reset_index()
    armsizes.sort_values(by=[config['chrom_cname']], key=particular_sort, inplace=True)
    armsizes[config['chrom_start']] = np.cumsum(np.concatenate([[0], armsizes[config['chrom_end_cname']].values[:-1]]))
    startDic = armsizes.set_index(config['chrom_cname'])[config['chrom_start']].to_dict()
    grdata[config['x_coverage']] = [p + startDic[a[:-1]] for p, a in zip(grdata[config['pos_cname']].tolist(), grdata[config['arm_cname']].tolist())]

    if config['rai_format']:
        columns_to_drop = ['group_tr', 'lcv', 'Pos', 'cv']
        grdata = grdata.drop(columns=columns_to_drop)

    grdata.rename(columns={config['arm_cname']: config['arm_coverage']}, inplace=True)
    return grdata

def getVaf(data_i, config):
    data_no_hol = data_i[data_i[config['outlier_cname']] != True]
    data_fil = data_no_hol[data_no_hol[config['cv_cname']] > 30]
    data = data_fil.dropna(subset=[config['vaf_cname'], config['pos_cname']])
    median_data = data.groupby([config['arm_cname'], config['group_cname']]).agg({
        config['vaf_cname']: np.median,
        config['pos_cname']: np.mean
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
                pattern = rf"{re.escape(config['arm_format'])}(\d+)([pq])"
                match = re.match(pattern, arm)
                if match:
                    number = int(match.group(1))
                    arm_type = match.group(2)
                    return number * 2 + (1 if arm_type == 'q' else 0)
        return -1

    median_data[config['arm_order']] = median_data[config['arm_cname']].apply(extract_chromosome)
    median_data = median_data.sort_values(by=[config['arm_order'], config['group_cname']]).reset_index(drop=True)
    median_data[config['x_vaf']] = range(len(median_data))
    median_data[config['y_vaf']] = median_data[config['vaf_cname']]

    if config['rai_format']:
        columns_to_drop = ['group_tr', 'v', 'Pos', 'arm_order']
        median_data = median_data.drop(columns=columns_to_drop)

    median_data.rename(columns={config['arm_cname']: config['arm_vaf']}, inplace=True)
    return median_data

def getAICN(kar_data, config):
    ai = kar_data[config['ai_cname']]
    dcn = kar_data[config['cn_cname']]
    arms_kar_data = kar_data[config['arm_cname']]
    new_data = pd.DataFrame({
        config['arm_CN_AI']: arms_kar_data,
        config['x_CN_AI']: dcn,
        config['y_CN_AI']: ai
    })
    return new_data

def process_files(cyto_folder, kars_folder, data_folder, output_folder, config, progress_var, max_progress):
    cyto_files = sorted([f for f in os.listdir(cyto_folder) if f.endswith('.tsv')])
    kars_files = sorted([f for f in os.listdir(kars_folder) if f.endswith('.tsv')])
    data_files = sorted([f for f in os.listdir(data_folder) if f.endswith('.tsv')])

    progress_var.set(0)
    max_progress.set(len(cyto_files))

    for i, (cyto_file, kars_file, data_file) in tqdm(enumerate(zip(cyto_files, kars_files, data_files)), total=len(cyto_files)):
        cyto_path = os.path.join(cyto_folder, cyto_file)
        kar_file_path = os.path.join(kars_folder, kars_file)
        data_file_path = os.path.join(data_folder, data_file)

        kar_data = pd.read_csv(kar_file_path, sep='\t')
        data_i = pd.read_csv(data_file_path, sep='\t')

        coverage_df = getCoverage(cyto_path, kar_data, data_i, config)
        vaf_df = getVaf(data_i, config)
        cn_ai_df = getAICN(kar_data, config)
        final_df = pd.concat([coverage_df, vaf_df, cn_ai_df], axis=1)

        output_path = os.path.join(output_folder, f'master{i}.tsv')
        final_df.to_csv(output_path, sep='\t', index=False)

        progress_var.set(i + 1)

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
    
    config = default_config.copy()
    config_path = config_file_entry.get()
    if config_path:
        config.update(pdConfig(config_path))

    try:
        process_files(cyto_folder, kars_folder, data_folder, output_folder, config, progress_var, max_progress)
        error_message_var.set("")
    except Exception as e:
        error_message_var.set("Please check that you have the correct folders, files, and file format. If you are not using rAI you might need to change the config file. If the issue persists please reach out to the dev's in the source code below.")
        print(f"Error: {e}")

def open_link(event):
    webbrowser.open_new("https://github.com/philip-hub/cnv-db-tools")

def open_interactive_plot_link(event):
    webbrowser.open_new("https://github.com/philip-hub/CNV-web-app-mock-up")

# GUI
root = Tk()
root.title("Red Panda - rAI 2 Interactive Plot")

# Set minimum size for the window
root.minsize(600, 400)

# Define a font style
font_style = font.Font(family="Helvetica", size=12)

progress_var = StringVar()
max_progress = StringVar()
error_message_var = StringVar()

Label(root, text="Cyto Folder:", font=font_style).grid(row=0, column=0, sticky='e', pady=(5, 0))
cyto_folder_entry = Entry(root, width=60, font=font_style)
cyto_folder_entry.grid(row=0, column=1, pady=(5, 0))
Button(root, text="Browse", command=lambda: browse_folder(cyto_folder_entry), font=font_style).grid(row=0, column=2, pady=(5, 0))

Label(root, text="Kars Folder:", font=font_style).grid(row=1, column=0, sticky='e', pady=(5, 0))
kars_folder_entry = Entry(root, width=60, font=font_style)
kars_folder_entry.grid(row=1, column=1, pady=(5, 0))
Button(root, text="Browse", command=lambda: browse_folder(kars_folder_entry), font=font_style).grid(row=1, column=2, pady=(5, 0))

Label(root, text="Data Folder:", font=font_style).grid(row=2, column=0, sticky='e', pady=(5, 0))
data_folder_entry = Entry(root, width=60, font=font_style)
data_folder_entry.grid(row=2, column=1, pady=(5, 0))
Button(root, text="Browse", command=lambda: browse_folder(data_folder_entry), font=font_style).grid(row=2, column=2, pady=(5, 0))

Label(root, text="Output Folder:", font=font_style).grid(row=3, column=0, sticky='e', pady=(5, 0))
output_folder_entry = Entry(root, width=60, font=font_style)
output_folder_entry.grid(row=3, column=1, pady=(5, 0))
Button(root, text="Browse", command=lambda: browse_folder(output_folder_entry), font=font_style).grid(row=3, column=2, pady=(5, 0))

Label(root, text="Config File (optional):", font=font_style).grid(row=4, column=0, sticky='e', pady=(5, 0))
config_file_entry = Entry(root, width=60, font=font_style)
config_file_entry.grid(row=4, column=1, pady=(5, 0))
Button(root, text="Browse", command=lambda: browse_file(config_file_entry), font=font_style).grid(row=4, column=2, pady=(5, 0))

interactive_plot_link_label = Label(root, text="Create Interactive CNV Plot", font=font_style, fg="blue", cursor="hand2")
interactive_plot_link_label.grid(row=5, column=0, columnspan=3, pady=(5, 0))
interactive_plot_link_label.bind("<Button-1>", open_interactive_plot_link)

Button(root, text="Start Processing", command=start_processing, font=font_style).grid(row=6, column=1, pady=(10, 0))

progress_bar = Progressbar(root, orient='horizontal', mode='determinate', length=500, maximum=100)
progress_bar.grid(row=7, column=0, columnspan=3, pady=(10, 10))

error_message_label = Label(root, textvariable=error_message_var, font=font_style, fg="red", wraplength=500, justify="center")
error_message_label.grid(row=8, column=0, columnspan=3, pady=(5, 0))

link_label = Label(root, text="Open Source Code ✨", font=font_style, fg="blue", cursor="hand2")
link_label.grid(row=9, column=0, columnspan=3)
link_label.bind("<Button-1>", open_link)

def update_progress(*args):
    try:
        progress_bar['value'] = int(progress_var.get()) / int(max_progress.get()) * 100
    except ValueError:
        progress_bar['value'] = 0

progress_var.trace("w", update_progress)

root.mainloop()
