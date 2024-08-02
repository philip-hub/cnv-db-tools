import pandas as pd
import numpy as np

vaf_file_path = "vaf/vaf_coverage_with_x_and_median.tsv"
vaf_df = pd.read_csv(vaf_file_path, sep="\t")

def sort_and_create_X_column(group):
    group = group.sort_values(by='v')
    group_size = len(group)
    group['Y4'] = np.linspace(0, 1, group_size)
    return group

sorted_vaf_df = vaf_df.groupby('arm', group_keys=False).apply(lambda group: sort_and_create_X_column(group.reset_index(drop=True)))
sorted_vaf_df = sorted_vaf_df.rename(columns={'v': 'X4', 'arm': 'arm4'})

# sorted_vaf_df['X4'] = sorted_vaf_df['X4'].apply(lambda x: round(x, 4 - int(np.floor(np.log10(abs(x)))) - 1) if x != 0 else 0)
# sorted_vaf_df['Y4'] = sorted_vaf_df['Y4'].apply(lambda x: round(x, 4 - int(np.floor(np.log10(abs(x)))) - 1) if x != 0 else 0)

final_vaf_df = sorted_vaf_df[['arm4', 'X4', 'Y4']]

final_vaf_file_path = "vaf/vaf_sorted_by_arm_and_v_with_X_final.csv"
final_vaf_df.to_csv(final_vaf_file_path, index=False)
