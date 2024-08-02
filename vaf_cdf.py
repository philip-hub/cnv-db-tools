import pandas as pd
import numpy as np

vaf_file_path = "vaf/vaf_coverage_with_x_and_median.tsv"
vaf_df = pd.read_csv(vaf_file_path, sep="\t")

#sort V within each arm and then apply the line space
def sort_and_create_X_column(group):
    group = group.sort_values(by='v')
    group_size = len(group)
    group['X'] = np.linspace(0, 1, group_size)
    return group

sorted_vaf_df = vaf_df.groupby('arm').apply(sort_and_create_X_column).reset_index(drop=True)

sorted_vaf_file_path = "vaf/vaf_sorted_by_arm_and_v_with_X.csv"
sorted_vaf_df.to_csv(sorted_vaf_file_path, index=False)