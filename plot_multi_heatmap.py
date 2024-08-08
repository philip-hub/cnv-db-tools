import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        super().__init__(vmin, vmax, clip)
        
    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint
        
        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between max value and min value.")
        elif vmin == vmax:
            return np.ma.zeros_like(value)
        
        rescaled_mid = 0.5 * (midpoint - vmin) / (vmax - vmin)
        result = np.ma.masked_array(np.interp(result, [vmin, midpoint, vmax], [0, rescaled_mid, 1]),
                                    mask=np.ma.getmask(result))
        if is_scalar:
            result = np.ma.masked
        return result

filepath = "multi-output/multi.tsv"

df = pd.read_csv(filepath, sep="\t")

def create_heatmap(data, value_column, title, cmap, midpoint):
    pivot_df = data.pivot_table(index="sample", columns="arm", values=value_column, aggfunc='mean')
    
    plt.figure(figsize=(10, 8))
    norm = MidpointNormalize(vmin=pivot_df.min().min(), vmax=pivot_df.max().max(), midpoint=midpoint)
    ax = sns.heatmap(pivot_df, cmap=cmap, annot=False, linewidths=.5, cbar_kws={'label': value_column}, norm=norm)
    plt.title(title)
    plt.xlabel('Arm')
    plt.ylabel('Sample Number')
    plt.show()

create_heatmap(df, 'cn', 'CN Heatmap', 'coolwarm', midpoint=2)

create_heatmap(df, 'ai', 'AI Heatmap', 'cool', midpoint=df['ai'].mean())  # Assuming the midpoint for AI is its mean value
