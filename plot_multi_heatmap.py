import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

filepath = "multi-output/multi.tsv"

df = pd.read_csv(filepath, sep="\t")

def create_heatmap(data, value_column, title, color):
    pivot_df = data.pivot_table(index="sample", columns="arm", values=value_column, aggfunc='mean')
    
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(pivot_df, cmap=color, annot=False, linewidths=.5, cbar_kws={'label': value_column})
    plt.title(title)
    plt.xlabel('Arm')
    plt.ylabel('Sample Number')
    plt.show()

create_heatmap(df, 'cn', 'CN Heatmap','coolwarm')

create_heatmap(df, 'ai', 'AI Heatmap','cool')
