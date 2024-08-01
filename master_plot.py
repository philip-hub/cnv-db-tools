import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df1 = pd.read_csv('master/master.tsv', sep="\t")
print(df1)


plt.figure(figsize=(20, 20))

# AI vs CN Plot
plt.subplot(3, 1, 1)
if 'X3' in df1.columns and 'Y3' in df1.columns:
    plt.scatter(df1['X3'], df1['Y3'], c='blue', alpha=0.5)
    plt.title('Scatter Plot of AI vs CN')
    plt.xlabel('CN')
    plt.ylabel('AI')
else:
    print("The required columns 'X3' and 'Y3' are not present in the CSV file.")

#Coverage Plot
plt.subplot(3, 1, 2)
if 'X1' in df1.columns and 'Y1' in df1.columns and 'arm1' in df1.columns:
    plt.scatter(df1['X1'].values, df1['Y1'].values, s=5, alpha=0.5, c='lightgreen')
    plt.title('Coverage Plot')
    plt.xlabel('position')
    plt.ylabel('log2(median/ref)')
    plt.grid(True)

    x_ticks = []
    x_labels = []
    previous_arm = None
    
    for i, (x, arm) in enumerate(zip(df1['X1'], df1['arm1'])):
        if arm != previous_arm:
            x_ticks.append(x)
            x_labels.append(arm)
            previous_arm = arm
    
    plt.xticks(ticks=x_ticks, labels=x_labels, rotation=90)
else:
    print("The required columns 'X1', 'Y1', and 'arm1' are not present in the CSV file.")

#Vaf Plot
plt.subplot(3, 1, 3)
if 'X2' in df1.columns and 'Y2' in df1.columns and 'arm2' in df1.columns:
    plt.scatter(df1['X2'], df1['Y2'], s=5, alpha=0.5, c='lightgreen')
    plt.title('Vaf Plot')
    plt.xlabel('position')
    plt.ylabel('Vaf')
    plt.grid(True)

    x_ticks = []
    x_labels = []
    previous_arm = None
    
    for i, (x, arm) in enumerate(zip(df1['X2'], df1['arm2'])):
        if arm != previous_arm:
            x_ticks.append(x)
            x_labels.append(arm)
            previous_arm = arm
    
    plt.xticks(ticks=x_ticks, labels=x_labels, rotation=90)
else:
    print("The required columns 'X2', 'Y2', and 'arm2' are not present in the CSV file.")


plt.tight_layout()
plt.savefig('plots/combined_plots.png')
plt.show()
