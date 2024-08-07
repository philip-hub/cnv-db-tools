import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

csv_file_path = "coverage/coverage_with_x_and_median.csv"
data = pd.read_csv(csv_file_path)

if 'X1' in data.columns and 'Y1' in data.columns and 'arm1' in data.columns:
    plt.figure(figsize=(20, 3))
    plt.scatter(data['X1'].values, data['Y1'].values, s=5, alpha=0.5, c='lightgreen')  # s parameter controls the size of the dots
    plt.title('Coverage Plot')
    plt.xlabel('position')
    plt.ylabel('log2(median/ref)')
    plt.grid(True)
    

    x_ticks = []
    x_labels = []
    previous_arm = None
    
    for i, (x, arm) in enumerate(zip(data['X1'], data['arm1'])):
        if arm != previous_arm:
            x_ticks.append(x)
            x_labels.append(arm)
            previous_arm = arm
    
    plt.xticks(ticks=x_ticks, labels=x_labels, rotation=90)
    
    # Show the plot
    plt.tight_layout()
    plt.show()
else:
    print("The required columns 'x', 'y', and 'arm' are not present in the CSV file.")