import pandas as pd
import matplotlib.pyplot as plt

# File path to the CSV file
csv_file_path = "coverage_with_x_and_median.csv"

# Read the CSV file into a DataFrame
data = pd.read_csv(csv_file_path)

# Check if the necessary columns exist
if 'x' in data.columns and 'y' in data.columns and 'arm' in data.columns:
    # Create a scatter plot with smaller dots
    plt.figure(figsize=(20, 4))
    plt.scatter(data['x'], data['y'], s=10, alpha=0.5)  # 's' parameter controls the size of the dots
    
    # Add title and labels
    plt.title('Coverage Plot')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

    # Add arm labels as tickers on the x-axis
    x_ticks = []
    x_labels = []
    previous_arm = None
    
    for i, (x, arm) in enumerate(zip(data['x'], data['arm'])):
        if arm != previous_arm:
            x_ticks.append(x)
            x_labels.append(arm)
            previous_arm = arm
    
    plt.xticks(ticks=x_ticks, labels=x_labels, rotation=90)  # Rotate labels for better visibility
    
    # Show the plot
    plt.tight_layout()
    plt.show()
else:
    print("The required columns 'x', 'y', and 'arm' are not present in the CSV file.")
