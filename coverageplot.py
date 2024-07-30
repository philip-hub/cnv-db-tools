import pandas as pd
import matplotlib.pyplot as plt

# File path to the CSV file
csv_file_path = "coverage_with_x_and_median.csv"

# Read the CSV file into a DataFrame
data = pd.read_csv(csv_file_path)

# Check if the necessary columns exist
if 'x' in data.columns and 'y' in data.columns:
    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(data['x'], data['y'], alpha=0.5)
    plt.title('Scatter Plot of x vs y')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
else:
    print("The required columns 'x' and 'y' are not present in the CSV file.")
