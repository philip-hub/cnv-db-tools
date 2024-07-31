import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV file
df = pd.read_csv('SJALL003310_D3_karyotype.tsv')

# Create the scatter plot
plt.figure(figsize=(6, 6))
plt.scatter(df['X2'], df['Y2'], c='blue', alpha=0.5)

# Add titles and labels
plt.title('Scatter Plot of AI vs CN')
plt.xlabel('CN')
plt.ylabel('AI')

# Set the limits for x and y axes
# plt.xlim(1, 3)
# plt.ylim(0, 0.5)

# Show the plot
plt.show()
