import pandas as pd
import matplotlib.pyplot as plt

csv_file_path = "coverage_with_x_and_median.csv"
data = pd.read_csv(csv_file_path)

plt.figure(figsize=(10, 6))
plt.scatter(data['x'], data['y'], s=1, alpha=0.5)
plt.title('Scatter Plot of x vs y')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.show()

