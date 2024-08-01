import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
df = pd.read_csv('CN_AI/karotype_graph.csv')

plt.figure(figsize=(6, 6))
plt.scatter(df['X2'], df['Y2'], c='blue', alpha=0.5)

plt.title('Scatter Plot of AI vs CN')
plt.xlabel('CN')
plt.ylabel('AI')
# plt.xlim(1, 3)
# plt.ylim(0, 0.5)

plt.show()