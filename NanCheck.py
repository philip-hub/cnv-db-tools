import pandas as pd

def check_nans(file_path):
    # Load the dataset into a DataFrame
    df = pd.read_csv(file_path, sep='\t')

    # Check for NaNs in each column
    nan_summary = df.isna().sum()

    # Print the count of NaNs for each column
    print("NaN counts per column:")
    print(nan_summary)

    # Optionally, you can also print rows with NaNs
    print("\nRows with NaN values:")
    print(df[df.isna().any(axis=1)])

# Example usage
file_path = "SJALL003310_D3.tsv"  # Replace with your TSV file path
check_nans(file_path)
