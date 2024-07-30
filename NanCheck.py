import pandas as pd

def check_nans(file_path):

    df = pd.read_csv(file_path, sep='\t')
    nan_summary = df.isna().sum()
    print("NaN counts per column:")
    print(nan_summary)
    print("\nRows with NaN values:")
    print(df[df.isna().any(axis=1)])

file_path = "SJALL003310_D3.tsv"  # Replace with your TSV file path
check_nans(file_path)
