import pandas as pd

def tsv_to_csv(tsv_file_path, csv_file_path):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv_file_path, sep='\t')
    
    # Write the DataFrame to a CSV file
    df.to_csv(csv_file_path, index=False)

# Example usage
tsv_file_path = "Data_D1_karyotype.tsv"  # Replace with your TSV file path
csv_file_path = "karyotype.csv"  # Replace with your desired CSV file path

tsv_to_csv(tsv_file_path, csv_file_path)
