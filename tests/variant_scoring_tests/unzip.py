import gzip
import pandas as pd

# Set path to one of your files
file_path = "/pool01/projects/abante_lab/genomic_llms/borzoi/proc_results/expression_vectors/caudate/caudate_expression_matrix_chr7.txt.gz"

# Open and load a preview
with gzip.open(file_path, 'rt') as f:
    # Read only the first few lines (for safety)
    for _ in range(10):
        print(f.readline())

# Optional: read into DataFrame
df = pd.read_csv(file_path, sep='\t', compression='gzip', nrows=5)
print("\nDataFrame Preview:")
print(df.head())
