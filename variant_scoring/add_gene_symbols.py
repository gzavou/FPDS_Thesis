import pandas as pd
import mygene
import os
import argparse

parser = argparse.ArgumentParser(description="Weight Borzoi logSED values.")
parser.add_argument("--chomosome", required=True, help="chromosome number to process (ex: 18)")
args = parser.parse_args()
chromosome = args.chomosome

# Set the working directory
os.chdir("/pool01/projects/abante_lab/genomic_llms/borzoi/")

df1 = pd.read_csv(f"test_outputs/sed_all_chr/chr{chromosome}.txt", sep="\t")

#  Combine to get all unique gene IDs 
combined = pd.concat([df1], ignore_index=True)
ensembl_ids = combined["gene"].unique().tolist()

# Query MyGene for gene symbol
mg = mygene.MyGeneInfo()
gene_info = mg.querymany(
    ensembl_ids,
    scopes="ensembl.gene",
    fields="symbol",
    species="human"
)

# Create mapping DataFrame 
gene_df = pd.DataFrame(gene_info)
if "symbol" in gene_df.columns:
    gene_map = gene_df[["query", "symbol"]].dropna()
else:
    gene_map = pd.DataFrame(columns=["query", "symbol"])

# Merge initial mapping into both files
df1 = df1.merge(gene_map, left_on="gene", right_on="query", how="left").drop(columns=["query"])

#  Save the outputs
df1.to_csv(f"test_outputs/chr{chromosome}_with_symbols.tsv", sep="\t", index=False)

