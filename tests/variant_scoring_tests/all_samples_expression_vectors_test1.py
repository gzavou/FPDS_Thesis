#!/usr/bin/env python3

import polars as pl
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Process Borzoi SED output for a given chromosome.")
parser.add_argument("--chrom", required=True, help="Number of chromosome (e.g., 18)")
args = parser.parse_args()
chrom = args.chrom
#%%
chrom = 10  # For testing

# File Paths
vcf_path = f"/pool01/projects/abante_lab/ao_prediction_enrollhd_2024/enroll_hd/regulatory_vcfs/gwa12345.mis4.9064.hg38.cisreg0.5.mad0.01.chr{chrom}.vcf" 
results_path = f"/pool01/projects/abante_lab/genomic_llms/borzoi/proc_results/weighted_logSED/chr{chrom}_weighted_predictions.tsv.gz" 


# Load VCF 
vcf = pl.read_csv(vcf_path, separator="\t", has_header=True, comment_prefix="##")
vcf = vcf.rename({"#CHROM": "CHROM"})

# Identify sample columns
meta_cols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
sample_names = [col for col in vcf.columns if col not in meta_cols]

# ALT genotypes of interest
alt_het_genos = {"0/1", "1/0"}

# Load predictions
predictions = pd.read_csv(results_path, sep="\t")
tissues = ["RNA-seq: Putamen", "RNA-seq: Caudate"]
#%%
# Initialize: {tissue: list of rows}
all_tissues = {}

sample_names = sample_names[:1]  # Limit to first 5 for testing

# Process each sample
for sample_name in sample_names:
    # Get ALT SNPs for this sample
    het_snps = vcf.filter(pl.col(sample_name).is_in(alt_het_genos)).select("ID").to_series().to_list()
    hom_snps = vcf.filter(pl.col(sample_name) == "1/1").select("ID").to_series().to_list()
    het_snps = set(str(s).strip() for s in het_snps)
    hom_snps = set(str(s).strip() for s in hom_snps)

    # Select SNPs 
    sample_preds = predictions[predictions["snp"].isin(het_snps | hom_snps)]
#%%
    # Double the logSED for alternative homozygous SNPs
    sample_preds[sample_preds["snp"].isin(hom_snps)]["weighted_logSED"] *= 2

    #  aggregate per gene
    summed = sample_preds.groupby("gene")["weighted_logSED"].sum().reset_index()
    summed["sample"] = sample_name
    all_tissues[sample_name] = [summed]

    #all_tissues[sample].append(summed)
# %%
# Write output for each tissue 
for tissue, rows in tissue_long_rows.items():
    if not rows:
        continue
    df_long = pd.concat(rows, ignore_index=True)
    matrix = df_long.pivot(index="sample", columns="gene", values="logSED").fillna(0)
    matrix.to_csv(f"expression_matrix2_chr18_{tissue}.tsv", sep="\t")



# %%
