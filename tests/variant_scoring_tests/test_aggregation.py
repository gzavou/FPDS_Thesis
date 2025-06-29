
import polars as pl
import pandas as pd
import os

# File Paths
vcf_path = "/pool01/projects/abante_lab/ao_prediction_enrollhd_2024/enroll_hd/regulatory_vcfs/gwa12345.mis4.9064.hg38.cisreg0.5.mad0.01.chr18.vcf" 
results_path = "/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/sed_all_chr/chr18.txt.gz" 

sample_names = ["gwa1.100_gwa1.100"]
tissues = ["RNA-seq: Putamen"]

# Load VCF 
vcf = pl.read_csv(vcf_path, separator="\t", has_header=True, comment_prefix="##")
vcf = vcf.rename({"#CHROM": "CHROM"})

# Identify ALT-genotype SNPs for sample 
meta_cols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
alt_genos = {"0/1", "1/0", "1/1"}

sample_snps = {}
for sample in sample_names:
    snps = vcf.filter(pl.col(sample).is_in(alt_genos)).select("ID").to_series().to_list()
    sample_snps[sample] = set(str(s).strip() for s in snps)


predictions = pd.read_csv(results_path, sep="\t")

# Process each tissue 
for tissue in tissues:
    tissue_preds = predictions[predictions["tissue"] == tissue]
    print(f"\nTissue: {tissue} | Prediction rows: {len(tissue_preds)}")

    long_rows = []
    for sample, snps in sample_snps.items():
        matched = tissue_preds[tissue_preds["snp"].isin(snps)]

        if matched.empty:
            print(f"No matching SNPs for sample {sample}.")
            continue

        print(f"\nMatching SNPs for sample '{sample}': {len(matched)}")

        # Debug Print first 5 genes and their SNPs/logSEDs
        printed = 0
        for gene, group in matched.groupby("gene"):
            logs = group["logSED"].tolist()
            gene_snps = group["snp"].tolist()
            total = sum(logs)

            print(f"\nGene: {gene}")
            print(f"  logSED values: {logs}")
            print(f"  SNPs used: {gene_snps}")
            print(f"  Sum: {total}")

            printed += 1
            if printed >= 5:
                break

        # Aggregate for matrix
        summed = matched.groupby("gene")["logSED"].sum().reset_index()
        summed["sample"] = sample
        long_rows.append(summed)

    if long_rows:
        df_long = pd.concat(long_rows, ignore_index=True)
        matrix = df_long.pivot(index="sample", columns="gene", values="logSED").fillna(0)

        # Save file
        matrix.to_csv(f"expression_matrix_chr18_{tissue}.tsv", sep="\t")
        print(f"\nSaved expression matrix for tissue '{tissue}'.")
    else:
        print(f"No data to save for tissue '{tissue}'.")