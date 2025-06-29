import polars as pl
import pandas as pd

vcf_path = "/pool01/projects/abante_lab/ao_prediction_enrollhd_2024/enroll_hd/regulatory_vcfs/gwa12345.mis4.9064.hg38.cisreg0.5.mad0.01.chr15.vcf" 
results_path = "/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_02.txt" 
sample_name = "gwa1.101_gwa1.101"  # one sample 

vcf = pl.read_csv(vcf_path, separator="\t", has_header=True, comment_prefix="##")
vcf = vcf.rename({"#CHROM": "CHROM"})

# Confirm sample exists
if sample_name not in vcf.columns:
    raise ValueError(f"Sample '{sample_name}' not found in VCF columns: {vcf.columns}")

# Filter to ALT genotypes
alt_genos = {"0/1", "1/0", "1/1"}
vcf = vcf.filter(pl.col(sample_name).is_in(alt_genos))

# Extract SNP IDs
vcf = vcf.with_columns([pl.col("ID").alias("snp")])
vcf_snps = vcf.select("snp").to_series().to_list()
vcf_snps = [str(s).strip() for s in vcf_snps]

# Load Predictions 
predictions = pd.read_csv(results_path, sep="\t")

# Filter relevant SNPs
filtered = predictions[predictions["snp"].isin(vcf_snps)]
print(f"\nTotal matching SNPs for sample '{sample_name}': {len(filtered)}\n")

# Group and print logSED values per gene
grouped = filtered.groupby("gene")

for gene, group in grouped:
    logsed_values = group["logSED"].tolist()
    print(f"Gene: {gene}")
    print(f"  logSED values: {logsed_values}")
    print(f"  Sum: {sum(logsed_values)}\n")

# Aggregate final expression vector
expression_vector = grouped["logSED"].sum().reset_index()
expression_vector.columns = ["gene", "predicted_logFC"]

# Save result
output_path = f"{sample_name}_expression_vector_chr15.tsv"
expression_vector.to_csv(output_path, sep="\t", index=False)

print(f"\nSaved expression vector to: {output_path}")
