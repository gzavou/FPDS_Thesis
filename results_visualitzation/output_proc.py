#!/usr/bin/env python3

#%%## score_variants output rearrangement ###
# run with enrollhd conda environment

import h5py
import pandas as pd

# Path to save results
outdir = "/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/"

# Path of borzoi output
output_path = "/pool01/code/src/genomic_llms/borzoi/tutorials/latest/score_variants/snp_sed/f0c0/sed.h5"

#%%## Open the .h5 file ###
sed_h5 = h5py.File(output_path,'r')

# Extract variables
snp_list = [s.decode() for s in sed_h5['snp']]
gene_list = [g.decode() for g in sed_h5['gene']]
snp_idx = sed_h5['si'][:]
logseds = sed_h5['logSED'][:]
target_labels = [t.decode() for t in sed_h5['target_labels']]
alt_alleles = [a.decode() for a in sed_h5['alt_allele']]

results = []

for row_ix, snp_gene_idx in enumerate(snp_idx):
    snp = snp_list[snp_gene_idx]
    gene = gene_list[row_ix]
    alt_allele = alt_alleles[snp_gene_idx]  
    for target_ix, target_label in enumerate(target_labels):
        results.append({
            "snp": snp,
            "gene": gene,
            "alt_allele": alt_allele,
            "tissue": target_label,
            "logSED": logseds[row_ix, target_ix]
        })

# Convert results to a pandas DataFrame
df = pd.DataFrame(results)

# Save the DataFrame to a tab separated file
df.to_csv(f"{outdir}sed.txt", sep="\t", index=False)