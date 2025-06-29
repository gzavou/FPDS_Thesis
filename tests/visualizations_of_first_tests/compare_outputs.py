import pandas as pd
import matplotlib.pyplot as plt

borzoi = pd.read_csv('/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_01.txt', sep='\t')
ours = pd.read_csv('/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_02.txt', sep='\t')

# Filter to FAN1 gene
borzoi_fan1 = borzoi[borzoi["gene"] == "ENSG00000198690"]
ours_fan1 = ours[ours["gene"] == "ENSG00000198690"]

# Rename columns for clarity
borzoi_fan1.rename(columns={"logSED": "logSED_K562"}, inplace=True)
ours_fan1.rename(columns={"logSED": "logSED_Putamen"}, inplace=True)

# Merge on SNP and allele
merged = pd.merge(borzoi_fan1, ours_fan1, on=["snp", "gene", "alt_allele"])

# Calculate difference
merged["diff"] = merged["logSED_Putamen"] - merged["logSED_K562"]

# Save the result
merged.to_csv("fan1_tissue_comparison.tsv", sep="\t", index=False)

plt.figure(figsize=(10, 6))
plt.hist(merged["logSED_K562"], bins=50, alpha=0.5, label='K562')
plt.hist(merged["logSED_Putamen"], bins=50, alpha=0.5, label='Putamen')
plt.hist(merged["diff"], bins=50, label='Difference')
plt.legend()
plt.show()

