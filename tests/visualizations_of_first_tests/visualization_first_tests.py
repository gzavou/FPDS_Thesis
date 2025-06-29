import pandas as pd
borzoi = pd.read_csv("/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_01.txt", sep="\t")
ours = pd.read_csv("/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_02.txt", sep="\t")

# Combine both sources into one DataFrame
df = pd.concat([borzoi, ours], ignore_index=True)

# Optional cleanup
df["tissue"] = df["tissue"].str.strip()

# Pivot to show SNPs across tissues
pivot = df.pivot_table(index="snp", columns="tissue", values="logSED")

# Calculate difference between tissues
pivot["diff"] = pivot.iloc[:, 0] - pivot.iloc[:, 1]  
pivot = pivot.reset_index()

# Save for inspection
pivot.to_csv("snp_tissue_comparison.tsv", sep="\t", index=False)


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load your merged tissue comparison file
df = pd.read_csv("fan1_tissue_comparison.tsv", sep="\t")

# Convert wide format to long format for plotting
melted = pd.melt(
    df,
    id_vars=["snp"],
    value_vars=["logSED_K562", "logSED_Putamen"],
    var_name="tissue",
    value_name="logSED"
)

# Make tissue names cleaner 
melted["tissue"] = melted["tissue"].str.replace("logSED_", "")

# Create the barplot
plt.figure(figsize=(12, 6))
sns.barplot(data=melted, x="snp", y="logSED", hue="tissue")

# Style and save
plt.title("Per-SNP Predicted Regulatory Effect (logSED) on FAN1")
plt.ylabel("logSED")
plt.xlabel("SNP")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("fan1_snp_by_tissue_barplot.png", dpi=300)



# Load the SNP comparison file
df = pd.read_csv("fan1_tissue_comparison.tsv", sep="\t")

# Create scatter plot of SNP logSED between tissues
plt.figure(figsize=(7, 6))
sns.scatterplot(data=df, x="logSED_K562", y="logSED_Putamen")

# Diagonal line = equal effect
plt.axline((0, 0), slope=1, color="red", linestyle="--", label="Equal Effect")

plt.title("FAN1 SNP Effects: K562 vs Putamen")
plt.xlabel("logSED in K562")
plt.ylabel("logSED in Putamen")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save plot as image
plt.savefig("fan1_cisregulatory_snp_scatter_k562_vs_putamen.png", dpi=300)


