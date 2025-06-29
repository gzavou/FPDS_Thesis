import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

borzoi = pd.read_csv("/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_01.txt", sep="\t")
ours = pd.read_csv("/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_02.txt", sep="\t")

# Combine both sources into one DataFrame
all_preds = pd.concat([borzoi, ours], ignore_index=True)

# Aggregate logSED per gene per tissue 
agg = all_preds.groupby(["gene", "tissue"])["logSED"].sum().reset_index()

# Save the output 
agg.to_csv("gene_tissue_logsed_summary.tsv", sep="\t", index=False)


# Create the plot
plt.figure(figsize=(6, 4))
sns.barplot(data=agg, x="tissue", y="logSED")
plt.title("Total Predicted Effect Across Tissues")
plt.ylabel("Sum of logSED")
plt.xticks(rotation=45)
plt.tight_layout()

# Save 
plt.savefig("fan1_cisregulatory_effect_by_tissue.png", dpi=300)



