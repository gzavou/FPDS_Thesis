import pandas as pd
import mygene
import matplotlib.pyplot as plt
import seaborn as sns


borzoi = pd.read_csv("/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test1_with_symbols.tsv", sep="\t")
ours = pd.read_csv("/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test2_with_symbols.tsv", sep="\t")

#  Combine both sources into one DataFrame 
all_preds = pd.concat([borzoi, ours], ignore_index=True)

#  Aggregate logSED per gene per tissue 
agg = all_preds.groupby(["symbol", "tissue"])["logSED"].sum().reset_index()

# Save the aggregate data
agg.to_csv("gene_tissue_logsed_summary.tsv", sep="\t", index=False)

top_genes = agg.sort_values("logSED", ascending=False).head(20)

plt.figure(figsize=(12, 6))
sns.barplot(data=top_genes, x="symbol", y="logSED", hue="tissue")
plt.xticks(rotation=45, ha='right')
plt.title("Top Genes by Sum of Predicted Regulatory Effect (logSED)")
plt.tight_layout()

plt.savefig("genes_logsed_by_tissue.png", dpi=300)
