import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("/pool01/projects/abante_lab/ao_prediction_enrollhd_2024/genes/gene_cisreg.txt", sep="\t")

assert "feature name" in df.columns, "Missing 'feature name' column"
assert "score" in df.columns, "Missing 'score' column"

df["score"] = pd.to_numeric(df["score"], errors="coerce")

# Plot enhancer score distribution
plt.figure(figsize=(10, 5))
sns.histplot(data=df[df["feature name"] == "Enhancer"], x="score", bins=30, kde=True)
plt.title("Enhancer Score Distribution")
plt.xlabel("Enhancer Score")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("enhancer_score_distribution.png", dpi=300)



