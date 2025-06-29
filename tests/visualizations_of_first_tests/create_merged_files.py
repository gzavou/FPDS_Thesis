import pandas as pd

# load files
borzoi = pd.read_csv('/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_01.txt', sep='\t')
ours = pd.read_csv('/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_02.txt', sep='\t')

# Filter to FAN1
fan1_borzoi = borzoi[borzoi["gene"] == "ENSG00000198690"].copy()
fan1_ours = ours[ours["gene"] == "ENSG00000198690"].copy()

fan1_borzoi["source"] = "borzoi_ref"
fan1_ours["source"] = "our_ref"

# concat both sets
fan1_all = pd.concat([fan1_borzoi, fan1_ours], ignore_index=True)

fan1_all["tissue"] = fan1_all["tissue"].str.strip().str.replace("-", "").str.replace("seq:", ":")

# Save combined file
fan1_all.to_csv("fan1_all_tissues_predictions.tsv", sep="\t", index=False)
