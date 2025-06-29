import os
import gzip
import pandas as pd
from glob import glob
from joblib import Parallel, delayed
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import linkage, leaves_list
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

def load_expression_matrix(folder_path):
    """
    Loads and merges all .txt.gz files from a folder into one big dataframe.
    Assumes rows = subjects, columns = genes.
    """
    all_files = sorted(glob(os.path.join(folder_path, "*.txt.gz")))
    print(f"Found {len(all_files)} files in {folder_path}")
    
    df_list = []
    for f in tqdm(all_files, desc="Reading files"):
        df = pd.read_csv(f, sep="\t", compression="gzip")
        df_list.append(df.set_index("sample"))  # first col is subject/sample ID

    # Concatenate along columns (genes)
    full_df = pd.concat(df_list, axis=1)
    return full_df

def compute_correlation_matrix(caudate_expr, putamen_expr, n_jobs=-1):
    """
    Computes Pearson correlation matrix between genes in caudate and putamen.
    Rows: putamen genes, Cols: caudate genes
    """
    # Transpose: genes x subjects
    caudate_expr = caudate_expr.T
    putamen_expr = putamen_expr.T

    # Keep only shared subjects
    shared_subjects = caudate_expr.columns.intersection(putamen_expr.columns)
    caudate_expr = caudate_expr[shared_subjects]
    putamen_expr = putamen_expr[shared_subjects]

    caudate_genes = caudate_expr.index
    putamen_genes = putamen_expr.index

    def correlate_gene_pair(gene_p):
        vec_p = putamen_expr.loc[gene_p]
        row = []
        for gene_c in caudate_genes:
            vec_c = caudate_expr.loc[gene_c]
            corr = pearsonr(vec_p, vec_c)[0]
            if pd.isna(corr): corr = 0
            row.append(corr)
        return row

    print("Computing correlations in parallel...")
    result = Parallel(n_jobs=n_jobs, batch_size=10)(
        delayed(correlate_gene_pair)(g) for g in tqdm(putamen_genes)
    )

    corr_df = pd.DataFrame(result, index=putamen_genes, columns=caudate_genes)
    return corr_df

def cluster_and_plot(corr_df, output="heatmap.png"):
    # Cluster rows and cols
    row_linkage = linkage(corr_df, method="average")
    col_linkage = linkage(corr_df.T, method="average")
    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)

    clustered = corr_df.iloc[row_order, col_order]

    # Plot
    plt.figure(figsize=(12, 10))
    sns.heatmap(clustered, cmap="vlag", xticklabels=False, yticklabels=False)
    plt.title("Gene Expression Correlation: Caudate vs Putamen")
    plt.xlabel("Caudate genes")
    plt.ylabel("Putamen genes")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.show()

# === MAIN SCRIPT ===
caudate_dir = "/pool01/projects/abante_lab/genomic_llms/borzoi/proc_results/expression_vectors/caudate"
putamen_dir = "/pool01/projects/abante_lab/genomic_llms/borzoi/proc_results/expression_vectors/putamen"

print("Loading expression data...")
caudate_df = load_expression_matrix(caudate_dir)
putamen_df = load_expression_matrix(putamen_dir)

print("Computing gene-gene correlation matrix...")
corr_matrix = compute_correlation_matrix(caudate_df, putamen_df)

print("Plotting clustered heatmap...")
cluster_and_plot(corr_matrix)
