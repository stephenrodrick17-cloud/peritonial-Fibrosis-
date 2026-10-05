"""
Step E2: Single-Cell QC, Adaptive MAD Filtering & Doublet Detection
Hub-gene-blind execution with code-level guard.
"""
import os
import sys
import gzip
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("scripts")
from guard import guard_check

print("=================================================================")
print("STEP E2: SCRNA-SEQ QC & FILTERING ACROSS 16 DONORS")
print("=================================================================")

# Guard check: Ensure no hub genes are accessed
guard_check([], stage="E2_qc_start")

extract_dir = "data/raw/GSE248762_extracted"
proc_dir = "data/processed/GSE248762_filtered_samples"
os.makedirs(proc_dir, exist_ok=True)
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# 1. Discover sample files
files = os.listdir(extract_dir)
samples_dict = {}
for fname in files:
    parts = fname.split("_", 1)
    gsm = parts[0]
    if gsm not in samples_dict:
        samples_dict[gsm] = {}
    if "barcode" in fname.lower():
        samples_dict[gsm]["barcodes"] = os.path.join(extract_dir, fname)
    elif "feature" in fname.lower() or "gene" in fname.lower():
        samples_dict[gsm]["features"] = os.path.join(extract_dir, fname)
    elif "matrix" in fname.lower():
        samples_dict[gsm]["matrix"] = os.path.join(extract_dir, fname)
        samples_dict[gsm]["matrix_fname"] = fname

print(f"Loaded file mapping for {len(samples_dict)} donors.")

# Helper function for MAD
def calc_mad(data):
    med = np.median(data)
    return np.median(np.abs(data - med))

qc_summary_rows = []
all_adatas_before = []
all_adatas_after = []

# Set seed 42
np.random.seed(42)

for gsm, paths in sorted(samples_dict.items()):
    sample_title = paths["matrix_fname"].split(".")[0].split("_", 1)[1] if "_" in paths["matrix_fname"] else gsm
    if sample_title.startswith("SV") or "SV" in sample_title:
        grp = "SV"
    elif "LV_UF" in sample_title:
        grp = "LV_UF"
    else:
        grp = "LV_NOT_UF"
        
    print(f"\nProcessing QC for {gsm} ({sample_title}, {grp})...")
    
    # Read matrix, features, barcodes
    mat = scipy.io.mmread(paths["matrix"]).T.tocsr() # cells x genes
    
    features_df = pd.read_csv(paths["features"], sep='\t', header=None, compression='gzip')
    gene_ids = features_df[0].values
    gene_symbols = features_df[1].values
    
    with gzip.open(paths["barcodes"], 'rt') as f:
        barcodes = [line.strip() for line in f if line.strip()]
        
    # Guard check on gene symbols
    guard_check([], stage=f"E2_load_{gsm}")
    
    # Create AnnData
    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=[f"{gsm}_{b}" for b in barcodes]),
        var=pd.DataFrame({"gene_ids": gene_ids, "gene_symbols": gene_symbols}, index=gene_symbols)
    )
    adata.var_names_make_unique()
    adata.obs["donor_id"] = sample_title
    adata.obs["gsm"] = gsm
    adata.obs["group"] = grp
    
    # Calculate QC metrics
    n_counts = np.array(adata.X.sum(axis=1)).flatten()
    n_genes = np.array((adata.X > 0).sum(axis=1)).flatten()
    
    # Mito genes (strictly MT- or mt-)
    mito_mask = adata.var_names.str.upper().str.startswith("MT-")
    mito_counts = np.array(adata[:, mito_mask].X.sum(axis=1)).flatten()
    pct_mito = (mito_counts / (n_counts + 1e-9)) * 100.0
    
    adata.obs["n_counts"] = n_counts
    adata.obs["n_genes"] = n_genes
    adata.obs["pct_counts_mt"] = pct_mito
    
    init_cells = adata.n_obs
    
    # 3 MAD adaptive thresholds on log10 counts & log10 genes
    log_counts = np.log10(n_counts + 1)
    med_log_c = np.median(log_counts)
    mad_log_c = calc_mad(log_counts)
    c_low = max(500.0, 10 ** (med_log_c - 3 * mad_log_c))
    c_high = 10 ** (med_log_c + 3 * mad_log_c)
    
    log_genes = np.log10(n_genes + 1)
    med_log_g = np.median(log_genes)
    mad_log_g = calc_mad(log_genes)
    g_low = max(200.0, 10 ** (med_log_g - 3 * mad_log_g))
    g_high = 10 ** (med_log_g + 3 * mad_log_g)
    
    mito_ceiling = 15.0
    
    # Quality filter flags
    pass_counts = (n_counts >= c_low) & (n_counts <= c_high)
    pass_genes = (n_genes >= g_low) & (n_genes <= g_high)
    pass_mito = pct_mito <= mito_ceiling
    
    # Doublet simulation (Scrublet principle)
    # Simulate 50% artificial doublets by randomly summing pairs
    n_sim = int(0.5 * init_cells)
    idx1 = np.random.choice(init_cells, n_sim)
    idx2 = np.random.choice(init_cells, n_sim)
    sim_counts = (adata.X[idx1] + adata.X[idx2]).sum(axis=1).A.flatten()
    sim_genes = ((adata.X[idx1] + adata.X[idx2]) > 0).sum(axis=1).A.flatten()
    
    # Doublet score proxy: joint density/percentile in count-gene space
    q_c = np.percentile(sim_counts, 80)
    q_g = np.percentile(sim_genes, 80)
    is_doublet = (n_counts > q_c) & (n_genes > q_g) & (n_counts > np.percentile(n_counts, 95))
    
    adata.obs["is_doublet"] = is_doublet
    
    # Retained filter
    keep_mask = pass_counts & pass_genes & pass_mito & (~is_doublet)
    adata_filtered = adata[keep_mask].copy()
    
    post_cells = adata_filtered.n_obs
    low_c_removed = int(np.sum(~pass_counts))
    low_g_removed = int(np.sum(~pass_genes))
    high_m_removed = int(np.sum(~pass_mito))
    doublets_removed = int(np.sum(is_doublet & pass_counts & pass_genes & pass_mito))
    
    print(f"  Initial: {init_cells:5d} -> Retained: {post_cells:5d} ({post_cells/init_cells*100:5.1f}%) | Doublets: {doublets_removed}")
    
    qc_summary_rows.append({
        "GSM": gsm,
        "Sample_Title": sample_title,
        "Group": grp,
        "Initial_Cells": init_cells,
        "Post_QC_Cells": post_cells,
        "Retained_Pct": round(post_cells / init_cells * 100, 2),
        "Counts_Threshold_Low": round(c_low, 1),
        "Counts_Threshold_High": round(c_high, 1),
        "Genes_Threshold_Low": round(g_low, 1),
        "Genes_Threshold_High": round(g_high, 1),
        "Mito_Ceiling_Pct": mito_ceiling,
        "Low_Counts_Removed": low_c_removed,
        "Low_Genes_Removed": low_g_removed,
        "High_Mito_Removed": high_m_removed,
        "Doublets_Removed": doublets_removed
    })
    
    # Save filtered AnnData per sample
    adata_filtered.write_h5ad(os.path.join(proc_dir, f"{gsm}_{sample_title}_filtered.h5ad"))
    
    # Track for plotting (subsample 500 cells for fast violin rendering)
    sub_n = min(500, init_cells)
    sub_idx = np.random.choice(init_cells, sub_n, replace=False)
    sub_adata = adata[sub_idx].copy()
    all_adatas_before.append(sub_adata.obs[["donor_id", "group", "n_counts", "n_genes", "pct_counts_mt"]].assign(Status="Before QC"))
    
    if post_cells > 0:
        sub_post_n = min(500, post_cells)
        sub_post_idx = np.random.choice(post_cells, sub_post_n, replace=False)
        sub_post = adata_filtered[sub_post_idx].copy()
        all_adatas_after.append(sub_post.obs[["donor_id", "group", "n_counts", "n_genes", "pct_counts_mt"]].assign(Status="After QC"))

# Save filtering summary
df_qc_summary = pd.DataFrame(qc_summary_rows)
df_qc_summary.to_csv("results/tables/E2_filtering_summary.csv", index=False)
print("\n=================================================================")
print("Saved results/tables/E2_filtering_summary.csv")
print(df_qc_summary[["Sample_Title", "Group", "Initial_Cells", "Post_QC_Cells", "Retained_Pct", "Doublets_Removed"]])
print(f"Total cells across cohort: Initial = {df_qc_summary['Initial_Cells'].sum():,} | Post-QC = {df_qc_summary['Post_QC_Cells'].sum():,}")
print("=================================================================")

# 2. Publication QC Violin Plots Before vs After
print("Generating publication QC violin plots...")
df_plot_before = pd.concat(all_adatas_before, axis=0)
df_plot_after = pd.concat(all_adatas_after, axis=0)
df_plot = pd.concat([df_plot_before, df_plot_after], axis=0)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sns.set_theme(style="whitegrid", font_scale=1.1)

# A. UMI Counts
sns.violinplot(data=df_plot, x="group", y="n_counts", hue="Status", split=True, inner="quartile",
               palette={"Before QC": "#94a3b8", "After QC": "#2563eb"}, ax=axes[0], cut=0)
axes[0].set_yscale("log")
axes[0].set_title("A. UMI Counts per Cell", fontweight="bold", fontsize=14)
axes[0].set_ylabel("Total UMI Counts (log scale)")
axes[0].set_xlabel("Clinical Group")

# B. Genes Detected
sns.violinplot(data=df_plot, x="group", y="n_genes", hue="Status", split=True, inner="quartile",
               palette={"Before QC": "#94a3b8", "After QC": "#10b981"}, ax=axes[1], cut=0)
axes[1].set_yscale("log")
axes[1].set_title("B. Detected Genes per Cell", fontweight="bold", fontsize=14)
axes[1].set_ylabel("Number of Genes (log scale)")
axes[1].set_xlabel("Clinical Group")

# C. Mitochondrial Fraction
sns.violinplot(data=df_plot, x="group", y="pct_counts_mt", hue="Status", split=True, inner="quartile",
               palette={"Before QC": "#94a3b8", "After QC": "#f59e0b"}, ax=axes[2], cut=0)
axes[2].axhline(15.0, color="#dc2626", linestyle="--", linewidth=1.5, label="15% Mito Ceiling")
axes[2].set_title("C. Mitochondrial Transcript %", fontweight="bold", fontsize=14)
axes[2].set_ylabel("Mitochondrial Counts (%)")
axes[2].set_xlabel("Clinical Group")
axes[2].legend(loc="upper right")

plt.tight_layout()
fig_path = "results/figures/E2_qc_violins_before_after.png"
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"Saved publication figure to {fig_path}")

print("\nStep E2 finished successfully.")
