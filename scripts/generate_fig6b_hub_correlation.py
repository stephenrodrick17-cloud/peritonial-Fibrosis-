"""
Generate Figure 6B: 11 Consensus Hub Gene Pearson Correlation Co-expression Heatmap.

Uses GSE62928 MaxMean collapsed full expression matrix (N=8 samples) to compute
pairwise Pearson correlation among the 11 ML consensus hub genes. Matches the
publication-quality aesthetic (seaborn vlag cmap, center at 0, publication-style
annotations, 300 DPI PNG + vector PDF) used throughout the pipeline.

Outputs:
  - results/figures/Hub_05_hub_gene_pearson_correlation_heatmap.png
  - results/figures/Hub_05_hub_gene_pearson_correlation_heatmap.pdf
  - results/figures/source_data/fig6b_hub_gene_correlation_source.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# ---- Directories ---------------------------------------------------------
TABLES_DIR = r"d:\Peritoneal Project\results\tables"
FIG_DIR = r"d:\Peritoneal Project\results\figures"
SRC_DIR = os.path.join(FIG_DIR, "source_data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(SRC_DIR, exist_ok=True)

# ---- Aesthetic parameters (match pipeline defaults) ----------------------
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#2B2B2B'
plt.rcParams['axes.linewidth'] = 0.9

# ---- Hub gene panel (canonical ordered from B_gse62928_recomputed_hub_genes) ---
HUB_GENES = [
    "ISM1", "FN1", "EDIL3", "VCAN", "COL3A1",
    "COMP", "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"
]

# ---- Load expression matrix ---------------------------------------------
expr_path = os.path.join(TABLES_DIR, "GSE62928_full_expression_matrix.csv")
df_expr = pd.read_csv(expr_path, index_col=0)
print(f"Loaded expression matrix: {df_expr.shape[0]} genes x {df_expr.shape[1]} samples")

missing = [g for g in HUB_GENES if g not in df_expr.index]
if missing:
    raise ValueError(f"Missing hub genes in expression matrix: {missing}")

# Extract hub genes x samples (N=8)
hub_expr = df_expr.loc[HUB_GENES]

# ---- Compute Pearson r matrix + two-sided p-values -----------------------
n_genes = len(HUB_GENES)
r_mat = pd.DataFrame(np.ones((n_genes, n_genes)), index=HUB_GENES, columns=HUB_GENES, dtype=float)
p_mat = pd.DataFrame(np.zeros((n_genes, n_genes)), index=HUB_GENES, columns=HUB_GENES, dtype=float)

for i, g1 in enumerate(HUB_GENES):
    for j, g2 in enumerate(HUB_GENES):
        if i == j:
            r_mat.iloc[i, j] = 1.0
            p_mat.iloc[i, j] = 0.0
        else:
            r, p = pearsonr(hub_expr.loc[g1], hub_expr.loc[g2])
            r_mat.iloc[i, j] = r
            p_mat.iloc[i, j] = p

# ---- Build annotation matrix: r with significance stars ------------------
annot_mat = np.empty((n_genes, n_genes), dtype=object)
for i in range(n_genes):
    for j in range(n_genes):
        r = r_mat.iloc[i, j]
        p = p_mat.iloc[i, j]
        if i == j:
            star_str = ""
        else:
            star_str = "**" if p < 0.01 else ("*" if p < 0.05 else "")
        annot_mat[i, j] = f"{r:.2f}{star_str}"

# ---- Save source data ----------------------------------------------------
src_path = os.path.join(SRC_DIR, "fig6b_hub_gene_correlation_source.csv")
source_out = r_mat.copy().reset_index().rename(columns={"index": "Gene"})
source_out.to_csv(src_path, index=False)
p_out_path = os.path.join(SRC_DIR, "fig6b_hub_gene_correlation_pvalues.csv")
p_out = p_mat.copy().reset_index().rename(columns={"index": "Gene"})
p_out.to_csv(p_out_path, index=False)
print(f"Saved correlation matrix source CSV -> {src_path}")
print(f"Saved p-values source CSV -> {p_out_path}")

# ---- Summary stats for caption ------------------------------------------
triu_mask = np.triu(np.ones((n_genes, n_genes), dtype=bool), k=1)
upper_r_vals = r_mat.values[triu_mask]
count_r_ge_085 = int(np.sum(np.abs(upper_r_vals) >= 0.85))
total_pairs = int(n_genes * (n_genes - 1) / 2)
count_p_lt_05 = int(np.sum(p_mat.values[triu_mask] < 0.05))
count_p_lt_01 = int(np.sum(p_mat.values[triu_mask] < 0.01))
median_r = float(np.median(upper_r_vals))
min_r = float(np.min(upper_r_vals))
max_r = float(np.max(upper_r_vals))

print("=" * 70)
print("11-HUB PEARSON CORRELATION SUMMARY (GSE62928 N=8)")
print("=" * 70)
print(f"  Total unique hub-gene pairs: {total_pairs}")
print(f"  Pairs with |r| >= 0.85:      {count_r_ge_085}")
print(f"  Pairs with P < 0.05:          {count_p_lt_05}")
print(f"  Pairs with P < 0.01:          {count_p_lt_01}")
print(f"  Median pairwise r:            {median_r:.3f}")
print(f"  Min / Max pairwise r:         {min_r:.3f} / {max_r:.3f}")
print("=" * 70)

# ---- Plot the heatmap ----------------------------------------------------
fig, ax = plt.subplots(figsize=(11.5, 9.5), dpi=300)

cbar_kws = {
    "label": "Pearson Correlation (r)\n(* P < 0.05, ** P < 0.01, unadjusted, N=8)",
    "shrink": 0.82
}

sns.heatmap(
    r_mat,
    annot=annot_mat,
    fmt="",
    cmap="vlag",
    center=0,
    vmin=-1.0,
    vmax=1.0,
    linewidths=0.7,
    linecolor="#E2E8F0",
    square=True,
    cbar_kws=cbar_kws,
    ax=ax,
    annot_kws={"size": 8.5, "weight": "normal"}
)

ax.set_title(
    "B. Pairwise Pearson Correlation Among 11 Consensus Hub Genes\n"
    f"(GSE62928 Discovery Tissue, N = 8: 4 EPS vs 4 Control; "
    f"P < 0.05 = {count_p_lt_05}/{total_pairs} pairs, P < 0.01 = {count_p_lt_01}/{total_pairs})",
    fontsize=12, fontweight="bold", color="#0F172A", pad=14
)
ax.set_xlabel("11 Consensus Pro-Fibrotic Hub Genes", fontsize=10.5, fontweight="bold", labelpad=8)
ax.set_ylabel("11 Consensus Pro-Fibrotic Hub Genes", fontsize=10.5, fontweight="bold", labelpad=8)
ax.tick_params(axis='x', labelsize=10)
plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
ax.tick_params(axis='y', labelsize=10, rotation=0)

# Footer callout (matches pipeline readability style)
stats_callout = (
    f"Median r = {median_r:.2f} | Range [{min_r:.2f}, {max_r:.2f}] | "
    f"|r| >= 0.85: {count_r_ge_085} pairs (self-pairs excluded)"
)
fig.text(
    0.5, 0.008, stats_callout,
    ha="center", va="bottom",
    fontsize=8.5, style="italic", color="#475569"
)

plt.tight_layout(rect=[0, 0.03, 1, 1])

# ---- Save PNG + PDF ------------------------------------------------------
png_path = os.path.join(FIG_DIR, "Hub_05_hub_gene_pearson_correlation_heatmap.png")
pdf_path = os.path.join(FIG_DIR, "Hub_05_hub_gene_pearson_correlation_heatmap.pdf")
fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"\nSaved PNG -> {png_path}")
print(f"Saved PDF -> {pdf_path}")
print("\nDone.")
