"""
Generate publication-quality heatmap of the 46 Matrisome DEGs
Comparing SPD (n=20) vs LPD (n=13) in GSE125498.
Matches style of the reference heatmap (C in paper):
- Top annotation bar showing groups (SPD vs LPD)
- Rows: 46 Matrisome DEGs (33 Upregulated, 13 Downregulated in LPD vs SPD)
- Z-score normalized expression per gene
- Blue-White-Pink / Red-White-Blue diverging color palette
- Dendrogram for gene clustering
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram

# Paths
out_dir = "results/revision/wgcna_paper_figures"
os.makedirs(out_dir, exist_ok=True)

# 1. Load data
degs_df = pd.read_csv("results/tables/reconciled_matrisome_DEGs.csv")
meta_df = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
expr_df = pd.read_csv("results/tables/GSE125498_full_expression_matrix.csv")

print(f"Total DEGs loaded: {len(degs_df)}")
up_genes = degs_df[degs_df["DEG_Status"] == "Upregulated"]["Gene"].tolist()
down_genes = degs_df[degs_df["DEG_Status"] == "Downregulated"]["Gene"].tolist()
print(f"Upregulated (in LPD vs SPD): {len(up_genes)}")
print(f"Downregulated (in LPD vs SPD): {len(down_genes)}")

# 2. Extract expression values for the 46 genes
mat_df = expr_df[expr_df["Gene"].isin(degs_df["Gene"])].copy()
mat_df = mat_df.set_index("Gene")

# Sort samples by group: SPD first, then LPD
spd_samples = meta_df[meta_df["group"] == "SPD"]["sample_id"].tolist()
lpd_samples = meta_df[meta_df["group"] == "LPD"]["sample_id"].tolist()
ordered_samples = spd_samples + lpd_samples

mat_df = mat_df[ordered_samples]

# 3. Compute row Z-scores (per gene across samples)
z_scores = mat_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Hierarchical clustering of genes (rows)
row_linkage = linkage(z_scores, method='average', metric='correlation')

# Color palette: reference uses cyan/blue -> white -> hotpink/magenta
# or let's create a custom colormap matching the image:
# -4 (cyan/blue) to 0 (white) to +4 (hot pink/rose)
cmap = mcolors.LinearSegmentedColormap.from_list(
    "custom_diverging",
    ["#0099cc", "#56b4e9", "#ffffff", "#f596c8", "#e63385"]
)

# Colors for sample groups:
# In reference: blue = Healthy/Control, pink/magenta = Disease
# For our data: SPD (short peritoneal dialysis / baseline) = #0099cc, LPD (long peritoneal dialysis / advanced) = #e63385
group_colors = {"SPD": "#0099cc", "LPD": "#e63385"}
col_colors = [group_colors[meta_df[meta_df["sample_id"] == s]["group"].values[0]] for s in ordered_samples]

# Plot clustermap with row clustering only (keeping sample columns ordered by group as in reference)
g = sns.clustermap(
    z_scores,
    row_linkage=row_linkage,
    col_cluster=False,
    col_colors=col_colors,
    cmap=cmap,
    center=0,
    vmin=-3,
    vmax=3,
    figsize=(12, 10),
    dendrogram_ratio=(0.18, 0.05),
    cbar_pos=(0.02, 0.15, 0.03, 0.25),
    yticklabels=True,
    xticklabels=False,
    linewidths=0.2,
    linecolor="#e0e0e0"
)

# Clean labels
g.ax_heatmap.set_ylabel("")
g.ax_heatmap.set_xlabel(f"Samples (SPD: n={len(spd_samples)}, LPD: n={len(lpd_samples)})", fontsize=11, labelpad=10)
g.ax_heatmap.tick_params(axis='y', labelsize=8)

# Add Legend for Group annotation
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=group_colors["SPD"], edgecolor="none", label="SPD (Short PD, n=20)"),
    Patch(facecolor=group_colors["LPD"], edgecolor="none", label="LPD (Long PD, n=13)")
]
g.ax_col_colors.legend(
    handles=legend_elements,
    title="Group",
    loc="upper right",
    bbox_to_anchor=(1.25, 2.5),
    frameon=False,
    fontsize=9,
    title_fontsize=10
)

# Title
g.fig.suptitle(
    "Expression Heatmap of 46 Matrisome DEGs\n(33 Upregulated, 13 Downregulated in LPD vs SPD)",
    fontsize=14,
    fontweight="bold",
    y=0.98
)

# Save figure in multiple formats
for ext in ["png", "pdf", "svg"]:
    out_file = os.path.join(out_dir, f"heatmap_46_matrisome_DEGs.{ext}")
    g.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_file}")

plt.close()
print("Done!")
