"""
Figure 7: Identifying Molecular Subtypes in Long-Term Peritoneal Dialysis Effluent
Replicating the exact figure layout and visual style of Sun et al. (Electron J Biotechnol 2026, Fig. 7).
Panels:
  A: Consensus clustering CDF curve (left), Relative change in area under CDF curve (middle), and Consensus matrix k=2 (right)
  B: Proportion of ambiguous clustering (PAC vs k, left) and PCA plot of LPD samples separated into C1 and C2 subtypes (right)
  C: Expression levels of diagnostic genes in C1 vs C2 molecular subtypes
  D: Differences in relative immune cell fractions between C1 and C2 molecular subtypes
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
from scipy.stats import mannwhitneyu

# Directory setup
OUT_DIR = os.path.join("results", "paper_style")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load Data
cdf_df = pd.read_csv(os.path.join(OUT_DIR, "Fig7A_consensus_CDF_source.csv"))
metrics_df = pd.read_csv(os.path.join(OUT_DIR, "Fig7A_consensus_metrics.csv"))
matrix_df = pd.read_csv(os.path.join(OUT_DIR, "Fig7A_best_k_consensus_matrix.csv"))
meta_df = pd.read_csv(os.path.join("results", "tables", "GSE125498_sample_metadata.csv"))
expr_df = pd.read_csv(os.path.join("results", "tables", "GSE125498_full_expression_matrix.csv"))
cell_scores = pd.read_csv(os.path.join(OUT_DIR, "Fig5E_saved_ssgsea_sample_scores.csv"))

# LPD discovery samples
lpd_samples = meta_df[meta_df["group"] == "LPD"]["sample_id"].tolist()
expr_lpd = expr_df.set_index("Gene")[lpd_samples]

# Reconstruct consensus matrix at k=2
piv_mat = matrix_df.pivot(index="Var1", columns="Var2", values="Freq")
samples_mat = piv_mat.index.tolist()
dist_mat = 1.0 - piv_mat.values
dist_mat = (dist_mat + dist_mat.T) / 2.0
np.fill_diagonal(dist_mat, 0.0)

# Hierarchical clustering on consensus distance
row_link = linkage(squareform(dist_mat), method="complete")
clusters_k2 = cut_tree(row_link, n_clusters=2).flatten() + 1
cluster_map = dict(zip(samples_mat, clusters_k2))

# Save cluster assignments at k=2
cluster_df = pd.DataFrame({"sample_id": samples_mat, "cluster": [cluster_map[s] for s in samples_mat]})
cluster_df.to_csv(os.path.join(OUT_DIR, "Fig7_cluster_assignments_k2.csv"), index=False)

# Order matrix by dendrogram
row_order = dendrogram(row_link, no_plot=True)["leaves"]
ordered_samples = [samples_mat[i] for i in row_order]
consensus_ordered = piv_mat.loc[ordered_samples, ordered_samples].values

# Diagnostic panel genes
panel_genes = ["TNFSF15", "ADAM19", "FLT3LG", "NAT6", "ECM1"]

# Subtype colors matching Sun et al. Fig 7
c_c1 = "#F43F5E"  # Rose / Pink for C1
c_c2 = "#0284C7"  # Cyan / Blue for C2

# Figure Setup (14.0 x 17.5 inches)
fig = plt.figure(figsize=(14.0, 17.5), facecolor="#FFFFFF")
gs = gridspec.GridSpec(3, 1, height_ratios=[1.1, 1.0, 1.15], hspace=0.34)

# ---------------------------------------------------------
# ROW 1: PANEL A (CDF left, Delta Area middle, Consensus Matrix right)
# ---------------------------------------------------------
gs_a = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], width_ratios=[1.0, 1.0, 1.15], wspace=0.28)
ax_a1 = fig.add_subplot(gs_a[0])
ax_a2 = fig.add_subplot(gs_a[1])
ax_a3 = fig.add_subplot(gs_a[2])

# A1: Consensus CDF
k_colors = {2: "#DC2626", 3: "#0284C7", 4: "#16A34A", 5: "#9333EA", 6: "#EA580C"}
for k_val in sorted(cdf_df["k"].unique()):
    sub_k = cdf_df[cdf_df["k"] == k_val].sort_values("consensus_value")
    col = k_colors.get(k_val, "#64748B")
    ax_a1.plot(sub_k["consensus_value"], sub_k["cumulative_fraction"], color=col, lw=1.6, label=f"k = {k_val}")

ax_a1.set_xlabel("Consensus Index", fontsize=9.5)
ax_a1.set_ylabel("CDF", fontsize=9.5)
ax_a1.set_title("Consensus CDF", fontsize=10.5, fontweight="bold", pad=8)
ax_a1.legend(loc="lower right", frameon=False, fontsize=8)
ax_a1.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
ax_a1.text(-0.20, 1.08, "A", transform=ax_a1.transAxes, fontsize=17, fontweight="bold", va="top")

# A2: Delta Area
k_vals = metrics_df["k"].values[1:]
delta_vals = metrics_df["Delta_area"].values[1:]
ax_a2.plot(k_vals, delta_vals, marker="o", color="#2563EB", lw=1.8, markersize=5.5)
ax_a2.set_xlabel("k", fontsize=9.5)
ax_a2.set_ylabel("Relative Change in Area under CDF", fontsize=9.5)
ax_a2.set_title("Delta Area", fontsize=10.5, fontweight="bold", pad=8)
ax_a2.set_xticks(k_vals)
ax_a2.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")

# A3: Consensus Matrix at k=2
gs_a3 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_a[2], height_ratios=[0.14, 0.05, 0.81], hspace=0.03)
ax_dend = fig.add_subplot(gs_a3[0])
dendrogram(row_link, orientation="top", ax=ax_dend, link_color_func=lambda k: "#475569", no_labels=True)
ax_dend.axis("off")

# Cluster annotation bar
ax_bar = fig.add_subplot(gs_a3[1])
clust_ordered = [cluster_map[s] for s in ordered_samples]
bar_arr = np.array([[1 if c == 1 else 0 for c in clust_ordered]])
cmap_c = plt.matplotlib.colors.ListedColormap([c_c2, c_c1])
ax_bar.imshow(bar_arr, aspect="auto", cmap=cmap_c)
ax_bar.set_xticks([])
ax_bar.set_yticks([])
ax_bar.spines[:].set_visible(False)

# Consensus matrix
ax_mat = fig.add_subplot(gs_a3[2])
im_mat = ax_mat.imshow(consensus_ordered, aspect="auto", cmap="Blues", vmin=0.0, vmax=1.0)
ax_mat.set_xticks([])
ax_mat.set_yticks([])
ax_mat.spines[:].set_visible(False)
ax_mat.set_title("Consensus Matrix k = 2", fontsize=10, fontweight="bold", pad=6)

# Colorbar for matrix
cb_mat = plt.colorbar(im_mat, ax=ax_mat, fraction=0.046, pad=0.04)
cb_mat.set_label("Consensus Value", fontsize=8)

# ---------------------------------------------------------
# ROW 2: PANEL B (PAC left, PCA right) & PANEL C (Expression by subtype)
# ---------------------------------------------------------
gs_bc = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], width_ratios=[0.85, 1.05, 1.1], wspace=0.28)
ax_b1 = fig.add_subplot(gs_bc[0])
ax_b2 = fig.add_subplot(gs_bc[1])
ax_c = fig.add_subplot(gs_bc[2])

# B1: PAC vs k
pacs = metrics_df["PAC_0_1"].values
ks = metrics_df["k"].values
ax_b1.plot(ks, pacs, marker="o", color="#0284C7", lw=1.8, markersize=5.5)
ax_b1.set_xlabel("Cluster Number K", fontsize=9.5)
ax_b1.set_ylabel("Proportion of Ambiguous Clustering", fontsize=9.0)
ax_b1.set_title("Ambiguous Clustering (PAC)", fontsize=10.5, fontweight="bold", pad=8)
ax_b1.set_xticks(ks)
ax_b1.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
ax_b1.text(-0.22, 1.08, "B", transform=ax_b1.transAxes, fontsize=17, fontweight="bold", va="top")

# B2: PCA of LPD samples
from sklearn.decomposition import PCA
pca_lpd = PCA(n_components=2, random_state=42)
lpd_matrix_scaled = (expr_lpd.loc[panel_genes].T - expr_lpd.loc[panel_genes].T.mean()) / expr_lpd.loc[panel_genes].T.std()
pca_coords_lpd = pca_lpd.fit_transform(lpd_matrix_scaled)
var1 = pca_lpd.explained_variance_ratio_[0] * 100
var2 = pca_lpd.explained_variance_ratio_[1] * 100

for sid, (px, py) in zip(lpd_samples, pca_coords_lpd):
    cl = cluster_map.get(sid, 1)
    col = c_c1 if cl == 1 else c_c2
    marker = "o" if cl == 1 else "^"
    ax_b2.scatter(px, py, color=col, marker=marker, s=38, alpha=0.85, edgecolors="#1E293B", lw=0.6)

# Confidence Ellipses
for cl, col in [(1, c_c1), (2, c_c2)]:
    pts = pca_coords_lpd[[cluster_map.get(s, 1) == cl for s in lpd_samples]]
    if len(pts) >= 3:
        center = np.mean(pts, axis=0)
        cov = np.cov(pts, rowvar=False)
        w, v = np.linalg.eigh(cov)
        angle = np.degrees(np.arctan2(v[1, 0], v[0, 0]))
        ell = patches.Ellipse(xy=center, width=2.2 * np.sqrt(max(w[0], 0.1)),
                              height=2.2 * np.sqrt(max(w[1], 0.1)), angle=angle,
                              facecolor=col, alpha=0.22, edgecolor=col, lw=1.4)
        ax_b2.add_patch(ell)

ax_b2.set_xlabel(f"PCA1 ({var1:.1f}%)", fontsize=9.5)
ax_b2.set_ylabel(f"PCA2 ({var2:.1f}%)", fontsize=9.5)
ax_b2.set_title("PCA of LPD Subtypes", fontsize=10.5, fontweight="bold", pad=8)
ax_b2.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")

# Legend for B2
h_c1 = ax_b2.scatter([], [], color=c_c1, marker="o", label="C1")
h_c2 = ax_b2.scatter([], [], color=c_c2, marker="^", label="C2")
ax_b2.legend(handles=[h_c1, h_c2], title="Cluster", loc="upper right", frameon=False, fontsize=8, title_fontsize=8.5)

# C: Panel Gene Expression by Subtype
for i, g in enumerate(panel_genes):
    v_c1 = [expr_lpd.loc[g, s] for s in lpd_samples if cluster_map.get(s, 1) == 1]
    v_c2 = [expr_lpd.loc[g, s] for s in lpd_samples if cluster_map.get(s, 1) == 2]
    
    pos_1 = i * 2.2 + 0.8
    pos_2 = i * 2.2 + 1.5
    
    ax_c.boxplot(v_c1, positions=[pos_1], widths=0.55, patch_artist=True,
                 boxprops=dict(facecolor=c_c1, color=c_c1, alpha=0.75),
                 medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    ax_c.boxplot(v_c2, positions=[pos_2], widths=0.55, patch_artist=True,
                 boxprops=dict(facecolor=c_c2, color=c_c2, alpha=0.75),
                 medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    
    ax_c.scatter(np.random.normal(pos_1, 0.05, len(v_c1)), v_c1, color=c_c1, s=8, alpha=0.7)
    ax_c.scatter(np.random.normal(pos_2, 0.05, len(v_c2)), v_c2, color=c_c2, s=8, alpha=0.7)
    
    # Wilcoxon test
    stat, pval = mannwhitneyu(v_c1, v_c2, alternative="two-sided")
    y_top = max(max(v_c1), max(v_c2)) + 0.25
    ast = "****" if pval < 0.0001 else ("***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns")))
    ax_c.text(i * 2.2 + 1.15, y_top, ast, ha="center", va="bottom", fontsize=8.5, fontweight="bold")

ax_c.set_xticks([i * 2.2 + 1.15 for i in range(len(panel_genes))])
ax_c.set_xticklabels(panel_genes, fontsize=8.5, fontweight="bold")
ax_c.set_ylabel("Normalized Expression", fontsize=9.5)
ax_c.set_title("Diagnostic Biomarkers in Subtypes", fontsize=10.5, fontweight="bold", pad=8)
ax_c.grid(True, linestyle=":", alpha=0.4, axis="y")

h1_patch = patches.Patch(facecolor=c_c1, label="C1")
h2_patch = patches.Patch(facecolor=c_c2, label="C2")
ax_c.legend(handles=[h1_patch, h2_patch], title="Cluster", loc="upper right", frameon=False, fontsize=8, title_fontsize=8.5)
ax_c.text(-0.16, 1.08, "C", transform=ax_c.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# ROW 3: PANEL D (Immune Differences between Subtypes)
# ---------------------------------------------------------
ax_d = fig.add_subplot(gs[2])

cell_types = ['CD4 T cell', 'CD8 T cell', 'NKT cell', 'NK cell', 'T cell (other)', 'B cell',
              'Plasma cell', 'Monocyte', 'Macrophage', 'Dendritic cell', 'Neutrophil', 'Eosinophil', 'Mast cell']

# Filter immune scores for LPD samples
cell_lpd = cell_scores[cell_scores["sample_id"].isin(lpd_samples)].copy()
cell_lpd["cluster"] = [cluster_map.get(s, 1) for s in cell_lpd["sample_id"]]

for i, ct in enumerate(cell_types):
    v_1 = cell_lpd[cell_lpd["cluster"] == 1][ct].values
    v_2 = cell_lpd[cell_lpd["cluster"] == 2][ct].values
    
    pos_1 = i * 2.2 + 0.8
    pos_2 = i * 2.2 + 1.5
    
    ax_d.boxplot(v_1, positions=[pos_1], widths=0.55, patch_artist=True,
                 boxprops=dict(facecolor=c_c1, color=c_c1, alpha=0.75),
                 medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    ax_d.boxplot(v_2, positions=[pos_2], widths=0.55, patch_artist=True,
                 boxprops=dict(facecolor=c_c2, color=c_c2, alpha=0.75),
                 medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    
    ax_d.scatter(np.random.normal(pos_1, 0.05, len(v_1)), v_1, color=c_c1, s=8, alpha=0.6)
    ax_d.scatter(np.random.normal(pos_2, 0.05, len(v_2)), v_2, color=c_c2, s=8, alpha=0.6)
    
    stat, pval = mannwhitneyu(v_1, v_2, alternative="two-sided")
    y_top = max(v_1.max(), v_2.max()) + 0.04
    if pval < 0.01:
        ax_d.text(i * 2.2 + 1.15, y_top, "**", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    elif pval < 0.05:
        ax_d.text(i * 2.2 + 1.15, y_top, "*", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

ax_d.set_xticks([i * 2.2 + 1.15 for i in range(len(cell_types))])
ax_d.set_xticklabels(cell_types, rotation=40, ha="right", fontsize=8.5, fontweight="medium")
ax_d.set_ylabel("ssGSEA Score", fontsize=9.5)
ax_d.set_title("Difference in Immune Cell Fractions between C1 and C2 Subtypes", fontsize=11, fontweight="bold", pad=8)
ax_d.grid(True, linestyle=":", alpha=0.4, axis="y")
ax_d.legend(handles=[h1_patch, h2_patch], title="Cluster", loc="upper right", frameon=False, fontsize=8.5, title_fontsize=9)
ax_d.text(-0.07, 1.08, "D", transform=ax_d.transAxes, fontsize=17, fontweight="bold", va="top")

# Caption at the very bottom
fig.text(0.5, 0.012,
         "Fig. 7. Identifying LPD-associated molecular subtypes. (A) LPD patients grouped into two clusters according to consensus clustering matrix (k = 2),\n"
         "along with CDF and delta area curves. (B) Proportion of ambiguous clustering (PAC) and PCA of LPD subtypes. (C) Expression levels of diagnostic genes\n"
         "in C1 vs C2 subtypes. (D) Difference in relative immune cell fractions between C1 and C2 subtypes. *p < 0.05, **p < 0.01, ***p < 0.001.",
         ha="center", va="bottom", fontsize=9.2, color="#0F172A", linespacing=1.35)

plt.subplots_adjust(bottom=0.065, top=0.96)

# Save high-res PNG and PDF
png_path = os.path.join(OUT_DIR, "Fig7_molecular_subtypes.png")
pdf_path = os.path.join(OUT_DIR, "Fig7_molecular_subtypes.pdf")
fig_png = os.path.join(OUT_DIR, "Fig7.png")
fig_pdf = os.path.join(OUT_DIR, "Fig7.pdf")

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
plt.savefig(fig_png, dpi=300, bbox_inches="tight")
plt.savefig(fig_pdf, bbox_inches="tight")
plt.close()

print(f"Successfully generated Figure 7: {png_path} and {pdf_path}")
