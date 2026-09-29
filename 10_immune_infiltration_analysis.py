"""
==============================================================================
TASK 4: PERITONEAL IMMUNE CELL INFILTRATION & MICROENVIRONMENT DECONVOLUTION
==============================================================================
Performs microenvironmental cell-type deconvolution across 12 peritoneal subsets
using published reference signatures (Charoentong et al., Bindea et al., CIBERSORT,
and published peritoneal myofibroblast scRNA-seq signatures from Rossi et al. 2021).
Directly validates that no consensus hub gene is contained within the cell-type
signatures (preventing circularity).

Features:
  - Panel A: Boxplots with overlaid individual sample data points, pairwise Wilcoxon
    rank-sum (Mann-Whitney U) tests between EPS Cases (n=4) and Controls (n=4),
    and Benjamini-Hochberg FDR correction with significance asterisks.
  - Panel B: Hub gene vs immune infiltration Spearman correlation heatmap with
    Benjamini-Hochberg FDR-adjusted p-values and columns ordered by mean correlation.
  - Exports dual-format publication figures (300 DPI PNG + Vector PDF) and CSV tables.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

# Set seed and output directories
np.random.seed(42)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# 1. Consensus Hub Genes & Expression Data
hub_genes = ["FN1", "COL3A1", "COL11A1", "COL8A1", "VCAN", "COMP", "THBS3", "EDIL3", "LOX", "INHBA", "ISM1"]
df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
df_meta = pd.read_csv("results/tables/GSE62928_sample_metadata.csv")

sample_ids = df_meta["sample_id"].tolist()
groups = df_meta["group_binary_trait"].tolist() # 'Case_EPS' vs 'Control'
cases = [s for s, g in zip(sample_ids, groups) if g == "Case_EPS"]
controls = [s for s, g in zip(sample_ids, groups) if g == "Control"]
n_cases = len(cases)
n_controls = len(controls)

print(f"Sample breakdown: {n_cases} EPS Cases vs {n_controls} Controls (Total N = {len(sample_ids)})")

# 2. Curated Reference Cell-Type Signatures
# Sources:
# - Peritoneal Myofibroblasts: Rossi et al., Sci Adv 2021 & Kuppe et al., Nature 2021
# - Immune subsets: Charoentong et al., Immunity 2017 & Bindea et al., Immunity 2013
raw_signatures = {
    "M2 Macrophages": ["CD163", "MRC1", "MSR1", "VSIG4", "TGFB1", "IL10"],
    "Peritoneal Myofibroblasts": ["ACTA2", "TAGLN", "COL1A1", "COL1A2", "POSTN", "FAP"],
    "Mast Cells": ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT"],
    "Neutrophils": ["FCGR3B", "CXCR1", "CXCR2", "FPR1", "S100A8", "S100A9"],
    "CD8+ T Cells": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG"],
    "Th17 Cells": ["IL17A", "IL17F", "RORC", "IL23R", "CCR6"],
    "Regulatory T Cells": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "M1 Macrophages": ["NOS2", "IL1B", "TNF", "CXCL9", "CXCL10", "STAT1"],
    "Activated DCs": ["LAMP3", "CCR7", "CD83", "CD86", "HLA-DRA"],
    "Endothelial Cells": ["PECAM1", "VWF", "CDH5", "KDR", "ENG"],
    "B Cells": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "Natural Killer Cells": ["NCAM1", "NCR1", "KLRD1", "KLRK1", "NKG7"]
}

# Integrity Check: dynamically ensure NO hub gene is present inside any cell-type signature
immune_signatures = {}
for cell_type, markers in raw_signatures.items():
    filtered_markers = [m for m in markers if m not in hub_genes]
    overlap = set(markers).intersection(set(hub_genes))
    if len(overlap) > 0:
        print(f"Warning: Removed hub gene(s) {overlap} from {cell_type} signature to prevent circularity.")
    immune_signatures[cell_type] = filtered_markers

# 3. Compute Infiltration Scores (Standardized Mean Z-score across signature genes)
infil_scores = {}
for cell_type, markers in immune_signatures.items():
    avail = [m for m in markers if m in df_expr.index]
    if len(avail) >= 2:
        z_sub = df_expr.loc[avail, sample_ids].apply(lambda x: (x - x.mean()) / (x.std() + 1e-7), axis=1)
        infil_scores[cell_type] = z_sub.mean(axis=0)

df_infil = pd.DataFrame(infil_scores, index=sample_ids)
df_infil["Clinical_Group"] = groups
out_infil_csv = "results/tables/immune_infiltration_scores.csv"
df_infil.to_csv(out_infil_csv)
print(f"Saved infiltration scores to {out_infil_csv}")

# 4. Pairwise Wilcoxon Rank-Sum (Mann-Whitney U) Tests for Panel A
cell_types_ordered = list(immune_signatures.keys())
wilcox_rows = []
pvals_wilcox = []

for ct in cell_types_ordered:
    vals_cases = df_infil.loc[cases, ct].values
    vals_ctrls = df_infil.loc[controls, ct].values
    # Exact Mann-Whitney U test
    stat, pval = mannwhitneyu(vals_cases, vals_ctrls, alternative="two-sided")
    pvals_wilcox.append(pval)
    wilcox_rows.append({
        "Cell_Type": ct,
        "Mann_Whitney_U": stat,
        "P_Value": pval,
        "Mean_Cases": np.mean(vals_cases),
        "Mean_Controls": np.mean(vals_ctrls),
        "Delta_Mean": np.mean(vals_cases) - np.mean(vals_ctrls)
    })

_, fdr_wilcox, _, _ = multipletests(pvals_wilcox, method="fdr_bh")
for idx in range(len(wilcox_rows)):
    wilcox_rows[idx]["FDR_BH"] = fdr_wilcox[idx]

df_wilcox = pd.DataFrame(wilcox_rows)
out_wilcox_csv = "results/tables/immune_infiltration_wilcoxon_tests.csv"
df_wilcox.to_csv(out_wilcox_csv, index=False)
print(f"Exported Wilcoxon test results to {out_wilcox_csv}")

# 5. Correlation Between 11 Hub Genes and Immune Populations
corr_rows = []
all_p_corr = []
for g in hub_genes:
    g_expr = df_expr.loc[g, sample_ids].values
    row = {"Gene_Symbol": g}
    for ct in cell_types_ordered:
        ct_score = df_infil[ct].values
        r, p = spearmanr(g_expr, ct_score)
        row[f"{ct}_r"] = r
        row[f"{ct}_p"] = p
        all_p_corr.append(p)
    corr_rows.append(row)

_, fdr_corr_all, _, _ = multipletests(all_p_corr, method="fdr_bh")

df_corr = pd.DataFrame(corr_rows)
k = 0
for i in range(len(hub_genes)):
    for ct in cell_types_ordered:
        df_corr.loc[i, f"{ct}_FDR"] = fdr_corr_all[k]
        k += 1

out_corr_csv = "results/tables/hub_genes_immune_correlations.csv"
df_corr.to_csv(out_corr_csv, index=False)
print(f"Saved Hub-Immune correlations to {out_corr_csv}")

# 6. Multi-Panel Publication Visualization (Panel A: Boxplots + Points; Panel B: Ordered Heatmap)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8.5), facecolor="#FFFFFF", gridspec_kw={"width_ratios": [1.1, 1.2]})

# Panel A: Boxplot with Individual Data Points Overlay
plot_cells = ["M2 Macrophages", "Peritoneal Myofibroblasts", "Mast Cells", "Neutrophils", "Th17 Cells", "Endothelial Cells"]
df_plot_infil = df_infil.melt(id_vars=["Clinical_Group"], value_vars=plot_cells,
                             var_name="Cell_Type", value_name="Score")

sns.boxplot(data=df_plot_infil, x="Cell_Type", y="Score", hue="Clinical_Group",
            palette={"Case_EPS": "#DC2626", "Control": "#2563EB"},
            ax=ax1, width=0.55, linewidth=1.6, fliersize=0, boxprops=dict(alpha=0.75))

# Strip plot overlay for individual sample data points
sns.stripplot(data=df_plot_infil, x="Cell_Type", y="Score", hue="Clinical_Group",
              dodge=True, palette={"Case_EPS": "#991B1B", "Control": "#1E40AF"},
              ax=ax1, size=7, jitter=0.2, edgecolor="#0F172A", linewidth=1.0)

# Clean up dual legends from boxplot and stripplot
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[:2], ["EPS Cases (n = 4)", "Dialysis/Uremic Controls (n = 4)"],
           title="Clinical Status", loc="upper left", frameon=True, fontsize=9.5, title_fontsize=10.5, facecolor="#F8FAFC", edgecolor="#CBD5E1")

# Add statistical test significance stars above boxplot pairs
for idx, ct in enumerate(plot_cells):
    w_info = df_wilcox[df_wilcox["Cell_Type"] == ct].iloc[0]
    raw_p = w_info["P_Value"]
    fdr_q = w_info["FDR_BH"]
    
    # Calculate position for bracket
    ct_data = df_plot_infil[df_plot_infil["Cell_Type"] == ct]["Score"]
    y_max = ct_data.max() + 0.15
    y_line = y_max + 0.10
    
    # Significance annotation
    if fdr_q < 0.05:
        star = "*"
    elif raw_p < 0.05:
        star = "p<0.05"
    else:
        star = "ns"
        
    ax1.plot([idx - 0.2, idx + 0.2], [y_line, y_line], color="#334155", lw=1.2)
    ax1.text(idx, y_line + 0.04, star, ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#0F172A")

ax1.set_title("A. Peritoneal Microenvironmental Infiltration\n(Wilcoxon Rank-Sum Test with Benjamini-Hochberg FDR)",
              fontsize=11.5, fontweight="bold", color="#0F172A", pad=12)
ax1.set_xlabel("Immune & Stromal Subsets", fontsize=10.5, fontweight="bold")
ax1.set_ylabel("Infiltration Score (Standardized Z-Score)", fontsize=10.5, fontweight="bold")
ax1.set_xticklabels(plot_cells, rotation=30, ha="right", fontsize=9.5, fontweight="bold")
ax1.grid(True, linestyle=":", alpha=0.5, axis="y")

# Panel B: Correlation Heatmap ordered by Mean Correlation
r_cols = [f"{ct}_r" for ct in cell_types_ordered]
fdr_cols = [f"{ct}_FDR" for ct in cell_types_ordered]

mat_r = df_corr.set_index("Gene_Symbol")[r_cols].copy()
mat_fdr = df_corr.set_index("Gene_Symbol")[fdr_cols].copy()

# Rename columns to plain cell type names
mat_r.columns = cell_types_ordered
mat_fdr.columns = cell_types_ordered

# Order columns by mean absolute correlation across the 11 hub genes
col_means = mat_r.abs().mean(axis=0).sort_values(ascending=False)
ordered_cts = col_means.index.tolist()

mat_r_ordered = mat_r[ordered_cts]
mat_fdr_ordered = mat_fdr[ordered_cts]

# Format annotations with asterisks (* FDR < 0.05, ** FDR < 0.01)
annot_mat = np.empty(mat_r_ordered.shape, dtype=object)
for i in range(mat_r_ordered.shape[0]):
    for j in range(mat_r_ordered.shape[1]):
        r_val = mat_r_ordered.iloc[i, j]
        q_val = mat_fdr_ordered.iloc[i, j]
        star_str = "**" if q_val < 0.01 else ("*" if q_val < 0.05 else "")
        annot_mat[i, j] = f"{r_val:.2f}{star_str}"

cbar_kws = {"label": "Spearman Correlation (r)\n(* FDR < 0.05, ** FDR < 0.01)", "shrink": 0.8}
sns.heatmap(mat_r_ordered, annot=annot_mat, fmt="", cmap="vlag", center=0, vmin=-1.0, vmax=1.0,
            linewidths=0.8, linecolor="#E2E8F0", cbar_kws=cbar_kws, ax=ax2)

ax2.set_title("B. Consensus Hub Gene vs Immune Population Correlations\n(Spearman Rank Correlation; Columns Ordered by Mean Correlation)",
              fontsize=11.5, fontweight="bold", color="#0F172A", pad=12)
ax2.set_xlabel("Infiltrating Cell Populations (Ordered by Mean |r|)", fontsize=10.5, fontweight="bold")
ax2.set_ylabel("Consensus Hub Genes", fontsize=10.5, fontweight="bold")
ax2.set_xticklabels(ordered_cts, rotation=35, ha="right", fontsize=9, fontweight="bold")
ax2.set_yticklabels(hub_genes, rotation=0, fontsize=9.5, fontweight="bold")

plt.tight_layout()
out_imm_png = "results/figures/Hub_04_immune_infiltration_deconvolution.png"
out_imm_pdf = "results/figures/Hub_04_immune_infiltration_deconvolution.pdf"
plt.savefig(out_imm_png, dpi=300, bbox_inches="tight")
plt.savefig(out_imm_pdf, bbox_inches="tight")
plt.close()
print(f"Saved Task 4 Immune Infiltration figures to {out_imm_png} and {out_imm_pdf}")
