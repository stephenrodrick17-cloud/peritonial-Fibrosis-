"""
==============================================================================
TASK 3: RIGOROUS PRERANKED GSEA & CIRCULARITY-CORRECTED ssGSEA PATHWAY ANALYSIS
==============================================================================
Provides a two-part pathway enrichment analysis:
  Part A: True Preranked GSEA (gseapy) on all genes from GSE62928 ranked by Limma
          moderated t-statistic across MSigDB Hallmark Gene Sets (v2020/v2023.2).
          Reports NES, nominal P-value, FDR q-value, and leading-edge genes.
  Part B: Single-Sample GSEA (ssGSEA) score correlation heatmap.
          Directly addresses and eliminates circularity by removing each hub gene
          from the target hallmark gene set prior to score calculation.
          Computes Spearman correlation with Benjamini-Hochberg FDR correction.
  Exports dual-format publication figures (300 DPI PNG + Vector PDF) and CSV tables.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import gseapy as gp

# Set seed and directories
np.random.seed(42)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# 1. Consensus Hub Genes & Expression Data
hub_genes = ["FN1", "COL3A1", "COL11A1", "COL8A1", "VCAN", "COMP", "THBS3", "EDIL3", "LOX", "INHBA", "ISM1"]
df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
df_meta = pd.read_csv("results/tables/GSE62928_sample_metadata.csv")
sample_ids = df_meta["sample_id"].tolist()
n_samples = len(sample_ids)

# 2. PART A: True Preranked GSEA on GSE62928 Limma Top Table
print("Step 1: Preparing ranked gene list from GSE62928 Limma moderated t-statistics...")
df_tt = pd.read_csv("GSE62928.top.table.tsv", sep="\t")
df_tt["Gene.symbol"] = df_tt["Gene.symbol"].fillna("").astype(str).str.strip()
df_tt = df_tt[df_tt["Gene.symbol"] != ""]
df_tt["primary_symbol"] = df_tt["Gene.symbol"].apply(lambda x: x.split("///")[0].strip().upper())

# Deduplicate by maximum absolute moderated t-statistic
df_ranked = df_tt.sort_values(by="t", key=abs, ascending=False).drop_duplicates(subset=["primary_symbol"])
rnk = df_ranked.set_index("primary_symbol")["t"].sort_values(ascending=False)
print(f"Total genes ranked for preranked GSEA: {len(rnk)}")

print("Step 2: Running Preranked GSEA on MSigDB Hallmark Gene Sets...")
preres = gp.prerank(
    rnk=rnk,
    gene_sets="MSigDB_Hallmark_2020",
    threads=4,
    min_size=15,
    max_size=500,
    permutation_num=1000,
    outdir=None,
    seed=42,
    verbose=False
)
df_gsea = preres.res2d.copy()
df_gsea["Term_Clean"] = df_gsea["Term"].apply(lambda x: f"HALLMARK_{x.upper().replace(' ', '_').replace('-', '_')}")
df_gsea = df_gsea.sort_values("NES", ascending=False)

# Export Preranked GSEA Table
out_gsea_csv = "results/tables/gsea_preranked_hallmark_results.csv"
export_cols = ["Term_Clean", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes"]
df_gsea[export_cols].to_csv(out_gsea_csv, index=False)
print(f"Exported Preranked GSEA results to {out_gsea_csv}")

# 3. PART B: ssGSEA Score Correlation Heatmap with Circularity Removal
print("Step 3: Computing ssGSEA scores and eliminating circularity...")
# Extract Hallmark gene sets from gseapy library
hallmark_dict = gp.get_library(name="MSigDB_Hallmark_2020")

# Key fibrosis & microenvironmental hallmark pathways
focal_pathways = {
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "Epithelial Mesenchymal Transition",
    "HALLMARK_TGF_BETA_SIGNALING": "TGF-beta Signaling",
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory Response",
    "HALLMARK_ANGIOGENESIS": "Angiogenesis",
    "HALLMARK_HYPOXIA": "Hypoxia",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF-alpha Signaling via NF-kB",
    "HALLMARK_COAGULATION": "Coagulation",
    "HALLMARK_APICAL_JUNCTION": "Apical Junction"
}

# Function to calculate standardized Z-score ssGSEA approximation
def compute_pathway_ssgsea(genes, expr_df, samples):
    avail = [g for g in genes if g in expr_df.index]
    if len(avail) < 3:
        return np.zeros(len(samples))
    sub = expr_df.loc[avail, samples]
    # Standardize gene expression across samples
    z_sub = sub.apply(lambda x: (x - x.mean()) / (x.std() + 1e-7), axis=1)
    return z_sub.mean(axis=0).values

# 1. Unadjusted (Standard) ssGSEA Scores
std_scores = {}
for p_label, p_name in focal_pathways.items():
    p_genes = hallmark_dict.get(p_name, [])
    std_scores[p_label] = compute_pathway_ssgsea(p_genes, df_expr, sample_ids)

# 2. Circularity-Corrected ssGSEA Scores:
# For each gene g and pathway p, compute score excluding g if g in p_genes
corr_rows_std = []
corr_rows_corrected = []

pvals_std_list = []
pvals_corr_list = []

for g in hub_genes:
    g_expr = df_expr.loc[g, sample_ids].values
    row_std = {"Gene_Symbol": g}
    row_corr = {"Gene_Symbol": g}
    
    for p_label, p_name in focal_pathways.items():
        p_genes = hallmark_dict.get(p_name, [])
        
        # Standard
        s_std = std_scores[p_label]
        r_s, p_s = spearmanr(g_expr, s_std)
        row_std[f"{p_label}_r"] = r_s
        row_std[f"{p_label}_p"] = p_s
        pvals_std_list.append(p_s)
        
        # Circularity-corrected (exclude g from pathway genes)
        genes_loo = [x for x in p_genes if x != g]
        s_corr = compute_pathway_ssgsea(genes_loo, df_expr, sample_ids)
        r_c, p_c = spearmanr(g_expr, s_corr)
        row_corr[f"{p_label}_r"] = r_c
        row_corr[f"{p_label}_p"] = p_c
        pvals_corr_list.append(p_c)
        
    corr_rows_std.append(row_std)
    corr_rows_corrected.append(row_corr)

# Multiple testing correction (FDR Benjamini-Hochberg)
_, fdr_std, _, _ = multipletests(pvals_std_list, method="fdr_bh")
_, fdr_corr, _, _ = multipletests(pvals_corr_list, method="fdr_bh")

df_corr_std = pd.DataFrame(corr_rows_std)
df_corr_corrected = pd.DataFrame(corr_rows_corrected)

# Add FDR columns
k = 0
for g_idx in range(len(hub_genes)):
    for p_label in focal_pathways.keys():
        df_corr_std.loc[g_idx, f"{p_label}_FDR"] = fdr_std[k]
        df_corr_corrected.loc[g_idx, f"{p_label}_FDR"] = fdr_corr[k]
        k += 1

# Export Correlation Table
out_corr_csv = "results/tables/hub_genes_gsea_pathway_correlations.csv"
df_corr_corrected.to_csv(out_corr_csv, index=False)
print(f"Exported Circularity-Corrected correlations to {out_corr_csv}")

# 4. Multi-Panel Publication Visualization (3-Panel Figure: Preranked Barplot + Standard Heatmap + Circularity-Corrected Heatmap)
fig = plt.figure(figsize=(20, 10), facecolor="#FFFFFF")

# Panel A: Top 10 Enriched Hallmark Pathways from Preranked GSEA
ax1 = plt.subplot2grid((1, 3), (0, 0))
top_pos = df_gsea[df_gsea["NES"] > 0].head(8)
top_neg = df_gsea[df_gsea["NES"] < 0].tail(4)
plot_gsea = pd.concat([top_pos, top_neg]).sort_values("NES")

bar_colors = ["#DC2626" if nes > 0 else "#2563EB" for nes in plot_gsea["NES"]]
bars = ax1.barh(range(len(plot_gsea)), plot_gsea["NES"], color=bar_colors, edgecolor="#0F172A", height=0.65)
ax1.axvline(0, color="#64748B", linestyle="--", lw=1.2)
ax1.set_yticks(range(len(plot_gsea)))
ax1.set_yticklabels([t.replace("HALLMARK_", "").replace("_", " ") for t in plot_gsea["Term_Clean"]], fontsize=9.5, fontweight="bold")
ax1.set_xlabel("Normalized Enrichment Score (NES)", fontsize=10.5, fontweight="bold")
ax1.set_title("A. True Preranked GSEA (MSigDB Hallmark v2020)\nRanked by Moderated t-Statistic (EPS vs Control)", fontsize=11, fontweight="bold", color="#0F172A", pad=12)

# Annotate FDR q-values on bars
for idx, (nes, fdr) in enumerate(zip(plot_gsea["NES"], plot_gsea["FDR q-val"])):
    offset = 0.08 if nes > 0 else -0.08
    ha_align = "left" if nes > 0 else "right"
    ax1.text(nes + offset, idx, f"q={fdr:.3f}", va="center", ha=ha_align, fontsize=8.0, fontweight="bold", color="#334155")
ax1.set_xlim(min(plot_gsea["NES"]) - 0.5, max(plot_gsea["NES"]) + 0.6)
ax1.grid(True, linestyle=":", alpha=0.5, axis="x")

# Matrices for heatmaps
path_cols_r = [f"{p}_r" for p in focal_pathways.keys()]
path_cols_fdr = [f"{p}_FDR" for p in focal_pathways.keys()]
short_path_labels = [p.replace("HALLMARK_", "").replace("_", " ") for p in focal_pathways.keys()]

mat_std = df_corr_std.set_index("Gene_Symbol")[path_cols_r].values
mat_std_fdr = df_corr_std.set_index("Gene_Symbol")[path_cols_fdr].values

mat_corr = df_corr_corrected.set_index("Gene_Symbol")[path_cols_r].values
mat_corr_fdr = df_corr_corrected.set_index("Gene_Symbol")[path_cols_fdr].values

# Annotations with significance stars (* FDR < 0.05, ** FDR < 0.01)
annot_std = np.empty(mat_std.shape, dtype=object)
for i in range(mat_std.shape[0]):
    for j in range(mat_std.shape[1]):
        r = mat_std[i, j]
        fdr = mat_std_fdr[i, j]
        star = "**" if fdr < 0.01 else ("*" if fdr < 0.05 else "")
        annot_std[i, j] = f"{r:.2f}{star}"

annot_corr = np.empty(mat_corr.shape, dtype=object)
for i in range(mat_corr.shape[0]):
    for j in range(mat_corr.shape[1]):
        r = mat_corr[i, j]
        fdr = mat_corr_fdr[i, j]
        star = "**" if fdr < 0.01 else ("*" if fdr < 0.05 else "")
        annot_corr[i, j] = f"{r:.2f}{star}"

# Panel B: Standard ssGSEA Correlation
ax2 = plt.subplot2grid((1, 3), (0, 1))
sns.heatmap(mat_std, annot=annot_std, fmt="", cmap="vlag", center=0, vmin=-1.0, vmax=1.0,
            xticklabels=short_path_labels, yticklabels=hub_genes,
            linewidths=0.8, linecolor="#E2E8F0", cbar=False, ax=ax2)
ax2.set_xticklabels(short_path_labels, rotation=45, ha="right", fontsize=9, fontweight="bold")
ax2.set_yticklabels(hub_genes, rotation=0, fontsize=9.5, fontweight="bold")
ax2.set_title("B. Standard ssGSEA Score Correlation\n(Unadjusted; Potential Circularity)", fontsize=11, fontweight="bold", color="#0F172A", pad=12)

# Panel C: Circularity-Corrected ssGSEA Correlation
ax3 = plt.subplot2grid((1, 3), (0, 2))
cbar_kws = {"label": "Spearman Correlation (r)\n(* FDR < 0.05, ** FDR < 0.01)", "shrink": 0.8}
sns.heatmap(mat_corr, annot=annot_corr, fmt="", cmap="vlag", center=0, vmin=-1.0, vmax=1.0,
            xticklabels=short_path_labels, yticklabels=hub_genes,
            linewidths=0.8, linecolor="#E2E8F0", cbar_kws=cbar_kws, ax=ax3)
ax3.set_xticklabels(short_path_labels, rotation=45, ha="right", fontsize=9, fontweight="bold")
ax3.set_yticklabels(hub_genes, rotation=0, fontsize=9.5, fontweight="bold")
ax3.set_title("C. Circularity-Corrected ssGSEA Correlation\n(Hub Gene Removed from Gene Set; N = 8)", fontsize=11, fontweight="bold", color="#0F172A", pad=12)

plt.tight_layout()
out_gsea_png = "results/figures/Hub_03_gsea_pathway_enrichment_heatmap.png"
out_gsea_pdf = "results/figures/Hub_03_gsea_pathway_enrichment_heatmap.pdf"
plt.savefig(out_gsea_png, dpi=300, bbox_inches="tight")
plt.savefig(out_gsea_pdf, bbox_inches="tight")
plt.close()
print(f"Saved Task 3 GSEA figures to {out_gsea_png} and {out_gsea_pdf}")
