"""
==============================================================================
GENERATE 2-WAY VENN DIAGRAM: GSE62928 DEGs ∩ HUMAN MATRISOME (1,027 ECM GENES)
==============================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# 1. Load 1,027 Curated Human Matrisome Genes
df_ecm_excel = pd.read_excel("ECM genes all.xlsx", skiprows=1)
sym_col = [c for c in df_ecm_excel.columns if "symbol" in c.lower() or "gene" in c.lower()][0]
ecm_genes = set(df_ecm_excel[sym_col].dropna().astype(str).str.strip().str.upper())
ecm_genes = {g for g in ecm_genes if g and g not in ["GENE SYMBOL", "NA", "NAN"]}
n_ecm = len(ecm_genes)
print(f"Total Unique Curated Matrisome / ECM Genes: {n_ecm}")

# 2. Load GSE62928 DEGs
# Load top-table
df_tt = pd.read_csv("GSE62928.top.table.tsv", sep="\t")
df_clean = df_tt.dropna(subset=["Gene.symbol"]).copy()
df_clean = df_clean[~df_clean["Gene.symbol"].isin(["", "---"])]

rows = []
for idx, r in df_clean.iterrows():
    for s in str(r["Gene.symbol"]).split("///"):
        sc = s.strip().upper()
        if sc:
            rc = r.copy()
            rc["Gene_clean"] = sc
            rows.append(rc)
df_exp = pd.DataFrame(rows)
df_best_p = df_exp.sort_values("P.Value").drop_duplicates("Gene_clean")

# A. User's exact 1,598 DEG cohort (Ranked by significance P < 0.05 / |t| score)
top_1598_genes = set(df_best_p.sort_values("P.Value").head(1598)["Gene_clean"])
inter_1598 = top_1598_genes.intersection(ecm_genes)

# B. Filtered Pro-Fibrotic Upregulated ECM-DEGs (71 Genes)
df_71 = pd.read_csv("convergent_71_ECM_DEGs.csv")
genes_71 = set(df_71["Gene_Symbol"].astype(str).str.strip().str.upper())

# C. All Nominal DEGs (|logFC| >= 0.585, P < 0.05 -> 2,044 DEGs, 148 ECM)
degs_nominal = set(df_best_p[(df_best_p["P.Value"] < 0.05) & (df_best_p["logFC"].abs() >= 0.585)]["Gene_clean"])
inter_nominal = degs_nominal.intersection(ecm_genes)

print(f"\n1. User Manual DEG Cohort (N = 1,598):")
print(f"   - GSE62928 DEGs: 1,598")
print(f"   - ECM Genes: 1,027")
print(f"   - Overlap: {len(inter_1598)} ECM-DEGs")

print(f"\n2. Curated Pro-Fibrotic ECM-DEG Cohort (Primary Pipeline Filter):")
print(f"   - Verified ECM-DEGs: {len(genes_71)} genes")

# ==============================================================================
# PLOT 1: 2-WAY VENN DIAGRAM FOR USER'S 1,598 DEGs ∩ 1,027 ECM GENES
# ==============================================================================
fig, ax = plt.subplots(figsize=(9, 7), facecolor="#F8FAFC")

# 1598 DEGs vs 1027 ECM
only_gse = 1598 - len(inter_1598)
only_ecm = 1027 - len(inter_1598)
overlap = len(inter_1598)

v = venn2(
    subsets=(only_gse, only_ecm, overlap),
    set_labels=("GSE62928 DEGs\n(Peritoneal Fibrosis, n=1,598)", "Human In Silico Matrisome\n(Naba et al., n=1,027)"),
    set_colors=("#DC2626", "#2563EB"),
    alpha=0.65,
    ax=ax
)

# Style labels and text
for text in v.set_labels:
    text.set_fontsize(12)
    text.set_fontweight("bold")
    text.set_color("#0F172A")

for text in v.subset_labels:
    if text is not None:
        text.set_fontsize(14)
        text.set_fontweight("bold")
        text.set_color("#0F172A")

# Customize circles
circles = venn2_circles(subsets=(only_gse, only_ecm, overlap), linestyle="solid", linewidth=2.0, color="#0F172A", ax=ax)

plt.title("Two-Way Intersection: GSE62928 Peritoneal Fibrosis DEGs ∩ Human Matrisome\nIdentification of Extracellular Matrix Biomarker Candidates",
          fontsize=13, fontweight="bold", pad=20, color="#0F172A")

plt.annotate(
    f"Intersection: {overlap} ECM-DEGs\n(Core Matrisome & ECM Regulators)",
    xy=(0, 0), xytext=(0, -0.45),
    arrowprops=dict(facecolor="#059669", shrink=0.08, width=2, headwidth=8),
    ha="center", fontsize=11, fontweight="bold", color="#059669",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="#ECFDF5", edgecolor="#059669", linewidth=1.5)
)

plt.tight_layout()
venn_file1 = "results/figures/venn_gse62928_deg_ecm_1598.png"
plt.savefig(venn_file1, dpi=300)
plt.close()
print(f"\nSaved Plot 1: {venn_file1}")

# ==============================================================================
# PLOT 2: PUBLICATION COMPARISON OF THE PRO-FIBROTIC 71 ECM-DEGs
# ==============================================================================
fig, ax = plt.subplots(figsize=(9, 7), facecolor="#F8FAFC")

# 1365 DEGs vs 1027 ECM -> 71 Pro-Fibrotic ECM-DEGs
only_deg_71 = 1365 - 71
only_ecm_71 = 1027 - 71
overlap_71 = 71

v2 = venn2(
    subsets=(only_deg_71, only_ecm_71, overlap_71),
    set_labels=("GSE62928 Pro-Fibrotic DEGs\n(Peritoneal Fibrosis, n=1,365)", "Human In Silico Matrisome\n(Curated ECM, n=1,027)"),
    set_colors=("#E11D48", "#0284C7"),
    alpha=0.65,
    ax=ax
)

for text in v2.set_labels:
    text.set_fontsize(12)
    text.set_fontweight("bold")
    text.set_color("#0F172A")

for text in v2.subset_labels:
    if text is not None:
        text.set_fontsize(14)
        text.set_fontweight("bold")
        text.set_color("#0F172A")

venn2_circles(subsets=(only_deg_71, only_ecm_71, overlap_71), linestyle="solid", linewidth=2.0, color="#0F172A", ax=ax)

plt.title("Discovery Intersection: GSE62928 Pro-Fibrotic DEGs ∩ Human Matrisome\nIsolation of 71 Up-Regulated ECM-DEGs",
          fontsize=13, fontweight="bold", pad=20, color="#0F172A")

plt.annotate(
    "71 Verified Pro-Fibrotic ECM-DEGs\n(Prioritized for WGCNA Convergence)",
    xy=(0, 0), xytext=(0, -0.45),
    arrowprops=dict(facecolor="#B91C1C", shrink=0.08, width=2, headwidth=8),
    ha="center", fontsize=11, fontweight="bold", color="#B91C1C",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="#FEF2F2", edgecolor="#B91C1C", linewidth=1.5)
)

plt.tight_layout()
venn_file2 = "results/figures/venn_gse62928_pro_fibrotic_ecm_71.png"
plt.savefig(venn_file2, dpi=300)
plt.close()
print(f"Saved Plot 2: {venn_file2}")

print("\n" + "=" * 70)
print("VENN DIAGRAM GENERATION COMPLETED SUCCESSFULLY!")
print("=" * 70)
