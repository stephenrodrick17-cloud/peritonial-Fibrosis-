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

# 2. Load GSE62928 Pro-Fibrotic Up-Regulated DEGs (P < 0.05, log2FC >= 0.585 [1.5-fold])
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

degs_profibrotic = set(df_best_p[(df_best_p["P.Value"] < 0.05) & (df_best_p["logFC"] >= 0.585)]["Gene_clean"])
n_degs = len(degs_profibrotic)

# Filtered Pro-Fibrotic Upregulated ECM-DEGs (81 Genes)
df_81 = pd.read_csv("convergent_81_ECM_DEGs.csv")
genes_81 = set(df_81["Gene_Symbol"].astype(str).str.strip().str.upper())
overlap_81 = len(genes_81)

print(f"\nCurated Pro-Fibrotic ECM-DEG Cohort (Primary Pipeline Filter):")
print(f"   - GSE62928 Pro-Fibrotic DEGs (log2FC >= 0.585 [1.5-fold], P < 0.05): {n_degs}")
print(f"   - Human Matrisome Genes: {n_ecm}")
print(f"   - Overlap: {overlap_81} ECM-DEGs")

# ==============================================================================
# 2-WAY VENN DIAGRAM: 534 PRO-FIBROTIC DEGs ∩ 1,027 MATRISOME GENES
# ==============================================================================
fig, ax = plt.subplots(figsize=(9, 7), facecolor="#F8FAFC")

only_deg_81 = n_degs - overlap_81  # 534 - 81 = 453
only_ecm_81 = n_ecm - overlap_81   # 1027 - 81 = 946

v2 = venn2(
    subsets=(only_deg_81, only_ecm_81, overlap_81),
    set_labels=(f"GSE62928 Pro-Fibrotic DEGs\n(Peritoneal Fibrosis, n={n_degs})", f"Human In Silico Matrisome\n(Curated ECM, n={n_ecm})"),
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

venn2_circles(subsets=(only_deg_81, only_ecm_81, overlap_81), linestyle="solid", linewidth=2.0, color="#0F172A", ax=ax)

plt.title("Discovery Intersection: GSE62928 Pro-Fibrotic DEGs ∩ Human Matrisome\nIsolation of 81 Up-Regulated ECM-DEGs",
          fontsize=13, fontweight="bold", pad=20, color="#0F172A")

plt.annotate(
    "81 Verified Pro-Fibrotic ECM-DEGs\n(Prioritized for WGCNA Convergence)",
    xy=(0, 0), xytext=(0, -0.45),
    arrowprops=dict(facecolor="#B91C1C", shrink=0.08, width=2, headwidth=8),
    ha="center", fontsize=11, fontweight="bold", color="#B91C1C",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="#FEF2F2", edgecolor="#B91C1C", linewidth=1.5)
)

plt.tight_layout()
venn_file2 = "results/figures/venn_gse62928_pro_fibrotic_ecm_81.png"
plt.savefig(venn_file2, dpi=300)
plt.close()
print(f"Saved 2-Way Venn: {venn_file2}")

print("\n" + "=" * 70)
print("VENN DIAGRAM GENERATION COMPLETED SUCCESSFULLY!")
print("=" * 70)
