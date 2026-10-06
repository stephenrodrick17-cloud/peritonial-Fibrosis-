# scripts/audit_checks/regenerate_e2_violins.py
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

print("=== REGENERATING E2 QC VIOLINS WITH UPDATED LEGEND ===")
h5ad_path = "data/processed/GSE248762_hubblind_allcells_qc.h5ad"
print(f"Loading {h5ad_path}...")
adata = ad.read_h5ad(h5ad_path)
obs = adata.obs.copy()

MITO_CEIL = 15.0

# ---------------------------------------------------------------------------
# 1. Per-sample QC violins: E2_qc_violins_per_sample.png
# ---------------------------------------------------------------------------
order = sorted(obs["donor_id"].unique())
# Match sample order: LV_UF-1..4, LV_NOT_UF-1..6, SV-1..6
uf = [f"LV_UF-{i}" for i in range(1, 5)]
not_uf = [f"LV_NOT_UF-{i}" for i in range(1, 7)]
sv = [f"SV-{i}" for i in range(1, 7)]
order = [d for d in uf + not_uf + sv if d in obs["donor_id"].unique()]

before = obs.assign(Status="Deposited barcodes (>= ~500 UMI)")
after = obs[obs["keep"]].assign(Status="After QC + doublet removal")
dfp = pd.concat([before, after])
dfp["donor_id"] = pd.Categorical(dfp["donor_id"], categories=order, ordered=True)
thr = obs.groupby("donor_id", observed=True)[["umi_lo", "umi_hi", "genes_lo", "genes_hi"]].first().reindex(order)

fig, axes = plt.subplots(3, 1, figsize=(24, 17))
pal = {"Deposited barcodes (>= ~500 UMI)": "#94a3b8", "After QC + doublet removal": "#2563eb"}
for ax, col, ylab, logy in [(axes[0], "n_counts", "UMI per cell (log10)", True),
                            (axes[1], "n_genes", "Genes per cell (log10)", True),
                            (axes[2], "pct_counts_mt", "Mitochondrial %", False)]:
    sns.violinplot(data=dfp, x="donor_id", y=col, hue="Status", split=True, inner="quartile",
                   cut=0, density_norm="width", palette=pal, ax=ax, linewidth=0.8)
    if logy:
        ax.set_yscale("log")
    for i, d in enumerate(order):
        if col == "n_counts":
            ax.hlines([thr.loc[d, "umi_lo"], thr.loc[d, "umi_hi"]], i - 0.45, i + 0.45, colors="#dc2626", lw=1.4)
        elif col == "n_genes":
            ax.hlines([thr.loc[d, "genes_lo"], thr.loc[d, "genes_hi"]], i - 0.45, i + 0.45, colors="#dc2626", lw=1.4)
    if col == "pct_counts_mt":
        ax.axhline(MITO_CEIL, color="#dc2626", lw=1.4, ls="--")
    ax.set_ylabel(ylab, fontsize=13)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=12, rotation=30)
    ax.legend(loc="upper right", fontsize=10)
axes[0].set_title("Per-sample QC (not pooled). Red bars = that sample's 3-MAD limits (UMI floor 500, gene floor 200); dashed = 15% mito ceiling",
                  fontsize=14, fontweight="bold")
plt.tight_layout()
out_sample = "results/figures/E2_qc_violins_per_sample.png"
plt.savefig(out_sample, dpi=130)
plt.close()
print(f"Saved: {out_sample}")

# ---------------------------------------------------------------------------
# 2. Before/After QC violins by group: E2_qc_violins_before_after.png
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
group_order = ["LV_UF", "LV_NOT_UF", "SV"]
dfp["group"] = pd.Categorical(dfp["group"], categories=group_order, ordered=True)

# A. UMI
sns.violinplot(data=dfp, x="group", y="n_counts", hue="Status", split=True, inner="quartile",
               palette=pal, ax=axes[0], cut=0)
axes[0].set_yscale("log")
axes[0].set_title("A. UMI Counts per Cell", fontweight="bold", fontsize=14)
axes[0].set_ylabel("Total UMI Counts (log scale)")
axes[0].set_xlabel("Clinical Group")
axes[0].legend(loc="upper right", fontsize=9)

# B. Genes
sns.violinplot(data=dfp, x="group", y="n_genes", hue="Status", split=True, inner="quartile",
               palette=pal, ax=axes[1], cut=0)
axes[1].set_yscale("log")
axes[1].set_title("B. Detected Genes per Cell", fontweight="bold", fontsize=14)
axes[1].set_ylabel("Number of Genes (log scale)")
axes[1].set_xlabel("Clinical Group")
axes[1].legend(loc="upper right", fontsize=9)

# C. Mitochondrial Fraction
sns.violinplot(data=dfp, x="group", y="pct_counts_mt", hue="Status", split=True, inner="quartile",
               palette=pal, ax=axes[2], cut=0)
axes[2].axhline(MITO_CEIL, color="#dc2626", linestyle="--", linewidth=1.5, label=f"{MITO_CEIL:.0f}% Mito Ceiling")
axes[2].set_title("C. Mitochondrial Transcript %", fontweight="bold", fontsize=14)
axes[2].set_ylabel("Mitochondrial Counts (%)")
axes[2].set_xlabel("Clinical Group")
axes[2].legend(loc="upper right", fontsize=9)

plt.tight_layout()
out_before_after = "results/figures/E2_qc_violins_before_after.png"
plt.savefig(out_before_after, dpi=300)
plt.close()
print(f"Saved: {out_before_after}")

print("Both violin plots regenerated successfully.")
