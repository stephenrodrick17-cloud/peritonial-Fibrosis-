"""
Script 36: Stage 4 Step 1 (Unseal) & Step 2 (Descriptive Profiling)
Rules:
- Verify SHA256 matches sealing hash
- Merge hub counts by barcode into annotated object (95,247 cells)
- Log guard release to audit/guard_audit.log
- Descriptive stats: % cells expressing, mean log1p(CP10k), number of donors with >=1 count (overall and per group)
- Flag pairs with <2 donors per group as INSUFFICIENT DATA
- Generate unscaled and scaled dotplots and violins
Seed 42, no causal language.
"""
import os
import hashlib
import json
import datetime
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(42)

print("=" * 80)
print("STAGE 4 - STEP 1: UNSEAL & STEP 2: DESCRIPTIVE PROFILING")
print("=" * 80)

# ---------------------------------------------------------------------------
# Step 1: Unseal
# ---------------------------------------------------------------------------
print("\n--- STEP 1: UNSEAL VERIFICATION ---")
sealed_file = "sealed/hub_counts.h5ad"
with open(sealed_file, "rb") as f:
    current_sha = hashlib.sha256(f.read()).hexdigest()

with open("provenance/seal_record.json") as f:
    seal_rec = json.load(f)

expected_sha = seal_rec["original_sealing_sha256"]
print(f"Current SHA256:  {current_sha}")
print(f"Expected SHA256: {expected_sha}")
assert current_sha == expected_sha, "CRITICAL ERROR: Sealed file checksum mismatch!"
print("Checksum verified: SEAL INTEGRITY CONFIRMED.")

# Log guard release
timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
guard_release_msg = f"[GUARD RELEASE - AUTHORIZED] {timestamp} | Stage 4 scRNA Validation | 11 Hub genes unsealed for Stage 4 descriptive & pseudobulk analysis only."
with open("audit/guard_audit.log", "a") as f:
    f.write(guard_release_msg + "\n")
print("Logged guard release to audit/guard_audit.log.")

# Load hub counts
print(f"Loading {sealed_file}...")
adata_hubs = ad.read_h5ad(sealed_file)
print(f"Loaded hub AnnData: {adata_hubs.shape[0]} cells x {adata_hubs.shape[1]} genes")
hub_genes = list(adata_hubs.var_names)
print(f"Hub genes ({len(hub_genes)}): {hub_genes}")

# Load annotated obs object
print("Loading data/processed/GSE248762_harmony_annotated_obs.h5ad...")
adata_annot = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad")
print(f"Annotated cells: {adata_annot.n_obs}")

# Barcode alignment
if "_index" in adata_annot.obs.columns:
    cell_barcodes = adata_annot.obs["_index"].values
else:
    cell_barcodes = adata_annot.obs.index.values

# Verify alignment
assert len(cell_barcodes) == adata_annot.n_obs, "Barcode length mismatch!"
hubs_sub = adata_hubs[cell_barcodes, :].copy()
assert hubs_sub.n_obs == adata_annot.n_obs, f"Cell count mismatch: {hubs_sub.n_obs} != {adata_annot.n_obs}"
print(f"Successfully aligned hub counts for all {hubs_sub.n_obs:,} kept annotated cells.")

# Create merged AnnData for hub genes with full obs metadata
hubs_merged = ad.AnnData(
    X=hubs_sub.X.copy(),
    obs=adata_annot.obs.copy(),
    var=pd.DataFrame(index=hub_genes)
)
if "_index" in hubs_merged.obs.columns:
    hubs_merged.obs.index = hubs_merged.obs["_index"]

# Log-normalize hub counts (CP10k + log1p)
# Total counts per cell from full transcriptomes
total_counts = hubs_merged.obs["n_counts"].values.astype(float)
total_counts[total_counts == 0] = 1.0

# Store raw counts in .layers["raw_counts"]
if sparse.issparse(hubs_merged.X):
    raw_mat = hubs_merged.X.tocsr()
else:
    raw_mat = sparse.csr_matrix(hubs_merged.X)
hubs_merged.layers["raw_counts"] = raw_mat.copy()

# Normalized expression
norm_mat = raw_mat.multiply(10000.0 / total_counts[:, None]).tocsr()
norm_mat.data = np.log1p(norm_mat.data)
hubs_merged.X = norm_mat

# Exclude donor-dominated clusters (4: Neutrophils, 12: Mono LV_UF-3)
hubs_clean = hubs_merged[~hubs_merged.obs["leiden"].isin(["4", "12"])].copy()
print(f"Filtered cohort (excluding clusters 4 & 12): {hubs_clean.n_obs:,} cells across {hubs_clean.obs['donor_id'].nunique()} donors.")

# ---------------------------------------------------------------------------
# Step 2: Descriptive Profiling (No Inference)
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("STEP 2: DESCRIPTIVE PROFILING OF 11 HUB GENES")
print("=" * 80)

CELL_TYPES = [
    "stromal / mesothelial-lineage (unresolved)",
    "Monocyte / macrophage",
    "cDC",
    "T cell",
    "NK cell",
    "B cell",
    "Mixed: pDC / Plasma cell"
]

desc_rows = []
raw_dense = hubs_clean.layers["raw_counts"].toarray()
norm_dense = hubs_clean.X.toarray()

for ig, g in enumerate(hub_genes):
    g_raw = raw_dense[:, ig]
    g_norm = norm_dense[:, ig]
    
    for ct in CELL_TYPES:
        ct_mask = (hubs_clean.obs["cell_type"] == ct).values
        n_ct_cells = int(ct_mask.sum())
        if n_ct_cells == 0:
            continue
            
        sub_raw = g_raw[ct_mask]
        sub_norm = g_norm[ct_mask]
        sub_obs = hubs_clean.obs.iloc[ct_mask]
        
        pct_expr = round(float((sub_raw > 0).mean() * 100), 2)
        mean_norm = round(float(sub_norm.mean()), 3)
        mean_expr_pos = round(float(sub_norm[sub_raw > 0].mean()), 3) if (sub_raw > 0).sum() > 0 else 0.0
        
        # Donors with >= 1 count
        donors_with_count = sub_obs.loc[sub_raw > 0, "donor_id"].nunique()
        tot_donors = sub_obs["donor_id"].nunique()
        
        # Donors with >= 1 count by group
        counts_by_grp = {}
        for grp in ["LV_UF", "LV_NOT_UF", "SV"]:
            grp_mask = (sub_obs["group"] == grp)
            tot_grp_donors = sub_obs.loc[grp_mask, "donor_id"].nunique()
            det_grp_donors = sub_obs.loc[grp_mask & (sub_raw > 0), "donor_id"].nunique()
            counts_by_grp[grp] = f"{det_grp_donors}/{tot_grp_donors}"
            
        # Evidence rule check: < 2 donors in ANY group with counts -> INSUFFICIENT DATA
        min_donors_in_any_group = min(
            sub_obs.loc[(sub_obs["group"] == grp) & (sub_raw > 0), "donor_id"].nunique()
            for grp in ["LV_UF", "LV_NOT_UF", "SV"]
        )
        status = "PASS" if min_donors_in_any_group >= 2 else "INSUFFICIENT DATA"
        
        desc_rows.append({
            "Gene": g,
            "Cell_Type": ct,
            "N_cells": n_ct_cells,
            "Pct_Expressing_%": pct_expr,
            "Mean_log1p_CP10k": mean_norm,
            "Mean_if_pos": mean_expr_pos,
            "Donors_ge1_count": f"{donors_with_count}/{tot_donors}",
            "LV_UF_donors": counts_by_grp["LV_UF"],
            "LV_NOT_UF_donors": counts_by_grp["LV_NOT_UF"],
            "SV_donors": counts_by_grp["SV"],
            "Evidence_Status": status
        })

df_desc = pd.DataFrame(desc_rows)
os.makedirs("results/tables", exist_ok=True)
df_desc.to_csv("results/tables/stage4_hub_genes_descriptive_summary.csv", index=False)

print("\n--- Descriptive Summary for 11 Hub Genes in Stromal & Immune Lineages ---")
# Show stromal first
print(df_desc[df_desc["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)"].to_string(index=False))

print("\n--- Evidence Status Summary across All Tested Lineages ---")
evidence_counts = df_desc.groupby(["Cell_Type", "Evidence_Status"]).size().unstack(fill_value=0)
print(evidence_counts.to_string())

# ---------------------------------------------------------------------------
# Generate Dotplots (Unscaled and Scaled) and Violins
# ---------------------------------------------------------------------------
print("\n--- Generating Dot Plots and Violins for Hub Genes ---")
os.makedirs("results/figures", exist_ok=True)

# 1. Unscaled Dotplot
test_cts = ["stromal / mesothelial-lineage (unresolved)", "Monocyte / macrophage", "cDC", "T cell", "NK cell"]
hubs_plot = hubs_clean[hubs_clean.obs["cell_type"].isin(test_cts)].copy()
hubs_plot.obs["cell_type"] = hubs_plot.obs["cell_type"].astype(str)

sc_dp_unscaled = sc.pl.dotplot(
    hubs_plot,
    var_names=hub_genes,
    groupby="cell_type",
    standard_scale=None,
    return_fig=True,
    colorbar_title="Mean log1p(CP10k)"
) if "sc" in locals() else None

# Build clean matplotlib dotplot directly
df_dp = []
for g in hub_genes:
    ig = hub_genes.index(g)
    for ct in test_cts:
        sub = hubs_plot[hubs_plot.obs["cell_type"] == ct]
        raw_vals = sub.layers["raw_counts"][:, ig].toarray().flatten()
        norm_vals = sub.X[:, ig].toarray().flatten()
        pct = float((raw_vals > 0).mean() * 100)
        mean_exp = float(norm_vals.mean())
        df_dp.append({"Gene": g, "Cell_Type": ct, "Pct": pct, "Mean": mean_exp})

df_dp = pd.DataFrame(df_dp)

# Unscaled Dot Plot
fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
ct_order = test_cts
gene_order = hub_genes
norm_mean = df_dp["Mean"]

scatter = ax.scatter(
    x=[gene_order.index(g) for g in df_dp["Gene"]],
    y=[ct_order.index(ct) for ct in df_dp["Cell_Type"]],
    s=df_dp["Pct"] * 3.0 + 5,
    c=df_dp["Mean"],
    cmap="Reds",
    edgecolors="black",
    linewidths=0.5
)
ax.set_xticks(range(len(gene_order)))
ax.set_xticklabels(gene_order, rotation=45, ha="right", fontweight="bold")
ax.set_yticks(range(len(ct_order)))
ax.set_yticklabels(ct_order, fontweight="bold")
ax.set_title("11 Hub Genes Expression Across Lineages (Unscaled Mean log1p(CP10k))", fontsize=12, fontweight="bold", pad=12)
cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
cbar.set_label("Mean log1p(CP10k)")
plt.tight_layout()
fig.savefig("results/figures/stage4_hub_genes_dotplot_unscaled.png")
plt.close(fig)

# Scaled Dot Plot (Min-Max Scaled 0-1 per gene)
df_dp["Mean_scaled"] = df_dp.groupby("Gene")["Mean"].transform(lambda x: (x - x.min()) / (x.max() - x.min() + 1e-9))
fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
scatter = ax.scatter(
    x=[gene_order.index(g) for g in df_dp["Gene"]],
    y=[ct_order.index(ct) for ct in df_dp["Cell_Type"]],
    s=df_dp["Pct"] * 3.0 + 5,
    c=df_dp["Mean_scaled"],
    cmap="viridis",
    edgecolors="black",
    linewidths=0.5
)
ax.set_xticks(range(len(gene_order)))
ax.set_xticklabels(gene_order, rotation=45, ha="right", fontweight="bold")
ax.set_yticks(range(len(ct_order)))
ax.set_yticklabels(ct_order, fontweight="bold")
ax.set_title("11 Hub Genes Expression Across Lineages (Per-Gene Min-Max Scaled 0-1)", fontsize=12, fontweight="bold", pad=12)
cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
cbar.set_label("Scaled Expression (0-1)")
plt.tight_layout()
fig.savefig("results/figures/stage4_hub_genes_dotplot_scaled.png")
plt.close(fig)

# Violins for Stromal Cells Across Groups (for genes passing evidence)
st_cells = hubs_plot[hubs_plot.obs["cell_type"] == "stromal / mesothelial-lineage (unresolved)"].copy()
st_df_expr = pd.DataFrame(st_cells.X.toarray(), columns=hub_genes, index=st_cells.obs_names)
st_df_expr["Group"] = st_cells.obs["group"].values
st_df_expr["Donor"] = st_cells.obs["donor_id"].values

fig, axes = plt.subplots(3, 4, figsize=(16, 11), dpi=150)
axes = axes.flatten()
for ig, g in enumerate(hub_genes):
    ax = axes[ig]
    sns.violinplot(data=st_df_expr, x="Group", y=g, order=["SV", "LV_NOT_UF", "LV_UF"],
                   palette=["#2b83ba", "#abdda4", "#d7191c"], ax=ax, inner="quartile", cut=0)
    ax.set_title(g, fontweight="bold", fontsize=11)
    ax.set_ylabel("log1p(CP10k)")
    ax.set_xlabel("")
# Hide 12th subplot
fig.delaxes(axes[11])
plt.suptitle("Stromal / Mesothelial-Lineage: Hub Gene Expression by PD Group", fontsize=14, fontweight="bold", y=0.99)
plt.tight_layout()
fig.savefig("results/figures/stage4_hub_genes_stromal_violins.png")
plt.close(fig)

# Save unsealed merged object for downstream pseudobulk
if "_index" in hubs_clean.obs.columns:
    hubs_clean.obs = hubs_clean.obs.drop(columns=["_index"])
hubs_clean.write_h5ad("data/processed/GSE248762_stage4_unsealed_hubs_clean.h5ad")
print("Saved data/processed/GSE248762_stage4_unsealed_hubs_clean.h5ad.")
print("Saved figures:")
print("  - results/figures/stage4_hub_genes_dotplot_unscaled.png")
print("  - results/figures/stage4_hub_genes_dotplot_scaled.png")
print("  - results/figures/stage4_hub_genes_stromal_violins.png")
print("Step 1 & Step 2 completed successfully.")
print("=" * 80)
