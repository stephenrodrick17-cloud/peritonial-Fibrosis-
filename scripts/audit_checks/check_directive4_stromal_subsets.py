# scripts/audit_checks/check_directive4_stromal_subsets.py
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats

print("=== CHECK DIRECTIVE 4: Stromal Heterogeneity Analysis ===")

# Load harmony obs to get stromal barcodes and metadata
obs_annot = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad").obs
st_mask = obs_annot["cell_type"] == "stromal / mesothelial-lineage (unresolved)"
obs_st = obs_annot[st_mask].copy()
st_barcodes = obs_st["_index"].values
print(f"Total stromal cells: {len(st_barcodes)}")

# Load full QC matrix to get raw/normalized gene expression
adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad")
adata_st = adata_qc[st_barcodes, :].copy()
adata_st.obs["donor_id"] = obs_st["donor_id"].values
adata_st.obs["group"] = obs_st["group"].values

# Check candidate genes
meso_genes = ["WT1", "MSLN", "UPK3B", "KRT19"]
fibro_genes = ["DCN", "PDGFRA", "LUM", "COL1A2"]

print("Mesothelial marker genes in var:", [g for g in meso_genes if g in adata_st.var_names])
print("Fibroblast marker genes in var:", [g for g in fibro_genes if g in adata_st.var_names])

# Normalize per cell (log1p CPM / 10,000)
sc.pp.normalize_total(adata_st, target_sum=1e4)
sc.pp.log1p(adata_st)

# Calculate signature scores using scanpy score_genes
sc.tl.score_genes(adata_st, meso_genes, score_name="meso_score")
sc.tl.score_genes(adata_st, fibro_genes, score_name="fibro_score")

# Assignment: cell is mesothelial-skewed if meso_score > fibro_score, else fibroblast-skewed
adata_st.obs["skew"] = np.where(adata_st.obs["meso_score"] > adata_st.obs["fibro_score"], "Mesothelial", "Fibroblast")

# Per-donor metrics
donor_stats = []
donors = sorted(adata_st.obs["donor_id"].unique())
for d in donors:
    sub = adata_st.obs[adata_st.obs["donor_id"] == d]
    n_tot = len(sub)
    n_meso = (sub["skew"] == "Mesothelial").sum()
    n_fibro = (sub["skew"] == "Fibroblast").sum()
    pct_meso = n_meso / n_tot * 100
    pct_fibro = n_fibro / n_tot * 100
    grp = sub["group"].iloc[0]
    mean_umi = sub["n_counts"].mean()
    mean_umi_meso = sub.loc[sub["skew"] == "Mesothelial", "n_counts"].mean() if n_meso > 0 else np.nan
    mean_umi_fibro = sub.loc[sub["skew"] == "Fibroblast", "n_counts"].mean() if n_fibro > 0 else np.nan
    donor_stats.append({
        "Donor": d,
        "Group": grp,
        "Total_Stromal": n_tot,
        "N_Mesothelial": n_meso,
        "Pct_Mesothelial": pct_meso,
        "N_Fibroblast": n_fibro,
        "Pct_Fibroblast": pct_fibro,
        "Mean_UMI_All": mean_umi,
        "Mean_UMI_Mesothelial": mean_umi_meso,
        "Mean_UMI_Fibroblast": mean_umi_fibro
    })

df_donor = pd.DataFrame(donor_stats)
print("\n=== PER-DONOR STROMAL SUBSET BREAKDOWN ===")
print(df_donor.to_string(index=False))

# Group summaries
print("\n=== GROUP COMPARISON (Mean ± SD) ===")
for g in ["LV_UF", "LV_NOT_UF", "SV"]:
    sub = df_donor[df_donor["Group"] == g]
    print(f"Group {g} (N={len(sub)} donors):")
    print(f"  Pct Mesothelial: {sub['Pct_Mesothelial'].mean():.2f}% ± {sub['Pct_Mesothelial'].std():.2f}%")
    print(f"  Mean UMI (All):  {sub['Mean_UMI_All'].mean():.1f} ± {sub['Mean_UMI_All'].std():.1f}")

# Statistical testing across the 3 groups
kw_pct = stats.kruskal(*[df_donor.loc[df_donor["Group"] == g, "Pct_Mesothelial"] for g in ["LV_UF", "LV_NOT_UF", "SV"]])
kw_umi = stats.kruskal(*[df_donor.loc[df_donor["Group"] == g, "Mean_UMI_All"] for g in ["LV_UF", "LV_NOT_UF", "SV"]])

print("\n=== STATISTICAL COMPARISON ACROSS GROUPS ===")
print(f"Kruskal-Wallis test on Pct Mesothelial across 3 groups: H = {kw_pct.statistic:.4f}, P = {kw_pct.pvalue:.4f}")
print(f"Kruskal-Wallis test on Mean UMI across 3 groups:        H = {kw_umi.statistic:.4f}, P = {kw_umi.pvalue:.4f}")

# Contrast: LV_UF vs LV_NOT_UF
u_pct = stats.mannwhitneyu(df_donor.loc[df_donor["Group"] == "LV_UF", "Pct_Mesothelial"],
                           df_donor.loc[df_donor["Group"] == "LV_NOT_UF", "Pct_Mesothelial"])
u_umi = stats.mannwhitneyu(df_donor.loc[df_donor["Group"] == "LV_UF", "Mean_UMI_All"],
                           df_donor.loc[df_donor["Group"] == "LV_NOT_UF", "Mean_UMI_All"])

print(f"\nMann-Whitney U test (LV_UF vs LV_NOT_UF):")
print(f"  Pct Mesothelial: U = {u_pct.statistic}, P = {u_pct.pvalue:.4f}")
print(f"  Mean UMI:        U = {u_umi.statistic}, P = {u_umi.pvalue:.4f}")
