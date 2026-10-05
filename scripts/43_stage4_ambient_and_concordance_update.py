"""
Script 43: Stage 4 Ambient Profile & Cross-Dataset Concordance Relabeling
Rules:
- Numbered script
- No numeric literals hardcoded for data/results
- Raw console output
- Seed 42
- No causal language
"""
import os
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse, stats

np.random.seed(42)

print("=" * 80)
print("STAGE 4: AMBIENT-RNA PROFILING & CONCORDANCE RELABELING")
print("=" * 80)

# ---------------------------------------------------------------------------
# 1. Load Unsealed AnnData Object
# ---------------------------------------------------------------------------
h5ad_clean = "data/processed/GSE248762_stage4_unsealed_hubs_clean.h5ad"
print(f"Loading {h5ad_clean}...")
adata = ad.read_h5ad(h5ad_clean)

hub_genes = ["COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", "EDIL3", "LOX", "INHBA", "ISM1", "COMP"]

# Exclude donor-dominated clusters 4 and 12
adata_clean = adata[~adata.obs["leiden"].isin(["4", "12"])].copy()

# Target cell populations
st_mask = (adata_clean.obs["cell_type"] == "stromal / mesothelial-lineage (unresolved)")
t_mask = (adata_clean.obs["cell_type"] == "T cell")
mo_mask = (adata_clean.obs["cell_type"] == "Monocyte / macrophage")

st_cells = adata_clean[st_mask].copy()
t_cells = adata_clean[t_mask].copy()
mo_cells = adata_clean[mo_mask].copy()

# ---------------------------------------------------------------------------
# 2. Ambient RNA Table Generation
# ---------------------------------------------------------------------------
print("\n>>> 5. AMBIENT REFERENCE PROFILE & WITHIN-STROMAL SPEARMAN CORRELATION:")

# Function to get mean CP10k and percent expressing
def get_expr_stats(subset_adata, genes):
    counts = subset_adata.X
    sums = np.array(counts.sum(axis=1)).flatten()
    sums[sums == 0] = 1.0
    cp10k = counts.multiply(10000.0 / sums[:, None]).tocsr()
    
    mean_cp10k = np.array(cp10k.mean(axis=0)).flatten()
    pct_expr = np.array((counts > 0).mean(axis=0)).flatten() * 100.0
    
    gene_map = {g: i for i, g in enumerate(subset_adata.var_names)}
    res = {}
    for g in genes:
        if g in gene_map:
            idx = gene_map[g]
            res[g] = (mean_cp10k[idx], pct_expr[idx])
        else:
            res[g] = (np.nan, np.nan)
    return res

st_stats = get_expr_stats(st_cells, hub_genes)
t_stats = get_expr_stats(t_cells, hub_genes)
mo_stats = get_expr_stats(mo_cells, hub_genes)

# Compute within-stromal Spearman correlation between log1p(CP10k) and total UMI (n_counts)
st_counts = st_cells.X
st_sums = np.array(st_counts.sum(axis=1)).flatten()
st_sums_safe = st_sums.copy()
st_sums_safe[st_sums_safe == 0] = 1.0
st_cp10k = st_counts.multiply(10000.0 / st_sums_safe[:, None]).tocsr()
st_log_expr = st_cp10k.copy()
st_log_expr.data = np.log1p(st_log_expr.data)

st_var_map = {g: i for i, g in enumerate(st_cells.var_names)}

ambient_rows = []
for g in hub_genes:
    st_mean, st_pct = st_stats[g]
    t_mean, t_pct = t_stats[g]
    mo_mean, mo_pct = mo_stats[g]
    
    # Spearman correlation in stromal cells
    if g in st_var_map:
        g_idx = st_var_map[g]
        g_vec = np.array(st_log_expr[:, g_idx].toarray()).flatten()
        spear_r, spear_p = stats.spearmanr(g_vec, st_sums)
    else:
        spear_r, spear_p = np.nan, np.nan
        
    ambient_rows.append({
        "Gene": g,
        "Stromal_Mean_CP10k": round(float(st_mean), 4),
        "Stromal_Pct_Expr_%": round(float(st_pct), 2),
        "Tcell_Mean_CP10k": round(float(t_mean), 4),
        "Tcell_Pct_Expr_%": round(float(t_pct), 2),
        "Mono_Macro_Mean_CP10k": round(float(mo_mean), 4),
        "Mono_Macro_Pct_Expr_%": round(float(mo_pct), 2),
        "Stromal_Spearman_rho_vs_TotalUMI": round(float(spear_r), 4),
        "Stromal_Spearman_P": spear_p
    })

df_ambient = pd.DataFrame(ambient_rows)
out_amb_csv = "results/tables/stage4_ambient_reference_comparison.csv"
df_ambient.to_csv(out_amb_csv, index=False)

print("\n--- Ambient RNA Reference Comparison Table (Verbatim from CSV) ---")
print(df_ambient.to_string(index=False))

print("\nMethodological note on Ambient RNA:")
print("  Because only cellranger-filtered matrices (and not raw unfiltered feature-barcode matrices) were deposited")
print("  in GEO GSE248762, empty-droplet profiling algorithms (such as SoupX or CellBender) could not be executed.")
print("  Expression in CD3+ T cells (a lineage where ECM structural collagen and matrisome transcripts are not")
print("  physiologically expected) is provided strictly as an empirical reference point.")
print("  Prior statements asserting that signal 'cannot be accounted for by ambient contamination alone' are rescinded,")
print("  as ambient RNA cannot be rigorously modeled without raw droplet background data.")

# ---------------------------------------------------------------------------
# 3. Cross-Dataset Concordance Relabeling
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print(">>> 6. CROSS-DATASET CONCORDANCE RELABELING & RECOUNT:")
print("=" * 80)

# Load existing concordance table and primary DE table
df_tissue = pd.read_csv("results/tables/stage4_cross_dataset_concordance.csv")
df_prim = pd.read_csv("results/tables/stage4_primary_edger_pseudobulk.csv")

# Filter primary stromal results for LV_UF_vs_LV_NOT_UF
df_prim_st = df_prim[
    (df_prim["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
    (df_prim["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].set_index("Gene")

relabelled_rows = []
for idx, row in df_tissue.iterrows():
    g = row["Gene"]
    tissue_lfc = row["Tissue_EPS_vs_Pool_log2FC"]
    tissue_p = row["Tissue_EPS_vs_Pool_P"]
    tissue_dir = row["Tissue_Dir"]
    gse125498_lfc = row["GSE125498_Effluent_log2FC"]
    gse125498_dir = row["GSE125498_Dir"]
    
    if g in df_prim_st.index:
        p_row = df_prim_st.loc[g]
        st_lfc = p_row["log2FC"]
        st_pval = p_row["PValue"]
        st_fdr = p_row["BH_FDR"]
        ev = p_row["Evidence_Status"]
    else:
        st_lfc, st_pval, st_fdr, ev = np.nan, np.nan, np.nan, "INSUFFICIENT DATA"
        
    # Relabeling logic:
    # Categories:
    # 1. concordant (P<0.05)
    # 2. same direction (not significant)
    # 3. opposite direction (P<0.05, FDR>=0.05)
    # 4. opposite direction (not significant)
    # 5. INSUFFICIENT DATA
    if ev == "INSUFFICIENT DATA" or pd.isna(st_lfc):
        category = "INSUFFICIENT DATA"
    elif tissue_lfc > 0 and st_lfc > 0:
        if st_pval < 0.05:
            category = "concordant (P<0.05)"
        else:
            category = "same direction (not significant)"
    elif tissue_lfc > 0 and st_lfc < 0:
        if st_pval < 0.05:
            category = "opposite direction (P<0.05, FDR>=0.05)"
        else:
            category = "opposite direction (not significant)"
    else:
        category = "unclassified"
        
    relabelled_rows.append({
        "Gene": g,
        "Tissue_EPS_vs_Pool_log2FC": tissue_lfc,
        "Tissue_EPS_vs_Pool_P": tissue_p,
        "Tissue_Dir": tissue_dir,
        "Stromal_Effluent_LV_UF_vs_LV_NOT_UF_log2FC": st_lfc,
        "Stromal_Effluent_PValue": st_pval,
        "Stromal_Effluent_BH_FDR": st_fdr,
        "GSE125498_Effluent_log2FC": gse125498_lfc,
        "GSE125498_Dir": gse125498_dir,
        "Evidence_Status": ev,
        "Concordance_Tissue_vs_Stromal_Effluent": category
    })

df_concordance = pd.DataFrame(relabelled_rows)
out_conc_csv = "results/tables/stage4_cross_dataset_concordance_relabelled.csv"
df_concordance.to_csv(out_conc_csv, index=False)

# Update the main concordance file as well
df_concordance.to_csv("results/tables/stage4_cross_dataset_concordance.csv", index=False)

print("\nTABLE HEADER NOTE:")
print("  Compartment Distinction: Tissue = EPS vs PD/uremic (4 vs 4 pooled peritoneal tissue biopsies, GSE62928);")
print("  Effluent scRNA = ultrafiltration failure vs no ultrafiltration failure (peritoneal dialysis effluent dialysate cells, GSE248762).")
print("  These biological compartments and clinical comparisons are not equivalent.\n")

print("--- Relabelled Cross-Dataset Concordance Table (Verbatim from CSV) ---")
print(df_concordance.to_string(index=False))

print("\n--- Summary Recount by Category ---")
cat_counts = df_concordance["Concordance_Tissue_vs_Stromal_Effluent"].value_counts()
print(cat_counts.to_string())

print("\nAmbient profiling and concordance relabeling completed successfully.")
print("=" * 80)
