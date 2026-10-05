"""
Script 37a: Prepare pseudobulk matrices across all 37,487 genes (37,476 hub-blind + 11 unsealed hubs)
Creates pseudobulk raw counts per donor for primary analysis and all sensitivity subsets:
- Primary: all kept cells (clusters 4 & 12 excluded)
- Sensitivity (b): stromal without LV_UF-3
- Sensitivity (c): after removing cells flagged by scDblFinder
- Sensitivity (d): stromal using only cells with >= 700 genes
- Sensitivity (e): stromal with donors >= 50 stromal cells
- Leave-one-donor-out (a): 16 subsets per cell type
Seed 42.
"""
import os
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse

np.random.seed(42)

print("=" * 80)
print("STAGE 4: PREPARING PSEUDOBULK MATRICES ACROSS ALL 37,487 GENES")
print("=" * 80)

# 1. Load full transcriptome (37,476 genes) and unsealed hubs (11 genes)
print("Loading data/processed/GSE248762_hubblind_allcells_qc.h5ad...")
adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad")

print("Loading sealed/hub_counts.h5ad...")
adata_hubs = ad.read_h5ad("sealed/hub_counts.h5ad")

# Load cell annotations
adata_annot = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad")
if "_index" in adata_annot.obs.columns:
    kept_barcodes = adata_annot.obs["_index"].values
else:
    kept_barcodes = adata_annot.obs.index.values

print(f"Total kept annotated cells: {len(kept_barcodes)}")

# Align barcodes
adata_qc_sub = adata_qc[kept_barcodes, :].copy()
adata_hubs_sub = adata_hubs[kept_barcodes, :].copy()

# Add cell metadata
for col in ["donor_id", "group", "cell_type", "leiden", "n_counts", "n_genes", "pct_counts_mt"]:
    adata_qc_sub.obs[col] = adata_annot.obs[col].values
    adata_hubs_sub.obs[col] = adata_annot.obs[col].values

# Load scDblFinder calls
df_scd = pd.read_csv("results/tables/scDblFinder_per_barcode_calls.csv")
scd_map = dict(zip(df_scd["barcode"], (df_scd["scDblFinder_class"] == "doublet")))
adata_qc_sub.obs["is_doublet_scDblFinder"] = [scd_map.get(b, False) for b in kept_barcodes]
adata_hubs_sub.obs["is_doublet_scDblFinder"] = adata_qc_sub.obs["is_doublet_scDblFinder"].values

# Combine matrices horizontally: 37,476 + 11 = 37,487 genes
print("Concatenating 37,476 blind genes + 11 unsealed hub genes...")
combined_X = sparse.hstack([adata_qc_sub.X, adata_hubs_sub.X]).tocsr()
combined_var_names = list(adata_qc_sub.var_names) + list(adata_hubs_sub.var_names)

adata_full = ad.AnnData(
    X=combined_X,
    obs=adata_qc_sub.obs.copy(),
    var=pd.DataFrame(index=combined_var_names)
)
print(f"Full transcriptome cohort: {adata_full.shape[0]} cells x {adata_full.shape[1]} genes")

# Exclude donor-dominated clusters 4 and 12
adata_primary = adata_full[~adata_full.obs["leiden"].isin(["4", "12"])].copy()
print(f"Primary clean cohort: {adata_primary.shape[0]} cells x {adata_primary.shape[1]} genes")

# Save pseudobulk output directory
pb_dir = "data/processed/pseudobulk"
os.makedirs(pb_dir, exist_ok=True)

# Helper function to compute pseudobulk raw counts
def aggregate_pseudobulk(adata_sub, cell_type_name):
    ct_mask = (adata_sub.obs["cell_type"] == cell_type_name).values
    adata_ct = adata_sub[ct_mask].copy()
    if adata_ct.n_obs == 0:
        return None, None
    
    donors = sorted(adata_ct.obs["donor_id"].unique())
    donor_map = {d: i for i, d in enumerate(donors)}
    donor_indices = adata_ct.obs["donor_id"].map(donor_map).values
    
    n_donors = len(donors)
    n_cells = adata_ct.n_obs
    
    # Indicator matrix
    S = sparse.csr_matrix((np.ones(n_cells), (donor_indices, np.arange(n_cells))), shape=(n_donors, n_cells))
    pb_counts = S @ adata_ct.X
    
    # Metadata per donor
    donor_meta = adata_ct.obs.groupby("donor_id").agg(
        group=("group", "first"),
        n_cells=("group", "size"),
        mean_umi=("n_counts", "mean"),
        mean_genes=("n_genes", "mean")
    ).loc[donors]
    
    df_pb = pd.DataFrame(pb_counts.toarray(), index=donors, columns=adata_ct.var_names)
    return df_pb, donor_meta

# Target cell types to test
TEST_CTS = [
    "stromal / mesothelial-lineage (unresolved)",
    "Monocyte / macrophage",
    "cDC",
    "T cell",
    "NK cell"
]

print("\n--- Generating Primary Pseudobulk Matrices ---")
for ct in TEST_CTS:
    clean_ct_name = ct.replace(" / ", "_").replace(" ", "_").replace("(", "").replace(")", "")
    df_pb, meta = aggregate_pseudobulk(adata_primary, ct)
    df_pb.to_csv(f"{pb_dir}/pb_primary_{clean_ct_name}_counts.csv")
    meta.to_csv(f"{pb_dir}/pb_primary_{clean_ct_name}_metadata.csv")
    print(f"  {ct}: {len(meta)} donors, {meta['n_cells'].sum()} cells. Saved {clean_ct_name}")

# Sensitivity (b): stromal without LV_UF-3
print("\n--- Generating Sensitivity (b): Stromal without LV_UF-3 ---")
sub_b = adata_primary[adata_primary.obs["donor_id"] != "LV_UF-3"]
df_pb_b, meta_b = aggregate_pseudobulk(sub_b, "stromal / mesothelial-lineage (unresolved)")
df_pb_b.to_csv(f"{pb_dir}/pb_sens_b_no_LV_UF3_counts.csv")
meta_b.to_csv(f"{pb_dir}/pb_sens_b_no_LV_UF3_metadata.csv")
print(f"  Sens (b): {len(meta_b)} donors, {meta_b['n_cells'].sum()} cells.")

# Sensitivity (c): after removing cells flagged by scDblFinder
print("\n--- Generating Sensitivity (c): After removing scDblFinder doublets ---")
sub_c = adata_primary[~adata_primary.obs["is_doublet_scDblFinder"]]
for ct in TEST_CTS:
    clean_ct_name = ct.replace(" / ", "_").replace(" ", "_").replace("(", "").replace(")", "")
    df_pb_c, meta_c = aggregate_pseudobulk(sub_c, ct)
    df_pb_c.to_csv(f"{pb_dir}/pb_sens_c_no_scDblFinder_{clean_ct_name}_counts.csv")
    meta_c.to_csv(f"{pb_dir}/pb_sens_c_no_scDblFinder_{clean_ct_name}_metadata.csv")
print("  Sens (c) completed for all test cell types.")

# Sensitivity (d): stromal using only cells with >= 700 genes
print("\n--- Generating Sensitivity (d): Stromal with genes >= 700 ---")
sub_d = adata_primary[(adata_primary.obs["n_genes"] >= 700)]
df_pb_d, meta_d = aggregate_pseudobulk(sub_d, "stromal / mesothelial-lineage (unresolved)")
df_pb_d.to_csv(f"{pb_dir}/pb_sens_d_ge700genes_counts.csv")
meta_d.to_csv(f"{pb_dir}/pb_sens_d_ge700genes_metadata.csv")
print(f"  Sens (d): {len(meta_d)} donors, {meta_d['n_cells'].sum()} cells.")

# Sensitivity (e): stromal pseudobulk only for donors with >= 50 stromal cells
print("\n--- Generating Sensitivity (e): Stromal donors with >= 50 cells ---")
df_pb_st, meta_st = aggregate_pseudobulk(adata_primary, "stromal / mesothelial-lineage (unresolved)")
keep_donors_e = meta_st[meta_st["n_cells"] >= 50].index
df_pb_e = df_pb_st.loc[keep_donors_e]
meta_e = meta_st.loc[keep_donors_e]
df_pb_e.to_csv(f"{pb_dir}/pb_sens_e_ge50cells_counts.csv")
meta_e.to_csv(f"{pb_dir}/pb_sens_e_ge50cells_metadata.csv")
print(f"  Sens (e): {len(meta_e)} donors (donors dropped: {len(meta_st) - len(meta_e)}), {meta_e['n_cells'].sum()} cells.")

# Sensitivity (a): Leave-one-donor-out (16 runs)
print("\n--- Generating Sensitivity (a): Leave-One-Donor-Out Subsets ---")
all_donors = sorted(adata_primary.obs["donor_id"].unique())
for donor_out in all_donors:
    sub_lodo = adata_primary[adata_primary.obs["donor_id"] != donor_out]
    df_pb_lodo, meta_lodo = aggregate_pseudobulk(sub_lodo, "stromal / mesothelial-lineage (unresolved)")
    df_pb_lodo.to_csv(f"{pb_dir}/pb_sens_lodo_st_{donor_out}_counts.csv")
    meta_lodo.to_csv(f"{pb_dir}/pb_sens_lodo_st_{donor_out}_metadata.csv")
print("  Sens (a) completed for all 16 leave-one-donor-out iterations.")

print("\nAll pseudobulk matrices successfully created and saved to data/processed/pseudobulk/.")
