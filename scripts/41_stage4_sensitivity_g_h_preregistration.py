"""
Script 41: Pre-Registration and Preparation of Sensitivity Analyses (g) and (h)
Rules:
- Pre-register in provenance/stage4_addendum.json BEFORE any execution
- Print filesystem mtime of the addendum
- Rescued stromal ceiling cells via hub-blind module scores
- Remove Scrublet (>0.4390) OR scDblFinder doublets
- Prepare pseudobulk matrices across all 37,487 genes
- Seed 42, no causal language.
"""
import os
import json
import datetime
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

np.random.seed(42)

print("=" * 80)
print("STAGE 4: PRE-SPECIFIED (SELF-DOCUMENTED) & PREPARATION FOR SENSITIVITY (g) & (h)")
print("=" * 80)

# ---------------------------------------------------------------------------
# 1. Pre-Specify Sensitivity (g) and (h) in provenance/stage4_addendum.json
# ---------------------------------------------------------------------------
print("\n>>> 1. PRE-SPECIFYING (SELF-DOCUMENTED) SENSITIVITY (g) & (h) IN provenance/stage4_addendum.json...")
addendum_path = "provenance/stage4_addendum.json"

CUTOFF_SCRUBLET = 0.4390
MIN_CELLS_H = 50
MIN_MEDIAN_GENES_H = 1000
MITO_CEILING = 15.0

addendum = {
    "timestamp_addendum_frozen_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "dataset": "GSE248762 (PD effluent scRNA-seq, 16 donors)",
    "status": "pre-specified (self-documented) BEFORE running sensitivity analyses (g) and (h)",
    "sensitivity_g": {
        "title": "Ceiling-rescued stromal/mesothelial lineage pseudobulk",
        "description": (
            "Re-run stromal pseudobulk for all 3 contrasts including cells excluded by upper "
            "UMI/gene ceilings that are assigned to the stromal / mesothelial lineage by "
            "hub-blind module scoring (WT1, MSLN, CALB2, DCN, LUM, PDGFRA; strictly 0 hub genes). "
            "Doublet filtering: remove any cell flagged by Scrublet (> 0.4390) OR scDblFinder ('doublet'). "
            "Mitochondrial ceiling: <= 15%."
        ),
        "scrublet_cutoff": CUTOFF_SCRUBLET,
        "mito_ceiling_pct": MITO_CEILING,
        "contrasts": ["LV_UF_vs_LV_NOT_UF", "LV_UF_vs_SV", "LV_NOT_UF_vs_SV"],
        "method": "edgeR quasi-likelihood across all 37,487 genes (filterByExpr + TMM over all genes)"
    },
    "sensitivity_h": {
        "title": "High-confidence donors only (>=50 stromal cells and median genes >=1000)",
        "description": (
            "Re-run stromal pseudobulk restricting to donors with >= 50 stromal cells AND "
            "median genes per cell >= 1000 in the primary stromal population. "
            "Evaluates whether low-yield or low-complexity donor libraries drive the primary result."
        ),
        "min_cells": MIN_CELLS_H,
        "min_median_genes": MIN_MEDIAN_GENES_H,
        "contrasts": ["LV_UF_vs_LV_NOT_UF", "LV_UF_vs_SV", "LV_NOT_UF_vs_SV"],
        "method": "edgeR quasi-likelihood across all 37,487 genes (filterByExpr + TMM over all genes)"
    },
    "protocol_constraints": {
        "primary_status": "Primary analysis results remain strictly unchanged and are labeled PRIMARY",
        "sensitivity_status": "Analyses (g) and (h) are exploratory sensitivity analyses and are labeled SENSITIVITY",
        "causal_language": "No causal language allowed; strictly associative phrasing"
    }
}

os.makedirs("provenance", exist_ok=True)
with open(addendum_path, "w", encoding="utf-8") as f:
    json.dump(addendum, f, indent=2)

st_add = os.stat(addendum_path)
mtime_add_utc = datetime.datetime.fromtimestamp(st_add.st_mtime, tz=datetime.timezone.utc).isoformat()
print(f"Pre-specification (self-documented) addendum successfully written to: {addendum_path}")
print(f"Filesystem mtime (os.stat UTC): {mtime_add_utc}")

# ---------------------------------------------------------------------------
# 2. Build Population for Sensitivity (g)
# ---------------------------------------------------------------------------
print("\n>>> 2. BUILDING SENSITIVITY (g) STROMAL POPULATION...")

# Load all-cells QC object (37,476 genes) and hub object (11 genes)
adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad")
adata_hubs = ad.read_h5ad("sealed/hub_counts.h5ad")

# Load scDblFinder calls
df_sc = pd.read_csv("results/tables/scDblFinder_per_barcode_calls.csv", index_col="barcode")
adata_qc.obs["scDblFinder_class"] = df_sc.loc[adata_qc.obs.index, "scDblFinder_class"]

# Load primary annotated obs
adata_annot = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad")
primary_obs = adata_annot.obs.copy()
primary_obs.index = primary_obs["_index"].values if "_index" in primary_obs.columns else primary_obs.index.values

# Primary stromal cells (excluding clusters 4 and 12)
st_primary_barcodes = primary_obs[
    (primary_obs["cell_type"] == "stromal / mesothelial-lineage (unresolved)") &
    (~primary_obs["leiden"].isin(["4", "12"]))
].index.values

# Apply doublet filter: remove if Scrublet > 0.4390 OR scDblFinder == 'doublet'
# (Primary cells already have Scrublet <= 0.4390; we now filter on scDblFinder)
st_primary_sub = adata_qc.obs.loc[st_primary_barcodes]
st_prim_keep_mask = (st_primary_sub["doublet_score"] <= CUTOFF_SCRUBLET) & (st_primary_sub["scDblFinder_class"] != "doublet")
st_prim_pass_barcodes = st_primary_sub[st_prim_keep_mask].index.values
print(f"Primary clean stromal cells: {len(st_primary_barcodes)}")
print(f"Primary stromal cells passing joint doublet filter: {len(st_prim_pass_barcodes)}")

# Identify ceiling-excluded cells
obs_all = adata_qc.obs
ceiling_mask = (obs_all["fail_umi_high"] | obs_all["fail_genes_high"]).values
ceiling_adata = adata_qc[ceiling_mask].copy()

# Score ceiling cells with hub-blind lineage modules
MODULES = {
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "Monocyte_macrophage": ["CD14", "FCGR3A", "CD68", "CD163"],
    "cDC": ["CD1C", "CLEC10A", "FCER1A", "CLEC9A"],
    "NK_cell": ["NCAM1", "NKG7", "GNLY", "KLRD1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B"],
    "Plasma_cell": ["MZB1", "SDC1", "JCHAIN"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "MNDA", "G0S2"],
    "Stromal_mesothelial": ["WT1", "MSLN", "CALB2", "DCN", "LUM", "PDGFRA"]
}

sc_counts = ceiling_adata.X
sc_sums = np.array(sc_counts.sum(axis=1)).flatten()
sc_sums[sc_sums == 0] = 1.0
norm_X = sc_counts.multiply(10000.0 / sc_sums[:, None]).tocsr()
norm_X.data = np.log1p(norm_X.data)

var_dict = {g: i for i, g in enumerate(ceiling_adata.var_names)}
scores = {}
for mod_name, genes in MODULES.items():
    present_genes = [g for g in genes if g in var_dict]
    if present_genes:
        col_indices = [var_dict[g] for g in present_genes]
        scores[mod_name] = np.array(norm_X[:, col_indices].mean(axis=1)).flatten()
    else:
        scores[mod_name] = np.zeros(ceiling_adata.n_obs)

df_scores = pd.DataFrame(scores, index=ceiling_adata.obs_names)
top_module = df_scores.idxmax(axis=1)
top_module[df_scores.max(axis=1) == 0] = "Unassigned"

# Filter ceiling cells assigned to stromal/mesothelial lineage
ceiling_st_mask = (top_module == "Stromal_mesothelial").values
ceiling_st_obs = ceiling_adata.obs[ceiling_st_mask].copy()

# Doublet filter on ceiling cells: Scrublet <= 0.4390 AND scDblFinder != 'doublet' AND mito <= 15%
ceiling_st_pass = ceiling_st_obs[
    (ceiling_st_obs["doublet_score"] <= CUTOFF_SCRUBLET) &
    (ceiling_st_obs["scDblFinder_class"] != "doublet") &
    (ceiling_st_obs["pct_counts_mt"] <= MITO_CEILING)
]
ceiling_rescued_barcodes = ceiling_st_pass.index.values

print(f"Ceiling-excluded cells scored as Stromal_mesothelial: {len(ceiling_st_obs)}")
print(f"Ceiling-excluded stromal cells passing doublet & mito filters: {len(ceiling_rescued_barcodes)}")

# Combine primary stromal singlets + ceiling rescued stromal singlets
all_sens_g_barcodes = np.concatenate([st_prim_pass_barcodes, ceiling_rescued_barcodes])
print(f"Total cells in Sensitivity (g) stromal pool: {len(all_sens_g_barcodes)}")

# Extract counts across all 37,487 genes (37,476 non-hub + 11 hubs)
adata_qc_g = adata_qc[all_sens_g_barcodes, :].copy()
adata_hubs_g = adata_hubs[all_sens_g_barcodes, :].copy()
full_mat_g = sparse.hstack([adata_qc_g.X, adata_hubs_g.X]).tocsr()
all_var_names = list(adata_qc_g.var_names) + list(adata_hubs_g.var_names)

adata_g = ad.AnnData(X=full_mat_g, obs=adata_qc_g.obs.copy())
adata_g.var_names = all_var_names

# Aggregate to pseudobulk per donor
donors_g = sorted(adata_g.obs["donor_id"].unique())
donor_map_g = {d: i for i, d in enumerate(donors_g)}
donor_idx_g = adata_g.obs["donor_id"].map(donor_map_g).values
n_donors_g = len(donors_g)
n_cells_g = adata_g.n_obs

S_g = sparse.csr_matrix((np.ones(n_cells_g), (donor_idx_g, np.arange(n_cells_g))), shape=(n_donors_g, n_cells_g))
pb_counts_g = S_g @ adata_g.X

donor_meta_g = adata_g.obs.groupby("donor_id").agg(
    group=("group", "first"),
    n_cells=("group", "size"),
    mean_umi=("n_counts", "mean"),
    mean_genes=("n_genes", "mean")
).loc[donors_g]

df_pb_g = pd.DataFrame(pb_counts_g.toarray(), index=donors_g, columns=adata_g.var_names)

pb_dir = "data/processed/pseudobulk"
df_pb_g.to_csv(f"{pb_dir}/pb_sens_g_ceiling_rescued_stromal_counts.csv")
donor_meta_g.to_csv(f"{pb_dir}/pb_sens_g_ceiling_rescued_stromal_metadata.csv")
print(f"Saved Sensitivity (g) pseudobulk: {len(donor_meta_g)} donors, {donor_meta_g['n_cells'].sum()} total cells.")

# ---------------------------------------------------------------------------
# 3. Build Population for Sensitivity (h)
# ---------------------------------------------------------------------------
print("\n>>> 3. BUILDING SENSITIVITY (h) STROMAL POPULATION...")
# Primary stromal donors with >= 50 cells AND median genes >= 1000
# Load primary stromal metadata and primary stromal cell obs
adata_primary_st = adata_qc[st_primary_barcodes, :].copy()
adata_primary_hubs_st = adata_hubs[st_primary_barcodes, :].copy()
full_mat_st = sparse.hstack([adata_primary_st.X, adata_primary_hubs_st.X]).tocsr()
adata_prim_st = ad.AnnData(X=full_mat_st, obs=primary_obs.loc[st_primary_barcodes].copy())
adata_prim_st.var_names = all_var_names

donor_stats_h = adata_prim_st.obs.groupby("donor_id").agg(
    group=("group", "first"),
    n_cells=("n_counts", "count"),
    median_genes=("n_genes", "median")
)

qualifying_donors_h = donor_stats_h[
    (donor_stats_h["n_cells"] >= MIN_CELLS_H) &
    (donor_stats_h["median_genes"] >= MIN_MEDIAN_GENES_H)
].index.tolist()

print("\nDonor statistics for Sensitivity (h) qualification:")
print(donor_stats_h.to_string())
print(f"\nQualifying donors ({len(qualifying_donors_h)}): {qualifying_donors_h}")
qual_meta = donor_stats_h.loc[qualifying_donors_h]
print("Breakdown by group:")
print(qual_meta["group"].value_counts().to_string())

# Filter cells to qualifying donors
adata_h = adata_prim_st[adata_prim_st.obs["donor_id"].isin(qualifying_donors_h)].copy()

donors_h = sorted(adata_h.obs["donor_id"].unique())
donor_map_h = {d: i for i, d in enumerate(donors_h)}
donor_idx_h = adata_h.obs["donor_id"].map(donor_map_h).values
n_donors_h = len(donors_h)
n_cells_h = adata_h.n_obs

S_h = sparse.csr_matrix((np.ones(n_cells_h), (donor_idx_h, np.arange(n_cells_h))), shape=(n_donors_h, n_cells_h))
pb_counts_h = S_h @ adata_h.X

donor_meta_h = adata_h.obs.groupby("donor_id").agg(
    group=("group", "first"),
    n_cells=("group", "size"),
    mean_umi=("n_counts", "mean"),
    mean_genes=("n_genes", "mean")
).loc[donors_h]

df_pb_h = pd.DataFrame(pb_counts_h.toarray(), index=donors_h, columns=adata_h.var_names)

df_pb_h.to_csv(f"{pb_dir}/pb_sens_h_ge50cells_med1000genes_counts.csv")
donor_meta_h.to_csv(f"{pb_dir}/pb_sens_h_ge50cells_med1000genes_metadata.csv")
print(f"Saved Sensitivity (h) pseudobulk: {len(donor_meta_h)} donors, {donor_meta_h['n_cells'].sum()} total cells.")

# ---------------------------------------------------------------------------
# 4. Print Per-Donor Library Sizes (Cells and Total Counts)
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print(">>> 4. PER-DONOR PSEUDOBULK LIBRARY SIZES FOR (g) and (h):")
print("=" * 80)

def print_lib_sizes(label, counts_df, meta_df):
    tot_counts = counts_df.sum(axis=1)
    df_sizes = pd.DataFrame({
        "Donor": meta_df.index,
        "Group": meta_df["group"],
        "Cells": meta_df["n_cells"],
        "Total_UMI_Counts": tot_counts.loc[meta_df.index]
    })
    print(f"\n[{label}]:")
    print(f"Total donors: {len(df_sizes)}, Total cells: {df_sizes['Cells'].sum()}, Total counts: {df_sizes['Total_UMI_Counts'].sum():,}")
    print(df_sizes.to_string(index=False))

print_lib_sizes("Sensitivity (g) Ceiling-Rescued Stromal", df_pb_g, donor_meta_g)
print_lib_sizes("Sensitivity (h) High-Confidence Donors (>=50 cells & med genes >=1000)", df_pb_h, donor_meta_h)

print("\nPre-specification (self-documented) and preparation of (g) and (h) completed successfully.")
print("=" * 80)
