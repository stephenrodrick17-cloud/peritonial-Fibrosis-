"""
scripts/26_stage3_repair_hub_blinding_and_harmony.py
Stage 3 Repair: Strict Hub-Gene Sealing, Pre-QC Scrublet, Harmony Integration,
and Refined Lineage Annotation.
"""
import os
import sys
import gzip
import json
import time
import hashlib
import subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import anndata as ad
import scanpy as sc
from sklearn.metrics import silhouette_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("scripts")
from guard import guard_check

print("=================================================================")
print("SCRIPT 26: STAGE 3 REPAIR - REAL BLINDING & HARMONY PIPELINE")
print("=================================================================")

# Pre-defined 11 Hub Genes to be permanently sealed
HUB_GENES = [
    "ISM1", "FN1", "EDIL3", "VCAN", "COL3A1",
    "COMP", "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"
]

os.makedirs("sealed", exist_ok=True)
os.makedirs("data/processed/GSE248762_blinded_samples", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# Guard check: Ensure guard is active
guard_check([], stage="repair_start")

# =====================================================================
# STEP A1 & B2: RAW MTX LOAD, HUB SEALING, AND PRE-QC SCRUBLET
# =====================================================================
raw_dir = "data/raw/GSE248762_extracted"
with open("provenance/GSE248762_acquisition_summary.json", "r") as f:
    acq_data = json.load(f)
manifest_samples = {s["GSM"]: s for s in acq_data["samples"]}

barcode_files = sorted(glob.glob(os.path.join(raw_dir, "*barcodes.tsv.gz"))) if "glob" in locals() else []
import glob
barcode_files = sorted(glob.glob(os.path.join(raw_dir, "*barcodes.tsv.gz")))

# Helper function for MAD
def calc_mad(data):
    med = np.median(data)
    return np.median(np.abs(data - med))

sealed_hub_adatas = []
filtered_blind_adatas = []
scrublet_summary_rows = []
qc_summary_rows = []

np.random.seed(42)
sc.settings.seed = 42

# Doublet score plots
fig_scrub, axes_scrub = plt.subplots(4, 4, figsize=(20, 16))
axes_scrub = axes_scrub.flatten()

print(f"\nProcessing {len(barcode_files)} raw donor libraries...")

for idx, bf in enumerate(barcode_files):
    fname = os.path.basename(bf)
    gsm = fname.split("_")[0]
    sample_title = fname.split(".")[0].split("_", 1)[1]
    
    m_info = manifest_samples[gsm]
    grp = m_info["Group"].split(" ")[0]
    
    feat_file = os.path.join(raw_dir, m_info["Files"]["features"]["filename"])
    mtx_file = os.path.join(raw_dir, m_info["Files"]["matrix"]["filename"])
    
    # 1. Load barcodes, features, matrix
    with gzip.open(bf, 'rt') as f:
        barcodes = [line.strip() for line in f if line.strip()]
        
    features_df = pd.read_csv(feat_file, sep='\t', header=None, compression='gzip')
    gene_ids = features_df[0].values
    gene_symbols = features_df[1].values
    
    mat = scipy.io.mmread(mtx_file).T.tocsr() # cells x genes
    
    # Check gene symbols for hub genes
    hub_mask = np.isin(gene_symbols, HUB_GENES)
    non_hub_mask = ~hub_mask
    
    # Extract hub counts to sealed object
    hub_mat = mat[:, hub_mask]
    hub_symbols = gene_symbols[hub_mask]
    hub_ids = gene_ids[hub_mask]
    
    adata_hub = ad.AnnData(
        X=hub_mat,
        obs=pd.DataFrame(index=[f"{gsm}_{b}" for b in barcodes]),
        var=pd.DataFrame({"gene_ids": hub_ids, "gene_symbols": hub_symbols}, index=hub_symbols)
    )
    adata_hub.obs["gsm"] = gsm
    adata_hub.obs["donor_id"] = sample_title
    adata_hub.obs["group"] = grp
    sealed_hub_adatas.append(adata_hub)
    
    # Create Hub-Blind AnnData
    blind_mat = mat[:, non_hub_mask]
    blind_symbols = gene_symbols[non_hub_mask]
    blind_ids = gene_ids[non_hub_mask]
    
    adata_blind = ad.AnnData(
        X=blind_mat,
        obs=pd.DataFrame(index=[f"{gsm}_{b}" for b in barcodes]),
        var=pd.DataFrame({"gene_ids": blind_ids, "gene_symbols": blind_symbols}, index=blind_symbols)
    )
    adata_blind.var_names_make_unique()
    adata_blind.obs["gsm"] = gsm
    adata_blind.obs["donor_id"] = sample_title
    adata_blind.obs["group"] = grp
    
    raw_cells = adata_blind.n_obs
    
    # 2. Run Scrublet on RAW cells BEFORE QC filtering (per plan)
    # Expected doublet rate 0.06
    sc.pp.scrublet(adata_blind, expected_doublet_rate=0.06, random_state=42)
    
    d_scores = adata_blind.obs["doublet_score"].values
    d_pred = adata_blind.obs["predicted_doublet"].values
    n_doublets = int(d_pred.sum())
    
    # Scrublet threshold
    threshold = getattr(adata_blind.uns.get("scrublet", {}), "threshold_", None)
    if threshold is None and "threshold" in adata_blind.uns.get("scrublet", {}):
        threshold = adata_blind.uns["scrublet"]["threshold"]
    if threshold is None:
        threshold = float(np.min(d_scores[d_pred])) if n_doublets > 0 else 0.25
        
    scrublet_summary_rows.append({
        "GSM": gsm,
        "Sample_Title": sample_title,
        "Group": grp,
        "Raw_Cells": raw_cells,
        "Scrublet_Threshold": round(float(threshold), 4),
        "Predicted_Doublets": n_doublets,
        "Predicted_Doublet_Pct": round(n_doublets / raw_cells * 100, 2),
        "Threshold_Rule": "Scanpy scrublet automatic bimodal simulation threshold"
    })
    
    # Plot histogram on ax
    ax = axes_scrub[idx]
    ax.hist(d_scores, bins=40, color="#3b82f6", alpha=0.7, density=True, label="Observed")
    ax.axvline(threshold, color="#dc2626", linestyle="--", linewidth=1.5, label=f"Cut: {threshold:.3f}")
    ax.set_title(f"{sample_title} ({n_doublets} doublets)", fontsize=10, fontweight="bold")
    if idx == 0:
        ax.legend(fontsize=8)
        
    # 3. Calculate Cell Quality Metrics on Hub-Blind Object
    n_counts = np.array(adata_blind.X.sum(axis=1)).flatten()
    n_genes = np.array((adata_blind.X > 0).sum(axis=1)).flatten()
    
    mito_mask = adata_blind.var_names.str.upper().str.startswith("MT-")
    mito_counts = np.array(adata_blind[:, mito_mask].X.sum(axis=1)).flatten()
    pct_mito = (mito_counts / (n_counts + 1e-9)) * 100.0
    
    adata_blind.obs["n_counts"] = n_counts
    adata_blind.obs["n_genes"] = n_genes
    adata_blind.obs["pct_counts_mt"] = pct_mito
    
    # 3 MAD adaptive thresholds on log10 counts & log10 genes
    log_c = np.log10(n_counts + 1)
    med_c = np.median(log_c)
    mad_c = calc_mad(log_c)
    c_low = max(500.0, 10 ** (med_c - 3 * mad_c))
    c_high = 10 ** (med_c + 3 * mad_c)
    
    log_g = np.log10(n_genes + 1)
    med_g = np.median(log_g)
    mad_g = calc_mad(log_g)
    g_low = max(200.0, 10 ** (med_g - 3 * mad_g))
    g_high = 10 ** (med_g + 3 * mad_g)
    
    mito_ceiling = 15.0
    
    pass_qc = (n_counts >= c_low) & (n_counts <= c_high) & (n_genes >= g_low) & (n_genes <= g_high) & (pct_mito <= mito_ceiling)
    keep_cell = pass_qc & (~d_pred)
    
    adata_filtered = adata_blind[keep_cell].copy()
    post_cells = adata_filtered.n_obs
    
    qc_summary_rows.append({
        "GSM": gsm,
        "Sample_Title": sample_title,
        "Group": grp,
        "Raw_Cells": raw_cells,
        "Post_QC_Cells": post_cells,
        "Retained_Pct": round(post_cells / raw_cells * 100, 2),
        "UMI_Low_MAD": round(float(c_low), 4),
        "UMI_High_MAD": round(float(c_high), 4),
        "Gene_Low_MAD": round(float(g_low), 4),
        "Gene_High_MAD": round(float(g_high), 4),
        "Mito_Ceiling_Pct": mito_ceiling,
        "Doublets_Removed_PreQC": n_doublets,
        "Doublets_Removed_In_Filtered": int((pass_qc & d_pred).sum())
    })
    
    filtered_blind_adatas.append(adata_filtered)
    print(f"  {gsm} ({sample_title}): Raw {raw_cells:5d} -> Doublets {n_doublets:4d} -> Post-QC {post_cells:5d} ({post_cells/raw_cells*100:5.1f}%)")

plt.tight_layout()
fig_scrub_path = "results/figures/B2_scrublet_doublet_score_histograms.png"
plt.savefig(fig_scrub_path, dpi=300)
plt.close()
print(f"Saved Scrublet doublet score histograms to {fig_scrub_path}")

# Save Scrublet Summary Table
df_scrublet = pd.DataFrame(scrublet_summary_rows)
df_scrublet.to_csv("results/tables/B2_scrublet_summary.csv", index=False)
print("\n--- B2 SCRUBLET SUMMARY TABLE ---")
print(df_scrublet.to_string(index=False))

# Save E2 Filtering Summary Table
df_qc_summary = pd.DataFrame(qc_summary_rows)
df_qc_summary.to_csv("results/tables/E2_filtering_summary.csv", index=False)
print("\n--- E2 FILTERING SUMMARY TABLE ---")
print(df_qc_summary.to_string(index=False))

# Save Sealed Hub Genes
print("\nSealing 11 Hub Genes into sealed/hub_counts.h5ad...")
adata_sealed_cohort = ad.concat(sealed_hub_adatas, join="outer", fill_value=0.0)
adata_sealed_cohort.write_h5ad("sealed/hub_counts.h5ad")

with open("sealed/hub_counts.h5ad", "rb") as f:
    sealed_sha = hashlib.sha256(f.read()).hexdigest()
sealed_mtime = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(os.path.getmtime("sealed/hub_counts.h5ad")))

print(f"Sealed file created: sealed/hub_counts.h5ad")
print(f"  Shape: {adata_sealed_cohort.shape} (Cells x 11 Hub Genes)")
print(f"  SHA256: {sealed_sha}")
print(f"  Timestamp: {sealed_mtime}")

# Concatenate Hub-Blind Cohort
print("\nConcatenating Hub-Blind Cohort...")
adata = ad.concat(filtered_blind_adatas, join="outer", fill_value=0.0)
adata.var_names_make_unique()
print(f"Hub-Blind Cohort: {adata.n_obs:,} cells x {adata.n_vars:,} genes across {adata.obs['donor_id'].nunique()} donors.")

# Verify zero hub genes present
guard_check(adata.var_names.tolist(), stage="verify_blinded_anndata")
print("Confirmed: Exactly 0 of the 11 hub genes are present in the processing AnnData object.")

# =====================================================================
# STEP B1: NORMALIZATION, HVG, PCA & HARMONY INTEGRATION
# =====================================================================
print("\n--- STEP B1: Normalization, HVG, PCA & Harmony Batch Integration ---")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

# HVG selection with batch balancing
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="donor_id", subset=False)
guard_check(adata.var_names[adata.var["highly_variable"]].tolist(), stage="verify_blinded_hvg")
print(f"HVG selection complete: {adata.var['highly_variable'].sum()} HVGs selected (0 hub genes).")

# PCA (30 PCs)
sc.pp.pca(adata, n_comps=30, use_highly_variable=True, random_state=42)

# Save unintegrated PCA to temp CSV for Harmony in R
pca_unintegrated = adata.obsm["X_pca"]
df_pca = pd.DataFrame(pca_unintegrated, index=adata.obs_names)
df_pca.to_csv("sealed/temp_pca_unintegrated.csv")
adata.obs[["donor_id"]].to_csv("sealed/temp_pca_metadata.csv")

# Run Harmony via R's harmony package (harmony::HarmonyMatrix)
r_harmony_script = """
library(harmony)
pca_mat <- as.matrix(read.csv('sealed/temp_pca_unintegrated.csv', row.names=1))
meta <- read.csv('sealed/temp_pca_metadata.csv', row.names=1)
set.seed(42)
harmony_embed <- HarmonyMatrix(
  data_mat = pca_mat,
  meta_data = meta,
  vars_use = 'donor_id',
  do_pca = FALSE,
  max.iter.harmony = 20,
  verbose = TRUE
)
write.csv(harmony_embed, 'sealed/temp_pca_harmony.csv')
"""
with open("sealed/run_harmony.R", "w") as f:
    f.write(r_harmony_script)

print("Executing Harmony matrix batch correction via R harmony::HarmonyMatrix...")
r_cmd = ["C:\\Program Files\\R\\R-4.4.2\\bin\\Rscript.exe", "sealed/run_harmony.R"]
res = subprocess.run(r_cmd, capture_output=True, text=True)
print(res.stdout)
if res.returncode != 0:
    print("R Harmony error:", res.stderr)
    raise RuntimeError("Harmony integration failed in R.")

harmony_df = pd.read_csv("sealed/temp_pca_harmony.csv", index=1)
adata.obsm["X_pca_harmony"] = harmony_df.values

# Clean temp files
for tf in ["sealed/temp_pca_unintegrated.csv", "sealed/temp_pca_metadata.csv", "sealed/run_harmony.R", "sealed/temp_pca_harmony.csv"]:
    if os.path.exists(tf): os.remove(tf)

# Neighborhood graph & UMAP on Harmony embedding
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony", random_state=42)
sc.tl.umap(adata, random_state=42)
print("Harmony integration complete and stored in adata.obsm['X_pca_harmony'].")

# =====================================================================
# STEP B3: LEIDEN RESOLUTION SELECTION (COMPUTED FROM DATA)
# =====================================================================
print("\n--- STEP B3: Leiden Resolution Selection via Silhouette Scores ---")
resolutions = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
res_data = []

# Subsample 5,000 cells for silhouette evaluation with fixed seed 42
eval_idx = np.random.RandomState(42).choice(adata.n_obs, min(5000, adata.n_obs), replace=False)
emb_eval = adata.obsm["X_pca_harmony"][eval_idx]

best_res = 0.3
best_sil = -1.0

for r in resolutions:
    key = f"leiden_res_{r:.1f}"
    sc.tl.leiden(adata, resolution=r, key_added=key, random_state=42, flavor="igraph", n_iterations=2, directed=False)
    n_cl = int(adata.obs[key].nunique())
    labels_eval = adata.obs[key].values[eval_idx]
    sil = float(silhouette_score(emb_eval, labels_eval))
    res_data.append({
        "Resolution": r,
        "N_Clusters": n_cl,
        "Silhouette_Score": round(sil, 4),
        "Subsample_Size": len(eval_idx),
        "Evaluation_Seed": 42
    })
    print(f"  Resolution {r:.1f} -> {n_cl:2d} clusters | Silhouette: {sil:+.4f}")
    if sil > best_sil:
        best_sil = sil
        best_res = r

df_leiden_eval = pd.DataFrame(res_data)
df_leiden_eval.to_csv("results/tables/B3_leiden_resolution_selection.csv", index=False)
print(f"\nOptimal resolution chosen: {best_res:.1f} (Silhouette: {best_sil:+.4f})")
adata.obs["leiden"] = adata.obs[f"leiden_res_{best_res:.1f}"]

# =====================================================================
# STEP B4: HUB-BLIND ANNOTATION, WILCOXON RANK GENES, EXPANDED MODULES
# =====================================================================
print("\n--- STEP B4: Differential Expression & Lineage Annotation ---")
# 1. Wilcoxon rank_genes_groups (top 15 genes per cluster) on hub-blind object
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="rank_genes_wilcoxon")

top_markers = pd.DataFrame(adata.uns["rank_genes_wilcoxon"]["names"]).head(15)
top_markers.to_csv("results/tables/B4_top15_cluster_markers_wilcoxon.csv")
print("Saved results/tables/B4_top15_cluster_markers_wilcoxon.csv (Top 15 markers per cluster)")

# 2. Lineage module definitions
lineage_modules = {
    "T_Cell": ["CD3D", "CD3E", "CD4", "CD8A", "TRAC"],
    "NK_Cell": ["NCAM1", "NKG7", "GNLY", "KLRD1"],
    "B_Cell": ["CD19", "MS4A1", "CD79A"],
    "Plasma_Cell": ["JCHAIN", "MZB1", "SDC1"],
    "Monocyte_Macrophage": ["C1QA", "C1QB", "MARCO", "CD163", "CD68", "CSF1R"],
    "Conventional_DC": ["CD1C", "CLEC10A", "FCER1A", "CLEC9A"],
    "Neutrophil": ["FCGR3B", "S100A8", "S100A9", "CSF3R", "CXCR2"],
    "Stromal_Mesothelial": ["DCN", "LUM", "PDGFRA", "COL1A2", "ACTA2", "TAGLN", "WT1", "MSLN", "CALB2", "UPK3B", "KRT19", "UPK1B"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"]
}

# Guard check on lineage module genes
all_module_genes = [g for sub in lineage_modules.values() for g in sub]
guard_check(all_module_genes, stage="verify_lineage_module_genes")

# Co-expression test for EMT vs separate mesothelial lineage in stromal cluster
meso_markers = ["WT1", "MSLN", "UPK3B", "CALB2", "KRT19"]
fibro_markers = ["PDGFRA", "DCN", "LUM", "COL1A2"]

present_meso = [g for g in meso_markers if g in adata.raw.var_names]
present_fibro = [g for g in fibro_markers if g in adata.raw.var_names]

meso_expr = np.array(adata.raw[:, present_meso].X.mean(axis=1)).flatten()
fibro_expr = np.array(adata.raw[:, present_fibro].X.mean(axis=1)).flatten()

co_expr_mask = (meso_expr > 0.1) & (fibro_expr > 0.1)
co_expr_frac = float(np.sum(co_expr_mask) / adata.n_obs)
print(f"Stromal/Mesothelial Co-expression test: {np.sum(co_expr_mask)} / {adata.n_obs} cells ({co_expr_frac*100:.2f}%) co-express both modules.")

# Cluster annotation mapping
cluster_annotations = {}
cluster_rows = []
clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))

for cl in clusters:
    cl_mask = adata.obs["leiden"] == cl
    n_cl = int(cl_mask.sum())
    
    # Check for small clusters (<50 cells)
    if n_cl < 50:
        donor_counts = adata.obs[cl_mask]["donor_id"].value_counts().to_dict()
        top_donor, top_d_cnt = list(donor_counts.items())[0]
        assigned_type = f"Small/Unassigned (Cluster {cl}, N={n_cl}, {top_d_cnt}/{n_cl} from {top_donor})"
        top_mod = "None"
        top_score = 0.0
    else:
        scores = {}
        for mod_name, mod_genes in lineage_modules.items():
            valid_g = [g for g in mod_genes if g in adata.raw.var_names]
            if valid_g:
                sub_mat = adata.raw[cl_mask, valid_g].X
                scores[mod_name] = float(np.mean(sub_mat.toarray() if sp.issparse(sub_mat) else sub_mat))
            else:
                scores[mod_name] = 0.0
                
        sorted_s = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        top_mod, top_score = sorted_s[0]
        
        if top_mod == "T_Cell":
            assigned_type = "T Cell"
        elif top_mod == "NK_Cell":
            assigned_type = "NK Cell"
        elif top_mod == "B_Cell":
            assigned_type = "B Cell"
        elif top_mod == "Plasma_Cell":
            assigned_type = "Plasma Cell"
        elif top_mod == "Monocyte_Macrophage":
            assigned_type = "Monocyte / Macrophage"
        elif top_mod == "Conventional_DC":
            assigned_type = "Conventional Dendritic Cell"
        elif top_mod == "Neutrophil":
            assigned_type = "Neutrophil"
        elif top_mod == "Stromal_Mesothelial":
            assigned_type = "Stromal / mesothelial-lineage (unresolved)"
        elif top_mod == "Endothelial":
            assigned_type = "Endothelial Cell"
        else:
            assigned_type = "Unassigned"
            
    cluster_annotations[cl] = assigned_type
    cluster_rows.append({
        "Cluster": cl,
        "N_Cells": n_cl,
        "Assigned_Lineage": assigned_type,
        "Top_Module": top_mod,
        "Top_Module_Score": round(top_score, 4)
    })
    print(f"  Cluster {cl:>2s} ({n_cl:5d} cells) -> {assigned_type}")

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_annotations)

df_comp = pd.DataFrame(cluster_rows)
df_comp.to_csv("results/tables/E3_cell_type_composition.csv", index=False)

# =====================================================================
# STEP B5: BATCH MIXING METRICS & CLUSTER X DONOR DIVERSITY
# =====================================================================
print("\n--- STEP B5: Batch Mixing Metrics Before & After Harmony ---")
sil_donor_unint = float(silhouette_score(adata.obsm["X_pca"][eval_idx], adata.obs["donor_id"].values[eval_idx]))
sil_donor_harm = float(silhouette_score(adata.obsm["X_pca_harmony"][eval_idx], adata.obs["donor_id"].values[eval_idx]))

sil_ct_unint = float(silhouette_score(adata.obsm["X_pca"][eval_idx], adata.obs["cell_type"].values[eval_idx]))
sil_ct_harm = float(silhouette_score(adata.obsm["X_pca_harmony"][eval_idx], adata.obs["cell_type"].values[eval_idx]))

df_batch_metrics = pd.DataFrame([
    {"Metric": "Donor Silhouette (Batch Effect, lower is better)", "Unintegrated_PCA": round(sil_donor_unint, 4), "Harmony_Integrated": round(sil_donor_harm, 4), "Interpretation": "Negative indicates donor mixing across clusters"},
    {"Metric": "Cell Type Silhouette (Lineage separation, higher is better)", "Unintegrated_PCA": round(sil_ct_unint, 4), "Harmony_Integrated": round(sil_ct_harm, 4), "Interpretation": "Positive indicates distinct lineage clusters"}
])
df_batch_metrics.to_csv("results/tables/E3_batch_mixing_metrics.csv", index=False)
print(df_batch_metrics.to_string(index=False))

# Cluster x Donor diversity table
cl_donor_ct = pd.crosstab(adata.obs["leiden"], adata.obs["donor_id"])
diversity_rows = []
for cl in clusters:
    cl_counts = cl_donor_ct.loc[cl]
    total_cl = cl_counts.sum()
    fracs = cl_counts / total_cl
    n_ge5 = int((fracs >= 0.05).sum())
    top_d = fracs.idxmax()
    max_f = float(fracs.max())
    diversity_rows.append({
        "Cluster": cl,
        "N_Cells": total_cl,
        "Assigned_Lineage": cluster_annotations[cl],
        "Donors_ge_5pct": f"{n_ge5}/16",
        "Top_Contributing_Donor": top_d,
        "Max_Donor_Fraction": round(max_f, 4)
    })

df_cluster_donor_div = pd.DataFrame(diversity_rows)
df_cluster_donor_div.to_csv("results/tables/B5_cluster_donor_diversity.csv", index=False)
print("\n--- B5 CLUSTER X DONOR DIVERSITY TABLE ---")
print(df_cluster_donor_div.to_string(index=False))

# =====================================================================
# STEP B6: DONOR X CELL TYPE TABLE & MINIMUM EVIDENCE TESTABILITY
# =====================================================================
print("\n--- STEP B6: Donor x Cell Type Contingency Table & Pre-Specified (Self-Documented) Testability ---")
ct_donor = pd.crosstab(adata.obs["donor_id"], adata.obs["cell_type"])
ct_donor["Group"] = ct_donor.index.map(adata.obs.groupby("donor_id")["group"].first().to_dict())
df_donor_ct_final = ct_donor.reset_index()
df_donor_ct_final.to_csv("results/tables/E4_donor_cell_counts_by_type.csv", index=False)
print("\n--- E4 DONOR X CELL TYPE TABLE (RAW DATAFRAME) ---")
print(df_donor_ct_final.to_string(index=False))

# Evaluate testability rule:
# Minimum evidence rule: >=20 cells in >=3 donors in target group AND >=2 donors in comparator group
contrasts = [
    {"name": "LV_UF_vs_LV_NOT_UF", "target": "LV_UF", "comp": "LV_NOT_UF", "min_t": 3, "min_c": 2},
    {"name": "LV_UF_vs_SV", "target": "LV_UF", "comp": "SV", "min_t": 3, "min_c": 2},
    {"name": "LV_NOT_UF_vs_SV", "target": "LV_NOT_UF", "comp": "SV", "min_t": 3, "min_c": 2}
]

eval_cell_types = [c for c in ct_donor.columns if c not in ["donor_id", "Group"] and not c.startswith("Small")]

testability_rows = []
for c_info in contrasts:
    c_name = c_info["name"]
    tgt_g = c_info["target"]
    cmp_g = c_info["comp"]
    
    tgt_df = ct_donor[ct_donor["Group"] == tgt_g]
    cmp_df = ct_donor[ct_donor["Group"] == cmp_g]
    
    for ct in eval_cell_types:
        t_ge20 = int((tgt_df[ct] >= 20).sum())
        c_ge20 = int((cmp_df[ct] >= 20).sum())
        
        is_testable = (t_ge20 >= c_info["min_t"]) and (c_ge20 >= c_info["min_c"])
        if is_testable:
            status = "PASS (Testable)"
        else:
            reasons = []
            if t_ge20 < c_info["min_t"]:
                reasons.append(f"Target {tgt_g} donors ge20: {t_ge20}/{len(tgt_df)} (need {c_info['min_t']})")
            if c_ge20 < c_info["min_c"]:
                reasons.append(f"Comparator {cmp_g} donors ge20: {c_ge20}/{len(cmp_df)} (need {c_info['min_c']})")
            status = f"FAIL ({'; '.join(reasons)})"
            
        testability_rows.append({
            "Contrast": c_name,
            "Cell_Type": ct,
            "Target_Group": tgt_g,
            "Target_Donors_ge20": f"{t_ge20}/{len(tgt_df)}",
            "Comparator_Group": cmp_g,
            "Comparator_Donors_ge20": f"{c_ge20}/{len(cmp_df)}",
            "Is_Testable": is_testable,
            "Testability_Status": status
        })

df_testability = pd.DataFrame(testability_rows)
df_testability.to_csv("results/tables/E4_celltype_testability_table.csv", index=False)
print("\n--- E4 TESTABILITY TABLE (RAW DATAFRAME) ---")
print(df_testability.to_string(index=False))

# Update plan file with dated deviations
plan_path = "provenance/analysis_plan_scRNA.json"
with open(plan_path, "r") as f:
    plan = json.load(f)

plan["repair_stage3_deviations"] = {
    "timestamp_utc": time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    "harmony_integration": "Harmony integration completed via R harmony::HarmonyMatrix on PCA (30 PCs, vars_use=donor_id, seed=42)",
    "scrublet_timing": "Scrublet run on raw counts per sample BEFORE QC filtering",
    "annotation_refinements": [
        "Added Plasma Cell module (JCHAIN, MZB1, SDC1)",
        "Separated Monocyte/Macrophage (C1QA, C1QB, MARCO, CD163, CD68, CSF1R) from cDC (CD1C, CLEC10A, FCER1A, CLEC9A)",
        "Relabeled stromal cluster as 'Stromal / mesothelial-lineage (unresolved)'",
        "Labeled clusters with <50 cells as 'Small/Unassigned' and excluded from downstream testing"
    ]
}
with open(plan_path, "w") as f:
    json.dump(plan, f, indent=2)

# Save repaired annotated AnnData
adata.write_h5ad("data/processed/GSE248762_annotated.h5ad")
print("\nSaved repaired cohort to data/processed/GSE248762_annotated.h5ad")
print("Script 26 execution complete.")
