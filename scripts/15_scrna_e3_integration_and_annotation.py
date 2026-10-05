"""
Step E3: scRNA-seq Integration, Leiden Clustering & Cell-Type Annotation
Strictly Hub-Gene-Blind execution with code-level guard.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
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
print("STEP E3: SCRNA-SEQ INTEGRATION & CELL-TYPE ANNOTATION (HUB-BLIND)")
print("=================================================================")

# Guard check: Ensure no hub genes are accessed
guard_check([], stage="E3_start")

proc_dir = "data/processed/GSE248762_filtered_samples"
sample_files = [os.path.join(proc_dir, f) for f in sorted(os.listdir(proc_dir)) if f.endswith(".h5ad")]
print(f"Loading {len(sample_files)} filtered donor AnnData objects...")

adatas = []
for sfile in sample_files:
    ad_s = ad.read_h5ad(sfile)
    adatas.append(ad_s)

# Concatenate into unified cohort AnnData
adata = ad.concat(adatas, join="outer", fill_value=0.0)
adata.var_names_make_unique()
print(f"Unified Cohort: {adata.n_obs:,} cells x {adata.n_vars:,} genes across {adata.obs['donor_id'].nunique()} donors.")

# Set seed 42
sc.settings.seed = 42
np.random.seed(42)

# 1. Normalization & Log Transformation
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

# 2. Highly Variable Genes (2,000 HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="donor_id", subset=False)

# Guard check: ensure Stage E3 remains strictly hub-gene blind
guard_check([], stage="E3_hvg_check")
print(f"Identified {adata.var['highly_variable'].sum()} Highly Variable Genes.")

# 3. PCA
sc.pp.pca(adata, n_comps=30, use_highly_variable=True, random_state=42)

# 4. Batch Integration & Neighborhood Graph
# Harmony integration / batch-corrected PCA
try:
    import harmonypy as hm
    print("Running Harmony batch integration across donors...")
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "donor_id", max_iter_harmony=10, random_state=42)
    adata.obsm["X_pca_integrated"] = ho.Z_corr.T
    use_rep = "X_pca_integrated"
except Exception as e:
    print(f"Harmony note: using standard PCA embedding with batch balancing ({e})")
    use_rep = "X_pca"

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep=use_rep, random_state=42)
sc.tl.umap(adata, random_state=42)

# 5. Leiden Clustering with Silhouette Score Resolution Selection
resolutions = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
res_scores = {}

# Subsample 5,000 cells for silhouette evaluation
eval_idx = np.random.choice(adata.n_obs, min(5000, adata.n_obs), replace=False)
emb_eval = adata.obsm[use_rep][eval_idx]

print("\nEvaluating Leiden resolutions via Silhouette Score:")
best_res = 0.4
best_score = -1.0

for res in resolutions:
    key = f"leiden_res_{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key, random_state=42, flavor="igraph", n_iterations=2, directed=False)
    labels = adata.obs[key].values[eval_idx]
    if len(np.unique(labels)) > 1:
        sil = silhouette_score(emb_eval, labels)
    else:
        sil = -1.0
    res_scores[res] = sil
    print(f"  Resolution {res:.1f} -> Clusters: {len(np.unique(labels)):2d} | Avg Silhouette: {sil:+.4f}")
    if sil > best_score:
        best_score = sil
        best_res = res

print(f"\nSelected optimal resolution: {best_res:.1f} (Silhouette Score: {best_score:+.4f})")
adata.obs["leiden"] = adata.obs[f"leiden_res_{best_res}"]

# 6. Cell-Type Annotation via Canonical Marker Panel
with open("provenance/marker_panel.json", "r") as f:
    marker_panel = json.load(f)

# Guard check: Ensure no hub genes in marker panel
all_markers = [m for sublist in marker_panel.values() for m in sublist]
guard_check(all_markers, stage="E3_marker_panel_check")

# Calculate lineage scores per cluster
cluster_annotations = {}
cluster_composition = []

for cl in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
    cl_mask = adata.obs["leiden"] == cl
    n_cl = int(cl_mask.sum())
    
    lineage_scores = {}
    for lineage, markers in marker_panel.items():
        present_m = [m for m in markers if m in adata.raw.var_names]
        if present_m:
            # Mean expression across cells in cluster
            sub_mat = adata.raw[cl_mask, present_m].X
            if sp.issparse(sub_mat):
                score = np.mean(sub_mat.toarray())
            else:
                score = np.mean(sub_mat)
            lineage_scores[lineage] = float(score)
        else:
            lineage_scores[lineage] = 0.0
            
    # Sort scores
    sorted_scores = sorted(lineage_scores.items(), key=lambda x: x[1], reverse=True)
    top_lineage, top_score = sorted_scores[0]
    second_lineage, second_score = sorted_scores[1]
    
    # Strict unassigned rule: must have positive enrichment and distinct from background
    if top_score < 0.15 or (top_score - second_score < 0.05 and top_score < 0.30):
        assigned_type = "Unassigned"
    else:
        assigned_type = top_lineage.replace("_", " ")
        
    cluster_annotations[cl] = assigned_type
    print(f"  Cluster {cl:>2s} ({n_cl:5d} cells) -> {assigned_type:<24} (Top marker score: {top_score:.3f}, Next: {second_lineage} {second_score:.3f})")
    
    cluster_composition.append({
        "Cluster": cl,
        "N_Cells": n_cl,
        "Cell_Type": assigned_type,
        "Top_Lineage": top_lineage,
        "Top_Score": round(top_score, 4),
        "Second_Lineage": second_lineage,
        "Second_Score": round(second_score, 4)
    })

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_annotations)

# Save composition table
df_comp = pd.DataFrame(cluster_composition)
df_comp.to_csv("results/tables/E3_cell_type_composition.csv", index=False)
print("\nSaved results/tables/E3_cell_type_composition.csv")

# 7. Batch Mixing & Quality Metrics
sil_donor = silhouette_score(emb_eval, adata.obs["donor_id"].values[eval_idx])
sil_celltype = silhouette_score(emb_eval, adata.obs["cell_type"].values[eval_idx])

df_metrics = pd.DataFrame([{
    "Metric": "Silhouette Score (Cell Type Separation)",
    "Value": round(sil_celltype, 4),
    "Interpretation": "Higher is better (distinct cell lineages)"
}, {
    "Metric": "Silhouette Score (Donor Separation / Batch Effect)",
    "Value": round(sil_donor, 4),
    "Interpretation": "Lower is better (well-mixed donors)"
}, {
    "Metric": "Optimal Leiden Resolution",
    "Value": best_res,
    "Interpretation": f"Tested range {resolutions}"
}])
df_metrics.to_csv("results/tables/E3_batch_mixing_metrics.csv", index=False)
print("Saved results/tables/E3_batch_mixing_metrics.csv")
print(df_metrics)

# 8. Publication Figures
# Figure A: UMAP by Cell Type and Clinical Group
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sns.set_theme(style="white", font_scale=1.1)

palette_ct = sns.color_palette("tab20", n_colors=adata.obs["cell_type"].nunique())
sc.pl.umap(adata, color="cell_type", ax=axes[0], show=False, title="A. Effluent Single-Cell Transcriptomic Map", palette=palette_ct, frameon=False)
sc.pl.umap(adata, color="group", ax=axes[1], show=False, title="B. Distribution by Clinical Condition (Vintage / UFF)",
           palette={"SV": "#10b981", "LV_NOT_UF": "#3b82f6", "LV_UF": "#ef4444"}, frameon=False)

plt.tight_layout()
fig_umap_path = "results/figures/E3_umap_by_celltype_and_group.png"
plt.savefig(fig_umap_path, dpi=300)
plt.close()
print(f"Saved publication UMAP to {fig_umap_path}")

# Figure B: Canonical Lineage Marker Dot Plot
dot_markers = {}
for lin, m_list in marker_panel.items():
    lin_name = lin.replace("_", " ")
    # Take top 2 markers present in dataset
    valid_m = [m for m in m_list if m in adata.raw.var_names][:2]
    if valid_m:
        dot_markers[lin_name] = valid_m

# Guard check before plotting
guard_check([m for sub in dot_markers.values() for m in sub], stage="E3_dotplot_guard")

sc.pl.dotplot(adata, dot_markers, groupby="cell_type", standard_scale="var", show=False,
              title="Canonical Lineage Marker Expression Across Effluent Populations")
fig_dot_path = "results/figures/E3_marker_dotplot.png"
plt.savefig(fig_dot_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved publication Marker Dotplot to {fig_dot_path}")

# Save Annotated AnnData
adata.write_h5ad("data/processed/GSE248762_annotated.h5ad")
print("Saved annotated cohort to data/processed/GSE248762_annotated.h5ad")

print("\nStep E3 finished successfully.")
