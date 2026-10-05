"""
scripts/30b_scrna_harmony_annotation_testability.py
Stage 3 repair (round 2), part B. Input: data/processed/GSE248762_hubblind_allcells_qc.h5ad (script 29).

  1. Regenerates E2_filtering_summary.csv FROM the h5ad obs (single source of truth).
  2. LV_UF-3: gene/UMI counts by cell type; marker make-up of QC-removed vs kept barcodes.
  3. Normalise -> HVG (batch-aware) -> PCA(30) -> UMAP (unintegrated) ->
     Harmony by donor (harmonypy 0.0.10) -> neighbours/UMAP -> Leiden (resolution chosen by the
     frozen plan rule: max silhouette over 0.3-0.8 on the Harmony embedding).
  4. Annotation with revised modules (Mono/Mac vs cDC separated, plasma, pDC, mast; endothelial
     = VWF/CDH5/KDR, PECAM1 dropped; stromal renamed "stromal / mesothelial-lineage (unresolved)").
  5. Neutrophil confirmation (FCGR3B, CSF3R, CXCR2, MNDA, G0S2).
  6. Per-cluster donor fractions; donor-dominated flag (largest donor fraction >= 0.50);
     testability (>=20 cells in >=3 target donors and >=2 comparator donors) before / after removal.
  7. Figures: UMAPs by donor & group before/after Harmony; dot plots (unscaled + scaled).
Hub genes are absent from the object (removed at load in script 29) and every gene list is guard-checked.
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import harmonypy
from sklearn.metrics import silhouette_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append("scripts")
from guard import guard_check, BLOCKED_HUB_GENES

pd.set_option("display.width", 260)
pd.set_option("display.max_columns", 60)
pd.set_option("display.max_rows", 500)
pd.set_option("display.max_colwidth", 120)
np.random.seed(0)
sc.settings.verbosity = 1
T0 = time.time()


def tick(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


IN = "data/processed/GSE248762_hubblind_allcells_qc.h5ad"
print("=" * 78)
print("SCRIPT 30b: E2 from h5ad, Harmony, annotation, donor fractions, testability")
print("=" * 78)
guard_check([], stage="30b_start")

full = ad.read_h5ad(IN)
guard_check(full.var_names.tolist(), stage="30b_load")
assert len(set(full.var_names.str.upper()) & BLOCKED_HUB_GENES) == 0
ORDER = list(dict.fromkeys(full.obs.sort_values("gsm")["donor_id"]))
GROUP_OF = full.obs.groupby("donor_id", observed=True)["group"].first().to_dict()
tick(f"loaded {full.n_obs} barcodes x {full.n_vars} genes; doublet rule: {full.uns['doublet_rule']}")

# ---------------------------------------------------------------------------
# 1. E2 regenerated from the h5ad
# ---------------------------------------------------------------------------
o = full.obs
e2 = o.groupby("donor_id", observed=True).agg(
    GSM=("gsm", "first"), Group=("group", "first"), Barcodes=("keep", "size"),
    UMI_low=("umi_lo", "first"), UMI_low_MAD_unfloored=("umi_lo_mad_unfloored", "first"), UMI_high=("umi_hi", "first"),
    Genes_low=("genes_lo", "first"), Genes_low_MAD_unfloored=("genes_lo_mad_unfloored", "first"), Genes_high=("genes_hi", "first"),
    Fail_UMI_low=("fail_umi_low", "sum"), Fail_UMI_high=("fail_umi_high", "sum"),
    Fail_genes_low=("fail_genes_low", "sum"), Fail_genes_high=("fail_genes_high", "sum"),
    Fail_mito_gt15=("fail_mito", "sum"), Pass_QC=("pass_qc", "sum"),
    Doublet_cutoff=("doublet_cutoff_used", "first"), Doublets_all=("predicted_doublet", "sum"),
    Kept=("keep", "sum")).reindex(ORDER)
e2["Doublets_among_pass_QC"] = o[o["pass_qc"]].groupby("donor_id", observed=True)["predicted_doublet"].sum().reindex(ORDER).fillna(0).astype(int)
e2["Doublet_rate_pct"] = (e2["Doublets_all"] / e2["Barcodes"] * 100).round(2)
e2["Retained_pct"] = (e2["Kept"] / e2["Barcodes"] * 100).round(2)
for c in ["UMI_low", "UMI_low_MAD_unfloored", "UMI_high", "Genes_low", "Genes_low_MAD_unfloored", "Genes_high"]:
    e2[c] = e2[c].round(1)
e2 = e2.reset_index().rename(columns={"donor_id": "Sample"})
e2.to_csv("results/tables/E2_filtering_summary.csv", index=False)
print("\n--- E2 FILTERING SUMMARY (regenerated from the h5ad obs; source: script 29 run) ---")
print(e2.to_string(index=False))
print(f"TOTAL: barcodes={e2.Barcodes.sum()} pass_QC={e2.Pass_QC.sum()} doublets_all={e2.Doublets_all.sum()} "
      f"doublets_among_pass_QC={e2.Doublets_among_pass_QC.sum()} kept={e2.Kept.sum()}")

# ---------------------------------------------------------------------------
# 2a. LV_UF-3: marker make-up of removed vs kept barcodes (counts > 0)
# ---------------------------------------------------------------------------
NEUT = ["FCGR3B", "CSF3R", "CXCR2", "MNDA", "G0S2"]
guard_check(NEUT, stage="30b_neut_markers")
m3 = (full.obs["donor_id"] == "LV_UF-3").values
sub3 = full[m3]
pos = pd.DataFrame((sub3[:, NEUT].X > 0).toarray(), columns=NEUT, index=sub3.obs_names)
pos["any_FCGR3B_CSF3R_CXCR2"] = pos[["FCGR3B", "CSF3R", "CXCR2"]].any(axis=1)
pos["status"] = np.where(sub3.obs["keep"], "kept",
                         np.where(sub3.obs["predicted_doublet"] & sub3.obs["pass_qc"], "removed: doublet only",
                                  np.where(sub3.obs["fail_genes_low"] | sub3.obs["fail_umi_low"], "removed: low UMI/genes",
                                           "removed: other QC (high/mito)")))
pos["n_genes"] = sub3.obs["n_genes"].values
pos["n_counts"] = sub3.obs["n_counts"].values
print("\n--- LV_UF-3: neutrophil-marker positivity (% of barcodes with counts>0) by QC fate ---")
print(pos.groupby("status").agg(N=("n_genes", "size"), median_genes=("n_genes", "median"), median_UMI=("n_counts", "median"),
                                **{f"pct_{g}": (g, lambda x: round(x.mean() * 100, 1)) for g in NEUT + ["any_FCGR3B_CSF3R_CXCR2"]}).to_string())
del sub3

# ---------------------------------------------------------------------------
# 3. Kept cells: normalise, HVG, PCA, UMAP (unintegrated), Harmony, UMAP, Leiden
# ---------------------------------------------------------------------------
adata = full[full.obs["keep"].values].copy()
del full
sc.pp.filter_genes(adata, min_cells=3)
adata.obs["donor_id"] = pd.Categorical(adata.obs["donor_id"], categories=ORDER)
adata.obs["group"] = pd.Categorical(adata.obs["group"], categories=["SV", "LV_NOT_UF", "LV_UF"])
tick(f"kept cells: {adata.n_obs} x {adata.n_vars} genes (min_cells=3)")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="donor_id")
guard_check(adata.var_names[adata.var["highly_variable"]].tolist(), stage="30b_hvg")
sc.pp.pca(adata, n_comps=30, mask_var="highly_variable", random_state=0)
tick("PCA done")

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca", random_state=0)
sc.tl.umap(adata, random_state=0)
adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_unintegrated", flavor="igraph", n_iterations=2, directed=False, random_state=0)
tick("unintegrated UMAP + Leiden(0.5) done")

ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs[["donor_id"]], "donor_id", max_iter_harmony=20, random_state=0)
Z = np.asarray(ho.Z_corr)
adata.obsm["X_pca_harmony"] = Z.T if Z.shape[0] != adata.n_obs else Z
tick(f"Harmony done (harmonypy 0.0.10), embedding {adata.obsm['X_pca_harmony'].shape}")

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony", random_state=0)
sc.tl.umap(adata, random_state=0)
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()

eval_idx = np.random.RandomState(42).choice(adata.n_obs, min(5000, adata.n_obs), replace=False)
res_rows, best = [], (None, -np.inf)
for r in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
    k = f"leiden_h_{r:.1f}"
    sc.tl.leiden(adata, resolution=r, key_added=k, flavor="igraph", n_iterations=2, directed=False, random_state=0)
    s = silhouette_score(adata.obsm["X_pca_harmony"][eval_idx], adata.obs[k].values[eval_idx])
    res_rows.append({"resolution": r, "n_clusters": adata.obs[k].nunique(), "silhouette": round(float(s), 4)})
    if s > best[1] + 1e-12:
        best = (r, s)
print("\n--- Leiden resolution selection on Harmony embedding (frozen rule: max silhouette, ties -> lower) ---")
print(pd.DataFrame(res_rows).to_string(index=False))
RES = best[0]
adata.obs["leiden"] = adata.obs[f"leiden_h_{RES:.1f}"]
tick(f"chosen resolution {RES} -> {adata.obs['leiden'].nunique()} clusters")

sil = {"donor_unintegrated": silhouette_score(adata.obsm["X_pca"][eval_idx], adata.obs["donor_id"].values[eval_idx]),
       "donor_harmony": silhouette_score(adata.obsm["X_pca_harmony"][eval_idx], adata.obs["donor_id"].values[eval_idx])}
print(f"Donor silhouette (lower = better mixed): unintegrated PCA {sil['donor_unintegrated']:+.4f} -> Harmony {sil['donor_harmony']:+.4f}")

# ---------------------------------------------------------------------------
# 4. Annotation
# ---------------------------------------------------------------------------
STROMAL = "stromal / mesothelial-lineage (unresolved)"
MODULES = {
    "T cell": ["CD3D", "CD3E", "CD2", "TRAC"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    "B cell": ["MS4A1", "CD79A", "CD19"],
    "Plasma cell": ["JCHAIN", "MZB1", "SDC1"],
    "Monocyte / macrophage": ["C1QA", "C1QB", "CD163", "MARCO", "CD14"],
    "cDC": ["CD1C", "CLEC10A", "FCER1A", "CLEC9A"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA"],
    "Neutrophil": NEUT,
    "Mast cell": ["TPSAB1", "CPA3", "KIT"],
    STROMAL: ["WT1", "MSLN", "UPK3B", "CALB2", "KRT19", "DCN", "PDGFRA", "LUM", "COL1A2"],
    "Endothelial": ["VWF", "CDH5", "KDR"],
}
guard_check([g for v in MODULES.values() for g in v], stage="30b_modules")
missing = {k: [g for g in v if g not in adata.var_names] for k, v in MODULES.items()}
print("Module genes missing from object:", {k: v for k, v in missing.items() if v})
for k, v in MODULES.items():
    sc.tl.score_genes(adata, [g for g in v if g in adata.var_names], score_name=f"score::{k}", random_state=0)
score_cols = [f"score::{k}" for k in MODULES]
cl_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
cl_scores.columns = list(MODULES)
cl_n = adata.obs["leiden"].value_counts()

labels = {}
for cl, row in cl_scores.iterrows():
    srt = row.sort_values(ascending=False)
    if cl_n[cl] < 50:
        labels[cl] = "Small (<50 cells)"
    elif srt.iloc[0] < 0.10:
        labels[cl] = "Unassigned (no module >= 0.10)"
    elif srt.iloc[0] - srt.iloc[1] < 0.05:
        labels[cl] = f"Mixed: {srt.index[0]} / {srt.index[1]}"
    else:
        labels[cl] = srt.index[0]
adata.obs["cell_type"] = adata.obs["leiden"].map(labels).astype(str)
out = cl_scores.round(3)
out.insert(0, "N", cl_n.reindex(out.index).values)
out["label"] = out.index.map(labels)
print("\n--- Cluster mean module scores (scanpy score_genes) and assigned label "
      "(top module; Unassigned if top<0.10; Mixed if margin<0.05; Small if <50 cells) ---")
print(out.to_string())


def expr_table(genes, groupby="leiden"):
    g = [x for x in genes if x in adata.var_names]
    X = adata[:, g].X
    X = X.toarray() if sp.issparse(X) else X
    df = pd.DataFrame(X > 0, columns=g, index=adata.obs_names)
    pct = df.groupby(adata.obs[groupby].values).mean().mul(100).round(1)
    pct.index.name = groupby
    return pct


print("\n--- Mono/Mac vs cDC markers: % cells with counts>0 per cluster ---")
mdc = expr_table(["CD1C", "CLEC10A", "FCER1A", "CLEC9A", "C1QA", "C1QB", "CD163", "MARCO", "CD14", "LILRA4"])
mdc.insert(0, "label", mdc.index.map(labels))
print(mdc.to_string())
Xd = {g: (adata[:, g].X > 0).toarray().ravel() for g in ["CD1C", "CLEC10A", "FCER1A", "CLEC9A", "C1QA", "C1QB"] if g in adata.var_names}
cdc2 = Xd["CD1C"] & (Xd["CLEC10A"] | Xd["FCER1A"]) & ~Xd["C1QA"] & ~Xd["C1QB"]
cdc1 = Xd["CLEC9A"] & ~Xd["C1QA"] & ~Xd["C1QB"]
adata.obs["cDC2_like_cell"], adata.obs["cDC1_like_cell"] = cdc2, cdc1
print(f"Cell-level: cDC2-like (CD1C+ & (CLEC10A+|FCER1A+) & C1QA-/C1QB-) = {cdc2.sum()} cells; "
      f"cDC1-like (CLEC9A+ & C1QA-/C1QB-) = {cdc1.sum()} cells")
print("cDC2-like cells by cluster:", adata.obs.loc[cdc2, "leiden"].value_counts().head(8).to_dict())

print("\n--- Plasma cell markers: % cells with counts>0 per cluster ---")
pl = expr_table(["JCHAIN", "MZB1", "SDC1", "MS4A1", "CD79A"])
pl.insert(0, "label", pl.index.map(labels))
print(pl.to_string())
pc_cell = ((adata[:, "JCHAIN"].X > 0).toarray().ravel() & (adata[:, "MZB1"].X > 0).toarray().ravel())
print(f"Cell-level JCHAIN+ & MZB1+ = {pc_cell.sum()} cells; by cluster:",
      adata.obs.loc[pc_cell, "leiden"].value_counts().head(8).to_dict())
print("JCHAIN+MZB1+ cells by donor:", adata.obs.loc[pc_cell, "donor_id"].value_counts().to_dict())

print("\n--- Endothelial (VWF, CDH5, KDR; PECAM1 not used): % cells with counts>0 per cluster (top 6 by VWF) ---")
en = expr_table(["VWF", "CDH5", "KDR", "PECAM1"])
en.insert(0, "label", en.index.map(labels))
print(en.sort_values("VWF", ascending=False).head(6).to_string())

# ---------------------------------------------------------------------------
# 5. Neutrophil confirmation
# ---------------------------------------------------------------------------
print("\n--- Neutrophil confirmation: % cells with counts>0 and mean log-norm expression, per cluster ---")
nt = expr_table(NEUT)
mean_n = pd.DataFrame(adata[:, NEUT].X.toarray(), columns=NEUT).groupby(adata.obs["leiden"].values).mean().round(2)
mean_n.columns = [f"mean_{c}" for c in NEUT]
nt = nt.join(mean_n)
nt.insert(0, "label", nt.index.map(labels))
nt.insert(1, "N", cl_n.reindex(nt.index).values)
print(nt.to_string())

# ---------------------------------------------------------------------------
# 2b. LV_UF-3 gene counts by cell type
# ---------------------------------------------------------------------------
print("\n--- LV_UF-3 (kept cells): genes / UMI per cell by cell type ---")
d3 = adata.obs[adata.obs["donor_id"] == "LV_UF-3"]
print(d3.groupby("cell_type").agg(N=("n_genes", "size"),
                                  genes_median=("n_genes", "median"), genes_q25=("n_genes", lambda x: x.quantile(.25)),
                                  genes_q75=("n_genes", lambda x: x.quantile(.75)), genes_min=("n_genes", "min"),
                                  UMI_median=("n_counts", "median"), UMI_min=("n_counts", "min"),
                                  mito_median=("pct_counts_mt", "median")).sort_values("N", ascending=False).round(1).to_string())
print("Same cell types, all OTHER donors (median genes):",
      adata.obs[adata.obs["donor_id"] != "LV_UF-3"].groupby("cell_type")["n_genes"].median().round(0).to_dict())

# ---------------------------------------------------------------------------
# Stromal co-expression
# ---------------------------------------------------------------------------
st_mask = (adata.obs["cell_type"] == STROMAL).values
cnt = {g: (adata[st_mask, g].X > 0).toarray().ravel() for g in ["WT1", "MSLN", "DCN", "PDGFRA"]}
n_st = int(st_mask.sum())
meso_any, fib_any = cnt["WT1"] | cnt["MSLN"], cnt["DCN"] | cnt["PDGFRA"]
rows = [("WT1+", cnt["WT1"]), ("MSLN+", cnt["MSLN"]), ("DCN+", cnt["DCN"]), ("PDGFRA+", cnt["PDGFRA"]),
        ("(WT1+ or MSLN+)", meso_any), ("(DCN+ or PDGFRA+)", fib_any),
        ("(WT1+ or MSLN+) AND (DCN+ or PDGFRA+)", meso_any & fib_any),
        ("WT1+ AND DCN+", cnt["WT1"] & cnt["DCN"]), ("MSLN+ AND PDGFRA+", cnt["MSLN"] & cnt["PDGFRA"]),
        ("WT1+ MSLN+ DCN+ PDGFRA+ (all four)", cnt["WT1"] & cnt["MSLN"] & cnt["DCN"] & cnt["PDGFRA"]),
        ("(WT1+ or MSLN+) only, DCN-/PDGFRA-", meso_any & ~fib_any), ("(DCN+ or PDGFRA+) only, WT1-/MSLN-", fib_any & ~meso_any)]
print(f"\n--- '{STROMAL}': co-expression (counts>0) among its {n_st} cells ---")
print(pd.DataFrame([(k, int(v.sum()), round(v.mean() * 100, 2) if n_st else np.nan) for k, v in rows],
                   columns=["pattern", "cells", "pct_of_stromal"]).to_string(index=False))
print("Stromal cells by donor:", adata.obs.loc[st_mask, "donor_id"].value_counts().to_dict())

# ---------------------------------------------------------------------------
# 6. Donor fractions, flags, testability
# ---------------------------------------------------------------------------
ct = pd.crosstab(adata.obs["leiden"], adata.obs["donor_id"])
frac = ct.div(ct.sum(axis=1), axis=0)
rows = []
for cl in ct.index:
    f = frac.loc[cl].sort_values(ascending=False)
    gfr = adata.obs.loc[adata.obs["leiden"] == cl, "group"].value_counts(normalize=True)
    rows.append({"cluster": cl, "label": labels[cl], "N": int(ct.loc[cl].sum()),
                 "largest_donor": f.index[0], "largest_frac": round(f.iloc[0], 3),
                 "n_donors_ge5pct": int((f >= 0.05).sum()),
                 "donors_ge5pct": ", ".join(f"{d}:{v*100:.0f}%" for d, v in f[f >= 0.05].items()),
                 "SV/LVnot/LVuf_%": "/".join(f"{gfr.get(g, 0)*100:.0f}" for g in ["SV", "LV_NOT_UF", "LV_UF"]),
                 "DONOR_DOMINATED(>=0.50)": bool(f.iloc[0] >= 0.50)})
dfr = pd.DataFrame(rows)
dfr.to_csv("results/tables/E3_cluster_donor_fractions_harmony.csv", index=False)
print("\n--- Per-cluster donor fractions after Harmony (flag: largest donor fraction >= 0.50) ---")
print(dfr.to_string(index=False))
ctu = pd.crosstab(adata.obs["leiden_unintegrated"], adata.obs["donor_id"])
mu = ctu.max(axis=1) / ctu.sum(axis=1)
print(f"Before Harmony (Leiden 0.5 on PCA): {len(mu)} clusters, {(mu >= 0.5).sum()} with largest donor >= 0.50; "
      f"median largest-donor fraction {mu.median():.3f}")
print(f"After Harmony  (Leiden {RES} on Harmony): {len(dfr)} clusters, {dfr['DONOR_DOMINATED(>=0.50)'].sum()} flagged; "
      f"median largest-donor fraction {dfr['largest_frac'].median():.3f}")

flagged = set(dfr.loc[dfr["DONOR_DOMINATED(>=0.50)"], "cluster"])
CONTRASTS = [("LV_UF", "LV_NOT_UF"), ("LV_UF", "SV"), ("LV_NOT_UF", "SV")]


def testability(obs, tag):
    obs = obs[~obs["cell_type"].str.startswith(("Small", "Unassigned", "Mixed"))]
    tab = pd.crosstab(obs["donor_id"], obs["cell_type"]).reindex(ORDER).fillna(0).astype(int)
    out = []
    for ctp in tab.columns:
        r = {"cell_type": ctp}
        for t, c in CONTRASTS:
            td = [d for d in ORDER if GROUP_OF[d] == t]; cd = [d for d in ORDER if GROUP_OF[d] == c]
            nt_, nc_ = int((tab.loc[td, ctp] >= 20).sum()), int((tab.loc[cd, ctp] >= 20).sum())
            r[f"{t}_vs_{c}"] = f"{'PASS' if (nt_ >= 3 and nc_ >= 2) else 'fail'} ({nt_}/{len(td)} vs {nc_}/{len(cd)})"
        out.append(r)
    res = pd.DataFrame(out)
    print(f"\n--- Donor x cell type counts [{tag}] ---")
    t2 = tab.copy(); t2.insert(0, "group", [GROUP_OF[d] for d in t2.index])
    print(t2.to_string())
    print(f"\n--- Testability [{tag}] (PASS = >=20 cells in >=3 target donors AND >=2 comparator donors; "
          f"shown as target donors>=20 / n vs comparator donors>=20 / n) ---")
    print(res.to_string(index=False))
    return res


ta = testability(adata.obs, "all clusters")
tb = testability(adata.obs[~adata.obs["leiden"].isin(flagged)], f"after removing donor-dominated clusters {sorted(flagged, key=int)}")
ta.to_csv("results/tables/E4_testability_all_clusters.csv", index=False)
tb.to_csv("results/tables/E4_testability_after_removing_donor_dominated.csv", index=False)

# ---------------------------------------------------------------------------
# 7. Figures
# ---------------------------------------------------------------------------
donor_cols = dict(zip(ORDER, plt.get_cmap("tab20").colors[:16]))
group_cols = {"SV": "#10b981", "LV_NOT_UF": "#3b82f6", "LV_UF": "#ef4444"}
perm = np.random.RandomState(0).permutation(adata.n_obs)
fig, axes = plt.subplots(2, 2, figsize=(20, 18))
for i, (emb, nm) in enumerate([("X_umap_unintegrated", "Before Harmony (PCA)"), ("X_umap_harmony", "After Harmony (by donor)")]):
    E = adata.obsm[emb][perm]
    for j, (key, cols) in enumerate([("donor_id", donor_cols), ("group", group_cols)]):
        ax = axes[i, j]
        v = adata.obs[key].astype(str).values[perm]
        ax.scatter(E[:, 0], E[:, 1], c=[cols[x] for x in v], s=0.4, rasterized=True, linewidths=0)
        ax.set_title(f"{nm} - coloured by {key}", fontsize=14, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(handles=[plt.Line2D([], [], ls="", marker="o", ms=7, color=c, label=k) for k, c in cols.items()],
                  fontsize=8 if key == "donor_id" else 11, loc="upper right", ncol=2 if key == "donor_id" else 1, markerscale=1)
plt.tight_layout()
plt.savefig("results/figures/E3_umap_donor_group_before_after_harmony.png", dpi=130)
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(22, 10))
E = adata.obsm["X_umap_harmony"]
sc.pl.embedding(adata, basis="X_umap_harmony", color="leiden", legend_loc="on data", ax=axes[0], show=False,
                title=f"Harmony UMAP - Leiden {RES} (donor-dominated: {sorted(flagged, key=int)})", size=2)
sc.pl.embedding(adata, basis="X_umap_harmony", color="cell_type", ax=axes[1], show=False, title="Harmony UMAP - cell type", size=2)
plt.tight_layout()
plt.savefig("results/figures/E3_umap_harmony_clusters_celltypes.png", dpi=130, bbox_inches="tight")
plt.close()

DOT = {"T": MODULES["T cell"][:3], "NK": MODULES["NK cell"][:3], "B": MODULES["B cell"], "Plasma": MODULES["Plasma cell"],
       "Mono/Mac": ["CD14", "C1QA", "C1QB", "CD163", "MARCO"], "cDC": MODULES["cDC"], "pDC": MODULES["pDC"],
       "Neut": NEUT, "Mast": MODULES["Mast cell"], "Strom/Meso": ["WT1", "MSLN", "UPK3B", "KRT19", "DCN", "PDGFRA", "LUM", "COL1A2"],
       "Endo": MODULES["Endothelial"]}
DOT = {k: [g for g in v if g in adata.var_names] for k, v in DOT.items()}
guard_check([g for v in DOT.values() for g in v], stage="30b_dotplot")
for scale, fn, ttl in [(None, "E3_marker_dotplot_unscaled.png", "mean log-normalised expression (no per-gene scaling)"),
                       ("var", "E3_marker_dotplot_scaled.png", "per-gene min-max scaled mean expression (0-1)")]:
    dp = sc.pl.dotplot(adata, DOT, groupby="cell_type", standard_scale=scale, return_fig=True, var_group_rotation=0,
                       colorbar_title="mean log1p(CP10k)" if scale is None else "scaled mean\n(per gene, 0-1)")
    dp.make_figure()
    fig = dp.fig
    fig.suptitle(f"Canonical lineage markers by annotated cell type - {ttl}", fontsize=14, fontweight="bold", y=1.06)
    fig.savefig(f"results/figures/{fn}", dpi=150, bbox_inches="tight")
    plt.close(fig)

keep_obs = ["gsm", "donor_id", "group", "n_counts", "n_genes", "pct_counts_mt", "doublet_score", "leiden",
            "leiden_unintegrated", "cell_type", "cDC2_like_cell", "cDC1_like_cell"] + score_cols
lite = ad.AnnData(X=sp.csr_matrix((adata.n_obs, 0), dtype=np.float32), obs=adata.obs[keep_obs].copy())
for k in ["X_pca", "X_pca_harmony", "X_umap_unintegrated", "X_umap_harmony"]:
    lite.obsm[k] = adata.obsm[k]
lite.uns["doublet_rule"] = adata.uns["doublet_rule"]
lite.uns["leiden_resolution"] = RES
lite.uns["harmony"] = "harmonypy 0.0.10 run_harmony(X_pca 30 PCs, vars_use=donor_id, max_iter_harmony=20, random_state=0)"
lite.obs.columns = [c.replace("::", "__").replace("/", "_").replace(" ", "_") for c in lite.obs.columns]
lite.write_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad")
tick("Saved figures + data/processed/GSE248762_harmony_annotated_obs.h5ad (obs/embeddings only). Script 30b complete.")
