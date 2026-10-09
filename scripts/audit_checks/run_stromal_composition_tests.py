# scripts/audit_checks/run_stromal_composition_tests.py
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
import subprocess
import os

print("=== SCRIPT: STROMAL COMPOSITION EDGE R PSEUDOBULK SUITE ===")

# 1. Load annotated obs to get stromal cell barcodes
obs_annot = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad").obs
st_mask = obs_annot["cell_type"] == "stromal / mesothelial-lineage (unresolved)"
obs_st = obs_annot[st_mask].copy()
st_barcodes = obs_st["_index"].values if "_index" in obs_st.columns else obs_st.index.values

# Exclude donor-dominated clusters 4 and 12 if any stromal cells are in them (check)
# In script 37a: adata_primary = adata_full[~adata_full.obs["leiden"].isin(["4", "12"])]
st_barcodes_clean = [b for b in st_barcodes if obs_st.loc[obs_st["_index"] == b if "_index" in obs_st.columns else b, "leiden"].values[0] not in ["4", "12"]]
print(f"Total stromal barcodes: {len(st_barcodes)}, Clean (excl C4/C12): {len(st_barcodes_clean)}")

# 2. Load QC matrix and Hub counts
adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad")
adata_hubs = ad.read_h5ad("sealed/hub_counts.h5ad")

# Align
st_qc = adata_qc[st_barcodes_clean, :].copy()
st_hubs = adata_hubs[st_barcodes_clean, :].copy()

# Add metadata
st_qc.obs["donor_id"] = [obs_st.loc[obs_st["_index"] == b if "_index" in obs_st.columns else b, "donor_id"].values[0] for b in st_barcodes_clean]
st_qc.obs["group"] = [obs_st.loc[obs_st["_index"] == b if "_index" in obs_st.columns else b, "group"].values[0] for b in st_barcodes_clean]
st_qc.obs["n_counts"] = st_qc.X.sum(axis=1).A1 if sparse.issparse(st_qc.X) else st_qc.X.sum(axis=1)

# Compute mesothelial and fibroblast scores
meso_genes = ["WT1", "MSLN", "UPK3B", "KRT19"]
fibro_genes = ["DCN", "PDGFRA", "LUM", "COL1A2"]

sc_norm = st_qc.copy()
sc.pp.normalize_total(sc_norm, target_sum=1e4)
sc.pp.log1p(sc_norm)
sc.tl.score_genes(sc_norm, meso_genes, score_name="meso_score")
sc.tl.score_genes(sc_norm, fibro_genes, score_name="fibro_score")

st_qc.obs["meso_score"] = sc_norm.obs["meso_score"].values
st_qc.obs["fibro_score"] = sc_norm.obs["fibro_score"].values
st_qc.obs["sublineage"] = np.where(st_qc.obs["meso_score"] > st_qc.obs["fibro_score"], "Mesothelial", "Fibroblast")

# Combine transcriptome + hub genes: 37,476 + 11 = 37,487
combined_X = sparse.hstack([st_qc.X, st_hubs.X]).tocsr()
combined_var_names = list(st_qc.var_names) + list(st_hubs.var_names)
st_full = ad.AnnData(X=combined_X, obs=st_qc.obs.copy(), var=pd.DataFrame(index=combined_var_names))

# Compute per-donor mesothelial fraction
donor_meso_frac = st_full.obs.groupby("donor_id").apply(lambda x: (x["sublineage"] == "Mesothelial").mean()).to_dict()

hub_genes = ["COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", "EDIL3", "LOX", "INHBA", "ISM1", "COMP"]

# Function to aggregate pseudobulk
def aggregate_pb(adata_sub):
    donors = sorted(adata_sub.obs["donor_id"].unique())
    donor_map = {d: i for i, d in enumerate(donors)}
    donor_idx = adata_sub.obs["donor_id"].map(donor_map).values
    n_donors = len(donors)
    n_cells = adata_sub.n_obs
    S = sparse.csr_matrix((np.ones(n_cells), (donor_idx, np.arange(n_cells))), shape=(n_donors, n_cells))
    pb = S @ adata_sub.X
    meta = adata_sub.obs.groupby("donor_id").agg(
        group=("group", "first"),
        n_cells=("group", "size")
    ).loc[donors]
    meta["meso_frac"] = [donor_meso_frac[d] for d in donors]
    df_counts = pd.DataFrame(pb.toarray(), index=donors, columns=adata_sub.var_names)
    return df_counts, meta

# 1. Full stromal
df_counts_all, meta_all = aggregate_pb(st_full)
df_counts_all.to_csv("data/processed/pseudobulk/pb_audit_stromal_all_counts.csv")
meta_all.to_csv("data/processed/pseudobulk/pb_audit_stromal_all_meta.csv")

# 2. Mesothelial subset with donors >= 20 cells
meso_cells = st_full[st_full.obs["sublineage"] == "Mesothelial"].copy()
meso_counts_per_donor = meso_cells.obs["donor_id"].value_counts()
meso_keep_donors = meso_counts_per_donor[meso_counts_per_donor >= 20].index.tolist()
meso_cells_filt = meso_cells[meso_cells.obs["donor_id"].isin(meso_keep_donors)].copy()
df_counts_meso, meta_meso = aggregate_pb(meso_cells_filt)
df_counts_meso.to_csv("data/processed/pseudobulk/pb_audit_stromal_meso_counts.csv")
meta_meso.to_csv("data/processed/pseudobulk/pb_audit_stromal_meso_meta.csv")

# 3. Fibroblast subset with donors >= 20 cells
fibro_cells = st_full[st_full.obs["sublineage"] == "Fibroblast"].copy()
fibro_counts_per_donor = fibro_cells.obs["donor_id"].value_counts()
fibro_keep_donors = fibro_counts_per_donor[fibro_counts_per_donor >= 20].index.tolist()
fibro_cells_filt = fibro_cells[fibro_cells.obs["donor_id"].isin(fibro_keep_donors)].copy()
df_counts_fibro, meta_fibro = aggregate_pb(fibro_cells_filt)
df_counts_fibro.to_csv("data/processed/pseudobulk/pb_audit_stromal_fibro_counts.csv")
meta_fibro.to_csv("data/processed/pseudobulk/pb_audit_stromal_fibro_meta.csv")

print(f"Full stromal: {len(meta_all)} donors")
print(f"Mesothelial subset (>=20 cells): {len(meta_meso)} donors ({dict(meta_meso['group'].value_counts())})")
print(f"Fibroblast subset (>=20 cells): {len(meta_fibro)} donors ({dict(meta_fibro['group'].value_counts())})")

# Write R script to run edgeR for the three analyses
r_script_content = """
suppressPackageStartupMessages({
  library(edgeR)
})

hub_genes <- c("COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
               "EDIL3", "LOX", "INHBA", "ISM1", "COMP")

run_edger <- function(counts_file, meta_file, add_covariate=FALSE) {
  dt <- data.table::fread(counts_file, data.table = FALSE)
  rownames(dt) <- dt[[1]]
  counts <- dt[, -1, drop = FALSE]
  meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
  
  common <- intersect(rownames(meta), rownames(counts))
  meta <- meta[common, ]
  counts <- counts[common, ]
  
  mat <- t(as.matrix(counts))
  groups <- factor(meta$group, levels = c("LV_NOT_UF", "LV_UF", "SV"))
  
  if (add_covariate) {
    cov_val <- meta$meso_frac
    design <- model.matrix(~ 0 + groups + cov_val)
    colnames(design) <- c("LV_NOT_UF", "LV_UF", "SV", "meso_frac")
  } else {
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- levels(groups)
  }
  
  dge <- DGEList(counts = mat, group = groups)
  keep <- filterByExpr(dge, design)
  hub_in_mat <- intersect(hub_genes, rownames(mat))
  keep[hub_in_mat] <- TRUE
  dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
  dge_filt <- calcNormFactors(dge_filt, method = "TMM")
  dge_filt <- estimateDisp(dge_filt, design, robust = TRUE)
  fit <- glmQLFit(dge_filt, design, robust = TRUE)
  
  contrast_uf <- makeContrasts(LV_UF - LV_NOT_UF, levels = design)
  qlf <- glmQLFTest(fit, contrast = contrast_uf)
  tab <- topTags(qlf, n = Inf)$table
  
  res <- list()
  for (g in hub_genes) {
    if (g %in% rownames(tab)) {
      res[[g]] <- data.frame(Gene = g, log2FC = tab[g, "logFC"], PValue = tab[g, "PValue"])
    } else {
      res[[g]] <- data.frame(Gene = g, log2FC = NA, PValue = NA)
    }
  }
  do.call(rbind, res)
}

cat("=== 1. Stromal with Mesothelial Fraction Covariate ===\\n")
res_cov <- run_edger("data/processed/pseudobulk/pb_audit_stromal_all_counts.csv", 
                     "data/processed/pseudobulk/pb_audit_stromal_all_meta.csv", add_covariate=TRUE)
print(res_cov)

cat("\\n=== 2. Mesothelial-Scored Cells Subset (>=20 cells/donor) ===\\n")
res_meso <- run_edger("data/processed/pseudobulk/pb_audit_stromal_meso_counts.csv", 
                      "data/processed/pseudobulk/pb_audit_stromal_meso_meta.csv", add_covariate=FALSE)
print(res_meso)

cat("\\n=== 3. Fibroblast-Scored Cells Subset (>=20 cells/donor) ===\\n")
res_fibro <- run_edger("data/processed/pseudobulk/pb_audit_stromal_fibro_counts.csv", 
                       "data/processed/pseudobulk/pb_audit_stromal_fibro_meta.csv", add_covariate=FALSE)
print(res_fibro)
"""

with open("scripts/audit_checks/run_stromal_composition_edger.R", "w") as f:
    f.write(r_script_content)

print("Saved scripts/audit_checks/run_stromal_composition_edger.R. Ready for Rscript execution.")
