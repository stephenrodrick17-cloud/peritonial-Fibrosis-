"""
Script 47: Pre-Registration, Data Preparation, and Analysis of QC Sensitivities (i) and (j)
Rules:
- Register provenance/stage4_addendum2.json BEFORE running
- Print filesystem mtime of the addendum
- (i) No upper UMI/gene ceiling (low floors, mito 15%, remove Scrublet & scDblFinder doublets)
- (j) One pooled fixed ceiling (99.5th percentile of pooled kept-cell UMI and genes)
- Print per-donor stromal cells, UMI per cell, genes per cell
- Run edgeR quasi-likelihood for LV_UF vs LV_NOT_UF and LV_UF vs SV
- Print side-by-side comparison of log2FC: Primary vs (g) vs (i) vs (j) with percent change
- Apply pre-set rule: if |log2FC| shrinks by >50% or flips sign, label QC-SENSITIVE
- Seed 42, no causal language.
"""
import os
import json
import hashlib
import datetime
import subprocess
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

np.random.seed(42)

print("=" * 80)
print("STAGE 4: PRE-SPECIFIED (SELF-DOCUMENTED) & EXECUTION OF QC SENSITIVITIES (i) AND (j)")
print("=" * 80)

# Helper function to get file provenance
def get_file_provenance(fpath):
    st = os.stat(fpath)
    mtime_utc = datetime.datetime.fromtimestamp(st.st_mtime, tz=datetime.timezone.utc).isoformat()
    with open(fpath, "rb") as f:
        sha256 = hashlib.sha256(f.read()).hexdigest()
    return {
        "Path": fpath,
        "Size_Bytes": st.st_size,
        "Mtime_UTC": mtime_utc,
        "SHA256": sha256
    }

# ---------------------------------------------------------------------------
# 1. Pre-Specify in provenance/stage4_addendum2.json
# ---------------------------------------------------------------------------
print("\n>>> 1. PRE-SPECIFYING (SELF-DOCUMENTED) SENSITIVITIES (i) & (j) IN provenance/stage4_addendum2.json...")
addendum2_path = "provenance/stage4_addendum2.json"

CUTOFF_SCRUBLET = 0.4390
FLOOR_UMI = 500.0
FLOOR_GENES = 200.0
MITO_CEILING = 15.0

addendum2 = {
    "timestamp_addendum2_frozen_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "dataset": "GSE248762 (PD effluent scRNA-seq, 16 donors)",
    "status": "pre-specified (self-documented) BEFORE running sensitivity analyses (i) and (j)",
    "sensitivity_i": {
        "title": "No upper UMI/gene ceilings",
        "description": (
            "Re-run stromal pseudobulk without upper adaptive UMI or gene ceilings. "
            "Keep low floors (UMI >= 500, genes >= 200, mito <= 15%). "
            "Filter doublets using both Scrublet (<= 0.4390) and scDblFinder (singlets only). "
            "Includes primary stromal cells and ceiling-excluded stromal/mesothelial lineage cells."
        ),
        "floor_umi": FLOOR_UMI,
        "floor_genes": FLOOR_GENES,
        "mito_ceiling_pct": MITO_CEILING,
        "scrublet_cutoff": CUTOFF_SCRUBLET,
        "contrasts": ["LV_UF_vs_LV_NOT_UF", "LV_UF_vs_SV"],
        "method": "edgeR quasi-likelihood across all 37,487 genes (filterByExpr + TMM over all genes)"
    },
    "sensitivity_j": {
        "title": "One pooled fixed upper ceiling across all donors",
        "description": (
            "Re-run stromal pseudobulk applying a single fixed upper ceiling defined by the "
            "99.5th percentile of pooled kept-cell UMI and 99.5th percentile of pooled kept-cell genes. "
            "Evaluates whether donor-specific adaptive 3-MAD ceilings introduced compositional bias."
        ),
        "percentile_ceiling": 99.5,
        "contrasts": ["LV_UF_vs_LV_NOT_UF", "LV_UF_vs_SV"],
        "method": "edgeR quasi-likelihood across all 37,487 genes (filterByExpr + TMM over all genes)"
    },
    "pre_set_interpretation_rule": (
        "If for any hub gene |log2FC| shrinks by more than 50% (|log2FC_sens| < 0.5 * |log2FC_primary|) "
        "or flips sign in sensitivity (g), (i), or (j), label it QC-SENSITIVE and exclude it from "
        "any 'supported' claim."
    )
}

os.makedirs("provenance", exist_ok=True)
with open(addendum2_path, "w", encoding="utf-8") as f:
    json.dump(addendum2, f, indent=2)

prov_add2 = get_file_provenance(addendum2_path)
print(f"Pre-specification (self-documented) addendum successfully written to: {prov_add2['Path']}")
print(f"Filesystem mtime (UTC): {prov_add2['Mtime_UTC']}")
print(f"SHA256:                 {prov_add2['SHA256']}")

# ---------------------------------------------------------------------------
# 2. Prepare Population for Sensitivity (i) and (j)
# ---------------------------------------------------------------------------
print("\n>>> 2. PREPARING POPULATIONS FOR SENSITIVITY (i) AND (j)...")

# Load full-transcriptome QC object (37,476 genes) and hub object (11 genes)
adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad")
adata_hubs = ad.read_h5ad("sealed/hub_counts.h5ad")

# Load scDblFinder calls
df_sc = pd.read_csv("results/tables/scDblFinder_per_barcode_calls.csv", index_col="barcode")
adata_qc.obs["scDblFinder_class"] = df_sc.loc[adata_qc.obs.index, "scDblFinder_class"]

# Load primary annotated obs
adata_annot = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad")
primary_obs = adata_annot.obs.copy()
primary_obs.index = primary_obs["_index"].values if "_index" in primary_obs.columns else primary_obs.index.values

# Compute 99.5th percentiles of kept cells
p995_umi = float(np.percentile(primary_obs["n_counts"], 99.5))
p995_genes = float(np.percentile(primary_obs["n_genes"], 99.5))
print(f"Pooled Kept-Cell 99.5th Percentiles: UMI = {p995_umi:.1f}, Genes = {p995_genes:.1f}")

# Primary stromal cells passing doublet filter
st_primary_barcodes = primary_obs[
    (primary_obs["cell_type"] == "stromal / mesothelial-lineage (unresolved)") &
    (~primary_obs["leiden"].isin(["4", "12"]))
].index.values

st_primary_sub = adata_qc.obs.loc[st_primary_barcodes]
st_prim_pass_mask = (st_primary_sub["doublet_score"] <= CUTOFF_SCRUBLET) & (st_primary_sub["scDblFinder_class"] != "doublet")
st_prim_pass_barcodes = st_primary_sub[st_prim_pass_mask].index.values

# Ceiling-excluded stromal cells
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

ceiling_st_obs = ceiling_adata.obs[(top_module == "Stromal_mesothelial").values].copy()
ceiling_st_pass = ceiling_st_obs[
    (ceiling_st_obs["doublet_score"] <= CUTOFF_SCRUBLET) &
    (ceiling_st_obs["scDblFinder_class"] != "doublet") &
    (ceiling_st_obs["pct_counts_mt"] <= MITO_CEILING)
]
ceiling_rescued_barcodes = ceiling_st_pass.index.values

# Population (i): No upper ceilings
barcodes_i = np.concatenate([st_prim_pass_barcodes, ceiling_rescued_barcodes])
obs_i = adata_qc.obs.loc[barcodes_i]

# Population (j): Apply pooled 99.5th percentile fixed ceiling
pass_j_mask = (obs_i["n_counts"] <= p995_umi) & (obs_i["n_genes"] <= p995_genes)
barcodes_j = obs_i[pass_j_mask].index.values
obs_j = obs_i.loc[barcodes_j]

print(f"Total Stromal Cells in Sensitivity (i) [No Upper Ceiling]: {len(barcodes_i)}")
print(f"Total Stromal Cells in Sensitivity (j) [Pooled 99.5% Ceiling]: {len(barcodes_j)}")

# Helper to aggregate pseudobulk
all_var_names = list(adata_qc.var_names) + list(adata_hubs.var_names)

def create_pseudobulk(barcodes, label, prefix):
    sub_qc = adata_qc[barcodes, :].copy()
    sub_hubs = adata_hubs[barcodes, :].copy()
    mat = sparse.hstack([sub_qc.X, sub_hubs.X]).tocsr()
    
    adata_sub = ad.AnnData(X=mat, obs=sub_qc.obs.copy())
    adata_sub.var_names = all_var_names
    
    donors = sorted(adata_sub.obs["donor_id"].unique())
    donor_map = {d: i for i, d in enumerate(donors)}
    d_idx = adata_sub.obs["donor_id"].map(donor_map).values
    n_don = len(donors)
    n_c = adata_sub.n_obs
    
    S = sparse.csr_matrix((np.ones(n_c), (d_idx, np.arange(n_c))), shape=(n_don, n_c))
    pb_mat = S @ adata_sub.X
    
    meta = adata_sub.obs.groupby("donor_id").agg(
        group=("group", "first"),
        n_cells=("group", "size"),
        mean_umi=("n_counts", "mean"),
        mean_genes=("n_genes", "mean")
    ).loc[donors]
    
    df_counts = pd.DataFrame(pb_mat.toarray(), index=donors, columns=adata_sub.var_names)
    
    pb_dir = "data/processed/pseudobulk"
    cpath = f"{pb_dir}/pb_sens_{prefix}_counts.csv"
    mpath = f"{pb_dir}/pb_sens_{prefix}_metadata.csv"
    df_counts.to_csv(cpath)
    meta.to_csv(mpath)
    
    print(f"\n[{label}]:")
    print(f"  Saved {cpath} and {mpath}")
    print(f"  Per-donor library sizes:")
    print(meta.to_string())
    return cpath, mpath, meta

c_i, m_i, meta_i = create_pseudobulk(barcodes_i, "Sensitivity (i) No Upper Ceilings", "i_no_upper_ceilings")
c_j, m_j, meta_j = create_pseudobulk(barcodes_j, "Sensitivity (j) Pooled Fixed Ceiling", "j_pooled_fixed_ceiling")

# ---------------------------------------------------------------------------
# 3. Execute edgeR on (i) and (j) via Rscript
# ---------------------------------------------------------------------------
print("\n>>> 3. EXECUTING edgeR QUASI-LIKELIHOOD FOR (i) AND (j)...")

r_code = f"""
suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
  library(data.table)
}})

set.seed(42)

pb_dir <- "data/processed/pseudobulk"
hub_genes <- c("COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
               "EDIL3", "LOX", "INHBA", "ISM1", "COMP")

run_de <- function(cfile, mfile, ana_name) {{
  dt <- data.table::fread(cfile, data.table = FALSE)
  rownames(dt) <- dt[[1]]
  counts <- dt[, -1, drop = FALSE]
  meta <- read.csv(mfile, row.names = 1, stringsAsFactors = FALSE)
  
  common_donors <- intersect(rownames(meta), rownames(counts))
  meta <- meta[common_donors, ]
  counts <- counts[common_donors, ]
  
  mat <- t(as.matrix(counts))
  groups <- factor(meta$group, levels = c("LV_NOT_UF", "LV_UF", "SV"))
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  dge <- DGEList(counts = mat, group = groups)
  keep <- filterByExpr(dge, design)
  hub_in_mat <- intersect(hub_genes, rownames(mat))
  keep[hub_in_mat] <- TRUE
  
  dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
  dge_filt <- calcNormFactors(dge_filt, method = "TMM")
  dge_filt <- estimateDisp(dge_filt, design, robust = TRUE)
  fit <- glmQLFit(dge_filt, design, robust = TRUE)
  
  contrasts_list <- list(
    LV_UF_vs_LV_NOT_UF = makeContrasts(LV_UF - LV_NOT_UF, levels = design),
    LV_UF_vs_SV = makeContrasts(LV_UF - SV, levels = design)
  )
  
  res_list <- list()
  for (cname in names(contrasts_list)) {{
    qlf <- glmQLFTest(fit, contrast = contrasts_list[[cname]])
    tab <- topTags(qlf, n = Inf)$table
    df_denom <- qlf$df.total[1]
    t_crit <- qt(0.975, df = df_denom)
    
    n_uf <- sum(meta$group == "LV_UF")
    n_not_uf <- sum(meta$group == "LV_NOT_UF")
    n_sv <- sum(meta$group == "SV")
    
    for (g in hub_genes) {{
      if (g %in% rownames(tab)) {{
        row <- tab[g, ]
        lfc <- row$logFC
        pval <- row$PValue
        f_stat <- row$F
        se <- if (f_stat > 0) abs(lfc) / sqrt(f_stat) else NA
        ci_lo <- lfc - t_crit * se
        ci_hi <- lfc + t_crit * se
      }} else {{
        lfc <- NA; pval <- NA; se <- NA; ci_lo <- NA; ci_hi <- NA
      }}
      res_list[[paste(cname, g, sep = "__")]] <- data.frame(
        Analysis = ana_name,
        Cell_Type = "stromal / mesothelial-lineage (unresolved)",
        Contrast = cname,
        Gene = g,
        log2FC = round(lfc, 4),
        SE = round(se, 4),
        CI_95_low = round(ci_lo, 4),
        CI_95_high = round(ci_hi, 4),
        PValue = pval,
        QL_Denom_df = round(df_denom, 2),
        Donors_LV_UF = n_uf,
        Donors_LV_NOT_UF = n_not_uf,
        Donors_SV = n_sv,
        Cells_total = sum(meta$n_cells),
        stringsAsFactors = FALSE
      )
    }}
  }}
  res_df <- do.call(rbind, res_list)
  res_df$BH_FDR <- p.adjust(res_df$PValue, method = "BH")
  res_df
}}

res_i <- run_de("{c_i}", "{m_i}", "sens_i_no_upper_ceilings")
res_j <- run_de("{c_j}", "{m_j}", "sens_j_pooled_fixed_ceiling")

res_ij <- rbind(res_i, res_j)
write.csv(res_ij, "results/tables/stage4_sensitivity_i_j_edger_pseudobulk.csv", row.names = FALSE)
cat("Sensitivity (i) and (j) edgeR differential testing complete.\\n")
"""

r_script_path = "scripts/47b_run_edger_sens_ij.R"
with open(r_script_path, "w", encoding="utf-8") as f:
    f.write(r_code)

res_cmd = subprocess.run(
    ["C:\\Program Files\\R\\R-4.4.2\\bin\\Rscript.exe", r_script_path],
    capture_output=True,
    text=True,
    check=False
)
print("Rscript Output:")
print(res_cmd.stdout)
if res_cmd.stderr:
    print("Rscript Messages/Warnings:")
    print(res_cmd.stderr)

# ---------------------------------------------------------------------------
# 4. Read Back Results & Compare Primary vs (g) vs (i) vs (j)
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print(">>> 4. READBACK AUDIT: PRIMARY vs (g) vs (i) vs (j) LOG2FC COMPARISON")
print("=" * 80)

ij_csv = "results/tables/stage4_sensitivity_i_j_edger_pseudobulk.csv"
prov_ij = get_file_provenance(ij_csv)
print(f"Target CSV:   {prov_ij['Path']}")
print(f"Mtime (UTC):  {prov_ij['Mtime_UTC']}")
print(f"SHA256:       {prov_ij['SHA256']}")
print(f"Size (Bytes): {prov_ij['Size_Bytes']}")

df_ij = pd.read_csv(ij_csv)

# Load Primary & (g)
df_prim = pd.read_csv("results/tables/stage4_primary_edger_pseudobulk.csv")
prim_st = df_prim[
    (df_prim["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
    (df_prim["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].set_index("Gene")

df_gh = pd.read_csv("results/tables/stage4_sensitivity_g_h_edger_pseudobulk.csv")
g_st = df_gh[
    (df_gh["Analysis"] == "sens_g_ceiling_rescued") &
    (df_gh["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].set_index("Gene")

i_st = df_ij[
    (df_ij["Analysis"] == "sens_i_no_upper_ceilings") &
    (df_ij["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].set_index("Gene")

j_st = df_ij[
    (df_ij["Analysis"] == "sens_j_pooled_fixed_ceiling") &
    (df_ij["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].set_index("Gene")

eval_genes = ["THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", "EDIL3", "LOX", "INHBA", "COMP"]

comp_rows = []
for g in eval_genes:
    p_lfc = prim_st.loc[g, "log2FC"]
    g_lfc = g_st.loc[g, "log2FC"]
    i_lfc = i_st.loc[g, "log2FC"]
    j_lfc = j_st.loc[g, "log2FC"]
    
    pct_chg_g = ((g_lfc - p_lfc) / abs(p_lfc)) * 100.0
    pct_chg_i = ((i_lfc - p_lfc) / abs(p_lfc)) * 100.0
    pct_chg_j = ((j_lfc - p_lfc) / abs(p_lfc)) * 100.0
    
    # Pre-set rule: |log2FC| shrinks by >50% or flips sign
    # Shrunk by > 50%: |sens_lfc| < 0.5 * |p_lfc|
    # Sign flip: (sens_lfc * p_lfc) < 0
    def is_qc_sensitive(sens_lfc, p_lfc):
        if (sens_lfc * p_lfc) < 0:
            return True # sign flip
        if abs(sens_lfc) < 0.5 * abs(p_lfc):
            return True # shrunk by > 50%
        return False
    
    flag_g = is_qc_sensitive(g_lfc, p_lfc)
    flag_i = is_qc_sensitive(i_lfc, p_lfc)
    flag_j = is_qc_sensitive(j_lfc, p_lfc)
    
    any_qc_sens = flag_g or flag_i or flag_j
    qc_label = "QC-SENSITIVE" if any_qc_sens else "ROBUST_TO_QC"
    
    comp_rows.append({
        "Gene": g,
        "Primary_log2FC": p_lfc,
        "Sens_g_log2FC": g_lfc,
        "Pct_Change_g": round(pct_chg_g, 1),
        "Sens_i_log2FC": i_lfc,
        "Pct_Change_i": round(pct_chg_i, 1),
        "Sens_j_log2FC": j_lfc,
        "Pct_Change_j": round(pct_chg_j, 1),
        "Classification": qc_label
    })

df_qc_comp = pd.DataFrame(comp_rows)
out_qc_csv = "results/tables/stage4_qc_ceiling_sensitivity_comparison.csv"
df_qc_comp.to_csv(out_qc_csv, index=False)
prov_qc_comp = get_file_provenance(out_qc_csv)

print(f"\nSaved QC comparison to: {prov_qc_comp['Path']}")
print(f"Mtime (UTC):            {prov_qc_comp['Mtime_UTC']}")
print(f"SHA256:                 {prov_qc_comp['SHA256']}")

print("\n--- Side-by-Side Comparison of Primary vs (g) vs (i) vs (j) [LV_UF vs LV_NOT_UF] ---")
print(df_qc_comp.to_string(index=False))

print("\nSummary Classification under Pre-Set Rule:")
print(df_qc_comp["Classification"].value_counts().to_string())

print("\nQC Sensitivity Analysis (i) and (j) completed successfully.")
print("=" * 80)
