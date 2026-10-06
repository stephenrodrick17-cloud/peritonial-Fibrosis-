# scripts/audit_checks/check_directive5_lineage_covariates.py
import gzip
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import anndata as ad
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

print("=== CHECK DIRECTIVE 5: LINEAGE SIGNATURE COVARIATES IN GSE125498 ===")

# 1. Load GPL10558 probe mappings
annot_file = "data/raw/GPL10558.annot.gz"
probe_gene_map = {}
with gzip.open(annot_file, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!platform_table_begin"):
            break
    header = f.readline().strip().split("\t")
    id_idx = header.index("ID")
    symbol_idx = header.index("Gene symbol")
    for line in f:
        if line.startswith("!platform_table_end"):
            break
        parts = line.strip().split("\t")
        if len(parts) > symbol_idx:
            pid = parts[id_idx]
            sym = parts[symbol_idx]
            if sym:
                for s in sym.split("///"):
                    s = s.strip()
                    if s:
                        probe_gene_map.setdefault(s, []).append(pid)

gpl_genes = set(probe_gene_map.keys())
print(f"Total unique genes mapped on GPL10558: {len(gpl_genes)}")

# 2. Check FCGR3A vs FCGR3B annotation in GPL10558
print("\n--- FCGR3A vs FCGR3B GPL10558 Probes & Sequence Specificity ---")
with gzip.open(annot_file, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!platform_table_begin"):
            break
    header = f.readline().strip().split("\t")
    id_idx = header.index("ID")
    symbol_idx = header.index("Gene symbol")
    title_idx = header.index("Gene title")
    seq_idx = header.index("Platform_SEQUENCE") if "Platform_SEQUENCE" in header else -1
    acc_idx = header.index("GenBank Accession") if "GenBank Accession" in header else -1
    for line in f:
        if line.startswith("!platform_table_end"):
            break
        parts = line.strip().split("\t")
        if len(parts) > symbol_idx:
            sym = parts[symbol_idx]
            if any(target in sym for target in ["FCGR3A", "FCGR3B"]):
                pid = parts[id_idx]
                title = parts[title_idx] if len(parts) > title_idx else ""
                acc = parts[acc_idx] if acc_idx != -1 and len(parts) > acc_idx else ""
                seq = parts[seq_idx] if seq_idx != -1 and len(parts) > seq_idx else ""
                print(f"Probe: {pid:<14} | Symbol: {sym:<12} | Title: {title:<35} | Acc: {acc:<12} | Seq: {seq[:30]}...")

# 3. Load GSE248762 to identify top 20 specific marker genes per lineage
print("\n--- Identifying Top 20 Lineage Specific Marker Genes from GSE248762 ---")
adata_annot = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad")
adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad")

# Subsample or take cell types of interest
lineages = ["T cell", "NK cell", "Monocyte / macrophage", "cDC", "Neutrophil"]

# Barcode alignment
kept_barcodes = adata_annot.obs["_index"].values if "_index" in adata_annot.obs.columns else adata_annot.obs.index.values
sub_qc = adata_qc[kept_barcodes, :].copy()
sub_qc.obs["cell_type"] = adata_annot.obs["cell_type"].values

# Filter genes to those present on GPL10558
valid_genes = [g for g in sub_qc.var_names if g in gpl_genes]
sub_qc = sub_qc[:, valid_genes].copy()
print(f"Matrix for marker scoring: {sub_qc.shape[0]} cells x {sub_qc.shape[1]} genes on GPL10558")

# Normalize and log1p
sc.pp.normalize_total(sub_qc, target_sum=1e4)
sc.pp.log1p(sub_qc)

# Filter out low-expressing genes
sc.pp.filter_genes(sub_qc, min_cells=50)

# Rank genes groups
sc.tl.rank_genes_groups(sub_qc, groupby="cell_type", method="wilcoxon", groups=lineages, reference="rest")

top_markers = {}
for ct in lineages:
    df_markers = sc.get.rank_genes_groups_df(sub_qc, group=ct)
    # top 20 genes with highest logfoldchanges and score > 0
    top20 = df_markers.head(20)["names"].tolist()
    top_markers[ct] = top20
    print(f"  {ct} Top 20: {', '.join(top20[:10])}...")

# 4. Load GSE125498 expression matrix
series_file = "data/raw/GSE125498_series_matrix.txt.gz"
with gzip.open(series_file, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
    header = f.readline().strip().split("\t")
    sample_ids = [s.strip('"') for s in header[1:]]
    expr_dict = {}
    for line in f:
        if line.startswith("!series_matrix_table_end"):
            break
        parts = line.strip().split("\t")
        pid = parts[0].strip('"')
        expr_dict[pid] = [float(x) for x in parts[1:]]

df_gse = pd.DataFrame(expr_dict, index=sample_ids).T # probes x samples

# Map each gene to its highest-mean probe in GSE125498
gene_to_probe = {}
all_needed_genes = set([g for sublist in top_markers.values() for g in sublist] + ["VCAN", "FN1", "COL8A1"])
for g in all_needed_genes:
    pids = [p for p in probe_gene_map.get(g, []) if p in df_gse.index]
    if pids:
        # choose probe with highest mean expression
        means = [df_gse.loc[p].mean() for p in pids]
        best_p = pids[np.argmax(means)]
        gene_to_probe[g] = best_p

# Load metadata
df_meta = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
meta_map = dict(zip(df_meta["sample_id"], df_meta["stage_binary"]))
y_group = np.array([meta_map[s] for s in df_gse.columns]) # 1 = LPD (13), 0 = SPD (20)

# 5. Score samples by mean z-score of top 20 genes
sig_scores = {}
for ct in lineages:
    ct_genes = [g for g in top_markers[ct] if g in gene_to_probe]
    probes = [gene_to_probe[g] for g in ct_genes]
    mat = df_gse.loc[probes].values # genes x samples
    # z-score per gene across samples
    z_mat = stats.zscore(mat, axis=1)
    mean_z = np.nanmean(z_mat, axis=0)
    sig_scores[ct] = mean_z

df_scores = pd.DataFrame(sig_scores, index=df_gse.columns)
df_scores["group"] = y_group

print("\n--- Correlation of Lineage Scores with Group (LPD vs SPD) ---")
for ct in lineages:
    score = df_scores[ct].values
    r_pear, p_pear = stats.pearsonr(score, y_group)
    r_spear, p_spear = stats.spearmanr(score, y_group)
    mean_lpd = score[y_group == 1].mean()
    mean_spd = score[y_group == 0].mean()
    tt = stats.ttest_ind(score[y_group == 1], score[y_group == 0])
    print(f"{ct:<25}: LPD mean = {mean_lpd:+.3f}, SPD mean = {mean_spd:+.3f} | Pearson r = {r_pear:+.3f} (P = {p_pear:.4f}) | Welch P = {tt.pvalue:.4f}")

# 6. Regress VCAN, FN1, COL8A1 adjusting for the 5 scores
print("\n--- Target Gene Adjustment Models (OLS / Limma equivalent) ---")
X_covs = df_scores[lineages].values
X_all = sm.add_constant(np.column_stack([y_group, X_covs]))

targets = [
    ("VCAN", gene_to_probe["VCAN"]),
    ("FN1 (MaxMean)", gene_to_probe["FN1"]),
    ("FN1 (Lowest-P)", "ILMN_1778237"),
    ("COL8A1 (MaxMean)", gene_to_probe["COL8A1"]),
    ("COL8A1 (Lowest-P)", "ILMN_2402392")
]

for label, pid in targets:
    y = df_gse.loc[pid].values
    # Unadjusted
    model_unadj = sm.OLS(y, sm.add_constant(y_group)).fit()
    lfc_unadj = model_unadj.params[1]
    p_unadj = model_unadj.pvalues[1]
    
    # Adjusted for all 5 lineage signatures
    model_adj = sm.OLS(y, X_all).fit()
    lfc_adj = model_adj.params[1]
    p_adj = model_adj.pvalues[1]
    t_adj = model_adj.tvalues[1]
    
    print(f"\nTarget: {label} ({pid})")
    print(f"  Unadjusted: log2FC = {lfc_unadj:+.4f}, P = {p_unadj:.5f}")
    print(f"  Adjusted for 5 Lineage Scores: log2FC = {lfc_adj:+.4f}, t = {t_adj:+.4f}, P = {p_adj:.5f}")
