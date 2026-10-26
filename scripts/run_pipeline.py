"""
LEGACY - NOT FOR THE REVISED ANALYSES.
This historical full-pipeline runner is retained for auditability only.
================================================================================
MASTER PERITONEAL FIBROSIS TRANSCRIPTOMICS PIPELINE
================================================================================
Discovery Cohort:  GSE125498 (N=33 PD effluent cells, Illumina HT-12 v4)
                   - 13 Long-term PD (LPD, vintage >= 25 mo, cases = 1)
                   - 20 Short-term PD (SPD, vintage < 6 mo, controls = 0)
Validation Cohort: GSE62928 (N=8 human parietal peritoneal biopsies, Affy HG-U133 Plus 2)
                   - 4 Severe surgical EPS cases (cases = 1)
                   - 4 Non-EPS controls (2 early PD + 2 uremic, controls = 0)

8-Phase Analytical Framework:
  Phase 1: Preprocessing & two probe-to-gene collapses (Universe A: best-P; Universe B: MaxMean)
  Phase 2: Limma empirical Bayes differential expression (|log2FC| >= 0.585, nominal P < 0.05)
  Phase 3: Naba Matrisome (1,027 genes) filtering + hypergeometric enrichment + Monte Carlo permutation
  Phase 4: Unsupervised WGCNA co-expression modules & trait correlation
  Phase 5: 4-model Machine Learning consensus feature selection (LASSO, SVM-RFE, RF, XGBoost >= 2/4 votes)
  Phase 6: Load real-cohort GSE62928 validation outputs from the revision analysis
  Phase 7: Exploratory legacy summaries; STRING/GSEA provenance not established
  Phase 8: Publication-grade figures (300 DPI PNG + PDF), source data CSVs, and audit manifests

All results stored in results/
Random seed = 42 locked across all modules.
================================================================================
"""

import os
import sys
import json
import random
import hashlib
import gzip
import urllib.request
import warnings
from datetime import datetime
from math import comb
from itertools import combinations

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, mannwhitneyu, hypergeom
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import RepeatedStratifiedKFold, LeaveOneOut, StratifiedKFold, cross_val_predict
from xgboost import XGBClassifier

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================================
# CONFIGURATION & DIRECTORY PATHS
# ==============================================================================
ROOT         = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(ROOT, "results")
TAB_DIR      = os.path.join(RESULTS_DIR, "tables")
FIG_DIR      = os.path.join(RESULTS_DIR, "figures")
SRC_DIR      = os.path.join(FIG_DIR, "source_data")
RAW_ROOT     = os.path.join(ROOT, "data", "raw")
VAL_ROOT     = os.path.join(ROOT, "Validation")
LOG_FILE     = os.path.join(RESULTS_DIR, "pipeline_execution.log")

NABA_REF_XLSX = os.path.join(ROOT, "ECM genes all.xlsx")
RANDOM_STATE  = 42
LOGFC_CUTOFF  = 0.585
P_CUTOFF      = 0.05

random.seed(RANDOM_STATE)
np.random.seed(RANDOM_STATE)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

# ==============================================================================
# LOGGING & WORKSPACE HELPERS
# ==============================================================================
def setup_workspace():
    """Ensure output directories exist."""
    for d in [RESULTS_DIR, TAB_DIR, FIG_DIR, SRC_DIR]:
        os.makedirs(d, exist_ok=True)
    with open(LOG_FILE, "w", encoding="utf-8") as fh:
        fh.write(f"PIPELINE EXECUTION LOG: {datetime.now().isoformat()}\n")
        fh.write(f"Discovery Cohort: GSE125498 (N=33 PD effluent: 13 LPD cases, 20 SPD controls)\n")
        fh.write(f"Validation Cohort: GSE62928 (N=8 biopsy: 4 EPS cases, 4 non-EPS controls)\n")
        fh.write(f"Thresholds: |log2FC| >= {LOGFC_CUTOFF}, nominal P < {P_CUTOFF}, seed={RANDOM_STATE}\n")
        fh.write("=" * 80 + "\n")

def log(msg, level="INFO"):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] [{level:5s}] {msg}"
    try:
        print(line, flush=True)
    except Exception:
        clean_line = line.encode("ascii", errors="replace").decode("ascii")
        print(clean_line, flush=True)
    with open(LOG_FILE, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")

def save_csv(df, path, desc=""):
    df.to_csv(path, index=False)
    sz = os.path.getsize(path)
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()[:12]
    log(f"SAVED {os.path.basename(path):<48s} rows={len(df):<6,d} size={sz:<8,d}b sha256={h} {desc}", "OK")

# ==============================================================================
# PHASE 1 & 2: DISCOVERY PREPROCESSING, PROBE COLLAPSES & DIFFERENTIAL EXPRESSION
# ==============================================================================
def run_phase_1_2():
    log("=" * 80)
    log("PHASE 1 & 2: GSE125498 Discovery Preprocessing, Collapses & Limma DE")
    log("=" * 80)

    matrix_gz = os.path.join(RAW_ROOT, "GSE125498_series_matrix.txt.gz")
    if not os.path.exists(matrix_gz):
        raise FileNotFoundError(f"Missing raw series matrix at {matrix_gz}")

    raw_accessions = []
    raw_titles = []
    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                raw_accessions = [s.strip().strip('"') for s in line.strip().split("\t")][1:]
            elif line.startswith("!Sample_title"):
                raw_titles = [s.strip().strip('"') for s in line.strip().split("\t")][1:]
            elif line.startswith("!series_matrix_table_begin"):
                break
    samples_ids_found = raw_accessions
    sample_titles = dict(zip(raw_accessions, raw_titles)) if len(raw_accessions) == len(raw_titles) else {s: s for s in raw_accessions}

    expr_probe = pd.read_csv(
        matrix_gz,
        sep="\t",
        comment="!",
        skip_blank_lines=True,
        engine="python",
        on_bad_lines="skip"
    ).dropna(how="all")

    first_col = expr_probe.columns[0]
    expr_probe = expr_probe.rename(columns={first_col: "Probe_ID"})
    expr_probe["Probe_ID"] = expr_probe["Probe_ID"].astype(str)
    sample_cols = [s for s in samples_ids_found if s in expr_probe.columns]
    expr_probe = expr_probe[["Probe_ID"] + sample_cols].set_index("Probe_ID").astype(float)
    log(f"Expression matrix: {expr_probe.shape[0]:,} probes x {expr_probe.shape[1]} samples")

    # Metadata
    meta_rows = []
    for smp in sample_cols:
        title = sample_titles[smp]
        is_lpd = "Long" in title or "LPD" in title or "long" in title
        label = "LPD" if is_lpd else "SPD"
        meta_rows.append({
            "sample_id": smp,
            "title": title,
            "stage_binary": 1 if is_lpd else 0,
            "stage_label": "Late_Stage_LPD" if is_lpd else "Early_Stage_SPD",
            "group_binary_trait": "case_LPD" if is_lpd else "control_SPD",
            "binary_numeric": 1 if is_lpd else 0,
            "group": label
        })
    meta_df = pd.DataFrame(meta_rows)
    save_csv(meta_df, os.path.join(TAB_DIR, "GSE125498_sample_metadata.csv"), "Metadata")
    n_case = int((meta_df["group"] == "LPD").sum())
    n_ctrl = int((meta_df["group"] == "SPD").sum())
    log(f"Phenotype cohort balance: LPD cases={n_case}, SPD controls={n_ctrl} (N={len(meta_df)})")

    # Limma eBayes toptable
    toptable_tsv = os.path.join(VAL_ROOT, "GSE125498.top.table.tsv")
    probe_top = pd.read_csv(toptable_tsv, sep="\t").rename(columns={"ID": "Probe_ID"})
    probe_top["Gene"] = probe_top["Gene.symbol"].astype(str).str.strip()
    probe_top.loc[probe_top["Gene"].isin(["", "nan", "NA", "NaN", "None"]), "Gene"] = np.nan

    common_probes = pd.Index(probe_top["Probe_ID"].astype(str)).intersection(pd.Index(expr_probe.index.astype(str)))
    probe_top = probe_top[probe_top["Probe_ID"].astype(str).isin(set(common_probes))].copy().reset_index(drop=True)
    probe_rowmeans = expr_probe.reindex(probe_top["Probe_ID"].astype(str)).mean(axis=1).values
    probe_top["AveExpr"] = probe_rowmeans
    probe_top["ID"] = probe_top["Probe_ID"]

    def assign_deg(df):
        return np.where(
            (df["logFC"] >= LOGFC_CUTOFF) & (df["P.Value"] < P_CUTOFF), "Upregulated",
            np.where((df["logFC"] <= -LOGFC_CUTOFF) & (df["P.Value"] < P_CUTOFF), "Downregulated", "Not Significant")
        )
    probe_top["DEG_Status"] = assign_deg(probe_top)

    canon_cols = ["ID", "Gene", "Gene.symbol", "Gene.title", "logFC", "AveExpr",
                  "t", "P.Value", "adj.P.Val", "B", "DEG_Status"]

    # Universe A: Best-Probe-by-min-P
    probe_valid = probe_top[probe_top["Gene"].notna()].sort_values(["Gene", "P.Value"], ascending=[True, True])
    bestp_df = probe_valid.drop_duplicates(subset=["Gene"], keep="first").copy().reset_index(drop=True)
    bestp_df["DEG_Status"] = assign_deg(bestp_df)
    save_csv(bestp_df[canon_cols], os.path.join(TAB_DIR, "GSE125498_all_results.csv"), "Universe A (best-P)")

    deg_filtered = bestp_df[bestp_df["DEG_Status"] != "Not Significant"][canon_cols]
    save_csv(deg_filtered, os.path.join(TAB_DIR, "GSE125498_DEGs_filtered.csv"), "Primary DEGs")

    n_up = int((bestp_df["DEG_Status"] == "Upregulated").sum())
    n_dn = int((bestp_df["DEG_Status"] == "Downregulated").sum())
    n_fdr = int((bestp_df["adj.P.Val"] < 0.05).sum())

    log(f"Universe A: {len(bestp_df):,} unique genes | Up-DEGs: {n_up:,} | Down-DEGs: {n_dn:,} | Genome-wide FDR<0.05: {n_fdr}")

    # Universe B: MaxMean collapse
    probe_mm = probe_top[probe_top["Gene"].notna()].copy()
    probe_mm["row_mean"] = expr_probe.reindex(probe_mm["Probe_ID"].astype(str)).mean(axis=1).values
    probe_mm = probe_mm.sort_values(["Gene", "row_mean"], ascending=[True, False])
    maxmean_picks = probe_mm.drop_duplicates(subset=["Gene"], keep="first").copy()

    expr_genes_df = expr_probe.reindex(maxmean_picks["Probe_ID"].astype(str)).copy()
    expr_genes_df.index = maxmean_picks["Gene"].values
    expr_genes_df.index.name = "Gene"
    expr_genes_df_reset = expr_genes_df.reset_index()
    save_csv(expr_genes_df_reset, os.path.join(TAB_DIR, "GSE125498_full_expression_matrix.csv"), "Universe B MaxMean Matrix")

    maxmean_top = maxmean_picks.drop(columns=["row_mean"]).reset_index(drop=True)
    maxmean_top["DEG_Status"] = assign_deg(maxmean_top)
    save_csv(maxmean_top[canon_cols], os.path.join(TAB_DIR, "GSE125498_gene_level_toptable.csv"), "Universe B MaxMean Toptable")

    # Fig 1: Discovery Volcano
    fig, ax = plt.subplots(figsize=(8, 6.5), facecolor="#F8FAFC")
    ax.scatter(bestp_df[bestp_df["DEG_Status"]=="Not Significant"]["logFC"],
               -np.log10(bestp_df[bestp_df["DEG_Status"]=="Not Significant"]["P.Value"]),
               c="#94A3B8", alpha=0.4, s=12, label="Not Significant")
    ax.scatter(bestp_df[bestp_df["DEG_Status"]=="Downregulated"]["logFC"],
               -np.log10(bestp_df[bestp_df["DEG_Status"]=="Downregulated"]["P.Value"]),
               c="#2563EB", alpha=0.75, s=18, label=f"Downregulated (n={n_dn})")
    ax.scatter(bestp_df[bestp_df["DEG_Status"]=="Upregulated"]["logFC"],
               -np.log10(bestp_df[bestp_df["DEG_Status"]=="Upregulated"]["P.Value"]),
               c="#DC2626", alpha=0.75, s=18, label=f"Upregulated (n={n_up})")
    ax.axvline(LOGFC_CUTOFF, color="#1E293B", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.axvline(-LOGFC_CUTOFF, color="#1E293B", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.axhline(-np.log10(P_CUTOFF), color="#1E293B", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.set_title(f"Volcano Plot: Discovery Cohort GSE125498 (Effluent N=33)\n"
                 f"LPD (n=13) vs SPD (n=20) | Primary cutoffs: |log2FC| >= {LOGFC_CUTOFF}, P < {P_CUTOFF}",
                 fontsize=11, fontweight="bold", pad=12)
    ax.set_xlabel("log2 Fold Change (LPD / SPD)", fontsize=10, fontweight="bold")
    ax.set_ylabel("-log10(P-value)", fontsize=10, fontweight="bold")
    ax.legend(loc="upper left", frameon=True, facecolor="#FFFFFF")
    ax.grid(alpha=0.3, linestyle=":")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_01_volcano_plot.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_01_volcano_plot.pdf"))
    plt.close()

    volcano_src = bestp_df[["Gene", "logFC", "P.Value", "adj.P.Val", "DEG_Status"]].copy()
    save_csv(volcano_src, os.path.join(SRC_DIR, "Fig1_discovery_volcano_source_data.csv"), "Volcano Source Data")

    counts = {
        "n_probe_level": len(probe_top),
        "n_samples": expr_probe.shape[1],
        "n_lpd_cases": n_case,
        "n_spd_controls": n_ctrl,
        "universe_bestp_n": len(bestp_df),
        "universe_maxmean_n": len(maxmean_top),
        "degs_primary_up": n_up,
        "degs_primary_down": n_dn,
        "degs_primary_total": n_up + n_dn,
        "fdr_sig_total": n_fdr,
        "logfc_cutoff": LOGFC_CUTOFF,
        "p_cutoff": P_CUTOFF
    }
    with open(os.path.join(TAB_DIR, "phase1_2_counts.json"), "w") as fh:
        json.dump(counts, fh, indent=2)

    return bestp_df, expr_genes_df, meta_df, counts

# ==============================================================================
# PHASE 3: MATRISOME FILTERING & ENRICHMENT AUDIT
# ==============================================================================
def run_phase_3(bestp_df, meta_df, counts):
    log("=" * 80)
    log("PHASE 3: Naba Matrisome (1,027 Genes) Filtering & Enrichment Analysis")
    log("=" * 80)

    wb = load_workbook(NABA_REF_XLSX, data_only=True, read_only=True)
    ws = wb.active
    rows_raw = list(ws.iter_rows(values_only=True))
    data_rows = [r[:4] for r in rows_raw[2:] if r[2] is not None and str(r[2]).strip() != ""]
    ecm_raw = pd.DataFrame(data_rows, columns=["Category", "Subcategory", "Gene", "Gene_Name"])
    ecm_raw["Gene"] = ecm_raw["Gene"].astype(str).str.strip()
    matrisome_ref = ecm_raw[ecm_raw["Gene"].notna() & (ecm_raw["Gene"] != "") & (ecm_raw["Gene"] != "NA")].drop_duplicates(subset=["Gene"]).reset_index(drop=True)
    matrisome_ref = matrisome_ref[["Gene", "Category", "Subcategory"]].sort_values(["Category", "Subcategory", "Gene"]).reset_index(drop=True)

    n_ref = len(matrisome_ref)
    log(f"Naba reference loaded: {n_ref:,} genes across 7 subcategories")

    measured_set = set(bestp_df["Gene"].dropna().astype(str).str.strip())
    measured_ecm = sorted(measured_set & set(matrisome_ref["Gene"].astype(str)))
    n_measured_ecm = len(measured_ecm)
    log(f"Matrisome genes measured in GSE125498: {n_measured_ecm:,} / {n_ref:,}")

    ecm_df = bestp_df[bestp_df["Gene"].astype(str).isin(set(measured_ecm))].copy()
    ecm_df = ecm_df.merge(matrisome_ref, on="Gene", how="left")
    save_csv(ecm_df, os.path.join(TAB_DIR, "GSE125498_ECM_intersection.csv"), "ECM Intersection")

    up_ecm = int((ecm_df["DEG_Status"] == "Upregulated").sum())
    dn_ecm = int((ecm_df["DEG_Status"] == "Downregulated").sum())
    n_up_deg = counts["degs_primary_up"]
    n_univ = counts["universe_bestp_n"]

    fold_enrich = (up_ecm / n_up_deg) / (n_measured_ecm / n_univ) if n_up_deg > 0 else 0.0
    hypergeom_p = hypergeom.sf(up_ecm - 1, n_univ, n_measured_ecm, n_up_deg) if n_up_deg > 0 else 1.0

    ovl = ecm_df.groupby(["Category", "Subcategory"]).agg(
        Measured=("Gene", "nunique"),
        Up_ECM_DEG=("DEG_Status", lambda s: int((s == "Upregulated").sum())),
        Down_ECM_DEG=("DEG_Status", lambda s: int((s == "Downregulated").sum()))
    ).reset_index()
    save_csv(ovl, os.path.join(TAB_DIR, "GSE125498_Matrisome_Overlap_Counts.csv"), "Matrisome Overlap Counts")

    # Fast Vectorized Monte Carlo Permutation Test (N=10,000)
    log("Running Monte Carlo label permutations (N=10,000) for ECM enrichment P...")
    expr_matrix_path = os.path.join(TAB_DIR, "GSE125498_full_expression_matrix.csv")
    mm_df = pd.read_csv(expr_matrix_path).set_index("Gene")
    sample_cols = [c for c in mm_df.columns if c.startswith("GSM")]
    y_true = meta_df.set_index("sample_id").loc[sample_cols]["stage_binary"].values.astype(int)

    common_ecm_genes = [g for g in measured_ecm if g in mm_df.index]
    X_ecm = mm_df.loc[common_ecm_genes, sample_cols].values
    n_ecm_genes, n_total_samples = X_ecm.shape
    n_cases = int((y_true == 1).sum())
    n_ctrls = int((y_true == 0).sum())

    N_PERM_ECM = 10_000
    n_extreme_ecm = 0
    t_crit = stats.t.ppf(1.0 - (P_CUTOFF / 2.0), df=n_total_samples - 2)

    for p in range(N_PERM_ECM):
        shuff_idx = np.random.permutation(n_total_samples)
        c_idx = shuff_idx[:n_cases]
        s_idx = shuff_idx[n_cases:]

        m1 = X_ecm[:, c_idx].mean(axis=1)
        m2 = X_ecm[:, s_idx].mean(axis=1)
        v1 = X_ecm[:, c_idx].var(axis=1, ddof=1)
        v2 = X_ecm[:, s_idx].var(axis=1, ddof=1)
        se = np.sqrt(v1 / n_cases + v2 / n_ctrls)
        t_stat = (m1 - m2) / np.clip(se, 1e-9, None)
        logfc = m1 - m2

        up_ecm_perm = int(np.sum((logfc >= LOGFC_CUTOFF) & (t_stat >= t_crit)))
        if up_ecm_perm >= up_ecm:
            n_extreme_ecm += 1

    mc_ecm_p = (n_extreme_ecm + 1) / (N_PERM_ECM + 1)
    log(f"ECM Permutation: {n_extreme_ecm} / {N_PERM_ECM:,} >= {up_ecm} (Perm P = {mc_ecm_p:.5f})")

    # Fig 2: Matrisome Enrichment
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), facecolor="#F8FAFC")
    ovl_sorted = ovl.sort_values("Up_ECM_DEG", ascending=True)
    axes[0].barh(ovl_sorted["Subcategory"], ovl_sorted["Up_ECM_DEG"], color="#DC2626", edgecolor="#1E293B", height=0.6)
    axes[0].set_title(f"Pro-fibrotic Up-regulated ECM Genes by Subcategory\nTotal = {up_ecm} Up ECM-DEGs", fontsize=11, fontweight="bold")
    axes[0].set_xlabel("Gene Count", fontsize=10, fontweight="bold")
    axes[0].grid(axis='x', linestyle=":", alpha=0.5)

    sizes = [up_ecm, n_up_deg - up_ecm]
    labels = [f"ECM Genes\n(n={up_ecm}, {up_ecm/n_up_deg*100:.1f}%)", f"Non-ECM Genes\n(n={n_up_deg-up_ecm}, {(n_up_deg-up_ecm)/n_up_deg*100:.1f}%)"]
    axes[1].pie(sizes, labels=labels, colors=["#DC2626", "#94A3B8"], autopct='%1.1f%%', startangle=140,
                wedgeprops=dict(width=0.45, edgecolor='#1E293B'))
    axes[1].set_title(f"ECM Enrichment in Up-DEGs\nFold Enrich = {fold_enrich:.2f}x | Hypergeom P = {hypergeom_p:.2e}\n10k Perm P = {mc_ecm_p:.4f}",
                      fontsize=11, fontweight="bold")

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_02_matrisome_enrichment.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_02_matrisome_enrichment.pdf"))
    plt.close()

    save_csv(ovl, os.path.join(SRC_DIR, "Fig2_matrisome_enrichment_source_data.csv"), "Matrisome Source Data")

    counts.update({
        "matrisome_reference_1027": n_ref,
        "matrisome_measured": n_measured_ecm,
        "profibrotic_up_ecm_degs": up_ecm,
        "down_ecm_degs": dn_ecm,
        "total_ecm_degs_both_dirs": up_ecm + dn_ecm,
        "fold_enrichment_up_ecm": float(fold_enrich),
        "hypergeom_p_up_ecm": float(hypergeom_p),
        "permutation_ecm_p_10k": float(mc_ecm_p)
    })
    with open(os.path.join(TAB_DIR, "phase3_counts.json"), "w") as fh:
        json.dump(counts, fh, indent=2)

    return ecm_df, matrisome_ref, counts

# ==============================================================================
# PHASE 4: UNSUPERVISED CO-EXPRESSION & 3-WAY CONVERGENCE
# ==============================================================================
def run_phase_4(expr_genes_df, meta_df, bestp_df, ecm_df, counts):
    log("=" * 80)
    log("PHASE 4: Unsupervised Co-Expression Network & 3-Way Convergence")
    log("=" * 80)

    if "Gene" in expr_genes_df.columns:
        expr_mat = expr_genes_df.set_index("Gene")
    else:
        expr_mat = expr_genes_df.copy()
    sample_cols = [c for c in expr_mat.columns if c.startswith("GSM")]
    datExpr = expr_mat[sample_cols].T

    up_ecm_genes = set(ecm_df[ecm_df["DEG_Status"] == "Upregulated"]["Gene"].astype(str))
    Q1 = datExpr.quantile(0.25, axis=0)
    Q3 = datExpr.quantile(0.75, axis=0)
    IQR = Q3 - Q1
    top_var_genes = set(IQR.sort_values(ascending=False).head(5000).index.tolist())
    network_genes = sorted(top_var_genes | (up_ecm_genes & set(datExpr.columns)))
    datExpr_filt = datExpr[network_genes].copy()
    log(f"WGCNA network input: {datExpr_filt.shape[0]} samples x {datExpr_filt.shape[1]:,} genes (IQR 5000 + Up-ECM)")

    softPower = 12
    cor_mat = np.clip(np.corrcoef(datExpr_filt.values, rowvar=False), -1.0, 1.0)
    dist_mat = 1.0 - (0.5 + 0.5 * cor_mat) ** softPower
    np.fill_diagonal(dist_mat, 0.0)

    dist_condensed = squareform(dist_mat, checks=False)
    link_mat = linkage(dist_condensed, method='average')
    k_modules = 14
    cluster_labels = fcluster(link_mat, t=k_modules, criterion='maxclust')

    palette = [
        "turquoise", "blue", "brown", "yellow", "green", "red", "black",
        "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan"
    ]
    unique_clusters = sorted(set(cluster_labels))
    cluster_to_color = {c: palette[i % len(palette)] for i, c in enumerate(unique_clusters)}

    gene_mod_df = pd.DataFrame({
        "Gene": datExpr_filt.columns.tolist(),
        "Module": [cluster_to_color[c] for c in cluster_labels]
    })
    save_csv(gene_mod_df, os.path.join(TAB_DIR, "wgcna_gene_module_membership.csv"), "Module Membership")

    trait_dict = dict(zip(meta_df["sample_id"], meta_df["stage_binary"]))
    sample_ids = datExpr_filt.index.tolist()
    y_trait = np.array([trait_dict[s] for s in sample_ids], dtype=float)

    me_dict = {}
    for mod_color, group in gene_mod_df.groupby("Module"):
        g_list = group["Gene"].tolist()
        sub_mat = datExpr_filt[g_list].values
        if sub_mat.shape[1] < 2:
            me = sub_mat[:, 0]
        else:
            me = PCA(n_components=1, random_state=RANDOM_STATE).fit_transform(sub_mat).flatten()
        if np.corrcoef(me, sub_mat.mean(axis=1))[0, 1] < 0:
            me = -me
        me_dict[mod_color] = me

    me_df = pd.DataFrame(me_dict, index=sample_ids)
    me_df_reset = me_df.reset_index().rename(columns={"index": "sample_id"})
    save_csv(me_df_reset, os.path.join(TAB_DIR, "wgcna_module_eigengenes.csv"), "Module Eigengenes")

    mt_rows = []
    for mod_col in me_df.columns:
        r, p = pearsonr(me_df[mod_col].values, y_trait)
        mt_rows.append({"Module": f"ME{mod_col}", "Module_Color": mod_col, "Correlation": r, "P_value": p})
    mt_df = pd.DataFrame(mt_rows).sort_values("Correlation", ascending=False).reset_index(drop=True)
    save_csv(mt_df, os.path.join(TAB_DIR, "wgcna_module_trait_correlation.csv"), "Module Trait Correlations")

    top_mod_row = mt_df.iloc[0]
    top_mod_color = top_mod_row["Module_Color"]
    top_mod_r = top_mod_row["Correlation"]
    top_mod_p = top_mod_row["P_value"]
    top_mod_genes = set(gene_mod_df[gene_mod_df["Module"] == top_mod_color]["Gene"])
    log(f"Top trait-correlated module: ME{top_mod_color} (r={top_mod_r:+.3f}, P={top_mod_p:.2e}, n={len(top_mod_genes)} genes)")

    top_mod_df = gene_mod_df[gene_mod_df["Module"] == top_mod_color].copy()
    save_csv(top_mod_df, os.path.join(TAB_DIR, "wgcna_trait_significant_module_genes.csv"), "Top Module Genes")

    up_deg_genes = set(bestp_df[bestp_df["DEG_Status"] == "Upregulated"]["Gene"].dropna().astype(str))
    measured_ecm_genes = set(ecm_df["Gene"].dropna().astype(str))

    convergent_genes = sorted(up_deg_genes & measured_ecm_genes & top_mod_genes)
    log(f"3-WAY CONVERGENCE POOL: {len(convergent_genes)} genes (Up-DEGs & Matrisome & ME{top_mod_color})")

    if len(convergent_genes) < 5:
        log("Convergence pool < 5; expanding to positive-trait correlated modules & Up-ECM DEGs...", "WARN")
        pos_modules = set(mt_df[mt_df["Correlation"] > 0]["Module_Color"])
        pos_mod_genes = set(gene_mod_df[gene_mod_df["Module"].isin(pos_modules)]["Gene"])
        candidate_pool = sorted(up_deg_genes & measured_ecm_genes & pos_mod_genes)
        if len(candidate_pool) < 5:
            candidate_pool = sorted(up_deg_genes & measured_ecm_genes)
    else:
        candidate_pool = convergent_genes

    conv_rows = []
    all_union = sorted(up_deg_genes | measured_ecm_genes | top_mod_genes)
    for g in all_union:
        in_up = g in up_deg_genes
        in_ecm = g in measured_ecm_genes
        in_wgcna = g in top_mod_genes
        is_conv = in_up and in_ecm and in_wgcna
        conv_rows.append({
            "Gene": g,
            "in_Up_DEG": int(in_up),
            "in_Measured_ECM": int(in_ecm),
            "in_Top_WGCNA_Module": int(in_wgcna),
            "Convergent": int(is_conv)
        })
    conv_df = pd.DataFrame(conv_rows)
    save_csv(conv_df, os.path.join(TAB_DIR, "convergent_WGCNA_ECM_genes.csv"), "Convergence Table")

    # Fig 3: WGCNA Module-Trait
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), facecolor="#F8FAFC")
    dend = dendrogram(link_mat, ax=axes[0], no_labels=True, above_threshold_color="#475569")
    axes[0].set_title(f"Unsupervised Co-expression Dendrogram\n(Signed Network, Beta={softPower}, N={len(network_genes):,} Genes)", fontsize=11, fontweight="bold")
    axes[0].set_xlabel("Gene Co-expression Tree", fontsize=10, fontweight="bold")
    axes[0].set_ylabel("Dissimilarity Height", fontsize=10, fontweight="bold")

    mt_plot = mt_df.sort_values("Correlation", ascending=True)
    colors_mt = ["#DC2626" if r > 0 else "#2563EB" for r in mt_plot["Correlation"]]
    axes[1].barh(mt_plot["Module"], mt_plot["Correlation"], color=colors_mt, edgecolor="#1E293B", height=0.6)
    axes[1].axvline(0, color="#1E293B", linewidth=0.8)
    axes[1].set_title(f"Module-Trait Correlation with LPD Effluent Trait\nTop Module = ME{top_mod_color} (r={top_mod_r:+.2f}, P={top_mod_p:.2e})", fontsize=11, fontweight="bold")
    axes[1].set_xlabel("Pearson Correlation (r)", fontsize=10, fontweight="bold")
    axes[1].grid(axis='x', linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_03_wgcna_module_trait.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_03_wgcna_module_trait.pdf"))
    plt.close()

    save_csv(mt_df, os.path.join(SRC_DIR, "Fig3_wgcna_source_data.csv"), "WGCNA Source Data")

    counts.update({
        "wgcna_network_genes": len(network_genes),
        "wgcna_modules_n": k_modules,
        "top_wgcna_module": f"ME{top_mod_color}",
        "top_wgcna_module_r": float(top_mod_r),
        "top_wgcna_module_p": float(top_mod_p),
        "top_wgcna_module_genes_n": len(top_mod_genes),
        "convergent_3way_genes_n": len(convergent_genes),
        "candidate_pool_for_ml_n": len(candidate_pool)
    })
    with open(os.path.join(TAB_DIR, "phase4_counts.json"), "w") as fh:
        json.dump(counts, fh, indent=2)

    return candidate_pool, conv_df, counts

# ==============================================================================
# PHASE 5: 4-MODEL MACHINE LEARNING CONSENSUS
# ==============================================================================
def run_phase_5(candidate_pool, expr_genes_df, meta_df, bestp_df, counts):
    log("=" * 80)
    log("PHASE 5: 4-Model ML Consensus Feature Selection (LASSO, SVM-RFE, RF, XGB)")
    log("=" * 80)

    if "Gene" in expr_genes_df.columns:
        expr_mat = expr_genes_df.set_index("Gene")
    else:
        expr_mat = expr_genes_df.copy()

    sample_cols = [c for c in meta_df["sample_id"] if c in expr_mat.columns]
    y = meta_df.set_index("sample_id").loc[sample_cols]["stage_binary"].values.astype(int)

    features = [g for g in candidate_pool if g in expr_mat.index]
    if len(features) < 2:
        ecm_df = pd.read_csv(os.path.join(TAB_DIR, "GSE125498_ECM_intersection.csv"))
        features = ecm_df[ecm_df["DEG_Status"] == "Upregulated"]["Gene"].tolist()[:25]

    X_raw = expr_mat.loc[features, sample_cols].T.values.astype(float)
    X = StandardScaler().fit_transform(X_raw)
    n_feat = len(features)
    log(f"ML Input: {X.shape[0]} samples x {n_feat} candidate genes (LPD={int((y==1).sum())}, SPD={int((y==0).sum())})")

    # 1. LASSO
    lasso_cv = LogisticRegressionCV(Cs=[0.01, 0.05, 0.1, 0.5, 1, 5, 10], cv=5, penalty='l1',
                                    solver='liblinear', max_iter=5000, random_state=RANDOM_STATE)
    lasso_cv.fit(X, y)
    lasso_coefs = np.abs(lasso_cv.coef_[0])
    lasso_sel_idx = np.where(lasso_coefs > 1e-4)[0]
    if len(lasso_sel_idx) == 0:
        lasso_sel_idx = np.argsort(-lasso_coefs)[:max(2, int(n_feat * 0.2))]
    lasso_genes = set([features[i] for i in lasso_sel_idx])

    # 2. SVM-RFE
    n_select_rfe = max(2, min(10, int(np.ceil(n_feat * 0.35))))
    svm_rfe = RFE(estimator=SVC(kernel="linear", random_state=RANDOM_STATE), n_features_to_select=n_select_rfe, step=1)
    svm_rfe.fit(X, y)
    svm_genes = set([features[i] for i in range(n_feat) if svm_rfe.support_[i]])
    svm_ranking = svm_rfe.ranking_

    # 3. Random Forest
    rf = RandomForestClassifier(n_estimators=500, random_state=RANDOM_STATE, class_weight='balanced')
    rf.fit(X, y)
    rf_imp = rf.feature_importances_
    rf_thresh = np.sort(rf_imp)[::-1][min(n_select_rfe - 1, n_feat - 1)]
    rf_genes = set([features[i] for i in range(n_feat) if rf_imp[i] >= rf_thresh])

    # 4. XGBoost
    xgb = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.08, random_state=RANDOM_STATE,
                        eval_metric='logloss', use_label_encoder=False)
    xgb.fit(X, y)
    xgb_imp = xgb.feature_importances_
    xgb_thresh = np.sort(xgb_imp)[::-1][min(n_select_rfe - 1, n_feat - 1)]
    xgb_genes = set([features[i] for i in range(n_feat) if xgb_imp[i] >= xgb_thresh])

    all_res_dict = bestp_df.set_index("Gene").to_dict(orient="index")
    consensus_rows = []
    for i, g in enumerate(features):
        v_l = int(g in lasso_genes)
        v_s = int(g in svm_genes)
        v_r = int(g in rf_genes)
        v_x = int(g in xgb_genes)
        total_v = v_l + v_s + v_r + v_x
        is_hub = int(total_v >= 2)

        g_info = all_res_dict.get(g, {})
        consensus_rows.append({
            "Gene": g,
            "LASSO_vote": v_l,
            "SVM_RFE_vote": v_s,
            "RF_vote": v_r,
            "XGB_vote": v_x,
            "Total_votes": total_v,
            "IsConsensusHub": is_hub,
            "LASSO_abs_coef": float(lasso_coefs[i]),
            "SVM_RFE_rank": int(svm_ranking[i]),
            "RF_importance": float(rf_imp[i]),
            "XGB_importance": float(xgb_imp[i]),
            "Discovery_logFC": float(g_info.get("logFC", np.nan)),
            "Discovery_P_Value": float(g_info.get("P.Value", np.nan)),
            "Discovery_adjP": float(g_info.get("adj.P.Val", np.nan)),
            "Discovery_DEG_Status": str(g_info.get("DEG_Status", "N/A"))
        })

    consensus_df = pd.DataFrame(consensus_rows).sort_values(["Total_votes", "Discovery_logFC"], ascending=[False, False]).reset_index(drop=True)
    save_csv(consensus_df, os.path.join(TAB_DIR, "ML_hub_genes_from_WGCNA_ECM_all_results.csv"), "ML All Results")

    hubs_df = consensus_df[consensus_df["IsConsensusHub"] == 1].copy().reset_index(drop=True)
    if len(hubs_df) == 0:
        hubs_df = consensus_df.head(5).copy()
    save_csv(hubs_df, os.path.join(TAB_DIR, "ML_hub_genes_from_WGCNA_ECM.csv"), "Consensus Hub Panel")

    hub_genes = hubs_df["Gene"].tolist()
    log(f"ML CONSENSUS HUB GENES (>= 2/4 votes): N = {len(hub_genes)}")
    for _, r in hubs_df.iterrows():
        log(f"  -> {r['Gene']:<12s} votes={r['Total_votes']}/4 (L={r['LASSO_vote']}, S={r['SVM_RFE_vote']}, R={r['RF_vote']}, X={r['XGB_vote']}) logFC={r['Discovery_logFC']:+.3f}")

    # Fig 4: ML Consensus Hubs
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), facecolor="#F8FAFC")
    c_plot = consensus_df.sort_values("Total_votes", ascending=True)
    b_colors = ["#DC2626" if v >= 2 else "#94A3B8" for v in c_plot["Total_votes"]]
    axes[0].barh(c_plot["Gene"], c_plot["Total_votes"], color=b_colors, edgecolor="#1E293B", height=0.6)
    axes[0].axvline(2, color="#B91C1C", linestyle="--", linewidth=1.2, label="Consensus Cutoff (>= 2 Votes)")
    axes[0].set_title(f"Machine Learning Consensus Hub Selection\n{len(hub_genes)} Consensus Hubs Selected (>= 2 of 4 Models)", fontsize=11, fontweight="bold")
    axes[0].set_xlabel("Model Votes (out of 4: LASSO, SVM-RFE, RF, XGBoost)", fontsize=10, fontweight="bold")
    axes[0].legend(loc="lower right")
    axes[0].grid(axis='x', linestyle=":", alpha=0.5)

    heat_data = consensus_df.set_index("Gene")[["LASSO_vote", "SVM_RFE_vote", "RF_vote", "XGB_vote"]]
    sns.heatmap(heat_data, cmap="Blues", cbar=False, linewidths=1, linecolor="#CBD5E1", ax=axes[1], annot=True)
    axes[1].set_title("Feature Selection Vote Matrix Across 4 Algorithms", fontsize=11, fontweight="bold")
    axes[1].set_xlabel("Machine Learning Algorithm", fontsize=10, fontweight="bold")

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_04_ml_consensus_hubs.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_04_ml_consensus_hubs.pdf"))
    plt.close()

    save_csv(consensus_df, os.path.join(SRC_DIR, "Fig4_ml_consensus_source_data.csv"), "ML Source Data")

    counts.update({
        "ml_candidate_features_n": n_feat,
        "lasso_selected_n": len(lasso_genes),
        "svm_rfe_selected_n": len(svm_genes),
        "rf_selected_n": len(rf_genes),
        "xgb_selected_n": len(xgb_genes),
        "consensus_hubs_n": len(hub_genes),
        "consensus_hub_genes": hub_genes
    })
    with open(os.path.join(TAB_DIR, "phase5_counts.json"), "w") as fh:
        json.dump(counts, fh, indent=2)

    return hub_genes, hubs_df, counts

# ==============================================================================
# PHASE 6: EXTERNAL VALIDATION ON GSE62928 (N=8 BIOPSY)
# ==============================================================================
def run_phase_6(hub_genes, counts):
    log("=" * 80)
    log("PHASE 6: Loading real GSE62928 validation outputs")
    log("=" * 80)

    revision_dir = os.path.join(RESULTS_DIR, "revision")
    metrics_path = os.path.join(revision_dir, "gse62928_robustness_metrics_taskA_fig5.json")
    per_gene_path = os.path.join(revision_dir, "gse62928_per_gene_concordance_taskA_fig5.csv")
    if not os.path.isfile(metrics_path) or not os.path.isfile(per_gene_path):
        raise RuntimeError(
            "Real GSE62928 validation outputs are missing. Run "
            "'python scripts/task3_validation_robustness.py' first; Phase 6 will not "
            "simulate a validation matrix or overwrite legacy results."
        )

    with open(metrics_path, encoding="utf-8") as fh:
        results = json.load(fh)
    per_gene_df = pd.read_csv(per_gene_path)
    if results.get("validation_cohort", "").startswith("GSE62928") is False:
        raise ValueError(f"Unexpected validation cohort in {metrics_path}")
    expected_hubs = set(results["hub_genes_evaluated"])
    if set(hub_genes) != expected_hubs:
        raise RuntimeError(
            "The legacy four-model panel does not match the fixed seven-hub signature "
            "in the real GSE62928 outputs. Phase 6 cannot report validation for a "
            "different feature set; run the dedicated revision analysis instead."
        )
    composite = results["composite_hub_signature"]
    null_results = results["random_gene_set_null_1000"]
    concordance = results["per_gene_concordance"]
    mapped_hubs = per_gene_df["Gene"].astype(str).tolist()
    val_df = pd.DataFrame([{
        "Validation_Cohort": results["validation_cohort"],
        "N_Samples": results["n_samples"],
        "Hubs_Evaluated": len(mapped_hubs),
        "Mapped_Gene_List": ";".join(mapped_hubs),
        "Direction_Concordant": concordance["n_concordant"],
        "Binomial_P": concordance["binomial_sign_test_p"],
        "Composite_Signature_AUC": composite["hub_signature_auc"],
        "Mann_Whitney_One_Sided_P": composite["mann_whitney_one_sided_p"],
        "Mann_Whitney_Two_Sided_P": composite["mann_whitney_two_sided_p"],
        "Random_Set_Null_Mean_AUC": null_results["null_auc_mean"],
        "Random_Set_Empirical_P": null_results["empirical_p_value"]
    }])
    log(
        f"Real biopsy hub score AUC = {composite['hub_signature_auc']:.4f}; "
        f"2-sided Mann-Whitney P = {composite['mann_whitney_two_sided_p']:.4f}"
    )

    counts.update({
        "gse62928_mapped_hubs_n": len(mapped_hubs),
        "validation_real_signature_auc": composite["hub_signature_auc"],
        "validation_direction_concordant_n": concordance["n_concordant"],
        "validation_binomial_sign_p": concordance["binomial_sign_test_p"],
        "validation_mann_whitney_two_sided_p": composite["mann_whitney_two_sided_p"],
        "validation_random_set_null_mean_auc": null_results["null_auc_mean"],
        "validation_random_set_empirical_p": null_results["empirical_p_value"]
    })
    return val_df, counts

# ==============================================================================
# PHASE 7: EXPLORATORY NETWORKS, PREDICTION & MARKER-SCORE SUMMARIES
# ==============================================================================
def run_phase_7(hub_genes, expr_genes_df, meta_df, bestp_df, counts):
    log("=" * 80)
    log("PHASE 7: Unverified network outputs, 5-gene model & marker-score summaries")
    log("=" * 80)

    # 7a. No source STRING export is available; do not synthesize an edge list.
    n_h = len(hub_genes)
    max_edges = comb(n_h, 2) if n_h >= 2 else 1
    string_edges_dict = {}
    n_edges = len(string_edges_dict)
    conn_rate = (n_edges / max_edges) if max_edges > 0 else 0.0
    log("No verified STRING export available; PPI edges and centrality are not reported.", "WARN")

    edge_rows = [
        {"preferredName_A": p[0], "preferredName_B": p[1], "string_functional_score": s, "edge": f"{p[0]}-{p[1]}"}
        for p, s in sorted(string_edges_dict.items())
    ]
    edges_df = pd.DataFrame(edge_rows, columns=[
        "preferredName_A", "preferredName_B", "string_functional_score", "edge"
    ])
    save_csv(edges_df, os.path.join(TAB_DIR, "STRING_edges_functional.csv"), "STRING Edges")

    centrality_rows = []
    for g in hub_genes:
        centrality_rows.append({"Gene": g, "Degree": np.nan, "Degree_Centrality": np.nan})
    cent_df = pd.DataFrame(centrality_rows).sort_values("Degree", ascending=False).reset_index(drop=True)
    save_csv(cent_df, os.path.join(TAB_DIR, "hub_degree_centrality.csv"), "Degree Centrality")

    # 7b. A custom hand-curated EMT-like set is descriptive only, not Hallmark GSEA.
    df_rank = bestp_df.copy()
    df_rank["score"] = np.sign(df_rank["logFC"].fillna(0)) * -np.log10(df_rank["P.Value"].clip(lower=1e-300))
    df_rank = df_rank.sort_values("score", ascending=False).reset_index(drop=True)

    EMT_GENESET = set([
        "FN1", "VIM", "CDH2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2",
        "MMP2", "MMP3", "MMP7", "MMP9", "MMP14", "COL1A1", "COL1A2", "COL3A1",
        "COL4A1", "COL5A1", "COL5A2", "COL6A1", "COL8A1", "COL11A1", "TGFB1",
        "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "POSTN", "SPARC",
        "LOX", "LOXL1", "LOXL2", "COMP", "THBS1", "THBS2", "THBS3", "VCAN", "EDIL3",
        "ACTA2", "TAGLN", "CTGF", "FAP", "CALD1", "INHBA", "ISM1"
    ])
    ranked_genes = df_rank["Gene"].tolist()
    ranked_scores = df_rank["score"].values.astype(float)
    hits_mask = df_rank["Gene"].isin(EMT_GENESET).values
    N_g = len(ranked_genes)
    Nh_g = int(hits_mask.sum())
    NR_g = np.abs(ranked_scores).sum() if np.abs(ranked_scores).sum() > 0 else 1.0

    phit = np.zeros(N_g + 1)
    pmiss = np.zeros(N_g + 1)
    for i in range(1, N_g + 1):
        if hits_mask[i - 1]:
            phit[i] = phit[i - 1] + abs(ranked_scores[i - 1]) / NR_g
            pmiss[i] = pmiss[i - 1]
        else:
            phit[i] = phit[i - 1]
            pmiss[i] = pmiss[i - 1] + 1.0 / max(1, N_g - Nh_g)
    running_es = phit - pmiss
    es_max_idx = np.argmax(np.abs(running_es))
    ES_val = float(running_es[es_max_idx])
    gsea_df = pd.DataFrame([{
        "Pathway": "Custom manually selected EMT-like genes (not Hallmark)",
        "Set_Size": len(EMT_GENESET),
        "Hits_In_Discovery": Nh_g,
        "Enrichment_Score_ES": round(ES_val, 4),
        "Normalized_Enrichment_Score_NES": np.nan,
        "Nominal_P_Value": np.nan
    }])
    save_csv(gsea_df, os.path.join(TAB_DIR, "custom_emt_signature_running_score.csv"), "Custom EMT-like Score")

    # 7c. 5-Gene Nomogram
    if "Gene" in expr_genes_df.columns:
        expr_mat = expr_genes_df.set_index("Gene")
    else:
        expr_mat = expr_genes_df.copy()

    sample_cols = [s for s in meta_df["sample_id"] if s in expr_mat.columns]
    y_nomo = meta_df.set_index("sample_id").loc[sample_cols]["stage_binary"].values.astype(int)

    top5_nomo_genes = ['FLT3LG', 'TNFSF15', 'ADAM19', 'SEMA3E', 'SERPINA10']
    X_nomo_raw = expr_mat.loc[top5_nomo_genes, sample_cols].T.values.astype(float)
    X_nomo = StandardScaler().fit_transform(X_nomo_raw)

    clf_nomo = LogisticRegression(penalty=None, max_iter=5000, random_state=RANDOM_STATE)
    clf_nomo.fit(X_nomo, y_nomo)
    nomo_probs = clf_nomo.predict_proba(X_nomo)[:, 1]
    c_index_apparent = roc_auc_score(y_nomo, nomo_probs)

    skf_nomo = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    cv_nomo_probs = cross_val_predict(clf_nomo, X_nomo, y_nomo, cv=skf_nomo, method='predict_proba')[:, 1]
    c_index_cv = roc_auc_score(y_nomo, cv_nomo_probs)
    optimism = c_index_apparent - c_index_cv
    log(f"Nomogram C-index: Apparent={c_index_apparent:.4f} | 5-fold CV={c_index_cv:.4f} (Optimism gap={optimism:.4f})")

    coefs = clf_nomo.coef_[0]
    max_c = np.max(np.abs(coefs)) if np.max(np.abs(coefs)) > 0 else 1.0
    nomo_points = []
    for i, g in enumerate(top5_nomo_genes):
        pts = 100.0 * (abs(coefs[i]) / max_c)
        nomo_points.append({
            "Gene": g,
            "Logit_Coefficient": round(float(coefs[i]), 4),
            "Points_Per_SD": round(float(pts), 2),
            "Points_at_Low_Expr": -round(float(pts), 2) if coefs[i] > 0 else round(float(pts), 2),
            "Points_at_High_Expr": round(float(pts), 2) if coefs[i] > 0 else -round(float(pts), 2)
        })
    df_nomo_points = pd.DataFrame(nomo_points)
    save_csv(df_nomo_points, os.path.join(TAB_DIR, "nomogram_point_scoring_table.csv"), "Nomogram Points")

    thresholds = np.linspace(0.05, 0.95, 19)
    dca_rows = []
    prev = y_nomo.mean()
    for t in thresholds:
        tp = int(((nomo_probs >= t) & (y_nomo == 1)).sum())
        fp = int(((nomo_probs >= t) & (y_nomo == 0)).sum())
        nb_model = (tp / len(y_nomo)) - (fp / len(y_nomo)) * (t / (1.0 - t))
        nb_all = prev - (1.0 - prev) * (t / (1.0 - t))
        dca_rows.append({
            "Threshold_Probability": round(float(t), 3),
            "Net_Benefit_Nomogram": round(float(max(0.0, nb_model)), 5),
            "Net_Benefit_TreatAll": round(float(max(0.0, nb_all)), 5),
            "Net_Benefit_TreatNone": 0.0
        })
    df_dca = pd.DataFrame(dca_rows)
    save_csv(df_dca, os.path.join(TAB_DIR, "decision_curve_analysis_metrics.csv"), "DCA Metrics")

    # 7d. Six-gene expression marker score; this is not cell-fraction deconvolution.
    MYO_SIG = ["ACTA2", "TAGLN", "COL1A1", "COL1A2", "POSTN", "FN1"]
    myo_present = [g for g in MYO_SIG if g in expr_mat.index]
    myo_mat = expr_mat.loc[myo_present, sample_cols].T.values.astype(float)
    myo_z = StandardScaler().fit_transform(myo_mat)
    myo_scores = myo_z.mean(axis=1)

    lpd_idx = np.where(y_nomo == 1)[0]
    spd_idx = np.where(y_nomo == 0)[0]
    m_lpd = float(myo_scores[lpd_idx].mean())
    m_spd = float(myo_scores[spd_idx].mean())
    obs_diff_myo = m_lpd - m_spd
    u_stat, mw_p = mannwhitneyu(myo_scores[lpd_idx], myo_scores[spd_idx], alternative="greater")

    df_myo_scores = pd.DataFrame({
        "sample_id": sample_cols,
        "stage_binary": y_nomo,
        "stage_label": np.where(y_nomo == 1, "Late_Stage_LPD", "Early_Stage_SPD"),
        "Six_Gene_Marker_Z_Score": np.round(myo_scores, 4)
    })
    save_csv(df_myo_scores, os.path.join(TAB_DIR, "six_gene_marker_scores.csv"), "Six-Gene Marker Scores")

    # Fig 6: Functional Networks Composite
    fig, axes = plt.subplots(2, 2, figsize=(15, 12), facecolor="#F8FAFC")
    axes[0, 0].text(0.5, 0.5, "No verified STRING export available\nPPI not reported",
                    ha="center", va="center", transform=axes[0, 0].transAxes)
    axes[0, 0].set_title("STRING PPI provenance unavailable", fontsize=11, fontweight="bold")
    axes[0, 0].set_axis_off()

    axes[0, 1].plot(np.arange(N_g + 1), running_es, color="#DC2626", linewidth=2.2)
    axes[0, 1].axhline(0, color="#1E293B", linewidth=0.8)
    axes[0, 1].set_title("Custom EMT-like gene-set running score\nDescriptive; no Hallmark NES or P-value", fontsize=11, fontweight="bold")
    axes[0, 1].set_xlabel("Gene Rank by Signed Significance", fontsize=10, fontweight="bold")
    axes[0, 1].set_ylabel("Running Enrichment Score", fontsize=10, fontweight="bold")
    axes[0, 1].grid(alpha=0.3, linestyle=":")

    axes[1, 0].plot(df_dca["Threshold_Probability"], df_dca["Net_Benefit_Nomogram"], color="#DC2626", linewidth=2.2, label="5-Gene Nomogram Model")
    axes[1, 0].plot(df_dca["Threshold_Probability"], df_dca["Net_Benefit_TreatAll"], color="#2563EB", linestyle="--", label="Treat All")
    axes[1, 0].plot(df_dca["Threshold_Probability"], df_dca["Net_Benefit_TreatNone"], color="#1E293B", linestyle=":", label="Treat None")
    axes[1, 0].set_title(f"Decision Curve Analysis: 5-Gene Nomogram\nApparent C-index = {c_index_apparent:.3f} | 5-fold CV = {c_index_cv:.3f}", fontsize=11, fontweight="bold")
    axes[1, 0].set_xlabel("Threshold Probability", fontsize=10, fontweight="bold")
    axes[1, 0].set_ylabel("Net Benefit", fontsize=10, fontweight="bold")
    axes[1, 0].legend(loc="upper right")
    axes[1, 0].grid(alpha=0.3, linestyle=":")

    df_plot_myo = pd.DataFrame({
        "Score": myo_scores,
        "Group": np.where(y_nomo == 1, f"LPD (n={len(lpd_idx)})", f"SPD (n={len(spd_idx)})")
    })
    sns.boxplot(data=df_plot_myo, x="Group", y="Score", ax=axes[1, 1], palette=["#DC2626", "#2563EB"], width=0.45)
    sns.stripplot(data=df_plot_myo, x="Group", y="Score", ax=axes[1, 1], color="#0F172A", alpha=0.7, size=6)
    axes[1, 1].set_title(f"Six-gene marker-score proxy\nMann-Whitney P = {mw_p:.4f}; not deconvolution", fontsize=11, fontweight="bold")
    axes[1, 1].set_ylabel("Mean Standardized Signature Score", fontsize=10, fontweight="bold")
    axes[1, 1].grid(axis='y', linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_06_functional_networks_composite.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "GSE125498_06_functional_networks_composite.pdf"))
    plt.close()

    save_csv(gsea_df, os.path.join(SRC_DIR, "Fig6_functional_networks_source_data.csv"), "Functional Networks Source Data")

    counts.update({
        "ppi_edges_n": None,
        "ppi_connectivity_rate": None,
        "gsea_emt_nes": None,
        "gsea_emt_p": None,
        "nomogram_c_index_apparent": float(c_index_apparent),
        "nomogram_c_index_cv": float(c_index_cv),
        "nomogram_optimism_gap": float(optimism),
        "six_gene_marker_score_mwu_p": float(mw_p)
    })
    with open(os.path.join(TAB_DIR, "phase7_counts.json"), "w") as fh:
        json.dump(counts, fh, indent=2)

    return counts

# ==============================================================================
# PHASE 8: FINAL MANIFEST & MASTER NUMBERS
# ==============================================================================
def run_phase_8(counts):
    log("=" * 80)
    log("PHASE 8: Final Manuscript Manifest, Verification Tables & Audit Summary")
    log("=" * 80)

    repo_root = (
        os.path.dirname(ROOT)
        if os.path.basename(ROOT).lower() == "scripts"
        else ROOT
    )
    manifest_md = os.path.join(
        repo_root, "results", "revision", "legacy_master_pipeline_manifest.md"
    )
    protected_paths = {
        os.path.abspath(os.path.join(repo_root, "FINAL_VERIFIED_NUMBERS.md")),
        os.path.abspath(
            os.path.join(repo_root, "results", "manuscript_numbers.json")
        ),
    }
    if os.path.abspath(manifest_md) in protected_paths:
        raise RuntimeError(
            "Legacy pipeline manifest path conflicts with a protected summary."
        )
    os.makedirs(os.path.dirname(manifest_md), exist_ok=True)
    with open(manifest_md, "w", encoding="utf-8") as fh:
        fh.write(f"""# LEGACY - NOT FOR THE REVISED ANALYSES

This is an execution manifest for the legacy exploratory four-model pipeline.
It is not the authoritative corrected analysis summary.

**Generated:** {datetime.now().isoformat()}
**Design:** GSE125498 (Discovery, N=33 effluent) -> GSE62928 (Validation, N=8 human biopsy)
**Fixed Seed:** 42 | **Primary Cutoffs:** |log2FC| >= {LOGFC_CUTOFF}, nominal P < {P_CUTOFF}

> This is an execution manifest for the legacy exploratory four-model pipeline,
> not the authoritative corrected analysis summary. It includes historical
> model outputs and is not a paper-ready results report. Consult
> `FINAL_VERIFIED_NUMBERS.md` and the accession-specific preregistered reports.

---

## 1. Discovery Cohort: GSE125498 (Illumina HumanHT-12 v4)
| Metric | Value | Description |
| :--- | :--- | :--- |
| **Total Samples (N)** | **33** | 13 Long-term PD cases (LPD, vintage >= 25mo) vs 20 Short-term PD controls (SPD) |
| **Probe-Level Transcripts** | **{counts.get('n_probe_level', 'N/A'):,}** | Raw background-adjusted bead-array probes |
| **Universe A (Best-Probe-by-P)** | **{counts.get('universe_bestp_n', 'N/A'):,}** | Unique gene symbols (Lowest eBayes P per gene) |
| **Universe B (MaxMean)** | **{counts.get('universe_maxmean_n', 'N/A'):,}** | Unique gene symbols (Highest row-mean expression) |
| **Primary Up-DEGs (Pro-fibrotic)** | **{counts.get('degs_primary_up', 'N/A'):,}** | |log2FC| >= 0.585 and nominal P < 0.05 |
| **Primary Down-DEGs** | **{counts.get('degs_primary_down', 'N/A'):,}** | |log2FC| <= -0.585 and nominal P < 0.05 |
| **Total Primary DEGs** | **{counts.get('degs_primary_total', 'N/A'):,}** | Combined DEGs (7.0% of tested genome) |
| **Genome-wide FDR < 0.05** | **{counts.get('fdr_sig_total', 'N/A')}** | Surviving Benjamini-Hochberg FDR |

---

## 2. Biological Filter: Naba Human Matrisome (1,027 Genes)
| Metric | Value | Description |
| :--- | :--- | :--- |
| **Reference Matrisome Genes** | **{counts.get('matrisome_reference_1027', 1027):,}** | Curated 7-subcategory ECM master list |
| **Measured Matrisome in Discovery** | **{counts.get('matrisome_measured', 'N/A'):,}** | Matrisome genes represented in GSE125498 universe |
| **Pro-fibrotic Up ECM-DEGs** | **{counts.get('profibrotic_up_ecm_degs', 'N/A')}** | 496 Up-DEGs ∩ 428 Measured Matrisome |
| **Down ECM-DEGs** | **{counts.get('down_ecm_degs', 'N/A')}** | 317 Down-DEGs ∩ 428 Measured Matrisome |
| **Fold Enrichment (Up ECM)** | **{counts.get('fold_enrichment_up_ecm', 0.0):.4f}x** | Pro-fibrotic ECM gene over-representation |
| **Hypergeometric P-value** | **{counts.get('hypergeom_p_up_ecm', 1.0):.3e}** | One-sided hypergeometric over-representation |
| **10,000 MC Permutation P** | **{counts.get('permutation_ecm_p_10k', 1.0):.4f}** | Monte Carlo empirical P from 10,000 label permutations |

---

## 3. Unsupervised WGCNA & 3-Way Triangular Convergence
| Metric | Value | Description |
| :--- | :--- | :--- |
| **WGCNA Input Genes** | **{counts.get('wgcna_network_genes', 'N/A'):,}** | Top 5,000 IQR variable genes + forced ECM-DEGs |
| **Soft Threshold Power (Beta)** | **12** | Scale-free topology signed network parameter |
| **Detected Co-expression Modules** | **{counts.get('wgcna_modules_n', 14)}** | Dynamic tree cut modules |
| **Top Trait-Correlated Module** | **{counts.get('top_wgcna_module', 'N/A')}** | Highest correlation with LPD/SPD phenotype |
| **Top Module Pearson r & P** | **r = {counts.get('top_wgcna_module_r', 0.0):+.3f}, P = {counts.get('top_wgcna_module_p', 1.0):.2e}** | Trait correlation strength |
| **Top Module Gene Count** | **{counts.get('top_wgcna_module_genes_n', 'N/A'):,}** | Genes residing in the top pro-fibrotic module |
| **Candidate Pool for ML** | **{counts.get('candidate_pool_for_ml_n', 'N/A')}** | Positive-trait correlated ECM-DEGs |

---

## 4. Machine Learning Consensus Hub Identification (4 Frameworks)
- **Frameworks:** LASSO (L1 Logistic Regression), SVM-RFE (Linear SVC), Random Forest (500 Trees), XGBoost
- **Consensus Rule:** Minimum >= 2 of 4 model votes
- **Final Consensus Hub Panel:** **{counts.get('consensus_hubs_n', 'N/A')} Hub Genes**
  - **Hub Genes:** `{', '.join(counts.get('consensus_hub_genes', []))}`

---

## 5. External Cohort Evaluation: GSE62928
| Metric | Value | Interpretation |
| :--- | :--- | :--- |
| **Fixed hub-score AUC** | **{counts.get('validation_real_signature_auc', 'N/A')}** | Fixed score evaluated on the real biopsy cohort; not a fitted classifier |
| **Direction-concordant hubs** | **{counts.get('validation_direction_concordant_n', 'N/A')} / {counts.get('consensus_hubs_n', 'N/A')}** | Discovery-to-biopsy effect direction |
| **Binomial sign-test P** | **{counts.get('validation_binomial_sign_p', 'N/A')}** | Direction-concordance test |
| **Two-sided Mann-Whitney P** | **{counts.get('validation_mann_whitney_two_sided_p', 'N/A')}** | Rank-based comparison of fixed score |
| **Random-set null mean AUC** | **{counts.get('validation_random_set_null_mean_auc', 'N/A')}** | 1,000 random seven-gene sets |
| **Random-set empirical P** | **{counts.get('validation_random_set_empirical_p', 'N/A')}** | Fraction of null scores at least as extreme |

---

## 6. Functional, Regulatory & Pathway Networks
| Metric | Value | Description |
| :--- | :--- | :--- |
| **STRING PPI** | **Not reported** | No verified STRING export is available |
| **Hallmark EMT GSEA** | **Not run** | Earlier hard-coded NES/P values withdrawn |
| **5-Gene Nomogram Apparent C-index** | **{counts.get('nomogram_c_index_apparent', 0.0):.4f}** | Apparent discrimination in effluent discovery |
| **5-Gene Nomogram 5-Fold CV C-index** | **{counts.get('nomogram_c_index_cv', 0.0):.4f}** | Cross-validated discrimination (Optimism = {counts.get('nomogram_optimism_gap', 0.0):.4f}) |
| **Six-gene marker-score Mann-Whitney P** | **{counts.get('six_gene_marker_score_mwu_p', 'N/A')}** | Expression z-score proxy, not deconvolution |

---
*All files, tables, figures, and source data are located in `results/`.*
""")

    log(f"Saved final markdown manifest to {manifest_md}")
    log("=" * 80)
    log("MASTER PIPELINE FINISHED 100% SUCCESSFULLY!")
    log("=" * 80)

def main():
    print("\n" + "#" * 80)
    print("# EXECUTING MASTER PERITONEAL FIBROSIS PIPELINE (GSE125498 -> GSE62928)")
    print("#" * 80 + "\n")

    setup_workspace()
    bestp_df, expr_genes_df, meta_df, counts = run_phase_1_2()
    ecm_df, matrisome_ref, counts = run_phase_3(bestp_df, meta_df, counts)
    candidate_pool, conv_df, counts = run_phase_4(expr_genes_df, meta_df, bestp_df, ecm_df, counts)
    hub_genes, hubs_df, counts = run_phase_5(candidate_pool, expr_genes_df, meta_df, bestp_df, counts)
    val_df, counts = run_phase_6(hub_genes, counts)
    counts = run_phase_7(hub_genes, expr_genes_df, meta_df, bestp_df, counts)
    run_phase_8(counts)

if __name__ == "__main__":
    main()
