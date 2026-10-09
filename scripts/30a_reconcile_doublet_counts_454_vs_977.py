"""
scripts/30a_reconcile_doublet_counts_454_vs_977.py
Names which run produced each doublet total:
  * 977 -> script 14 (14_scrna_e2_qc_and_filtering.py): NOT Scrublet. A count/gene percentile
           heuristic. Its output is the 94,517-cell object data/processed/GSE248762_annotated.h5ad
           and audit/stage3_audit_summary.json ("Total_Doublets_Removed": 977, counted AFTER QC).
  * 454 -> script 26 (26_stage3_repair_hub_blinding_and_harmony.py): scanpy Scrublet with the
           automatic per-sample threshold (expected rate 0.06), counted BEFORE QC. Wrote
           B2_scrublet_summary.csv and overwrote E2_filtering_summary.csv, then crashed
           (pd.read_csv(..., index=1) / R Harmony) before writing any h5ad.
This script re-executes script 14's heuristic exactly (same RNG seed and call order) to get the
per-sample 977 breakdown, and checks it against the per-donor cell counts stored in the old h5ad.
Only per-cell totals are used (no gene-level values), so no hub gene is read.
"""
import os, sys, gzip
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad

sys.path.append("scripts")
from guard import guard_check
guard_check([], stage="30a_start")
pd.set_option("display.width", 250)

extract_dir = "data/raw/GSE248762_extracted"
files = os.listdir(extract_dir)
sd = {}
for fname in files:
    gsm = fname.split("_", 1)[0]
    sd.setdefault(gsm, {})
    if "barcode" in fname.lower():
        sd[gsm]["barcodes"] = os.path.join(extract_dir, fname)
    elif "feature" in fname.lower() or "gene" in fname.lower():
        sd[gsm]["features"] = os.path.join(extract_dir, fname)
    elif "matrix" in fname.lower():
        sd[gsm]["matrix"] = os.path.join(extract_dir, fname)
        sd[gsm]["matrix_fname"] = fname


def calc_mad(x):
    return np.median(np.abs(x - np.median(x)))


np.random.seed(42)
rows = []
for gsm, p in sorted(sd.items()):
    title = p["matrix_fname"].split(".")[0].split("_", 1)[1]
    mat = scipy.io.mmread(p["matrix"]).T.tocsr()
    feats = pd.read_csv(p["features"], sep="\t", header=None, compression="gzip")
    mito = feats[1].astype(str).str.upper().str.startswith("MT-").values
    n = mat.shape[0]
    nc = np.asarray(mat.sum(axis=1)).ravel()
    ng = np.asarray((mat > 0).sum(axis=1)).ravel()
    pm = np.asarray(mat[:, mito].sum(axis=1)).ravel() / (nc + 1e-9) * 100
    lc, lg = np.log10(nc + 1), np.log10(ng + 1)
    c_lo = max(500.0, 10 ** (np.median(lc) - 3 * calc_mad(lc))); c_hi = 10 ** (np.median(lc) + 3 * calc_mad(lc))
    g_lo = max(200.0, 10 ** (np.median(lg) - 3 * calc_mad(lg))); g_hi = 10 ** (np.median(lg) + 3 * calc_mad(lg))
    pc = (nc >= c_lo) & (nc <= c_hi); pg = (ng >= g_lo) & (ng <= g_hi); pmt = pm <= 15.0
    n_sim = int(0.5 * n)
    i1 = np.random.choice(n, n_sim); i2 = np.random.choice(n, n_sim)
    s = mat[i1] + mat[i2]
    sim_c = np.asarray(s.sum(axis=1)).ravel(); sim_g = np.asarray((s > 0).sum(axis=1)).ravel()
    is_d = (nc > np.percentile(sim_c, 80)) & (ng > np.percentile(sim_g, 80)) & (nc > np.percentile(nc, 95))
    keep = pc & pg & pmt & ~is_d
    # replicate script 14's violin subsampling RNG calls so later samples see the same RNG state
    np.random.choice(n, min(500, n), replace=False)
    if keep.sum() > 0:
        np.random.choice(int(keep.sum()), min(500, int(keep.sum())), replace=False)
    rows.append({"GSM": gsm, "Sample": title, "Barcodes": n,
                 "Heuristic_flagged_all": int(is_d.sum()),
                 "Heuristic_doublets_counted_after_QC": int((is_d & pc & pg & pmt).sum()),
                 "Post_QC_cells": int(keep.sum())})
    print(f"  {gsm} {title:12s} flagged={is_d.sum():4d} counted(after QC)={(is_d & pc & pg & pmt).sum():4d} kept={keep.sum():5d}")

df = pd.DataFrame(rows)
old = ad.read_h5ad("data/processed/GSE248762_annotated.h5ad", backed="r")
old_counts = old.obs["donor_id"].value_counts()
df["Cells_in_old_annotated_h5ad"] = df["Sample"].map(old_counts).fillna(0).astype(int)
df["Match"] = df["Post_QC_cells"] == df["Cells_in_old_annotated_h5ad"]
b2 = pd.read_csv("results/tables/B2_scrublet_summary.csv").set_index("Sample_Title")["Predicted_Doublets"]
df["Script26_scrublet_auto_preQC"] = df["Sample"].map(b2)
print("\n--- RECONCILIATION: 977 (script 14 heuristic) vs 454 (script 26 Scrublet auto-threshold) ---")
print(df.to_string(index=False))
print(f"\nScript 14 heuristic, counted after QC (stage3_audit_summary 'Total_Doublets_Removed'): {df['Heuristic_doublets_counted_after_QC'].sum()}")
print(f"Script 14 heuristic, all flagged barcodes: {df['Heuristic_flagged_all'].sum()}")
print(f"Script 14 post-QC cells: {df['Post_QC_cells'].sum()} | old annotated h5ad n_obs: {old.n_obs} | per-donor all match: {bool(df['Match'].all())}")
print(f"Script 26 Scrublet auto-threshold, pre-QC (B2 table): {int(df['Script26_scrublet_auto_preQC'].sum())}")
print(f"Old annotated h5ad obs columns: {list(old.obs.columns)}")
df.to_csv("results/tables/E2_doublet_reconciliation_454_vs_977.csv", index=False)
old.file.close()
