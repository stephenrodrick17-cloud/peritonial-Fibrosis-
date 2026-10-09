"""
scripts/29_scrna_doublets_qc_rerun.py
Stage 3 repair (round 2), part A: doublets + QC on GSE248762, hub-gene-blind.

ONE doublet rule for all 16 samples (pre-specified before any call was made):
  * Scrublet (scanpy.pp.scrublet), run per sample on all deposited barcodes, BEFORE QC.
  * Fixed parameters for every sample: expected_doublet_rate=0.06, sim_doublet_ratio=2.0,
    n_prin_comps=30, log_transform=False (scanpy default), random_state=0.
  * ONE score cutoff for all samples, set from the simulated-doublet histogram:
    the scores of simulated doublets from all 16 samples are pooled, a Gaussian KDE is fit,
    and the cutoff is the density minimum (valley) between the low "embedded" mode and the
    high "neotypic" mode. Every sample uses that same cutoff.
  * Each sample's own valley is computed and printed for transparency only (not used).

QC (unchanged from the frozen plan / earlier scripts): per sample, 3 MAD on log10 UMI and
log10 genes with hard floors (UMI >= 500, genes >= 200), mito ceiling 15%.

Outputs
  data/processed/GSE248762_hubblind_allcells_qc.h5ad   (ALL raw barcodes, hub genes removed,
                                                        obs carries QC + doublet flags)
  data/processed/scrublet_scores/<GSM>.npz             (observed + simulated scores)
  results/tables/E2_doublets_scrublet_fixedrule.csv
  results/tables/E2_filtering_summary.csv              (regenerated FROM the h5ad, see script 29b)
  results/figures/E2_doublet_hist_observed_vs_simulated.png
  results/figures/E2_doublet_pooled_simulated_kde.png
  results/figures/E2_qc_violins_per_sample.png
"""
import os
import sys
import gzip
import json
import time
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("scripts")
from guard import guard_check, BLOCKED_HUB_GENES

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
pd.set_option("display.max_rows", 500)

RAW_DIR = "data/raw/GSE248762_extracted"
CACHE = "data/processed/scrublet_scores"
OUT_H5AD = "data/processed/GSE248762_hubblind_allcells_qc.h5ad"
os.makedirs(CACHE, exist_ok=True)
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

EXPECTED_RATE = 0.06
SIM_RATIO = 2.0
N_PCS_SCRUB = 30
SEED = 0
UMI_FLOOR, GENE_FLOOR, MITO_CEIL = 500.0, 200.0, 15.0

print("=" * 78)
print("SCRIPT 29: DOUBLETS (one Scrublet rule) + PER-SAMPLE QC, hub-gene-blind")
print("=" * 78)
guard_check([], stage="29_start")

with open("provenance/GSE248762_acquisition_summary.json") as f:
    acq = json.load(f)
samples = sorted(acq["samples"], key=lambda s: s["GSM"])


def mad(x):
    return np.median(np.abs(x - np.median(x)))


adatas = []
SIMS = {}
raw_min_rows = []
for s in samples:
    gsm, title = s["GSM"], s["Sample_Title"]
    grp = s["Group"].split(" ")[0]
    files = s["Files"]
    with gzip.open(os.path.join(RAW_DIR, files["barcodes"]["filename"]), "rt") as f:
        barcodes = [l.strip() for l in f if l.strip()]
    feats = pd.read_csv(os.path.join(RAW_DIR, files["features"]["filename"]), sep="\t",
                        header=None, compression="gzip")
    mat = scipy.io.mmread(os.path.join(RAW_DIR, files["matrix"]["filename"])).T.tocsr().astype(np.float32)

    # deposited-matrix floor check on the FULL gene set (before hub removal); only totals are used
    full_umi = np.asarray(mat.sum(axis=1)).ravel()
    full_genes = np.asarray((mat > 0).sum(axis=1)).ravel()

    keep_g = ~feats[1].str.upper().isin(BLOCKED_HUB_GENES).values
    mat = mat[:, keep_g]
    var = pd.DataFrame({"gene_ids": feats[0].values[keep_g]}, index=feats[1].values[keep_g])
    a = ad.AnnData(X=mat, obs=pd.DataFrame(index=[f"{gsm}_{b}" for b in barcodes]), var=var)
    a.var_names_make_unique()
    guard_check(a.var_names.tolist(), stage=f"29_load_{gsm}")
    a.obs["gsm"], a.obs["donor_id"], a.obs["group"] = gsm, title, grp

    n_counts = np.asarray(a.X.sum(axis=1)).ravel()
    n_genes = np.asarray((a.X > 0).sum(axis=1)).ravel()
    mt = a.var_names.str.upper().str.startswith("MT-")
    pct_mt = np.asarray(a[:, mt].X.sum(axis=1)).ravel() / np.maximum(n_counts, 1) * 100
    a.obs["n_counts"], a.obs["n_genes"], a.obs["pct_counts_mt"] = n_counts, n_genes, pct_mt

    raw_min_rows.append({"Sample": title, "Barcodes": a.n_obs,
                         "Min_UMI_all_genes": int(full_umi.min()),
                         "Min_genes_all_genes": int(full_genes.min()),
                         "Barcodes_UMI_lt_500": int((full_umi < 500).sum()),
                         "Barcodes_UMI_lt_1000": int((full_umi < 1000).sum()),
                         "P1_UMI": round(float(np.percentile(full_umi, 1)), 1)})

    # ---- Scrublet (cached) ----
    cpath = os.path.join(CACHE, f"{gsm}.npz")
    if os.path.exists(cpath):
        z = np.load(cpath)
        obs_s, sim_s = z["obs"], z["sim"]
        auto_thr = float(z["auto_thr"])
    else:
        tmp = ad.AnnData(X=a.X.copy(), obs=pd.DataFrame(index=a.obs_names), var=pd.DataFrame(index=a.var_names))
        sc.pp.scrublet(tmp, expected_doublet_rate=EXPECTED_RATE, sim_doublet_ratio=SIM_RATIO,
                       n_prin_comps=N_PCS_SCRUB, random_state=SEED, verbose=False)
        obs_s = tmp.obs["doublet_score"].values.astype(float)
        sim_s = np.asarray(tmp.uns["scrublet"]["doublet_scores_sim"], dtype=float)
        auto_thr = float(tmp.uns["scrublet"].get("threshold", np.nan))
        np.savez_compressed(cpath, obs=obs_s, sim=sim_s, auto_thr=auto_thr)
        del tmp
    a.obs["doublet_score"] = obs_s
    SIMS[gsm] = sim_s
    a.obs["scanpy_auto_threshold"] = auto_thr
    print(f"  {gsm} {title:12s} {grp:9s} barcodes={a.n_obs:6d}  scrublet obs/sim={len(obs_s)}/{len(sim_s)}  scanpy-auto-thr={auto_thr:.4f}")
    adatas.append(a)

# ---------------------------------------------------------------------------
# ONE cutoff from the pooled simulated-doublet histogram
# ---------------------------------------------------------------------------
GRID = np.linspace(0.0, 1.0, 1001)


def kde_valley(sim):
    dens = gaussian_kde(sim)(GRID)
    maxima = argrelextrema(dens, np.greater, order=10)[0]
    minima = argrelextrema(dens, np.less, order=10)[0]
    if len(maxima) < 2:
        return np.nan, dens
    # low mode = highest peak; high mode = highest peak to its right
    p1 = maxima[np.argmax(dens[maxima])]
    right = maxima[maxima > p1]
    if len(right) == 0:
        return np.nan, dens
    p2 = right[np.argmax(dens[right])]
    mins = minima[(minima > p1) & (minima < p2)]
    if len(mins) == 0:
        return np.nan, dens
    v = mins[np.argmin(dens[mins])]
    return float(GRID[v]), dens


pooled_sim = np.concatenate([SIMS[a.obs["gsm"].iloc[0]] for a in adatas])
rng = np.random.RandomState(SEED)
kde_sample = pooled_sim if len(pooled_sim) <= 60000 else rng.choice(pooled_sim, 60000, replace=False)
CUTOFF, pooled_dens = kde_valley(kde_sample)
print(f"\nPooled simulated doublets: n={len(pooled_sim)} (KDE fit on {len(kde_sample)})")
print(f"ONE GLOBAL SCRUBLET CUTOFF (valley of pooled simulated-doublet KDE) = {CUTOFF:.4f}")
if not np.isfinite(CUTOFF):
    raise RuntimeError("Pooled simulated-doublet distribution is not bimodal; rule undefined. Stop.")

# ---------------------------------------------------------------------------
# Per-sample calls, QC flags
# ---------------------------------------------------------------------------
dbl_rows = []
for a in adatas:
    title = a.obs["donor_id"].iloc[0]
    own_valley, _ = kde_valley(SIMS[a.obs["gsm"].iloc[0]])
    pred = a.obs["doublet_score"].values >= CUTOFF
    a.obs["predicted_doublet"] = pred
    a.obs["doublet_cutoff_used"] = CUTOFF
    sim_above = float((SIMS[a.obs["gsm"].iloc[0]] >= CUTOFF).mean())
    dbl_rows.append({
        "GSM": a.obs["gsm"].iloc[0], "Sample": title, "Group": a.obs["group"].iloc[0],
        "Barcodes": a.n_obs, "Expected_rate": EXPECTED_RATE,
        "Threshold_used": round(CUTOFF, 4),
        "Doublets": int(pred.sum()), "Doublet_rate_pct": round(pred.mean() * 100, 2),
        "Sim_doublets_above_cut_pct": round(sim_above * 100, 1),
        "Est_detectable_doublet_rate_pct": round(pred.mean() / max(sim_above, 1e-9) * 100, 2),
        "Own_sim_valley_(diagnostic)": round(own_valley, 4) if np.isfinite(own_valley) else np.nan,
        "Scanpy_auto_threshold_(old run)": round(float(a.obs["scanpy_auto_threshold"].iloc[0]), 4),
    })

    log_c, log_g = np.log10(a.obs["n_counts"].values + 1), np.log10(a.obs["n_genes"].values + 1)
    c_lo = max(UMI_FLOOR, 10 ** (np.median(log_c) - 3 * mad(log_c)))
    c_hi = 10 ** (np.median(log_c) + 3 * mad(log_c))
    g_lo = max(GENE_FLOOR, 10 ** (np.median(log_g) - 3 * mad(log_g)))
    g_hi = 10 ** (np.median(log_g) + 3 * mad(log_g))
    nc, ng, pm = a.obs["n_counts"].values, a.obs["n_genes"].values, a.obs["pct_counts_mt"].values
    a.obs["fail_umi_low"] = nc < c_lo
    a.obs["fail_umi_high"] = nc > c_hi
    a.obs["fail_genes_low"] = ng < g_lo
    a.obs["fail_genes_high"] = ng > g_hi
    a.obs["fail_mito"] = pm > MITO_CEIL
    a.obs["pass_qc"] = ~(a.obs[["fail_umi_low", "fail_umi_high", "fail_genes_low", "fail_genes_high", "fail_mito"]].any(axis=1))
    a.obs["keep"] = a.obs["pass_qc"] & ~a.obs["predicted_doublet"]
    a.obs["umi_lo"], a.obs["umi_hi"], a.obs["genes_lo"], a.obs["genes_hi"] = c_lo, c_hi, g_lo, g_hi
    a.obs["umi_lo_mad_unfloored"] = 10 ** (np.median(log_c) - 3 * mad(log_c))
    a.obs["genes_lo_mad_unfloored"] = 10 ** (np.median(log_g) - 3 * mad(log_g))

df_dbl = pd.DataFrame(dbl_rows)
df_dbl.to_csv("results/tables/E2_doublets_scrublet_fixedrule.csv", index=False)
print("\n--- DOUBLETS: ONE RULE (Scrublet, expected rate 0.06, single cutoff from pooled simulated histogram) ---")
print(df_dbl.to_string(index=False))
print(f"TOTAL doublets called = {df_dbl['Doublets'].sum()} / {df_dbl['Barcodes'].sum()} barcodes "
      f"({df_dbl['Doublets'].sum() / df_dbl['Barcodes'].sum() * 100:.2f}%)")

print("\n--- DEPOSITED-MATRIX FLOOR CHECK (all genes, before any QC) ---")
print(pd.DataFrame(raw_min_rows).to_string(index=False))

# ---------------------------------------------------------------------------
# Figures: observed AND simulated histograms
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(4, 4, figsize=(22, 17))
bins = np.linspace(0, 1, 81)
for ax, a, row in zip(axes.ravel(), adatas, dbl_rows):
    ax.hist(SIMS[a.obs["gsm"].iloc[0]], bins=bins, density=True, color="#f59e0b", alpha=0.55, label="Simulated doublets")
    ax.hist(a.obs["doublet_score"].values, bins=bins, density=True, color="#2563eb", alpha=0.55, label="Observed cells")
    ax.axvline(CUTOFF, color="#dc2626", ls="--", lw=1.6, label=f"Global cutoff {CUTOFF:.3f}")
    if np.isfinite(row["Own_sim_valley_(diagnostic)"]):
        ax.axvline(row["Own_sim_valley_(diagnostic)"], color="#6b7280", ls=":", lw=1.2, label="Own sim valley (diag.)")
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.set_title(f"{row['Sample']}  thr={CUTOFF:.3f}  n={row['Doublets']} ({row['Doublet_rate_pct']:.2f}%)",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel("Scrublet doublet score")
    ax.set_ylabel("density (log)")
axes.ravel()[0].legend(fontsize=8, loc="upper right")
fig.suptitle("Scrublet doublet scores per sample: observed cells (blue) vs simulated doublets (orange)\n"
             f"One rule: expected rate {EXPECTED_RATE}, single cutoff {CUTOFF:.3f} = valley of pooled simulated-doublet KDE",
             fontsize=15, fontweight="bold", y=1.0)
plt.tight_layout()
plt.savefig("results/figures/E2_doublet_hist_observed_vs_simulated.png", dpi=150, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(1, 2, figsize=(15, 5))
ax[0].hist(pooled_sim, bins=bins, density=True, color="#f59e0b", alpha=0.6, label="Pooled simulated doublets")
ax[0].plot(GRID, pooled_dens, color="k", lw=1.5, label="Gaussian KDE")
ax[0].axvline(CUTOFF, color="#dc2626", ls="--", lw=1.8, label=f"Valley = {CUTOFF:.3f}")
ax[0].set_title("Pooled simulated-doublet scores (16 samples)", fontweight="bold")
ax[0].set_xlabel("Scrublet doublet score"); ax[0].legend()
pooled_obs = np.concatenate([a.obs["doublet_score"].values for a in adatas])
ax[1].hist(pooled_obs, bins=bins, color="#2563eb", alpha=0.7, label="Pooled observed cells")
ax[1].axvline(CUTOFF, color="#dc2626", ls="--", lw=1.8)
ax[1].set_yscale("log"); ax[1].set_title("Pooled observed-cell scores", fontweight="bold")
ax[1].set_xlabel("Scrublet doublet score"); ax[1].legend()
plt.tight_layout()
plt.savefig("results/figures/E2_doublet_pooled_simulated_kde.png", dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# Concatenate and save ALL barcodes with flags
# ---------------------------------------------------------------------------
for a in adatas:
    pass
adata = ad.concat(adatas, join="inner", merge="same")
guard_check(adata.var_names.tolist(), stage="29_concat")
adata.uns["doublet_rule"] = {
    "method": "scanpy.pp.scrublet per sample, pre-QC", "expected_doublet_rate": EXPECTED_RATE,
    "sim_doublet_ratio": SIM_RATIO, "n_prin_comps": N_PCS_SCRUB, "random_state": SEED,
    "cutoff": CUTOFF, "cutoff_rule": "valley of Gaussian KDE of pooled simulated-doublet scores (all 16 samples)",
    "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
for c in ["fail_umi_low", "fail_umi_high", "fail_genes_low", "fail_genes_high", "fail_mito", "pass_qc", "keep", "predicted_doublet"]:
    adata.obs[c] = adata.obs[c].astype(bool)
adata.write_h5ad(OUT_H5AD)
print(f"\nSaved {OUT_H5AD}: {adata.n_obs} barcodes x {adata.n_vars} genes (hub genes absent: "
      f"{len(set(adata.var_names.str.upper()) & BLOCKED_HUB_GENES) == 0})")

# ---------------------------------------------------------------------------
# Per-sample QC violins (before vs after), thresholds drawn per sample
# ---------------------------------------------------------------------------
obs = adata.obs.copy()
order = [s["Sample_Title"] for s in samples]
before = obs.assign(Status="Deposited barcodes (>= ~500 UMI)")
after = obs[obs["keep"]].assign(Status="After QC + doublet removal")
dfp = pd.concat([before, after])
dfp["donor_id"] = pd.Categorical(dfp["donor_id"], categories=order, ordered=True)
thr = obs.groupby("donor_id", observed=True)[["umi_lo", "umi_hi", "genes_lo", "genes_hi"]].first().reindex(order)

fig, axes = plt.subplots(3, 1, figsize=(24, 17))
pal = {"Deposited barcodes (>= ~500 UMI)": "#94a3b8", "After QC + doublet removal": "#2563eb"}
for ax, col, ylab, logy in [(axes[0], "n_counts", "UMI per cell (log10)", True),
                            (axes[1], "n_genes", "Genes per cell (log10)", True),
                            (axes[2], "pct_counts_mt", "Mitochondrial %", False)]:
    sns.violinplot(data=dfp, x="donor_id", y=col, hue="Status", split=True, inner="quartile",
                   cut=0, density_norm="width", palette=pal, ax=ax, linewidth=0.8)
    if logy:
        ax.set_yscale("log")
    for i, d in enumerate(order):
        if col == "n_counts":
            ax.hlines([thr.loc[d, "umi_lo"], thr.loc[d, "umi_hi"]], i - 0.45, i + 0.45, colors="#dc2626", lw=1.4)
        elif col == "n_genes":
            ax.hlines([thr.loc[d, "genes_lo"], thr.loc[d, "genes_hi"]], i - 0.45, i + 0.45, colors="#dc2626", lw=1.4)
    if col == "pct_counts_mt":
        ax.axhline(MITO_CEIL, color="#dc2626", lw=1.4, ls="--")
    ax.set_ylabel(ylab, fontsize=13)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=12, rotation=30)
    ax.legend(loc="upper right", fontsize=10)
axes[0].set_title("Per-sample QC (not pooled). Red bars = that sample's 3-MAD limits (UMI floor 500, gene floor 200); dashed = 15% mito ceiling",
                  fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("results/figures/E2_qc_violins_per_sample.png", dpi=130)
plt.close()
print("Saved figures: E2_doublet_hist_observed_vs_simulated.png, E2_doublet_pooled_simulated_kde.png, E2_qc_violins_per_sample.png")
print("Script 29 complete.")
