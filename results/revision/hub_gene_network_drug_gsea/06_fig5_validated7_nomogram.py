"""
Figure 5: Nomogram construction and association between immune cell infiltration and validated hub genes.
Replicating the layout of Sun et al. (Electron J Biotechnol 2026, Fig. 5) using the seven validated hub genes.

Panels:
  A: Ridge logistic regression nomogram (TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, CST7).
  B: Calibration curve - apparent + bootstrap-corrected, five equal-count bins, exact binomial CIs.
  C: Decision curve analysis (DCA).
  D: Clinical impact curve (counts scaled to 1000 from n = 33).
  E: ssGSEA immune cell enrichment differences (LPD vs SPD), Mann-Whitney U, BH-FDR stars.
  F: Immune-immune Spearman correlation heatmap.
  G: Hub gene-immune Spearman correlation heatmap.

Method notes recorded here for the figure legend:
  * Nomogram gene axes are cohort-standardised z-scores; the grey text beneath each gene is the
    raw (log2) discovery-cohort mean +/- SD, provided as a read-out aid for the z-scale.
  * Panel B uses FIVE equal-count bins (n = 7, 6, 7, 6, 7). Bins 1 and 2 both contain zero events,
    which is why only four visually distinct points appear.
  * Panel E stars are BH-adjusted across the 13 cell types.
  * Panel G stars are BH-adjusted across all 91 gene-cell-type tests.
  * Discrimination: apparent AUC 0.969; leave-one-out CV AUC 0.892; 5-fold CV AUC 0.919.
"""

from __future__ import annotations

import os
import textwrap
import warnings

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import beta, mannwhitneyu, spearmanr, ttest_ind
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

matplotlib.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 0.9})

# ----------------------------------------------------------------------------- paths
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
EXPR_PATH = os.path.join(ROOT, "results", "tables", "GSE125498_full_expression_matrix.csv")
META_PATH = os.path.join(ROOT, "results", "tables", "GSE125498_sample_metadata.csv")
SSGSEA_PATH = os.path.join(ROOT, "results", "paper_style", "Fig5E_saved_ssgsea_sample_scores.csv")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

GENES = ["TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "ECM1", "CST7"]
CELL_TYPES = ["CD4 T cell", "CD8 T cell", "NKT cell", "NK cell", "T cell (other)", "B cell",
              "Plasma cell", "Monocyte", "Macrophage", "Dendritic cell", "Neutrophil",
              "Eosinophil", "Mast cell"]
N_BINS = 5
N_BOOT = 1000
SEED = 2026

C_LPD = "#F43F5E"
C_SPD = "#0284C7"
C_MODEL = "#C2185B"
C_ALL = "#94A3B8"


def stars(p: float) -> str:
    """Return the conventional asterisk code for a P value or BH-adjusted q value."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


# ----------------------------------------------------------------------------- data
expr = pd.read_csv(EXPR_PATH)
meta = pd.read_csv(META_PATH)
ssgsea = pd.read_csv(SSGSEA_PATH)

samples = meta["sample_id"].tolist()
y = meta["stage_binary"].values
X_raw = expr.set_index("Gene")[samples].loc[GENES].T
X_z = pd.DataFrame(StandardScaler().fit_transform(X_raw), columns=GENES, index=samples)

model = LogisticRegression(C=1.0, solver="liblinear", random_state=SEED).fit(X_z, y)
probs = model.predict_proba(X_z)[:, 1]
coefs = dict(zip(GENES, model.coef_[0]))
intercept = float(model.intercept_[0])

# ----------------------------------------------------------------------------- A: nomogram
z_min = {g: float(X_z[g].min()) for g in GENES}
z_max = {g: float(X_z[g].max()) for g in GENES}
lp_span = {g: coefs[g] * (z_max[g] - z_min[g]) for g in GENES}
max_span = max(lp_span.values())
points_per_lp = 100.0 / max_span
gene_points = {g: lp_span[g] * points_per_lp for g in GENES}
total_points_max = float(sum(gene_points.values()))
base_lp = intercept + sum(coefs[g] * z_min[g] for g in GENES)

nomogram_source = pd.DataFrame(
    {
        "gene": GENES,
        "coefficient_z": [coefs[g] for g in GENES],
        "z_min": [z_min[g] for g in GENES],
        "z_max": [z_max[g] for g in GENES],
        "points_at_z_max": [gene_points[g] for g in GENES],
        "raw_mean_log2": [X_raw[g].mean() for g in GENES],
        "raw_sd_log2": [X_raw[g].std() for g in GENES],
    }
).sort_values("points_at_z_max", ascending=False)
nomogram_source.to_csv(os.path.join(OUT_DIR, "Fig5A_nomogram_points_source.csv"), index=False)

# ----------------------------------------------------------------------------- B: calibration
cal = pd.DataFrame({"prob": probs, "y": y}).sort_values("prob").reset_index(drop=True)
cal["bin"] = pd.qcut(cal["prob"], N_BINS, labels=False, duplicates="drop")
bin_counts = cal.groupby("bin").size().tolist()


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact binomial (Clopper-Pearson) confidence interval."""
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


cal_summary = cal.groupby("bin").agg(mean_pred=("prob", "mean"), obs=("y", "mean"), n=("y", "count"))
cal_summary["k"] = (cal_summary["obs"] * cal_summary["n"]).round().astype(int)
ci = [clopper_pearson(r.k, r.n) for r in cal_summary.itertuples()]
cal_summary["ci_low"] = [c[0] for c in ci]
cal_summary["ci_high"] = [c[1] for c in ci]
cal_summary["bin_n"] = bin_counts

# Bootstrap-corrected calibration: for each resample, refit, bin the ORIGINAL cohort by the
# bootstrap model's predicted probabilities, and record mean predicted vs observed per bin.
rng = np.random.default_rng(SEED)
boot_curves = []
for _ in range(N_BOOT):
    idx = rng.integers(0, len(y), len(y))
    if len(np.unique(y[idx])) < 2:
        continue
    mb = LogisticRegression(C=1.0, solver="liblinear").fit(X_z.iloc[idx], y[idx])
    pb = mb.predict_proba(X_z)[:, 1]
    tmp = pd.DataFrame({"prob": pb, "y": y})
    tmp["b"] = pd.qcut(tmp["prob"], N_BINS, labels=False, duplicates="drop")
    g = tmp.groupby("b").agg(mp=("prob", "mean"), ob=("y", "mean"))
    if len(g) == N_BINS:
        boot_curves.append(g[["mp", "ob"]].to_numpy())
boot_arr = np.array(boot_curves)
cal_summary["boot_mean_pred"] = boot_arr[:, :, 0].mean(axis=0)
cal_summary["boot_obs"] = boot_arr[:, :, 1].mean(axis=0)
cal_summary.to_csv(os.path.join(OUT_DIR, "Fig5B_calibration_source.csv"))

# Discrimination with optimism control
apparent_auc = roc_auc_score(y, probs)
loo_auc = roc_auc_score(
    y, cross_val_predict(LogisticRegression(C=1.0, solver="liblinear"), X_z, y,
                         cv=LeaveOneOut(), method="predict_proba")[:, 1])
cv_auc = roc_auc_score(
    y, cross_val_predict(LogisticRegression(C=1.0, solver="liblinear"), X_z, y,
                         cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                         method="predict_proba")[:, 1])

# ----------------------------------------------------------------------------- C/D: DCA
thresholds = np.linspace(0.01, 0.99, 99)
n = len(y)
nb_model, nb_all = [], []
for t in thresholds:
    tp = int(((probs >= t) & (y == 1)).sum())
    fp = int(((probs >= t) & (y == 0)).sum())
    nb_model.append(tp / n - (fp / n) * (t / (1 - t)))
    nb_all.append(int((y == 1).sum()) / n - int((y == 0).sum()) / n * (t / (1 - t)))
nb_model, nb_all = np.array(nb_model), np.array(nb_all)
pd.DataFrame({"Threshold": thresholds, "Net_benefit_nomogram": nb_model,
              "Net_benefit_treat_all": nb_all,
              "Net_benefit_treat_none": 0.0}).to_csv(
    os.path.join(OUT_DIR, "Fig5C_dca_source.csv"), index=False)

impact_n = np.array([int((probs >= t).sum()) for t in thresholds]) * (1000.0 / n)
impact_true = np.array([int(((probs >= t) & (y == 1)).sum()) for t in thresholds]) * (1000.0 / n)
pd.DataFrame({"Threshold": thresholds, "Number_high_risk_per_1000": impact_n,
              "Number_high_risk_with_event_per_1000": impact_true}).to_csv(
    os.path.join(OUT_DIR, "Fig5D_clinical_impact_source.csv"), index=False)

# ----------------------------------------------------------------------------- E: ssGSEA stats
rows = []
for ct in CELL_TYPES:
    lpd = ssgsea.loc[ssgsea["group"] == "LPD", ct].values
    spd = ssgsea.loc[ssgsea["group"] == "SPD", ct].values
    _, p_raw = mannwhitneyu(lpd, spd, alternative="two-sided")
    rows.append({"cell_type": ct, "median_LPD": np.median(lpd), "median_SPD": np.median(spd),
                 "delta_median": np.median(lpd) - np.median(spd), "P_raw": p_raw})
ssgsea_stats = pd.DataFrame(rows)
ssgsea_stats["P_FDR_BH"] = multipletests(ssgsea_stats["P_raw"], method="fdr_bh")[1]
ssgsea_stats["stars_raw"] = ssgsea_stats["P_raw"].map(stars)
ssgsea_stats["stars_fdr"] = ssgsea_stats["P_FDR_BH"].map(stars)
ssgsea_stats["survives_FDR"] = ssgsea_stats["P_FDR_BH"] < 0.05
ssgsea_stats.to_csv(os.path.join(OUT_DIR, "Fig5E_ssgsea_stats_source.csv"), index=False)

# ----------------------------------------------------------------------------- F: immune-immune
imm_cor = ssgsea[CELL_TYPES].corr(method="spearman")
imm_long = imm_cor.stack().rename("Spearman_rho").reset_index()
imm_long.columns = ["Cell_type_1", "Cell_type_2", "Spearman_rho"]
imm_long.to_csv(os.path.join(OUT_DIR, "Fig5F_immune_immune_corr_source.csv"), index=False)

# ----------------------------------------------------------------------------- G: gene-immune
gene_rho = pd.DataFrame(index=GENES, columns=CELL_TYPES, dtype=float)
gene_p = pd.DataFrame(index=GENES, columns=CELL_TYPES, dtype=float)
for g in GENES:
    for ct in CELL_TYPES:
        rho, p = spearmanr(X_raw[g], ssgsea.set_index("sample_id").loc[samples, ct])
        gene_rho.loc[g, ct], gene_p.loc[g, ct] = rho, p
gene_rho.to_csv(os.path.join(OUT_DIR, "Fig5G_gene_immune_rho_source.csv"))
gene_p.to_csv(os.path.join(OUT_DIR, "Fig5G_gene_immune_p_source.csv"))
gene_fdr = pd.DataFrame(
    multipletests(gene_p.to_numpy().ravel(), method="fdr_bh")[1].reshape(gene_p.shape),
    index=GENES,
    columns=CELL_TYPES,
)
gene_fdr.to_csv(os.path.join(OUT_DIR, "Fig5G_gene_immune_fdr_source.csv"))

# ============================================================================= FIGURE
fig = plt.figure(figsize=(15.5, 19.5), facecolor="white")
gs = gridspec.GridSpec(4, 2, height_ratios=[1.25, 1.0, 1.05, 1.2], hspace=0.50, wspace=0.26)

# ------------------------------------------------------------------- Panel A: nomogram
ax_a = fig.add_subplot(gs[0, 0])
ax_a.axis("off")
ax_a.set_xlim(-42, 112)
ax_a.set_ylim(-2.1, 11.2)
ax_a.set_title("Constructed Nomogram (ridge logistic regression)", fontsize=12,
               fontweight="bold", pad=14)

row_pitch = 1.25
row_y = {"points": 10.1}
row_y.update({g: 10.1 - row_pitch * (i + 1) for i, g in enumerate(GENES)})
row_y["total"] = row_y[GENES[-1]] - row_pitch - 0.45
row_y["risk"] = row_y["total"] - row_pitch

# Points ruler
ax_a.text(-6, row_y["points"], "Points", ha="right", va="center", fontsize=9.5, fontweight="bold")
ax_a.plot([0, 100], [row_y["points"]] * 2, color="#0F172A", lw=1.5)
for v in range(0, 101, 10):
    ax_a.plot([v, v], [row_y["points"], row_y["points"] + 0.13], color="#0F172A", lw=0.9)
    ax_a.text(v, row_y["points"] + 0.20, str(v), ha="center", va="bottom", fontsize=7.2)

# Gene axes - tick density scales with axis length so short axes never crowd
for g in GENES:
    yy = row_y[g]
    p_span = gene_points[g]
    ax_a.text(-6, yy, g, ha="right", va="center", fontsize=9.5, fontweight="bold", color="#1D4ED8")
    ax_a.plot([0, p_span], [yy, yy], color="#2563EB", lw=1.4)
    n_ticks = int(np.clip(round(p_span / 14) + 1, 2, 6))
    for z in np.linspace(z_min[g], z_max[g], n_ticks):
        px = coefs[g] * (z - z_min[g]) * points_per_lp
        ax_a.plot([px, px], [yy - 0.12, yy + 0.12], color="#2563EB", lw=1.0)
        ax_a.text(px, yy - 0.26, f"{z:+.1f}", ha="center", va="top", fontsize=7.2, color="#334155")
    mu, sd = X_raw[g].mean(), X_raw[g].std()
    ax_a.text(-6, yy - 0.64, f"{mu:.2f} ± {sd:.2f}", ha="right", va="top",
              fontsize=6.8, color="#6B7280", style="italic")

ax_a.text(-6, row_y["points"] + 0.60, "grey: cohort mean ± SD (log₂)",
          ha="right", va="center", fontsize=6.8, color="#6B7280", style="italic")

# Total points ruler
ax_a.text(-6, row_y["total"], "Total Points", ha="right", va="center", fontsize=9.5, fontweight="bold")
ax_a.plot([0, 100], [row_y["total"]] * 2, color="#0F172A", lw=1.5)
tp_ticks = np.linspace(0, np.floor(total_points_max / 20) * 20, 6)
for v in tp_ticks:
    sx = v / total_points_max * 100
    ax_a.plot([sx, sx], [row_y["total"], row_y["total"] + 0.13], color="#0F172A", lw=0.9)
    ax_a.text(sx, row_y["total"] + 0.20, f"{int(v)}", ha="center", va="bottom", fontsize=7.2)

# Risk axis - exact mapping: LP = base_lp + total_points / points_per_lp.
# The risk band occupies only ~x = 41-74 of the total-points scale, so the ticks are chosen
# for even horizontal spacing (0.1/0.3/0.6/0.9) rather than even probability spacing.
ax_a.text(-6, row_y["risk"], "Risk of LPD", ha="right", va="center", fontsize=9.5,
          fontweight="bold", color="#DC2626")
ax_a.plot([0, 100], [row_y["risk"]] * 2, color="#DC2626", lw=1.0)
risk_ticks = [0.1, 0.3, 0.6, 0.9]
risk_x = []
for pr in risk_ticks:
    tp = (np.log(pr / (1 - pr)) - base_lp) * points_per_lp
    risk_x.append(tp / total_points_max * 100)
ax_a.plot([min(risk_x), max(risk_x)], [row_y["risk"]] * 2, color="#DC2626", lw=2.2)
for pr, sx in zip(risk_ticks, risk_x):
    ax_a.plot([sx, sx], [row_y["risk"] - 0.15, row_y["risk"] + 0.15], color="#DC2626", lw=1.2)
    ax_a.text(sx, row_y["risk"] - 0.28, f"{pr:g}", ha="center", va="top",
              fontsize=8.0, color="#B91C1C", fontweight="bold")

# ------------------------------------------------------------------- Panel B: calibration
ax_b = fig.add_subplot(gs[0, 1])
ax_b.plot([0, 1], [0, 1], ls="--", color="#475569", lw=1.3, label="Ideal")
ax_b.plot(cal_summary["boot_mean_pred"], cal_summary["boot_obs"], ls=":", color="#7C3AED",
          lw=1.6, marker="s", ms=4.5, label="Bootstrap-corrected")
ax_b.errorbar(cal_summary["mean_pred"], cal_summary["obs"],
              yerr=[cal_summary["obs"] - cal_summary["ci_low"],
                    cal_summary["ci_high"] - cal_summary["obs"]],
              fmt="o-", color=C_MODEL, ecolor="#F472B6", elinewidth=1.5, capsize=4,
              ms=6.5, lw=1.9, label="Apparent", zorder=3)
ax_b.set_xlim(-0.03, 1.03)
ax_b.set_ylim(-0.05, 1.05)
ax_b.set_xlabel("Predicted probability of LPD", fontsize=9.5)
ax_b.set_ylabel("Observed proportion with LPD", fontsize=9.5)
ax_b.set_title("Calibration Curve", fontsize=12, fontweight="bold", pad=14)
ax_b.legend(loc="lower right", frameon=False, fontsize=8.2)
ax_b.grid(ls=":", alpha=0.5, color="#CBD5E1")
ax_b.text(0.03, 0.90,
         f"Five equal-count bins (n = {', '.join(map(str, bin_counts))})\n"
         f"Exact binomial 95% CI; {N_BOOT:,} bootstrap resamples\n"
         f"Apparent AUC {apparent_auc:.3f} | LOO-CV AUC {loo_auc:.3f} | 5-fold CV AUC {cv_auc:.3f}",
         fontsize=7.4, color="#334155", va="top", linespacing=1.5)

# ------------------------------------------------------------------- Panel C: DCA
ax_c = fig.add_subplot(gs[1, 0])
ax_c.plot(thresholds, nb_model, color=C_MODEL, lw=2.2, label="Nomogram")
ax_c.plot(thresholds, nb_all, color=C_ALL, lw=1.7, ls="--", label="Treat all")
ax_c.axhline(0, color="#111827", lw=1.2, label="Treat none")
ax_c.fill_between(thresholds, 0, nb_model, where=nb_model > np.maximum(nb_all, 0),
                  color=C_MODEL, alpha=0.10)
ax_c.set_xlim(0, 1)
ax_c.set_ylim(-0.12, 0.45)
ax_c.set_xlabel("High-risk threshold\n(upper row: cost:benefit ratio)", fontsize=9.2, labelpad=26)
ax_c.set_ylabel("Net benefit", fontsize=9.5)
ax_c.set_title("Decision Curve Analysis", fontsize=12, fontweight="bold", pad=16)
ax_c.legend(loc="upper right", frameon=False, fontsize=8.5)
ax_c.grid(ls=":", alpha=0.5, color="#CBD5E1")
for pos, lab in zip([0.05, 0.2, 0.4, 0.6, 0.8, 0.95], ["1:20", "1:4", "2:3", "3:2", "4:1", "20:1"]):
    ax_c.text(pos, -0.115, lab, transform=ax_c.get_xaxis_transform(), fontsize=7.2,
              ha="center", va="top", color="#64748B")

# ------------------------------------------------------------------- Panel D: impact
ax_d = fig.add_subplot(gs[1, 1])
ax_d.plot(thresholds, impact_n, color=C_MODEL, lw=2.0, label="Number high risk")
ax_d.plot(thresholds, impact_true, color="#0891B2", lw=1.9, ls=":", label="Number high risk with event")
ax_d.set_xlim(0, 1)
ax_d.set_ylim(0, 1080)
ax_d.set_xlabel("High-risk threshold\n(upper row: cost:benefit ratio)", fontsize=9.2, labelpad=26)
ax_d.set_ylabel("Number per 1,000 patients", fontsize=9.5)
ax_d.set_title("Clinical Impact Curve", fontsize=12, fontweight="bold", pad=16)
ax_d.legend(loc="upper right", frameon=False, fontsize=8.5)
ax_d.grid(ls=":", alpha=0.5, color="#CBD5E1")
for pos, lab in zip([0.05, 0.2, 0.4, 0.6, 0.8, 0.95], ["1:20", "1:4", "2:3", "3:2", "4:1", "20:1"]):
    ax_d.text(pos, -0.115, lab, transform=ax_d.get_xaxis_transform(), fontsize=7.2,
              ha="center", va="top", color="#64748B")

# ------------------------------------------------------------------- Panel E: ssGSEA
ax_e = fig.add_subplot(gs[2, :])
pos = np.arange(len(CELL_TYPES))
lpd_data = [ssgsea.loc[ssgsea["group"] == "LPD", ct].values for ct in CELL_TYPES]
spd_data = [ssgsea.loc[ssgsea["group"] == "SPD", ct].values for ct in CELL_TYPES]
bp1 = ax_e.boxplot(lpd_data, positions=pos - 0.19, widths=0.33, patch_artist=True,
                   boxprops=dict(facecolor=C_LPD, edgecolor=C_LPD, alpha=0.60),
                   medianprops=dict(color="black", lw=1.3), showfliers=False)
bp2 = ax_e.boxplot(spd_data, positions=pos + 0.19, widths=0.33, patch_artist=True,
                   boxprops=dict(facecolor=C_SPD, edgecolor=C_SPD, alpha=0.60),
                   medianprops=dict(color="black", lw=1.3), showfliers=False)
rs = np.random.default_rng(SEED)
for i in range(len(CELL_TYPES)):
    ax_e.scatter(rs.normal(i - 0.19, 0.045, len(lpd_data[i])), lpd_data[i], s=9,
                 color=C_LPD, alpha=0.55, zorder=3)
    ax_e.scatter(rs.normal(i + 0.19, 0.045, len(spd_data[i])), spd_data[i], s=9,
                 color=C_SPD, alpha=0.55, zorder=3)
ymax = max(max(v) for v in lpd_data + spd_data)
# Significance is reported on a single, consistent convention: BH-FDR across the 13 cell types.
# Raw P values are retained in Fig5E_ssgsea_stats_source.csv.
for i, r in ssgsea_stats.iterrows():
    if not r["stars_fdr"]:
        continue
    top = max(max(lpd_data[i]), max(spd_data[i]))
    ax_e.text(i, top + 0.04, r["stars_fdr"], ha="center", va="bottom", fontsize=10,
              fontweight="bold", color="#111827")
ax_e.set_xticks(pos)
ax_e.set_xticklabels(CELL_TYPES, rotation=35, ha="right", fontsize=9)
ax_e.set_ylim(min(min(v) for v in lpd_data + spd_data) - 0.06, ymax + 0.16)
ax_e.set_ylabel("ssGSEA enrichment score", fontsize=9.5)
ax_e.set_title("Immune cell enrichment in LPD versus SPD", fontsize=12, fontweight="bold", pad=16)
ax_e.grid(ls=":", alpha=0.4, axis="y")
ax_e.legend([bp1["boxes"][0], bp2["boxes"][0], Patch(facecolor="none", edgecolor="none")],
            ["LPD", "SPD", "Mann-Whitney U, BH-FDR across 13 cell types"],
            loc="upper right", frameon=False, fontsize=8.5)

# ------------------------------------------------------------------- Panel F: immune-immune
ax_f = fig.add_subplot(gs[3, 0])
mat = imm_cor.loc[CELL_TYPES, CELL_TYPES].to_numpy()
mask = np.triu(np.ones_like(mat, dtype=bool), k=1)
im_f = ax_f.imshow(np.ma.masked_array(mat, mask=mask), cmap="RdBu_r", vmin=-1, vmax=1)
ax_f.set_xticks(range(len(CELL_TYPES)))
ax_f.set_yticks(range(len(CELL_TYPES)))
ax_f.set_xticklabels(CELL_TYPES, rotation=90, fontsize=7.4)
ax_f.set_yticklabels(CELL_TYPES, fontsize=7.4)
ax_f.set_title("Correlation between immune cell signatures", fontsize=11, fontweight="bold", pad=18)
cb_f = plt.colorbar(im_f, ax=ax_f, fraction=0.046, pad=0.03)
cb_f.set_label("Spearman ρ", fontsize=8)
cb_f.ax.tick_params(labelsize=7)
ax_f.spines[:].set_visible(False)

# ------------------------------------------------------------------- Panel G: gene-immune
ax_g = fig.add_subplot(gs[3, 1])
gc = ax_g.imshow(gene_rho.loc[GENES, CELL_TYPES].to_numpy(), cmap="RdBu_r", vmin=-0.8, vmax=0.8,
                 aspect="auto")
ax_g.set_xticks(range(len(CELL_TYPES)))
ax_g.set_yticks(range(len(GENES)))
ax_g.set_xticklabels(CELL_TYPES, rotation=90, fontsize=7.4)
ax_g.set_yticklabels(GENES, fontsize=8.5, fontweight="bold")
for i, g in enumerate(GENES):
    for j, ct in enumerate(CELL_TYPES):
        s = stars(gene_fdr.loc[g, ct])
        if s:
            ax_g.text(j, i, s, ha="center", va="center", fontsize=7.4, fontweight="bold",
                      color="white" if abs(gene_rho.loc[g, ct]) > 0.6 else "black")
ax_g.set_title("Correlation between hub genes and immune cells", fontsize=11,
               fontweight="bold", pad=18)
cb_g = plt.colorbar(gc, ax=ax_g, fraction=0.046, pad=0.03)
cb_g.set_label("Spearman ρ", fontsize=8)
cb_g.ax.tick_params(labelsize=7)

for ax, lab in zip([ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g], "ABCDEFG"):
    ax.text(-0.135, 1.145, lab, transform=ax.transAxes, fontsize=19, fontweight="bold",
            va="bottom", ha="left")

caption = textwrap.fill(
    "Figure 5. Nomogram based on the seven validated hub genes and its immune context. "
    "(A) Nomogram from a ridge logistic regression model fitted on GSE125498 (n = 33); gene axes are "
    "cohort z-scores and the italic grey text is the raw log₂ cohort mean ± SD. "
    "(B) Apparent and bootstrap-corrected calibration in five equal-count bins with exact binomial "
    "95% CIs, using 1,000 bootstrap resamples; apparent, LOOCV, and 5-fold CV AUCs were "
    f"{apparent_auc:.3f}, {loo_auc:.3f}, and {cv_auc:.3f}, respectively. "
    "(C) Decision curve analysis. (D) Clinical impact curve. "
    "(E) ssGSEA scores in LPD (red) versus SPD (blue); Mann-Whitney U, BH-FDR adjusted across the 13 cell types. "
    "(F) Spearman correlation among immune signatures. "
    "(G) Spearman correlation between hub genes and immune signatures; BH-FDR adjusted across all "
    "91 gene-cell-type tests. "
    "*FDR < 0.05, **FDR < 0.01, ***FDR < 0.001.",
    width=148,
)
fig.text(0.5, 0.010, caption, ha="center", va="bottom", fontsize=12, color="#0F172A",
         linespacing=1.6)

fig.subplots_adjust(left=0.075, right=0.975, top=0.955, bottom=0.155)
png = os.path.join(OUT_DIR, "Fig5_nomogram_immune_validated7.png")
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(png.replace(".png", ".pdf"), bbox_inches="tight")
plt.close(fig)

print(f"Figure 5 written: {png}")
print("\nPanel E statistics (raw and BH-FDR):")
print(ssgsea_stats[["cell_type", "median_LPD", "median_SPD", "delta_median",
                    "P_raw", "stars_raw", "P_FDR_BH", "survives_FDR"]].round(4).to_string(index=False))
print(f"\nDiscrimination: apparent AUC {apparent_auc:.4f} | LOO-CV {loo_auc:.4f} | 5-fold CV {cv_auc:.4f}")
print(f"Calibration bins: {bin_counts}")
print("\nNomogram ranking by point span:")
print(nomogram_source[["gene", "coefficient_z", "points_at_z_max"]].round(3).to_string(index=False))
