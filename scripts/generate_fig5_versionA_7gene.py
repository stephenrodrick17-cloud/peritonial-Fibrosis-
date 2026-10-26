"""
Figure 5 (Version A): Nomogram Construction and Association between Immune
Infiltration and Candidate Diagnostic Genes — based on the 7 candidate genes
with concordant direction in GSE248762 from the ml_39 (Version A) analysis.

9 hub genes initially identified:
  TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, CST7, NAT6, CCBE1.
NAT6 and CCBE1 excluded from the nomogram:
  In GSE248762 pseudobulk, NAT6 was detected (all-cell AUC = 0.500, P = 1.0)
  and CCBE1 was detected and inverted (AUC = 0.067, P = 0.003).

Panels:
  A: Nomogram point display (7-gene ridge model with discovery mean ± SD annotations)
  B: Calibration curve with 5 equal-count bins, exact binomial CIs, and bootstrap calibration slope [95% CI]
  C: Apparent Decision Curve Analysis (DCA)
  D: Apparent Clinical Impact Curve (scaled to 1,000 from n=33)
  E: ssGSEA enrichment score differences between LPD and SPD (displaying BH-adjusted q-values)
  F: Correlation between immune cell types (Spearman rho)
  G: Correlation between 7 candidate genes and immune cells (displaying BH-adjusted q-values across 91 tests)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.special import logit
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score

# ---------------------------------------------------------------------------
# PATHS & OUTPUT SETUP
# ---------------------------------------------------------------------------
ROOT = r"d:/Peritoneal Project"
OUT_DIR = os.path.join(ROOT, "results", "revision", "ml_39_matrisome_wgcna")
os.makedirs(OUT_DIR, exist_ok=True)

# Colors
C_LPD = "#F43F5E"     # Rose / Pink
C_SPD = "#0284C7"     # Cyan / Blue
C_ACCENT = "#EC4899"  # Accent Pink

GENES7 = ["TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "ECM1", "CST7"]
CELL_TYPES = [
    "CD4 T cell", "CD8 T cell", "NKT cell", "NK cell", "T cell (other)",
    "B cell", "Plasma cell", "Monocyte", "Macrophage", "Dendritic cell",
    "Neutrophil", "Eosinophil", "Mast cell"
]

def plabel(ax, label):
    ax.text(-0.12, 1.05, label, transform=ax.transAxes,
            fontsize=16, fontweight="bold", va="top", ha="right")

def qstar(q):
    if q != q:
        return "ns"
    if q < 1e-4:
        return "****"
    if q < 1e-3:
        return "***"
    if q < 1e-2:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"

# ---------------------------------------------------------------------------
# 1. LOAD DISCOVERY DATA (GSE125498)
# ---------------------------------------------------------------------------
print("Loading discovery data ...")
expr_df = pd.read_csv(os.path.join(ROOT, "results", "tables", "GSE125498_full_expression_matrix.csv"))
expr = expr_df.set_index("Gene")
meta = pd.read_csv(os.path.join(ROOT, "results", "tables", "GSE125498_sample_metadata.csv"))
meta["sample_id"] = meta["sample_id"].astype(str)
expr.columns = expr.columns.astype(str)

sc = [s for s in meta["sample_id"] if s in expr.columns]
y = np.array([int(meta.loc[meta["sample_id"] == s, "stage_binary"].iloc[0]) for s in sc])
n_samples = len(y)
n_events = int(y.sum())
n_nonevents = n_samples - n_events

X_raw = expr.loc[GENES7, sc].T.values.astype(float)
gene_means = X_raw.mean(axis=0)
gene_sds = X_raw.std(axis=0, ddof=0)
scaler = StandardScaler()
X = scaler.fit_transform(X_raw)

# ---------------------------------------------------------------------------
# 2. FIT 7-GENE RIDGE MODEL & IN-SAMPLE METRICS
# ---------------------------------------------------------------------------
print("Fitting ridge logistic regression (5-fold CV, neg_log_loss) ...")
lr_cv = LogisticRegressionCV(Cs=10, cv=5, penalty="l2", scoring="neg_log_loss",
                             random_state=42, max_iter=2000)
lr_cv.fit(X, y)
best_C = lr_cv.C_[0]
coefs = dict(zip(GENES7, lr_cv.coef_[0]))
intercept = float(lr_cv.intercept_[0])
pred = lr_cv.predict_proba(X)[:, 1]
apparent_auc = roc_auc_score(y, pred)

# LOOCV
loo = LeaveOneOut()
loo_preds = np.zeros(len(y))
for tr, te in loo.split(X, y):
    sc_i = StandardScaler().fit(X_raw[tr])
    X_tr = sc_i.transform(X_raw[tr])
    X_te = sc_i.transform(X_raw[te])
    lr_i = LogisticRegression(C=best_C, penalty="l2", solver="lbfgs", max_iter=10000)
    lr_i.fit(X_tr, y[tr])
    loo_preds[te] = lr_i.predict_proba(X_te)[:, 1]
loocv_auc = roc_auc_score(y, loo_preds)

# Bootstrap optimism and calibration slope (B=500)
rng = np.random.RandomState(42)
n_boot = 500
optimisms = []
boot_cal_slopes = []
boot_cal_intercepts = []
for b in range(n_boot):
    idx = rng.choice(len(y), len(y), replace=True)
    Xb = X_raw[idx]; yb = y[idx]
    if len(np.unique(yb)) < 2:
        continue
    sc_b = StandardScaler().fit(Xb)
    Xb_z = sc_b.transform(Xb)
    lr_b = LogisticRegression(C=best_C, penalty="l2", solver="lbfgs", max_iter=10000)
    lr_b.fit(Xb_z, yb)
    auc_boot = roc_auc_score(yb, lr_b.predict_proba(Xb_z)[:, 1])
    Xorig_z = sc_b.transform(X_raw)
    auc_test = roc_auc_score(y, lr_b.predict_proba(Xorig_z)[:, 1])
    optimisms.append(auc_boot - auc_test)
    
    lp_orig = lr_b.decision_function(Xorig_z)
    try:
        cal_lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=10000)
        cal_lr.fit(lp_orig.reshape(-1, 1), y)
        boot_cal_slopes.append(cal_lr.coef_[0][0])
        boot_cal_intercepts.append(cal_lr.intercept_[0])
    except Exception:
        pass

mean_optimism = np.mean(optimisms)
corrected_auc = apparent_auc - mean_optimism
mean_cal_slope = np.mean(boot_cal_slopes)
ci_slope_low = np.percentile(boot_cal_slopes, 2.5)
ci_slope_high = np.percentile(boot_cal_slopes, 97.5)
mean_cal_intercept = np.mean(boot_cal_intercepts)

print(f"Apparent AUC: {apparent_auc:.4f}, LOOCV AUC: {loocv_auc:.4f}, Optimism-corrected AUC: {corrected_auc:.4f}")
print(f"Bootstrap calibration slope: {mean_cal_slope:.2f} [95% CI: {ci_slope_low:.2f}-{ci_slope_high:.2f}]")

# ---------------------------------------------------------------------------
# 3. EXTERNAL VALIDATION IN GSE130888 (n = 10; 4 LV, 6 SV)
# ---------------------------------------------------------------------------
print("Evaluating frozen model on GSE130888 ...")
mapping = pd.read_csv(os.path.join(ROOT, "results/revision/single_cell_context/gene_id_symbol_map.csv"))
sym2ens = dict(zip(mapping["gene_symbol"], mapping["feature_id"]))
g130_cpm = pd.read_csv(os.path.join(ROOT, "results/revision/fibrosis_context/GSE130888/edger_logCPM.csv")).set_index("Unnamed: 0")
g130_samples = pd.read_csv(os.path.join(ROOT, "results/revision/ml_8_consensus_hub_validation/GSE130888/sample_group_mapping.csv"))

s_ids_130 = g130_samples["Sample_ID"].values
y_130 = np.array([1 if g == "LV" else 0 for g in g130_samples["Group"]])
ens_ids_7 = [sym2ens[g] for g in GENES7]
X_130_raw = g130_cpm.loc[ens_ids_7, s_ids_130].values.T

# Standardize GSE130888 using discovery parameters vs cohort-internal
X_130_disc_scale = scaler.transform(X_130_raw)
prob_130_disc = lr_cv.predict_proba(X_130_disc_scale)[:, 1]
auc_130_disc = roc_auc_score(y_130, prob_130_disc)

s_130_internal = StandardScaler()
X_130_internal = s_130_internal.fit_transform(X_130_raw)
prob_130_internal = lr_cv.predict_proba(X_130_internal)[:, 1]
auc_130_internal = roc_auc_score(y_130, prob_130_internal)

print(f"GSE130888 External AUC: {auc_130_disc:.4f} (scale caveat: discovery log-intensity vs RNA-seq logCPM)")

# ---------------------------------------------------------------------------
# 4. CALIBRATION BINS (PANEL B)
# ---------------------------------------------------------------------------
bins = np.quantile(pred, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
bins[0] -= 1e-4; bins[-1] += 1e-4
bin_indices = np.digitize(pred, bins) - 1
cal_pred = [float(pred[bin_indices == k].mean()) for k in range(5)]
cal_actual = [float(y[bin_indices == k].mean()) for k in range(5)]
cal_counts = [int((bin_indices == k).sum()) for k in range(5)]
cal_events = [int(y[bin_indices == k].sum()) for k in range(5)]

ci_low, ci_high = [], []
for k in range(5):
    lo, hi = stats.binom.interval(0.95, cal_counts[k], cal_actual[k])
    ci_low.append(max(0.0, lo / cal_counts[k]))
    ci_high.append(min(1.0, hi / cal_counts[k]))

mae_cal = float(np.mean(np.abs(np.array(cal_pred) - np.array(cal_actual))))

# ---------------------------------------------------------------------------
# 5. DCA & CLINICAL IMPACT (PANELS C & D)
# ---------------------------------------------------------------------------
thresholds = np.linspace(0.05, 0.95, 37)
N_dca = len(y)
dca_records = []
impact_records = []

for pt in thresholds:
    pred_pos = (pred >= pt).astype(int)
    tp = np.sum((y == 1) & (pred_pos == 1))
    fp = np.sum((y == 0) & (pred_pos == 1))
    nb_model = (tp / N_dca) - (fp / N_dca) * (pt / (1.0 - pt))
    
    tp_all = np.sum(y == 1)
    fp_all = np.sum(y == 0)
    nb_all = (tp_all / N_dca) - (fp_all / N_dca) * (pt / (1.0 - pt))
    
    dca_records.append({'Threshold': pt, 'Strategy': 'Seven-gene model', 'Net_benefit': nb_model})
    dca_records.append({'Threshold': pt, 'Strategy': 'Treat all', 'Net_benefit': nb_all})
    dca_records.append({'Threshold': pt, 'Strategy': 'Treat none', 'Net_benefit': 0.0})
    
    impact_records.append({'Threshold': pt, 'High_risk_n': int(pred_pos.sum()), 'True_LPD_in_high_risk': int(tp)})

dca_df = pd.DataFrame(dca_records)
impact_df = pd.DataFrame(impact_records)

# ---------------------------------------------------------------------------
# 6. IMMUNE DATA & MULTIPLE TESTING ADJUSTMENTS (PANELS E, F, G)
# ---------------------------------------------------------------------------
imm_df = pd.read_csv(os.path.join(ROOT, "results", "paper_style", "Fig5E_saved_ssgsea_sample_scores.csv")).set_index("sample_id")
imm_matched = imm_df.loc[sc]
imm_scores = imm_matched[CELL_TYPES].values.astype(float)

# Panel E: Mann-Whitney tests + BH FDR correction across 13 tests
immune_pvals = []
for i, ct in enumerate(CELL_TYPES):
    v_l = imm_scores[y == 1, i]
    v_s = imm_scores[y == 0, i]
    _, p = stats.mannwhitneyu(v_l, v_s, alternative="two-sided")
    immune_pvals.append(p)

rej_e, immune_qvals, _, _ = multipletests(immune_pvals, method="fdr_bh")
immune_q_dict = dict(zip(CELL_TYPES, immune_qvals))
immune_p_dict = dict(zip(CELL_TYPES, immune_pvals))
n_sig_q = int(rej_e.sum())
print(f"Panel E: {n_sig_q}/13 immune cell types significant at BH q < 0.05")

# Panel F: Immune-immune correlation
imm_corr = np.zeros((len(CELL_TYPES), len(CELL_TYPES)))
imm_pmat = np.zeros((len(CELL_TYPES), len(CELL_TYPES)))
for i in range(len(CELL_TYPES)):
    for j in range(len(CELL_TYPES)):
        r, p = stats.spearmanr(imm_scores[:, i], imm_scores[:, j])
        imm_corr[i, j] = r
        imm_pmat[i, j] = p

# Panel G: Gene-immune correlation (91 tests) + BH FDR correction
gene_imm_corr = np.zeros((len(GENES7), len(CELL_TYPES)))
gene_imm_p = np.zeros((len(GENES7), len(CELL_TYPES)))
for i, g in enumerate(GENES7):
    g_vals = expr.loc[g, sc].values.astype(float)
    for j in range(len(CELL_TYPES)):
        r, p = stats.spearmanr(g_vals, imm_scores[:, j])
        gene_imm_corr[i, j] = r
        gene_imm_p[i, j] = p

p_g_flat = gene_imm_p.flatten()
rej_g, q_g_flat, _, _ = multipletests(p_g_flat, method="fdr_bh")
gene_imm_q = q_g_flat.reshape(len(GENES7), len(CELL_TYPES))
n_sig_g_fdr = int(rej_g.sum())
print(f"Panel G: {n_sig_g_fdr}/91 gene-immune correlations significant at BH q < 0.05 (raw p < 0.05: {(p_g_flat < 0.05).sum()})")

# ---------------------------------------------------------------------------
# 7. GENERATE COMPOSITE FIGURE
# ---------------------------------------------------------------------------
print("Generating composite figure ...")
fig = plt.figure(figsize=(14.0, 20.5), facecolor="#FFFFFF")
gs = gridspec.GridSpec(4, 2, height_ratios=[1.25, 0.95, 1.10, 1.15],
                       hspace=0.38, wspace=0.25)

# ---------------- PANEL A: Nomogram Display ----------------
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_xlim(-22, 115)
ax_a.set_ylim(-2.0, 10.0)
ax_a.axis("off")

z_spans = {g: (float(X[:, i].min()), float(X[:, i].max())) for i, g in enumerate(GENES7)}
lp_spans = {g: abs(coefs[g]) * (z_spans[g][1] - z_spans[g][0]) for g in GENES7}
max_span = max(lp_spans.values())
pts_per_lp = 100.0 / max_span

total_pts_min = sum(coefs[g] * z_spans[g][0] for g in GENES7 if coefs[g] > 0)
total_pts_max = sum(coefs[g] * z_spans[g][1] for g in GENES7 if coefs[g] > 0)
total_max_points = sum(lp_spans.values()) * pts_per_lp

y_pos = {
    "Points": 9.2,
    "TNFSF15": 8.0,
    "FLT3LG": 6.9,
    "EBI3": 5.8,
    "LTB": 4.7,
    "ADAM19": 3.6,
    "ECM1": 2.5,
    "CST7": 1.4,
    "Total Points": -0.1,
    "Risk": -1.3
}

# Points axis
ax_a.plot([0, 100], [y_pos["Points"], y_pos["Points"]], color="#000000", lw=1.2)
for pt in [0, 20, 40, 60, 80, 100]:
    ax_a.plot([pt, pt], [y_pos["Points"], y_pos["Points"] + 0.14], color="#000000", lw=1.0)
    ax_a.text(pt, y_pos["Points"] + 0.22, str(pt), ha="center", va="bottom", fontsize=7.2)
ax_a.text(-2, y_pos["Points"], "Points", ha="right", va="center", fontsize=8.8, fontweight="bold")

# Gene axes
for g in GENES7:
    yp = y_pos[g]
    span = lp_spans[g] * pts_per_lp
    ax_a.plot([0, span], [yp, yp], color="#0284C7", lw=1.2)
    z_min, z_max = z_spans[g]
    g_idx = GENES7.index(g)
    m_val = gene_means[g_idx]
    s_val = gene_sds[g_idx]
    
    ticks = np.linspace(z_min, z_max, 5)
    for tk in ticks:
        x_tk = (tk - z_min) / (z_max - z_min) * span
        ax_a.plot([x_tk, x_tk], [yp, yp - 0.12], color="#0284C7", lw=1.0)
        ax_a.text(x_tk, yp - 0.20, f"{tk:+.1f}", ha="center", va="top", fontsize=6.2)
    
    label_text = f"{g}\n({m_val:.1f}±{s_val:.1f})"
    ax_a.text(-2, yp, label_text, ha="right", va="center", fontsize=7.8,
              fontweight="bold", color="#0369A1")

# Total points axis
tp_max_disp = 400
ax_a.plot([0, 100], [y_pos["Total Points"], y_pos["Total Points"]], color="#000000", lw=1.2)
for tp in np.arange(0, tp_max_disp + 1, 40):
    x_tp = (tp / tp_max_disp) * 100.0
    ax_a.plot([x_tp, x_tp], [y_pos["Total Points"], y_pos["Total Points"] + 0.14], color="#000000", lw=1.0)
    ax_a.text(x_tp, y_pos["Total Points"] + 0.22, str(tp), ha="center", va="bottom", fontsize=6.8)
ax_a.text(-2, y_pos["Total Points"], "Total Points", ha="right", va="center", fontsize=8.8, fontweight="bold")

# Risk axis
ax_a.plot([0, 100], [y_pos["Risk"], y_pos["Risk"]], color="#DC2626", lw=1.4)
ax_a.text(-2, y_pos["Risk"], "Risk of LPD", ha="right", va="center",
          fontsize=8.8, fontweight="bold", color="#DC2626")

risk_probs = [0.1, 0.3, 0.5, 0.7, 0.9]
for pr in risk_probs:
    lp_val = np.log(pr / (1.0 - pr))
    pts = (lp_val - intercept) * pts_per_lp / max_span
    x_pr = (pts / tp_max_disp) * 100.0
    if 0 <= x_pr <= 100:
        ax_a.plot([x_pr, x_pr], [y_pos["Risk"], y_pos["Risk"] - 0.14], color="#DC2626", lw=1.2)
        ax_a.text(x_pr, y_pos["Risk"] - 0.22, f"{pr:.1f}", ha="center", va="top",
                  fontsize=7.0, color="#DC2626", fontweight="bold")

ax_a.set_title("Constructed Nomogram", fontsize=11.5, fontweight="bold", pad=8)
plabel(ax_a, "A")

# ---------------- PANEL B: Apparent Calibration Curve ----------------
ax_b = fig.add_subplot(gs[0, 1])
ax_b.plot([0, 1], [0, 1], linestyle="--", color="#64748B", lw=1.2, label="Ideal")
ax_b.errorbar(cal_pred, cal_actual,
              yerr=[np.array(cal_actual) - np.array(ci_low),
                    np.array(ci_high) - np.array(cal_actual)],
              fmt="o-", color=C_ACCENT, ecolor=C_ACCENT, elinewidth=1.2,
              capsize=3, markersize=6, lw=1.8, label="Apparent (in-sample)")
ax_b.set_xlim(-0.02, 1.02)
ax_b.set_ylim(-0.02, 1.02)
ax_b.set_xlabel("Predicted Pr (LPD = 1)", fontsize=9.5)
ax_b.set_ylabel("Actual Probability", fontsize=9.5)
ax_b.set_title("Apparent Calibration Curve", fontsize=11.5, fontweight="bold", pad=8)
ax_b.legend(loc="lower right", frameon=False, fontsize=8)
ax_b.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")

stat_text = (
    f"Discovery: n=33 (13 LPD, 20 SPD)\n"
    f"EPV = {n_events/len(GENES7):.1f}\n"
    f"Apparent AUC = {apparent_auc:.3f}\n"
    f"LOOCV AUC = {loocv_auc:.3f}\n"
    f"Optimism-corrected AUC = {corrected_auc:.3f}\n"
    f"Bootstrap slope = {mean_cal_slope:.2f} [{ci_slope_low:.2f}-{ci_slope_high:.2f}]\n"
    f"External GSE130888 AUC = {auc_130_disc:.3f}*"
)
ax_b.text(0.04, 0.96, stat_text, transform=ax_b.transAxes, fontsize=7.2,
          bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8FAFC", edgecolor="#E2E8F0"),
          verticalalignment="top")
plabel(ax_b, "B")

# ---------------- PANEL C: DCA ----------------
ax_c = fig.add_subplot(gs[1, 0])
dca_m = dca_df[dca_df["Strategy"] == "Seven-gene model"].sort_values("Threshold")
dca_a = dca_df[dca_df["Strategy"] == "Treat all"].sort_values("Threshold")
dca_n = dca_df[dca_df["Strategy"] == "Treat none"].sort_values("Threshold")
ax_c.plot(dca_m["Threshold"], dca_m["Net_benefit"], color=C_ACCENT, lw=2.2,
          label="Nomogram (apparent)")
ax_c.plot(dca_a["Threshold"], dca_a["Net_benefit"], color="#94A3B8", lw=1.6, label="All")
ax_c.plot(dca_n["Threshold"], dca_n["Net_benefit"], color="#000000", lw=1.2, label="None")
ax_c.fill_between(dca_m["Threshold"].values, 0, dca_m["Net_benefit"].values,
                  where=dca_m["Net_benefit"].values > np.maximum(dca_a["Net_benefit"].values, 0),
                  color=C_ACCENT, alpha=0.12, interpolate=True)
ax_c.set_xlim(0.0, 1.0)
ax_c.set_ylim(-0.05, max(0.45, dca_m["Net_benefit"].max() + 0.05))
ax_c.set_xlabel("High Risk Threshold", fontsize=9.5)
ax_c.set_ylabel("Net Benefit", fontsize=9.5)
ax_c.set_title("Apparent Decision Curve Analysis (DCA)", fontsize=11.5, fontweight="bold", pad=8)
ax_c.legend(loc="upper right", frameon=False, fontsize=8)
ax_c.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
plabel(ax_c, "C")

for xv, lbl in [(0.0, "1:100"), (0.2, "1:4"), (0.4, "2:3"), (0.6, "3:2"),
                (0.8, "4:1"), (1.0, "100:1")]:
    ax_c.text(xv, -0.16, lbl, transform=ax_c.get_xaxis_transform(),
              fontsize=7.0, ha="center")
ax_c.text(0.5, -0.25, "Cost:Benefit Ratio", transform=ax_c.get_xaxis_transform(),
          fontsize=8.0, ha="center")

# ---------------- PANEL D: Clinical Impact Curve ----------------
ax_d = fig.add_subplot(gs[1, 1])
th_vals = impact_df["Threshold"].values
high_n = impact_df["High_risk_n"].values * (1000.0 / y.shape[0])
high_true = impact_df["True_LPD_in_high_risk"].values * (1000.0 / y.shape[0])
ax_d.plot(th_vals, high_n, color=C_ACCENT, lw=2.0, label="Number high risk")
ax_d.plot(th_vals, high_true, color="#06B6D4", lw=1.8, linestyle=":",
          label="Number high risk with event")
ax_d.set_xlim(0.0, 1.0)
ax_d.set_ylim(0, 1050)
ax_d.set_xlabel("High Risk Threshold", fontsize=9.5)
ax_d.set_ylabel("Number high risk (per 1000)", fontsize=9.5)
ax_d.set_title("Apparent Clinical Impact Curve (scaled from n=33)",
               fontsize=11.5, fontweight="bold", pad=8)
ax_d.legend(loc="upper right", frameon=False, fontsize=8)
ax_d.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
plabel(ax_d, "D")

for xv, lbl in [(0.0, "1:100"), (0.2, "1:4"), (0.4, "2:3"), (0.6, "3:2"),
                (0.8, "4:1"), (1.0, "100:1")]:
    ax_d.text(xv, -0.16, lbl, transform=ax_d.get_xaxis_transform(),
              fontsize=7.0, ha="center")
ax_d.text(0.5, -0.25, "Cost:Benefit Ratio", transform=ax_d.get_xaxis_transform(),
          fontsize=8.0, ha="center")

# ---------------- PANEL E: ssGSEA Enrichment Score Differences ----------------
ax_e = fig.add_subplot(gs[2, :])
np.random.seed(42)
for i, ct in enumerate(CELL_TYPES):
    v_l = imm_scores[y == 1, i]
    v_s = imm_scores[y == 0, i]
    pos_l = i * 2.2 + 0.8
    pos_s = i * 2.2 + 1.5
    bp_l = ax_e.boxplot(v_l, positions=[pos_l], widths=0.55, patch_artist=True,
                        boxprops=dict(facecolor=C_LPD, color=C_LPD, alpha=0.75),
                        medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    bp_s = ax_e.boxplot(v_s, positions=[pos_s], widths=0.55, patch_artist=True,
                        boxprops=dict(facecolor=C_SPD, color=C_SPD, alpha=0.75),
                        medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    ax_e.scatter(np.random.normal(pos_l, 0.06, len(v_l)), v_l, color=C_LPD, s=8, alpha=0.6)
    ax_e.scatter(np.random.normal(pos_s, 0.06, len(v_s)), v_s, color=C_SPD, s=8, alpha=0.6)
    y_top = max(float(v_l.max()), float(v_s.max())) + 0.03
    star = qstar(immune_q_dict[ct])
    if star != "ns":
        ax_e.text(i * 2.2 + 1.15, y_top, star, ha="center", va="bottom",
                  fontsize=9, fontweight="bold")

ax_e.set_xticks([i * 2.2 + 1.15 for i in range(len(CELL_TYPES))])
ax_e.set_xticklabels(CELL_TYPES, rotation=35, ha="right", fontsize=8.2, fontweight="medium")
ax_e.set_ylabel("ssGSEA Enrichment Score", fontsize=9.5)
ax_e.set_title("ssGSEA Enrichment Score Differences between LPD and SPD (BH-adjusted q)",
               fontsize=11.5, fontweight="bold", pad=8)
ax_e.grid(True, linestyle=":", alpha=0.4, axis="y")
h_l = patches.Patch(facecolor=C_LPD, label="LPD")
h_s = patches.Patch(facecolor=C_SPD, label="SPD")
ax_e.legend(handles=[h_l, h_s], loc="upper right", frameon=False, fontsize=8.5)
plabel(ax_e, "E")

# ---------------- PANEL F: Immune-Immune Correlation ----------------
ax_f = fig.add_subplot(gs[3, 0])
mask_tri = np.triu(np.ones_like(imm_corr, dtype=bool), k=1)
imm_masked = np.ma.masked_array(imm_corr, mask=mask_tri)
im_f = ax_f.imshow(imm_masked, cmap="coolwarm", vmin=-1.0, vmax=1.0, aspect="auto")
ax_f.set_xticks(range(len(CELL_TYPES)))
ax_f.set_xticklabels(CELL_TYPES, rotation=60, ha="right", fontsize=6.8)
ax_f.set_yticks(range(len(CELL_TYPES)))
ax_f.set_yticklabels(CELL_TYPES, fontsize=6.8)
ax_f.set_title("Correlation Between Immune Cells", fontsize=10.5, fontweight="bold", pad=8)
ax_f.spines[:].set_visible(False)
for i in range(len(CELL_TYPES)):
    for j in range(i):
        p = imm_pmat[i, j]
        if p < 0.001:
            ax_f.text(j, i, "***", ha="center", va="center", color="white",
                      fontsize=5.5, fontweight="bold")
        elif p < 0.01:
            ax_f.text(j, i, "**", ha="center", va="center", color="white",
                      fontsize=5.5, fontweight="bold")
        elif p < 0.05:
            ax_f.text(j, i, "*", ha="center", va="center", color="white",
                      fontsize=5.5, fontweight="bold")
cb_f = plt.colorbar(im_f, ax=ax_f, fraction=0.046, pad=0.04)
cb_f.set_label("Spearman Rho", fontsize=8)
plabel(ax_f, "F")

# ---------------- PANEL G: Gene-Immune Correlation (BH FDR Annotated) ----------------
ax_g = fig.add_subplot(gs[3, 1])
im_g = ax_g.imshow(gene_imm_corr, cmap="coolwarm", vmin=-0.8, vmax=0.8, aspect="auto")
ax_g.set_xticks(range(len(CELL_TYPES)))
ax_g.set_xticklabels(CELL_TYPES, rotation=60, ha="right", fontsize=6.8)
ax_g.set_yticks(range(len(GENES7)))
ax_g.set_yticklabels(GENES7, fontsize=8.5, fontweight="bold")
for i in range(len(GENES7)):
    for j in range(len(CELL_TYPES)):
        q = gene_imm_q[i, j]
        star = qstar(q)
        if star != "ns":
            ax_g.text(j, i, star, ha="center", va="center", color="#000000",
                      fontsize=6.8, fontweight="bold")

ax_g.set_title("Correlation Between Candidate Genes and Immune Cells (BH q)",
               fontsize=10.5, fontweight="bold", pad=8)
cb_g = plt.colorbar(im_g, ax=ax_g, fraction=0.046, pad=0.04)
cb_g.set_label("Spearman Rho", fontsize=8)
plabel(ax_g, "G")

# ---------------- CAPTION ----------------
fig.text(0.5, 0.008,
         "Fig. 5. Nomogram construction and association between immune cell infiltration and candidate diagnostic genes\n"
         "(TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, CST7 \u2014 the 7 candidate genes concordant in GSE248762 out of 9 Version A hub genes;\n"
         "NAT6 was detected [all-cell AUC = 0.500, P = 1.0] and CCBE1 was detected and inverted [AUC = 0.067, P = 0.003], and both were dropped after viewing verification data).\n"
         "(A) Constructed nomogram (7-gene ridge model; gene labels display discovery mean \u00b1 SD). "
         "(B) Apparent (in-sample) calibration curve (5 equal-count bins with exact binomial CIs; bootstrap calibration slope = 1.96 [95% CI 1.28–3.41]; *external frozen model AUC on GSE130888 = 0.542, subject to platform/scale differences).\n"
         "(C) Apparent DCA evaluating clinical net benefit vs treat-all/treat-none. "
         "(D) Apparent clinical impact curve (scaled to 1,000 from discovery n = 33).\n"
         "(E) ssGSEA enrichment score differences between LPD and SPD (6 of 13 cell types survive BH adjustment at q < 0.05: CD4 T cell, CD8 T cell, NKT cell, NK cell, T cell other, B cell; plasma cell q = 0.079, dendritic cell q = 0.069). "
         "(F) Spearman correlation between immune cell types.\n"
         "(G) Spearman correlation between 7 candidate genes and immune cell scores (38 of 91 tests survive BH FDR adjustment across the panel; note that 6 of 7 candidate genes appear as markers in CellMarker signatures, e.g. CST7 and LTB in lymphocyte sets).\n"
         f"All panels: apparent (in-sample) discovery estimates; n = {n_samples} ({n_events} LPD, {n_nonevents} SPD); 7 predictors. "
         "*q < 0.05, **q < 0.01, ***q < 0.001, ****q < 0.0001 (Benjamini–Hochberg adjusted).",
         ha="center", va="bottom", fontsize=8.2, color="#0F172A", linespacing=1.35)

plt.subplots_adjust(bottom=0.09, top=0.965)

png_out = os.path.join(OUT_DIR, "Fig5_versionA_7gene_nomogram_immune.png")
pdf_out = os.path.join(OUT_DIR, "Fig5_versionA_7gene_nomogram_immune.pdf")
fig.savefig(png_out, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"\nSaved PNG: {png_out}")
print(f"Saved PDF: {pdf_out}")
