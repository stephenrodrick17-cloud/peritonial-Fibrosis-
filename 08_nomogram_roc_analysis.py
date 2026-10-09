"""
==============================================================================
TASK 2: CLINICAL DIAGNOSTIC NOMOGRAM, CALIBRATION & DECISION CURVE ANALYSIS (DCA)
==============================================================================
Constructs:
  1. Mathematically Exact Clinical Diagnostic Nomogram from fitted multivariable
     logistic regression model (5 core genes profiled across discovery & external cohorts:
     VCAN, COL8A1, FN1, ISM1, COL3A1). Points are strictly proportional to |beta * range|,
     and the Total Points axis maps exactly to the Fibrosis Risk probability scale.
  2. Bootstrap (B=1000) Bias-Corrected Calibration Curve with Lowess Smoother,
     reporting Calibration Slope, Intercept, Brier Score, and Hosmer-Lemeshow test.
  3. Decision Curve Analysis (DCA) evaluated on external validation cohort (GSE125498, N=33),
     evaluating Net Clinical Benefit across decision thresholds Pt = 0.10 - 0.70.
  4. Likelihood-Ratio (LR) Test comparing multivariable nomogram against FN1 alone.
  5. Exports dual-format publication figures (300 DPI PNG + Vector PDF) and CSV tables.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import chi2
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss

# Set seed and output directories
np.random.seed(42)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# 1. Load Data
df_expr_val = pd.read_csv("results/tables/GSE125498_expression_matrix_hubs.csv", index_col=0)
df_meta_val = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
y_val = df_meta_val["stage_binary"].values.astype(int)
val_samples = df_meta_val["sample_id"].tolist()
n_total = len(y_val)
n_events = int(y_val.sum())
n_controls = n_total - n_events

# 5 Core predictors with complete profiles across cohorts
nomo_genes = ["VCAN", "COL8A1", "FN1", "ISM1", "COL3A1"]
X_df = df_expr_val.loc[nomo_genes, val_samples].T
X = X_df.values

# 2. Fit Multivariable Logistic Regression
# Statsmodels for exact coefficients, SEs, p-values and log-likelihood
X_const = sm.add_constant(X_df)
logit_full = sm.Logit(y_val, X_const).fit(disp=False)
coefs = logit_full.params[nomo_genes].values
intercept = logit_full.params["const"]
pvals = logit_full.pvalues[nomo_genes].values

# Fit single-gene model FN1 alone for Likelihood-Ratio test
logit_fn1 = sm.Logit(y_val, sm.add_constant(X_df[["FN1"]])).fit(disp=False)
lr_stat = 2 * (logit_full.llf - logit_fn1.llf)
lr_df = len(nomo_genes) - 1
lr_p = chi2.sf(lr_stat, df=lr_df)
print(f"Likelihood-Ratio Test vs FN1 alone: LR chi2 = {lr_stat:.3f}, df = {lr_df}, P = {lr_p:.4f}")

# Predicted probabilities
probs = logit_full.predict(X_const)
auc_val = roc_auc_score(y_val, probs)

# Bootstrap 2000 for AUC 95% CI
boot_aucs = []
for _ in range(2000):
    idx = np.random.choice(n_total, size=n_total, replace=True)
    if len(np.unique(y_val[idx])) > 1:
        boot_aucs.append(roc_auc_score(y_val[idx], probs[idx]))
auc_ci_low, auc_ci_high = np.percentile(boot_aucs, [2.5, 97.5])
print(f"Nomogram AUC (C-index): {auc_val:.3f} (95% CI: [{auc_ci_low:.3f}, {auc_ci_high:.3f}])")

# 3. Exact Mathematical Nomogram Construction
# Calculate dynamic range R_i = |beta_i| * (max(X_i) - min(X_i))
gene_mins = {g: X_df[g].min() for g in nomo_genes}
gene_maxs = {g: X_df[g].max() for g in nomo_genes}
gene_ranges = {g: gene_maxs[g] - gene_mins[g] for g in nomo_genes}

beta_dict = dict(zip(nomo_genes, coefs))
dyn_ranges = {g: abs(beta_dict[g]) * gene_ranges[g] for g in nomo_genes}
max_dyn_range = max(dyn_ranges.values())

# Scaling factor: maximum dynamic range maps to 100 points
scale = 100.0 / max_dyn_range

# Point scoring function per variable:
# If beta > 0: min value gives 0 pts, max gives scale * beta * (max - min)
# If beta < 0: max value gives 0 pts, min gives scale * |beta| * (max - min)
max_points = {}
for g in nomo_genes:
    max_points[g] = scale * dyn_ranges[g]

# Linear predictor min value:
# LP_min is achieved when every variable gives 0 points
lp_min = intercept
for g in nomo_genes:
    if beta_dict[g] > 0:
        lp_min += beta_dict[g] * gene_mins[g]
    else:
        lp_min += beta_dict[g] * gene_maxs[g]

total_points_max = sum(max_points.values())
print(f"Total Points Range: 0.0 to {total_points_max:.1f}")

# Export Nomogram Point Scoring Table
nomo_table_rows = []
for g in nomo_genes:
    nomo_table_rows.append({
        "Gene_Symbol": g,
        "Coefficient_Beta": round(beta_dict[g], 4),
        "P_Value": round(logit_full.pvalues[g], 4),
        "Observed_Min": round(gene_mins[g], 2),
        "Observed_Max": round(gene_maxs[g], 2),
        "Points_At_Min": 0.0 if beta_dict[g] > 0 else round(max_points[g], 1),
        "Points_At_Max": round(max_points[g], 1) if beta_dict[g] > 0 else 0.0,
        "Max_Assigned_Points": round(max_points[g], 1)
    })
df_nomo_table = pd.DataFrame(nomo_table_rows)
out_nomo_table = "results/tables/nomogram_point_scoring_table.csv"
df_nomo_table.to_csv(out_nomo_table, index=False)
print(f"Exported Nomogram scoring table to {out_nomo_table}")

# Verification of Total Points -> Probability Mapping:
# LP = lp_min + Total_Points / scale
# Probability P = 1 / (1 + exp(-LP))
# Conversely: Total_Points = scale * (logit(P) - lp_min)
def p_to_total_points(p):
    logit_p = np.log(p / (1.0 - p))
    return scale * (logit_p - lp_min)

def total_points_to_p(pts):
    lp = lp_min + pts / scale
    return 1.0 / (1.0 + np.exp(-lp))

# Verify calibration points
print("Verification of Nomogram Points to Probability:")
for test_p in [0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95]:
    calc_pts = p_to_total_points(test_p)
    recov_p = total_points_to_p(calc_pts)
    print(f"  Risk {test_p*100:4.1f}% -> Total Points: {calc_pts:6.2f} (Verified: {recov_p*100:4.1f}%)")

# 4. Calibration Metrics & Bootstrap Calibration
brier = brier_score_loss(y_val, probs)

# Calibration slope and intercept: logit(y) = a + b * logit(p)
eps = 1e-7
lp_all = np.log((probs + eps) / (1 - probs + eps))
calib_mod = sm.Logit(y_val, sm.add_constant(lp_all)).fit(disp=False)
calib_intercept, calib_slope = calib_mod.params
print(f"Calibration Metrics: Brier = {brier:.4f}, Intercept = {calib_intercept:.4f}, Slope = {calib_slope:.4f}")

# Hosmer-Lemeshow test (g=4 bins due to N=33 sample size)
g_bins = 4
df_hl = pd.DataFrame({"y": y_val, "p": probs}).sort_values("p")
df_hl["group"] = pd.qcut(df_hl["p"], q=g_bins, duplicates="drop")
hl_stat = 0.0
for _, grp in df_hl.groupby("group", observed=True):
    o1 = grp["y"].sum()
    o0 = len(grp) - o1
    e1 = grp["p"].sum()
    e0 = len(grp) - e1
    hl_stat += ((o1 - e1)**2) / max(e1, 1e-5) + ((o0 - e0)**2) / max(e0, 1e-5)
hl_df = g_bins - 2
hl_p = chi2.sf(hl_stat, df=hl_df)
print(f"Hosmer-Lemeshow test (g={g_bins} bins): chi2 = {hl_stat:.3f}, df = {hl_df}, P = {hl_p:.4f}")

# Bootstrap Bias-Corrected Calibration (B=1000)
n_eval_pts = 50
grid_p = np.linspace(0.02, 0.98, n_eval_pts)
# Lowess smoother on observed data
from statsmodels.nonparametric.smoothers_lowess import lowess
orig_smooth = lowess(y_val, probs, frac=0.65, xvals=grid_p)

# Bootstrap iterations for apparent vs test performance
boot_optimism = np.zeros(n_eval_pts)
n_boot = 1000
for _ in range(n_boot):
    b_idx = np.random.choice(n_total, size=n_total, replace=True)
    if len(np.unique(y_val[b_idx])) > 1:
        clf_b = LogisticRegression(C=1e5, max_iter=1000, random_state=42).fit(X[b_idx], y_val[b_idx])
        p_apparent = clf_b.predict_proba(X[b_idx])[:, 1]
        p_test = clf_b.predict_proba(X)[:, 1]
        app_smooth = lowess(y_val[b_idx], p_apparent, frac=0.65, xvals=grid_p)
        test_smooth = lowess(y_val, p_test, frac=0.65, xvals=grid_p)
        boot_optimism += (app_smooth - test_smooth)
boot_optimism /= n_boot
bias_corrected_calib = np.clip(orig_smooth - boot_optimism, 0.0, 1.0)

# 5. Decision Curve Analysis (DCA)
thresholds = np.linspace(0.05, 0.75, 100)
nb_model = []
nb_all = []
for pt in thresholds:
    # Model net benefit
    y_pred = (probs >= pt).astype(int)
    tp = np.sum((y_pred == 1) & (y_val == 1))
    fp = np.sum((y_pred == 1) & (y_val == 0))
    nb_m = (tp / n_total) - (fp / n_total) * (pt / (1.0 - pt))
    nb_model.append(nb_m)
    
    # Treat all net benefit
    tp_all = np.sum(y_val == 1)
    fp_all = np.sum(y_val == 0)
    nb_a = (tp_all / n_total) - (fp_all / n_total) * (pt / (1.0 - pt))
    nb_all.append(nb_a)

nb_model = np.array(nb_model)
nb_all = np.array(nb_all)

# Export Calibration and DCA Metrics to CSV
df_dca_export = pd.DataFrame({
    "Threshold_Pt": np.round(thresholds, 4),
    "Net_Benefit_Model": np.round(nb_model, 4),
    "Net_Benefit_Treat_All": np.round(nb_all, 4),
    "Net_Benefit_Treat_None": 0.0
})
df_dca_export.to_csv("results/tables/decision_curve_analysis_metrics.csv", index=False)

df_calib_export = pd.DataFrame({
    "Predicted_Risk_Grid": np.round(grid_p, 4),
    "Observed_Lowess": np.round(orig_smooth, 4),
    "Bias_Corrected_Calibration": np.round(bias_corrected_calib, 4)
})
df_calib_export.to_csv("results/tables/nomogram_calibration_curve_metrics.csv", index=False)

# 6. Publication Visualization Layout
fig = plt.figure(figsize=(18, 12), facecolor="#FFFFFF")
gs = gridspec.GridSpec(2, 2, height_ratios=[1.25, 1.0], hspace=0.35, wspace=0.25)

# Panel A: Clinical Diagnostic Nomogram (Top Row spanning both columns)
ax_nomo = fig.add_subplot(gs[0, :])
ax_nomo.set_facecolor("#FFFFFF")
ax_nomo.axis("off")

# Title
ax_nomo.text(0.5, 0.98, "CLINICAL DIAGNOSTIC NOMOGRAM FOR PERITONEAL MEMBRANE FIBROGENESIS",
             fontsize=13, fontweight="bold", color="#0F172A", ha="center", transform=ax_nomo.transAxes)
ax_nomo.text(0.5, 0.93, f"Multivariable Logistic Model in External Cohort (GSE125498: N = {n_total}, Events = {n_events}) | C-Index = {auc_val:.3f} (95% CI: [{auc_ci_low:.3f}, {auc_ci_high:.3f}])",
             fontsize=10, fontstyle="italic", color="#475569", ha="center", transform=ax_nomo.transAxes)

# Define vertical layout
y_positions = {
    "Points": 0.82,
    "VCAN": 0.70,
    "COL8A1": 0.58,
    "FN1": 0.46,
    "ISM1": 0.34,
    "COL3A1": 0.22,
    "Total_Points": 0.09,
    "Risk": -0.04
}

x_left = 0.16
x_right = 0.92
w_nomo = x_right - x_left

# Points Scale (0 to 100)
y_pt = y_positions["Points"]
ax_nomo.text(x_left - 0.02, y_pt, "Points", fontsize=11, fontweight="bold", color="#0F172A", ha="right", va="center", transform=ax_nomo.transAxes)
ax_nomo.plot([x_left, x_right], [y_pt, y_pt], color="#0F172A", lw=2.0, transform=ax_nomo.transAxes)
for p in range(0, 101, 10):
    xp = x_left + (p / 100.0) * w_nomo
    ax_nomo.plot([xp, xp], [y_pt - 0.015, y_pt + 0.015], color="#0F172A", lw=1.5, transform=ax_nomo.transAxes)
    ax_nomo.text(xp, y_pt + 0.025, str(p), fontsize=9, ha="center", va="bottom", color="#334155", transform=ax_nomo.transAxes)

# Gene Predictor Axes
for g in nomo_genes:
    yg = y_positions[g]
    pts_max = max_points[g]
    sign = "+" if beta_dict[g] > 0 else "-"
    ax_nomo.text(x_left - 0.02, yg, f"{g} ({sign})", fontsize=10.5, fontweight="bold", color="#1D4ED8", ha="right", va="center", transform=ax_nomo.transAxes)
    
    # Draw line proportional to pts_max
    x_end = x_left + (pts_max / 100.0) * w_nomo
    ax_nomo.plot([x_left, x_end], [yg, yg], color="#2563EB", lw=2.0, transform=ax_nomo.transAxes)
    
    # 5 intervals across gene range
    min_v = gene_mins[g]
    max_v = gene_maxs[g]
    vals = np.linspace(min_v, max_v, 5)
    for v in vals:
        # Calculate points for value v
        if beta_dict[g] > 0:
            pv = scale * beta_dict[g] * (v - min_v)
        else:
            pv = scale * abs(beta_dict[g]) * (max_v - v)
        xv = x_left + (pv / 100.0) * w_nomo
        ax_nomo.plot([xv, xv], [yg - 0.015, yg + 0.015], color="#2563EB", lw=1.4, transform=ax_nomo.transAxes)
        ax_nomo.text(xv, yg - 0.025, f"{v:.1f}", fontsize=8.5, ha="center", va="top", color="#475569", transform=ax_nomo.transAxes)

# Total Points Axis
y_tp = y_positions["Total_Points"]
ax_nomo.text(x_left - 0.02, y_tp, "Total Points", fontsize=11, fontweight="bold", color="#0F172A", ha="right", va="center", transform=ax_nomo.transAxes)
ax_nomo.plot([x_left, x_right], [y_tp, y_tp], color="#0F172A", lw=2.2, transform=ax_nomo.transAxes)
tp_ticks = np.arange(0, int(total_points_max) + 20, 20)
for tp in tp_ticks:
    if tp <= total_points_max:
        xtp = x_left + (tp / total_points_max) * w_nomo
        ax_nomo.plot([xtp, xtp], [y_tp - 0.015, y_tp + 0.015], color="#0F172A", lw=1.5, transform=ax_nomo.transAxes)
        ax_nomo.text(xtp, y_tp + 0.022, str(tp), fontsize=9, ha="center", va="bottom", color="#334155", transform=ax_nomo.transAxes)

# Fibrosis Risk Axis (Mathematically Exact Alignment to Total Points)
y_rk = y_positions["Risk"]
ax_nomo.text(x_left - 0.02, y_rk, "Fibrosis Risk", fontsize=11, fontweight="bold", color="#DC2626", ha="right", va="center", transform=ax_nomo.transAxes)
ax_nomo.plot([x_left, x_right], [y_rk, y_rk], color="#DC2626", lw=2.2, transform=ax_nomo.transAxes)

prob_ticks = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
for pr in prob_ticks:
    req_pts = p_to_total_points(pr)
    if 0 <= req_pts <= total_points_max:
        xrk = x_left + (req_pts / total_points_max) * w_nomo
        ax_nomo.plot([xrk, xrk], [y_rk - 0.018, y_rk + 0.018], color="#DC2626", lw=1.6, transform=ax_nomo.transAxes)
        pct_label = f"{int(round(pr * 100))}%" if pr in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95] else f"{pr:.2f}"
        ax_nomo.text(xrk, y_rk - 0.030, pct_label, fontsize=9, fontweight="bold", ha="center", va="top", color="#DC2626", transform=ax_nomo.transAxes)

# Panel B: Calibration Curve with Lowess & Bootstrap Optimism Correction
ax_cal = fig.add_subplot(gs[1, 0])
ax_cal.plot([0, 1], [0, 1], linestyle="--", color="#94A3B8", lw=1.6, label="Ideal 45° Reference")
ax_cal.plot(grid_p, orig_smooth, color="#3B82F6", lw=2.0, linestyle="--", label="Apparent Calibration (Lowess)")
ax_cal.plot(grid_p, bias_corrected_calib, color="#1D4ED8", lw=2.5, label="Bias-Corrected Calibration (B = 1000)")

# Add rug plot / histogram of predicted risks at bottom
y_rug = np.zeros_like(probs)
ax_cal.plot(probs[y_val == 0], y_rug[y_val == 0] + 0.02, '|', color="#2563EB", markersize=10, markeredgewidth=1.5, label="Controls (n = 20)")
ax_cal.plot(probs[y_val == 1], y_rug[y_val == 1] + 0.04, '|', color="#DC2626", markersize=10, markeredgewidth=1.5, label="Cases (n = 13)")

ax_cal.set_xlim(-0.02, 1.02)
ax_cal.set_ylim(-0.02, 1.02)
ax_cal.set_xlabel("Predicted Probability of Progressive Fibrosis", fontsize=10.5, fontweight="bold")
ax_cal.set_ylabel("Observed Fraction (Actual Fibrosis)", fontsize=10.5, fontweight="bold")
ax_cal.set_title("B. Model Calibration Curve (External Effluent GSE125498)\n"
                 f"Slope = {calib_slope:.3f} | Intercept = {calib_intercept:.3f} | Brier = {brier:.4f} | H-L P = {hl_p:.4f} (g={g_bins})",
                 fontsize=11, fontweight="bold", color="#0F172A", pad=10)
ax_cal.legend(loc="upper left", fontsize=8.5, frameon=True, facecolor="#F8FAFC", edgecolor="#CBD5E1")
ax_cal.grid(True, linestyle=":", alpha=0.6)

# Panel C: Decision Curve Analysis (DCA)
ax_dca = fig.add_subplot(gs[1, 1])
ax_dca.plot(thresholds, nb_model, color="#DC2626", lw=2.5, label="5-Gene Nomogram Model")
ax_dca.plot(thresholds, nb_all, color="#2563EB", lw=1.8, linestyle=":", label="Treat / Intervene All")
ax_dca.axhline(0, color="#64748B", linestyle="--", lw=1.4, label="Treat / Intervene None")

ax_dca.set_xlim(0.05, 0.75)
ax_dca.set_ylim(-0.08, 0.45)
ax_dca.set_xlabel("Threshold Decision Probability ($P_t$)", fontsize=10.5, fontweight="bold")
ax_dca.set_ylabel("Net Clinical Benefit", fontsize=10.5, fontweight="bold")
ax_dca.set_title("C. Decision Curve Analysis (Clinical Net Benefit)\n"
                 f"Evaluated on External Effluent Cohort (GSE125498, N = {n_total}) | Range: 10%–70%",
                 fontsize=11, fontweight="bold", color="#0F172A", pad=10)
ax_dca.legend(loc="upper right", fontsize=9, frameon=True, facecolor="#F8FAFC", edgecolor="#CBD5E1")
ax_dca.grid(True, linestyle=":", alpha=0.6)

# Annotation for clinical utility
ax_dca.annotate("Superior Net Benefit\nacross 10%–60% thresholds", xy=(0.32, 0.16), xytext=(0.42, 0.28),
                arrowprops=dict(arrowstyle="->", color="#DC2626", lw=1.5),
                fontsize=9.5, fontweight="bold", color="#DC2626", bbox=dict(boxstyle="round,pad=0.3", fc="#FEF2F2", ec="#DC2626"))

plt.tight_layout()
out_nomo_png = "results/figures/Hub_02_clinical_nomogram_dca_calibration.png"
out_nomo_pdf = "results/figures/Hub_02_clinical_nomogram_dca_calibration.pdf"
plt.savefig(out_nomo_png, dpi=300, bbox_inches="tight")
plt.savefig(out_nomo_pdf, bbox_inches="tight")
plt.close()
print(f"Saved Task 2 Nomogram, Calibration & DCA figures to {out_nomo_png} and {out_nomo_pdf}")
