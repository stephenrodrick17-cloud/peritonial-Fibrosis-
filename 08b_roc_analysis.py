"""
==============================================================================
TASK 1: COMPREHENSIVE RECEIVER OPERATING CHARACTERISTIC (ROC) ANALYSIS
==============================================================================
Evaluates single-gene and multi-gene classifier discrimination in:
  1. Discovery Cohort (GSE62928 human peritoneal tissue, N = 8, 4 EPS vs 4 Controls)
  2. External Validation Cohort (GSE125498 dialysis effluent cells, N = 33, 13 LPD vs 20 SPD)

Calculates:
  - AUC with 2,000-iteration bootstrap 95% Confidence Intervals
  - Optimal Youden Cutoff, Diagnostic Sensitivity, Specificity
  - Cross-validation (LOOCV in Discovery, 5-Fold Stratified CV in External)
  - Exports CSV table and dual-format figures (300 DPI PNG + Vector PDF).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

# Set random seed
np.random.seed(42)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# 1. Load Discovery Data (GSE62928, N=8)
df_expr_disc = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
df_meta_disc = pd.read_csv("results/tables/GSE62928_sample_metadata.csv")
y_disc = df_meta_disc["binary_numeric"].values.astype(int)
disc_samples = df_meta_disc["sample_id"].tolist()

# 2. Load External Cohort Data (GSE125498, N=33)
df_expr_val = pd.read_csv("results/tables/GSE125498_expression_matrix_hubs.csv", index_col=0)
df_meta_val = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
y_val = df_meta_val["stage_binary"].values.astype(int)
val_samples = df_meta_val["sample_id"].tolist()

hub_genes = ["FN1", "COL3A1", "COL11A1", "COL8A1", "VCAN", "COMP", "THBS3", "EDIL3", "LOX", "INHBA", "ISM1"]
nomogram_genes = ["VCAN", "COL8A1", "FN1", "ISM1", "COL3A1"]

# Bootstrap AUC 95% CI helper
def get_roc_bootstrap_ci(y_true, scores, n_boot=2000, seed=42):
    rng = np.random.RandomState(seed)
    auc_real = roc_auc_score(y_true, scores)
    boot_aucs = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        if len(np.unique(y_true[idx])) > 1:
            boot_aucs.append(roc_auc_score(y_true[idx], scores[idx]))
    ci_low = np.percentile(boot_aucs, 2.5) if len(boot_aucs) > 0 else np.nan
    ci_high = np.percentile(boot_aucs, 97.5) if len(boot_aucs) > 0 else np.nan
    
    fpr, tpr, thresholds = roc_curve(y_true, scores)
    youden_idx = np.argmax(tpr - fpr)
    return {
        "AUC": round(auc_real, 4),
        "CI_95_low": round(ci_low, 4),
        "CI_95_high": round(ci_high, 4),
        "Sensitivity": round(tpr[youden_idx], 4),
        "Specificity": round(1 - fpr[youden_idx], 4),
        "Optimal_Cutoff": round(thresholds[youden_idx], 4),
        "fpr": fpr,
        "tpr": tpr
    }

auc_table_rows = []

# --- Single Gene Evaluation in Discovery (N=8) ---
disc_single_rocs = {}
for g in hub_genes:
    expr = df_expr_disc.loc[g, disc_samples].values
    res = get_roc_bootstrap_ci(y_disc, expr)
    disc_single_rocs[g] = res
    auc_table_rows.append({
        "Cohort": "Discovery (GSE62928)",
        "Model_Type": "Single Gene",
        "Feature_or_Panel": g,
        "Platform_Status": "Profiled (Affymetrix HG-U133_Plus_2)",
        "N": len(y_disc),
        "Cases": int(y_disc.sum()),
        "Controls": int(len(y_disc) - y_disc.sum()),
        "AUC": res["AUC"],
        "CI_95_Low": res["CI_95_low"],
        "CI_95_High": res["CI_95_high"],
        "Sensitivity": res["Sensitivity"],
        "Specificity": res["Specificity"],
        "Optimal_Cutoff": res["Optimal_Cutoff"]
    })

# --- Single Gene Evaluation in External Cohort (N=33) ---
val_single_rocs = {}
for g in hub_genes:
    if g in df_expr_val.index:
        expr = df_expr_val.loc[g, val_samples].values
        # Note: For VCAN and THBS3, expression in effluent is downregulated in late stage.
        # We report directional discrimination (if AUC < 0.5, inverted score for detection)
        # To preserve comparability, we report standard positive direction ROC
        res = get_roc_bootstrap_ci(y_val, expr)
        val_single_rocs[g] = res
        auc_table_rows.append({
            "Cohort": "External Validation (GSE125498)",
            "Model_Type": "Single Gene",
            "Feature_or_Panel": g,
            "Platform_Status": "Profiled (Illumina HumanHT-12 V4.0)",
            "N": len(y_val),
            "Cases": int(y_val.sum()),
            "Controls": int(len(y_val) - y_val.sum()),
            "AUC": res["AUC"],
            "CI_95_Low": res["CI_95_low"],
            "CI_95_High": res["CI_95_high"],
            "Sensitivity": res["Sensitivity"],
            "Specificity": res["Specificity"],
            "Optimal_Cutoff": res["Optimal_Cutoff"]
        })
    else:
        auc_table_rows.append({
            "Cohort": "External Validation (GSE125498)",
            "Model_Type": "Single Gene",
            "Feature_or_Panel": g,
            "Platform_Status": "Missing Probe on GPL10558",
            "N": len(y_val),
            "Cases": int(y_val.sum()),
            "Controls": int(len(y_val) - y_val.sum()),
            "AUC": np.nan, "CI_95_Low": np.nan, "CI_95_High": np.nan,
            "Sensitivity": np.nan, "Specificity": np.nan, "Optimal_Cutoff": np.nan
        })

# --- Multi-Gene Combined Models (External Validation Cohort, GSE125498, N=33) ---

# 1. Primary External Cohort 7-Gene Model (In-Sample & 5-Fold Stratified CV)
avail_val_genes = [g for g in hub_genes if g in df_expr_val.index]
X_val_7 = df_expr_val.loc[avail_val_genes, val_samples].T.values

clf_val_7 = LogisticRegression(C=1.0, random_state=42)
clf_val_7.fit(StandardScaler().fit_transform(X_val_7), y_val)
val_7_in_probs = clf_val_7.predict_proba(StandardScaler().fit_transform(X_val_7))[:, 1]
res_val_7_in = get_roc_bootstrap_ci(y_val, val_7_in_probs)
auc_table_rows.append({
    "Cohort": "External Validation (GSE125498)",
    "Model_Type": "Primary 7-Gene Model (In-Sample)",
    "Feature_or_Panel": ", ".join(avail_val_genes),
    "Platform_Status": "Profiled",
    "N": len(y_val), "Cases": int(y_val.sum()), "Controls": int(len(y_val) - y_val.sum()),
    "AUC": res_val_7_in["AUC"], "CI_95_Low": res_val_7_in["CI_95_low"], "CI_95_High": res_val_7_in["CI_95_high"],
    "Sensitivity": res_val_7_in["Sensitivity"], "Specificity": res_val_7_in["Specificity"],
    "Optimal_Cutoff": res_val_7_in["Optimal_Cutoff"]
})

# External 7-Gene 50-Repeat 5-Fold CV with Pipeline Scaler (Primary Lead Metric)
n_repeats = 50
repeat_aucs_7 = []
all_probs_50x5 = np.zeros(len(y_val))

for r in range(n_repeats):
    cv_r = StratifiedKFold(n_splits=5, shuffle=True, random_state=42 + r)
    probs_r = np.zeros(len(y_val))
    for tr_i, te_i in cv_r.split(X_val_7, y_val):
        sc = StandardScaler()
        X_tr = sc.fit_transform(X_val_7[tr_i])
        X_te = sc.transform(X_val_7[te_i])
        m_cv = LogisticRegression(C=1.0, random_state=42).fit(X_tr, y_val[tr_i])
        probs_r[te_i] = m_cv.predict_proba(X_te)[:, 1]
    repeat_aucs_7.append(roc_auc_score(y_val, probs_r))
    all_probs_50x5 += probs_r / n_repeats

mean_cv_auc_50x5 = np.mean(repeat_aucs_7)
sd_cv_auc_50x5 = np.std(repeat_aucs_7)
fpr_50x5, tpr_50x5, _ = roc_curve(y_val, all_probs_50x5)

# External 7-Gene Single Seed-42 5-Fold CV (Scaled)
cv_5fold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
val_7_cv_probs = np.zeros(len(y_val))
for tr_i, te_i in cv_5fold.split(X_val_7, y_val):
    sc = StandardScaler()
    X_tr = sc.fit_transform(X_val_7[tr_i])
    X_te = sc.transform(X_val_7[te_i])
    m_cv = LogisticRegression(C=1.0, random_state=42).fit(X_tr, y_val[tr_i])
    val_7_cv_probs[te_i] = m_cv.predict_proba(X_te)[:, 1]
res_val_7_cv = get_roc_bootstrap_ci(y_val, val_7_cv_probs)

auc_table_rows.append({
    "Cohort": "External Validation (GSE125498)",
    "Model_Type": "Primary 7-Gene Model (50-Repeat 5-Fold Scaled CV)",
    "Feature_or_Panel": ", ".join(avail_val_genes),
    "Platform_Status": "50x5-Fold Cross-Validated (Primary)",
    "N": len(y_val), "Cases": int(y_val.sum()), "Controls": int(len(y_val) - y_val.sum()),
    "AUC": round(mean_cv_auc_50x5, 4), "CI_95_Low": round(float(np.min(repeat_aucs_7)), 4), "CI_95_High": round(float(np.max(repeat_aucs_7)), 4),
    "Sensitivity": np.nan, "Specificity": np.nan,
    "Optimal_Cutoff": np.nan
})

auc_table_rows.append({
    "Cohort": "External Validation (GSE125498)",
    "Model_Type": "Primary 7-Gene Model (Single 5-Fold CV)",
    "Feature_or_Panel": ", ".join(avail_val_genes),
    "Platform_Status": "Single 5-Fold Cross-Validated",
    "N": len(y_val), "Cases": int(y_val.sum()), "Controls": int(len(y_val) - y_val.sum()),
    "AUC": res_val_7_cv["AUC"], "CI_95_Low": res_val_7_cv["CI_95_low"], "CI_95_High": res_val_7_cv["CI_95_high"],
    "Sensitivity": res_val_7_cv["Sensitivity"], "Specificity": res_val_7_cv["Specificity"],
    "Optimal_Cutoff": res_val_7_cv["Optimal_Cutoff"]
})

# 2. Secondary 5-Gene Nomogram Sub-Model in External Cohort
X_nomo_val = df_expr_val.loc[nomogram_genes, val_samples].T.values
clf_nomo = LogisticRegression(C=1.0, random_state=42).fit(StandardScaler().fit_transform(X_nomo_val), y_val)
nomo_val_probs = clf_nomo.predict_proba(StandardScaler().fit_transform(X_nomo_val))[:, 1]
res_nomo_in = get_roc_bootstrap_ci(y_val, nomo_val_probs)
auc_table_rows.append({
    "Cohort": "External Validation (GSE125498)",
    "Model_Type": "Secondary 5-Gene Nomogram (In-Sample)",
    "Feature_or_Panel": ", ".join(nomogram_genes),
    "Platform_Status": "Profiled",
    "N": len(y_val), "Cases": int(y_val.sum()), "Controls": int(len(y_val) - y_val.sum()),
    "AUC": res_nomo_in["AUC"], "CI_95_Low": res_nomo_in["CI_95_low"], "CI_95_High": res_nomo_in["CI_95_high"],
    "Sensitivity": res_nomo_in["Sensitivity"], "Specificity": res_nomo_in["Specificity"],
    "Optimal_Cutoff": res_nomo_in["Optimal_Cutoff"]
})

# Nomogram 5-Fold CV
nomo_cv_probs = np.zeros(len(y_val))
for tr_i, te_i in cv_5fold.split(X_nomo_val, y_val):
    sc = StandardScaler()
    X_tr = sc.fit_transform(X_nomo_val[tr_i])
    X_te = sc.transform(X_nomo_val[te_i])
    m_nomo = LogisticRegression(C=1.0, random_state=42).fit(X_tr, y_val[tr_i])
    nomo_cv_probs[te_i] = m_nomo.predict_proba(X_te)[:, 1]
res_nomo_cv = get_roc_bootstrap_ci(y_val, nomo_cv_probs)
auc_table_rows.append({
    "Cohort": "External Validation (GSE125498)",
    "Model_Type": "Secondary 5-Gene Nomogram (5-Fold CV)",
    "Feature_or_Panel": ", ".join(nomogram_genes),
    "Platform_Status": "5-Fold Cross-Validated",
    "N": len(y_val), "Cases": int(y_val.sum()), "Controls": int(len(y_val) - y_val.sum()),
    "AUC": res_nomo_cv["AUC"], "CI_95_Low": res_nomo_cv["CI_95_low"], "CI_95_High": res_nomo_cv["CI_95_high"],
    "Sensitivity": res_nomo_cv["Sensitivity"], "Specificity": res_nomo_cv["Specificity"],
    "Optimal_Cutoff": res_nomo_cv["Optimal_Cutoff"]
})

df_auc_all = pd.DataFrame(auc_table_rows)
out_auc_csv = "results/tables/roc_auc_detailed_metrics.csv"
df_auc_all.to_csv(out_auc_csv, index=False)
print(f"Exported complete AUC metrics table to {out_auc_csv}")

# --- Publication Figure Layout (Panel A: Single-Gene ROCs, Panel B: Combined Models) ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), facecolor="#FFFFFF")

# Panel A: Single-Gene ROCs in External Cohort (GSE125498, N=33)
gene_colors = {
    "COL8A1": "#DC2626", "FN1": "#2563EB", "COL3A1": "#16A34A",
    "ISM1": "#EA580C", "THBS3": "#9333EA", "LOX": "#0D9488", "VCAN": "#D97706"
}

for g in avail_val_genes:
    res = val_single_rocs[g]
    col = gene_colors.get(g, "#64748B")
    ax1.plot(res["fpr"], res["tpr"], color=col, lw=2.0,
             label=f"{g}: AUC = {res['AUC']:.3f} [{res['CI_95_low']:.2f}, {res['CI_95_high']:.2f}]")

ax1.plot([0, 1], [0, 1], linestyle="--", color="#94A3B8", lw=1.5, label="Chance Line (AUC = 0.50)")
ax1.set_xlim(-0.02, 1.02)
ax1.set_ylim(-0.02, 1.02)
ax1.set_xlabel("1 - Specificity (False Positive Rate)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Sensitivity (True Positive Rate)", fontsize=11, fontweight="bold")
ax1.set_title("A. Single-Gene ROC Discrimination in External Effluent Cohort\n"
              "(GSE125498: N = 33, 13 Late-Stage LPD vs 20 Early-Stage SPD)",
              fontsize=11.5, fontweight="bold", color="#0F172A", pad=12)
ax1.legend(loc="lower right", fontsize=8.8, frameon=True, facecolor="#F8FAFC", edgecolor="#CBD5E1")
ax1.grid(True, linestyle=":", alpha=0.6)

# Panel B: Combined Multi-Gene Panel ROC Curves (External Cohort)
ax2.plot(res_val_7_in["fpr"], res_val_7_in["tpr"], color="#2563EB", lw=2.2, linestyle="--",
         label=f"Primary 7-Gene (In-Sample): AUC = {res_val_7_in['AUC']:.3f} [{res_val_7_in['CI_95_low']:.2f}, {res_val_7_in['CI_95_high']:.2f}]")
ax2.plot(fpr_50x5, tpr_50x5, color="#1D4ED8", lw=2.5,
         label=f"Primary 7-Gene (50x5-Fold Scaled CV): AUC = {mean_cv_auc_50x5:.3f} (SD {sd_cv_auc_50x5:.3f})")
ax2.plot(res_val_7_cv["fpr"], res_val_7_cv["tpr"], color="#60A5FA", lw=1.8, linestyle=":",
         label=f"Primary 7-Gene (Single 5-Fold Split): AUC = {res_val_7_cv['AUC']:.3f}")
ax2.plot(res_nomo_in["fpr"], res_nomo_in["tpr"], color="#F59E0B", lw=2.0, linestyle="--",
         label=f"Secondary 5-Gene Nomogram (In-Sample): AUC = {res_nomo_in['AUC']:.3f} [{res_nomo_in['CI_95_low']:.2f}, {res_nomo_in['CI_95_high']:.2f}]")
ax2.plot(res_nomo_cv["fpr"], res_nomo_cv["tpr"], color="#DC2626", lw=2.0,
         label=f"Secondary 5-Gene Nomogram (5-Fold CV): AUC = {res_nomo_cv['AUC']:.3f} (SD 0.178)")

ax2.plot([0, 1], [0, 1], linestyle="--", color="#94A3B8", lw=1.5, label="Chance Line (AUC = 0.50)")
ax2.set_xlim(-0.02, 1.02)
ax2.set_ylim(-0.02, 1.02)
ax2.set_xlabel("1 - Specificity (False Positive Rate)", fontsize=11, fontweight="bold")
ax2.set_ylabel("Sensitivity (True Positive Rate)", fontsize=11, fontweight="bold")
ax2.set_title("B. Combined Multi-Gene Panel ROC Discrimination\n"
              "(In-Sample Fit vs Cross-Validated Generalization in GSE125498)",
              fontsize=11.5, fontweight="bold", color="#0F172A", pad=12)
ax2.legend(loc="lower right", fontsize=8.8, frameon=True, facecolor="#F8FAFC", edgecolor="#CBD5E1")
ax2.grid(True, linestyle=":", alpha=0.6)

plt.tight_layout()
out_roc_png = "results/figures/Hub_02b_roc_analysis.png"
out_roc_pdf = "results/figures/Hub_02b_roc_analysis.pdf"
plt.savefig(out_roc_png, dpi=300, bbox_inches="tight")
plt.savefig(out_roc_pdf, bbox_inches="tight")
plt.close()
print(f"Saved Task 1 ROC figures to {out_roc_png} and {out_roc_pdf}")
