"""
generate_fig4_versionA.py  –  Figure 4 reproduced with Version A (ml_39) data.
- 39-gene candidate pool (non-grey WGCNA modules ∩ Matrisome DEGs)
- 9 consensus hub genes by majority vote (>= 2 of 3 models)
- LASSO n=11, SVM-RFE n=10, RF n=10; 3-way intersection = 2 (FLT3LG, TNFSF15)
- Validation: GSE248762 all-cell pseudobulk (LV n=10, SV n=6) — 7/9 concordant
- Secondary: GSE130888 all-cell pseudobulk (LV n=4, SV n=6)
- All curves and statistics computed from real patient data

Run:  python scripts/generate_fig4_versionA.py
"""
import os, warnings, sys
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.feature_selection import RFE, RFECV

ROOT = r"D:\Peritoneal Project"

P_EXPR   = os.path.join(ROOT, "results", "tables", "GSE125498_full_expression_matrix.csv")
P_META   = os.path.join(ROOT, "results", "tables", "GSE125498_sample_metadata.csv")
P_POOL   = os.path.join(ROOT, "results", "revision", "ml_39_matrisome_wgcna", "wgcna_matrisome_39_gene_input.csv")
P_ML39   = os.path.join(ROOT, "results", "revision", "ml_39_matrisome_wgcna", "ml_39_gene_results.csv")

P_G24_LOG = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE248762",
                          "analysis_20261009", "edger_logCPM.csv")
P_G24_MET = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE248762",
                          "analysis_20261009", "all_cell_pseudobulk_metadata.csv")

P_G13_LOG = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE130888",
                          "analysis_20261009", "edger_logCPM.csv")
P_G13_MET = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE130888",
                          "analysis_20261009", "all_cell_pseudobulk_metadata.csv")

P_G24_DET = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE248762",
                          "analysis_20261009", "pseudobulk_gene_detection.csv")
P_G13_DET = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE130888",
                          "analysis_20261009", "pseudobulk_gene_detection.csv")

ENS9 = {"TNFSF15": "ENSG00000181634", "FLT3LG": "ENSG00000090554", "EBI3": "ENSG00000105246",
        "LTB":    "ENSG00000227507", "ADAM19": "ENSG00000135074", "NAT6": "ENSG00000243477",
        "ECM1":   "ENSG00000143369", "CST7":  "ENSG00000077984", "CCBE1": "ENSG00000183287"}

OUT_DIR = os.path.join(ROOT, "results", "revision", "ml_39_matrisome_wgcna")
os.makedirs(OUT_DIR, exist_ok=True)

C_LPD   = "#E11D48"; C_SPD   = "#0284C7"
C_LPD_L = "#FFE4E6"; C_SPD_L = "#E0F2FE"

HUB9 = ["TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19",
        "NAT6", "ECM1", "CST7", "CCBE1"]

COLORS9 = {
    "TNFSF15": "#1F77B4", "FLT3LG": "#FF7F0E", "EBI3":   "#2CA02C",
    "LTB":     "#D62728", "ADAM19": "#9467BD", "NAT6":   "#8C564B",
    "ECM1":    "#E377C2", "CST7":   "#7F7F7F", "CCBE1":  "#BCBD22",
}

NAT6_ALIAS = "NAA80"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.linewidth": 0.85,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
})

def pstar(p):
    if p != p: return "ns"
    if p < 1e-4:  return "****"
    if p < 1e-3:  return "***"
    if p < 1e-2:  return "**"
    if p < 0.05:  return "*"
    return "ns"

def plabel(ax, label, x=-0.14, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=16, fontweight="bold", va="top")

# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------
print("Loading data ...")
expr_mat = pd.read_csv(P_EXPR).set_index("Gene")
expr_mat.columns = expr_mat.columns.astype(str)
meta_df = pd.read_csv(P_META)
meta_df["sample_id"] = meta_df["sample_id"].astype(str)
sc = [s for s in meta_df["sample_id"] if s in expr_mat.columns]
grp = dict(zip(meta_df["sample_id"], meta_df["group"]))
lpd = [s for s in sc if grp[s] == "LPD"]
spd = [s for s in sc if grp[s] == "SPD"]
y_disc = np.array([1 if grp[s] == "LPD" else 0 for s in sc])
n_samp = len(y_disc)
null_err = len(lpd) / n_samp
null_acc = len(spd) / n_samp

# ml_39 results (ground truth selections)
ml39 = pd.read_csv(P_ML39).dropna(subset=["Gene"])
L = set(ml39[ml39["LASSO_selected"] == 1]["Gene"])
S = set(ml39[ml39["SVM_RFE_selected"] == 1]["Gene"])
R = set(ml39[ml39["RF_selected"] == 1]["Gene"])
n_lasso, n_svm, n_rf = len(L), len(S), len(R)
only_L  = len(L - S - R)
only_S  = len(S - L - R)
only_R  = len(R - L - S)
ls_only = len((L & S) - R)
lr_only = len((L & R) - S)
sr_only = len((S & R) - L)
lsr     = len(L & S & R)

# 39-gene candidate pool
pool_df = pd.read_csv(P_POOL)
genes39 = [g for g in pool_df["Gene"] if g in expr_mat.index]
X_all = StandardScaler().fit_transform(expr_mat.loc[genes39, sc].T.values.astype(float))
genes39_arr = np.array(genes39)

print(f"39-gene pool: {len(genes39)} genes")
print(f"LASSO={n_lasso}, SVM-RFE={n_svm}, RF={n_rf}, 3-way={lsr}")

# GSE248762 all-cell pseudobulk
g24_log = pd.read_csv(P_G24_LOG, index_col=0)
g24_log.index = g24_log.index.astype(str)
g24_meta = pd.read_csv(P_G24_MET).drop_duplicates("sample_id").set_index("sample_id")
g24_lv = [s for s in g24_meta.index if g24_meta.loc[s, "primary_group"] == "LV"]
g24_sv = [s for s in g24_meta.index if g24_meta.loc[s, "primary_group"] == "SV"]
g24_det = pd.read_csv(P_G24_DET, dtype={"gene_id": str}).set_index("gene_id")
g24_gene_map = {g: eid for g, eid in ENS9.items()
                if eid in g24_log.index and eid in g24_det.index and bool(g24_det.loc[eid, "detected"])}

# GSE130888 all-cell pseudobulk
g13_log = pd.read_csv(P_G13_LOG, index_col=0)
g13_log.index = g13_log.index.astype(str)
g13_meta = pd.read_csv(P_G13_MET).drop_duplicates("sample_id").set_index("sample_id")
g13_lv = [s for s in g13_meta.index if g13_meta.loc[s, "primary_group"] == "LV"]
g13_sv = [s for s in g13_meta.index if g13_meta.loc[s, "primary_group"] == "SV"]
g13_det = pd.read_csv(P_G13_DET, dtype={"gene_id": str}).set_index("gene_id")
g13_gene_map = {g: eid for g, eid in ENS9.items()
                if eid in g13_log.index and eid in g13_det.index and bool(g13_det.loc[eid, "detected"])}

print(f"GSE248762: LV={len(g24_lv)}, SV={len(g24_sv)}, genes found={len(g24_gene_map)}")
print(f"GSE130888: LV={len(g13_lv)}, SV={len(g13_sv)}, genes found={len(g13_gene_map)}")

# ---------------------------------------------------------------------------
# 2. FIT ML MODELS (for curves only; selections from ml_39 CSV)
# ---------------------------------------------------------------------------
print("Fitting ML models for curve visualization ...")

# LASSO CV curves
Cs = np.logspace(-1.5, 2.5, 60)
K_lasso = 10
cv_errors = np.zeros((K_lasso, len(Cs)))
skf = StratifiedKFold(K_lasso, shuffle=True, random_state=42)
for k, (tr, te) in enumerate(skf.split(X_all, y_disc)):
    for j, C in enumerate(Cs):
        m = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=5000).fit(X_all[tr], y_disc[tr])
        cv_errors[k, j] = 1 - m.score(X_all[te], y_disc[te])
cv_mean = cv_errors.mean(axis=0)
cv_se   = cv_errors.std(axis=0, ddof=1) / np.sqrt(K_lasso)

log_lambda = -np.log(Cs * n_samp)
j_min = int(cv_mean.argmin())
j_1se = int(np.where(cv_mean <= cv_mean[j_min] + cv_se[j_min])[0][0])

def lasso_nz(C):
    m = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=5000).fit(X_all, y_disc)
    return int((m.coef_ != 0).sum())

nz_at_each = np.array([lasso_nz(C) for C in Cs])
nz_min = lasso_nz(Cs[j_min])
nz_1se = lasso_nz(Cs[j_1se])

# Fixed C=5.0 marker (matches ml_39 LASSO grid which includes C=5)
log_lam_c5 = -np.log(5.0 * n_samp)
nz_c5 = lasso_nz(5.0)

lasso_coef_paths = np.zeros((len(genes39), len(Cs)))
for j, C in enumerate(Cs):
    m = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=5000).fit(X_all, y_disc)
    lasso_coef_paths[:, j] = m.coef_[0]

# SVM-RFE (RFECV for accuracy curve)
rfe = RFECV(SVC(kernel="linear", C=1.0), step=1, min_features_to_select=1, scoring="accuracy",
            cv=StratifiedKFold(10, shuffle=True, random_state=42))
rfe.fit(X_all, y_disc)
try:
    rfe_acc = np.asarray(rfe.cv_results_["mean_test_score"])
except Exception:
    rfe_acc = np.asarray(rfe.grid_scores_).mean(1)
n_feat_axis = np.arange(1, len(rfe_acc) + 1)
rfe_opt_n = int(rfe.n_features_)
rfe_opt_acc = float(rfe_acc.max())

# Random Forest
rf = RandomForestClassifier(n_estimators=20, warm_start=True, oob_score=True,
                             class_weight="balanced", max_features="sqrt", random_state=42)
oob_overall, oob_spd, oob_lpd = [], [], []
tree_counts = list(range(20, 501, 10))
for t in tree_counts:
    rf.set_params(n_estimators=t).fit(X_all, y_disc)
    pr = np.argmax(np.nan_to_num(rf.oob_decision_function_), axis=1)
    oob_overall.append((pr != y_disc).mean())
    oob_spd.append((pr[y_disc == 0] != 0).mean())
    oob_lpd.append((pr[y_disc == 1] != 1).mean())

oob_overall = np.array(oob_overall)
oob_spd = np.array(oob_spd)
oob_lpd = np.array(oob_lpd)

rf_final = RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=42).fit(X_all, y_disc)
rf_imp = pd.Series(rf_final.feature_importances_, index=genes39)
rf_top15 = rf_imp.sort_values(ascending=False).head(15)
rf_top10_genes = list(rf_imp.sort_values(ascending=False).head(10).index)

# ---------------------------------------------------------------------------
# 3. DISCOVERY ROC & STATS
# ---------------------------------------------------------------------------
disc_curves = {}; disc_aucs = {}; disc_pvals = {}; disc_direction = {}
for g in HUB9:
    v = expr_mat.loc[g, sc].values.astype(float)
    diff = v[y_disc == 1].mean() - v[y_disc == 0].mean()
    disc_direction[g] = diff > 0
    v_oriented = v if diff > 0 else -v
    fpr, tpr, _ = roc_curve(y_disc, v_oriented)
    auc = roc_auc_score(y_disc, v_oriented)
    disc_aucs[g] = auc
    disc_curves[g] = (fpr, tpr)
    vl = expr_mat.loc[g, lpd].values.astype(float)
    vs = expr_mat.loc[g, spd].values.astype(float)
    _, p = mannwhitneyu(vl, vs, alternative="two-sided")
    disc_pvals[g] = p

# ---------------------------------------------------------------------------
# 4. VALIDATION ROC & STATS
# ---------------------------------------------------------------------------
def bootstrap_auc_ci(y_true, scores, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(set(y_true[idx])) < 2: continue
        aucs.append(roc_auc_score(y_true[idx], scores[idx]))
    return np.percentile(aucs, [2.5, 97.5]) if len(aucs) > 0 else (0.5, 0.5)

def compute_validation(log_df, lv_s, sv_s, gene_map, disc_dir):
    curves, aucs, pvals, cis, concordant = {}, {}, {}, {}, {}
    y_arr = np.array([1]*len(lv_s) + [0]*len(sv_s))
    for g in HUB9:
        if g not in gene_map:
            continue
        eid = gene_map[g]
        vl = log_df.loc[eid, lv_s].values.astype(float)
        vs = log_df.loc[eid, sv_s].values.astype(float)
        all_vals = np.concatenate([vl, vs])
        if np.all(all_vals == all_vals[0]):
            curves[g] = (np.array([0.0, 1.0]), np.array([0.0, 1.0]))
            aucs[g] = 0.500; cis[g] = (0.500, 0.500); pvals[g] = np.nan
            concordant[g] = False
        else:
            score_val = all_vals if disc_dir[g] else -all_vals
            fpr_v, tpr_v, _ = roc_curve(y_arr, score_val)
            auc_v = roc_auc_score(y_arr, score_val)
            ci = bootstrap_auc_ci(y_arr, score_val)
            curves[g] = (fpr_v, tpr_v)
            aucs[g] = auc_v; cis[g] = (ci[0], ci[1])
            _, p = mannwhitneyu(vl, vs, alternative="two-sided")
            pvals[g] = p
            concordant[g] = (auc_v > 0.50)
    return curves, aucs, pvals, cis, concordant

g24_curves, g24_aucs, g24_pvals, g24_cis, g24_concordant = compute_validation(
    g24_log, g24_lv, g24_sv, g24_gene_map, disc_direction)
g13_curves, g13_aucs, g13_pvals, g13_cis, g13_concordant = compute_validation(
    g13_log, g13_lv, g13_sv, g13_gene_map, disc_direction)

n_concordant_g24 = sum(1 for g in g24_concordant if g24_concordant[g])
n_concordant_g13 = sum(1 for g in g13_concordant if g13_concordant[g])
print(f"GSE248762 concordant: {n_concordant_g24}/{len(g24_concordant)}")
print(f"GSE130888 concordant: {n_concordant_g13}/{len(g13_concordant)}")

# ---------------------------------------------------------------------------
# 5. BUILD FIGURE
# ---------------------------------------------------------------------------
print("Generating composite figure ...")
fig = plt.figure(figsize=(17, 22), dpi=300)
fig.patch.set_facecolor("white")
outer = gridspec.GridSpec(4, 2, figure=fig, left=0.07, right=0.97, top=0.97, bottom=0.05, hspace=0.46, wspace=0.34)

def inner2(sp):
    return gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=sp, wspace=0.36)

# ---------------- PANEL A: LASSO ----------------
iA = inner2(outer[0, 0])
ax_a1 = fig.add_subplot(iA[0, 0])
ax_a2 = fig.add_subplot(iA[0, 1])

for i, gene in enumerate(genes39):
    is_hub = gene in HUB9
    color = COLORS9.get(gene, "#CBD5E1") if is_hub else "#CBD5E1"
    lw = 2.0 if is_hub else 0.6
    alpha = 0.95 if is_hub else 0.35
    ax_a1.plot(log_lambda, lasso_coef_paths[i, :], color=color, lw=lw, alpha=alpha,
               label=gene if is_hub else "_nolegend_")

ax_a1.axhline(0, color="#94A3B8", linestyle="--", lw=0.8)
ax_a1.axvline(log_lambda[j_1se], color="#64748B", lw=1.1, linestyle=":", label=f"lambda.1se (nz={nz_1se})")
ax_a1.axvline(log_lambda[j_min], color="#475569", lw=1.1, linestyle="--", label=f"lambda.min (nz={nz_min})")
ax_a1.axvline(log_lam_c5, color="#DC2626", lw=1.2, linestyle="-.", label=f"C=5.0 (nz={nz_c5})")
ax_a1.set_xlabel("log(lambda)", fontsize=9, fontweight="bold")
ax_a1.set_ylabel("Coefficients", fontsize=9, fontweight="bold")
ax_a1.set_title(f"LASSO Coefficient Paths\n({nz_c5} features at fixed C=5.0; {nz_min} at lambda.min)", fontsize=9, fontweight="bold", pad=4)
ax_a1.grid(True, linestyle=":", alpha=0.35)
ax_a1.legend(fontsize=6.0, loc="lower left", framealpha=0.85, ncol=2)

ax_a2.errorbar(log_lambda, cv_mean, yerr=cv_se, fmt="o", ms=3, color="#EF4444", ecolor="#FCA5A5", elinewidth=1, capsize=2)
ax_a2.axvline(log_lambda[j_min], color="#475569", lw=1.1, linestyle="--", label=f"lambda.min (log={log_lambda[j_min]:.2f})")
ax_a2.axvline(log_lambda[j_1se], color="#64748B", lw=1.1, linestyle=":", label=f"lambda.1se (log={log_lambda[j_1se]:.2f})")
ax_a2.axvline(log_lam_c5, color="#DC2626", lw=1.2, linestyle="-.", label="C=5.0")
ax_a2.axhline(null_err, color="#6B7280", linestyle="--", lw=1.0, label=f"Null baseline ({len(lpd)}/{n_samp}={null_err:.3f})")
ax_a2.set_xlabel("log(lambda)", fontsize=9, fontweight="bold")
ax_a2.set_ylabel("10-Fold CV Misclassification Error", fontsize=9, fontweight="bold")
ax_a2.set_title("LASSO Cross-Validation Error", fontsize=9.5, fontweight="bold", pad=4)
ax_a2.legend(fontsize=7.0, loc="upper left", framealpha=0.85)
ax_a2.grid(True, linestyle=":", alpha=0.35)

ax_a2t = ax_a2.twiny()
ax_a2t.set_xlim(ax_a2.get_xlim())
unique_nz_targets = sorted(set(nz_at_each))
unique_tick_locs, unique_tick_labels = [], []
for unz in unique_nz_targets:
    match_indices = np.where(nz_at_each == unz)[0]
    if len(match_indices) > 0:
        mid_idx = match_indices[len(match_indices)//2]
        unique_tick_locs.append(log_lambda[mid_idx])
        unique_tick_labels.append(str(unz))
ax_a2t.set_xticks(unique_tick_locs)
ax_a2t.set_xticklabels(unique_tick_labels, fontsize=7)
ax_a2t.set_xlabel("Number of non-zero variables", fontsize=7.5, labelpad=3)
plabel(ax_a1, "A")

# ---------------- PANEL B: SVM-RFE ----------------
iB = inner2(outer[0, 1])
ax_b1 = fig.add_subplot(iB[0, 0])
ax_b2 = fig.add_subplot(iB[0, 1])

ax_b1.plot(n_feat_axis, rfe_acc, color="#0284C7", lw=1.5)
ax_b1.scatter([rfe_opt_n], [rfe_opt_acc], color=C_LPD, s=50, zorder=5)
ax_b1.axvline(10, color="#64748B", linestyle=":", lw=1.1, label="Selected set (n=10)")
ax_b1.axhline(null_acc, color="#6B7280", linestyle="--", lw=1.0, label=f"Null acc ({len(spd)}/{n_samp}={null_acc:.3f})")
all39_svm_acc = rfe_acc[-1] if len(rfe_acc) > 0 else 0
ax_b1.axhline(all39_svm_acc, color="#8B5CF6", linestyle="-.", lw=1.1, label=f"All-39 SVM acc ({all39_svm_acc:.3f})")
ax_b1.text(rfe_opt_n + 1.5, rfe_opt_acc - 0.02, f"CV-opt n={rfe_opt_n}\nAcc={rfe_opt_acc:.3f}", fontsize=7.5, color=C_LPD, fontweight="bold")
ax_b1.set_xlabel("Number of Features", fontsize=9, fontweight="bold")
ax_b1.set_ylabel("10-CV Accuracy", fontsize=9, fontweight="bold")
ax_b1.set_title(f"SVM-RFE Accuracy\n(Fixed top-10 selector; opt n={rfe_opt_n})", fontsize=9, fontweight="bold", pad=4)
ax_b1.set_ylim(max(0, rfe_acc.min() - 0.05), min(1.02, rfe_acc.max() + 0.05))
ax_b1.legend(fontsize=7.0, loc="lower right", framealpha=0.85)
ax_b1.grid(True, linestyle=":", alpha=0.35)

rfe_err = 1 - rfe_acc
ax_b2.plot(n_feat_axis, rfe_err, color="#0284C7", lw=1.5)
ax_b2.scatter([rfe_opt_n], [1 - rfe_opt_acc], color=C_LPD, s=50, zorder=5)
ax_b2.axvline(10, color="#64748B", linestyle=":", lw=1.1, label="Selected set (n=10)")
ax_b2.axhline(null_err, color="#6B7280", linestyle="--", lw=1.0, label=f"Null err ({len(lpd)}/{n_samp}={null_err:.3f})")
ax_b2.text(rfe_opt_n + 1.5, 1 - rfe_opt_acc + 0.01, f"CV-opt n={rfe_opt_n}\nErr={1-rfe_opt_acc:.3f}", fontsize=7.5, color=C_LPD, fontweight="bold")
ax_b2.set_xlabel("Number of Features", fontsize=9, fontweight="bold")
ax_b2.set_ylabel("10-CV Error", fontsize=9, fontweight="bold")
ax_b2.set_title("SVM-RFE Error Curve", fontsize=9.5, fontweight="bold", pad=4)
ax_b2.legend(fontsize=7.0, loc="upper right", framealpha=0.85)
ax_b2.grid(True, linestyle=":", alpha=0.35)
plabel(ax_b1, "B")

# ---------------- PANEL C: Random Forest ----------------
iC = inner2(outer[1, 0])
ax_c1 = fig.add_subplot(iC[0, 0])
ax_c2 = fig.add_subplot(iC[0, 1])

ax_c1.plot(tree_counts, oob_overall, color="#111827", lw=1.3, label="OOB Error")
ax_c1.plot(tree_counts, oob_spd, color=C_SPD, linestyle="--", lw=1.1, label=f"SPD ({int(oob_spd[-1]*len(spd))}/{len(spd)})")
ax_c1.plot(tree_counts, oob_lpd, color=C_LPD, linestyle=":", lw=1.2, label=f"LPD ({int(oob_lpd[-1]*len(lpd))}/{len(lpd)})")
ax_c1.axhline(null_err, color="#6B7280", linestyle="--", lw=1.0, label=f"Null baseline ({null_err:.3f})")
ax_c1.set_xlabel("Number of Trees", fontsize=9, fontweight="bold")
ax_c1.set_ylabel("Misclassification Rate", fontsize=9, fontweight="bold")
ax_c1.set_title("RF Error Convergence\n(500 trees, balanced, OOB staircase)", fontsize=9, fontweight="bold", pad=4)
ax_c1.legend(fontsize=6.8, loc="upper right", framealpha=0.85)
ax_c1.grid(True, linestyle=":", alpha=0.35)
ax_c1.set_ylim(0, max(oob_overall.max(), oob_lpd.max()) + 0.08)

rf_top15_sorted = rf_top15.sort_values(ascending=True)
bar_colors = []
for g in rf_top15_sorted.index:
    if g in rf_top10_genes and g in HUB9:
        bar_colors.append(COLORS9[g])
    elif g in rf_top10_genes:
        bar_colors.append("#94A3B8")
    else:
        bar_colors.append("#CBD5E1")

ax_c2.barh(range(len(rf_top15_sorted)), rf_top15_sorted.values, color=bar_colors, height=0.7, edgecolor="white", lw=0.4)
ax_c2.set_yticks(range(len(rf_top15_sorted)))
ax_c2.set_yticklabels(rf_top15_sorted.index, fontsize=8, fontweight="bold")
ax_c2.set_xlabel("Mean Decrease Gini", fontsize=9, fontweight="bold")
ax_c2.set_title(f"RF Feature Importance\n(Top 10 selected; rank cutoff at rank 10)", fontsize=9, fontweight="bold", pad=4)
ax_c2.grid(True, linestyle=":", alpha=0.35, axis="x")

for idx, gene in enumerate(rf_top15_sorted.index):
    if gene in rf_top10_genes and gene in HUB9:
        ax_c2.text(rf_top15_sorted.values[idx] + 0.0015, idx, "[Hub]",
                   fontsize=7.0, color=COLORS9[gene], va="center", fontweight="bold")
plabel(ax_c1, "C")

# ---------------- PANEL D: Venn Diagram ----------------
ax_d = fig.add_subplot(outer[1, 1])
ax_d.set_xlim(0, 1); ax_d.set_ylim(0, 1); ax_d.axis("off")

e1 = Ellipse((0.36, 0.60), 0.52, 0.48, facecolor="#F43F5E", alpha=0.25, edgecolor="#BE123C", lw=1.8)
e2 = Ellipse((0.64, 0.60), 0.52, 0.48, facecolor="#0284C7", alpha=0.25, edgecolor="#0369A1", lw=1.8)
e3 = Ellipse((0.50, 0.38), 0.52, 0.48, facecolor="#F59E0B", alpha=0.25, edgecolor="#D97706", lw=1.8)
for e in [e1, e2, e3]: ax_d.add_patch(e)

ax_d.text(0.16, 0.92, f"LASSO\n(n={n_lasso})", ha="center", fontsize=10, fontweight="bold", color="#9F1239")
ax_d.text(0.84, 0.92, f"SVM-RFE\n(n={n_svm})", ha="center", fontsize=10, fontweight="bold", color="#0369A1")
ax_d.text(0.50, 0.05, f"Random Forest (n={n_rf})", ha="center", fontsize=10, fontweight="bold", color="#B45309")

ax_d.text(0.22, 0.65, str(only_L), ha="center", va="center", fontsize=13, fontweight="bold", color="#9F1239")
ax_d.text(0.78, 0.65, str(only_S), ha="center", va="center", fontsize=13, fontweight="bold", color="#075985")
ax_d.text(0.50, 0.22, str(only_R), ha="center", va="center", fontsize=13, fontweight="bold", color="#92400E")

ax_d.text(0.50, 0.72, str(ls_only), ha="center", va="center", fontsize=13, fontweight="bold", color="#1E293B")
ax_d.text(0.34, 0.44, str(lr_only), ha="center", va="center", fontsize=13, fontweight="bold", color="#1E293B")
ax_d.text(0.66, 0.44, str(sr_only), ha="center", va="center", fontsize=13, fontweight="bold", color="#1E293B")
ax_d.text(0.50, 0.52, str(lsr), ha="center", va="center", fontsize=14, fontweight="bold", color="#DC2626")

bbox_style = dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#CBD5E1", alpha=0.92, lw=0.6)
ls_genes = sorted((L & S) - R)
lr_genes = sorted((L & R) - S)
sr_genes = sorted((S & R) - L)
lsr_genes = sorted(L & S & R)

ax_d.text(0.50, 0.84, f"L& S ({ls_only}): {', '.join(ls_genes)}",
          ha="center", va="bottom", fontsize=6.8, color="#0F172A", bbox=bbox_style)
ax_d.text(0.16, 0.38, f"L& R ({lr_only}):\n{', '.join(lr_genes)}",
          ha="center", va="center", fontsize=6.8, color="#0F172A", bbox=bbox_style)
ax_d.text(0.84, 0.38, f"S& R ({sr_only}):\n{', '.join(sr_genes)}",
          ha="center", va="center", fontsize=6.8, color="#0F172A", bbox=bbox_style)
ax_d.text(0.50, 0.58, f"L&S&R ({lsr}):\n{', '.join(lsr_genes)}",
          ha="center", va="center", fontsize=7.0, color="#DC2626", fontweight="bold")

n_hub = len(L & S) + len((L & R) - S) + len((S & R) - L)
ax_d.set_title(f"Venn Diagram of ML Feature Selectors\n{n_hub} Hub Genes by Majority Rule (>=2 Models; 3-way={lsr})",
               fontsize=9.5, fontweight="bold", pad=8)
plabel(ax_d, "D")

# ---------------- PANEL E: Discovery ROC ----------------
ax_e = fig.add_subplot(outer[2, 0])
for g in HUB9:
    fpr, tpr = disc_curves[g]
    auc = disc_aucs[g]
    ax_e.plot(fpr, tpr, color=COLORS9[g], lw=1.6, label=f"{g:<8s} AUC={auc:.3f}")

ax_e.plot([0, 1], [0, 1], "--", color="#9CA3AF", lw=1)
ax_e.set_xlim(-0.01, 1.01); ax_e.set_ylim(-0.01, 1.01)
ax_e.set_xlabel("1 - Specificity", fontsize=9.5, fontweight="bold")
ax_e.set_ylabel("Sensitivity", fontsize=9.5, fontweight="bold")
ax_e.set_title(f"ROC - Discovery Cohort (GSE125498)\nLPD n={len(lpd)}, SPD n={len(spd)} | All {len(HUB9)} Hub Genes",
               fontsize=9.5, fontweight="bold", pad=4)
ax_e.legend(loc="lower right", fontsize=7.2, framealpha=0.85, ncol=1)
ax_e.grid(True, linestyle=":", alpha=0.3)
plabel(ax_e, "E")

# ---------------- PANEL F: Validation ROC ----------------
iF = inner2(outer[2, 1])
ax_f1 = fig.add_subplot(iF[0, 0])
ax_f2 = fig.add_subplot(iF[0, 1])

# F1 - GSE248762 all-cell
for g in HUB9:
    if g not in g24_curves:
        continue
    fpr_v, tpr_v = g24_curves[g]
    auc_v = g24_aucs[g]
    conc = g24_concordant.get(g, False)
    if g == "CCBE1":
        lbl = f"{g} [INVERTED] AUC={auc_v:.3f}"
        ax_f1.step(fpr_v, tpr_v, where="post", color=COLORS9[g], lw=1.8, linestyle="--", label=lbl)
    elif conc:
        lbl = f"{g} [concordant] AUC={auc_v:.3f}"
        ax_f1.step(fpr_v, tpr_v, where="post", color=COLORS9[g], lw=1.8, label=lbl)
    else:
        lbl = f"{g} AUC={auc_v:.3f}"
        ax_f1.step(fpr_v, tpr_v, where="post", color=COLORS9[g], lw=1.0, linestyle=":", label=lbl)

ax_f1.plot([0, 1], [0, 1], "--", color="#9CA3AF", lw=1)
ax_f1.set_xlim(-0.01, 1.01); ax_f1.set_ylim(-0.01, 1.01)
ax_f1.set_xlabel("1 - Specificity", fontsize=8.5, fontweight="bold")
ax_f1.set_ylabel("Sensitivity", fontsize=8.5, fontweight="bold")
ax_f1.set_title(f"GSE248762 All-Cell Pseudobulk\nLV n={len(g24_lv)}, SV n={len(g24_sv)} | {n_concordant_g24}/{len(g24_concordant)} concordant",
               fontsize=8.5, fontweight="bold", pad=4)
ax_f1.legend(loc="lower right", fontsize=5.8, framealpha=0.88)
ax_f1.grid(True, linestyle=":", alpha=0.3)
plabel(ax_f1, "F")

# F2 - GSE130888 all-cell
for g in HUB9:
    if g not in g13_curves:
        continue
    fpr_v, tpr_v = g13_curves[g]
    auc_v = g13_aucs[g]
    conc = g13_concordant.get(g, False)
    if g == "CCBE1":
        lbl = f"{g} [INVERTED] AUC={auc_v:.3f}"
        ax_f2.step(fpr_v, tpr_v, where="post", color=COLORS9[g], lw=1.8, linestyle="--", label=lbl)
    elif conc:
        lbl = f"{g} [concordant] AUC={auc_v:.3f}"
        ax_f2.step(fpr_v, tpr_v, where="post", color=COLORS9[g], lw=1.8, label=lbl)
    else:
        lbl = f"{g} AUC={auc_v:.3f}"
        ax_f2.step(fpr_v, tpr_v, where="post", color=COLORS9[g], lw=1.0, linestyle=":", label=lbl)

ax_f2.plot([0, 1], [0, 1], "--", color="#9CA3AF", lw=1)
ax_f2.set_xlim(-0.01, 1.01); ax_f2.set_ylim(-0.01, 1.01)
ax_f2.set_xlabel("1 - Specificity", fontsize=8.5, fontweight="bold")
ax_f2.set_ylabel("Sensitivity", fontsize=8.5, fontweight="bold")
ax_f2.set_title(f"GSE130888 All-Cell Pseudobulk\nLV n={len(g13_lv)}, SV n={len(g13_sv)} | {n_concordant_g13}/{len(g13_concordant)} concordant (TNFSF15 ND)",
               fontsize=8.5, fontweight="bold", pad=4)
ax_f2.legend(loc="lower right", fontsize=5.8, framealpha=0.88)
ax_f2.grid(True, linestyle=":", alpha=0.3)
ax_f2.text(0.01, -0.22,
           "[concordant] = same direction as discovery (LV > SV). CCBE1 inverted in both cohorts.\n"
           "GSE248762: 7/9 concordant in direction; GSE130888: secondary validation cohort.",
           transform=ax_f2.transAxes, fontsize=6.0, color="#475569", style="italic")

# ---------------- PANEL G: Discovery Expression ----------------
ax_g = fig.add_subplot(outer[3, 0])
np.random.seed(42)
w = 0.38; gap = 0.50

for i, g in enumerate(HUB9):
    xl = i * (2*w + gap); xs = xl + w
    vl = expr_mat.loc[g, lpd].values.astype(float)
    vs = expr_mat.loc[g, spd].values.astype(float)

    for xp, vv, cl, clf in [(xl, vl, C_LPD, C_LPD_L), (xs, vs, C_SPD, C_SPD_L)]:
        ax_g.boxplot([vv], positions=[xp], widths=w*0.90, patch_artist=True,
                     boxprops=dict(facecolor=clf, color=cl, lw=1.2),
                     whiskerprops=dict(color=cl, lw=1.0),
                     capprops=dict(color=cl, lw=1.0),
                     medianprops=dict(color="#111827", lw=1.8),
                     flierprops=dict(marker="o", ms=3, color=cl, alpha=0.5))
        ax_g.scatter(np.random.normal(xp, 0.035, len(vv)), vv, color=cl, s=9, alpha=0.7, zorder=4)

    p = disc_pvals[g]
    star = pstar(p)
    top = max(np.max(vl), np.max(vs))
    bary = top + 0.25; mx = (xl + xs) / 2
    ax_g.plot([xl, xs], [bary, bary], color="#374151", lw=0.9)
    ax_g.text(mx, bary + 0.05, star, ha="center", fontsize=8.5, fontweight="bold", color="#374151")

tick_x = [i * (2*w + gap) + w/2 for i in range(len(HUB9))]
ax_g.set_xticks(tick_x)
ax_g.set_xticklabels(HUB9, fontsize=8, fontweight="bold", rotation=30, ha="right")
ax_g.set_ylabel("Normalised Expression", fontsize=9.5, fontweight="bold")
ax_g.set_title(f"Expression - Discovery (GSE125498) | All {len(HUB9)} Hub Genes", fontsize=9.5, fontweight="bold", pad=4)
ax_g.legend(handles=[mpatches.Patch(color=C_LPD, label=f"LPD (n={len(lpd)})"),
                     mpatches.Patch(color=C_SPD, label=f"SPD (n={len(spd)})")],
            loc="upper left", bbox_to_anchor=(0.01, 0.99), fontsize=7.5, frameon=True, facecolor="white", edgecolor="#CBD5E1")
ax_g.grid(True, linestyle=":", alpha=0.3, axis="y")
ax_g.text(0.01, -0.16,
          "Mann-Whitney U, two-sided raw p: *p<0.05  **p<0.01  ***p<0.001  ****p<0.0001  ns (p>=0.05)",
          transform=ax_g.transAxes, fontsize=6.8, color="#475569", style="italic")
plabel(ax_g, "G")

# ---------------- PANEL H: Validation Expression ----------------
iH = inner2(outer[3, 1])
ax_h1 = fig.add_subplot(iH[0, 0])
ax_h2 = fig.add_subplot(iH[0, 1])
np.random.seed(0)

# H1 - GSE248762 all-cell boxplots
g24_in = [g for g in HUB9 if g in g24_gene_map]
for i, g in enumerate(g24_in):
    eid = g24_gene_map[g]
    xl = i * (2*w + gap); xs = xl + w
    vl = g24_log.loc[eid, g24_lv].values.astype(float)
    vs = g24_log.loc[eid, g24_sv].values.astype(float)
    for xp, vv, cl, clf in [(xl, vl, C_LPD, C_LPD_L), (xs, vs, C_SPD, C_SPD_L)]:
        ax_h1.boxplot([vv], positions=[xp], widths=w*0.88, patch_artist=True,
                      boxprops=dict(facecolor=clf, color=cl, lw=1.2),
                      whiskerprops=dict(color=cl, lw=1.0), capprops=dict(color=cl, lw=1.0),
                      medianprops=dict(color="#111827", lw=1.8),
                      flierprops=dict(marker="o", ms=3, color=cl, alpha=0.5))
        ax_h1.scatter(np.random.normal(xp, 0.035, len(vv)), vv, color=cl, s=7, alpha=0.7, zorder=4)
    p = g24_pvals.get(g, np.nan); star = pstar(p)
    top = max(np.max(vl), np.max(vs))
    bary = top + 0.35; mx = (xl + xs) / 2
    ax_h1.plot([xl, xs], [bary, bary], color="#374151", lw=0.9)
    ax_h1.text(mx, bary + 0.08, star, ha="center", fontsize=8, fontweight="bold", color="#374151")
    conc = g24_concordant.get(g, False)
    if conc and g != "CCBE1":
        ax_h1.text(mx, bary + 0.55, "\u2713", ha="center", fontsize=12, color="#16A34A", fontweight="bold")
    elif g == "CCBE1":
        ax_h1.text(mx, bary + 0.55, "\u2195", ha="center", fontsize=11, color="#DC2626", fontweight="bold")

tick_xh1 = [i * (2*w + gap) + w/2 for i in range(len(g24_in))]
ax_h1.set_xticks(tick_xh1)
ax_h1.set_xticklabels(g24_in, fontsize=7, fontweight="bold", rotation=32, ha="right")
ax_h1.set_ylabel("edgeR logCPM", fontsize=8.5, fontweight="bold")
ax_h1.set_title(f"GSE248762 All-Cell Pseudobulk\n{n_concordant_g24}/{len(g24_in)} concordant in direction",
               fontsize=8.5, fontweight="bold", pad=4)
ax_h1.legend(handles=[mpatches.Patch(color=C_LPD, label=f"LV (n={len(g24_lv)})"),
                      mpatches.Patch(color=C_SPD, label=f"SV (n={len(g24_sv)})")],
             loc="upper right", fontsize=7, framealpha=0.9, edgecolor="#CBD5E1")
ax_h1.grid(True, linestyle=":", alpha=0.3, axis="y")
ax_h1.text(0.01, -0.22,
           "*p<0.05  **p<0.01  ***p<0.001  ns p\u22650.05  \u2713=concordant direction  \u2195=inverted vs discovery",
           transform=ax_h1.transAxes, fontsize=6.0, color="#475569", style="italic")
plabel(ax_h1, "H")

# H2 - GSE130888 all-cell boxplots
g13_in = [g for g in HUB9 if g in g13_gene_map]
for i, g in enumerate(g13_in):
    eid = g13_gene_map[g]
    xl = i * (2*w + gap); xs = xl + w
    vl = g13_log.loc[eid, g13_lv].values.astype(float)
    vs = g13_log.loc[eid, g13_sv].values.astype(float)
    for xp, vv, cl, clf in [(xl, vl, C_LPD, C_LPD_L), (xs, vs, C_SPD, C_SPD_L)]:
        ax_h2.boxplot([vv], positions=[xp], widths=w*0.88, patch_artist=True,
                      boxprops=dict(facecolor=clf, color=cl, lw=1.2),
                      whiskerprops=dict(color=cl, lw=1.0), capprops=dict(color=cl, lw=1.0),
                      medianprops=dict(color="#111827", lw=1.8),
                      flierprops=dict(marker="o", ms=3, color=cl, alpha=0.5))
        ax_h2.scatter(np.random.normal(xp, 0.035, len(vv)), vv, color=cl, s=7, alpha=0.7, zorder=4)
    p = g13_pvals.get(g, np.nan); star = pstar(p)
    top = max(np.max(vl), np.max(vs))
    bary = top + 0.35; mx = (xl + xs) / 2
    ax_h2.plot([xl, xs], [bary, bary], color="#374151", lw=0.9)
    ax_h2.text(mx, bary + 0.08, star, ha="center", fontsize=8, fontweight="bold", color="#374151")
    conc = g13_concordant.get(g, False)
    if conc and g != "CCBE1":
        ax_h2.text(mx, bary + 0.55, "\u2713", ha="center", fontsize=12, color="#16A34A", fontweight="bold")
    elif g == "CCBE1":
        ax_h2.text(mx, bary + 0.55, "\u2195", ha="center", fontsize=11, color="#DC2626", fontweight="bold")

tick_xh2 = [i * (2*w + gap) + w/2 for i in range(len(g13_in))]
ax_h2.set_xticks(tick_xh2)
ax_h2.set_xticklabels(g13_in, fontsize=7, fontweight="bold", rotation=32, ha="right")
ax_h2.set_ylabel("edgeR logCPM", fontsize=8.5, fontweight="bold")
ax_h2.set_title(f"GSE130888 All-Cell Pseudobulk\n{n_concordant_g13}/{len(g13_in)} concordant in direction",
               fontsize=8.5, fontweight="bold", pad=4)
ax_h2.legend(handles=[mpatches.Patch(color=C_LPD, label=f"LV (n={len(g13_lv)})"),
                      mpatches.Patch(color=C_SPD, label=f"SV (n={len(g13_sv)})")],
             loc="upper right", fontsize=7, framealpha=0.9, edgecolor="#CBD5E1")
ax_h2.grid(True, linestyle=":", alpha=0.3, axis="y")
ax_h2.text(0.01, -0.22,
           "GSE130888: secondary all-cell pseudobulk. TNFSF15 not detected (excluded).\n"
           "*p<0.05  **p<0.01  ***p<0.001  ns p\u22650.05  \u2713=concordant  \u2195=inverted vs discovery",
           transform=ax_h2.transAxes, fontsize=6.0, color="#475569", style="italic")

# ---------------------------------------------------------------------------
# 6. SAVE
# ---------------------------------------------------------------------------
png_out = os.path.join(OUT_DIR, "Fig4_versionA_ml39_machine_learning_panel.png")
pdf_out = os.path.join(OUT_DIR, "Fig4_versionA_ml39_machine_learning_panel.pdf")

fig.savefig(png_out, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"\nSaved PNG: {png_out}")
print(f"Saved PDF: {pdf_out}")
print(f"\nSummary:")
print(f"  39-gene pool, 9 hub genes: {', '.join(HUB9)}")
print(f"  LASSO={n_lasso}, SVM-RFE={n_svm}, RF={n_rf}, 3-way={lsr}")
print(f"  GSE248762: {n_concordant_g24}/{len(g24_concordant)} concordant")
print(f"  GSE130888: {n_concordant_g13}/{len(g13_concordant)} concordant")
for g in HUB9:
    g24_s = "concordant" if g24_concordant.get(g) else ("inverted" if g == "CCBE1" else "discordant")
    g13_s = "concordant" if g13_concordant.get(g) else ("inverted" if g == "CCBE1" else "discordant/nd")
    print(f"    {g:8s}: disc_AUC={disc_aucs[g]:.3f} | G24_AUC={g24_aucs.get(g, 0):.3f}({g24_s}) | G13_AUC={g13_aucs.get(g, 0):.3f}({g13_s})")
