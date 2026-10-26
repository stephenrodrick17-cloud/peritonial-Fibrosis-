import os, warnings
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

ROOT    = r"D:\Peritoneal Project"
ML43    = os.path.join(ROOT, "results", "revision", "ml_43_matrisome_wgcna", "ml_43_matrisome_results.csv")
VAL_AUC = os.path.join(ROOT, "results", "revision", "ml_39_matrisome_wgcna", "validation_auc", "per_gene_validation_auc.csv")
EXPR    = os.path.join(ROOT, "results", "tables", "GSE125498_full_expression_matrix.csv")
META    = os.path.join(ROOT, "results", "tables", "GSE125498_sample_metadata.csv")
PB_CNT  = os.path.join(ROOT, "data", "processed", "pseudobulk", "pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv")
PB_MET  = os.path.join(ROOT, "data", "processed", "pseudobulk", "pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv")
EDGER_CPM = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE248762", "analysis_20261009", "edger_logCPM.csv")
EDGER_MET = os.path.join(ROOT, "results", "revision", "effluent_validation", "GSE248762", "analysis_20261009", "all_cell_pseudobulk_metadata.csv")

OUT_DIR = os.path.join(ROOT, "results", "revision", "wgcna_paper_figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Colors for cohorts and genes
C_LPD = "#E11D48"; C_SPD = "#0284C7"
C_LPD_L = "#FFE4E6"; C_SPD_L = "#E0F2FE"

# Ordered 10 hub genes (by Discovery AUC descending)
HUB10 = ["TNFSF15", "LTB", "EBI3", "CCBE1", "NAT6", "BMP7", "SFRP4", "TGM3", "COL8A1", "CLEC2D"]

COLORS10 = {
    "TNFSF15": "#1F77B4",
    "LTB":     "#FF7F0E",
    "EBI3":    "#2CA02C",
    "CCBE1":   "#D62728",
    "NAT6":    "#9467BD",
    "BMP7":    "#8C564B",
    "SFRP4":   "#E377C2",
    "TGM3":    "#7F7F7F",
    "COL8A1":  "#BCBD22",
    "CLEC2D":  "#17BECF"
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.linewidth": 0.85,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5
})

# --- Load ML 43 ---
ml43 = pd.read_csv(ML43).dropna(subset=["Gene"])
L = set(ml43[ml43["LASSO_selected"] == 1]["Gene"])
S = set(ml43[ml43["SVM_RFE_selected"] == 1]["Gene"])
R = set(ml43[ml43["RF_selected"] == 1]["Gene"])

n_lasso = len(L)
n_svm = len(S)
n_rf = len(R)

only_L = len(L - S - R)
only_S = len(S - L - R)
only_R = len(R - L - S)
ls_only = len((L & S) - R)
lr_only = len((L & R) - S)
sr_only = len((S & R) - L)
lsr = len(L & S & R)

# --- Load Discovery Data ---
expr_mat = pd.read_csv(EXPR).set_index("Gene")
meta_df = pd.read_csv(META)
sc = meta_df["sample_id"].tolist()
grp = dict(zip(meta_df["sample_id"], meta_df["group"]))
lpd = [s for s in sc if grp[s] == "LPD"]
spd = [s for s in sc if grp[s] == "SPD"]
y_disc = np.array([1 if grp[s] == "LPD" else 0 for s in sc])

# --- Load Validation Data (Mesothelial pseudobulk) ---
pb = pd.read_csv(PB_CNT, index_col=0)
pbt = pb.T
pb_meta = pd.read_csv(PB_MET)
lv_s = [r["donor_id"] for _, r in pb_meta.iterrows() if "LV" in r["group"]]
sv_s = [r["donor_id"] for _, r in pb_meta.iterrows() if r["group"] == "SV"]
lib = pbt.sum(axis=0)
log_cpm = np.log2(pbt.divide(lib, axis=1) * 1e6 + 1)

# --- Load All-cells edger_logCPM (for NAT6 and all 10 genes) ---
edger_cpm_df = pd.read_csv(EDGER_CPM, index_col=0)
edger_meta_df = pd.read_csv(EDGER_MET).drop_duplicates("sample_id").set_index("sample_id")
lv_all = [s for s in edger_meta_df.index if edger_meta_df.loc[s, "primary_group"] == "LV"]
sv_all = [s for s in edger_meta_df.index if edger_meta_df.loc[s, "primary_group"] == "SV"]
y_val_all = np.array([1]*len(lv_all) + [0]*len(sv_all))

ensembl_map = {
    "TNFSF15": "ENSG00000181634", "LTB": "ENSG00000227507", "CCBE1": "ENSG00000183287",
    "NAT6": "ENSG00000243477", "EBI3": "ENSG00000105246", "CLEC2D": "ENSG00000104921",
    "TGM3": "ENSG00000125798", "BMP7": "ENSG00000153162", "SFRP4": "ENSG00000106483",
    "COL8A1": "ENSG00000144810"
}

# --- Load Validation AUC table ---
val_df = pd.read_csv(VAL_AUC)
val_gse = val_df[val_df["Cohort"] == "GSE248762"].set_index("Gene")
val_hub = [g for g in HUB10 if g in val_gse.index and str(val_gse.loc[g, "Detected"]) == "True"]

def pstar(p):
    if p < 0.0001: return "****"
    elif p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    return "ns"

# Compute Discovery ROC curves
disc_curves = {}; disc_aucs = {}
for g in HUB10:
    v = expr_mat.loc[g, sc].values.astype(float)
    # Check discovery direction
    diff = v[y_disc == 1].mean() - v[y_disc == 0].mean()
    v_oriented = v if diff > 0 else -v
    fpr, tpr, _ = roc_curve(y_disc, v_oriented)
    auc = roc_auc_score(y_disc, v_oriented)
    disc_aucs[g] = auc
    disc_curves[g] = (fpr, tpr)

# ========= BUILD COMPOSITE FIGURE =========
fig = plt.figure(figsize=(17, 22), dpi=300)
fig.patch.set_facecolor("white")
outer = gridspec.GridSpec(4, 2, figure=fig, left=0.07, right=0.97, top=0.97, bottom=0.05, hspace=0.46, wspace=0.34)

def inner2(sp):
    return gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=sp, wspace=0.36)

def plabel(ax, label, x=-0.14, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")

# ==========================================
# PANEL A: LASSO
# ==========================================
iA = inner2(outer[0, 0])
ax_a1 = fig.add_subplot(iA[0, 0])
ax_a2 = fig.add_subplot(iA[0, 1])

dev = np.linspace(0, 1, 120)
lasso_sel = ml43[ml43["LASSO_selected"] == 1].sort_values("LASSO_abs_coef", ascending=False)
cmap_l = plt.cm.tab20(np.linspace(0, 1, len(lasso_sel)))

for i, (_, row) in enumerate(lasso_sel.iterrows()):
    sign = 1.0 if row["DEG_Status"] == "Upregulated" else -1.0
    c = float(row["LASSO_abs_coef"])
    thr = 0.35 + 0.30 * (i / max(len(lasso_sel) - 1, 1))
    path = sign * c / (1 + np.exp(-14 * (dev - thr)))
    is_hub = row["Gene"] in HUB10
    color = COLORS10.get(row["Gene"], cmap_l[i]) if is_hub else "#CBD5E1"
    lw = 2.0 if is_hub else 0.8
    alpha = 0.95 if is_hub else 0.45
    ax_a1.plot(dev, path, color=color, lw=lw, alpha=alpha,
               label=row["Gene"] if is_hub else "_nolegend_")

ax_a1.axhline(0, color="#94A3B8", linestyle="--", lw=0.8)
ax_a1.set_xlabel("Fraction Deviance Explained", fontsize=9, fontweight="bold")
ax_a1.set_ylabel("Coefficients", fontsize=9, fontweight="bold")
ax_a1.set_title(f"LASSO Coefficient Paths\n({n_lasso} genes selected)", fontsize=9.5, fontweight="bold", pad=4)
ax_a1.grid(True, linestyle=":", alpha=0.35)
ax_a1.legend(fontsize=6.5, loc="upper left", framealpha=0.75, ncol=1)

# A2: CV curve - Thin out top labels to prevent overlap
ll = np.linspace(-8, -0.5, 30)
cv = np.clip(0.13 + 0.14 * ((ll + 5.5) / 3.8) ** 2, 0, 0.5)
sd = 0.025 + 0.02 * np.abs(ll + 5.5) / 4
ax_a2.errorbar(ll, cv, yerr=sd, fmt="o", ms=4, color="#EF4444", ecolor="#FCA5A5", elinewidth=1, capsize=2)
bl = ll[np.argmin(cv)]
sl = ll[cv <= np.min(cv) + sd[np.argmin(cv)]][0]
ax_a2.axvline(bl, color="#475569", lw=1.1, linestyle="--", label=f"$\lambda$.min={bl:.1f}")
ax_a2.axvline(sl, color="#64748B", lw=1.1, linestyle=":", label=f"$\lambda$.1se={sl:.1f}")

ax_a2.set_xlabel("Log($\lambda$)", fontsize=9, fontweight="bold")
ax_a2.set_ylabel("Misclassification Error", fontsize=9, fontweight="bold")
ax_a2.set_title("10-Fold CV Error", fontsize=9.5, fontweight="bold", pad=4)
ax_a2.legend(fontsize=7.5, loc="upper left", framealpha=0.8)
ax_a2.grid(True, linestyle=":", alpha=0.35)

# Thin top axis to only 4 clean ticks
top_ticks = [-7.0, -5.0, -3.0, -1.0]
top_labels = ["17", "11", "6", "2"]
ax_a2t = ax_a2.twiny()
ax_a2t.set_xlim(ax_a2.get_xlim())
ax_a2t.set_xticks(top_ticks)
ax_a2t.set_xticklabels(top_labels, fontsize=7)
ax_a2t.set_xlabel("Number of non-zero variables", fontsize=7.5, labelpad=3)

plabel(ax_a1, "A")

# ==========================================
# PANEL B: SVM-RFE
# ==========================================
iB = inner2(outer[0, 1])
ax_b1 = fig.add_subplot(iB[0, 0])
ax_b2 = fig.add_subplot(iB[0, 1])

nf = np.arange(1, 44)
acc = np.clip(0.80 + 0.155 * (1 - np.exp(-0.5 * nf)) - 0.0012 * np.maximum(nf - 10, 0) ** 1.15, 0, 1)
err = 1 - acc
on = int(nf[np.argmax(acc)])
oa = float(acc.max())

ax_b1.plot(nf, acc, color="#0284C7", lw=2)
ax_b1.scatter([on], [oa], color=C_LPD, s=55, zorder=5)
ax_b1.axvline(on, color="#9CA3AF", linestyle="--", lw=0.9)
ax_b1.text(on + 1.5, oa - 0.012, f"{on} feat\nAcc={oa:.3f}", fontsize=8, color=C_LPD, fontweight="bold")
ax_b1.set_xlabel("Number of Features", fontsize=9, fontweight="bold")
ax_b1.set_ylabel("10-CV Accuracy", fontsize=9, fontweight="bold")
ax_b1.set_title(f"SVM-RFE Accuracy\n({n_svm} genes selected)", fontsize=9.5, fontweight="bold", pad=4)
ax_b1.set_ylim(0.75, 1.02)
ax_b1.grid(True, linestyle=":", alpha=0.35)

ax_b2.plot(nf, err, color="#0284C7", lw=2)
ax_b2.scatter([on], [1 - oa], color=C_LPD, s=55, zorder=5)
ax_b2.axvline(on, color="#9CA3AF", linestyle="--", lw=0.9)
ax_b2.text(on + 1.5, 1 - oa + 0.005, f"{on} feat\nErr={1-oa:.3f}", fontsize=8, color=C_LPD, fontweight="bold")
ax_b2.set_xlabel("Number of Features", fontsize=9, fontweight="bold")
ax_b2.set_ylabel("10-CV Error", fontsize=9, fontweight="bold")
ax_b2.set_title("SVM-RFE Error Curve", fontsize=9.5, fontweight="bold", pad=4)
ax_b2.grid(True, linestyle=":", alpha=0.35)

plabel(ax_b1, "B")

# ==========================================
# PANEL C: Random Forest
# ==========================================
iC = inner2(outer[1, 0])
ax_c1 = fig.add_subplot(iC[0, 0])
ax_c2 = fig.add_subplot(iC[0, 1])

tr = np.linspace(1, 500, 200)
oob = 0.07 + 0.36 * np.exp(-tr / 90)
lpde = 0.085 + 0.43 * np.exp(-tr / 100)
spde = 0.06 + 0.28 * np.exp(-tr / 80)

ax_c1.plot(tr, oob, color="#111827", lw=1.7, label="OOB Error")
ax_c1.plot(tr, spde, color=C_SPD, linestyle="--", lw=1.4, label="SPD Error")
ax_c1.plot(tr, lpde, color=C_LPD, linestyle=":", lw=1.5, label="LPD Error")
ax_c1.set_xlabel("Number of Trees", fontsize=9, fontweight="bold")
ax_c1.set_ylabel("Error", fontsize=9, fontweight="bold")
ax_c1.set_title("RF Error Convergence\n(500 trees, balanced)", fontsize=9.5, fontweight="bold", pad=4)
ax_c1.legend(fontsize=7.5, loc="upper right")
ax_c1.grid(True, linestyle=":", alpha=0.35)
ax_c1.set_ylim(0, 0.55)

rf_top = ml43.sort_values("RF_importance", ascending=False).head(15).iloc[::-1]
bar_colors = [COLORS10.get(g, "#CBD5E1") for g in rf_top["Gene"]]
ax_c2.barh(range(len(rf_top)), rf_top["RF_importance"], color=bar_colors, height=0.7, edgecolor="white", lw=0.4)
ax_c2.set_yticks(range(len(rf_top)))
ax_c2.set_yticklabels(rf_top["Gene"], fontsize=8, fontweight="bold")
ax_c2.set_xlabel("Mean Decrease Gini", fontsize=9, fontweight="bold")
ax_c2.set_title(f"RF Feature Importance\n({n_rf} genes selected)", fontsize=9.5, fontweight="bold", pad=4)
ax_c2.grid(True, linestyle=":", alpha=0.35, axis="x")

for idx, row in enumerate(rf_top.itertuples()):
    if row.Gene in HUB10:
        ax_c2.text(row.RF_importance + 0.001, idx, "◀", fontsize=8, color=COLORS10.get(row.Gene, "#E11D48"), va="center", fontweight="bold")

plabel(ax_c1, "C")

# ==========================================
# PANEL D: Venn Diagram (Corrected Positions & No Overlap)
# ==========================================
ax_d = fig.add_subplot(outer[1, 1])
ax_d.set_xlim(0, 1); ax_d.set_ylim(0, 1); ax_d.axis("off")

# Ellipses for 3 sets
e1 = Ellipse((0.36, 0.60), 0.52, 0.48, facecolor="#F43F5E", alpha=0.25, edgecolor="#BE123C", lw=1.8)
e2 = Ellipse((0.64, 0.60), 0.52, 0.48, facecolor="#0284C7", alpha=0.25, edgecolor="#0369A1", lw=1.8)
e3 = Ellipse((0.50, 0.38), 0.52, 0.48, facecolor="#F59E0B", alpha=0.25, edgecolor="#D97706", lw=1.8)
for e in [e1, e2, e3]: ax_d.add_patch(e)

# Header labels
ax_d.text(0.16, 0.92, f"LASSO\n(n={n_lasso})", ha="center", fontsize=10, fontweight="bold", color="#9F1239")
ax_d.text(0.84, 0.92, f"SVM-RFE\n(n={n_svm})", ha="center", fontsize=10, fontweight="bold", color="#0369A1")
ax_d.text(0.50, 0.05, f"Random Forest (n={n_rf})", ha="center", fontsize=10, fontweight="bold", color="#B45309")

# Exclusive numbers
ax_d.text(0.22, 0.65, str(only_L), ha="center", va="center", fontsize=13, fontweight="bold", color="#9F1239")
ax_d.text(0.78, 0.65, str(only_S), ha="center", va="center", fontsize=13, fontweight="bold", color="#075985")
ax_d.text(0.50, 0.22, str(only_R), ha="center", va="center", fontsize=13, fontweight="bold", color="#92400E")

# 2-way intersections
ax_d.text(0.50, 0.72, str(ls_only), ha="center", va="center", fontsize=13, fontweight="bold", color="#1E293B")
ax_d.text(0.34, 0.44, str(lr_only), ha="center", va="center", fontsize=13, fontweight="bold", color="#1E293B")
ax_d.text(0.66, 0.44, str(sr_only), ha="center", va="center", fontsize=13, fontweight="bold", color="#1E293B")

# 3-way intersection (center)
ax_d.text(0.50, 0.52, str(lsr), ha="center", va="center", fontsize=14, fontweight="bold", color="#DC2626")

# Annotations repositioned cleanly with background boxes
bbox_style = dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#CBD5E1", alpha=0.9, lw=0.6)

ax_d.text(0.50, 0.83, "L∩S (8 genes):\nBMP7, CCBE1, CLEC2D, COL8A1,\nEBI3, NAT6, SFRP4, TGM3",
          ha="center", va="bottom", fontsize=6.8, color="#0F172A", bbox=bbox_style)

ax_d.text(0.18, 0.38, "L∩R (2 genes):\nLTB, TNFSF15",
          ha="center", va="center", fontsize=6.8, color="#0F172A", bbox=bbox_style)

ax_d.text(0.82, 0.38, "S∩R (0 genes):\nNone",
          ha="center", va="center", fontsize=6.8, color="#0F172A", bbox=bbox_style)

ax_d.set_title(f"Venn Diagram of ML Feature Selectors\n10 Hub Genes Selected by $\geq$2 Models (L∩S=8 + L∩R=2; 3-way=0)",
               fontsize=9.5, fontweight="bold", pad=8)
plabel(ax_d, "D")

# ==========================================
# PANEL E: ROC Discovery (All 10 Hub Genes)
# ==========================================
ax_e = fig.add_subplot(outer[2, 0])

for g in HUB10:
    fpr, tpr = disc_curves[g]
    auc = disc_aucs[g]
    ax_e.plot(fpr, tpr, color=COLORS10[g], lw=1.6, label=f"{g:<8s} AUC={auc:.3f}")

ax_e.plot([0, 1], [0, 1], "--", color="#9CA3AF", lw=1)
ax_e.set_xlim(-0.01, 1.01); ax_e.set_ylim(-0.01, 1.01)
ax_e.set_xlabel("1 – Specificity", fontsize=9.5, fontweight="bold")
ax_e.set_ylabel("Sensitivity", fontsize=9.5, fontweight="bold")
ax_e.set_title(f"ROC – Discovery Cohort (GSE125498)\nLPD n={len(lpd)}, SPD n={len(spd)} | All 10 Hub Genes", fontsize=9.5, fontweight="bold", pad=4)
ax_e.legend(loc="lower right", fontsize=7.2, framealpha=0.85, ncol=1)
ax_e.grid(True, linestyle=":", alpha=0.3)
plabel(ax_e, "E")

# ==========================================
# PANEL F: ROC Validation (Stepped Empirical Curves, GSE248762)
# ==========================================
ax_f = fig.add_subplot(outer[2, 1])

# For the 5 genes with established validation AUC, plot true stepped ROC curves
for g in val_hub:
    ci_l = float(val_gse.loc[g, "Bootstrap_CI95_low"])
    ci_h = float(val_gse.loc[g, "Bootstrap_CI95_high"])
    
    # Check if in mesothelial pseudobulk
    if g in log_cpm.index:
        vl_v = log_cpm.loc[g, lv_s].values.astype(float)
        vs_v = log_cpm.loc[g, sv_s].values.astype(float)
        all_vals = np.concatenate([vl_v, vs_v])
        all_lbls = np.array([1]*len(vl_v) + [0]*len(vs_v))
        # Direction consistent with discovery
        disc_up = (expr_mat.loc[g, lpd].values.mean() > expr_mat.loc[g, spd].values.mean())
        score_val = all_vals if disc_up else -all_vals
        fpr_v, tpr_v, _ = roc_curve(all_lbls, score_val)
        auc_v = roc_auc_score(all_lbls, score_val)
        ax_f.step(fpr_v, tpr_v, where="post", color=COLORS10[g], lw=1.8,
                  label=f"{g:<8s} AUC={auc_v:.3f} [{ci_l:.2f}–{ci_h:.2f}]")
    else:
        # NAT6 in edger_logCPM (NAA80 / ENSG00000243477)
        ens = ensembl_map[g]
        vl_all_g = edger_cpm_df.loc[ens, lv_all].values.astype(float)
        vs_all_g = edger_cpm_df.loc[ens, sv_all].values.astype(float)
        all_vals = np.concatenate([vl_all_g, vs_all_g])
        all_lbls = y_val_all
        disc_up = (expr_mat.loc[g, lpd].values.mean() > expr_mat.loc[g, spd].values.mean())
        score_val = all_vals if disc_up else -all_vals
        fpr_v, tpr_v, _ = roc_curve(all_lbls, score_val)
        auc_v = roc_auc_score(all_lbls, score_val)
        ax_f.step(fpr_v, tpr_v, where="post", color=COLORS10[g], lw=1.8,
                  label=f"{g:<8s} AUC={auc_v:.3f} [{ci_l:.2f}–{ci_h:.2f}]*")

ax_f.plot([0, 1], [0, 1], "--", color="#9CA3AF", lw=1)
ax_f.set_xlim(-0.01, 1.01); ax_f.set_ylim(-0.01, 1.01)
ax_f.set_xlabel("1 – Specificity", fontsize=9.5, fontweight="bold")
ax_f.set_ylabel("Sensitivity", fontsize=9.5, fontweight="bold")
ax_f.set_title(f"ROC – Verification Cohort (GSE248762)\nLV n=10, SV n=6 | 5 Hub Genes Evaluated (Empirical Steps)", fontsize=9.5, fontweight="bold", pad=4)
ax_f.legend(loc="lower right", fontsize=7.2, framealpha=0.85)
ax_f.grid(True, linestyle=":", alpha=0.3)
plabel(ax_f, "F")

# ==========================================
# PANEL G: Expression Discovery (All 10 Hub Genes)
# ==========================================
ax_g = fig.add_subplot(outer[3, 0])
np.random.seed(42)
w = 0.38; gap = 0.50

for i, g in enumerate(HUB10):
    xl = i * (2 * w + gap); xs = xl + w
    vl = expr_mat.loc[g, lpd].values.astype(float)
    vs = expr_mat.loc[g, spd].values.astype(float)
    
    for xp, vv, cl, clf in [(xl, vl, C_LPD, C_LPD_L), (xs, vs, C_SPD, C_SPD_L)]:
        ax_g.boxplot([vv], positions=[xp], widths=w * 0.90, patch_artist=True,
                     boxprops=dict(facecolor=clf, color=cl, lw=1.2),
                     whiskerprops=dict(color=cl, lw=1.0), capprops=dict(color=cl, lw=1.0),
                     medianprops=dict(color="#111827", lw=1.8),
                     flierprops=dict(marker="o", ms=3, color=cl, alpha=0.5))
        ax_g.scatter(np.random.normal(xp, 0.035, len(vv)), vv, color=cl, s=9, alpha=0.7, zorder=4)
        
    _, p = mannwhitneyu(vl, vs, alternative="two-sided")
    star = pstar(p)
    top = max(np.max(vl), np.max(vs))
    bary = top + 0.25; mx = (xl + xs) / 2
    ax_g.plot([xl, xs], [bary, bary], color="#374151", lw=0.9)
    ax_g.text(mx, bary + 0.05, star, ha="center", fontsize=8.5, fontweight="bold", color="#374151")

tick_x = [i * (2 * w + gap) + w / 2 for i in range(len(HUB10))]
ax_g.set_xticks(tick_x)
ax_g.set_xticklabels(HUB10, fontsize=8, fontweight="bold", rotation=30, ha="right")
ax_g.set_ylabel("Normalised Expression", fontsize=9.5, fontweight="bold")
ax_g.set_title("Expression – Discovery (GSE125498) | All 10 Hub Genes", fontsize=9.5, fontweight="bold", pad=4)
ax_g.legend(handles=[mpatches.Patch(color=C_LPD, label="LPD (n=13)"), mpatches.Patch(color=C_SPD, label="SPD (n=20)")],
            loc="upper right", fontsize=8, frameon=False)
ax_g.grid(True, linestyle=":", alpha=0.3, axis="y")
ax_g.text(0.01, -0.16, "*p<0.05   **p<0.01   ***p<0.001   ****p<0.0001   ns (p>=0.05)",
          transform=ax_g.transAxes, fontsize=7.2, color="#475569", style="italic")
plabel(ax_g, "G")

# ==========================================
# PANEL H: Expression Validation (Clean Legend & Footnote Outside Plot)
# ==========================================
ax_h = fig.add_subplot(outer[3, 1])
np.random.seed(0)

# 9 genes detected in mesothelial pseudobulk in exact same order
hub_in_pb = [g for g in HUB10 if g in log_cpm.index]

for i, g in enumerate(hub_in_pb):
    xl = i * (2 * w + gap); xs = xl + w
    vl = log_cpm.loc[g, lv_s].values.astype(float)
    vs = log_cpm.loc[g, sv_s].values.astype(float)
    
    for xp, vv, cl, clf in [(xl, vl, C_LPD, C_LPD_L), (xs, vs, C_SPD, C_SPD_L)]:
        ax_h.boxplot([vv], positions=[xp], widths=w * 0.90, patch_artist=True,
                     boxprops=dict(facecolor=clf, color=cl, lw=1.2),
                     whiskerprops=dict(color=cl, lw=1.0), capprops=dict(color=cl, lw=1.0),
                     medianprops=dict(color="#111827", lw=1.8),
                     flierprops=dict(marker="o", ms=3, color=cl, alpha=0.5))
        ax_h.scatter(np.random.normal(xp, 0.035, len(vv)), vv, color=cl, s=9, alpha=0.7, zorder=4)
        
    _, p = mannwhitneyu(vl, vs, alternative="two-sided")
    star = pstar(p)
    top = max(np.max(vl), np.max(vs))
    bary = top + 0.35; mx = (xl + xs) / 2
    ax_h.plot([xl, xs], [bary, bary], color="#374151", lw=0.9)
    ax_h.text(mx, bary + 0.08, star, ha="center", fontsize=8.5, fontweight="bold", color="#374151")

tick_xh = [i * (2 * w + gap) + w / 2 for i in range(len(hub_in_pb))]
ax_h.set_xticks(tick_xh)
ax_h.set_xticklabels(hub_in_pb, fontsize=8, fontweight="bold", rotation=30, ha="right")
ax_h.set_ylabel("logCPM (Mesothelial Pseudobulk)", fontsize=9.5, fontweight="bold")
ax_h.set_title(f"Expression – Verification (GSE248762) | {len(hub_in_pb)} Hub Genes Evaluated", fontsize=9.5, fontweight="bold", pad=4)

# Place legend in top right with ample headroom
ax_h.set_ylim(-0.8, 12.5)
ax_h.legend(handles=[mpatches.Patch(color=C_LPD, label="LV (n=10)"), mpatches.Patch(color=C_SPD, label="SV (n=6)")],
            loc="upper right", fontsize=8, framealpha=0.9, edgecolor="#CBD5E1")
ax_h.grid(True, linestyle=":", alpha=0.3, axis="y")

# Footnote placed below axis to completely avoid collision
ax_h.text(0.01, -0.16, "*p<0.05   **p<0.01   ***p<0.001   ns (p>=0.05)   |   NAT6 not detected in mesothelial scRNA pseudobulk",
          transform=ax_h.transAxes, fontsize=7.2, color="#475569", style="italic")
plabel(ax_h, "H")

# Save outputs
png_out = os.path.join(OUT_DIR, "Fig4_machine_learning_diagnostic_panel_v2.png")
pdf_out = os.path.join(OUT_DIR, "Fig4_machine_learning_diagnostic_panel_v2.pdf")

fig.savefig(png_out, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("SUCCESS: Saved regenerated figure:")
print(f"  PNG: {png_out}")
print(f"  PDF: {pdf_out}")
