"""
==============================================================================
SCRIPT 05b: MACHINE LEARNING CONSENSUS HUB GENE IDENTIFICATION (WGCNA-ECM)
Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
Purpose: Feature selection on convergent WGCNA ∩ 71 ECM-DEGs using 4 consensus
         ML models: LASSO, SVM-RFE, Random Forest, and XGBoost
==============================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

# Set seed
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

print("=" * 70)
print("STEP 1: LOADING CONVERGENT WGCNA-ECM GENES & GSE62928 EXPRESSION")
print("=" * 70)

# 1. Load convergent WGCNA-ECM candidate genes
wgcna_ecm_file = "convergent_WGCNA_ECM_genes.csv"
if not os.path.exists(wgcna_ecm_file):
    wgcna_ecm_file = "results/tables/convergent_WGCNA_ECM_genes.csv"

df_candidates = pd.read_csv(wgcna_ecm_file)
candidate_genes = df_candidates["gene_symbol"].tolist()
print(f"Loaded {len(candidate_genes)} convergent WGCNA-ECM genes.")

# 2. Load full GSE62928 expression matrix and metadata
expr_file = "results/tables/GSE62928_full_expression_matrix.csv"
meta_file = "results/tables/GSE62928_sample_metadata.csv"

if not os.path.exists(expr_file) or not os.path.exists(meta_file):
    raise FileNotFoundError("Missing expression matrix or metadata! Please run 00_fetch_gse62928_matrix.R first.")

df_expr = pd.read_csv(expr_file, index_col=0)
df_meta = pd.read_csv(meta_file)

# Match available candidate genes in expression matrix
available_genes = [g for g in candidate_genes if g in df_expr.index]
print(f"Matched {len(available_genes)} / {len(candidate_genes)} candidate genes in full expression matrix.")

# Prepare sample X and y
sample_ids = df_meta["sample_id"].tolist()
# Extract binary label y
if "binary_numeric" in df_meta.columns:
    y = df_meta["binary_numeric"].values.astype(int)
elif "group_binary_trait" in df_meta.columns:
    y = np.array([1 if "case" in str(g).lower() or "eps" in str(g).lower() else 0 for g in df_meta["group_binary_trait"]])
elif "group" in df_meta.columns:
    y = np.array([1 if "case" in str(g).lower() or "eps" in str(g).lower() else 0 for g in df_meta["group"]])
else:
    y = np.array([1 if "eps" in str(s).lower() else 0 for s in df_meta["sample_id"]])

X_df = df_expr.loc[available_genes, sample_ids].T
features = available_genes

print(f"Feature matrix shape: {X_df.shape} ({len(sample_ids)} samples x {len(features)} genes)")
print(f"Class distribution: {np.sum(y == 1)} Peritoneal Fibrosis (Case) vs {np.sum(y == 0)} Controls")

# Standardize feature matrix
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_df)
X_scaled_df = pd.DataFrame(X_scaled, index=X_df.index, columns=features)

# ==============================================================================
# STEP 2: RUN 4 MACHINE LEARNING FEATURE SELECTION MODELS
# ==============================================================================
print("\n" + "=" * 70)
print("STEP 2: RUNNING 4 MACHINE LEARNING FEATURE SELECTION MODELS")
print("=" * 70)

# ------------------------------------------------------------------------------
# Model 1: LASSO Logistic Regression (L1-regularization)
# ------------------------------------------------------------------------------
print("\n[Model 1/4] Running LASSO Logistic Regression (L1 Penalty)...")
try:
    lasso_cv = LogisticRegressionCV(
        Cs=20, cv=4, penalty="l1", solver="liblinear",
        random_state=RANDOM_STATE, scoring="accuracy"
    )
    lasso_cv.fit(X_scaled, y)
    lasso_coefs = lasso_cv.coef_[0]
except Exception:
    clf_lasso = LogisticRegression(penalty="l1", C=0.5, solver="liblinear", random_state=RANDOM_STATE)
    clf_lasso.fit(X_scaled, y)
    lasso_coefs = clf_lasso.coef_[0]

lasso_selected = np.abs(lasso_coefs) > 1e-4

if np.sum(lasso_selected) < 2:
    clf_lasso = LogisticRegression(penalty="l1", C=0.8, solver="liblinear", random_state=RANDOM_STATE)
    clf_lasso.fit(X_scaled, y)
    lasso_coefs = clf_lasso.coef_[0]
    lasso_selected = np.abs(lasso_coefs) > 1e-4

lasso_score_dict = dict(zip(features, np.abs(lasso_coefs)))
print(f"  LASSO selected {np.sum(lasso_selected)} genes: {[f for f, s in zip(features, lasso_selected) if s]}")

# ------------------------------------------------------------------------------
# Model 2: SVM-RFE (Support Vector Machine Recursive Feature Elimination)
# ------------------------------------------------------------------------------
print("\n[Model 2/4] Running SVM-RFE (Linear Kernel)...")
svm = SVC(kernel="linear", C=1.0, random_state=RANDOM_STATE)
n_rfe_select = max(5, int(len(features) * 0.35))
rfe = RFE(estimator=svm, n_features_to_select=n_rfe_select, step=1)
rfe.fit(X_scaled, y)
svm_selected = rfe.support_
svm_score_dict = dict(zip(features, 1.0 / rfe.ranking_))
print(f"  SVM-RFE selected {np.sum(svm_selected)} genes: {[f for f, s in zip(features, svm_selected) if s]}")

# ------------------------------------------------------------------------------
# Model 3: Random Forest Classifier (Gini Feature Importance)
# ------------------------------------------------------------------------------
print("\n[Model 3/4] Running Random Forest Classifier...")
rf = RandomForestClassifier(n_estimators=500, max_depth=3, random_state=RANDOM_STATE)
rf.fit(X_scaled, y)
rf_importances = rf.feature_importances_
rf_threshold = np.mean(rf_importances)
rf_selected = rf_importances >= rf_threshold
rf_score_dict = dict(zip(features, rf_importances))
print(f"  Random Forest threshold: {rf_threshold:.4f}")
print(f"  Random Forest selected {np.sum(rf_selected)} genes: {[f for f, s in zip(features, rf_selected) if s]}")

# ------------------------------------------------------------------------------
# Model 4: XGBoost Classifier (Gradient Boosted Trees)
# ------------------------------------------------------------------------------
print("\n[Model 4/4] Running XGBoost Classifier...")
xgb = XGBClassifier(
    n_estimators=100, max_depth=2, learning_rate=0.08,
    random_state=RANDOM_STATE, eval_metric="logloss"
)
xgb.fit(X_scaled, y)
xgb_importances = xgb.feature_importances_
xgb_threshold = np.mean(xgb_importances)
xgb_selected = xgb_importances >= xgb_threshold
xgb_score_dict = dict(zip(features, xgb_importances))
print(f"  XGBoost threshold: {xgb_threshold:.4f}")
print(f"  XGBoost selected {np.sum(xgb_selected)} genes: {[f for f, s in zip(features, xgb_selected) if s]}")

# ==============================================================================
# STEP 3: CONSENSUS INTEGRATION & HUB GENE DISCOVERY
# ==============================================================================
print("\n" + "=" * 70)
print("STEP 3: CONSENSUS VOTING INTEGRATION (>= 2 OF 4 MODELS)")
print("=" * 70)

summary_rows = []
for gene in features:
    l_sel = 1 if gene in [f for f, s in zip(features, lasso_selected) if s] else 0
    s_sel = 1 if gene in [f for f, s in zip(features, svm_selected) if s] else 0
    r_sel = 1 if gene in [f for f, s in zip(features, rf_selected) if s] else 0
    x_sel = 1 if gene in [f for f, s in zip(features, xgb_selected) if s] else 0
    votes = l_sel + s_sel + r_sel + x_sel
    
    g_info = df_candidates[df_candidates["gene_symbol"] == gene].iloc[0]
    
    # Check which models selected this gene
    sel_models = []
    if l_sel: sel_models.append("LASSO")
    if s_sel: sel_models.append("SVM-RFE")
    if r_sel: sel_models.append("RF")
    if x_sel: sel_models.append("XGBoost")
    
    summary_rows.append({
        "Gene_Symbol": gene,
        "Module_Color": g_info["module_color"],
        "WGCNA_MM": g_info["MM"],
        "WGCNA_GS": g_info["GS"],
        "Peritoneal_logFC": g_info["logFC"],
        "Peritoneal_Pval": g_info["P.Value"],
        "Peritoneal_adjPval": g_info["adj.P.Val"],
        "Matrisome_Division": g_info["Matrisome_Division"],
        "Matrisome_Category": g_info["Matrisome_Category"],
        "LASSO_Selected": l_sel,
        "LASSO_Score": lasso_score_dict[gene],
        "SVMRFE_Selected": s_sel,
        "SVMRFE_Score": svm_score_dict[gene],
        "RandomForest_Selected": r_sel,
        "RandomForest_Score": rf_score_dict[gene],
        "XGBoost_Selected": x_sel,
        "XGBoost_Score": xgb_score_dict[gene],
        "Votes": votes,
        "Selecting_Models": " + ".join(sel_models) if sel_models else "None",
        "Is_Hub_Gene": "YES" if votes >= 2 else "NO"
    })

df_consensus = pd.DataFrame(summary_rows).sort_values(
    ["Votes", "Peritoneal_Pval"], ascending=[False, True]
)

hub_genes_df = df_consensus[df_consensus["Is_Hub_Gene"] == "YES"].copy()
print(f"\nIDENTIFIED {len(hub_genes_df)} CONSENSUS HUB GENES (Votes >= 2):")
for idx, r in hub_genes_df.iterrows():
    print(f"  * {r['Gene_Symbol']:<10} | Votes: {r['Votes']}/4 ({r['Selecting_Models']}) | logFC: {r['Peritoneal_logFC']:+.2f} | MM: {r['WGCNA_MM']:.3f} | GS: {r['WGCNA_GS']:.3f}")

# Save CSVs
df_consensus.to_csv("results/tables/ML_hub_genes_from_WGCNA_ECM_all_results.csv", index=False)
hub_genes_df.to_csv("results/tables/ML_hub_genes_from_WGCNA_ECM.csv", index=False)
print("\nSaved final Hub Genes list to: results/tables/ML_hub_genes_from_WGCNA_ECM.csv")

# ==============================================================================
# STEP 4: VISUALIZATION SUITE
# ==============================================================================
print("\n" + "=" * 70)
print("STEP 4: GENERATING PUBLICATION-GRADE VISUALIZATIONS")
print("=" * 70)

# ------------------------------------------------------------------------------
# Plot 1: Consensus Votes Bar Chart
# ------------------------------------------------------------------------------
plt.figure(figsize=(12, 6.5), facecolor="#F8FAFC")
colors = ["#DC2626" if v >= 2 else "#94A3B8" for v in df_consensus["Votes"]]
bars = plt.bar(df_consensus["Gene_Symbol"], df_consensus["Votes"], color=colors, edgecolor="#1E293B", width=0.65)
plt.axhline(2, color="#B91C1C", linestyle="--", linewidth=1.5, label="Hub Selection Cutoff (>= 2 Models)")

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.08, int(yval), ha="center", va="bottom", fontsize=9.5, fontweight="bold")

plt.title("WGCNA-ECM Machine Learning Consensus Feature Selection (4 Models)\nHub Biomarkers Identified by >= 2 Algorithms", fontsize=13, fontweight="bold", pad=15)
plt.xlabel("Convergent WGCNA ∩ ECM-DEGs", fontsize=11, fontweight="bold")
plt.ylabel("Consensus Votes (out of 4)", fontsize=11, fontweight="bold")
plt.ylim(0, 4.6)
plt.xticks(rotation=60, ha="right", fontsize=9.5, fontweight="medium")
plt.legend(loc="upper right", frameon=True, facecolor="#EFF6FF")
plt.grid(axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig("results/figures/WGCNA_ML_01_consensus_votes_barchart.png", dpi=300)
plt.close()
print("Saved Plot 1: results/figures/WGCNA_ML_01_consensus_votes_barchart.png")

# ------------------------------------------------------------------------------
# Plot 2: Model Selection Heatmap
# ------------------------------------------------------------------------------
df_binary = df_consensus[["Gene_Symbol", "LASSO_Selected", "SVMRFE_Selected", "RandomForest_Selected", "XGBoost_Selected"]].set_index("Gene_Symbol")
df_binary.columns = ["LASSO", "SVM-RFE", "Random Forest", "XGBoost"]

plt.figure(figsize=(9, 11), facecolor="#F8FAFC")
sns.heatmap(
    df_binary, cmap=["#F1F5F9", "#DC2626"], cbar=False,
    linewidths=0.5, linecolor="#CBD5E1", annot=True, fmt="d",
    annot_kws={"fontsize": 9, "fontweight": "bold"}
)
plt.title("Binary Feature Selection Matrix Across 4 ML Models (WGCNA-ECM)\n(1 = Selected, 0 = Excluded)", fontsize=12, fontweight="bold", pad=15)
plt.xlabel("Machine Learning Algorithm", fontsize=11, fontweight="bold")
plt.ylabel("Convergent WGCNA ∩ ECM-DEGs", fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("results/figures/WGCNA_ML_02_model_selection_heatmap.png", dpi=300)
plt.close()
print("Saved Plot 2: results/figures/WGCNA_ML_02_model_selection_heatmap.png")

# ------------------------------------------------------------------------------
# Plot 3: 2x2 Feature Importance Panel
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(15, 12), facecolor="#F8FAFC")
top_k = min(15, len(features))

# 1. LASSO
df_lasso = df_consensus.sort_values("LASSO_Score", ascending=False).head(top_k)
axes[0, 0].barh(df_lasso["Gene_Symbol"][::-1], df_lasso["LASSO_Score"][::-1], color="#2563EB", edgecolor="#1E293B", height=0.65)
axes[0, 0].set_title("LASSO Feature Coefficients (|Beta|)", fontsize=11, fontweight="bold")
axes[0, 0].grid(axis="x", linestyle=":", alpha=0.6)

# 2. SVM-RFE
df_svm = df_consensus.sort_values("SVMRFE_Score", ascending=False).head(top_k)
axes[0, 1].barh(df_svm["Gene_Symbol"][::-1], df_svm["SVMRFE_Score"][::-1], color="#059669", edgecolor="#1E293B", height=0.65)
axes[0, 1].set_title("SVM-RFE Inverse Rank Score", fontsize=11, fontweight="bold")
axes[0, 1].grid(axis="x", linestyle=":", alpha=0.6)

# 3. Random Forest
df_rf = df_consensus.sort_values("RandomForest_Score", ascending=False).head(top_k)
axes[1, 0].barh(df_rf["Gene_Symbol"][::-1], df_rf["RandomForest_Score"][::-1], color="#D97706", edgecolor="#1E293B", height=0.65)
axes[1, 0].axvline(rf_threshold, color="#B91C1C", linestyle="--", label=f"Mean Cutoff ({rf_threshold:.4f})")
axes[1, 0].set_title("Random Forest Gini Importance", fontsize=11, fontweight="bold")
axes[1, 0].legend(loc="lower right")
axes[1, 0].grid(axis="x", linestyle=":", alpha=0.6)

# 4. XGBoost
df_xgb = df_consensus.sort_values("XGBoost_Score", ascending=False).head(top_k)
axes[1, 1].barh(df_xgb["Gene_Symbol"][::-1], df_xgb["XGBoost_Score"][::-1], color="#7C3AED", edgecolor="#1E293B", height=0.65)
axes[1, 1].axvline(xgb_threshold, color="#B91C1C", linestyle="--", label=f"Mean Cutoff ({xgb_threshold:.4f})")
axes[1, 1].set_title("XGBoost Gain Feature Importance", fontsize=11, fontweight="bold")
axes[1, 1].legend(loc="lower right")
axes[1, 1].grid(axis="x", linestyle=":", alpha=0.6)

plt.suptitle("Comparative Feature Importance Across 4 Machine Learning Frameworks (WGCNA-ECM)", fontsize=13, fontweight="bold", y=0.98)
plt.tight_layout()
plt.savefig("results/figures/WGCNA_ML_03_per_model_importance_2x2.png", dpi=300)
plt.close()
print("Saved Plot 3: results/figures/WGCNA_ML_03_per_model_importance_2x2.png")

# ------------------------------------------------------------------------------
# Plot 4: Hub Genes Expression Heatmap
# ------------------------------------------------------------------------------
hub_symbols = hub_genes_df["Gene_Symbol"].tolist()
hub_expr = X_df[hub_symbols].T

# Z-score normalize
hub_expr_z = hub_expr.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

col_colors = ["#EF4444" if yi == 1 else "#3B82F6" for yi in y]
col_colors_series = pd.Series(col_colors, index=sample_ids, name="Group")

cg = sns.clustermap(
    hub_expr_z, cmap="vlag", figsize=(10, 8),
    col_cluster=True, row_cluster=True,
    col_colors=col_colors_series,
    dendrogram_ratio=(0.18, 0.12),
    cbar_kws={"label": "Expression Z-Score"},
    linewidths=0.5, linecolor="#CBD5E1"
)
cg.figure.subplots_adjust(top=0.90)
cg.figure.suptitle("Expression Profile of Consensus WGCNA-ECM Hub Genes in GSE62928\n(Blue = Control, Red = Peritoneal Fibrosis Case)",
                  fontsize=12, fontweight="bold", y=0.96)
plt.savefig("results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png", dpi=300)
plt.close()
print("Saved Plot 4: results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png")

print("\n" + "=" * 70)
print("TASK 3 COMPLETED SUCCESSFULLY!")
print("=" * 70)
