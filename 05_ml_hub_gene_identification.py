"""
==============================================================================
SCRIPT 05: MACHINE LEARNING FEATURE SELECTION & HUB GENE IDENTIFICATION
Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
Models: LASSO, SVM-RFE, Random Forest, XGBoost
Input: 15 Convergent Causal ECM Genes across GSE62928 (EPS Fibrosis vs Control)
Consensus Rule: Selected by >= 2 of the 4 Machine Learning models
==============================================================================
"""

import os
import GEOparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

# Set random seed
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# Create output directories
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)
os.makedirs("data", exist_ok=True)

print("=" * 70)
print("STEP 1: LOADING 15 CONVERGENT CAUSAL ECM GENES & GSE62928 EXPRESSION")
print("=" * 70)

# 1. Load 15 convergent causal genes
df_15 = pd.read_csv("convergent_15_causal_genes.csv")
genes_15 = df_15["GeneSymbol"].tolist()
print(f"Loaded {len(genes_15)} candidate genes: {', '.join(genes_15)}")

# 2. Load GSE62928 SOFT data
gse = GEOparse.get_GEO(filepath="data/GSE62928_family.soft.gz")
piv = gse.pivot_samples("VALUE")
print(f"Platform probe matrix shape: {piv.shape}")

# 3. Best probe mapping from limma top table
df_top = pd.read_csv("GSE62928.top.table.tsv", sep="\t")

def parse_syms(s):
    if pd.isna(s) or str(s).strip() in ["", "---", "NAN"]:
        return []
    return [g.strip().upper() for g in str(s).replace(";", "///").split("///") if g.strip()]

rows = []
for idx, r in df_top.iterrows():
    for g in parse_syms(r["Gene.symbol"]):
        if g in genes_15:
            rows.append({
                "Gene": g,
                "ID": r["ID"],
                "P.Value": r["P.Value"],
                "logFC": r["logFC"]
            })

df_map = pd.DataFrame(rows).sort_values("P.Value").drop_duplicates("Gene")
print(f"Mapped {len(df_map)} genes to high-affinity probes.")

# Extract expression for the 15 genes
expr_df = piv.loc[df_map["ID"].tolist()].copy()
expr_df.index = df_map["Gene"].tolist()

# Define sample classes
# GSM1536406 - GSM1536409: EPS (Peritoneal Fibrosis / Sclerosis, Class 1)
# GSM1536410 - GSM1536413: Controls (PD + Uremic, Class 0)
sample_classes = {
    "GSM1536406": 1, "GSM1536407": 1, "GSM1536408": 1, "GSM1536409": 1,
    "GSM1536410": 0, "GSM1536411": 0, "GSM1536412": 0, "GSM1536413": 0
}

X_df = expr_df[list(sample_classes.keys())].T
y = np.array([sample_classes[s] for s in X_df.index])
features = X_df.columns.tolist()

print(f"Feature matrix X: {X_df.shape} (8 samples x {len(features)} genes)")
print(f"Class distribution: {np.sum(y == 1)} Peritoneal Fibrosis (EPS) vs {np.sum(y == 0)} Controls")

# Standardize feature matrix
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_df)
X_scaled_df = pd.DataFrame(X_scaled, index=X_df.index, columns=features)

# ==============================================================================
# STEP 2: RUN 4 MACHINE LEARNING MODELS
# ==============================================================================
print("\n" + "=" * 70)
print("STEP 2: RUNNING 4 MACHINE LEARNING FEATURE SELECTION MODELS")
print("=" * 70)

# ------------------------------------------------------------------------------
# Model 1: LASSO Logistic Regression (L1-regularization)
# ------------------------------------------------------------------------------
print("\n[Model 1/4] Running LASSO Logistic Regression (L1 Penalty)...")
lasso_cv = LogisticRegressionCV(
    Cs=20, cv=4, penalty="l1", solver="liblinear",
    random_state=RANDOM_STATE, scoring="accuracy"
)
lasso_cv.fit(X_scaled, y)
best_C = lasso_cv.C_[0]
lasso_coefs = lasso_cv.coef_[0]
lasso_selected = np.abs(lasso_coefs) > 1e-4

# Fallback check: if penalty was too aggressive, calibrate with soft threshold
if np.sum(lasso_selected) < 2:
    clf_lasso = LogisticRegression(penalty="l1", C=0.5, solver="liblinear", random_state=RANDOM_STATE)
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
# Select top 7 features (~half of candidate panel)
rfe = RFE(estimator=svm, n_features_to_select=7, step=1)
rfe.fit(X_scaled, y)
svm_selected = rfe.support_
# RFE ranking: 1 is top selected, convert to score (inverse rank)
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
    
    # Gene metadata
    g_info = df_15[df_15["GeneSymbol"] == gene].iloc[0]
    
    summary_rows.append({
        "GeneSymbol": gene,
        "Direction": g_info["Direction"],
        "Peritoneal_logFC": g_info["logFC"],
        "Peritoneal_Pval": g_info["P.Value"],
        "TWMR_Beta": g_info["MR_Beta"],
        "TWMR_Pval": g_info["MR_Pvalue"],
        "TWMR_FDR": g_info["MR_FDR"],
        "LASSO_Selected": l_sel,
        "LASSO_Score": lasso_score_dict[gene],
        "SVMRFE_Selected": s_sel,
        "SVMRFE_Score": svm_score_dict[gene],
        "RandomForest_Selected": r_sel,
        "RandomForest_Score": rf_score_dict[gene],
        "XGBoost_Selected": x_sel,
        "XGBoost_Score": xgb_score_dict[gene],
        "Consensus_Votes": votes,
        "Is_Hub_Gene": "YES" if votes >= 2 else "NO"
    })

df_consensus = pd.DataFrame(summary_rows).sort_values(
    ["Consensus_Votes", "Peritoneal_Pval"], ascending=[False, True]
)

hub_genes_df = df_consensus[df_consensus["Is_Hub_Gene"] == "YES"]
print(f"\nIDENTIFIED {len(hub_genes_df)} CONSENSUS HUB GENES (Votes >= 2):")
for idx, r in hub_genes_df.iterrows():
    print(f"  * {r['GeneSymbol']:<10} | Votes: {r['Consensus_Votes']}/4 | Direction: {r['Direction']:<4} | logFC: {r['Peritoneal_logFC']:+.2f} | MR P: {r['TWMR_Pval']:.3e} | MR FDR: {r['TWMR_FDR']:.3e}")

# Save CSVs
df_consensus.to_csv("results/tables/ML_hub_genes_consensus_results.csv", index=False)
hub_genes_df.to_csv("results/tables/ML_hub_genes_final_list.csv", index=False)
print("\nSaved consensus results to results/tables/ML_hub_genes_consensus_results.csv")

# ==============================================================================
# STEP 4: PUBLICATION-GRADE VISUALIZATION
# ==============================================================================
print("\n" + "=" * 70)
print("STEP 4: GENERATING PUBLICATION-GRADE VISUALIZATIONS")
print("=" * 70)

# ------------------------------------------------------------------------------
# Plot 1: Consensus Votes Bar Chart
# ------------------------------------------------------------------------------
plt.figure(figsize=(10, 6), facecolor="#F8FAFC")
colors = ["#DC2626" if v >= 2 else "#94A3B8" for v in df_consensus["Consensus_Votes"]]
bars = plt.bar(df_consensus["GeneSymbol"], df_consensus["Consensus_Votes"], color=colors, edgecolor="#1E293B", width=0.65)
plt.axhline(2, color="#B91C1C", linestyle="--", linewidth=1.5, label="Hub Selection Cutoff (>= 2 Models)")

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.08, int(yval), ha="center", va="bottom", fontsize=11, fontweight="bold")

plt.title("Machine Learning Consensus Feature Selection (4 Models)\nHub Biomarkers Identified by >= 2 Algorithms", fontsize=13, fontweight="bold", pad=15)
plt.xlabel("Candidate Causal ECM Genes", fontsize=11, fontweight="bold")
plt.ylabel("Consensus Votes (out of 4)", fontsize=11, fontweight="bold")
plt.ylim(0, 4.6)
plt.xticks(rotation=45, ha="right", fontsize=10, fontweight="medium")
plt.legend(loc="upper right", frameon=True, facecolor="#EFF6FF")
plt.grid(axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig("results/figures/ML_01_consensus_votes_barchart.png", dpi=300)
plt.close()
print("Saved Plot 1: results/figures/ML_01_consensus_votes_barchart.png")

# ------------------------------------------------------------------------------
# Plot 2: Model Selection Heatmap
# ------------------------------------------------------------------------------
heatmap_data = df_consensus.set_index("GeneSymbol")[
    ["LASSO_Selected", "SVMRFE_Selected", "RandomForest_Selected", "XGBoost_Selected"]
].T
heatmap_data.index = ["LASSO", "SVM-RFE", "Random Forest", "XGBoost"]

plt.figure(figsize=(12, 4.5), facecolor="#F8FAFC")
sns.heatmap(
    heatmap_data, cmap=["#F1F5F9", "#2563EB"], cbar=False,
    linewidths=1.2, linecolor="#CBD5E1", annot=True, fmt="d",
    annot_kws={"fontsize": 11, "fontweight": "bold"}
)
plt.title("Algorithm-Specific Selection Matrix Across 15 Convergent Causal ECM Genes", fontsize=13, fontweight="bold", pad=15)
plt.xlabel("Gene Symbol", fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("results/figures/ML_02_model_selection_heatmap.png", dpi=300)
plt.close()
print("Saved Plot 2: results/figures/ML_02_model_selection_heatmap.png")

# ------------------------------------------------------------------------------
# Plot 3: 2x2 Per-Model Feature Importance
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor="#F8FAFC")

# Subplot A: LASSO
order_lasso = df_consensus.sort_values("LASSO_Score", ascending=True)
axes[0, 0].barh(order_lasso["GeneSymbol"], order_lasso["LASSO_Score"], color="#3B82F6", edgecolor="#1E293B")
axes[0, 0].set_title("A. LASSO Logistic Regression (|Beta|)", fontsize=11, fontweight="bold")
axes[0, 0].set_xlabel("Absolute Coefficient Magnitude")

# Subplot B: SVM-RFE
order_svm = df_consensus.sort_values("SVMRFE_Score", ascending=True)
axes[0, 1].barh(order_svm["GeneSymbol"], order_svm["SVMRFE_Score"], color="#10B981", edgecolor="#1E293B")
axes[0, 1].set_title("B. SVM-RFE (Relative Ranking Score)", fontsize=11, fontweight="bold")
axes[0, 1].set_xlabel("1 / RFE Rank")

# Subplot C: Random Forest
order_rf = df_consensus.sort_values("RandomForest_Score", ascending=True)
axes[1, 0].barh(order_rf["GeneSymbol"], order_rf["RandomForest_Score"], color="#F59E0B", edgecolor="#1E293B")
axes[1, 0].axvline(rf_threshold, color="#B91C1C", linestyle="--", linewidth=1.2, label=f"Mean ({rf_threshold:.3f})")
axes[1, 0].set_title("C. Random Forest (Gini Impurity Importance)", fontsize=11, fontweight="bold")
axes[1, 0].set_xlabel("Gini Importance")
axes[1, 0].legend(loc="lower right")

# Subplot D: XGBoost
order_xgb = df_consensus.sort_values("XGBoost_Score", ascending=True)
axes[1, 1].barh(order_xgb["GeneSymbol"], order_xgb["XGBoost_Score"], color="#8B5CF6", edgecolor="#1E293B")
axes[1, 1].axvline(xgb_threshold, color="#B91C1C", linestyle="--", linewidth=1.2, label=f"Mean ({xgb_threshold:.3f})")
axes[1, 1].set_title("D. XGBoost (Feature Weight Importance)", fontsize=11, fontweight="bold")
axes[1, 1].set_xlabel("Feature Importance")
axes[1, 1].legend(loc="lower right")

for ax in axes.flat:
    ax.grid(axis="x", linestyle=":", alpha=0.6)

plt.suptitle("Comparative Feature Importance Across 4 Machine Learning Architectures", fontsize=14, fontweight="bold", y=0.99)
plt.tight_layout()
plt.savefig("results/figures/ML_03_per_model_importance_2x2.png", dpi=300)
plt.close()
print("Saved Plot 3: results/figures/ML_03_per_model_importance_2x2.png")

# ------------------------------------------------------------------------------
# Plot 4: Clustered Expression Heatmap of Hub Genes
# ------------------------------------------------------------------------------
hub_gene_names = hub_genes_df["GeneSymbol"].tolist()
hub_expr = expr_df.loc[hub_gene_names].copy()

# Sample annotations
sample_labels = ["EPS1 (Fibrosis)", "EPS2 (Fibrosis)", "EPS3 (Fibrosis)", "EPS4 (Fibrosis)",
                 "PD1 (Control)", "PD2 (Control)", "UREMIC1 (Control)", "UREMIC2 (Control)"]
hub_expr.columns = sample_labels

# Z-score row normalization
hub_expr_z = hub_expr.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

cg = sns.clustermap(
    hub_expr_z, cmap="vlag", figsize=(9, 7),
    col_cluster=False, row_cluster=True,
    dendrogram_ratio=(0.2, 0.02),
    cbar_kws={"label": "Expression Z-Score"},
    linewidths=0.5, linecolor="#E2E8F0"
)
cg.ax_heatmap.set_title("Supervised Expression Profile of Consensus Hub Biomarkers\n(Peritoneal Fibrosis EPS vs Non-Fibrotic Controls)", fontsize=12, fontweight="bold", pad=20)
plt.savefig("results/figures/ML_04_hub_genes_expression_heatmap.png", dpi=300)
plt.close()
print("Saved Plot 4: results/figures/ML_04_hub_genes_expression_heatmap.png")

print("\n" + "=" * 70)
print("MACHINE LEARNING PIPELINE COMPLETED SUCCESSFULLY!")
print("=" * 70)
