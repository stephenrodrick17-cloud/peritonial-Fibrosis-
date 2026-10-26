"""
Run Machine Learning feature selection (LASSO, SVM-RFE, Random Forest)
directly on the 43 module-assigned matrisome DEGs.
Cohort: GSE125498 (33 samples: 13 LPD vs 20 SPD)
Output:
- Detailed selection tables per algorithm
- Consensus hub genes (intersection >= 2 models and strict 3 models)
- JSON summary and comparison to 38-gene model
"""

import json
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

SEED = 42
OUT_DIR = "results/revision/ml_43_matrisome_wgcna"
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load Data
expression = pd.read_csv("results/tables/GSE125498_full_expression_matrix.csv")
metadata = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
intersect = pd.read_csv("results/revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv")
all_results = pd.read_csv("results/tables/GSE125498_all_results.csv")

gcol = "Gene.symbol" if "Gene.symbol" in intersect.columns else "Annotated Gene"
genes_43 = sorted(intersect[intersect["Module"].notna()][gcol].tolist())
print(f"Loaded {len(genes_43)} input features (43 module-assigned matrisome DEGs).")

sample_ids = metadata["sample_id"].astype(str).tolist()
expr = expression.set_index("Gene")

x_raw = expr.loc[genes_43, sample_ids].to_numpy(dtype=float).T
y = metadata.set_index("sample_id").loc[sample_ids, "stage_binary"].to_numpy(dtype=int)

x = StandardScaler().fit_transform(x_raw)
n_select = 10  # top 10 features for SVM and RF

# 2. LASSO
lasso = LogisticRegressionCV(
    Cs=[0.01, 0.05, 0.1, 0.5, 1, 5, 10],
    cv=5,
    penalty="l1",
    solver="liblinear",
    max_iter=5000,
    random_state=SEED,
)
lasso.fit(x, y)
lasso_coef = np.abs(lasso.coef_[0])
lasso_selected_idx = np.flatnonzero(lasso_coef > 1e-4)
if lasso_selected_idx.size == 0:
    lasso_selected_idx = np.argsort(-lasso_coef)[:10]

lasso_set = {genes_43[i] for i in lasso_selected_idx}

# 3. SVM-RFE
svm_rfe = RFE(
    estimator=SVC(kernel="linear", random_state=SEED),
    n_features_to_select=n_select,
    step=1,
)
svm_rfe.fit(x, y)
svm_set = {genes_43[i] for i in np.flatnonzero(svm_rfe.support_)}

# 4. Random Forest
rf = RandomForestClassifier(
    n_estimators=500,
    random_state=SEED,
    class_weight="balanced",
)
rf.fit(x, y)
rf_importance = rf.feature_importances_
rf_threshold = np.sort(rf_importance)[::-1][n_select - 1]
rf_set = {genes_43[i] for i in np.flatnonzero(rf_importance >= rf_threshold)}

# 5. Assemble Results
module_lookup = intersect.set_index(gcol)["Module"].to_dict()
deg_lookup = intersect.set_index(gcol)["DEG_Status"].to_dict()
stats_lookup = all_results.set_index("Gene")

rows = []
for i, gene in enumerate(genes_43):
    in_lasso = int(gene in lasso_set)
    in_svm = int(gene in svm_set)
    in_rf = int(gene in rf_set)
    vote_count = in_lasso + in_svm + in_rf
    
    stat_row = stats_lookup.loc[gene] if gene in stats_lookup.index else {}
    
    rows.append({
        "Gene": gene,
        "WGCNA_Module": module_lookup.get(gene, "Unknown"),
        "DEG_Status": deg_lookup.get(gene, "Unknown"),
        "LASSO_selected": in_lasso,
        "SVM_RFE_selected": in_svm,
        "RF_selected": in_rf,
        "Model_Vote_Count": vote_count,
        "Majority_Consensus": int(vote_count >= 2),
        "Strict_Consensus": int(vote_count == 3),
        "LASSO_abs_coef": float(lasso_coef[i]),
        "SVM_RFE_rank": int(svm_rfe.ranking_[i]),
        "RF_importance": float(rf_importance[i]),
        "logFC": float(stat_row.get("logFC", np.nan)),
        "adj.P.Val": float(stat_row.get("adj.P.Val", np.nan))
    })

res_df = pd.DataFrame(rows).sort_values(
    by=["Model_Vote_Count", "RF_importance", "LASSO_abs_coef"],
    ascending=[False, False, False]
)

csv_out = os.path.join(OUT_DIR, "ml_43_matrisome_results.csv")
res_df.to_csv(csv_out, index=False)

majority_genes = res_df[res_df["Majority_Consensus"] == 1]["Gene"].tolist()
strict_genes = res_df[res_df["Strict_Consensus"] == 1]["Gene"].tolist()

summary = {
    "dataset": "GSE125498",
    "total_samples": len(sample_ids),
    "input_features": len(genes_43),
    "algorithms": ["LASSO (L1)", "SVM-RFE", "Random Forest"],
    "counts_selected": {
        "LASSO": len(lasso_set),
        "SVM_RFE": len(svm_set),
        "Random_Forest": len(rf_set)
    },
    "majority_consensus_count": len(majority_genes),
    "majority_consensus_genes": majority_genes,
    "strict_consensus_count": len(strict_genes),
    "strict_consensus_genes": strict_genes,
    "output_file": csv_out
}

with open(os.path.join(OUT_DIR, "ml_43_matrisome_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

print("\n=== MACHINE LEARNING ON 43 GENES RESULTS ===")
print(f"LASSO selected: {len(lasso_set)} genes")
print(f"SVM-RFE selected: {len(svm_set)} genes")
print(f"Random Forest selected: {len(rf_set)} genes")
print(f"Majority Consensus (>= 2 models): {len(majority_genes)} genes -> {majority_genes}")
print(f"Strict Consensus (3 of 3 models): {len(strict_genes)} genes -> {strict_genes}")
print(f"Results saved to: {csv_out}")
