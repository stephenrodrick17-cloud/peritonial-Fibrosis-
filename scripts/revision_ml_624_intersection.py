import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


SEED = 42
ROOT = Path(__file__).resolve().parents[1]
REVISION_DIR = ROOT / "results" / "revision"
REVISION_DIR.mkdir(parents=True, exist_ok=True)

expression = pd.read_csv(ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv")
metadata = pd.read_csv(ROOT / "results" / "tables" / "GSE125498_sample_metadata.csv")
modules = pd.read_csv(REVISION_DIR / "wgcna_modules_dynamic_taskB.csv")
degs = pd.read_csv(ROOT / "results" / "tables" / "GSE125498_DEGs_filtered.csv")
all_results = pd.read_csv(ROOT / "results" / "tables" / "GSE125498_all_results.csv")

module_genes = set(
    modules.loc[modules["Module"].astype(str).str.lower() != "grey", "Gene"]
    .dropna()
    .astype(str)
)
deg_genes = set(degs["Gene"].dropna().astype(str))
features = sorted(module_genes & deg_genes)
if len(features) != 624:
    raise ValueError(f"Expected 624 shared WGCNA-module/DEG genes; found {len(features)}.")

sample_ids = metadata["sample_id"].astype(str).tolist()
if not set(sample_ids).issubset(expression.columns):
    raise ValueError("Metadata contains sample IDs missing from the expression matrix.")
if not set(features).issubset(set(expression["Gene"].astype(str))):
    raise ValueError("Some selected genes are missing from the expression matrix.")

expr = expression.set_index("Gene")
x_raw = expr.loc[features, sample_ids].to_numpy(dtype=float).T
y = metadata.set_index("sample_id").loc[sample_ids, "stage_binary"].to_numpy(dtype=int)
if np.unique(y).size != 2:
    raise ValueError("The discovery cohort must contain both stage groups.")

# This is feature selection on the full discovery cohort, not an estimate of
# out-of-sample performance. Predictive performance requires nested validation.
x = StandardScaler().fit_transform(x_raw)
n_select = max(2, min(10, int(np.ceil(len(features) * 0.35))))

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
lasso_selected = np.flatnonzero(lasso_coef > 1e-4)
if lasso_selected.size == 0:
    lasso_selected = np.argsort(-lasso_coef)[: max(2, int(len(features) * 0.2))]

svm_rfe = RFE(
    estimator=SVC(kernel="linear", random_state=SEED),
    n_features_to_select=n_select,
    step=1,
)
svm_rfe.fit(x, y)

rf = RandomForestClassifier(
    n_estimators=500,
    random_state=SEED,
    class_weight="balanced",
)
rf.fit(x, y)
rf_importance = rf.feature_importances_
rf_threshold = np.sort(rf_importance)[::-1][n_select - 1]

lasso_set = {features[i] for i in lasso_selected}
svm_set = {features[i] for i in np.flatnonzero(svm_rfe.support_)}
rf_set = {features[i] for i in np.flatnonzero(rf_importance >= rf_threshold)}

stats = all_results.set_index("Gene")
module_lookup = modules.set_index("Gene")["Module"]
rows = []
for i, gene in enumerate(features):
    result = stats.loc[gene]
    votes = {
        "LASSO_selected": int(gene in lasso_set),
        "SVM_RFE_selected": int(gene in svm_set),
        "RF_selected": int(gene in rf_set),
    }
    rows.append(
        {
            "Gene": gene,
            "WGCNA_Module": str(module_lookup.loc[gene]),
            **votes,
            "Total_Votes_out_of_3": sum(votes.values()),
            "Majority_Selected_2of3": int(sum(votes.values()) >= 2),
            "LASSO_abs_coefficient": float(lasso_coef[i]),
            "SVM_RFE_rank": int(svm_rfe.ranking_[i]),
            "RF_importance": float(rf_importance[i]),
            "Discovery_log2FC": float(result["logFC"]),
            "Discovery_nominal_P": float(result["P.Value"]),
            "Discovery_BH_adjusted_P": float(result["adj.P.Val"]),
            "Discovery_DEG_Status": str(result["DEG_Status"]),
        }
    )

results = pd.DataFrame(rows).sort_values(
    ["Total_Votes_out_of_3", "LASSO_abs_coefficient", "RF_importance"],
    ascending=[False, False, False],
)
result_path = REVISION_DIR / "ml_624_intersection_results.csv"
results.to_csv(result_path, index=False)

selected = results.loc[results["Majority_Selected_2of3"] == 1]
summary = {
    "dataset": "GSE125498",
    "samples": int(len(sample_ids)),
    "case_samples": int(y.sum()),
    "control_samples": int(len(y) - y.sum()),
    "input_feature_count": int(len(features)),
    "feature_pool_definition": (
        "Intersection of 813 primary DEGs with non-grey genes assigned by dynamic WGCNA"
    ),
    "primary_deg_criteria": "|log2FC| >= 0.585 and nominal P < 0.05",
    "random_seed": SEED,
    "selection_method": {
        "LASSO": "LogisticRegressionCV; L1; liblinear; 5-fold CV; Cs=[0.01,0.05,0.1,0.5,1,5,10]; selected |coef| > 1e-4",
        "SVM_RFE": f"Linear SVC recursive elimination; top {n_select} features; step=1",
        "Random_Forest": f"500 trees; balanced class weights; top {n_select} importance threshold (ties included)",
    },
    "selected_counts": {
        "LASSO": int(len(lasso_set)),
        "SVM_RFE": int(len(svm_set)),
        "Random_Forest": int(len(rf_set)),
        "majority_consensus_at_least_2_of_3": int(len(selected)),
        "strict_consensus_3_of_3": int((results["Total_Votes_out_of_3"] == 3).sum()),
    },
    "majority_consensus_genes": selected["Gene"].tolist(),
    "performance_note": (
        "Feature selection used the full discovery cohort. These outputs do not estimate "
        "out-of-sample performance; that requires selection nested inside cross-validation."
    ),
    "output_csv": str(result_path.relative_to(ROOT)),
}
summary_path = REVISION_DIR / "ml_624_intersection_summary.json"
summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
