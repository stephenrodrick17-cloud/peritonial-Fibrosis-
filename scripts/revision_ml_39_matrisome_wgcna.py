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
OUTPUT_DIR = REVISION_DIR / "ml_39_matrisome_wgcna"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

expression = pd.read_csv(ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv")
metadata = pd.read_csv(ROOT / "results" / "tables" / "GSE125498_sample_metadata.csv")
modules = pd.read_csv(REVISION_DIR / "wgcna_modules_dynamic_taskB.csv")
matrisome = pd.read_csv(
    REVISION_DIR / "matrisome_volcano" / "GSE125498_matrisome_volcano_source_data.csv"
)
all_results = pd.read_csv(ROOT / "results" / "tables" / "GSE125498_all_results.csv")

required_matrisome_columns = {"Gene", "DEG_Status"}
if not required_matrisome_columns.issubset(matrisome.columns):
    raise ValueError("MatrisomeAnalyzeR source table is missing Gene or DEG_Status.")
if matrisome["Gene"].duplicated().any():
    raise ValueError("MatrisomeAnalyzeR source table contains duplicate gene symbols.")

significant_matrisome = matrisome.loc[
    matrisome["DEG_Status"].isin(["Upregulated", "Downregulated"])
].copy()
non_grey_modules = modules.loc[
    modules["Module"].astype(str).str.lower() != "grey"
].copy()
module_genes = set(non_grey_modules["Gene"].dropna().astype(str))
matrisome_genes = set(significant_matrisome["Gene"].dropna().astype(str))
features = sorted(module_genes & matrisome_genes)

if len(significant_matrisome) != 46 or len(features) != 39:
    raise ValueError(
        f"Expected 46 significant MatrisomeAnalyzeR genes and 39 non-grey "
        f"WGCNA overlaps; found {len(significant_matrisome)} and {len(features)}."
    )

input_genes = significant_matrisome.set_index("Gene").loc[features].reset_index()
direction_counts = input_genes["DEG_Status"].value_counts().to_dict()
if direction_counts != {"Upregulated": 28, "Downregulated": 11}:
    raise ValueError(f"Expected 28 upregulated and 11 downregulated genes; found {direction_counts}.")

sample_ids = metadata["sample_id"].astype(str).tolist()
if not set(sample_ids).issubset(expression.columns):
    raise ValueError("Metadata contains sample IDs missing from the expression matrix.")
expression_genes = set(expression["Gene"].astype(str))
if not set(features).issubset(expression_genes):
    raise ValueError("Some intersection genes are missing from the expression matrix.")

expr = expression.set_index("Gene")
x_raw = expr.loc[features, sample_ids].to_numpy(dtype=float).T
meta_by_sample = metadata.set_index("sample_id").loc[sample_ids]
y = meta_by_sample["stage_binary"].to_numpy(dtype=int)
if np.unique(y).size != 2:
    raise ValueError("The discovery cohort must contain both phenotype groups.")

# This full-cohort feature-selection run is exploratory, not an estimate of
# out-of-sample performance; that requires feature selection nested in CV.
x = StandardScaler().fit_transform(x_raw)
n_select = max(2, min(10, int(np.ceil(len(features) * 0.35))))

lasso = LogisticRegressionCV(
    Cs=[0.01, 0.05, 0.1, 0.5, 1, 5, 10],
    cv=5,
    scoring="accuracy",
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
deg_lookup = significant_matrisome.set_index("Gene")["DEG_Status"]
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
            "Matrisome_DEG_Status": str(deg_lookup.loc[gene]),
            **votes,
            "Total_Votes_out_of_3": sum(votes.values()),
            "Consensus_At_Least_2_of_3": int(sum(votes.values()) >= 2),
            "LASSO_abs_coefficient": float(lasso_coef[i]),
            "SVM_RFE_rank": int(svm_rfe.ranking_[i]),
            "RF_importance": float(rf_importance[i]),
            "Discovery_log2FC": float(result["logFC"]),
            "Discovery_nominal_P": float(result["P.Value"]),
            "Discovery_BH_adjusted_P": float(result["adj.P.Val"]),
        }
    )

results = pd.DataFrame(rows).sort_values(
    ["Total_Votes_out_of_3", "LASSO_abs_coefficient", "RF_importance"],
    ascending=[False, False, False],
)
selected = results.loc[results["Consensus_At_Least_2_of_3"] == 1].copy()

input_path = OUTPUT_DIR / "wgcna_matrisome_39_gene_input.csv"
results_path = OUTPUT_DIR / "ml_39_gene_results.csv"
summary_path = OUTPUT_DIR / "ml_39_gene_summary.json"
input_genes.to_csv(input_path, index=False)
results.to_csv(results_path, index=False)

summary = {
    "dataset": "GSE125498",
    "samples": int(len(sample_ids)),
    "long_term_samples": int(y.sum()),
    "short_term_samples": int(len(y) - y.sum()),
    "matrisome_deg_count": int(len(significant_matrisome)),
    "input_feature_count": int(len(features)),
    "input_direction_counts": {
        "upregulated": int(direction_counts["Upregulated"]),
        "downregulated": int(direction_counts["Downregulated"]),
    },
    "feature_pool_definition": (
        "Intersection of the 46 MatrisomeAnalyzeR DEGs "
        "(|log2FC| >= 0.585 and nominal P < 0.05) with genes assigned "
        "to non-grey modules in the current dynamic WGCNA result"
    ),
    "random_seed": SEED,
    "selection_method": {
        "LASSO": (
            "LogisticRegressionCV; L1; liblinear; 5-fold CV; "
            "accuracy scoring; Cs=[0.01,0.05,0.1,0.5,1,5,10]; "
            "selected |coefficient| > 1e-4; "
            "if none pass, retain the top max(2, floor(20% of input)) coefficients"
        ),
        "SVM_RFE": (
            f"Linear SVC recursive elimination; top {n_select} features; step=1; "
            "no separate RFE cross-validation"
        ),
        "Random_Forest": (
            f"500 trees; balanced class weights; top {n_select} importance "
            "threshold (ties included)"
        ),
    },
    "selected_counts": {
        "LASSO": int(len(lasso_set)),
        "SVM_RFE": int(len(svm_set)),
        "Random_Forest": int(len(rf_set)),
        "union_at_least_1_of_3": int((results["Total_Votes_out_of_3"] >= 1).sum()),
        "consensus_at_least_2_of_3": int(len(selected)),
        "strict_consensus_3_of_3": int((results["Total_Votes_out_of_3"] == 3).sum()),
    },
    "consensus_genes": selected["Gene"].tolist(),
    "performance_note": (
        "Feature selection used the full discovery cohort. These results do not "
        "estimate out-of-sample performance and do not constitute independent validation."
    ),
    "input_csv": str(input_path.relative_to(ROOT)),
    "results_csv": str(results_path.relative_to(ROOT)),
}
summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
