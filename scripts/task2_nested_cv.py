import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

os.makedirs("results/revision", exist_ok=True)
os.makedirs("scripts", exist_ok=True)

# Set random seed
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# 1. Load data
meta = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
expr_df = pd.read_csv("results/tables/GSE125498_full_expression_matrix.csv").set_index("Gene")
sample_cols = [c for c in meta["sample_id"] if c in expr_df.columns]
y = meta.set_index("sample_id").loc[sample_cols]["stage_binary"].values.astype(int)
genes = np.array(expr_df.index)
X_all = expr_df[sample_cols].values # shape: (11741, 33)

ecm_ref = pd.read_excel("ECM genes all.xlsx", sheet_name="Hs_ECM_Masterlist", header=1)
ecm_set = set(ecm_ref["Gene Symbol"].dropna().astype(str).str.strip())
is_ecm = np.array([g in ecm_set for g in genes])

print(f"Dataset: N={len(sample_cols)} samples (Cases LPD={int(y.sum())}, Controls SPD={int((1-y).sum())})")
print(f"Total Genes: {len(genes)}, Matrisome Genes Measured: {int(is_ecm.sum())}")

# 2. Setup Repeated Stratified 5-Fold CV (10 repeats, seed 42)
N_SPLITS = 5
N_REPEATS = 10
rskf = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS, random_state=RANDOM_STATE)

repeat_aucs = []
repeat_c_indices = []
fold_records = []

total_folds = N_SPLITS * N_REPEATS
current_repeat = 0
oof_preds = np.zeros(len(y))
fold_count = 0

print(f"Starting Nested Cross-Validation ({N_REPEATS} repeats of {N_SPLITS}-fold CV = {total_folds} total folds)...")

for fold_idx, (train_idx, test_idx) in enumerate(rskf.split(X_all.T, y)):
    rep_idx = fold_idx // N_SPLITS
    fold_in_rep = fold_idx % N_SPLITS
    
    # Train / Test split
    X_train = X_all[:, train_idx]
    y_train = y[train_idx]
    X_test = X_all[:, test_idx]
    y_test = y[test_idx]
    
    n_cases = int(y_train.sum())
    n_ctrls = int((1 - y_train).sum())
    
    # --- STEP A: INSIDE-FOLD DIFFERENTIAL EXPRESSION ---
    m1 = X_train[:, y_train == 1].mean(axis=1)
    m0 = X_train[:, y_train == 0].mean(axis=1)
    logfc = m1 - m0
    
    v1 = X_train[:, y_train == 1].var(axis=1, ddof=1)
    v0 = X_train[:, y_train == 0].var(axis=1, ddof=1)
    se = np.sqrt(v1 / n_cases + v0 / n_ctrls)
    t_stat = (m1 - m0) / np.clip(se, 1e-9, None)
    
    # Two-sided p-value
    df = n_cases + n_ctrls - 2
    p_vals = 2.0 * stats.t.sf(np.abs(t_stat), df=df)
    
    # Up-DEGs filter: log2FC >= 0.585 and nominal P < 0.05
    is_up_deg = (logfc >= 0.585) & (p_vals < 0.05)
    
    # Intersect with Matrisome
    is_candidate = is_up_deg & is_ecm
    cand_indices = np.where(is_candidate)[0]
    cand_genes = genes[cand_indices].tolist()
    
    if len(cand_genes) < 2:
        raise RuntimeError(
            f"Repeat {rep_idx + 1}, fold {fold_in_rep + 1} produced fewer than two "
            "upregulated Matrisome DEGs under the prespecified thresholds."
        )
    
    # Extract candidate feature matrix for train and test
    X_cand_tr = X_train[cand_indices, :].T
    X_cand_te = X_test[cand_indices, :].T
    
    scaler = StandardScaler()
    X_sc_tr = scaler.fit_transform(X_cand_tr)
    X_sc_te = scaler.transform(X_cand_te)
    
    n_cand = len(cand_genes)
    
    # --- STEP B: INSIDE-FOLD 3-MODEL ML SELECTION (LASSO, SVM-RFE, RF) ---
    # 1. LASSO
    lasso = LogisticRegressionCV(Cs=[0.01, 0.1, 1.0, 10.0], cv=3, penalty='l1',
                                 solver='liblinear', max_iter=2000, random_state=RANDOM_STATE)
    lasso.fit(X_sc_tr, y_train)
    lasso_coefs = np.abs(lasso.coef_[0])
    lasso_sel = set([cand_genes[i] for i in np.where(lasso_coefs > 1e-4)[0]])
    
    # 2. SVM-RFE
    n_rfe = max(2, min(7, int(np.ceil(n_cand * 0.35))))
    svm = RFE(estimator=SVC(kernel='linear', random_state=RANDOM_STATE), n_features_to_select=n_rfe, step=1)
    svm.fit(X_sc_tr, y_train)
    svm_sel = set([cand_genes[i] for i in range(n_cand) if svm.support_[i]])
        
    # 3. Random Forest
    rf = RandomForestClassifier(n_estimators=100, random_state=RANDOM_STATE, class_weight='balanced')
    rf.fit(X_sc_tr, y_train)
    rf_imp = rf.feature_importances_
    rf_thresh = np.sort(rf_imp)[::-1][min(n_rfe - 1, n_cand - 1)]
    rf_sel = set([cand_genes[i] for i in range(n_cand) if rf_imp[i] >= rf_thresh])
    
    # Consensus voting (>= 2 of 3 votes)
    votes = {g: (int(g in lasso_sel) + int(g in svm_sel) + int(g in rf_sel)) for g in cand_genes}
    selected_hubs = [g for g, v in votes.items() if v >= 2]
    if len(selected_hubs) < 2:
        sorted_by_votes = sorted(votes.items(), key=lambda x: x[1], reverse=True)
        selected_hubs = [x[0] for x in sorted_by_votes[:min(5, len(sorted_by_votes))]]
        
    # --- STEP C: INSIDE-FOLD FINAL MODEL FITTING & PREDICTION ---
    hub_indices = [cand_genes.index(g) for g in selected_hubs]
    X_hub_tr = X_sc_tr[:, hub_indices]
    X_hub_te = X_sc_te[:, hub_indices]
    
    final_clf = LogisticRegression(penalty='l2', C=1.0, max_iter=1000, random_state=RANDOM_STATE)
    final_clf.fit(X_hub_tr, y_train)
    pred_probs = final_clf.predict_proba(X_hub_te)[:, 1]
    
    oof_preds[test_idx] = pred_probs
    
    fold_records.append({
        "Repeat": rep_idx + 1,
        "Fold": fold_in_rep + 1,
        "N_Train": len(train_idx),
        "N_Test": len(test_idx),
        "N_Up_DEGs_Fold": int(is_up_deg.sum()),
        "N_Up_ECM_DEGs_Fold": n_cand,
        "N_Selected_Hubs_Fold": len(selected_hubs),
        "Selected_Hubs": ";".join(selected_hubs)
    })
    
    if fold_in_rep == N_SPLITS - 1:
        # Completed one full repeat
        rep_auc = roc_auc_score(y, oof_preds)
        # C-index for binary logistic regression equals AUC
        rep_c_index = rep_auc
        repeat_aucs.append(rep_auc)
        repeat_c_indices.append(rep_c_index)
        print(f"  -> Repeat {rep_idx + 1}/{N_REPEATS}: Out-of-Fold AUC = {rep_auc:.4f}")
        oof_preds = np.zeros(len(y))

# 3. Compute Summary Statistics across 10 Repeats
mean_auc = float(np.mean(repeat_aucs))
std_auc = float(np.std(repeat_aucs, ddof=1))
mean_c = float(np.mean(repeat_c_indices))
std_c = float(np.std(repeat_c_indices, ddof=1))

bootstrap_rng = np.random.default_rng(RANDOM_STATE)
bootstrap_auc_means = np.array([
    bootstrap_rng.choice(repeat_aucs, size=N_REPEATS, replace=True).mean()
    for _ in range(10_000)
])
ci_lower_auc, ci_upper_auc = np.percentile(bootstrap_auc_means, [2.5, 97.5])
ci_lower_c, ci_upper_c = ci_lower_auc, ci_upper_auc

print("\n" + "=" * 70)
print("NESTED CROSS-VALIDATION SUMMARY (10 Repeats x 5-Fold CV):")
print("=" * 70)
print(
    f"Mean Nested AUC:     {mean_auc:.4f} (SD across repeats: {std_auc:.4f}; "
    f"bootstrap 95% CI: [{ci_lower_auc:.4f}, {ci_upper_auc:.4f}])"
)
print(
    f"Mean Nested C-index: {mean_c:.4f} (SD across repeats: {std_c:.4f}; "
    f"bootstrap 95% CI: [{ci_lower_c:.4f}, {ci_upper_c:.4f}])"
)

# 4. Save Results
legacy_results_path = "results/revision/nested_cv_results.csv"
if not os.path.isfile(legacy_results_path):
    raise FileNotFoundError(
        f"Historical comparison table is required to retain non-nested values: {legacy_results_path}"
    )
legacy_rows = pd.read_csv(legacy_results_path).iloc[:2].to_dict(orient="records")
comparison_rows = legacy_rows + [
    {
        "Evaluation_Paradigm": "Nested CV (Repeated 5-Fold, 10 Repeats)",
        "Leakage_Status": "DEG filter + Matrisome + ML selection inside folds; WGCNA not recomputed",
        "Metric": "ROC AUC",
        "Estimate": round(mean_auc, 4),
        "CI_95_Lower": round(ci_lower_auc, 4),
        "CI_95_Upper": round(ci_upper_auc, 4),
        "SD": round(std_auc, 4),
        "Notes": "Mean of repeat-level OOF AUCs; percentile bootstrap resamples the 10 repeats"
    },
    {
        "Evaluation_Paradigm": "Nested CV (Repeated 5-Fold, 10 Repeats)",
        "Leakage_Status": "DEG filter + Matrisome + ML selection inside folds; WGCNA not recomputed",
        "Metric": "C-index",
        "Estimate": round(mean_c, 4),
        "CI_95_Lower": round(ci_lower_c, 4),
        "CI_95_Upper": round(ci_upper_c, 4),
        "SD": round(std_c, 4),
        "Notes": "Binary-model C-index equals AUC; percentile bootstrap resamples the 10 repeats"
    }
]

res_df = pd.DataFrame(comparison_rows)
res_df.to_csv("results/revision/nested_cv_results_taskE_v2.csv", index=False)

# Also save per-fold detailed log
folds_df = pd.DataFrame(fold_records)
folds_df.to_csv("results/revision/nested_cv_fold_log_taskE_v2.csv", index=False)

summary_json = {
    "cv_protocol": "Repeated Stratified 5-Fold CV (10 repeats, 50 total splits)",
    "seed": RANDOM_STATE,
    "fold_feature_selection": "Unmoderated two-sample t filter with Welch variance estimate, Naba Matrisome intersection, LASSO, SVM-RFE, RF",
    "wgcna_recomputed_within_folds": False,
    "ci_method": "Percentile bootstrap over 10 repeat-level out-of-fold AUC estimates; repeats are the bootstrap units",
    "bootstrap_replicates": 10000,
    "bootstrap_seed": RANDOM_STATE,
    "non_nested_apparent_auc": float(legacy_rows[0]["Estimate"]),
    "non_nested_cv_c_index": float(legacy_rows[1]["Estimate"]),
    "nested_cv_auc_mean": round(mean_auc, 4),
    "nested_cv_auc_sd": round(std_auc, 4),
    "nested_cv_auc_ci95": [round(ci_lower_auc, 4), round(ci_upper_auc, 4)],
    "nested_cv_c_index_mean": round(mean_c, 4),
    "nested_cv_c_index_sd": round(std_c, 4),
    "nested_cv_c_index_ci95": [round(ci_lower_c, 4), round(ci_upper_c, 4)],
    "repeat_aucs": [round(x, 4) for x in repeat_aucs],
    "repeat_c_indices": [round(x, 4) for x in repeat_c_indices]
}
with open("results/revision/nested_cv_summary_taskE_v2.json", "w") as f:
    json.dump(summary_json, f, indent=2)

print("Saved Task E nested CV outputs under results/revision with the taskE_v2 suffix.")
