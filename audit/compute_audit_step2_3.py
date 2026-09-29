"""
Comprehensive Statistical & Technical Audit Computation Script for Step 2 and Step 3
"""

import os
import itertools
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu, binomtest
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.metrics import roc_curve, auc, accuracy_score, recall_score

print("=" * 80)
print("COMPUTING ALL STEP 2 & STEP 3 STATISTICAL AUDIT REQUIREMENTS")
print("=" * 80)

# ==============================================================================
# STEP 3.1: DIRECTION CONCORDANCE (GSE62928 vs GSE125498)
# ==============================================================================
print("\n--- STEP 3.1: DIRECTION CONCORDANCE TABLE ---")
df_hubs = pd.read_csv("results/tables/ML_hub_genes_from_WGCNA_ECM.csv")
df_val = pd.read_csv("results/tables/GSE125498_wgcna_hub_validation_metrics.csv")

# Merge on Gene_Symbol
df_concordance = df_val.merge(df_hubs[["Gene_Symbol", "Peritoneal_logFC", "Peritoneal_Pval", "WGCNA_MM", "WGCNA_GS"]], on="Gene_Symbol")

df_concordance["GSE62928_Direction"] = np.where(df_concordance["Peritoneal_logFC"] > 0, "UP (Pro-fibrotic)", "DOWN")
df_concordance["GSE125498_Direction"] = np.where(df_concordance["log2_FC_GSE125498"] > 0, "UP (Late > Early)", "DOWN (Late < Early)")
df_concordance["Concordance_Status"] = np.where(
    np.sign(df_concordance["Peritoneal_logFC"]) == np.sign(df_concordance["log2_FC_GSE125498"]),
    "CONCORDANT", "DISCORDANT (Opposite Sign)"
)

cols_conc = [
    "Gene_Symbol", "Probe_ID", "Peritoneal_logFC", "GSE62928_Direction",
    "log2_FC_GSE125498", "GSE125498_Direction", "Concordance_Status",
    "Limma_P_Value", "Mann_Whitney_Pval", "ROC_AUC"
]
print(df_concordance[cols_conc].to_string(index=False))


# ==============================================================================
# STEP 3.2: PERMUTATION TEST WITH FAMILY-WISE MAX |r| ACROSS 14 MODULES
# ==============================================================================
print("\n--- STEP 3.2: PERMUTATION TESTS (SINGLE-HYPOTHESIS & FAMILY-WISE) ---")
df_mes = pd.read_csv("results/tables/wgcna_module_eigengenes.csv", index_col=0)
df_meta = pd.read_csv("results/tables/GSE62928_sample_metadata.csv")
df_mes = df_mes.loc[df_meta["sample_id"]]
me_salmon = df_mes["MEsalmon"].values
true_trait = df_meta["binary_numeric"].values.astype(int)

obs_r_salmon = pearsonr(me_salmon, true_trait)[0]
mod_cols = [c for c in df_mes.columns if c.startswith("ME")]

# Exact 70 permutations
indices = list(range(8))
all_combinations = list(itertools.combinations(indices, 4))
n_exact = len(all_combinations)

salmon_r_list = []
max_abs_r_all_14_list = []

for case_idx in all_combinations:
    perm_trait = np.zeros(8, dtype=int)
    perm_trait[list(case_idx)] = 1
    
    # Salmon r
    r_sal, _ = pearsonr(me_salmon, perm_trait)
    salmon_r_list.append(r_sal)
    
    # Max abs r across all 14 modules
    all_14_cors = [abs(pearsonr(df_mes[m].values, perm_trait)[0]) for m in mod_cols]
    max_abs_r_all_14_list.append(max(all_14_cors))

salmon_r_arr = np.array(salmon_r_list)
max_abs_r_arr = np.array(max_abs_r_all_14_list)

p_exact_one_sided = np.mean(salmon_r_arr >= obs_r_salmon - 1e-9)
p_exact_two_sided = np.mean(np.abs(salmon_r_arr) >= abs(obs_r_salmon) - 1e-9)
p_family_wise = np.mean(max_abs_r_arr >= abs(obs_r_salmon) - 1e-9)

print(f"Observed Salmon r: {obs_r_salmon:.5f}")
print(f"1. Single-Hypothesis Exact One-Sided P (r >= {obs_r_salmon:.3f}): {p_exact_one_sided:.4f} ({np.sum(salmon_r_arr >= obs_r_salmon - 1e-9)}/{n_exact})")
print(f"2. Single-Hypothesis Exact Two-Sided P (|r| >= {abs(obs_r_salmon):.3f}): {p_exact_two_sided:.4f} ({np.sum(np.abs(salmon_r_arr) >= abs(obs_r_salmon) - 1e-9)}/{n_exact})")
print(f"3. Family-Wise Max-|r| Permutation P (Max |r| across 14 modules >= {abs(obs_r_salmon):.3f}): {p_family_wise:.4f} ({np.sum(max_abs_r_arr >= abs(obs_r_salmon) - 1e-9)}/{n_exact})")


# ==============================================================================
# STEP 3.3: NEGATIVE CONTROL (Check 6 False Positive Rate)
# ==============================================================================
print("\n--- STEP 3.3: NEGATIVE CONTROL (RANDOM NOISE TRAITS) ---")
non_true_max_r = []
non_true_min_p = []

for case_idx in all_combinations:
    trait_i = np.zeros(8, dtype=int)
    trait_i[list(case_idx)] = 1
    if np.array_equal(trait_i, true_trait) or np.array_equal(trait_i, 1 - true_trait):
        continue
    
    cors_i = [abs(pearsonr(df_mes[m].values, trait_i)[0]) for m in mod_cols]
    pvals_i = [pearsonr(df_mes[m].values, trait_i)[1] for m in mod_cols]
    non_true_max_r.append(max(cors_i))
    non_true_min_p.append(min(pvals_i))

fpr_p05 = np.mean(np.array(non_true_min_p) < 0.05) * 100
fpr_r70 = np.mean(np.array(non_true_max_r) >= 0.70) * 100
print(f"Across {len(non_true_max_r)} non-true random 4v4 trait permutations:")
print(f"  * False Positive Rate (At least 1 module P < 0.05): {fpr_p05:.1f}%")
print(f"  * False Positive Rate (At least 1 module |r| >= 0.70): {fpr_r70:.1f}%")
print(f"  * Average Maximum |r| from pure noise: {np.mean(non_true_max_r):.4f}")
print(f"  * Maximum |r| from pure noise: {np.max(non_true_max_r):.4f} (P = {np.min(non_true_min_p):.4e})")


# ==============================================================================
# STEP 3.4: NESTED VALIDATION (LOOCV WITH FULL FEATURE SELECTION INSIDE EACH FOLD)
# ==============================================================================
print("\n--- STEP 3.4: NESTED LOOCV VALIDATION (FEATURE SELECTION INSIDE EACH FOLD) ---")
# Feature selection inside each fold from the 40 convergent candidate genes
df_candidates = pd.read_csv("results/tables/convergent_WGCNA_ECM_genes.csv")
candidate_genes = df_candidates["gene_symbol"].tolist()
df_expr62 = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
available_genes = [g for g in candidate_genes if g in df_expr62.index]
sample_ids = df_meta["sample_id"].tolist()
y = df_meta["binary_numeric"].values.astype(int)
X = df_expr62.loc[available_genes, sample_ids].T.values

loo = LeaveOneOut()
nested_preds_lasso = []
nested_preds_svm = []
nested_preds_rf = []
nested_preds_xgb = []
nested_preds_consensus = []

for fold, (train_idx, test_idx) in enumerate(loo.split(X)):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    # Internal Standardizing
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    
    # 1. LASSO inside fold
    clf_l = LogisticRegression(penalty="l1", C=0.8, solver="liblinear", random_state=42).fit(X_train_s, y_train)
    p_l = clf_l.predict(X_test_s)[0]
    nested_preds_lasso.append(p_l)
    
    # 2. SVM-RFE inside fold
    rfe = RFE(estimator=SVC(kernel="linear", C=1.0, random_state=42), n_features_to_select=max(5, int(len(available_genes)*0.35)), step=1).fit(X_train_s, y_train)
    p_s = rfe.predict(X_test_s)[0]
    nested_preds_svm.append(p_s)
    
    # 3. Random Forest inside fold
    rf = RandomForestClassifier(n_estimators=500, max_depth=3, random_state=42).fit(X_train_s, y_train)
    p_r = rf.predict(X_test_s)[0]
    nested_preds_rf.append(p_r)
    
    # 4. XGBoost inside fold
    xgb = XGBClassifier(n_estimators=100, max_depth=2, learning_rate=0.08, random_state=42, eval_metric="logloss").fit(X_train_s, y_train)
    p_x = xgb.predict(X_test_s)[0]
    nested_preds_xgb.append(p_x)
    
    # Consensus vote
    p_c = 1 if (p_l + p_s + p_r + p_x) >= 2 else 0
    nested_preds_consensus.append(p_c)

acc_cons = accuracy_score(y, nested_preds_consensus)
k_correct = int(np.sum(nested_preds_consensus == y))
n_total_loo = len(y)

# Exact Clopper-Pearson / Binomial Confidence Interval
b_test = binomtest(k_correct, n_total_loo, p=0.5)
ci_low, ci_high = b_test.proportion_ci(confidence_level=0.95, method="exact")

print(f"Nested LOOCV Consensus Accuracy on GSE62928 (N=8): {k_correct}/{n_total_loo} = {acc_cons*100:.1f}%")
print(f"Exact Binomial 95% Confidence Interval: [{ci_low*100:.1f}%, {ci_high*100:.1f}%]")
print(f"Binomial Test vs Null (50% chance): P = {b_test.pvalue:.4f}")
print(f"Per-Model Nested LOOCV Accuracies: LASSO = {accuracy_score(y, nested_preds_lasso)*100:.1f}%, SVM-RFE = {accuracy_score(y, nested_preds_svm)*100:.1f}%, RF = {accuracy_score(y, nested_preds_rf)*100:.1f}%, XGBoost = {accuracy_score(y, nested_preds_xgb)*100:.1f}%")


# ==============================================================================
# STEP 3.6: THRESHOLD ROBUSTNESS ACROSS VOTE CUTOFFS
# ==============================================================================
print("\n--- STEP 3.6: THRESHOLD ROBUSTNESS (VOTES >= 2, >= 3, 4/4) ---")
df_ml_all = pd.read_csv("results/tables/ML_hub_genes_from_WGCNA_ECM_all_results.csv")
v_4 = df_ml_all[df_ml_all["Votes"] == 4]["Gene_Symbol"].tolist()
v_3 = df_ml_all[df_ml_all["Votes"] >= 3]["Gene_Symbol"].tolist()
v_2 = df_ml_all[df_ml_all["Votes"] >= 2]["Gene_Symbol"].tolist()

print(f"Hub genes with 4/4 Votes ({len(v_4)} genes): {v_4 if v_4 else 'None'}")
print(f"Hub genes with >=3 Votes ({len(v_3)} genes): {v_3}")
print(f"Hub genes with >=2 Votes ({len(v_2)} genes): {v_2}")


# ==============================================================================
# STEP 3.7: COMPOSITE CLASSIFIER LOGISTIC REGRESSION COEFFICIENTS & CV +/- DISCORDANT GENES
# ==============================================================================
print("\n--- STEP 3.7: COMPOSITE CLASSIFIER IN GSE125498 (+/- VCAN & THBS3) ---")
import GEOparse
gse_val = GEOparse.get_GEO(filepath="data/GSE125498_family.soft.gz")
piv_val = gse_val.pivot_samples("VALUE")
df_tt_val = pd.read_csv("Validation/GSE125498.top.table.tsv", sep="\t")
df_tt_val["Gene_clean"] = df_tt_val["Gene.symbol"].fillna("").astype(str).str.strip().str.upper()

sample_stages = {}
for name, gsm in gse_val.gsms.items():
    title = gsm.metadata.get("title", [""])[0]
    chars = str(gsm.metadata.get("characteristics_ch1", []))
    is_long = "long-term" in title.lower() or "long-term" in chars.lower()
    sample_stages[name] = "Late_Stage_LPD" if is_long else "Early_Stage_SPD"

early_samples = [s for s, st in sample_stages.items() if st == "Early_Stage_SPD"]
late_samples = [s for s, st in sample_stages.items() if st == "Late_Stage_LPD"]
all_samples_val = early_samples + late_samples
y_val = np.array([1 if sample_stages[s] == "Late_Stage_LPD" else 0 for s in all_samples_val])

# Extract probe data for 7 available genes
gene_exprs_val = {}
for g in ["ISM1", "FN1", "VCAN", "COL3A1", "COL8A1", "THBS3", "LOX"]:
    tt_m = df_tt_val[df_tt_val["Gene_clean"].apply(lambda x: g == x or g in [s.strip() for s in x.split("///")])]
    best_p = tt_m.sort_values("P.Value").iloc[0]["ID"]
    gene_exprs_val[g] = piv_val.loc[best_p, all_samples_val].astype(float).values

X_all_7 = pd.DataFrame(gene_exprs_val, index=all_samples_val)

# Fit Logistic Regression on all 7 genes to get coefficients
clf_7 = LogisticRegression(random_state=42)
clf_7.fit(X_all_7, y_val)

print("Logistic Regression Model Coefficients (All 7 Profiled Genes in GSE125498):")
for g, coef in zip(X_all_7.columns, clf_7.coef_[0]):
    print(f"  * {g:<8}: Coef = {coef:+.4f} (Direction in Model: {'PRO-LATE' if coef > 0 else 'PRO-EARLY'})")

# Repeated Stratified CV (10 x 5-fold = 50 repeats) for:
# Model A: Full 7-Gene Panel
# Model B: Concordant 5-Gene Panel (Excluding VCAN and THBS3)
def eval_repeated_cv(X_mat, y_vec, n_repeats=10, n_splits=5):
    aucs = []
    for r in range(n_repeats):
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42 + r)
        oof_probs = np.zeros(len(y_vec))
        for tr, te in skf.split(X_mat, y_vec):
            clf = LogisticRegression(random_state=42)
            clf.fit(X_mat.iloc[tr], y_vec[tr])
            oof_probs[te] = clf.predict_proba(X_mat.iloc[te])[:, 1]
        fpr, tpr, _ = roc_curve(y_vec, oof_probs)
        aucs.append(auc(fpr, tpr))
    return np.mean(aucs), np.std(aucs)

mean_7, std_7 = eval_repeated_cv(X_all_7, y_val, n_repeats=10, n_splits=5)

# Concordant 5 (dropping VCAN and THBS3)
X_5_concordant = X_all_7[["ISM1", "FN1", "COL3A1", "COL8A1", "LOX"]]
mean_5, std_5 = eval_repeated_cv(X_5_concordant, y_val, n_repeats=10, n_splits=5)

# Single gene VCAN alone
mean_vcan, std_vcan = eval_repeated_cv(X_all_7[["VCAN"]], y_val, n_repeats=10, n_splits=5)

# Single gene COL8A1 alone
mean_col8, std_col8 = eval_repeated_cv(X_all_7[["COL8A1"]], y_val, n_repeats=10, n_splits=5)

print("\nRepeated Stratified 5-Fold Cross-Validation Performance (10 Repeats x 5 Folds):")
print(f"  * Full 7-Gene Panel (All Available):            Mean AUC = {mean_7:.3f} +/- {std_7:.3f}")
print(f"  * Concordant 5-Gene Panel (Excl. VCAN & THBS3): Mean AUC = {mean_5:.3f} +/- {std_5:.3f}")
print(f"  * VCAN Alone (Single Feature):                 Mean AUC = {mean_vcan:.3f} +/- {std_vcan:.3f}")
print(f"  * COL8A1 Alone (Single Feature):               Mean AUC = {mean_col8:.3f} +/- {std_col8:.3f}")
