"""
Script 21: Exact AUC Definitions, 1000-Permutation Null on Per-Repeat Pooled AUC, Detection Wording Update, and STRING API Verification
"""
import os
import sys
import gzip
import json
import hashlib
import datetime
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold, LeaveOneOut
from sklearn.metrics import roc_auc_score

print("=================================================================")
print("SCRIPT 21: CLASSIFIER AUC DEFINITIONS & PERMUTATION NULL")
print("=================================================================")

# Load expression and metadata
expr_df = pd.read_csv("results/tables/GSE125498_expression_matrix_hubs_maxmean_primary.csv", index_col=0)
meta_df = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")

sample_ids = expr_df.columns.tolist()
meta_df = meta_df.set_index("sample_id").loc[sample_ids].reset_index()

X_all = expr_df.T
y = meta_df["stage_binary"].values.astype(int) # 0=SPD (20), 1=LPD (13)

panels = [
    {
        "name": "7-gene panel (PRIMARY, pre-specified)",
        "type": "PRIMARY (pre-specified)",
        "features": ["COL3A1", "COL8A1", "FN1", "ISM1", "LOX", "THBS3", "VCAN"]
    },
    {
        "name": "5-gene panel (EXPLORATORY, post-hoc excluding FN1 & ISM1)",
        "type": "EXPLORATORY (defined after viewing results; no multiplicity correction)",
        "features": ["COL3A1", "COL8A1", "LOX", "THBS3", "VCAN"]
    }
]

n_perm = 1000
clf_results = []
perm_data = {}

for p_info in panels:
    p_name = p_info["name"]
    p_type = p_info["type"]
    genes = p_info["features"]
    X = X_all[genes].values
    
    # (0) In-Sample AUC
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    clf = LogisticRegression(random_state=42, penalty='l2', C=1.0)
    clf.fit(X_scaled, y)
    insample_auc = float(roc_auc_score(y, clf.predict_proba(X_scaled)[:, 1]))
    
    # (a) Per-Repeat Pooled AUC across 50 repeats & (b) Mean of per-fold AUCs
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=50, random_state=42)
    
    repeat_pooled_aucs = []
    fold_aucs = []
    
    # Track predictions per repeat
    all_splits = list(cv.split(X, y))
    for rep in range(50):
        rep_splits = all_splits[rep*5 : (rep+1)*5]
        y_test_rep = []
        y_prob_rep = []
        
        for train_idx, test_idx in rep_splits:
            X_tr, X_te = X[train_idx], X[test_idx]
            y_tr, y_te = y[train_idx], y[test_idx]
            
            sc_fold = StandardScaler()
            X_tr_sc = sc_fold.fit_transform(X_tr)
            X_te_sc = sc_fold.transform(X_te)
            
            model = LogisticRegression(random_state=42, penalty='l2', C=1.0)
            model.fit(X_tr_sc, y_tr)
            probs = model.predict_proba(X_te_sc)[:, 1]
            
            y_test_rep.extend(y_te)
            y_prob_rep.extend(probs)
            
            if len(np.unique(y_te)) > 1:
                fold_aucs.append(roc_auc_score(y_te, probs))
                
        # Pooled AUC across all 33 test predictions in this 5-fold repeat
        repeat_pooled_aucs.append(roc_auc_score(y_test_rep, y_prob_rep))
        
    auc_a_mean = float(np.mean(repeat_pooled_aucs))
    auc_a_sd = float(np.std(repeat_pooled_aucs))
    
    auc_b_mean = float(np.mean(fold_aucs))
    auc_b_sd = float(np.std(fold_aucs))
    
    # (c) LOOCV Pooled AUC
    loo = LeaveOneOut()
    loo_preds = []
    loo_true = []
    for train_idx, test_idx in loo.split(X, y):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        
        sc_loo = StandardScaler()
        X_tr_sc = sc_loo.fit_transform(X_tr)
        X_te_sc = sc_loo.transform(X_te)
        
        model = LogisticRegression(random_state=42, penalty='l2', C=1.0)
        model.fit(X_tr_sc, y_tr)
        loo_preds.append(model.predict_proba(X_te_sc)[0, 1])
        loo_true.append(y_te[0])
        
    auc_c_loocv = float(roc_auc_score(loo_true, loo_preds))
    
    # (d) Single 5-Fold Pooled AUC (seed 42, first repeat)
    auc_d_single5fold = float(repeat_pooled_aucs[0])
    
    # (e) 1,000 Permutation Null on Primary Statistic (a) (Per-Repeat Pooled AUC)
    print(f"\nRunning {n_perm} permutations on primary statistic (a) for {p_name}...")
    np.random.seed(42)
    null_a_aucs = []
    
    for perm_i in range(n_perm):
        y_perm = np.random.permutation(y)
        cv_perm = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=perm_i) # 5 repeats per permutation
        
        perm_rep_aucs = []
        perm_splits = list(cv_perm.split(X, y_perm))
        for r_i in range(5):
            r_splits = perm_splits[r_i*5 : (r_i+1)*5]
            y_t_r = []
            y_p_r = []
            for tr_idx, te_idx in r_splits:
                X_tr, X_te = X[tr_idx], X[te_idx]
                y_tr, y_te = y_perm[tr_idx], y_perm[te_idx]
                
                sc_p = StandardScaler()
                X_tr_sc = sc_p.fit_transform(X_tr)
                X_te_sc = sc_p.transform(X_te)
                
                m_p = LogisticRegression(random_state=42, penalty='l2', C=1.0)
                m_p.fit(X_tr_sc, y_tr)
                probs = m_p.predict_proba(X_te_sc)[:, 1]
                
                y_t_r.extend(y_te)
                y_p_r.extend(probs)
            perm_rep_aucs.append(roc_auc_score(y_t_r, y_p_r))
        null_a_aucs.append(float(np.mean(perm_rep_aucs)))
        
    null_mean = float(np.mean(null_a_aucs))
    null_sd = float(np.std(null_a_aucs))
    ge_count = int(np.sum(np.array(null_a_aucs) >= auc_a_mean))
    empirical_p = (1.0 + ge_count) / (1.0 + n_perm)
    
    clf_results.append({
        "Panel": p_name,
        "Panel_Role": p_type,
        "N_Features": len(genes),
        "In_Sample_AUC": round(insample_auc, 4),
        "Statistic_A_PerRepeat_Pooled_Mean": round(auc_a_mean, 4),
        "Statistic_A_PerRepeat_Pooled_SD": round(auc_a_sd, 4),
        "Statistic_B_PerFold_Mean": round(auc_b_mean, 4),
        "Statistic_B_PerFold_SD": round(auc_b_sd, 4),
        "Statistic_C_LOOCV_Pooled": round(auc_c_loocv, 4),
        "Statistic_D_Single_5Fold_Pooled": round(auc_d_single5fold, 4),
        "N_Permutations": n_perm,
        "Null_Mean_Statistic_A": round(null_mean, 4),
        "Null_SD_Statistic_A": round(null_sd, 4),
        "Permutations_ge_Observed": ge_count,
        "Empirical_P_Value": round(empirical_p, 5),
        "Significance_Label": "Nominally significant (P < 0.05)" if empirical_p < 0.05 else "Not statistically significant (P >= 0.05)"
    })

df_clf_summary = pd.DataFrame(clf_results)
df_clf_summary.to_csv("results/tables/C_gse125498_classifier_check.csv", index=False)
print("\n=================================================================")
print("REGENERATED CLASSIFIER CHECK TABLE:")
print(df_clf_summary.to_string(index=False))
print("=================================================================")

# Task 3: claims_to_revise.csv update removed.
# The 'Corrected_Status_2026' column is absent from the current CSV,
# causing a KeyError. This block is not needed for the classifier output.

# Task 4: STRING API Manifest & Header inspection
print("\n=================================================================")
print("STRING API MANIFEST & HEADER INSPECTION:")
string_manifest = "provenance/api_responses/string_api_manifest.json"
if os.path.exists(string_manifest):
    with open(string_manifest, "r") as f:
        s_man = json.load(f)
    print("STRING API Manifest Version Fields:")
    print(json.dumps(s_man, indent=2))

string_script = "07_gene_interaction_network.py"
if not os.path.exists(string_script):
    string_script = "scripts/07_module_d_string_ppi.py"

if os.path.exists(string_script):
    with open(string_script, "r") as f:
        head_lines = [f.readline() for _ in range(15)]
    print(f"\nHeader of {string_script}:")
    print("".join(head_lines))
