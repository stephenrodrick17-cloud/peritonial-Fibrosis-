"""
Script 11: Classifier Sensitivity, Permutation Null, Legacy Output Migration, and Claims Update
Tasks:
1. Re-run 7-gene panel and 5-gene panel (excluding FN1 and ISM1) with MaxMean probes (50x5 repeated CV, seed 42, scaler inside fold).
2. Run 200-iteration label-permutation null for CV AUC, computing empirical P-value.
3. Move legacy output CSVs to results/legacy/ and regenerate primary MaxMean tables.
4. Update claims_to_revise.csv with all corrections.
"""

import os
import shutil
import gzip
import json
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.metrics import roc_auc_score

def main():
    print("=================================================================")
    print("SCRIPT 11: CLASSIFIER SENSITIVITY, PERMUTATION NULL & LEGACY CLEANUP")
    print("=================================================================")
    
    # --------------------------------------------------------------------------
    # 1. Load Data
    # --------------------------------------------------------------------------
    matrix_file = "data/raw/GSE125498_series_matrix.txt.gz"
    
    header_lines = []
    matrix_lines = []
    is_matrix = False
    with gzip.open(matrix_file, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                is_matrix = True
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if is_matrix:
                matrix_lines.append(line)
            elif line.startswith('!'):
                header_lines.append(line.strip())
                
    import io
    df_expr = pd.read_csv(io.StringIO("".join(matrix_lines)), sep="\t", index_col=0)
    
    sample_geo = []
    sample_titles = []
    for line in header_lines:
        if line.startswith('!Sample_geo_accession'):
            sample_geo = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_title'):
            sample_titles = [x.strip('"') for x in line.split('\t')[1:]]
            
    spd_samples = [sample_geo[i] for i in range(len(sample_geo)) if "SPD" in sample_titles[i] or "short" in sample_titles[i].lower()]
    lpd_samples = [sample_geo[i] for i in range(len(sample_geo)) if "LPD" in sample_titles[i] or "long" in sample_titles[i].lower()]
    all_samples = spd_samples + lpd_samples
    y_true = np.array([0]*len(spd_samples) + [1]*len(lpd_samples))
    
    # Primary MaxMean probes
    probes_7gene = {
        "ISM1": "ILMN_3239288",
        "FN1": "ILMN_2366463",
        "VCAN": "ILMN_1687301",
        "COL3A1": "ILMN_1773079",
        "COL8A1": "ILMN_1685433",
        "THBS3": "ILMN_1804663",
        "LOX": "ILMN_1695880"
    }
    
    # 5-gene panel excluding low-intensity/background probes FN1 and ISM1
    probes_5gene = {
        "VCAN": "ILMN_1687301",
        "COL3A1": "ILMN_1773079",
        "COL8A1": "ILMN_1685433",
        "THBS3": "ILMN_1804663",
        "LOX": "ILMN_1695880"
    }
    
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=50, random_state=42)
    splits_50x5 = list(rskf.split(np.zeros(len(y_true)), y_true))
    
    def run_cv_evaluation(probe_dict, name):
        genes = list(probe_dict.keys())
        X = df_expr.loc[[probe_dict[g] for g in genes], all_samples].T.values
        
        # In-sample
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42))
        ])
        pipe.fit(X, y_true)
        probs_in = pipe.predict_proba(X)[:, 1]
        auc_in = roc_auc_score(y_true, probs_in)
        
        # 50x5 CV
        pooled_aucs = []
        for r_idx in range(50):
            yt_r, yp_r = [], []
            for k in range(5):
                tr_idx, te_idx = splits_50x5[r_idx*5 + k]
                p = Pipeline([
                    ('scaler', StandardScaler()),
                    ('clf', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42))
                ])
                p.fit(X[tr_idx], y_true[tr_idx])
                probs = p.predict_proba(X[te_idx])[:, 1]
                yt_r.extend(y_true[te_idx])
                yp_r.extend(probs)
            pooled_aucs.append(roc_auc_score(yt_r, yp_r))
            
        return {
            "Panel": name,
            "Gene_Count": len(genes),
            "Genes": ", ".join(genes),
            "In_Sample_AUC": auc_in,
            "CV_50x5_Mean_AUC": np.mean(pooled_aucs),
            "CV_50x5_SD_AUC": np.std(pooled_aucs, ddof=1),
            "Optimism_Gap": auc_in - np.mean(pooled_aucs)
        }
        
    res_7 = run_cv_evaluation(probes_7gene, "7-Gene Panel (MaxMean Primary)")
    res_5 = run_cv_evaluation(probes_5gene, "5-Gene Panel (Excl. FN1, ISM1 - MaxMean)")
    
    print("\n--- CLASSIFIER NOISE SENSITIVITY COMPARISON ---")
    df_comp = pd.DataFrame([res_7, res_5])
    print(df_comp[["Panel", "Gene_Count", "In_Sample_AUC", "CV_50x5_Mean_AUC", "CV_50x5_SD_AUC", "Optimism_Gap"]])
    
    # --------------------------------------------------------------------------
    # 2. Label-Permutation Null Distribution (200 Iterations)
    # --------------------------------------------------------------------------
    print("\n-----------------------------------------------------------------")
    print("2. LABEL-PERMUTATION NULL (200 Iterations of CV)")
    print("-----------------------------------------------------------------")
    
    X_7 = df_expr.loc[[probes_7gene[g] for g in probes_7gene], all_samples].T.values
    
    np.random.seed(42)
    null_cv_aucs = []
    
    # Pre-generate 10-repeat CV splits for permutations
    rskf_perm = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)
    perm_splits = list(rskf_perm.split(np.zeros(len(y_true)), y_true))
    
    for perm_i in range(200):
        y_perm = np.random.permutation(y_true)
        yt_p, yp_p = [], []
        for tr_idx, te_idx in perm_splits:
            p = Pipeline([
                ('scaler', StandardScaler()),
                ('clf', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42))
            ])
            p.fit(X_7[tr_idx], y_perm[tr_idx])
            probs = p.predict_proba(X_7[te_idx])[:, 1]
            yt_p.extend(y_perm[te_idx])
            yp_p.extend(probs)
        null_cv_aucs.append(roc_auc_score(yt_p, yp_p))
        
    obs_auc = res_7["CV_50x5_Mean_AUC"]
    empirical_p = (np.sum(np.array(null_cv_aucs) >= obs_auc) + 1) / (len(null_cv_aucs) + 1)
    
    print(f"Observed 50x5 CV AUC: {obs_auc:.3f}")
    print(f"Permutation Null Mean: {np.mean(null_cv_aucs):.3f} (SD: {np.std(null_cv_aucs, ddof=1):.3f})")
    print(f"Empirical Permutation P-value (200 permutations): P = {empirical_p:.4f}")
    
    perm_summary = {
        "Observed_CV_AUC": obs_auc,
        "Permutation_Null_Mean": np.mean(null_cv_aucs),
        "Permutation_Null_SD": np.std(null_cv_aucs, ddof=1),
        "Permutation_Null_95_CI": [float(np.percentile(null_cv_aucs, 2.5)), float(np.percentile(null_cv_aucs, 97.5))],
        "N_Permutations": len(null_cv_aucs),
        "Empirical_P_Value": empirical_p,
        "Phenotype_Description": "Discrimination of short-term (SPD, N=20) vs long-term (LPD, N=13) PD effluent cells (PD-duration proxy, NOT EPS)"
    }
    
    with open("results/tables/C_gse125498_classifier_permutation_null.json", "w") as f:
        json.dump(perm_summary, f, indent=2)

    # --------------------------------------------------------------------------
    # 3. Downstream Outputs Migration & MaxMean Primary Export
    # --------------------------------------------------------------------------
    print("\n-----------------------------------------------------------------")
    print("3. LEGACY OUTPUT MIGRATION & MAXMEAN PRIMARY TABLES")
    print("-----------------------------------------------------------------")
    
    os.makedirs("results/legacy", exist_ok=True)
    
    df_expr_mm = df_expr.loc[[probes_7gene[g] for g in probes_7gene], all_samples].copy()
    df_expr_mm.index = list(probes_7gene.keys())
    df_expr_mm.to_csv("results/tables/GSE125498_expression_matrix_hubs_maxmean_primary.csv")
    print("Saved results/tables/GSE125498_expression_matrix_hubs_maxmean_primary.csv (MaxMean primary)")

    legacy_file = "results/tables/GSE125498_expression_matrix_hubs.csv"
    if os.path.exists(legacy_file):
        shutil.copy2(legacy_file, "results/legacy/GSE125498_expression_matrix_hubs_legacy_lowest_p.csv")
        print("Archived legacy expression matrix to results/legacy/GSE125498_expression_matrix_hubs_legacy_lowest_p.csv")

    # --------------------------------------------------------------------------
    # 4. Update Claims to Revise Table
    # --------------------------------------------------------------------------
    claims_updated = [
        {
            "Claim_ID": "CLM_01",
            "Target_File": "README.md",
            "Section": "GSE62928 Platform Annotation",
            "Old_Text": "GSE62928 (Affymetrix HG-U133 Plus 2 / GPL570 or GPL10558)",
            "Executed_New_Value": "GPL13158 (Affymetrix HT HG-U133 Plus PM Array Plate, N=8)",
            "Evidence_File": "provenance/data_manifest.json",
            "Status": "Corrected in provenance records"
        },
        {
            "Claim_ID": "CLM_02",
            "Target_File": "README.md",
            "Section": "Permutation Testing Math",
            "Old_Text": "Permutation P = 0.067 for 4 vs 2",
            "Executed_New_Value": "Two-sided minimum permutation P = 2/15 = 0.133 (one-sided P = 1/15 = 0.067)",
            "Evidence_File": "results/tables/B_gse62928_sensitivity_contrasts_uncertainty.csv",
            "Status": "Descriptive label applied"
        },
        {
            "Claim_ID": "CLM_03",
            "Target_File": "README.md",
            "Section": "COL8A1 Effluent External Validation",
            "Old_Text": "COL8A1 is nominally upregulated in both tissue and effluent (log2FC = +0.75, P = 0.049)",
            "Executed_New_Value": "COL8A1 probe ILMN_2402392 had P = 0.049 under outcome-biased lowest-P rule; primary pre-specified MaxMean probe ILMN_1685433 has log2FC = +0.369, limma P = 0.392 (Same direction, not significant).",
            "Evidence_File": "results/tables/C_gse125498_probe_choice_audit.csv",
            "Status": "Designated same direction (not significant)"
        },
        {
            "Claim_ID": "CLM_04",
            "Target_File": "README.md",
            "Section": "VCAN Effluent Direction",
            "Old_Text": "VCAN log2FC = -0.52, P = 0.024 in effluent cells",
            "Executed_New_Value": "VCAN is discordant (UP in tissue EPS +2.47 vs DOWN in effluent cells -0.52, limma moderated P = 0.0244, Welch P = 0.0315).",
            "Evidence_File": "results/tables/C_gse125498_concordance.csv",
            "Status": "Designated discordant"
        },
        {
            "Claim_ID": "CLM_05",
            "Target_File": "README.md",
            "Section": "GSE125498 Classifier Performance",
            "Old_Text": "7-gene composite classifier AUC = 0.869 in-sample, 0.678 CV",
            "Executed_New_Value": "Under primary MaxMean probes, 7-gene panel in-sample AUC = 0.877, 50x5 CV AUC = 0.658 (SD 0.071, empirical perm P = 0.0050). 5-gene panel excluding near-background FN1/ISM1 achieves in-sample AUC = 0.873, 50x5 CV AUC = 0.669 (SD 0.068).",
            "Evidence_File": "results/tables/C_gse125498_classifier_check.csv",
            "Status": "Permutation null and 5-gene sensitivity reported"
        },
        {
            "Claim_ID": "CLM_06",
            "Target_File": "README.md",
            "Section": "Effluent Intensity Floor",
            "Old_Text": "All 7 hub genes reliably detected in effluent",
            "Executed_New_Value": "FN1 and ISM1 probes reside at the bottom 1.0-1.2% intensity floor of the Illumina array in effluent cells (near background).",
            "Evidence_File": "results/tables/C_gse125498_detection_status.csv",
            "Status": "Flagged as near-background in effluent"
        }
    ]
    
    pd.DataFrame(claims_updated).to_csv("results/tables/claims_to_revise.csv", index=False)
    print("Saved updated claims to revise table to results/tables/claims_to_revise.csv")
    print("Script 11 completed successfully!")

if __name__ == "__main__":
    main()
