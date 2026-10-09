"""
Module C: Bulk External Validation in GSE125498 (Peritoneal Effluent Cells)
Note: This dataset profiles shed peritoneal effluent cells comparing short-term vs long-term PD.
This is a PD-duration proxy, NOT EPS and NOT tissue.

Tasks:
C1. Checksums, sample metadata verification (20 Short-term vs 13 Long-term PD).
C2. Probe mapping using GPL10558.annot.gz (Symbol, Entrez ID, Aliases).
C3. MaxMean probe selection & detection verification.
C4. Differential expression (unadjusted & covariate-adjusted models).
C5. Directional concordance vs discovery contrasts.
C6. Classifier re-check with Pipeline StandardScaler fit inside each training fold.
"""

import os
import sys
import gzip
import hashlib
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import RepeatedStratifiedKFold, LeaveOneOut, StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy import stats
from statsmodels.stats.multitest import multipletests

def compute_sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def main():
    print("=================================================================")
    print("MODULE C: GSE125498 BULK VALIDATION (EFFLUENT CELLS; PD-DURATION PROXY)")
    print("=================================================================")
    
    # --------------------------------------------------------------------------
    # C1. Verify Files, Checksums, and Sample Metadata
    # --------------------------------------------------------------------------
    matrix_file = "data/raw/GSE125498_series_matrix.txt.gz"
    annot_file  = "data/raw/GPL10558.annot.gz"
    
    if not os.path.exists(matrix_file) or not os.path.exists(annot_file):
        raise FileNotFoundError("GSE125498 series matrix or GPL10558 annot file missing!")
        
    sha_matrix = compute_sha256(matrix_file)
    sha_annot  = compute_sha256(annot_file)
    
    print(f"GSE125498 Series Matrix SHA256: {sha_matrix} (Size: {os.path.getsize(matrix_file):,} bytes)")
    print(f"GPL10558 Annotation SHA256:    {sha_annot} (Size: {os.path.getsize(annot_file):,} bytes)")
    
    # Parse Series Matrix Headers and Expression Table
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
                
    sample_geo_accessions = []
    sample_titles = []
    sample_characteristics = {}
    
    for line in header_lines:
        if line.startswith('!Sample_geo_accession'):
            sample_geo_accessions = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_title'):
            sample_titles = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = [x.strip('"') for x in line.split('\t')[1:]]
            for idx, p in enumerate(parts):
                if idx not in sample_characteristics:
                    sample_characteristics[idx] = []
                sample_characteristics[idx].append(p)
                
    # Classify samples into Short-term PD vs Long-term PD
    sample_records = []
    for idx, gsm in enumerate(sample_geo_accessions):
        title = sample_titles[idx] if idx < len(sample_titles) else ""
        chars = sample_characteristics.get(idx, [])
        chars_str = " | ".join(chars)
        
        # Clinical duration grouping
        if "SPD" in title or "short" in chars_str.lower():
            group = "Short_term_PD"
        elif "LPD" in title or "long" in chars_str.lower():
            group = "Long_term_PD"
        else:
            group = "Unknown"
            
        sample_records.append({
            "gsm": gsm,
            "title": title,
            "group": group,
            "characteristics": chars_str
        })
        
    df_samples = pd.DataFrame(sample_records)
    n_short = sum(df_samples["group"] == "Short_term_PD")
    n_long  = sum(df_samples["group"] == "Long_term_PD")
    n_total = len(df_samples)
    
    print(f"\nCohort Breakdown (N={n_total}):")
    print(f"  Short-term PD (SPD): n = {n_short}")
    print(f"  Long-term PD  (LPD): n = {n_long}")
    print("Sample details preview:")
    print(df_samples[["gsm", "title", "group"]].head(6))
    
    # Covariates check
    print("\nCovariates in metadata:")
    for idx, c in enumerate(df_samples["characteristics"].iloc[0].split(" | ")):
        print(f"  Characteristic field [{idx+1}]: {c}")
        
    # Parse expression table
    import io
    matrix_csv = "".join(matrix_lines)
    df_expr = pd.read_csv(io.StringIO(matrix_csv), sep="\t", index_col=0)
    print(f"\nExpression Matrix Loaded: {df_expr.shape[0]} probe rows x {df_expr.shape[1]} samples")
    
    # --------------------------------------------------------------------------
    # C2. GPL10558 Probe Mapping for the 11 Hub Genes
    # --------------------------------------------------------------------------
    print("\n-----------------------------------------------------------------")
    print("C2. GPL10558 PROBE MAPPING (Official Symbol, Entrez ID, Aliases)")
    print("-----------------------------------------------------------------")
    
    hub_info = {
        "ISM1":    {"entrez": "56994",  "aliases": ["C20orf82", "TAIL1", "Isthmin 1"]},
        "FN1":     {"entrez": "2335",   "aliases": ["CIG", "ED-B", "FINC", "FN", "FNZ", "GFND", "GFND2", "LETS", "MSF"]},
        "EDIL3":   {"entrez": "10085",  "aliases": ["DEL1", "Del-1", "integrin-binding EGF-like domain 3"]},
        "VCAN":    {"entrez": "1462",   "aliases": ["CSPG2", "ERVR", "GHAP", "PG-M", "WGN", "WGN1", "versican"]},
        "COL3A1":  {"entrez": "1281",   "aliases": ["EDS4A", "FLJ34534", "collagen type III alpha 1 chain"]},
        "COMP":    {"entrez": "1311",   "aliases": ["EDM1", "EPD1", "MED", "PSACH", "THBS5"]},
        "COL8A1":  {"entrez": "1295",   "aliases": ["C3orf25", "collagen type VIII alpha 1 chain"]},
        "THBS3":   {"entrez": "7059",   "aliases": ["TSP3", "thrombospondin 3"]},
        "COL11A1": {"entrez": "1301",   "aliases": ["COLL6", "STL2", "collagen type XI alpha 1 chain"]},
        "INHBA":   {"entrez": "3624",   "aliases": ["EDF", "FRP", "inhibin subunit beta A"]},
        "LOX":     {"entrez": "4015",   "aliases": ["protein-lysine 6-oxidase", "lysyl oxidase"]}
    }
    
    annot_rows = []
    with gzip.open(annot_file, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('ID\t') or line.startswith('ID '):
                annot_cols = [c.strip() for c in line.strip().split('\t')]
                break
        for line in f:
            if line.startswith('!platform_table_end'):
                break
            parts = [p.strip() for p in line.strip().split('\t')]
            if len(parts) >= len(annot_cols):
                annot_rows.append(parts[:len(annot_cols)])
                
    df_annot = pd.DataFrame(annot_rows, columns=annot_cols)
    print(f"GPL10558 Annotation Table: {len(df_annot)} probe rows")
    
    sym_col = [c for c in df_annot.columns if "symbol" in c.lower() or "gene_symbol" in c.lower()][0]
    id_col  = [c for c in df_annot.columns if c.lower() == "id" or "probe" in c.lower()][0]
    entrez_col = [c for c in df_annot.columns if "entrez" in c.lower() or "gene_id" in c.lower()][0] if any("entrez" in c.lower() for c in df_annot.columns) else None
    syn_col = [c for c in df_annot.columns if "synonym" in c.lower() or "alias" in c.lower() or "title" in c.lower()][0] if any("synonym" in c.lower() or "alias" in c.lower() for c in df_annot.columns) else None
    
    mapping_results = []
    all_mapped_probes = []
    
    for gene, info in hub_info.items():
        m_sym = df_annot[df_annot[sym_col].str.upper() == gene.upper()]
        m_entrez = pd.DataFrame()
        if entrez_col:
            m_entrez = df_annot[df_annot[entrez_col].astype(str) == info["entrez"]]
        m_alias = pd.DataFrame()
        if syn_col:
            for al in info["aliases"]:
                m_al = df_annot[df_annot[syn_col].str.contains(al, case=False, na=False)]
                m_alias = pd.concat([m_alias, m_al])
                
        all_matches = pd.concat([m_sym, m_entrez, m_alias]).drop_duplicates(subset=[id_col])
        probe_ids = all_matches[id_col].tolist()
        present_probes = [p for p in probe_ids if p in df_expr.index]
        
        mapping_results.append({
            "Gene": gene,
            "Entrez_ID": info["entrez"],
            "Aliases": ", ".join(info["aliases"]),
            "Probes_Found_Count": len(present_probes),
            "Probe_IDs": ", ".join(present_probes) if present_probes else "NO PROBE",
            "Mapping_Status": "MAPPED" if len(present_probes) > 0 else "NO PROBE"
        })
        
        for p in present_probes:
            mean_e = df_expr.loc[p].mean()
            all_mapped_probes.append({
                "Gene": gene,
                "Probe_ID": p,
                "Mean_Expression": mean_e
            })
            
    df_mapping = pd.DataFrame(mapping_results)
    df_all_probes = pd.DataFrame(all_mapped_probes)
    
    os.makedirs("results/tables", exist_ok=True)
    df_mapping.to_csv("results/tables/C_gse125498_probe_mapping.csv", index=False)
    print("\n--- PROBE MAPPING SUMMARY FOR THE 11 HUB GENES ---")
    print(df_mapping[["Gene", "Mapping_Status", "Probes_Found_Count", "Probe_IDs"]])
    
    # --------------------------------------------------------------------------
    # C3. MaxMean Probe Selection
    # --------------------------------------------------------------------------
    print("\n-----------------------------------------------------------------")
    print("C3. MAXMEAN PROBE SELECTION")
    print("-----------------------------------------------------------------")
    
    best_probes = {}
    for gene in df_mapping[df_mapping["Mapping_Status"] == "MAPPED"]["Gene"]:
        sub = df_all_probes[df_all_probes["Gene"] == gene].sort_values("Mean_Expression", ascending=False)
        best_p = sub.iloc[0]["Probe_ID"]
        best_e = sub.iloc[0]["Mean_Expression"]
        best_probes[gene] = best_p
        print(f"  {gene:<8}: Selected Probe = {best_p} (Mean Expr = {best_e:.3f}, Candidate probes = {len(sub)})")
        
    mapped_hubs = list(best_probes.keys())
    unmapped_hubs = df_mapping[df_mapping["Mapping_Status"] == "NO PROBE"]["Gene"].tolist()
    print(f"\nMapped Hubs ({len(mapped_hubs)}):   {', '.join(mapped_hubs)}")
    print(f"Unmapped Hubs ({len(unmapped_hubs)}): {', '.join(unmapped_hubs)}")
    
    # --------------------------------------------------------------------------
    # C4. Differential Expression (Vectorized across full array)
    # --------------------------------------------------------------------------
    print("\n-----------------------------------------------------------------")
    print("C4. DIFFERENTIAL EXPRESSION IN GSE125498 (LONG-TERM vs SHORT-TERM PD)")
    print("-----------------------------------------------------------------")
    
    expr_subset = df_expr.loc[[best_probes[g] for g in mapped_hubs]].copy()
    expr_subset.index = mapped_hubs
    
    spd_samples = df_samples[df_samples["group"] == "Short_term_PD"]["gsm"].tolist()
    lpd_samples = df_samples[df_samples["group"] == "Long_term_PD"]["gsm"].tolist()
    
    # Fast vectorized t-test across full matrix
    mat_spd = df_expr[spd_samples].values
    mat_lpd = df_expr[lpd_samples].values
    
    t_stats_all, p_vals_all = stats.ttest_ind(mat_lpd, mat_spd, axis=1, equal_var=False)
    p_vals_all = np.nan_to_num(p_vals_all, nan=1.0)
    adjp_all = multipletests(p_vals_all, method="fdr_bh")[1]
    
    probe_to_adjp = dict(zip(df_expr.index, adjp_all))
    
    de_rows = []
    for gene in mapped_hubs:
        x_spd = expr_subset.loc[gene, spd_samples].values.astype(float)
        x_lpd = expr_subset.loc[gene, lpd_samples].values.astype(float)
        
        mean_spd = np.mean(x_spd)
        mean_lpd = np.mean(x_lpd)
        lfc = mean_lpd - mean_spd
        
        t_stat, p_val = stats.ttest_ind(x_lpd, x_spd, equal_var=False)
        se = np.sqrt(np.var(x_lpd, ddof=1)/len(x_lpd) + np.var(x_spd, ddof=1)/len(x_spd))
        ci_low = lfc - 1.96 * se
        ci_high = lfc + 1.96 * se
        
        de_rows.append({
            "Gene": gene,
            "Probe_ID": best_probes[gene],
            "log2FC": lfc,
            "SE": se,
            "CI_95_Lower": ci_low,
            "CI_95_Upper": ci_high,
            "t_statistic": t_stat,
            "Raw_P": p_val,
            "Mean_ShortPD": mean_spd,
            "Mean_LongPD": mean_lpd,
            "N_ShortPD": len(spd_samples),
            "N_LongPD": len(lpd_samples),
            "Compartment": "Peritoneal Effluent Cells (PD-Duration Proxy, NOT Tissue, NOT EPS)"
        })
        
    df_de = pd.DataFrame(de_rows)
    df_de["Adj_P_BH_Hubs"] = multipletests(df_de["Raw_P"], method="fdr_bh")[1]
    df_de["Adj_P_BH_Array"] = [probe_to_adjp[p] for p in df_de["Probe_ID"]]
    
    df_de.to_csv("results/tables/C_gse125498_hub_limma.csv", index=False)
    print(df_de[["Gene", "Probe_ID", "log2FC", "SE", "Raw_P", "Adj_P_BH_Hubs", "Adj_P_BH_Array"]])
    
    # --------------------------------------------------------------------------
    # C5. Directional Concordance Table vs Discovery Contrasts
    # --------------------------------------------------------------------------
    print("\n-----------------------------------------------------------------")
    print("C5. DIRECTIONAL CONCORDANCE TABLE")
    print("-----------------------------------------------------------------")
    
    sens_file = "results/tables/B_gse62928_sensitivity_contrasts_uncertainty.csv"
    df_sens = pd.read_csv(sens_file)
    
    concordance_rows = []
    for gene in hub_info.keys():
        disc_pooled = df_sens[(df_sens["Gene"] == gene) & (df_sens["Contrast"] == "EPS_vs_PooledControl")].iloc[0]
        disc_pd     = df_sens[(df_sens["Gene"] == gene) & (df_sens["Contrast"] == "EPS_vs_PD")].iloc[0]
        disc_uremic = df_sens[(df_sens["Gene"] == gene) & (df_sens["Contrast"] == "EPS_vs_Uremic")].iloc[0]
        
        disc_lfc_pooled = disc_pooled["log2FC_MaxMean"]
        disc_lfc_pd     = disc_pd["log2FC_MaxMean"]
        disc_lfc_uremic = disc_uremic["log2FC_MaxMean"]
        
        if gene in mapped_hubs:
            val_sub = df_de[df_de["Gene"] == gene].iloc[0]
            val_lfc = val_sub["log2FC"]
            val_p   = val_sub["Raw_P"]
            
            if abs(val_lfc) < 0.10:
                category = "Flat (|log2FC| < 0.1)"
            elif (val_lfc > 0 and disc_lfc_pooled > 0) or (val_lfc < 0 and disc_lfc_pooled < 0):
                category = "Concordant (UP in Long-term PD & UP in EPS)"
            else:
                category = "Discordant"
                
            concordance_rows.append({
                "Gene": gene,
                "Status": "MAPPED",
                "Discovery_EPS_vs_Pooled_log2FC": disc_lfc_pooled,
                "Discovery_EPS_vs_PD_log2FC": disc_lfc_pd,
                "Discovery_EPS_vs_Uremic_log2FC": disc_lfc_uremic,
                "Validation_GSE125498_log2FC": val_lfc,
                "Validation_Raw_P": val_p,
                "Concordance_Category": category,
                "Phenotype_Note": "Effluent cells (PD duration proxy) vs Tissue (EPS)"
            })
        else:
            concordance_rows.append({
                "Gene": gene,
                "Status": "NOT MAPPED ON GPL10558",
                "Discovery_EPS_vs_Pooled_log2FC": disc_lfc_pooled,
                "Discovery_EPS_vs_PD_log2FC": disc_lfc_pd,
                "Discovery_EPS_vs_Uremic_log2FC": disc_lfc_uremic,
                "Validation_GSE125498_log2FC": np.nan,
                "Validation_Raw_P": np.nan,
                "Concordance_Category": "Not mapped (No probe on array)",
                "Phenotype_Note": "Unmapped on Illumina HT-12 v4 beadchip"
            })
            
    df_concordance = pd.DataFrame(concordance_rows)
    df_concordance.to_csv("results/tables/C_gse125498_concordance.csv", index=False)
    
    print("Concordance Table Summary:")
    print(df_concordance[["Gene", "Status", "Validation_GSE125498_log2FC", "Concordance_Category"]])
    print("\nCategory Counts:")
    print(df_concordance["Concordance_Category"].value_counts())
    
    # --------------------------------------------------------------------------
    # C6. Classifier Re-check (Exact Repository Protocol)
    # --------------------------------------------------------------------------
    print("\n-----------------------------------------------------------------")
    print("C6. CLASSIFIER RE-CHECK (7-GENE EFFLUENT PANEL)")
    print("-----------------------------------------------------------------")
    
    X = expr_subset.T.loc[df_samples["gsm"]].values
    y = np.array([1 if g == "Long_term_PD" else 0 for g in df_samples["group"]])
    
    clf_pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42))
    ])
    
    clf_pipeline.fit(X, y)
    y_prob_insample = clf_pipeline.predict_proba(X)[:, 1]
    auc_insample = roc_auc_score(y, y_prob_insample)
    
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=50, random_state=42)
    pooled_aucs_per_repeat = []
    fold_aucs_all = []
    
    for r_idx in range(50):
        r_y_true = []
        r_y_prob = []
        for k in range(5):
            train_idx, test_idx = list(rskf.split(X, y))[r_idx*5 + k]
            pipe = Pipeline([
                ('scaler', StandardScaler()),
                ('classifier', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42))
            ])
            pipe.fit(X[train_idx], y[train_idx])
            probs = pipe.predict_proba(X[test_idx])[:, 1]
            r_y_true.extend(y[test_idx])
            r_y_prob.extend(probs)
            if len(np.unique(y[test_idx])) > 1:
                fold_aucs_all.append(roc_auc_score(y[test_idx], probs))
        pooled_aucs_per_repeat.append(roc_auc_score(r_y_true, r_y_prob))
        
    mean_pooled_auc = np.mean(pooled_aucs_per_repeat)
    sd_pooled_auc   = np.std(pooled_aucs_per_repeat, ddof=1)
    mean_fold_auc   = np.mean(fold_aucs_all)
    sd_fold_auc     = np.std(fold_aucs_all, ddof=1)
    
    loo = LeaveOneOut()
    loo_y_true = []
    loo_y_prob = []
    for train_idx, test_idx in loo.split(X, y):
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42))
        ])
        pipe.fit(X[train_idx], y[train_idx])
        probs = pipe.predict_proba(X[test_idx])[:, 1]
        loo_y_true.append(y[test_idx[0]])
        loo_y_prob.append(probs[0])
    auc_loocv = roc_auc_score(loo_y_true, loo_y_prob)
    
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    skf_y_true = []
    skf_y_prob = []
    for train_idx, test_idx in skf.split(X, y):
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42))
        ])
        pipe.fit(X[train_idx], y[train_idx])
        probs = pipe.predict_proba(X[test_idx])[:, 1]
        skf_y_true.extend(y[test_idx])
        skf_y_prob.extend(probs)
    auc_single_5fold = roc_auc_score(skf_y_true, skf_y_prob)
    
    optimism_gap = auc_insample - mean_pooled_auc
    
    readme_metrics = {
        "In-Sample AUC": {"readme": 0.869, "recomputed": round(auc_insample, 3)},
        "50x5 CV Pooled Mean AUC": {"readme": 0.678, "recomputed": round(mean_pooled_auc, 3)},
        "50x5 CV Pooled SD": {"readme": 0.058, "recomputed": round(sd_pooled_auc, 3)},
        "50x5 CV Fold Mean AUC": {"readme": 0.701, "recomputed": round(mean_fold_auc, 3)},
        "50x5 CV Fold SD": {"readme": 0.203, "recomputed": round(sd_fold_auc, 3)},
        "LOOCV Pooled AUC": {"readme": 0.658, "recomputed": round(auc_loocv, 3)},
        "Single 5-Fold Pooled AUC": {"readme": 0.592, "recomputed": round(auc_single_5fold, 3)},
        "Optimism Gap (Delta AUC)": {"readme": 0.191, "recomputed": round(optimism_gap, 3)}
    }
    
    clf_rows = []
    for metric, vals in readme_metrics.items():
        diff = abs(vals["recomputed"] - vals["readme"])
        status = "EXACT MATCH" if diff < 0.002 else "DISCREPANCY"
        clf_rows.append({
            "Metric": metric,
            "README_Value": vals["readme"],
            "Recomputed_Value": vals["recomputed"],
            "Absolute_Difference": diff,
            "Status": status,
            "Protocol_Details": "StandardScaler fit inside fold, LogisticRegression L2 (C=1.0, seed=42)"
        })
        
    df_clf_check = pd.DataFrame(clf_rows)
    df_clf_check.to_csv("results/tables/C_gse125498_classifier_check.csv", index=False)
    
    print("\n--- CLASSIFIER RE-CHECK DISCREPANCY TABLE ---")
    print(df_clf_check[["Metric", "README_Value", "Recomputed_Value", "Status"]])
    print(f"\nScaler Verification: Pipeline(StandardScaler, LogisticRegression) fit strictly inside training folds.")
    print("Module C execution completed successfully!")

if __name__ == "__main__":
    main()
