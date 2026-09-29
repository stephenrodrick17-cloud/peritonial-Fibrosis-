"""
Automated Repository Audit Script for Peritoneal Fibrosis Project
Usage: python audit/repo_audit.py .
"""

import sys
import os
import glob
import re
import pandas as pd
import numpy as np

def audit_repository(repo_dir):
    repo_dir = os.path.abspath(repo_dir)
    os.chdir(repo_dir)
    
    report_lines = []
    def log(msg="", status=None):
        if status:
            formatted = f"[{status}] {msg}"
        else:
            formatted = msg
        print(formatted)
        report_lines.append(formatted)

    log(f"=====================================================================")
    log(f"AUTOMATED REPOSITORY AUDIT: {repo_dir}")
    log(f"=====================================================================\n")

    fails = []
    warns = []

    # ------------------------------------------------------------------
    # 1. Pipeline Scripts Check
    # ------------------------------------------------------------------
    log("### 1. Core Pipeline Scripts Audit")
    expected_scripts = [
        "00_fetch_gse62928_matrix.R",
        "01_load_qc_preprocess.R",
        "02_differential_expression.R",
        "02b_wgcna_analysis.R",
        "03_matrisome_filtering.R",
        "04_functional_enrichment.R",
        "05b_ml_hub_gene_identification_wgcna.py",
        "06_external_validation_GSE125498.py",
        "analyze_convergence_wgcna.py",
        "statistical_rigor_audit.py",
        "audit_pipeline_errors.py",
        "generate_gse62928_ecm_venn.py"
    ]
    
    for s in expected_scripts:
        found = os.path.exists(s) or os.path.exists(os.path.join("audit", s))
        if found:
            loc = s if os.path.exists(s) else os.path.join("audit", s)
            log(f"Script exists: {loc}", "PASS")
        else:
            log(f"Missing required script: {s}", "FAIL")
            fails.append(f"Missing script: {s}")

    # ------------------------------------------------------------------
    # 2. Hardcoded Absolute Paths & Random Seeds
    # ------------------------------------------------------------------
    log("\n### 2. Code Inspection: Absolute Paths & Seeds")
    all_code_files = glob.glob("*.R") + glob.glob("*.py") + glob.glob("audit/*.py")
    for cf in all_code_files:
        with open(cf, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        
        # Check for hardcoded absolute paths
        abs_paths = re.findall(r'\b[A-Za-z]:[\\\/][A-Za-z0-9_\-\.\/]+', content)
        # filter out benign executable or library references if any
        abs_paths = [p for p in abs_paths if "Rscript" not in p and "python" not in p and not p.endswith(('.R', '.py'))]
        if abs_paths:
            log(f"{cf}: Found hard-coded absolute path(s): {abs_paths[:3]}", "WARN")
            warns.append(f"Hard-coded absolute path in {cf}: {abs_paths[:2]}")
        else:
            log(f"{cf}: No hard-coded absolute paths.", "PASS")

        # Check for random seeds
        has_seed = "set.seed" in content or "random_state" in content or "np.random.seed" in content or "Random.seed" in content
        if has_seed:
            log(f"{cf}: Random seed explicitly set.", "PASS")
        else:
            if "random" in content.lower() or "wgcna" in cf.lower() or "ml" in cf.lower():
                log(f"{cf}: No explicit random seed found in stochastic code.", "WARN")
                warns.append(f"Missing random seed in {cf}")
            else:
                log(f"{cf}: Deterministic script (seed not strictly required).", "PASS")

    # ------------------------------------------------------------------
    # 3. Mendelian Randomization (MR) Remnants Check
    # ------------------------------------------------------------------
    log("\n### 3. Deprecated Mendelian Randomization (MR) Remnants Audit")
    mr_files = [
        "20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt",
        "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
        "pilot_test_results",
        "transcriptome_wide_mr_results"
    ]
    for mf in mr_files:
        if os.path.exists(mf):
            log(f"Deprecated MR artifact present in workspace: {mf}", "FAIL")
            fails.append(f"MR artifact present in workspace: {mf}")
        else:
            log(f"MR artifact cleanly absent: {mf}", "PASS")

    # Check for MR references in active results/tables
    if os.path.exists("convergent_71_ECM_DEGs.csv"):
        df_71 = pd.read_csv("convergent_71_ECM_DEGs.csv")
        mr_cols = [c for c in df_71.columns if "mr" in c.lower()]
        if mr_cols:
            log(f"convergent_71_ECM_DEGs.csv contains leftover MR columns: {mr_cols}", "WARN")
            warns.append(f"Leftover MR columns in convergent_71_ECM_DEGs.csv: {mr_cols}")
        else:
            log("convergent_71_ECM_DEGs.csv has clean columns.", "PASS")

    # ------------------------------------------------------------------
    # 4. Data Tables & Gene Counts Verification (Source of Truth)
    # ------------------------------------------------------------------
    log("\n### 4. Data Table Dimensions & Gene Count Audit")
    
    # GSE62928 expression matrix
    if os.path.exists("results/tables/GSE62928_full_expression_matrix.csv"):
        df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
        n_genes_expr = len(df_expr)
        n_samples_expr = len(df_expr.columns)
        log(f"GSE62928 Expression Matrix: {n_genes_expr} genes x {n_samples_expr} samples", "PASS" if (n_genes_expr==20940 and n_samples_expr==8) else "WARN")
    else:
        log("Missing results/tables/GSE62928_full_expression_matrix.csv", "FAIL")
        fails.append("Missing GSE62928_full_expression_matrix.csv")

    # Salmon module
    if os.path.exists("results/tables/wgcna_trait_significant_module_genes.csv"):
        df_mod = pd.read_csv("results/tables/wgcna_trait_significant_module_genes.csv")
        n_salmon = len(df_mod)
        log(f"WGCNA Salmon module genes: {n_salmon} (Expected: 604)", "PASS" if n_salmon == 604 else "FAIL")
        if n_salmon != 604: fails.append(f"Salmon module count is {n_salmon}, expected 604")
    else:
        log("Missing results/tables/wgcna_trait_significant_module_genes.csv", "FAIL")
        fails.append("Missing wgcna_trait_significant_module_genes.csv")

    # Convergent 40
    if os.path.exists("results/tables/convergent_WGCNA_ECM_genes.csv"):
        df_conv = pd.read_csv("results/tables/convergent_WGCNA_ECM_genes.csv")
        n_conv = len(df_conv)
        log(f"Convergent WGCNA-ECM genes: {n_conv} (Expected: 40)", "PASS" if n_conv == 40 else "FAIL")
        if n_conv != 40: fails.append(f"Convergent gene count is {n_conv}, expected 40")
    else:
        log("Missing results/tables/convergent_WGCNA_ECM_genes.csv", "FAIL")
        fails.append("Missing convergent_WGCNA_ECM_genes.csv")

    # Hub 11
    if os.path.exists("results/tables/ML_hub_genes_from_WGCNA_ECM.csv"):
        df_hubs = pd.read_csv("results/tables/ML_hub_genes_from_WGCNA_ECM.csv")
        n_hubs = len(df_hubs)
        hub_list = df_hubs['Gene_Symbol'].tolist() if 'Gene_Symbol' in df_hubs.columns else []
        expected_11 = ['ISM1', 'FN1', 'EDIL3', 'VCAN', 'COL3A1', 'COMP', 'COL8A1', 'THBS3', 'COL11A1', 'INHBA', 'LOX']
        match_11 = set(hub_list) == set(expected_11)
        log(f"ML Consensus Hub Genes: {n_hubs} genes -> {hub_list}", "PASS" if (n_hubs == 11 and match_11) else "FAIL")
        if not (n_hubs == 11 and match_11): fails.append(f"Hub gene list mismatch: {hub_list}")
    else:
        log("Missing results/tables/ML_hub_genes_from_WGCNA_ECM.csv", "FAIL")
        fails.append("Missing ML_hub_genes_from_WGCNA_ECM.csv")

    # Validation table
    if os.path.exists("results/tables/GSE125498_wgcna_hub_validation_metrics.csv"):
        df_val = pd.read_csv("results/tables/GSE125498_wgcna_hub_validation_metrics.csv")
        n_val = len(df_val)
        log(f"Validation Table Profiled Genes: {n_val} genes (Expected: 7)", "PASS" if n_val == 7 else "FAIL")
        if n_val != 7: fails.append(f"Validation table row count is {n_val}, expected 7")
    else:
        log("Missing results/tables/GSE125498_wgcna_hub_validation_metrics.csv", "FAIL")
        fails.append("Missing GSE125498_wgcna_hub_validation_metrics.csv")

    # ------------------------------------------------------------------
    # 5. Figures Check
    # ------------------------------------------------------------------
    log("\n### 5. Publication Figures Audit")
    expected_figures = [
        "results/figures/graphical_abstract.jpg",
        "results/figures/WGCNA_00_sample_outlier_dendrogram.png",
        "results/figures/WGCNA_01_soft_threshold_selection.png",
        "results/figures/WGCNA_02_gene_dendrogram_modules.png",
        "results/figures/WGCNA_03_module_trait_heatmap.png",
        "results/figures/venn_wgcna_convergence.png",
        "results/figures/WGCNA_ML_01_consensus_votes_barchart.png",
        "results/figures/WGCNA_ML_02_model_selection_heatmap.png",
        "results/figures/WGCNA_ML_03_per_model_importance_2x2.png",
        "results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png",
        "results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png",
        "results/figures/Validation_02_roc_curves_early_vs_late.png",
        "results/figures/Validation_03_stage_progression_trajectories.png",
        "results/figures/Validation_04_patient_cohort_heatmap.png"
    ]
    for fig in expected_figures:
        if os.path.exists(fig):
            size_kb = os.path.getsize(fig) / 1024
            log(f"Figure present: {fig} ({size_kb:.1f} KB)", "PASS")
        else:
            log(f"Missing figure: {fig}", "FAIL")
            fails.append(f"Missing figure: {fig}")

    # ------------------------------------------------------------------
    # 6. Threshold Consistency across Scripts
    # ------------------------------------------------------------------
    log("\n### 6. Threshold Consistency Analysis Across Pipeline Scripts")
    # Check 02_differential_expression.R vs README claims
    with open("02_differential_expression.R", "r", encoding="utf-8", errors="ignore") as f:
        c02 = f.read()
    fc_02 = re.findall(r'logfc_cutoff\s*<-\s*([0-9.]+)', c02)
    pval_02 = re.findall(r'pval_cutoff\s*<-\s*([0-9.]+)', c02)
    log(f"02_differential_expression.R thresholds: logfc_cutoff = {fc_02}, pval_cutoff = {pval_02}")

    with open("02b_wgcna_analysis.R", "r", encoding="utf-8", errors="ignore") as f:
        c02b = f.read()
    if 'ecm_71_file <- "convergent_ECM_DEGs_nominal.csv"' in c02b:
        log("02b_wgcna_analysis.R line 107 loads 'convergent_ECM_DEGs_nominal.csv' (148 genes) into variable ecm_71_file", "WARN")
        warns.append("02b_wgcna_analysis.R variable ecm_71_file points to 148-gene nominal file")

    # ------------------------------------------------------------------
    # 7. Summary of Audit
    # ------------------------------------------------------------------
    log("\n=====================================================================")
    log("AUDIT SUMMARY RESULTS")
    log("=====================================================================")
    log(f"Total FAILs: {len(fails)}")
    for f in fails:
        log(f"  * {f}", "FAIL")
    log(f"Total WARNs: {len(warns)}")
    for w in warns:
        log(f"  * {w}", "WARN")

    # Write report file
    with open("repo_audit_report.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines) + "\n")
    log(f"\nWritten complete audit report to {os.path.join(repo_dir, 'repo_audit_report.md')}")

if __name__ == "__main__":
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    audit_repository(target_dir)
