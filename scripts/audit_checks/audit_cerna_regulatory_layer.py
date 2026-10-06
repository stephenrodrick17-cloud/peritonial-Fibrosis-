"""
audit_cerna_regulatory_layer.py
Audit script for the ceRNA regulatory section:
1. Quantifies all miRNA and lncRNA counts across databases and cohorts.
2. Documents all databases, API endpoints, assemblies, versions, and filters applied.
3. Traces each edge's provenance (lncRNA -> miRNA -> Hub Gene).
4. Systematically audits code for hand-written lists, hardcoded heuristics, synthetic identifiers, or unverifiable records.
"""

import os
import re
import glob
import pandas as pd
import numpy as np

ROOT = r"d:\Peritoneal Project"
TABLES_DIR = os.path.join(ROOT, "results", "tables")
SCRIPTS_DIR = os.path.join(ROOT, "scripts")
DATA_DIR = os.path.join(ROOT, "data")

def print_section(title):
    print("\n" + "=" * 80)
    print(f" {title}")
    print("=" * 80)

def main():
    print("STARTING COMPREHENSIVE ceRNA REGULATORY LAYER AUDIT\n")

    # =========================================================================
    # 1. COUNTS: miRNA LAYER (multiMiR, miRTarBase, TarBase, GSE182736, GSE130387)
    # =========================================================================
    print_section("1. miRNA QUANTIFICATION & EVIDENCE TIERS (multiMiR)")

    val_raw_path = os.path.join(TABLES_DIR, "F_multimir_validated_raw.csv")
    pred_raw_path = os.path.join(TABLES_DIR, "F_multimir_predicted_raw.csv")
    tiered_path = os.path.join(TABLES_DIR, "F_hub_to_mirna_tiered.csv")
    tier_cnt_path = os.path.join(TABLES_DIR, "F_hub_to_mirna_tier_counts.csv")

    if os.path.exists(val_raw_path):
        df_val = pd.read_csv(val_raw_path)
        mir_col = "mature_mirna_id" if "mature_mirna_id" in df_val.columns else "mature.mirna.id"
        tgt_col = "target_symbol" if "target_symbol" in df_val.columns else "target.symbol"
        print(f"Raw Validated multiMiR records: {len(df_val)} rows")
        print(f"  Databases present in validated raw: {df_val['database'].value_counts().to_dict()}")
        print(f"  Unique mature miRNAs in validated raw: {df_val[mir_col].nunique()}")
        print(f"  Target genes in validated raw: {df_val[tgt_col].unique().tolist()}")
    else:
        print("F_multimir_validated_raw.csv NOT FOUND")

    if os.path.exists(pred_raw_path):
        df_pred = pd.read_csv(pred_raw_path)
        mir_col_p = "mature_mirna_id" if "mature_mirna_id" in df_pred.columns else "mature.mirna.id"
        print(f"Raw Predicted multiMiR records: {len(df_pred)} rows")
        print(f"  Databases present in predicted raw: {df_pred['database'].value_counts().to_dict()}")
        print(f"  Unique mature miRNAs in predicted raw: {df_pred[mir_col_p].nunique()}")
    else:
        print("F_multimir_predicted_raw.csv NOT FOUND")

    if os.path.exists(tiered_path):
        df_tiered = pd.read_csv(tiered_path)
        print(f"\nUnified Tiered table: {len(df_tiered)} rows")
        print("  Breakdown by Evidence Tier:")
        for tier, grp in df_tiered.groupby("evidence_tier"):
            print(f"    - {tier}: {len(grp)} records, {grp['mature_mirna_id'].nunique()} unique miRNAs, {grp[['target_symbol', 'mature_mirna_id']].drop_duplicates().shape[0]} unique gene-miRNA pairs")
    
    if os.path.exists(tier_cnt_path):
        df_cnt = pd.read_csv(tier_cnt_path)
        print("\nPer-Hub Gene miRNA Counts (F_hub_to_mirna_tier_counts.csv):")
        print(df_cnt.to_string(index=False))

    # =========================================================================
    # 2. COUNTS: lncRNA LAYER (GENCODE v32 & GSE248762)
    # =========================================================================
    print_section("2. lncRNA QUANTIFICATION & SINGLE-CELL EXPRESSION (GSE248762)")

    lnc_det_path = os.path.join(TABLES_DIR, "G_lncRNA_detected_list.csv")
    lnc_sum_path = os.path.join(TABLES_DIR, "G_lncRNA_GSE248762_detection_summary.csv")
    lnc_sens_path = os.path.join(TABLES_DIR, "G_lncRNA_threshold_sensitivity.csv")

    if os.path.exists(lnc_sum_path):
        df_lsum = pd.read_csv(lnc_sum_path)
        print("lncRNA Detection Summary Across GSE248762 Cell Types:")
        print(df_lsum.to_string(index=False))

    if os.path.exists(lnc_det_path):
        df_ldet = pd.read_csv(lnc_det_path)
        print(f"\nDetected lncRNAs table: {len(df_ldet)} rows")
        if "credibility_status" in df_ldet.columns:
            print("  Credibility classification breakdown:")
            print(df_ldet["credibility_status"].value_counts().to_string())
        print("\nTop 10 Detected lncRNAs by Stromal Expression:")
        show_cols = [c for c in ["symbol", "ensembl_id", "pct_stromal_cells_expressing_%", "donors_expressing", "Stromal_total_counts", "credibility_status"] if c in df_ldet.columns]
        print(df_ldet[show_cols].head(10).to_string(index=False))

    if os.path.exists(lnc_sens_path):
        df_lsens = pd.read_csv(lnc_sens_path)
        print("\nThreshold Sensitivity Matrix (G_lncRNA_threshold_sensitivity.csv):")
        print(df_lsens.to_string(index=False))

    # =========================================================================
    # 3. ENCORI (starBase v3.0) QUERY & lncRNA-miRNA INTERACTIONS
    # =========================================================================
    print_section("3. ENCORI LIVE API RETRIEVAL AUDIT")

    encori_raw_path = os.path.join(TABLES_DIR, "G_encori_mirna_lncrna_raw.csv")
    if os.path.exists(encori_raw_path):
        df_encori = pd.read_csv(encori_raw_path)
        print(f"Raw ENCORI interactions: {len(df_encori)} rows")
        print(f"Columns: {df_encori.columns.tolist()}")
        print(f"Queried miRNAs: {df_encori['queried_miRNA'].value_counts().to_dict() if 'queried_miRNA' in df_encori.columns else 'N/A'}")
        print(f"Unique lncRNA symbols (geneName): {df_encori['geneName'].nunique() if 'geneName' in df_encori.columns else 'N/A'}")
        print(f"Unique lncRNA Ensembl IDs (geneID): {df_encori['geneID'].nunique() if 'geneID' in df_encori.columns else 'N/A'}")
        
        # Check distribution of CLIP evidence
        if "clipExpNum" in df_encori.columns:
            print("clipExpNum summary:")
            print(df_encori["clipExpNum"].describe().to_string())
            print(f"Interactions with clipExpNum >= 2: {(df_encori['clipExpNum'] >= 2).sum()}")
            print(f"Interactions with clipExpNum >= 5: {(df_encori['clipExpNum'] >= 5).sum()}")
    else:
        print("G_encori_mirna_lncrna_raw.csv NOT FOUND")

    # =========================================================================
    # 4. PREDICTED ceRNA REGULATORY NETWORK AUDIT
    # =========================================================================
    print_section("4. ceRNA PREDICTED NETWORK (lncRNA -> miRNA -> Hub Gene)")

    net_path = os.path.join(TABLES_DIR, "G_ceRNA_network_predicted.csv")
    sum_path = os.path.join(TABLES_DIR, "G_ceRNA_network_summary.csv")
    alt_path = os.path.join(TABLES_DIR, "G_ceRNA_alternative_top10_tierA.csv")

    if os.path.exists(net_path):
        df_net = pd.read_csv(net_path)
        print(f"Total predicted network rows: {len(df_net)}")
        print(f"Columns: {df_net.columns.tolist()}")
        print(f"Unique lncRNAs: {df_net['lncRNA'].nunique() if 'lncRNA' in df_net.columns else df_net['lncRNA_symbol'].nunique()}")
        print(f"Unique miRNAs: {df_net['miRNA'].nunique() if 'miRNA' in df_net.columns else df_net['miRNA_name'].nunique()}")
        print(f"Unique Hub Genes: {df_net['Hub_Gene'].nunique()}")
        
        if "Evidence_Classification" in df_net.columns:
            print("\nEvidence Classification Breakdown:")
            print(df_net["Evidence_Classification"].value_counts().to_string())
        
        if "Spearman_rho_stromal_pseudobulk" in df_net.columns:
            valid_rho = df_net["Spearman_rho_stromal_pseudobulk"].dropna()
            print(f"\nSpearman correlations computed: {len(valid_rho)} axes")
            print(f"  Positive rho (rho > 0): {(valid_rho > 0).sum()}")
            print(f"  Nominal significant (rho > 0 & P < 0.05): {((df_net['Spearman_rho_stromal_pseudobulk'] > 0) & (df_net['Spearman_pvalue'] < 0.05)).sum()}")
    
    if os.path.exists(sum_path):
        df_sum = pd.read_csv(sum_path)
        print("\nNetwork Summary Table (G_ceRNA_network_summary.csv):")
        print(df_sum.to_string(index=False))

    if os.path.exists(alt_path):
        df_alt = pd.read_csv(alt_path)
        print("\nAlternative Objective Top 10 miRNAs Ranked by Tier A Functional Targets:")
        print(df_alt.to_string(index=False))

    # =========================================================================
    # 5. CODE AUDIT: CHECK FOR HAND-WRITTEN / UNVERIFIABLE / FABRICATED ENTRIES
    # =========================================================================
    print_section("5. CODE AUDIT: HAND-WRITTEN LISTS & UNVERIFIABLE ENTRIES")

    stage5_scripts = glob.glob(os.path.join(SCRIPTS_DIR, "*stage5*.py")) + glob.glob(os.path.join(SCRIPTS_DIR, "*stage5*.R"))
    stage5_scripts += [os.path.join(SCRIPTS_DIR, "70_stage6_generate_figures.py")]
    
    print(f"Auditing {len(stage5_scripts)} scripts for static dictionaries, fake PMIDs, and hardcoded edges:\n")

    suspect_patterns = [
        (r'PMID[\s_]*[:=]\s*["\']\d+', "Hardcoded PMID assignment"),
        (r'HPA\d{6}', "HPA antibody catalog ID"),
        (r'CAB\d{6}', "CAB antibody catalog ID"),
        (r'fake|dummy|synthetic', "Synthetic marker"),
        (r'query_mirnas\s*=\s*\[', "Hardcoded miRNA query list"),
        (r'top_lncs\s*=\s*\[', "Hardcoded top lncRNA list"),
        (r'edges\s*=\s*\[\s*\(', "Static edge list tuple")
    ]

    for s_path in sorted(stage5_scripts):
        if not os.path.exists(s_path):
            continue
        rel_s = os.path.relpath(s_path, ROOT)
        with open(s_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        
        findings = []
        for line_num, line in enumerate(lines, 1):
            for pat, desc in suspect_patterns:
                if re.search(pat, line, re.IGNORECASE):
                    findings.append((line_num, desc, line.strip()))
        
        if findings:
            print(f"File: {rel_s} ({len(findings)} flagged lines)")
            for ln, desc, text in findings:
                print(f"  Line {ln:3d} [{desc}]: {text[:100]}")
            print()
        else:
            print(f"File: {rel_s} -> Clean (0 flagged patterns)")

    # =========================================================================
    # 6. PROVENANCE VERIFICATION OF ENCORI API CALLS
    # =========================================================================
    print_section("6. ENCORI API PROVENANCE VERIFICATION")

    prov_path = os.path.join(TABLES_DIR, "G_ceRNA_provenance.csv")
    if os.path.exists(prov_path):
        df_prov = pd.read_csv(prov_path)
        print("Provenance Table (G_ceRNA_provenance.csv):")
        print(df_prov.to_string(index=False))

    print("\nAUDIT COMPLETE.")

if __name__ == "__main__":
    main()
