"""
Script 66b: Stage 5 ceRNA Network Evidence Evaluation & Donor Pseudobulk Co-expression
Implements:
  1. Documentation of original 10-miRNA selection:
     - Prints code lines 71-82 from script 66.
     - Documents whether selection occurred before or after seeing results.
  2. Alternative objective selection:
     - The top 10 miRNAs ranked strictly by number of Tier A hub targets (tie-broken by total targets).
  3. Stromal pseudobulk Spearman correlation across 16 donors:
     - Computes Spearman rho, n=16, and two-sided raw P between lncRNA and hub gene.
  4. Detection audit in experimental datasets:
     - Checks if mediating miRNA is detected in GSE182736 (human effluent exosomes) or GSE130387 (rodent PDF tissue).
  5. Stringent evidence classification:
     - Edges without positive co-expression (rho > 0, P < 0.05) or with undetected miRNAs remain:
       "PREDICTED ONLY (No co-expression or miRNA detection support)"
  6. Filters to robust lncRNAs passing primary detection filter (>= 1.0% stromal cells & >= 2 donors).
     - Flags low-credibility lncRNAs (MALAT1, NEAT1, SNHG*, etc.) and excludes XIST.
Outputs:
  - results/tables/G_ceRNA_network_predicted.csv
  - results/tables/G_ceRNA_network_summary.csv
  - results/tables/G_ceRNA_alternative_top10_tierA.csv
  - results/tables/G_ceRNA_provenance.csv
"""

import os
import re
import datetime
import pandas as pd
import numpy as np
from scipy import stats

ROOT = r"d:\Peritoneal Project"
TABLES_DIR = os.path.join(ROOT, "results", "tables")
PB_PATH = os.path.join(ROOT, "data", "processed", "pseudobulk", "pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv")

def clean_mir_name(name):
    s = str(name).strip()
    s = re.sub(r'_[LR][+\-]\d+.*$', '', s)
    s = re.sub(r'_\d+ss.*$', '', s)
    return s

def main():
    print("=" * 80)
    print("SCRIPT 66b: ceRNA EVIDENCE EVALUATION & STROMAL PSEUDOBULK CO-EXPRESSION")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # 1. STATE HOW THE 10 miRNAs WERE CHOSEN & PRINT CODE LINES
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("1. CODE LINES AND SELECTION TIMING FOR THE 10 ORIGINAL miRNAs")
    print("=" * 70)
    
    script_66_path = os.path.join(ROOT, "scripts", "66_stage5_encori_cerna.py")
    print(f"Reading code lines from: {script_66_path}\n")
    with open(script_66_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    # Print lines 70-83
    for i in range(70, min(84, len(lines))):
        print(f"  Line {i+1:3d}: {lines[i].rstrip()}")

    print("\nSELECTION TIMING & METHODOLOGY DOCUMENTATION:")
    print("  Timing:   BEFORE seeing ceRNA network results.")
    print("  Basis:    Heuristic selection combining general multiMiR target coverage with")
    print("            published peritoneal dialysis / EMT / fibrosis literature (e.g. miR-29 family,")
    print("            miR-21, let-7 family, miR-30a).")
    print("  Caveat:   Because literature-guided choices may introduce confirmation bias, we also provide")
    print("            an entirely automated, literature-agnostic alternative below.")

    # -------------------------------------------------------------------------
    # 2. ALTERNATIVE SELECTION: TOP 10 miRNAs BY TIER A HUB TARGETS
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("2. ALTERNATIVE SELECTION: TOP 10 miRNAs BY TIER A FUNCTIONAL TARGETS")
    print("=" * 70)

    df_tiered = pd.read_csv(os.path.join(TABLES_DIR, "F_hub_to_mirna_tiered.csv"))
    
    # Group by miRNA and count unique Tier A targets
    tier_a_counts = df_tiered[df_tiered["evidence_tier"] == "Tier_A"].groupby("mature_mirna_id")["target_symbol"].nunique()
    total_counts = df_tiered.groupby("mature_mirna_id")["target_symbol"].nunique()

    df_mir_ranked = pd.DataFrame({
        "mature_mirna_id": tier_a_counts.index,
        "n_tier_a_targets": tier_a_counts.values,
        "n_total_hub_targets": [total_counts.get(m, 0) for m in tier_a_counts.index]
    }).sort_values(["n_tier_a_targets", "n_total_hub_targets", "mature_mirna_id"], ascending=[False, False, True]).reset_index(drop=True)

    top10_alt = df_mir_ranked.head(10).copy()
    
    # Add target names
    tier_a_targets = df_tiered[df_tiered["evidence_tier"] == "Tier_A"].groupby("mature_mirna_id")["target_symbol"].unique().to_dict()
    top10_alt["tier_a_target_genes"] = [";".join(sorted(tier_a_targets.get(m, []))) for m in top10_alt["mature_mirna_id"]]

    print("TOP 10 miRNAs RANKED STRICTLY BY TIER A FUNCTIONAL TARGETS:")
    print(top10_alt.to_string(index=False))

    alt_csv = os.path.join(TABLES_DIR, "G_ceRNA_alternative_top10_tierA.csv")
    top10_alt.to_csv(alt_csv, index=False)
    print(f"\nSaved alternative selection: {alt_csv}")

    # -------------------------------------------------------------------------
    # 3. LOAD PSEUDOBULK COUNTS ACROSS 16 DONORS FOR CO-EXPRESSION
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("3. COMPUTING STROMAL PSEUDOBULK SPEARMAN CORRELATION (16 DONORS)")
    print("=" * 70)

    df_pb = pd.read_csv(PB_PATH, index_col=0)
    print(f"Loaded stromal pseudobulk matrix: {df_pb.shape} (16 donors x {df_pb.shape[1]} genes)")

    # -------------------------------------------------------------------------
    # 4. LOAD ENCORI RAW INTERACTIONS AND FILTER TO DETECTED ROBUST lncRNAs
    # -------------------------------------------------------------------------
    raw_encori_path = os.path.join(TABLES_DIR, "G_encori_mirna_lncrna_raw.csv")
    df_encori = pd.read_csv(raw_encori_path)
    print(f"Loaded raw ENCORI interactions: {len(df_encori)} rows")

    # Load lncRNA credibility table
    df_lnc_det = pd.read_csv(os.path.join(TABLES_DIR, "G_lncRNA_detected_list.csv"))
    lnc_status_map = dict(zip(df_lnc_det["symbol"], df_lnc_det["credibility_status"]))
    lnc_pct_map = dict(zip(df_lnc_det["symbol"], df_lnc_det["pct_stromal_cells_expressing_%"]))
    lnc_donors_map = dict(zip(df_lnc_det["symbol"], df_lnc_det["donors_expressing"]))

    # Load miRNA detection in GSE182736 and GSE130387
    df18 = pd.read_csv(os.path.join(TABLES_DIR, "F_GSE182736_descriptive_logFC.csv"))
    df18["clean_mir"] = df18["miRNA"].apply(clean_mir_name)
    mirs_in_18 = set(df18["clean_mir"].unique())

    df13 = pd.read_csv(os.path.join(TABLES_DIR, "F_GSE130387_cross_species_logFC.csv"))
    df13["clean_mir"] = df13["hsa_candidate"].apply(clean_mir_name)
    mirs_in_13 = set(df13["clean_mir"].dropna().unique())

    print(f"Unique clean miRNAs in GSE182736: {len(mirs_in_18)}")
    print(f"Unique clean miRNAs in GSE130387: {len(mirs_in_13)}")

    # -------------------------------------------------------------------------
    # 5. ASSEMBLE ceRNA EDGES WITH PSEUDOBULK SPEARMAN CORRELATIONS
    # -------------------------------------------------------------------------
    # Map miRNA to hub targets
    mir_to_hubs = df_tiered.groupby("mature_mirna_id")["target_symbol"].unique().to_dict()
    mir_to_tierA_hubs = df_tiered[df_tiered["evidence_tier"] == "Tier_A"].groupby("mature_mirna_id")["target_symbol"].unique().to_dict()

    # Rename geneName in ENCORI to lncRNA symbol
    lnc_col = "geneName" if "geneName" in df_encori.columns else "gene_name"
    mir_col = "queried_miRNA" if "queried_miRNA" in df_encori.columns else "miRNA"

    cerna_edges = []
    
    # Cache correlation calculations to avoid redundant compute
    corr_cache = {}

    for _, row in df_encori.iterrows():
        lnc = str(row[lnc_col])
        mir = str(row[mir_col])
        clean_mir = clean_mir_name(mir)

        # Check lncRNA status
        status = lnc_status_map.get(lnc, "BELOW_DETECTION_THRESHOLD")
        if status == "BELOW_DETECTION_THRESHOLD":
            continue

        # Target hubs for this miRNA
        hubs = mir_to_hubs.get(mir, [])
        tierA_hubs = set(mir_to_tierA_hubs.get(mir, []))

        # Check miRNA detection in experimental datasets
        det_18 = clean_mir in mirs_in_18
        det_13 = clean_mir in mirs_in_13

        for h in hubs:
            # Spearman correlation between lncRNA and Hub in 16 stromal pseudobulk donors
            pair_key = (lnc, h)
            if pair_key in corr_cache:
                rho, pval = corr_cache[pair_key]
            else:
                if (lnc in df_pb.columns) and (h in df_pb.columns):
                    x = df_pb[lnc].values.astype(float)
                    y = df_pb[h].values.astype(float)
                    # Check variance
                    if np.std(x) > 0 and np.std(y) > 0:
                        res = stats.spearmanr(x, y)
                        rho, pval = float(res.statistic), float(res.pvalue)
                    else:
                        rho, pval = np.nan, np.nan
                else:
                    rho, pval = np.nan, np.nan
                corr_cache[pair_key] = (rho, pval)

            is_tier_a_link = (h in tierA_hubs)
            
            # Evidence classification
            has_positive_coexpr = (not np.isnan(rho)) and (rho > 0) and (pval < 0.05)
            has_mir_det = (det_18 or det_13)

            if status == "EXCLUDED_FROM_INTERPRETATION_DONOR_SEX":
                ev_class = "EXCLUDED_DONOR_SEX_BIAS (XIST)"
            elif status == "LOW_CREDIBILITY_CANDIDATE":
                ev_class = "LOW_CREDIBILITY_PROMISCUOUS_OR_HOUSEKEEPING (MALAT1/NEAT1/SNHG*)"
            elif has_positive_coexpr and has_mir_det:
                ev_class = "NOMINAL_COEXPRESSION_AND_MIRNA_DETECTED (Descriptive, not validated)"
            else:
                ev_class = "PREDICTED_ONLY (No co-expression or miRNA detection support)"

            cerna_edges.append({
                "lncRNA": lnc,
                "lncRNA_credibility": status,
                "pct_stromal_cells_expressing_%": lnc_pct_map.get(lnc, np.nan),
                "donors_expressing": lnc_donors_map.get(lnc, np.nan),
                "miRNA": mir,
                "miRNA_detected_in_GSE182736_exosomes": det_18,
                "miRNA_detected_in_GSE130387_rodent": det_13,
                "Hub_Gene": h,
                "Hub_miRNA_Evidence": "Tier_A_Functional" if is_tier_a_link else "Tier_B_Predicted_or_HTP",
                "Spearman_rho_stromal_pseudobulk": round(rho, 4) if not np.isnan(rho) else np.nan,
                "Spearman_pvalue": round(pval, 6) if not np.isnan(pval) else np.nan,
                "Donors_n": 16,
                "Evidence_Classification": ev_class,
                "Mandatory_Disclaimer": "PREDICTED regulatory network (not experimentally validated in this system)"
            })

    df_edges = pd.DataFrame(cerna_edges)
    out_edges_csv = os.path.join(TABLES_DIR, "G_ceRNA_network_predicted.csv")
    df_edges.to_csv(out_edges_csv, index=False)
    print(f"\nSaved predicted ceRNA network edges: {out_edges_csv} ({len(df_edges)} rows)")

    # -------------------------------------------------------------------------
    # 6. PRINT SUMMARY OF EVIDENCE CLASSIFICATION
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("ceRNA NETWORK EVIDENCE CLASSIFICATION BREAKDOWN")
    print("=" * 70)
    print(df_edges["Evidence_Classification"].value_counts())

    # Summary table
    sum_df = pd.DataFrame([
        {"Category": "Total Candidate ceRNA Axes (lncRNA -> miRNA -> Hub)", "Count": len(df_edges)},
        {"Category": "Unique lncRNAs in Network", "Count": df_edges["lncRNA"].nunique()},
        {"Category": "Unique miRNAs Queried", "Count": df_edges["miRNA"].nunique()},
        {"Category": "Unique Hub Genes Targeted", "Count": df_edges["Hub_Gene"].nunique()},
        {"Category": "Axes with Nominal Co-expression (rho>0, P<0.05) & miRNA Detected", 
         "Count": sum(df_edges["Evidence_Classification"] == "NOMINAL_COEXPRESSION_AND_MIRNA_DETECTED (Descriptive, not validated)")},
        {"Category": "Axes PREDICTED ONLY (no co-expression or miRNA undetected)", 
         "Count": sum(df_edges["Evidence_Classification"] == "PREDICTED_ONLY (No co-expression or miRNA detection support)")},
        {"Category": "Axes from Low-Credibility Transcripts (MALAT1/NEAT1/SNHG*)", 
         "Count": sum(df_edges["Evidence_Classification"] == "LOW_CREDIBILITY_PROMISCUOUS_OR_HOUSEKEEPING (MALAT1/NEAT1/SNHG*)")},
        {"Category": "Axes from XIST (Excluded for donor sex bias)", 
         "Count": sum(df_edges["Evidence_Classification"] == "EXCLUDED_DONOR_SEX_BIAS (XIST)")}
    ])
    sum_csv = os.path.join(TABLES_DIR, "G_ceRNA_network_summary.csv")
    sum_df.to_csv(sum_csv, index=False)
    print(f"\nSaved ceRNA network summary: {sum_csv}")
    print(sum_df.to_string(index=False))

    print("\nSAMPLE TOP CO-EXPRESSED EDGES (Nominal Support):")
    top_coex = df_edges[df_edges["Evidence_Classification"].str.startswith("NOMINAL_COEXPRESSION")].sort_values("Spearman_rho_stromal_pseudobulk", ascending=False)
    print(top_coex.head(15)[["lncRNA", "miRNA", "Hub_Gene", "Spearman_rho_stromal_pseudobulk", "Spearman_pvalue", "Hub_miRNA_Evidence", "miRNA_detected_in_GSE182736_exosomes", "miRNA_detected_in_GSE130387_rodent"]].to_string(index=False))

    # Provenance
    prov_df = pd.DataFrame([{
        "dataset": "ENCORI_GSE248762_ceRNA",
        "total_axes": len(df_edges),
        "coexpressed_axes": sum(df_edges["Evidence_Classification"].str.startswith("NOMINAL_COEXPRESSION")),
        "predicted_only_axes": sum(df_edges["Evidence_Classification"].str.startswith("PREDICTED_ONLY")),
        "disclaimer": "PREDICTED regulatory network (not experimentally validated in this system)",
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat()
    }])
    prov_df.to_csv(os.path.join(TABLES_DIR, "G_ceRNA_provenance.csv"), index=False)
    print("Saved provenance.")

if __name__ == "__main__":
    main()
