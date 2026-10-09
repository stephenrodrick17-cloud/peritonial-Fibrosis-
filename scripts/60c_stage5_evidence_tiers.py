"""
Script 60c: Stage 5 Evidence Tiers for Hub Gene - miRNA Interactions
Classifies every hub-miRNA record:
  - Tier A: miRTarBase functional evidence (reporter assay, luciferase, western blot, qPCR, qRT-PCR, immunoblot)
  - Tier B: everything else (TarBase high-throughput CLIP/sequencing, miRecords, predicted from multiMiR)
Outputs:
  - results/tables/F_hub_to_mirna_tiered.csv
  - results/tables/F_hub_to_mirna_tier_counts.csv
Prints counts per tier and per hub gene.
"""

import os
import pandas as pd
import numpy as np

ROOT = r"d:\Peritoneal Project"
TABLES_DIR = os.path.join(ROOT, "results", "tables")

def main():
    print("=" * 80)
    print("SCRIPT 60c: STAGE 5 EVIDENCE TIERS CLASSIFICATION")
    print("=" * 80)

    # 1. Load multiMiR validated and predicted tables
    val_path = os.path.join(TABLES_DIR, "F_multimir_validated_raw.csv")
    pred_path = os.path.join(TABLES_DIR, "F_multimir_predicted_raw.csv")

    df_val = pd.read_csv(val_path)
    df_pred = pd.read_csv(pred_path)

    print(f"Loaded validated raw records: {len(df_val)}")
    print(f"Loaded predicted raw records: {len(df_pred)}")

    # 2. Define functional evidence classifier for Tier A
    functional_keywords = [
        "reporter assay", "luciferase", "western blot", 
        "qpcr", "qrt-pcr", "immunoblot"
    ]

    def classify_record(row):
        db = str(row.get("database", "")).lower()
        exp = str(row.get("experiment", "")).lower()
        if db == "mirtarbase":
            if any(k in exp for k in functional_keywords):
                return "Tier_A"
        return "Tier_B"

    df_val["evidence_tier"] = df_val.apply(classify_record, axis=1)

    # All predicted records are Tier B
    df_pred["evidence_tier"] = "Tier_B"
    df_pred["experiment"] = "Bioinformatic prediction (multiMiR)"

    # Combine into unified tiered table
    cols = ["target_symbol", "mature_mirna_id", "database", "evidence_tier", "experiment", "type"]
    val_sub = df_val[cols].copy()
    pred_sub = df_pred[["target_symbol", "mature_mirna_id", "database", "evidence_tier", "experiment", "type"]].copy()

    df_combined = pd.concat([val_sub, pred_sub], ignore_index=True).drop_duplicates()

    # Save tiered table
    tiered_out_path = os.path.join(TABLES_DIR, "F_hub_to_mirna_tiered.csv")
    df_combined.to_csv(tiered_out_path, index=False)
    print(f"\nSaved tiered interactions: {tiered_out_path} ({len(df_combined)} rows)")

    # 3. Summarize counts per tier
    print("\n" + "=" * 60)
    print("EVIDENCE TIER SUMMARY")
    print("=" * 60)
    
    tier_summary = []
    for tier in ["Tier_A", "Tier_B"]:
        sub = df_combined[df_combined["evidence_tier"] == tier]
        unique_mirs = sub["mature_mirna_id"].nunique()
        unique_pairs = sub[["target_symbol", "mature_mirna_id"]].drop_duplicates().shape[0]
        total_records = len(sub)
        tier_summary.append({
            "Evidence_Tier": tier,
            "Description": "miRTarBase functional assays (reporter/luciferase/western/qPCR)" if tier == "Tier_A" else "High-throughput CLIP/microarray/predicted",
            "Unique_miRNAs": unique_mirs,
            "Unique_Hub_miRNA_Pairs": unique_pairs,
            "Total_Records": total_records
        })

    df_tier_sum = pd.DataFrame(tier_summary)
    print(df_tier_sum.to_string(index=False))

    # 4. Summarize counts per hub gene
    print("\n" + "=" * 60)
    print("COUNTS PER HUB GENE AND EVIDENCE TIER")
    print("=" * 60)
    
    hub_genes = ["ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP",
                 "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"]

    hub_rows = []
    for h in hub_genes:
        sub_h = df_combined[df_combined["target_symbol"] == h]
        tA_mirs = sub_h[sub_h["evidence_tier"] == "Tier_A"]["mature_mirna_id"].nunique()
        tB_mirs = sub_h[sub_h["evidence_tier"] == "Tier_B"]["mature_mirna_id"].nunique()
        total_mirs = sub_h["mature_mirna_id"].nunique()
        hub_rows.append({
            "Hub_Gene": h,
            "Tier_A_Functional_miRNAs": tA_mirs,
            "Tier_B_Other_miRNAs": tB_mirs,
            "Total_Unique_miRNAs": total_mirs
        })

    df_hub_counts = pd.DataFrame(hub_rows).sort_values("Tier_A_Functional_miRNAs", ascending=False)
    print(df_hub_counts.to_string(index=False))

    hub_counts_path = os.path.join(TABLES_DIR, "F_hub_to_mirna_tier_counts.csv")
    df_hub_counts.to_csv(hub_counts_path, index=False)
    print(f"\nSaved per-hub tier counts: {hub_counts_path}")

    # Also list the Tier A miRNAs and their validated targets
    print("\n" + "=" * 60)
    print("ALL TIER A FUNCTIONAL miRTarBase INTERACTIONS (n=24 pairs)")
    print("=" * 60)
    tier_a_pairs = df_combined[df_combined["evidence_tier"] == "Tier_A"][["target_symbol", "mature_mirna_id", "experiment"]].drop_duplicates()
    print(tier_a_pairs.sort_values(["target_symbol", "mature_mirna_id"]).to_string(index=False))

if __name__ == "__main__":
    main()
