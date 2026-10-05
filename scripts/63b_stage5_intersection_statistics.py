"""
Script 63b: Stage 5 Comprehensive miRNA Intersection Statistics and Directional Audit
Computes:
  - Universe N: miRNAs measured and passing expression filter.
  - K: Hub-targeting miRNAs in universe (Tier A primary; Tier B sensitivity).
  - n: DE miRNAs in universe (|log2FC| >= 1.0 primary; 0.58 and 2.0 sensitivities).
  - k: Observed overlap.
  - k_exp = n * K / N: Expected overlap under the null.
  - Enrichment ratio: k / k_exp.
  - Both one-sided P-values:
      * P_enrichment = P(X >= k) = hypergeom.sf(k - 1, N, K, n)
      * P_depletion  = P(X <= k) = hypergeom.cdf(k, N, K, n)
  - Direction analysis: Checks whether overlapping miRNA log2FC is opposite to hub-gene upregulation in tissue.
  - Explicit labeling: Every P-value labeled EXPLORATORY (n=3 vs 3, no variance model).
  - Replaces previous inaccurate statement that "no P-values exist".
Outputs:
  - results/tables/F_mirna_intersection_summary.csv
  - results/tables/F_mirna_hub_intersection.csv
"""

import os
import re
import pandas as pd
import numpy as np
from scipy import stats

ROOT = r"d:\Peritoneal Project"
TABLES_DIR = os.path.join(ROOT, "results", "tables")

def clean_mir_name(name):
    s = str(name).strip()
    s = re.sub(r'_[LR][+\-]\d+.*$', '', s)
    s = re.sub(r'_\d+ss.*$', '', s)
    return s

def main():
    print("=" * 80)
    print("SCRIPT 63b: STAGE 5 miRNA INTERSECTION STATISTICS & DIRECTIONAL AUDIT")
    print("=" * 80)

    # 1. Load Tiered Hub-miRNA Interactions
    tiered_path = os.path.join(TABLES_DIR, "F_hub_to_mirna_tiered.csv")
    df_tiered = pd.read_csv(tiered_path)

    # Tier A unique clean miRNAs
    tier_a_raw = df_tiered[df_tiered["evidence_tier"] == "Tier_A"]["mature_mirna_id"].unique()
    tier_a_mirs = set(clean_mir_name(m) for m in tier_a_raw)
    print(f"Tier A Functional Hub miRNAs: {len(tier_a_mirs)} unique clean miRNAs")

    # Tier B all miRNAs (validated + predicted)
    tier_b_raw = df_tiered["mature_mirna_id"].unique()
    tier_b_mirs = set(clean_mir_name(m) for m in tier_b_raw)
    print(f"Tier B (All Validated + Predicted) Hub miRNAs: {len(tier_b_mirs)} unique clean miRNAs")

    # Target map for reporting
    hub_targets_a = df_tiered[df_tiered["evidence_tier"] == "Tier_A"].groupby("mature_mirna_id")["target_symbol"].unique().to_dict()
    clean_target_map_a = {}
    for m, targets in hub_targets_a.items():
        cm = clean_mir_name(m)
        if cm not in clean_target_map_a:
            clean_target_map_a[cm] = set()
        clean_target_map_a[cm].update(targets)

    # 2. Load GSE182736 (Human Effluent Exosomes, Filtered)
    df18 = pd.read_csv(os.path.join(TABLES_DIR, "F_GSE182736_descriptive_logFC.csv"))
    df18["clean_miRNA"] = df18["miRNA"].apply(clean_mir_name)
    u18 = df18.groupby("clean_miRNA")["log2FC_UF_vs_nonUF"].mean()
    N_18 = len(u18)

    # 3. Load GSE130387 (Rodent PDF vs Saline Tissue, Filtered)
    df13 = pd.read_csv(os.path.join(TABLES_DIR, "F_GSE130387_cross_species_logFC.csv"))
    df13["clean_miRNA"] = df13["hsa_candidate"].apply(clean_mir_name)
    u13 = df13.groupby("clean_miRNA")["log2FC_PDF_vs_saline"].mean()
    N_13 = len(u13)

    print(f"\nUniverse sizes after expression filters:")
    print(f"  GSE182736 (Human Effluent Exosomes): N = {N_18} miRNAs")
    print(f"  GSE130387 (Rodent Tissue PDF):        N = {N_13} miRNAs")

    # 4. Compute Hypergeometric Statistics Across Datasets and Thresholds
    datasets = [
        {"name": "GSE182736_Human_Effluent_Exosomes", "series": u18, "N": N_18, "species": "Homo sapiens", "matrix": "Effluent Exosomes (UF fail vs non-UF)"},
        {"name": "GSE130387_Rodent_PDF_Tissue",        "series": u13, "N": N_13, "species": "CROSS-SPECIES (Mouse tissue / Rat probes)", "matrix": "Peritoneal Tissue (PDF vs Saline)"}
    ]

    fc_cutoffs = [1.0, 0.58, 2.0]
    tier_configs = [
        ("Tier_A_Primary", tier_a_mirs),
        ("Tier_B_Sensitivity", tier_b_mirs)
    ]

    summary_rows = []
    overlap_detail_rows = []

    for d in datasets:
        dname = d["name"]
        u_series = d["series"]
        N = d["N"]

        for tier_name, hub_mirs in tier_configs:
            # Universe intersection with hub list
            k_in_universe = set(u_series.index).intersection(hub_mirs)
            K = len(k_in_universe)

            for fc_cut in fc_cutoffs:
                de_mirs = set(u_series[u_series.abs() >= fc_cut].index)
                n = len(de_mirs)
                overlap = de_mirs.intersection(hub_mirs)
                k = len(overlap)

                k_exp = (n * K) / N if N > 0 else 0.0
                ratio = k / k_exp if k_exp > 0 else np.nan

                # One-sided tests
                p_enrich = stats.hypergeom.sf(k - 1, N, K, n) if (N > 0 and K > 0 and n > 0) else 1.0
                p_deplet = stats.hypergeom.cdf(k, N, K, n) if (N > 0 and K > 0 and n > 0) else 1.0

                summary_rows.append({
                    "Dataset": dname,
                    "Species": d["species"],
                    "Matrix": d["matrix"],
                    "Evidence_Tier": tier_name,
                    "DE_Threshold_log2FC": fc_cut,
                    "Universe_N": N,
                    "Hub_miRNAs_in_Universe_K": K,
                    "DE_miRNAs_n": n,
                    "Observed_Overlap_k": k,
                    "Expected_Overlap_k_exp": round(k_exp, 3),
                    "Enrichment_Ratio": round(ratio, 3) if not np.isnan(ratio) else np.nan,
                    "One_Sided_P_Enrichment": round(p_enrich, 6),
                    "One_Sided_P_Depletion": round(p_deplet, 6),
                    "Exploratory_P_Status": "EXPLORATORY_DESCRIPTIVE_ONLY (n=3 vs 3, no negative binomial dispersion model)",
                    "Overlapping_miRNAs": ";".join(sorted(overlap)) if k > 0 else "NONE"
                })

                # Record details for overlapping miRNAs
                if tier_name == "Tier_A_Primary" and fc_cut in [1.0, 0.58]:
                    for m in sorted(overlap):
                        lfc = u_series[m]
                        targs = clean_target_map_a.get(m, set())
                        # Hub genes are UPREGULATED in peritoneal fibrosis
                        # Repression expectation: miRNA should be DOWNREGULATED (opposite to hub)
                        is_opposite = (lfc < 0)
                        rep_status = "OPPOSITE_TO_HUB (Consistent with loss of repression)" if is_opposite else "SAME_DIRECTION_AS_HUB (Inconsistent with simple repression)"
                        overlap_detail_rows.append({
                            "Dataset": dname,
                            "Species": d["species"],
                            "Evidence_Tier": tier_name,
                            "DE_Cutoff_log2FC": fc_cut,
                            "miRNA": m,
                            "log2FC": round(lfc, 4),
                            "miRNA_Direction": "UP" if lfc > 0 else "DOWN",
                            "Hub_Direction_In_Tissue": "UP (Pro-fibrotic ECM)",
                            "Direction_Consistency": rep_status,
                            "Validated_Hub_Targets": ";".join(sorted(targs)),
                            "Inference_Label": "EXPLORATORY_ONLY (n=3 vs 3)"
                        })

    df_summary = pd.DataFrame(summary_rows)
    summary_path = os.path.join(TABLES_DIR, "F_mirna_intersection_summary.csv")
    df_summary.to_csv(summary_path, index=False)
    print(f"\nSaved intersection summary: {summary_path}")

    df_details = pd.DataFrame(overlap_detail_rows).drop_duplicates()
    details_path = os.path.join(TABLES_DIR, "F_mirna_hub_intersection.csv")
    df_details.to_csv(details_path, index=False)
    print(f"Saved intersection details: {details_path}")

    # 5. Print Primary Statistics
    print("\n" + "=" * 100)
    print("TABLE 3A: PRIMARY INTERSECTION STATISTICS (Tier A Functional miRTarBase)")
    print("=" * 100)
    sub_primary = df_summary[df_summary["Evidence_Tier"] == "Tier_A_Primary"]
    print(sub_primary[["Dataset", "DE_Threshold_log2FC", "Universe_N", "Hub_miRNAs_in_Universe_K", "DE_miRNAs_n", "Observed_Overlap_k", "Expected_Overlap_k_exp", "Enrichment_Ratio", "One_Sided_P_Enrichment", "One_Sided_P_Depletion", "Overlapping_miRNAs"]].to_string(index=False))

    print("\n" + "=" * 100)
    print("TABLE 3B: SENSITIVITY INTERSECTION STATISTICS (Tier B All Validated + Predicted)")
    print("=" * 100)
    sub_tier_b = df_summary[df_summary["Evidence_Tier"] == "Tier_B_Sensitivity"]
    print(sub_tier_b[["Dataset", "DE_Threshold_log2FC", "Universe_N", "Hub_miRNAs_in_Universe_K", "DE_miRNAs_n", "Observed_Overlap_k", "Expected_Overlap_k_exp", "Enrichment_Ratio", "One_Sided_P_Enrichment", "One_Sided_P_Depletion"]].to_string(index=False))

    print("\n" + "=" * 100)
    print("TABLE 3C: DIRECTIONAL CONCORDANCE AUDIT (Repression Expectation vs Observation)")
    print("=" * 100)
    print(df_details[["Dataset", "DE_Cutoff_log2FC", "miRNA", "log2FC", "miRNA_Direction", "Hub_Direction_In_Tissue", "Direction_Consistency", "Validated_Hub_Targets"]].to_string(index=False))

if __name__ == "__main__":
    main()
