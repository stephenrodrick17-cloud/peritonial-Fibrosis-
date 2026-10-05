"""
Script 61c: Stage 5 GSE182736 Exosome miRNA Processing with Strict Expression Filter
Applies:
  - Expression filter: raw count >= 10 in >= 3 samples across the 6 samples (n=3 non-UF vs n=3 UF failure).
  - Explicit pseudocount: +1.0 added to normalized expression mean before log2 ratio:
      log2FC = log2((mean_norm_UF + 1.0) / (mean_norm_nonUF + 1.0))
  - Explains the extreme values -10.30 and +6.58:
      * -10.30 (hsa-mir-7703-p5_1ss3GA): raw counts control=[45, 2813, 0], disease=[0, 0, 0].
        Only 2 samples have >= 10 counts; filtered OUT by the >=3 samples rule.
      * +6.58 (hsa-miR-892b): raw counts control=[0, 0, 0], disease=[319, 97, 30].
        All 3 disease samples have >= 10 counts; retained by filter.
Outputs:
  - results/tables/F_GSE182736_descriptive_logFC.csv
  - results/tables/F_GSE182736_sample_metadata.csv
  - results/tables/F_GSE182736_provenance.csv
"""

import os
import hashlib
import datetime
import pandas as pd
import numpy as np

ROOT = r"d:\Peritoneal Project"
RAW_DIR = os.path.join(ROOT, "data", "raw")
OUT_DIR = os.path.join(ROOT, "results", "tables")

def main():
    print("=" * 80)
    print("SCRIPT 61c: GSE182736 FILTERED miRNA LOG2FC COMPUTATION")
    print("=" * 80)

    xlsx_path = os.path.join(RAW_DIR, "GSE182736_All_Expressed_miRNA.xlsx")
    sha256_xlsx = hashlib.sha256(open(xlsx_path, "rb").read()).hexdigest()
    print(f"Loading: {xlsx_path}")
    print(f"SHA256:  {sha256_xlsx}")

    # Read sheet with header at row 51
    df = pd.read_excel(xlsx_path, sheet_name="Sheet1", header=51)
    print(f"Total rows in raw Excel: {len(df)}")

    # Restrict to human miRNAs (hsa-)
    mirna_col = "miR_name"
    hsa_mask = df[mirna_col].astype(str).str.startswith("hsa-")
    df_hsa = df[hsa_mask].copy()
    print(f"Total hsa- miRNAs: {len(df_hsa)}")

    raw_cols = ["control_1(raw)", "control_2(raw)", "control_3(raw)", 
                "disease_1(raw)", "disease_2(raw)", "disease_3(raw)"]
    norm_ctrl_cols = ["control_1(norm)", "control_2(norm)", "control_3(norm)"]
    norm_dis_cols  = ["disease_1(norm)", "disease_2(norm)", "disease_3(norm)"]

    for col in raw_cols + norm_ctrl_cols + norm_dis_cols:
        df_hsa[col] = pd.to_numeric(df_hsa[col], errors="coerce").fillna(0.0)

    # -------------------------------------------------------------------------
    # EXPLANATION OF -10.30 AND +6.58 VALUES
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("AUDIT: DETAILED INVESTIGATION OF -10.30 AND +6.58 EXTREME VALUES")
    print("=" * 70)

    # Check -10.30: hsa-mir-7703-p5_1ss3GA
    row_7703 = df_hsa[df_hsa[mirna_col].str.contains("7703")]
    if len(row_7703) > 0:
        r = row_7703.iloc[0]
        c_raw = [r["control_1(raw)"], r["control_2(raw)"], r["control_3(raw)"]]
        d_raw = [r["disease_1(raw)"], r["disease_2(raw)"], r["disease_3(raw)"]]
        c_norm = [r["control_1(norm)"], r["control_2(norm)"], r["control_3(norm)"]]
        d_norm = [r["disease_1(norm)"], r["disease_2(norm)"], r["disease_3(norm)"]]
        m_c_norm = np.mean(c_norm)
        m_d_norm = np.mean(d_norm)
        lfc_unfilt = np.log2((m_d_norm + 1.0) / (m_c_norm + 1.0))
        n_samples_ge10 = sum(x >= 10 for x in c_raw + d_raw)
        print(f"miRNA: {r[mirna_col]}")
        print(f"  Raw counts Control (non-UF):  {c_raw}")
        print(f"  Raw counts Disease (UF fail): {d_raw}")
        print(f"  Samples with raw >= 10 counts: {n_samples_ge10} of 6")
        print(f"  Mean norm: Disease={m_d_norm:.2f}, Control={m_c_norm:.2f}")
        print(f"  Unfiltered log2FC with +1 pseudocount: {lfc_unfilt:.4f} (~ -10.30)")
        print(f"  Filter status (>= 10 in >= 3 samples): {'PASS' if n_samples_ge10 >= 3 else 'FILTERED OUT (FAIL)'}")
        print(f"  Explanation: 7703 has count=0 in all 3 disease samples and 1 control sample.")
        print(f"               Only 2 control samples (45 and 2813) have counts. It fails the >= 3 samples threshold.")

    # Check +6.58: hsa-miR-892b
    row_892b = df_hsa[df_hsa[mirna_col] == "hsa-miR-892b"]
    if len(row_892b) > 0:
        r = row_892b.iloc[0]
        c_raw = [r["control_1(raw)"], r["control_2(raw)"], r["control_3(raw)"]]
        d_raw = [r["disease_1(raw)"], r["disease_2(raw)"], r["disease_3(raw)"]]
        c_norm = [r["control_1(norm)"], r["control_2(norm)"], r["control_3(norm)"]]
        d_norm = [r["disease_1(norm)"], r["disease_2(norm)"], r["disease_3(norm)"]]
        m_c_norm = np.mean(c_norm)
        m_d_norm = np.mean(d_norm)
        lfc_unfilt = np.log2((m_d_norm + 1.0) / (m_c_norm + 1.0))
        n_samples_ge10 = sum(x >= 10 for x in c_raw + d_raw)
        print(f"\nmiRNA: {r[mirna_col]}")
        print(f"  Raw counts Control (non-UF):  {c_raw}")
        print(f"  Raw counts Disease (UF fail): {d_raw}")
        print(f"  Samples with raw >= 10 counts: {n_samples_ge10} of 6")
        print(f"  Mean norm: Disease={m_d_norm:.2f}, Control={m_c_norm:.2f}")
        print(f"  log2FC with +1 pseudocount: {lfc_unfilt:.4f} (~ +6.58)")
        print(f"  Filter status (>= 10 in >= 3 samples): {'PASS' if n_samples_ge10 >= 3 else 'FILTERED OUT (FAIL)'}")
        print(f"  Explanation: 892b is detected robustly in all 3 UF failure samples (319, 97, 30 counts).")
        print(f"               It has zero counts in all 3 controls. With pseudocount +1.0: log2((94.39+1)/1) = +6.58.")

    # -------------------------------------------------------------------------
    # APPLY EXPRESSION FILTER: raw >= 10 in >= 3 samples
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("APPLYING EXPRESSION FILTER (raw >= 10 in >= 3 samples)")
    print("=" * 70)
    pass_filter = (df_hsa[raw_cols] >= 10.0).sum(axis=1) >= 3
    df_filt = df_hsa[pass_filter].copy()

    print(f"Total hsa- miRNAs before filter: {len(df_hsa)}")
    print(f"Total hsa- miRNAs passing filter: {len(df_filt)} (removed {len(df_hsa) - len(df_filt)} low-count miRNAs)")

    # -------------------------------------------------------------------------
    # COMPUTE LOG2 FOLD CHANGE WITH EXPLICIT PSEUDOCOUNT +1.0
    # -------------------------------------------------------------------------
    pseudocount = 1.0
    mean_ctrl_norm = df_filt[norm_ctrl_cols].mean(axis=1)
    mean_dis_norm  = df_filt[norm_dis_cols].mean(axis=1)
    mean_ctrl_raw  = df_filt[["control_1(raw)", "control_2(raw)", "control_3(raw)"]].mean(axis=1)
    mean_dis_raw   = df_filt[["disease_1(raw)", "disease_2(raw)", "disease_3(raw)"]].mean(axis=1)

    log2fc = np.log2((mean_dis_norm + pseudocount) / (mean_ctrl_norm + pseudocount))

    df_out = pd.DataFrame({
        "miRNA": df_filt[mirna_col].astype(str).values,
        "mean_raw_UF_failure": mean_dis_raw.round(2).values,
        "mean_raw_nonUF": mean_ctrl_raw.round(2).values,
        "mean_norm_UF_failure": mean_dis_norm.round(4).values,
        "mean_norm_nonUF": mean_ctrl_norm.round(4).values,
        "pseudocount": pseudocount,
        "log2FC_UF_vs_nonUF": log2fc.round(4).values,
        "direction": np.where(log2fc >= 1.0, "UP", np.where(log2fc <= -1.0, "DOWN", "UNCHANGED")),
        "dataset": "GSE182736",
        "sample_size": "n=3 UF vs n=3 non-UF",
        "analysis_note": "EXPLORATORY_DESCRIPTIVE_ONLY_n3vs3; pseudocount=1.0; filter: raw>=10 in >=3 samples"
    })

    # Sort by log2FC descending
    df_out = df_out.sort_values("log2FC_UF_vs_nonUF", ascending=False).reset_index(drop=True)

    out_csv = os.path.join(OUT_DIR, "F_GSE182736_descriptive_logFC.csv")
    df_out.to_csv(out_csv, index=False)
    print(f"\nSaved filtered table: {out_csv} ({len(df_out)} rows)")

    print("\nTOP 15 UPREGULATED miRNAs in UF FAILURE (passing filter):")
    print(df_out.head(15)[["miRNA", "log2FC_UF_vs_nonUF", "mean_raw_UF_failure", "mean_raw_nonUF", "mean_norm_UF_failure", "mean_norm_nonUF"]].to_string(index=False))

    print("\nTOP 15 DOWNREGULATED miRNAs in UF FAILURE (passing filter):")
    print(df_out.tail(15)[["miRNA", "log2FC_UF_vs_nonUF", "mean_raw_UF_failure", "mean_raw_nonUF", "mean_norm_UF_failure", "mean_norm_nonUF"]].to_string(index=False))

    # Provenance
    prov = pd.DataFrame([{
        "dataset": "GSE182736",
        "raw_file": "GSE182736_All_Expressed_miRNA.xlsx",
        "sha256": sha256_xlsx,
        "header_row": 51,
        "total_probes_in_file": len(df),
        "total_hsa_mirnas": len(df_hsa),
        "filtered_mirnas_universe": len(df_out),
        "filter_applied": "raw_count >= 10 in >= 3 samples",
        "pseudocount": pseudocount,
        "n_up_log2fc_ge1": int((df_out["log2FC_UF_vs_nonUF"] >= 1.0).sum()),
        "n_down_log2fc_le_neg1": int((df_out["log2FC_UF_vs_nonUF"] <= -1.0).sum()),
        "processing_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat()
    }])
    prov.to_csv(os.path.join(OUT_DIR, "F_GSE182736_provenance.csv"), index=False)
    print(f"Saved provenance: {os.path.join(OUT_DIR, 'F_GSE182736_provenance.csv')}")

if __name__ == "__main__":
    main()
