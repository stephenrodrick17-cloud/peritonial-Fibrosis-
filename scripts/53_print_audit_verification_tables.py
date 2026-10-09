"""
Script 53: Print Audit Verification Tables Verbatim via df.to_string()
Reads directly from manifest files:
1. BH comparison: computed vs statsmodels multipletests (max difference)
2. (i) and (j) count-matrix hashes, provenance, edgeR code, and SEs
3. Mathematical proof of ambient CP10k change (~3.6x in stromal/monocyte vs T/NK)
4. Non-stromal primary table (132 rows)
5. Stromal LODO table (144 rows)
6. Sensitivities (b)-(f) for all 11 genes (165 rows)
7. State which conclusions change
"""
import os
import json
import hashlib
import datetime
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def get_file_info(filepath):
    abs_path = os.path.abspath(filepath)
    stat = os.stat(abs_path)
    mtime_utc = datetime.datetime.fromtimestamp(stat.st_mtime, tz=datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
    h = hashlib.sha256()
    with open(abs_path, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return abs_path, mtime_utc, h.hexdigest(), stat.st_size

def main():
    print("=" * 100)
    print("STAGE 4 CONSOLIDATION AUDIT VERIFICATION & VERBATIM MANIFEST READBACK")
    print("=" * 100)
    
    # -------------------------------------------------------------------------
    # 1. BH Independent Calculation with Statsmodels
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("SECTION 1: INDEPENDENT BENJAMINI-HOCHBERG FDR VALIDATION (statsmodels vs edgeR)")
    print("=" * 100)
    bh_path = "results/tables/stage4_primary_96_bh_family_sorted.csv"
    p, mtime, sha, sz = get_file_info(bh_path)
    df_bh = pd.read_csv(bh_path)
    print(f"File: {p}")
    print(f"mtime: {mtime} | SHA256: {sha} | Shape: {df_bh.shape}")
    
    # Run statsmodels multipletests
    pvals = df_bh["PValue"].values
    reject, sm_fdr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    
    diff_computed_sm = np.abs(df_bh["computed_BH_FDR"].values - sm_fdr)
    diff_edger_sm = np.abs(df_bh["BH_FDR"].values - sm_fdr)
    
    print(f"\nMax difference between step-by-step computed_BH_FDR and statsmodels multipletests: {np.max(diff_computed_sm):.2e}")
    print(f"Max difference between edgeR stored BH_FDR and statsmodels multipletests:         {np.max(diff_edger_sm):.2e}")
    
    print("\nTop 15 Tests in 96-Test Family (sorted strictly by P-value):")
    cols_bh_show = ["rank", "Cell_Type", "Contrast", "Gene", "PValue", "raw_BH_ratio", "BH_FDR", "computed_BH_FDR", "statsmodels_BH_FDR"]
    print(df_bh[cols_bh_show].head(15).to_string())

    # -------------------------------------------------------------------------
    # 2. Sensitivities (i) and (j): Count Matrix Hashes & edgeR Fits
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("SECTION 2: SENSITIVITIES (i) AND (j) - COUNT MATRIX PROVENANCE & edgeR FITS")
    print("=" * 100)
    
    matrix_files = [
        "data/processed/pseudobulk/pb_sens_i_no_upper_ceilings_counts.csv",
        "data/processed/pseudobulk/pb_sens_i_no_upper_ceilings_metadata.csv",
        "data/processed/pseudobulk/pb_sens_j_pooled_fixed_ceiling_counts.csv",
        "data/processed/pseudobulk/pb_sens_j_pooled_fixed_ceiling_metadata.csv"
    ]
    for mf in matrix_files:
        mf_p, mf_mt, mf_sha, mf_sz = get_file_info(mf)
        df_tmp = pd.read_csv(mf)
        print(f"Matrix File: {mf}")
        print(f"  Path:   {mf_p}")
        print(f"  mtime:  {mf_mt} | Size: {mf_sz} bytes | Shape: {df_tmp.shape}")
        print(f"  SHA256: {mf_sha}\n")
    
    # Load sensitivity results
    sens_path = "results/tables/stage4_sensitivity_edger_pseudobulk.csv"
    sp, smt, ssha, ssz = get_file_info(sens_path)
    df_sens = pd.read_csv(sens_path)
    print(f"Results File: {sp}")
    print(f"mtime: {smt} | SHA256: {ssha} | Shape: {df_sens.shape}")
    
    eval_genes = ["THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", "EDIL3", "LOX", "INHBA", "COMP"]
    cols_fit = ["Analysis", "Contrast", "Gene", "log2FC", "SE", "CI_95_low", "CI_95_high", "t_stat", "PValue", "QL_Denom_df", "Donors_LV_UF", "Donors_LV_NOT_UF", "Cells_total"]
    
    print("\nSensitivity (i) [No Upper Ceilings, scrublet+scDblFinder removed] - LV_UF vs LV_NOT_UF:")
    df_i = df_sens[(df_sens["Analysis"] == "sens_i_no_upper_ceilings") & 
                   (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
                   (df_sens["Gene"].isin(eval_genes))]
    print(df_i[cols_fit].to_string(index=False))
    
    print("\nSensitivity (j) [Pooled 99.5th Percentile Ceiling] - LV_UF vs LV_NOT_UF:")
    df_j = df_sens[(df_sens["Analysis"] == "sens_j_pooled_fixed_ceiling") & 
                   (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
                   (df_sens["Gene"].isin(eval_genes))]
    print(df_j[cols_fit].to_string(index=False))

    # Side-by-side comparison table
    comp_path = "results/tables/stage4_qc_ceiling_sensitivity_comparison.csv"
    cp, cmt, csha, csz = get_file_info(comp_path)
    df_comp = pd.read_csv(comp_path)
    print(f"\nQC Ceiling Sensitivity Comparison Table ({cp}):")
    print(f"mtime: {cmt} | SHA256: {csha}")
    print(df_comp.to_string(index=False))

    # -------------------------------------------------------------------------
    # 3. Ambient Profiling Table Readback
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("SECTION 3: AMBIENT RNA PROFILING RECOMPUTATION (True Raw Counts / Total UMI)")
    print("=" * 100)
    amb_path = "results/tables/stage4_ambient_recomputed_true_cp10k.csv"
    ap, amt, asha, asz = get_file_info(amb_path)
    df_amb = pd.read_csv(amb_path)
    print(f"File: {ap}")
    print(f"mtime: {amt} | SHA256: {asha} | Shape: {df_amb.shape}")
    print(df_amb.to_string(index=False))

    # -------------------------------------------------------------------------
    # 4. Complete Non-Stromal Primary Table (132 Rows)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("SECTION 4: COMPLETE NON-STROMAL PRIMARY TABLE (132 Rows: 4 Cell Types x 3 Contrasts x 11 Genes)")
    print("=" * 100)
    prim_path = "results/tables/stage4_primary_edger_pseudobulk.csv"
    pp, pmt, psha, psz = get_file_info(prim_path)
    df_prim = pd.read_csv(prim_path)
    print(f"File: {pp}")
    print(f"mtime: {pmt} | SHA256: {psha} | Shape: {df_prim.shape}")
    
    df_non_stromal = df_prim[~df_prim["Cell_Type"].str.contains("stromal")].copy()
    print(f"Non-Stromal Subset Rows: {len(df_non_stromal)}")
    cols_prim = ["Cell_Type", "Contrast", "Gene", "log2FC", "SE", "CI_95_low", "CI_95_high", "t_stat", "PValue", "BH_FDR", "QL_Denom_df", "Evidence_Status"]
    print(df_non_stromal[cols_prim].to_string(index=False))

    # -------------------------------------------------------------------------
    # 5. Stromal Leave-One-Donor-Out (LODO) Table (144 Rows)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("SECTION 5: STROMAL LEAVE-ONE-DONOR-OUT (LODO) TABLE (144 Rows: 16 Dropped Donors x 9 Evaluable Genes)")
    print("=" * 100)
    df_lodo = df_sens[(df_sens["Analysis"].str.startswith("lodo_")) & 
                      (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
                      (df_sens["Gene"].isin(eval_genes))].copy()
    df_lodo = df_lodo.sort_values(["Analysis", "Gene"]).reset_index(drop=True)
    print(f"Stromal LODO Subset Rows: {len(df_lodo)}")
    cols_lodo = ["Analysis", "Contrast", "Gene", "log2FC", "SE", "CI_95_low", "CI_95_high", "t_stat", "PValue", "QL_Denom_df", "Donors_LV_UF", "Donors_LV_NOT_UF"]
    print(df_lodo[cols_lodo].to_string(index=False))

    # -------------------------------------------------------------------------
    # 6. Sensitivities (b)-(f) for All 11 Genes (165 Rows)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("SECTION 6: SENSITIVITY ANALYSES (b)-(f) ACROSS ALL 11 GENES (165 Rows)")
    print("=" * 100)
    bf_analyses = [
        "sens_b_no_LV_UF3", 
        "sens_c_no_scDblFinder", 
        "sens_d_ge700genes", 
        "sens_e_ge50cells", 
        "sens_f_top2000_hvg_tmm"
    ]
    df_bf = df_sens[df_sens["Analysis"].isin(bf_analyses)].copy()
    df_bf = df_bf.sort_values(["Analysis", "Contrast", "Gene"]).reset_index(drop=True)
    print(f"Sensitivities (b)-(f) Subset Rows: {len(df_bf)}")
    cols_bf = ["Analysis", "Contrast", "Gene", "log2FC", "SE", "CI_95_low", "CI_95_high", "t_stat", "PValue", "QL_Denom_df", "Donors_LV_UF", "Donors_LV_NOT_UF", "Donors_SV", "Evidence_Status"]
    print(df_bf[cols_bf].to_string(index=False))

if __name__ == "__main__":
    main()
