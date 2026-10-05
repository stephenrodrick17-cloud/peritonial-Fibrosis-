"""
Script 51: Stage 4 Tables, Ambient Recomputation, Strict Assertions, and Manifest Generation
Executes Step 2 (tables + manifest), Step 3 (all assertions), and Step 4 (readback printing).
Strictly reads all printed tables back from disk CSVs with path, SHA256, mtime, shape.
"""

import os
import sys
import json
import hashlib
import datetime
import pandas as pd
import numpy as np
from scipy import stats
import anndata as ad

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def get_file_info(filepath):
    abs_path = os.path.abspath(filepath)
    stat = os.stat(abs_path)
    mtime_utc = datetime.datetime.fromtimestamp(stat.st_mtime, tz=datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
    sha256 = compute_sha256(abs_path)
    return abs_path, mtime_utc, sha256

def main():
    print("=" * 100)
    print("SCRIPT 51: STAGE 4 CONSOLIDATION, AMBIENT PROFILING, ASSERTIONS & MANIFEST")
    print("=" * 100)
    start_time = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
    print(f"Timestamp: {start_time}")
    print(f"Working Directory: {os.getcwd()}")
    print("=" * 100)

    # -------------------------------------------------------------------------
    # 1. AMBIENT RNA RECOMPUTATION (True Raw Counts / True Total Transcriptome)
    # -------------------------------------------------------------------------
    print("\n--- 1. Recomputing Ambient RNA Profiling from Raw Matrix ---")
    adata_clean = ad.read_h5ad("data/processed/GSE248762_stage4_unsealed_hubs_clean.h5ad")
    # Total UMI per cell across all 37,487 genes is stored in obs['n_counts']
    target_cts = [
        "stromal / mesothelial-lineage (unresolved)",
        "Monocyte / macrophage",
        "cDC",
        "T cell",
        "NK cell"
    ]
    hub_genes = ["COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", "EDIL3", "LOX", "INHBA", "ISM1", "COMP"]

    ambient_rows = []
    # Use raw_counts layer from unsealed anndata for true transcriptome-referenced ambient CP10k
    raw_dense = adata_clean.layers["raw_counts"].toarray() if hasattr(adata_clean.layers["raw_counts"], "toarray") else adata_clean.layers["raw_counts"]
    gene_map = {g: i for i, g in enumerate(adata_clean.var_names)}
    n_counts = adata_clean.obs["n_counts"].values.astype(float)

    for g in hub_genes:
        g_idx = gene_map[g]
        g_counts = raw_dense[:, g_idx].astype(float)
        cp10k = (g_counts / n_counts) * 10000.0
        log1p_cp10k = np.log1p(cp10k)
        
        # Track lineage with highest expression
        ct_means = {}
        for ct in target_cts:
            ct_mask = (adata_clean.obs["cell_type"] == ct).values
            ct_means[ct] = np.mean(cp10k[ct_mask])
        
        highest_ct = max(ct_means, key=ct_means.get)
        highest_val = ct_means[highest_ct]

        for ct in target_cts:
            ct_mask = (adata_clean.obs["cell_type"] == ct).values
            tot_cells = int(ct_mask.sum())
            expr_mask = (g_counts[ct_mask] > 0)
            n_expr = int(expr_mask.sum())
            pct_expr = (n_expr / tot_cells) * 100.0 if tot_cells > 0 else 0.0
            mean_cp = float(np.mean(cp10k[ct_mask]))
            mean_log = float(np.mean(log1p_cp10k[ct_mask]))
            mean_in_expr = float(np.mean(cp10k[ct_mask][expr_mask])) if n_expr > 0 else 0.0

            ambient_rows.append({
                "Gene": g,
                "Cell_Type": ct,
                "Total_Cells": tot_cells,
                "Expressing_Cells": n_expr,
                "Percent_Expressing_%": round(pct_expr, 2),
                "Mean_CP10k": round(mean_cp, 4),
                "Mean_log1p_CP10k": round(mean_log, 4),
                "Mean_in_Expressing_CP10k": round(mean_in_expr, 4),
                "Lineage_With_Highest_Expression": f"{highest_ct} ({highest_val:.4f} CP10k)"
            })

    df_amb = pd.DataFrame(ambient_rows)
    amb_csv = "results/tables/stage4_ambient_recomputed_true_cp10k.csv"
    df_amb.to_csv(amb_csv, index=False)
    print(f"Saved true ambient table to: {amb_csv} ({len(df_amb)} rows)")

    # -------------------------------------------------------------------------
    # 2. BH FAMILY TABLE DERIVATION (m = 96 Tests)
    # -------------------------------------------------------------------------
    print("\n--- 2. Deriving Sorted 96-Test Benjamini-Hochberg Family Table ---")
    prim_csv = "results/tables/stage4_primary_edger_pseudobulk.csv"
    df_prim = pd.read_csv(prim_csv)

    pass_ev = df_prim[(df_prim["Evidence_Status"] == "PASS") & df_prim["PValue"].notna()].copy()
    m = len(pass_ev)
    assert m == 96, f"Expected 96 tests, found {m}"

    pass_ev = pass_ev.sort_values("PValue").reset_index(drop=True)
    pass_ev["rank"] = np.arange(1, m + 1)
    pass_ev["raw_BH_ratio"] = pass_ev["PValue"] * m / pass_ev["rank"]

    raw_ratios = pass_ev["raw_BH_ratio"].values
    monotonic_fdr = np.minimum.accumulate(raw_ratios[::-1])[::-1]
    monotonic_fdr = np.minimum(monotonic_fdr, 1.0)
    pass_ev["computed_BH_FDR"] = monotonic_fdr

    # Check match with edgeR BH_FDR
    diff = np.abs(pass_ev["computed_BH_FDR"].values - pass_ev["BH_FDR"].values)
    max_fdr_diff = float(np.max(diff))
    print(f"Maximum discrepancy between manual step-by-step BH and edgeR BH_FDR: {max_fdr_diff:.2e}")
    assert max_fdr_diff < 1e-6, "FDR reproduction assertion failed!"

    # Independent validation using statsmodels.stats.multitest.multipletests
    from statsmodels.stats.multitest import multipletests
    _, sm_fdr, _, _ = multipletests(pass_ev["PValue"].values, alpha=0.05, method="fdr_bh")
    pass_ev["statsmodels_BH_FDR"] = sm_fdr
    sm_diff = np.abs(pass_ev["computed_BH_FDR"].values - pass_ev["statsmodels_BH_FDR"].values)
    max_sm_diff = float(np.max(sm_diff))
    print(f"Maximum discrepancy between manual BH and statsmodels multipletests: {max_sm_diff:.2e}")
    assert max_sm_diff < 1e-12, "Statsmodels FDR reproduction assertion failed!"

    bh_csv = "results/tables/stage4_primary_96_bh_family_sorted.csv"
    pass_ev.to_csv(bh_csv, index=False)
    print(f"Saved sorted 96-test family table to: {bh_csv} ({len(pass_ev)} rows)")

    # -------------------------------------------------------------------------
    # 3. QC CEILING SENSITIVITY COMPARISON TABLE (Primary vs (i) vs (j))
    # -------------------------------------------------------------------------
    print("\n--- 3. Constructing Side-by-Side QC Ceiling Comparison Table ---")
    sens_csv = "results/tables/stage4_sensitivity_edger_pseudobulk.csv"
    df_sens = pd.read_csv(sens_csv)

    eval_genes = ["THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", "EDIL3", "LOX", "INHBA", "COMP"]
    st_prim_uf = df_prim[(df_prim["Cell_Type"].str.contains("stromal")) & 
                         (df_prim["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
                         (df_prim["Gene"].isin(eval_genes))].set_index("Gene")

    st_i_uf = df_sens[(df_sens["Analysis"] == "sens_i_no_upper_ceilings") & 
                      (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
                      (df_sens["Gene"].isin(eval_genes))].set_index("Gene")

    st_j_uf = df_sens[(df_sens["Analysis"] == "sens_j_pooled_fixed_ceiling") & 
                      (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
                      (df_sens["Gene"].isin(eval_genes))].set_index("Gene")

    comp_rows = []
    for g in eval_genes:
        lfc_p = float(st_prim_uf.loc[g, "log2FC"])
        lfc_i = float(st_i_uf.loc[g, "log2FC"])
        lfc_j = float(st_j_uf.loc[g, "log2FC"])

        pct_i = ((lfc_i - lfc_p) / abs(lfc_p)) * 100.0 if abs(lfc_p) > 0 else 0.0
        pct_j = ((lfc_j - lfc_p) / abs(lfc_p)) * 100.0 if abs(lfc_p) > 0 else 0.0

        # Pre-set rule: shrinks >50% (|lfc_sens| < 0.5 * |lfc_prim|) or flips sign
        is_i_sens = (abs(lfc_i) < 0.5 * abs(lfc_p)) or (np.sign(lfc_i) != np.sign(lfc_p))
        is_j_sens = (abs(lfc_j) < 0.5 * abs(lfc_p)) or (np.sign(lfc_j) != np.sign(lfc_p))
        is_qc_sens = is_i_sens or is_j_sens

        comp_rows.append({
            "Gene": g,
            "Primary_log2FC": round(lfc_p, 4),
            "Sens_i_log2FC": round(lfc_i, 4),
            "Pct_Change_i": round(pct_i, 1),
            "Sens_j_log2FC": round(lfc_j, 4),
            "Pct_Change_j": round(pct_j, 1),
            "Classification": "QC-SENSITIVE" if is_qc_sens else "ROBUST"
        })

    df_comp = pd.DataFrame(comp_rows)
    comp_csv = "results/tables/stage4_qc_ceiling_sensitivity_comparison.csv"
    df_comp.to_csv(comp_csv, index=False)
    print(f"Saved QC ceiling comparison table to: {comp_csv} ({len(df_comp)} rows)")

    # -------------------------------------------------------------------------
    # 4. EXECUTE ALL RIGOROUS ASSERTIONS (FAIL LOUDLY)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("STEP 3: RUNNING STRICT PROTOCOL ASSERTIONS")
    print("=" * 100)

    # Assertion 1: Unique (Analysis, Cell_Type, Contrast, Gene) keys
    print("Assertion 1: Testing key uniqueness...")
    assert not df_prim.duplicated(subset=["Analysis", "Cell_Type", "Contrast", "Gene"]).any(), "Duplicate keys in primary table!"
    assert not df_sens.duplicated(subset=["Analysis", "Cell_Type", "Contrast", "Gene"]).any(), "Duplicate keys in sensitivity table!"
    assert not pass_ev.duplicated(subset=["Cell_Type", "Contrast", "Gene"]).any(), "Duplicate keys in BH family table!"
    print("  PASS: All keys (Analysis, Cell_Type, Contrast, Gene) are 100% unique.")

    # Assertion 2: Mathematical reproduction of P-values via Student's t
    print("Assertion 2: Testing t = log2FC/SE and P = 2*pt(-|t|, df) reproduction...")
    # Test on primary estimable rows
    est_prim = df_prim[df_prim["PValue"].notna() & df_prim["SE"].notna() & (df_prim["SE"] > 0)].copy()
    est_prim["t_calc"] = est_prim["log2FC"] / est_prim["SE"]
    est_prim["p_reproduced"] = 2.0 * stats.t.sf(np.abs(est_prim["t_calc"]), df=est_prim["QL_Denom_df"])
    p_diff_prim = np.max(np.abs(est_prim["PValue"] - est_prim["p_reproduced"]))
    print(f"  Primary model max |P_stored - P_reproduced|:     {p_diff_prim:.2e}")
    assert p_diff_prim < 1e-4, f"P-value reproduction assertion failed for primary: {p_diff_prim}"

    # Test on sensitivity estimable rows
    est_sens = df_sens[df_sens["PValue"].notna() & df_sens["SE"].notna() & (df_sens["SE"] > 0)].copy()
    est_sens["t_calc"] = est_sens["log2FC"] / est_sens["SE"]
    est_sens["p_reproduced"] = 2.0 * stats.t.sf(np.abs(est_sens["t_calc"]), df=est_sens["QL_Denom_df"])
    p_diff_sens = np.max(np.abs(est_sens["PValue"] - est_sens["p_reproduced"]))
    print(f"  Sensitivity model max |P_stored - P_reproduced|: {p_diff_sens:.2e}")
    assert p_diff_sens < 1e-4, f"P-value reproduction assertion failed for sensitivities: {p_diff_sens}"
    print("  PASS: Mathematical identity t = log2FC/SE and P = 2*pt(-|t|, df) confirmed across all rows.")

    # Assertion 3: Donors per group follow definitions
    print("Assertion 3: Testing donor counts per group...")
    # Primary: 4 LV_UF, 6 LV_NOT_UF, 6 SV
    assert (df_prim["Donors_LV_UF"] == 4).all() and (df_prim["Donors_LV_NOT_UF"] == 6).all() and (df_prim["Donors_SV"] == 6).all()
    # (b): 3 LV_UF, 6 LV_NOT_UF, 6 SV
    b_rows = df_sens[df_sens["Analysis"] == "sens_b_no_LV_UF3"]
    assert (b_rows["Donors_LV_UF"] == 3).all() and (b_rows["Donors_LV_NOT_UF"] == 6).all() and (b_rows["Donors_SV"] == 6).all()
    # (e): 3 LV_UF, 2 LV_NOT_UF, 5 SV
    e_rows = df_sens[df_sens["Analysis"] == "sens_e_ge50cells"]
    assert (e_rows["Donors_LV_UF"] == 3).all() and (e_rows["Donors_LV_NOT_UF"] == 2).all() and (e_rows["Donors_SV"] == 5).all()
    # (h): 1 LV_UF, 2 LV_NOT_UF, 5 SV -> NOT ESTIMABLE
    h_rows = df_sens[df_sens["Analysis"] == "sens_h_ge50cells_med1000genes"]
    assert (h_rows["Donors_LV_UF"] == 1).all() and (h_rows["Evidence_Status"] == "NOT ESTIMABLE").all()
    assert h_rows["PValue"].isna().all()
    print("  PASS: Donors per group follow exact protocol definitions.")

    # Assertion 4: LODO "drop LV_UF-3" equals (b)
    print("Assertion 4: Testing LODO 'drop LV_UF-3' == (b)...")
    lodo_uf3 = df_sens[df_sens["Analysis"] == "lodo_out_LV_UF-3"].set_index(["Contrast", "Gene"])
    b_df = df_sens[df_sens["Analysis"] == "sens_b_no_LV_UF3"].set_index(["Contrast", "Gene"])
    lfc_diff = np.max(np.abs(lodo_uf3["log2FC"] - b_df["log2FC"]))
    se_diff = np.max(np.abs(lodo_uf3["SE"] - b_df["SE"]))
    pval_diff = np.max(np.abs(lodo_uf3["PValue"] - b_df["PValue"]))
    print(f"  Max difference between LODO drop LV_UF-3 and (b): log2FC diff={lfc_diff:.2e}, SE diff={se_diff:.2e}, P diff={pval_diff:.2e}")
    assert lfc_diff < 1e-10 and se_diff < 1e-10 and pval_diff < 1e-10, "LODO drop LV_UF-3 does not equal (b)!"
    print("  PASS: LODO 'drop LV_UF-3' equals sensitivity (b) identically across all columns.")

    # Assertion 5: No two LODO rows identical unless dropped libraries are identical
    print("Assertion 5: Testing LODO distinctness across dropped donors...")
    lodo_all = df_sens[df_sens["Analysis"].str.startswith("lodo_out_")].copy()
    donors = lodo_all["Analysis"].unique()
    assert len(donors) == 16, f"Expected 16 LODO donors, found {len(donors)}"
    for i in range(len(donors)):
        for j in range(i + 1, len(donors)):
            d1, d2 = donors[i], donors[j]
            s1 = lodo_all[lodo_all["Analysis"] == d1]["log2FC"].values
            s2 = lodo_all[lodo_all["Analysis"] == d2]["log2FC"].values
            assert not np.allclose(s1, s2, atol=1e-5), f"LODO runs for {d1} and {d2} are unexpectedly identical!"
    print("  PASS: No two distinct LODO runs produce identical log2FC profiles.")

    # Assertion 6: Ambient: percent expressing == expressing/total
    print("Assertion 6: Testing percent expressing == expressing/total...")
    pct_calc = (df_amb["Expressing_Cells"] / df_amb["Total_Cells"]) * 100.0
    pct_diff = np.max(np.abs(df_amb["Percent_Expressing_%"] - np.round(pct_calc, 2)))
    assert pct_diff < 1e-4, "Percent expressing assertion failed!"
    print("  PASS: Percent expressing equals expressing/total*100 across all cell types.")

    # Assertion 7: Per-cell-type totals equal obs-table sums
    print("Assertion 7: Testing per-cell-type totals against obs-table sums...")
    adata_obs_full = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad").obs
    n_mono_full = int((adata_obs_full["cell_type"] == "Monocyte / macrophage").sum())
    n_mono_no12 = int(((adata_obs_full["cell_type"] == "Monocyte / macrophage") & (adata_obs_full["leiden"] != "12")).sum())
    n_stromal = int((adata_obs_full["cell_type"] == "stromal / mesothelial-lineage (unresolved)").sum())

    print(f"  Obs table Monocyte/macrophage (incl. cluster 12): {n_mono_full} (expected: 22,925)")
    print(f"  Obs table Monocyte/macrophage (without cluster 12): {n_mono_no12} (expected: 22,655)")
    print(f"  Obs table Stromal cells:                            {n_stromal} (expected: 3,670)")
    assert n_mono_full == 22925, f"Expected 22,925 mono incl. 12, got {n_mono_full}"
    assert n_mono_no12 == 22655, f"Expected 22,655 mono without 12, got {n_mono_no12}"
    assert n_stromal == 3670, f"Expected 3,670 stromal, got {n_stromal}"
    print("  PASS: Per-cell-type totals equal obs-table sums exactly.")

    # -------------------------------------------------------------------------
    # 5. WRITE MANIFEST.JSON (LAST STEP BEFORE PRINTING)
    # -------------------------------------------------------------------------
    print("\n--- 5. Generating Cryptographic manifest.json as Last Pipeline Step ---")
    manifest_files = [
        "results/tables/stage4_primary_edger_pseudobulk.csv",
        "results/tables/stage4_primary_96_bh_family_sorted.csv",
        "results/tables/stage4_sensitivity_edger_pseudobulk.csv",
        "results/tables/stage4_qc_ceiling_sensitivity_comparison.csv",
        "results/tables/stage4_ambient_recomputed_true_cp10k.csv"
    ]

    manifest_entries = {}
    for fpath in manifest_files:
        p, mt, h = get_file_info(fpath)
        df_tmp = pd.read_csv(fpath)
        manifest_entries[os.path.basename(fpath)] = {
            "path": p,
            "mtime_utc": mt,
            "sha256": h,
            "n_rows": len(df_tmp),
            "shape": list(df_tmp.shape)
        }

    manifest_path = "results/tables/stage4_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest_entries, f, indent=2)

    p_man, mt_man, h_man = get_file_info(manifest_path)
    print(f"Manifest written successfully to: {p_man}")
    print(f"Manifest mtime: {mt_man}")
    print(f"Manifest SHA256: {h_man}")
    print(json.dumps(manifest_entries, indent=2))

    # -------------------------------------------------------------------------
    # 6. READ BACK AND PRINT VERBATIM FROM MANIFEST FILES
    # -------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("STEP 4: VERBATIM READBACK FROM MANIFEST CSV FILES")
    print("=" * 100)

    # Table 1: Primary Table (5 Cell Types x 3 Contrasts x 11 Genes)
    print("\n" + "#" * 90)
    print("TABLE 1: PRIMARY PSEUDOBULK edgeR RESULTS (5 Cell Types x 3 Contrasts x 11 Genes)")
    print("#" * 90)
    info1 = manifest_entries["stage4_primary_edger_pseudobulk.csv"]
    print(f"Path:   {info1['path']}\nmtime:  {info1['mtime_utc']}\nSHA256: {info1['sha256']}\nShape:  {info1['shape']}\n")
    df1 = pd.read_csv(info1["path"])
    print(df1.to_string(index=False))

    # Table 2: BH Family Table Sorted by P-value (m = 96 Tests)
    print("\n" + "#" * 90)
    print("TABLE 2: BENJAMINI-HOCHBERG MULTIPLE TESTING FAMILY (m = 96 TESTS)")
    print("#" * 90)
    info2 = manifest_entries["stage4_primary_96_bh_family_sorted.csv"]
    print(f"Path:   {info2['path']}\nmtime:  {info2['mtime_utc']}\nSHA256: {info2['sha256']}\nShape:  {info2['shape']}\n")
    df2 = pd.read_csv(info2["path"])
    print(df2.to_string(index=False))

    # Table 3: Stromal Sensitivity Analyses (b)-(f), (i), (j)
    print("\n" + "#" * 90)
    print("TABLE 3: STROMAL SENSITIVITY ANALYSES (b)-(f), (i), (j)")
    print("#" * 90)
    info3 = manifest_entries["stage4_sensitivity_edger_pseudobulk.csv"]
    print(f"Path:   {info3['path']}\nmtime:  {info3['mtime_utc']}\nSHA256: {info3['sha256']}\nShape:  {info3['shape']}\n")
    df3 = pd.read_csv(info3["path"])
    sub_bfij = df3[df3["Analysis"].isin([
        "sens_b_no_LV_UF3", "sens_c_no_scDblFinder", "sens_d_ge700genes",
        "sens_e_ge50cells", "sens_f_top2000_hvg_tmm", "sens_h_ge50cells_med1000genes",
        "sens_i_no_upper_ceilings", "sens_j_pooled_fixed_ceiling"
    ])]
    print(sub_bfij.to_string(index=False))

    # Table 4: Stromal Leave-One-Donor-Out (LODO) Runs
    print("\n" + "#" * 90)
    print("TABLE 4: STROMAL LEAVE-ONE-DONOR-OUT (LODO) RUNS (16 DONORS x 3 CONTRASTS x 11 GENES)")
    print("#" * 90)
    sub_lodo = df3[df3["Analysis"].str.startswith("lodo_out_")]
    print(f"Total LODO rows: {len(sub_lodo)}")
    print(sub_lodo.to_string(index=False))

    # Table 5: QC Ceiling Comparison Table Side-by-Side
    print("\n" + "#" * 90)
    print("TABLE 5: QC CEILING COMPARISON TABLE (Primary vs (i) vs (j))")
    print("#" * 90)
    info5 = manifest_entries["stage4_qc_ceiling_sensitivity_comparison.csv"]
    print(f"Path:   {info5['path']}\nmtime:  {info5['mtime_utc']}\nSHA256: {info5['sha256']}\nShape:  {info5['shape']}\n")
    df5 = pd.read_csv(info5["path"])
    print(df5.to_string(index=False))

    # Table 6: Ambient RNA Profiling Table
    print("\n" + "#" * 90)
    print("TABLE 6: AMBIENT RNA PROFILING (Recomputed True CP10k across 5 Lineages)")
    print("#" * 90)
    info6 = manifest_entries["stage4_ambient_recomputed_true_cp10k.csv"]
    print(f"Path:   {info6['path']}\nmtime:  {info6['mtime_utc']}\nSHA256: {info6['sha256']}\nShape:  {info6['shape']}\n")
    df6 = pd.read_csv(info6["path"])
    print(df6.to_string(index=False))

    print("\n" + "=" * 100)
    print("SCRIPT 51 COMPLETED ALL ASSERTIONS AND READBACK PRINTS SUCCESSFULLY.")
    print("=" * 100)

if __name__ == "__main__":
    main()
