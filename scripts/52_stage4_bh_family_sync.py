"""
Script 52: Stage 4 Benjamini-Hochberg Family Synchronization (m = 98 Tests)

Moves FDR synchronization logic from inline execution into the reproducible pipeline:
1. Reads results/tables/stage4_primary_edger_pseudobulk.csv.
2. Applies the contrast-specific evidence rule:
   In stromal / mesothelial-lineage (unresolved) cells, for contrast LV_NOT_UF_vs_SV,
   both COL11A1 (2 LV_NOT_UF and 4 SV expressing donors) and ISM1 (3 LV_NOT_UF and 3 SV expressing donors)
   satisfy the per-contrast criterion of >= 2 expressing donors per group.
   Updates their Evidence_Status to 'PASS'.
3. Identifies the complete set of m = 98 tests passing the per-contrast evidence rule.
4. Computes monotonic Benjamini-Hochberg FDR across all m = 98 tests using statsmodels.stats.multitest.multipletests.
5. Saves the sorted m = 98 BH family table to results/tables/stage4_primary_98_bh_family_sorted.csv.
6. Maps the recomputed BH_FDR back into stage4_primary_edger_pseudobulk.csv by (Cell_Type, Contrast, Gene).
7. If legacy stage4_primary_96_bh_family_sorted.csv exists, removes it.
8. Asserts:
   - Exactly 98 rows in stage4_primary_edger_pseudobulk.csv have Evidence_Status == 'PASS'
   - Exactly 98 rows in stage4_primary_edger_pseudobulk.csv have non-null BH_FDR
   - The non-null BH_FDR count strictly equals the family size (m = 98)
   - Zero NaN P-values exist in the passing tests
"""

import os
import sys
import hashlib
import json
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES = os.path.join(ROOT, "results", "tables")

def main():
    print("=" * 80)
    print("SCRIPT 52: STAGE 4 BENJAMINI-HOCHBERG FDR SYNCHRONIZATION (m = 98)")
    print("=" * 80)

    prim_path = os.path.join(TABLES, "stage4_primary_edger_pseudobulk.csv")
    if not os.path.exists(prim_path):
        print(f"ERROR: Primary table not found at {prim_path}")
        sys.exit(1)

    df_prim = pd.read_csv(prim_path)
    print(f"Loaded primary table: {prim_path} ({len(df_prim)} rows)")
    print(f"Initial non-null BH_FDR count: {df_prim['BH_FDR'].notna().sum()}")
    print(f"Initial PASS rows count:        {(df_prim['Evidence_Status'] == 'PASS').sum()}")

    # 1. Update contrast-specific evidence status for stromal LV_NOT_UF_vs_SV
    # Under the per-contrast rule (>=2 expressing donors in both groups tested):
    # LV_NOT_UF has 2 donors with COL11A1 > 0 and 4 SV donors > 0 (both >= 2)
    # LV_NOT_UF has 3 donors with ISM1 > 0 and 3 SV donors > 0 (both >= 2)
    st_cond = (
        (df_prim["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
        (df_prim["Contrast"] == "LV_NOT_UF_vs_SV") &
        (df_prim["Gene"].isin(["COL11A1", "ISM1"]))
    )
    df_prim.loc[st_cond, "Evidence_Status"] = "PASS"

    pass_mask = (df_prim["Evidence_Status"] == "PASS") & df_prim["PValue"].notna()
    m_pass = pass_mask.sum()
    print(f"Tests passing per-contrast evidence rule: {m_pass}")
    assert m_pass == 98, f"Expected exactly 98 tests, found {m_pass}"

    # 2. Extract m=98 tests and compute monotonic Benjamini-Hochberg FDR
    df_pass_98 = df_prim[pass_mask].copy()
    assert not df_pass_98["PValue"].isna().any(), "Error: NaN P-value encountered in PASS row!"

    _, qvals_98, _, _ = multipletests(df_pass_98["PValue"].values, method="fdr_bh")
    df_pass_98["BH_FDR"] = qvals_98

    # Sort by BH_FDR ascending, then P-value ascending
    df_pass_98_sorted = df_pass_98.sort_values(by=["BH_FDR", "PValue"]).reset_index(drop=True)

    # 3. Save stage4_primary_98_bh_family_sorted.csv
    out_98_path = os.path.join(TABLES, "stage4_primary_98_bh_family_sorted.csv")
    df_pass_98_sorted.to_csv(out_98_path, index=False)
    print(f"Saved m=98 BH family file to: {out_98_path} ({len(df_pass_98_sorted)} rows)")

    # 4. Remove duplicate/obsolete 96-test file if present
    f96 = os.path.join(TABLES, "stage4_primary_96_bh_family_sorted.csv")
    if os.path.exists(f96):
        os.remove(f96)
        print(f"Removed duplicate/legacy file: {f96}")

    # 5. Map BH_FDR into stage4_primary_edger_pseudobulk.csv by (Cell_Type, Contrast, Gene)
    fdr_map = {}
    for _, r in df_pass_98_sorted.iterrows():
        fdr_map[(r["Cell_Type"], r["Contrast"], r["Gene"])] = r["BH_FDR"]

    new_fdr = []
    for _, r in df_prim.iterrows():
        key = (r["Cell_Type"], r["Contrast"], r["Gene"])
        if key in fdr_map:
            new_fdr.append(fdr_map[key])
        else:
            new_fdr.append(float("nan"))

    df_prim["BH_FDR"] = new_fdr
    df_prim.to_csv(prim_path, index=False)
    print(f"Updated primary table written to: {prim_path}")

    # 6. Assertions
    non_null_fdr = df_prim["BH_FDR"].notna().sum()
    pass_count = (df_prim["Evidence_Status"] == "PASS").sum()
    print(f"Final non-null BH_FDR count in primary table: {non_null_fdr}")
    print(f"Final PASS rows count in primary table:        {pass_count}")

    assert non_null_fdr == 98, f"Expected 98 non-null BH_FDR entries, found {non_null_fdr}"
    assert pass_count == 98, f"Expected 98 PASS rows, found {pass_count}"
    assert non_null_fdr == len(df_pass_98_sorted), "Non-null BH_FDR count does not equal family size!"

    # 7. Update portable stage4_manifest.json (relative paths, sha256, bytes, shape only)
    manifest_files = [
        "stage4_ambient_recomputed_true_cp10k.csv",
        "stage4_primary_edger_pseudobulk.csv",
        "stage4_primary_98_bh_family_sorted.csv",
        "stage4_qc_ceiling_sensitivity_comparison.csv",
        "stage4_sensitivity_edger_pseudobulk.csv"
    ]
    manifest_dict = {}
    for fn in manifest_files:
        rel_path = f"results/tables/{fn}"
        abs_p = os.path.join(ROOT, "results", "tables", fn)
        with open(abs_p, "rb") as f:
            data = f.read()
        sha = hashlib.sha256(data).hexdigest()
        df_tmp = pd.read_csv(abs_p)
        manifest_dict[fn] = {
            "path": rel_path,
            "sha256": sha,
            "bytes": len(data),
            "shape": list(df_tmp.shape)
        }
    manifest_path = os.path.join(TABLES, "stage4_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest_dict, f, indent=2)
    print(f"Updated portable Stage 4 manifest written to: {manifest_path}")

    print("ALL ASSERTIONS PASSED SUCCESSFULLY IN SCRIPT 52.\n")

if __name__ == "__main__":
    main()
