"""
Script 52: Print Manifest Tables Verbatim
Reads each file from results/tables/stage4_manifest.json and prints df.to_string()
Strictly zero typed markdown tables.
"""

import json
import pandas as pd

def print_manifest_entry(manifest, key, filter_query=None, title=None):
    entry = manifest[key]
    print("=" * 100)
    print(f"FILE:   {key}")
    if title:
        print(f"TITLE:  {title}")
    print(f"PATH:   {entry['path']}")
    print(f"MTIME:  {entry['mtime_utc']}")
    print(f"SHA256: {entry['sha256']}")
    print(f"SHAPE:  {entry['shape']}")
    print("=" * 100)
    df = pd.read_csv(entry["path"])
    if filter_query is not None:
        df = df.query(filter_query)
        print(f"SUBSET: {filter_query} (Rows: {len(df)})")
    print(df.to_string(index=False))
    print()

def main():
    with open("results/tables/stage4_manifest.json", "r") as f:
        manifest = json.load(f)

    # 1. Primary table (5 cell types x 3 contrasts x 11 genes = 165 rows)
    print_manifest_entry(
        manifest, 
        "stage4_primary_edger_pseudobulk.csv", 
        title="PRIMARY PSEUDOBULK edgeR RESULTS (5 Cell Types x 3 Contrasts x 11 Genes = 165 Rows)"
    )

    # 2. BH Family (m = 96 tests)
    print_manifest_entry(
        manifest, 
        "stage4_primary_96_bh_family_sorted.csv", 
        title="BENJAMINI-HOCHBERG MULTIPLE TESTING FAMILY (m = 96 Tests Sorted by P-value)"
    )

    # 3. Sensitivities (b)-(f), (h), (i), (j)
    sens_filter = "Analysis in ['sens_b_no_LV_UF3', 'sens_c_no_scDblFinder', 'sens_d_ge700genes', 'sens_e_ge50cells', 'sens_f_top2000_hvg_tmm', 'sens_h_ge50cells_med1000genes', 'sens_i_no_upper_ceilings', 'sens_j_pooled_fixed_ceiling']"
    print_manifest_entry(
        manifest, 
        "stage4_sensitivity_edger_pseudobulk.csv", 
        filter_query=sens_filter,
        title="STROMAL SENSITIVITY ANALYSES (b)-(f), (h), (i), (j)"
    )

    # 4. Stromal LODO runs (16 donors x 3 contrasts x 11 genes = 528 rows)
    lodo_filter = "Analysis.str.startswith('lodo_out_')"
    print_manifest_entry(
        manifest, 
        "stage4_sensitivity_edger_pseudobulk.csv", 
        filter_query=lodo_filter,
        title="STROMAL LEAVE-ONE-DONOR-OUT (LODO) RUNS (16 Donors x 3 Contrasts x 11 Genes = 528 Rows)"
    )

    # 5. QC Ceiling Comparison
    print_manifest_entry(
        manifest,
        "stage4_qc_ceiling_sensitivity_comparison.csv",
        title="QC CEILING SENSITIVITY COMPARISON (Primary vs (i) vs (j) with Pre-set Rule)"
    )

    # 6. Ambient table
    print_manifest_entry(
        manifest, 
        "stage4_ambient_recomputed_true_cp10k.csv", 
        title="AMBIENT RNA PROFILING RECOMPUTATION (True Raw Counts across 5 Lineages = 55 Rows)"
    )

if __name__ == "__main__":
    main()
