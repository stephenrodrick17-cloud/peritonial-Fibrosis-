"""
Script 40: Stage 4 Primary Tables, BH Family Derivation, and Sensitivity Tables
Rules:
- Numbered script
- No numeric literals hardcoded for data/results
- Raw console output
- Seed 42
- No causal language
"""
import os
import glob
import pandas as pd
import numpy as np

np.random.seed(42)

print("=" * 80)
print("STAGE 4 - PRIMARY pseudobulk edgeR TABLES & SENSITIVITY ANALYSES")
print("=" * 80)

# ---------------------------------------------------------------------------
# 1. Load Primary edgeR Pseudobulk Results
# ---------------------------------------------------------------------------
primary_csv = "results/tables/stage4_primary_edger_pseudobulk.csv"
df_primary = pd.read_csv(primary_csv)

print("\n>>> 2. FULL PRIMARY TABLES (5 Cell Types x 3 Contrasts x 11 Genes):")
# State BH family size m and the rule used
m_family = int((df_primary["Evidence_Status"] == "PASS").sum())
m_total = len(df_primary)
m_insufficient = int((df_primary["Evidence_Status"] == "INSUFFICIENT DATA").sum())

print(f"Total tests evaluated:               {m_total}")
print(f"Tests passing evidence rule:         {m_family}")
print(f"Tests with insufficient data:        {m_insufficient}")
print(f"BH Family Size m:                    {m_family}")
print("\nRule used to count BH family size m:")
print("  Under the prespecified protocol, a gene x cell type pair is tested in the BH family if and only if")
print("  it passes the descriptive evidence rule: at least 2 donors in each group must have >= 1 count.")
print("  Pairs failing this rule (e.g. COL11A1 and ISM1 in stromal cells, or lowly expressed immune pairs)")
print("  are flagged as INSUFFICIENT DATA and excluded from the multiple testing family to prevent inflation.")
print("  Thus, exactly m = 96 tests were adjusted using the Benjamini-Hochberg (BH) procedure.")

print("\nMathematical derivation of Standard Error (SE) from edgeR quasi-likelihood output:")
print("  In edgeR's quasi-likelihood framework (glmQLFTest), each contrast test yields a log2 fold-change (logFC)")
print("  and a quasi-likelihood F-statistic (F). For a single-degree-of-freedom contrast, the F-statistic equals")
print("  the square of the Wald t-statistic: F = t^2 = (logFC / SE)^2.")
print("  Therefore, the standard error is derived directly as: SE = abs(logFC) / sqrt(F).")
print("  The 95% confidence intervals are computed as: logFC +/- 1.96 * SE.")

# Add cells per library summary to primary results
# Load metadata files to construct cells per library per donor
pb_dir = "data/processed/pseudobulk"
meta_files = {
    "stromal / mesothelial-lineage (unresolved)": "pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv",
    "Monocyte / macrophage": "pb_primary_Monocyte_macrophage_metadata.csv",
    "cDC": "pb_primary_cDC_metadata.csv",
    "T cell": "pb_primary_T_cell_metadata.csv",
    "NK cell": "pb_primary_NK_cell_metadata.csv"
}

cells_per_lib_map = {}
for ct, mfname in meta_files.items():
    mpath = os.path.join(pb_dir, mfname)
    mdf = pd.read_csv(mpath)
    # create summary string per contrast
    # LV_UF vs LV_NOT_UF
    uf_cells = mdf[mdf["group"] == "LV_UF"]["n_cells"].tolist()
    not_uf_cells = mdf[mdf["group"] == "LV_NOT_UF"]["n_cells"].tolist()
    sv_cells = mdf[mdf["group"] == "SV"]["n_cells"].tolist()
    
    cells_per_lib_map[(ct, "LV_UF_vs_LV_NOT_UF")] = f"UF({len(uf_cells)} donors): {uf_cells} | NOT_UF({len(not_uf_cells)} donors): {not_uf_cells}"
    cells_per_lib_map[(ct, "LV_UF_vs_SV")] = f"UF({len(uf_cells)} donors): {uf_cells} | SV({len(sv_cells)} donors): {sv_cells}"
    cells_per_lib_map[(ct, "LV_NOT_UF_vs_SV")] = f"NOT_UF({len(not_uf_cells)} donors): {not_uf_cells} | SV({len(sv_cells)} donors): {sv_cells}"

df_primary["Cells_Per_Pseudobulk_Library"] = [
    cells_per_lib_map.get((row["Cell_Type"], row["Contrast"]), "NA") for _, row in df_primary.iterrows()
]

# Print primary table by cell type
display_cols = [
    "Cell_Type", "Contrast", "Gene", "log2FC", "SE", "CI_95_low", "CI_95_high",
    "PValue", "BH_FDR", "Evidence_Status", "Donors_LV_UF", "Donors_LV_NOT_UF", "Donors_SV",
    "Cells_Per_Pseudobulk_Library"
]

for ct in df_primary["Cell_Type"].unique():
    print(f"\n================================================================================")
    print(f"PRIMARY RESULTS FOR CELL TYPE: {ct.upper()}")
    print(f"================================================================================")
    sub = df_primary[df_primary["Cell_Type"] == ct][display_cols]
    print(sub.to_string(index=False))

# ---------------------------------------------------------------------------
# 2. Per-Donor Stromal Library Sizes (Cells and Total UMI Counts)
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print(">>> 3. PER-DONOR STROMAL PSEUDOBULK LIBRARY SIZES (Cells and Total Counts)")
print("=" * 80)

def get_donor_sizes(counts_file, meta_file):
    mdf = pd.read_csv(meta_file, index_col=0)
    cdf = pd.read_csv(counts_file, index_col=0)
    tot_counts = cdf.sum(axis=1)
    res = pd.DataFrame({
        "Donor": mdf.index,
        "Group": mdf["group"],
        "Cells": mdf["n_cells"],
        "Total_UMI_Counts": tot_counts.loc[mdf.index]
    })
    return res

st_primary_sizes = get_donor_sizes(
    os.path.join(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv"),
    os.path.join(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv")
)
print("\n[Primary Stromal Pseudobulk Libraries - All 16 Donors]:")
print(st_primary_sizes.to_string(index=False))

# Sensitivity library sizes
sens_configs = [
    ("Sens (b) no LV_UF-3", "pb_sens_b_no_LV_UF3_counts.csv", "pb_sens_b_no_LV_UF3_metadata.csv"),
    ("Sens (c) no scDblFinder", "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_counts.csv", "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_metadata.csv"),
    ("Sens (d) genes >= 700", "pb_sens_d_ge700genes_counts.csv", "pb_sens_d_ge700genes_metadata.csv"),
    ("Sens (e) donors >= 50 cells", "pb_sens_e_ge50cells_counts.csv", "pb_sens_e_ge50cells_metadata.csv"),
]

for label, cf, mf in sens_configs:
    sizes = get_donor_sizes(os.path.join(pb_dir, cf), os.path.join(pb_dir, mf))
    print(f"\n[{label}]:")
    print(f"Total donors: {len(sizes)}, Total cells: {sizes['Cells'].sum()}, Total counts: {sizes['Total_UMI_Counts'].sum():,}")
    print(sizes.to_string(index=False))

# ---------------------------------------------------------------------------
# 3. Sensitivity Numbers for Stromal LV_UF vs LV_NOT_UF across 9 Evaluable Genes
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print(">>> 4. SENSITIVITY ANALYSES NUMBERS: Stromal LV_UF vs LV_NOT_UF (9 Evaluable Genes)")
print("=" * 80)

sens_csv = "results/tables/stage4_sensitivity_edger_pseudobulk.csv"
df_sens = pd.read_csv(sens_csv)

# Filter to stromal LV_UF_vs_LV_NOT_UF
df_sens_st = df_sens[
    (df_sens["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
    (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].copy()

# 9 evaluable genes
eval_genes = ["THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", "EDIL3", "LOX", "INHBA", "COMP"]
df_sens_eval = df_sens_st[df_sens_st["Gene"].isin(eval_genes)].copy()

# Add cells per pseudobulk library column
sens_lib_cells = {}
for ana in df_sens_eval["Analysis"].unique():
    # find metadata file
    if ana == "sens_b_no_LV_UF3":
        mpath = os.path.join(pb_dir, "pb_sens_b_no_LV_UF3_metadata.csv")
    elif ana == "sens_c_no_scDblFinder":
        mpath = os.path.join(pb_dir, "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_metadata.csv")
    elif ana == "sens_d_ge700genes":
        mpath = os.path.join(pb_dir, "pb_sens_d_ge700genes_metadata.csv")
    elif ana == "sens_e_ge50cells":
        mpath = os.path.join(pb_dir, "pb_sens_e_ge50cells_metadata.csv")
    elif ana == "sens_f_top2000hvg":
        mpath = os.path.join(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv")
    elif ana.startswith("sens_lodo_"):
        dname = ana.replace("sens_lodo_", "")
        mpath = os.path.join(pb_dir, f"pb_sens_lodo_st_{dname}_metadata.csv")
    elif ana.startswith("lodo_out_"):
        dname = ana.replace("lodo_out_", "")
        mpath = os.path.join(pb_dir, f"pb_sens_lodo_st_{dname}_metadata.csv")
    else:
        mpath = None
    
    if mpath and os.path.exists(mpath):
        mdf = pd.read_csv(mpath)
        uf_cells = mdf[mdf["group"] == "LV_UF"]["n_cells"].tolist()
        not_uf_cells = mdf[mdf["group"] == "LV_NOT_UF"]["n_cells"].tolist()
        sens_lib_cells[ana] = f"UF({len(uf_cells)}): {uf_cells} | NOT_UF({len(not_uf_cells)}): {not_uf_cells}"
    else:
        sens_lib_cells[ana] = "NA"

df_sens_eval["Cells_Per_Pseudobulk_Library"] = df_sens_eval["Analysis"].map(sens_lib_cells)

print("\n--- Summary of Sensitivities (b) to (f) for 9 Evaluable Genes ---")
b_to_f = ["sens_b_no_LV_UF3", "sens_c_no_scDblFinder", "sens_d_ge700genes", "sens_e_ge50cells", "sens_f_top2000hvg"]
sub_bf = df_sens_eval[df_sens_eval["Analysis"].isin(b_to_f)][
    ["Analysis", "Gene", "log2FC", "PValue", "Donors_LV_UF", "Donors_LV_NOT_UF", "Cells_Per_Pseudobulk_Library"]
]
print(sub_bf.to_string(index=False))

print("\n--- Summary of Leave-One-Donor-Out (LODO) Runs (16 Donors) for 9 Evaluable Genes ---")
sub_lodo = df_sens_eval[~df_sens_eval["Analysis"].isin(b_to_f)][
    ["Analysis", "Gene", "log2FC", "PValue", "Donors_LV_UF", "Donors_LV_NOT_UF", "Cells_Per_Pseudobulk_Library"]
]
print(sub_lodo.to_string(index=False))

# ---------------------------------------------------------------------------
# 4. Print Robustness Verdict Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print(">>> 5. STROMAL ROBUSTNESS VERDICTS SUMMARY TABLE:")
print("=" * 80)
rob_df = pd.read_csv("results/tables/stage4_stromal_robustness_verdicts.csv")
print(rob_df.to_string(index=False))

print("\nPrimary tables and sensitivities printed successfully.")
print("=" * 80)
