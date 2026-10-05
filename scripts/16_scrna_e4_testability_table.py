"""
Step E4: Minimum Evidence Testability Matrix across Cell Types and Donors
Strictly Hub-Gene-Blind execution with code-level guard.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import anndata as ad

sys.path.append("scripts")
from guard import guard_check

print("=================================================================")
print("STEP E4: MINIMUM EVIDENCE TESTABILITY MATRIX (HUB-BLIND)")
print("=================================================================")

# Guard check: Ensure no hub genes are accessed
guard_check([], stage="E4_start")

# Load annotated AnnData
adata_path = "data/processed/GSE248762_annotated.h5ad"
adata = ad.read_h5ad(adata_path)
print(f"Loaded annotated cohort: {adata.n_obs:,} cells across {adata.obs['donor_id'].nunique()} donors.")

# 1. Donor x Cell Type contingency table
ct_donor = pd.crosstab(adata.obs["donor_id"], adata.obs["cell_type"])
# Add group column
donor_to_group = adata.obs.groupby("donor_id")["group"].first().to_dict()
ct_donor["Group"] = ct_donor.index.map(donor_to_group)
ct_donor = ct_donor.reset_index()

ct_donor.to_csv("results/tables/E4_donor_cell_counts_by_type.csv", index=False)
print("Saved results/tables/E4_donor_cell_counts_by_type.csv")

# 2. Evaluate pre-specified testability rule for planned contrasts:
# Minimum Evidence Rule: >= 20 cells in >= 3 donors in target group, and >= 2 donors in comparator group.

contrasts = [
    {
        "name": "LV_UF_vs_LV_NOT_UF",
        "target_group": "LV_UF",
        "comp_group": "LV_NOT_UF",
        "min_target_donors": 3,
        "min_comp_donors": 2,
        "description": "Ultrafiltration failure effect in long vintage effluent"
    },
    {
        "name": "LV_UF_vs_SV",
        "target_group": "LV_UF",
        "comp_group": "SV",
        "min_target_donors": 3,
        "min_comp_donors": 2,
        "description": "UFF vs Short Vintage effluent"
    },
    {
        "name": "LV_NOT_UF_vs_SV",
        "target_group": "LV_NOT_UF",
        "comp_group": "SV",
        "min_target_donors": 3,
        "min_comp_donors": 2,
        "description": "PD duration vintage effect without UFF"
    }
]

testability_rows = []

all_cell_types = [c for c in ct_donor.columns if c not in ["donor_id", "Group"]]

for c_info in contrasts:
    c_name = c_info["name"]
    tgt_grp = c_info["target_group"]
    cmp_grp = c_info["comp_group"]
    
    tgt_donors = ct_donor[ct_donor["Group"] == tgt_grp]
    cmp_donors = ct_donor[ct_donor["Group"] == cmp_grp]
    
    n_tgt_total = len(tgt_donors)
    n_cmp_total = len(cmp_donors)
    
    for ct in all_cell_types:
        tgt_ge20 = int((tgt_donors[ct] >= 20).sum())
        cmp_ge20 = int((cmp_donors[ct] >= 20).sum())
        
        tgt_total_cells = int(tgt_donors[ct].sum())
        cmp_total_cells = int(cmp_donors[ct].sum())
        
        is_testable = (tgt_ge20 >= c_info["min_target_donors"]) and (cmp_ge20 >= c_info["min_comp_donors"])
        
        if is_testable:
            status = "TESTABLE (Sufficient cell & donor coverage)"
        else:
            reasons = []
            if tgt_ge20 < c_info["min_target_donors"]:
                reasons.append(f"Target {tgt_grp} has only {tgt_ge20}/{n_tgt_total} donors with >=20 cells (need {c_info['min_target_donors']})")
            if cmp_ge20 < c_info["min_comp_donors"]:
                reasons.append(f"Comparator {cmp_grp} has only {cmp_ge20}/{n_cmp_total} donors with >=20 cells (need {c_info['min_comp_donors']})")
            status = f"UNDERPOWERED ({'; '.join(reasons)})"
            
        testability_rows.append({
            "Contrast": c_name,
            "Description": c_info["description"],
            "Cell_Type": ct,
            "Target_Group": tgt_grp,
            "Target_Donors_ge20": f"{tgt_ge20}/{n_tgt_total}",
            "Target_Total_Cells": tgt_total_cells,
            "Comparator_Group": cmp_grp,
            "Comparator_Donors_ge20": f"{cmp_ge20}/{n_cmp_total}",
            "Comparator_Total_Cells": cmp_total_cells,
            "Testability_Status": status,
            "Is_Testable": is_testable
        })

df_testability = pd.DataFrame(testability_rows)
df_testability.to_csv("results/tables/E4_celltype_testability_table.csv", index=False)
print("\n=================================================================")
print("Saved results/tables/E4_celltype_testability_table.csv")
print(df_testability[["Contrast", "Cell_Type", "Target_Donors_ge20", "Comparator_Donors_ge20", "Testability_Status"]])
print("=================================================================")
print("\nStep E4 finished successfully.")
