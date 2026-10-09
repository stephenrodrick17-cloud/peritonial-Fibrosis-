"""
Step E: Stage 3 scRNA-seq Audit & Summary Assembly
Verifies all outputs, checks code-level guard logs for zero hub gene queries, and summarizes findings.
"""
import os
import sys
import json
import pandas as pd

sys.path.append("scripts")
from guard import guard_check

print("=================================================================")
print("STAGE 3 SUMMARY AUDIT & CODE-LEVEL GUARD VERIFICATION")
print("=================================================================")

# 1. Inspect Guard Log
guard_log = "audit/guard_audit.log"
guard_passed = False
guard_violations = 0
guard_checks_count = 0

if os.path.exists(guard_log):
    with open(guard_log, "r") as f:
        log_lines = f.readlines()
    for line in log_lines:
        if "[GUARD VIOLATION]" in line:
            guard_violations += 1
        elif "[GUARD PASS]" in line:
            guard_checks_count += 1
            
    if guard_violations == 0 and guard_checks_count > 0:
        guard_passed = True
        print(f"GUARD INTEGRITY PASS: {guard_checks_count} guard checks executed, 0 violations.")
    else:
        print(f"GUARD WARNING: {guard_violations} violations detected!")
else:
    print("Guard log not found.")

# 2. Verify Output Tables & Figures
required_files = [
    "provenance/analysis_plan_scRNA.json",
    "provenance/marker_panel.json",
    "provenance/GSE248762_acquisition_summary.json",
    "provenance/GSE130888_suppl_file_manifest.json",
    "results/tables/E1_gse248762_sample_summary.csv",
    "results/tables/E2_filtering_summary.csv",
    "results/tables/E3_cell_type_composition.csv",
    "results/tables/E3_batch_mixing_metrics.csv",
    "results/tables/E4_celltype_testability_table.csv",
    "results/tables/E4_donor_cell_counts_by_type.csv",
    "results/tables/E5_gse130888_metadata_listing.csv",
    "results/figures/E2_qc_violins_before_after.png",
    "results/figures/E3_umap_by_celltype_and_group.png",
    "results/figures/E3_marker_dotplot.png"
]

file_status = {}
all_exist = True
for fpath in required_files:
    exists = os.path.exists(fpath)
    size = os.path.getsize(fpath) if exists else 0
    file_status[fpath] = {"exists": exists, "size_bytes": size}
    if not exists:
        all_exist = False
    print(f"  [{'EXISTS' if exists else 'MISSING'}] {fpath} ({size:,} bytes)")

# 3. Load Key Summary Tables
summary_data = {}
if os.path.exists("results/tables/E2_filtering_summary.csv"):
    df_e2 = pd.read_csv("results/tables/E2_filtering_summary.csv")
    summary_data["E2_QC"] = {
        "Total_Initial_Cells": int(df_e2["Initial_Cells"].sum()),
        "Total_Post_QC_Cells": int(df_e2["Post_QC_Cells"].sum()),
        "Overall_Retained_Pct": round(df_e2["Post_QC_Cells"].sum() / df_e2["Initial_Cells"].sum() * 100, 2),
        "Total_Doublets_Removed": int(df_e2["Doublets_Removed"].sum())
    }

if os.path.exists("results/tables/E3_cell_type_composition.csv"):
    df_e3 = pd.read_csv("results/tables/E3_cell_type_composition.csv")
    summary_data["E3_Cell_Types"] = df_e3.to_dict(orient="records")

if os.path.exists("results/tables/E4_celltype_testability_table.csv"):
    df_e4 = pd.read_csv("results/tables/E4_celltype_testability_table.csv")
    summary_data["E4_Testability"] = df_e4.to_dict(orient="records")

audit_summary = {
    "stage": "Stage 3 (Module E, Steps E0-E5)",
    "hub_gene_blind_policy_enforced": guard_passed,
    "guard_checks_count": guard_checks_count,
    "guard_violations_count": guard_violations,
    "all_required_files_generated": all_exist,
    "files": file_status,
    "summary": summary_data
}

with open("audit/stage3_audit_summary.json", "w") as f:
    json.dump(audit_summary, f, indent=2)

print("\nSaved audit/stage3_audit_summary.json")
print("Stage 3 audit verification complete.")
