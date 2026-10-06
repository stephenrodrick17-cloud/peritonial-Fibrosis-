"""
Master Reproducible Pipeline Runner
Peritoneal Dialysis Fibrosis Multi-Omic & Single-Cell Study

Executes end-to-end pipeline stages:
  - Stage 4: Pseudobulk Single-Cell edgeR Modelling & Sensitivity Suite
  - Stage 5: Regulatory Layers (GSE130387 rodent miRNA, GSE121372 in vitro HPMC audit)
  - Integrity & Manifest Verification
"""

import os
import sys
import subprocess
import hashlib
import json
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
RESULTS = os.path.join(ROOT, "results")
TABLES = os.path.join(RESULTS, "tables")

RSCRIPT = r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
if not os.path.exists(RSCRIPT):
    RSCRIPT = "Rscript"

def print_header(title):
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80 + "\n")

def run_cmd(cmd_list, desc):
    print(f">>> Running: {desc}...")
    print(f"    Command: {' '.join(cmd_list)}")
    res = subprocess.run(cmd_list, cwd=ROOT, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"ERROR in {desc} (Exit code {res.returncode}):")
        print(res.stderr[:2000])
        sys.exit(res.returncode)
    else:
        print(f"    [SUCCESS] Completed {desc}.\n")

def check_file(path, expected_rows=None):
    if not os.path.exists(path):
        print(f"FAIL: Missing expected output: {path}")
        sys.exit(1)
    df = pd.read_csv(path)
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    sz = os.path.getsize(path)
    print(f"  Verified {os.path.basename(path)}: {len(df)} rows, {sz} bytes, SHA256: {h[:16]}...")
    if expected_rows is not None and len(df) != expected_rows:
        print(f"FAIL: Row count mismatch for {path}: got {len(df)}, expected {expected_rows}")
        sys.exit(1)

def main():
    print_header("PERITONEAL FIBROSIS MASTER REPRODUCIBLE PIPELINE")
    print(f"Root Directory: {ROOT}")
    print(f"Python:         {sys.version.split()[0]}")
    print(f"Rscript:        {RSCRIPT}")

    # Stage 4: Run Final edgeR Single-Cell Pseudobulk
    print_header("STAGE 4: SINGLE-CELL PSEUDOBULK edgeR & SENSITIVITY SUITE")
    run_cmd([RSCRIPT, os.path.join(SCRIPTS, "50_stage4_final.R")], "Stage 4 edgeR Final")
    run_cmd([sys.executable, os.path.join(SCRIPTS, "51_stage4_tables.py")], "Stage 4 Tables & Manifest")
    run_cmd([sys.executable, os.path.join(SCRIPTS, "52_stage4_bh_family_sync.py")], "Stage 4 BH Family Sync (m=98)")

    # Stage 5: Regulatory Layers
    print_header("STAGE 5: REGULATORY AUDITS & DESCRIPTIVE VALIDATION")
    run_cmd([sys.executable, os.path.join(SCRIPTS, "62b_stage5_gse130387_species_audit.py")], "GSE130387 Species Audit")
    run_cmd([sys.executable, os.path.join(SCRIPTS, "67b_stage5_gse121372_hpmc_audit.py")], "GSE121372 HPMC In Vitro Audit")
    run_cmd([sys.executable, os.path.join(SCRIPTS, "75_stage5_audit_master_resolution.py")], "Stage 5 Audit Master Resolution")
    run_cmd([sys.executable, os.path.join(SCRIPTS, "76_stage5_final_corrections_readback.py")], "Stage 5 Final Readbacks")
    run_cmd([sys.executable, os.path.join(SCRIPTS, "make_manifests_portable.py")], "Regenerate Portable Manifests")


    # Assertions & Verification
    print_header("FINAL CRYPTOGRAPHIC INTEGRITY VERIFICATION")
    check_file(os.path.join(TABLES, "stage4_primary_edger_pseudobulk.csv"), 165)
    check_file(os.path.join(TABLES, "stage4_primary_98_bh_family_sorted.csv"), 98)
    check_file(os.path.join(TABLES, "stage4_sensitivity_edger_pseudobulk.csv"), 792)
    check_file(os.path.join(TABLES, "F_GSE130387_cross_species_logFC.csv"), 479)
    check_file(os.path.join(TABLES, "H_GSE121372_hub_fold_changes.csv"), 13)

    print("\nAll pipeline assertions PASSED successfully!")

if __name__ == "__main__":
    main()
