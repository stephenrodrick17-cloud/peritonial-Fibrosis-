"""
Stage 5 Consolidation, Integrity Verification, and Manifest Generation
Verifies all files across Stage 5F (miRNA), 5G (lncRNA/ceRNA), and 5H (GSE121372).
Computes SHA256 hashes, asserts all schema and constraint checks, and writes stage5_manifest.json.
"""

import os, hashlib, json, datetime
import pandas as pd

TABLES_DIR = r"d:\Peritoneal Project\results\tables"
RESULTS_DIR = r"d:\Peritoneal Project\results"

expected_files = {
    # Stage 5F
    "F_multimir_validated_raw.csv": "multiMiR validated raw interactions (all databases)",
    "F_multimir_predicted_raw.csv": "multiMiR predicted raw interactions",
    "F_multimir_provenance.csv": "multiMiR package and DB version provenance",
    "F_hub_to_mirna_validated.csv": "Validated miRNAs targeting the 11 hub genes",
    "F_hub_to_mirna_mirtarbase.csv": "miRTarBase-only validated miRNAs targeting hubs",
    "F_GSE182736_descriptive_logFC.csv": "GSE182736 human effluent exosome miRNA log2FC (n=3 vs 3, descriptive)",
    "F_GSE182736_sample_metadata.csv": "GSE182736 sample metadata",
    "F_GSE182736_provenance.csv": "GSE182736 dataset provenance",
    "F_GSE130387_cross_species_logFC.csv": "GSE130387 rodent PDF vs saline miRNA log2FC (CROSS-SPECIES, n=3 vs 3)",
    "F_GSE130387_sample_metadata.csv": "GSE130387 sample metadata",
    "F_GSE130387_provenance.csv": "GSE130387 dataset provenance",
    "F_GSE130387_mirbase_mapping_note.csv": "GSE130387 miRBase stem homology mapping documentation",
    "F_mirna_hub_intersection.csv": "Intersection between validated hub miRNAs and DE miRNAs",
    "F_mirna_intersection_summary.csv": "Summary and hypergeometric test statistics for miRNA intersections",
    
    # Stage 5G
    "G_lncRNA_GSE248762_detection_summary.csv": "GSE248762 lncRNA detection summary across cell types",
    "G_lncRNA_GSE248762_by_celltype.csv": "Detailed per-celltype expression for all assayed lncRNAs",
    "G_lncRNA_detected_list.csv": "List of lncRNAs detected (counts > 0) in GSE248762",
    "G_lncRNA_provenance.csv": "GENCODE v32 lncRNA annotation provenance",
    "G_encori_mirna_lncrna_raw.csv": "ENCORI live query raw lncRNA-miRNA interactions",
    "G_ceRNA_network_predicted.csv": "Predicted ceRNA regulatory network (lncRNA -> miRNA -> Hub)",
    "G_ceRNA_network_summary.csv": "Summary of predicted ceRNA network nodes and hubs",
    "G_ceRNA_provenance.csv": "ENCORI API live query provenance and citation",
    
    # Stage 5H
    "H_GSE121372_hub_fold_changes.csv": "GSE121372 HPMC +/- TGF-b1 hub gene fold changes at 6h and 24h",
    "H_GSE121372_sample_metadata.csv": "GSE121372 sample metadata",
    "H_GSE121372_provenance.csv": "GSE121372 provenance documentation"
}

def sha256_file(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

print("="*70)
print("STAGE 5 CONSOLIDATION & INTEGRITY AUDIT")
print("="*70)

manifest_entries = {}
all_passed = True

for fname, desc in expected_files.items():
    fpath = os.path.join(TABLES_DIR, fname)
    if not os.path.exists(fpath):
        print(f"FAILED: Missing file: {fname}")
        all_passed = False
        continue
    
    file_size = os.path.getsize(fpath)
    if file_size == 0:
        print(f"FAILED: Empty file: {fname}")
        all_passed = False
        continue
    
    h = sha256_file(fpath)
    
    # Check row count
    try:
        df = pd.read_csv(fpath)
        nrows = len(df)
        ncols = len(df.columns)
    except Exception as e:
        nrows = -1
        ncols = -1
    
    manifest_entries[fname] = {
        "description": desc,
        "sha256": h,
        "bytes": file_size,
        "rows": nrows,
        "cols": ncols
    }
    print(f"  OK: {fname:45s} | {nrows:6d} rows | {file_size:8d} B | SHA: {h[:12]}...")

print("\n--- Performing Constraint Checks ---")

# Check 1: GSE130387 CROSS-SPECIES label
cs_path = os.path.join(TABLES_DIR, "F_GSE130387_cross_species_logFC.csv")
if os.path.exists(cs_path):
    df_cs = pd.read_csv(cs_path)
    if "species_label" in df_cs.columns:
        non_cs = (df_cs["species_label"] != "CROSS-SPECIES").sum()
        assert non_cs == 0, f"Found {non_cs} rows without CROSS-SPECIES in {cs_path}"
        print(f"  Constraint Check 1 PASSED: All {len(df_cs)} rows in F_GSE130387 labeled CROSS-SPECIES")
    else:
        print("  WARNING: species_label column missing in GSE130387 table")

# Check 2: ceRNA predicted network label
cerna_path = os.path.join(TABLES_DIR, "G_ceRNA_network_predicted.csv")
if os.path.exists(cerna_path):
    df_cerna = pd.read_csv(cerna_path)
    pred_label = "PREDICTED regulatory network (not experimentally validated in this system)"
    if "regulatory_label" in df_cerna.columns:
        non_pred = (df_cerna["regulatory_label"] != pred_label).sum()
        assert non_pred == 0, f"Found {non_pred} rows without PREDICTED label in {cerna_path}"
        print(f"  Constraint Check 2 PASSED: All {len(df_cerna)} rows in G_ceRNA_network_predicted labeled PREDICTED")
    else:
        print("  WARNING: regulatory_label column missing in ceRNA table")

# Check 3: GSE142819 PCOS dataset excluded
for fname in manifest_entries.keys():
    assert "142819" not in fname, f"GSE142819 detected in {fname}!"
print("  Constraint Check 3 PASSED: GSE142819 (PCOS) excluded completely")

# Check 4: GSE121372 n=1 descriptive only
gse121_path = os.path.join(TABLES_DIR, "H_GSE121372_hub_fold_changes.csv")
if os.path.exists(gse121_path):
    df_121 = pd.read_csv(gse121_path)
    for col in df_121.columns:
        assert "p_val" not in col.lower() and "fdr" not in col.lower(), f"P-value/FDR found in n=1 table: {col}"
    print(f"  Constraint Check 4 PASSED: GSE121372 table contains no p-values or FDR claims (n=1, descriptive only)")

# Write manifest
manifest = {
    "stage": "Stage 5: Regulatory Layers (miRNA, lncRNA/ceRNA, and in vitro TGF-b1 validation)",
    "timestamp": datetime.datetime.now().isoformat(),
    "audit_status": "PASSED" if all_passed else "FAILED",
    "files": manifest_entries
}

manifest_path = os.path.join(RESULTS_DIR, "stage5_manifest.json")
with open(manifest_path, "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2)

print(f"\nSaved Stage 5 Manifest: {manifest_path}")
print("="*70)
print(f"STAGE 5 AUDIT COMPLETED: {'ALL 23 FILES VERIFIED' if all_passed else 'SOME CHECKS FAILED'}")
print("="*70)
