"""
Script: Regenerate Stage 4 and Stage 5 Manifests with Portable Relative Paths

Formats manifests strictly with:
- Relative path (forward slashes: results/tables/<file>)
- SHA-256 hash (full 64-character lowercase hex)
- Byte count (bytes)
- Shape ([rows, cols])
Removes all absolute paths (e.g., D:\\...) and timestamps (mtimes).
"""

import os
import sys
import json
import hashlib
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES = os.path.join(ROOT, "results", "tables")

def compute_entry(rel_path):
    abs_path = os.path.join(ROOT, rel_path.replace("/", os.sep))
    if not os.path.exists(abs_path):
        raise FileNotFoundError(f"File not found: {abs_path}")
    
    with open(abs_path, "rb") as f:
        data = f.read()
    sha256 = hashlib.sha256(data).hexdigest()
    size = len(data)
    
    df = pd.read_csv(abs_path)
    shape = list(df.shape)
    
    return {
        "path": rel_path.replace("\\", "/"),
        "sha256": sha256,
        "bytes": size,
        "shape": shape
    }

def main():
    print("=" * 80)
    print("REGENERATING PORTABLE MANIFESTS (RELATIVE PATHS, SHA256, BYTES, SHAPES)")
    print("=" * 80)
    
    # 1. Stage 4 Manifest
    stage4_files = [
        "stage4_ambient_recomputed_true_cp10k.csv",
        "stage4_primary_edger_pseudobulk.csv",
        "stage4_primary_98_bh_family_sorted.csv",
        "stage4_qc_ceiling_sensitivity_comparison.csv",
        "stage4_sensitivity_edger_pseudobulk.csv"
    ]
    
    stage4_manifest = {}
    for fn in stage4_files:
        rel_p = f"results/tables/{fn}"
        stage4_manifest[fn] = compute_entry(rel_p)
    
    s4_path = os.path.join(TABLES, "stage4_manifest.json")
    with open(s4_path, "w", encoding="utf-8") as f:
        json.dump(stage4_manifest, f, indent=2)
    print(f"Saved portable Stage 4 manifest to: {s4_path} ({len(stage4_manifest)} files)")
    
    # 2. Stage 5 Manifest
    stage5_files = [
        "F_GSE130387_cross_species_logFC.csv",
        "F_GSE130387_mirbase_mapping_note.csv",
        "F_GSE130387_provenance.csv",
        "F_GSE130387_sample_metadata.csv",
        "F_GSE182736_descriptive_logFC.csv",
        "F_GSE182736_provenance.csv",
        "F_GSE182736_sample_metadata.csv",
        "F_hub_to_mirna_tier_counts.csv",
        "F_hub_to_mirna_tiered.csv",
        "F_mirna_hub_intersection.csv",
        "F_mirna_intersection_summary.csv",
        "F_multimir_predicted_raw.csv",
        "F_multimir_provenance.csv",
        "F_multimir_validated_raw.csv",
        "G_ceRNA_alternative_top10_tierA.csv",
        "G_ceRNA_network_predicted.csv",
        "G_ceRNA_network_summary.csv",
        "G_ceRNA_provenance.csv",
        "G_encori_mirna_lncrna_raw.csv",
        "G_lncRNA_GSE248762_detection_summary.csv",
        "G_lncRNA_detected_list.csv",
        "G_lncRNA_provenance.csv",
        "G_lncRNA_threshold_sensitivity.csv",
        "H_GSE121372_hub_fold_changes.csv",
        "H_GSE121372_provenance.csv",
        "H_GSE121372_sample_metadata.csv"
    ]
    
    stage5_manifest = {}
    for fn in stage5_files:
        rel_p = f"results/tables/{fn}"
        stage5_manifest[fn] = compute_entry(rel_p)
    
    s5_path = os.path.join(TABLES, "stage5_manifest.json")
    with open(s5_path, "w", encoding="utf-8") as f:
        json.dump(stage5_manifest, f, indent=2)
    print(f"Saved portable Stage 5 manifest to: {s5_path} ({len(stage5_manifest)} files)")

    # Also update results/stage5_manifest.json
    s5_legacy = os.path.join(ROOT, "results", "stage5_manifest.json")
    with open(s5_legacy, "w", encoding="utf-8") as f:
        json.dump(stage5_manifest, f, indent=2)
    print(f"Saved portable Stage 5 root manifest to: {s5_legacy}")
    
    print("\nManifest regeneration complete.")

if __name__ == "__main__":
    main()
