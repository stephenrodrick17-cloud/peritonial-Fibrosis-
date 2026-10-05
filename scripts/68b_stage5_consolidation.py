"""
Script 68b: Stage 5 Manifest Generation & Cryptographic Integrity Audit
Consolidates all verified Stage 5 tables across:
  - 5F: miRNA profiling, evidence tiers, and exploratory intersections
  - 5G: lncRNA detection thresholds, credibility audit, and ceRNA network
  - 5H: GSE121372 calibrated descriptive analysis
Computes SHA256 hashes, mtime, shape, and generates results/tables/stage5_manifest.json.
"""

import os
import json
import hashlib
import datetime
import pandas as pd

ROOT = r"d:\Peritoneal Project"
TABLES_DIR = os.path.join(ROOT, "results", "tables")

stage5_files = {
    # Stage 5F: miRNA
    "F_multimir_validated_raw.csv": "multiMiR validated raw interactions (miRTarBase, TarBase, miRecords)",
    "F_multimir_predicted_raw.csv": "multiMiR predicted raw interactions (TargetScan, miRanda, etc.)",
    "F_multimir_provenance.csv": "multiMiR query provenance and database versions",
    "F_hub_to_mirna_tiered.csv": "Hub-miRNA interactions classified into Tier A (functional) and Tier B",
    "F_hub_to_mirna_tier_counts.csv": "Summary counts of miRNAs per evidence tier and per hub gene",
    "F_GSE182736_descriptive_logFC.csv": "GSE182736 human effluent exosome miRNA log2FC (filtered: raw>=10 in >=3 samples)",
    "F_GSE182736_sample_metadata.csv": "GSE182736 sample metadata and group definitions",
    "F_GSE182736_provenance.csv": "GSE182736 provenance documentation and extreme value audit",
    "F_GSE130387_cross_species_logFC.csv": "GSE130387 rodent PDF vs saline miRNA log2FC (CROSS-SPECIES, filtered)",
    "F_GSE130387_sample_metadata.csv": "GSE130387 sample metadata",
    "F_GSE130387_provenance.csv": "GSE130387 provenance documentation and probe prefix audit",
    "F_GSE130387_mirbase_mapping_note.csv": "GSE130387 miRBase stem homology mapping documentation",
    "F_mirna_intersection_summary.csv": "Summary and hypergeometric statistics (Tier A primary, Tier B sensitivity)",
    "F_mirna_hub_intersection.csv": "Overlapping miRNAs with log2FC, hub targets, and directionality check",
    
    # Stage 5G: lncRNA & ceRNA
    "G_lncRNA_threshold_sensitivity.csv": "GSE248762 single-cell stromal lncRNA detection across multiple cell/donor thresholds",
    "G_lncRNA_detected_list.csv": "Complete list of GENCODE v32 lncRNAs with stromal cell/donor detection and credibility flags",
    "G_lncRNA_GSE248762_detection_summary.csv": "Summary counts of detected, robust, and flagged lncRNAs",
    "G_lncRNA_provenance.csv": "GENCODE v32 lncRNA detection provenance",
    "G_encori_mirna_lncrna_raw.csv": "ENCORI live query raw lncRNA-miRNA interactions",
    "G_ceRNA_alternative_top10_tierA.csv": "Alternative objective selection: top 10 miRNAs ranked by Tier A functional targets",
    "G_ceRNA_network_predicted.csv": "Predicted ceRNA network edges with 16-donor stromal pseudobulk Spearman correlation",
    "G_ceRNA_network_summary.csv": "Summary of predicted ceRNA network nodes, hubs, and evidence classifications",
    "G_ceRNA_provenance.csv": "ceRNA network provenance and evidence constraints",
    
    # Stage 5H: GSE121372 HPMC
    "H_GSE121372_hub_fold_changes.csv": "GSE121372 HPMC 6h, 24h, and time-in-culture drift log2FC (calibrated wording)",
    "H_GSE121372_sample_metadata.csv": "GSE121372 sample metadata and conditions",
    "H_GSE121372_provenance.csv": "GSE121372 provenance and scientific constraint documentation"
}

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def main():
    print("=" * 80)
    print("SCRIPT 68b: STAGE 5 CONSOLIDATION & MANIFEST GENERATION")
    print("=" * 80)

    manifest_entries = {}
    all_exist = True

    for fname, desc in stage5_files.items():
        fpath = os.path.join(TABLES_DIR, fname)
        if not os.path.exists(fpath):
            print(f"ERROR: Missing file: {fname}")
            all_exist = False
            continue

        size = os.path.getsize(fpath)
        if size == 0:
            print(f"ERROR: Empty file: {fname}")
            all_exist = False
            continue

        sha256 = compute_sha256(fpath)
        stat = os.stat(fpath)
        mtime_utc = datetime.datetime.fromtimestamp(stat.st_mtime, tz=datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        
        # Read shape
        df = pd.read_csv(fpath)
        shape = list(df.shape)

        manifest_entries[fname] = {
            "path": os.path.abspath(fpath),
            "mtime_utc": mtime_utc,
            "sha256": sha256,
            "n_rows": shape[0],
            "shape": shape,
            "description": desc
        }

    assert all_exist, "One or more Stage 5 files are missing or empty!"

    manifest_path = os.path.join(TABLES_DIR, "stage5_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest_entries, f, indent=2)

    print(f"\nSuccessfully generated: {manifest_path}")
    print(f"Total files cataloged: {len(manifest_entries)}")

    # Verify each file against disk
    print("\n" + "=" * 80)
    print("VERIFYING EACH MANIFEST SHA256 AGAINST DISK")
    print("=" * 80)
    all_verified = True
    for fname, meta in manifest_entries.items():
        disk_sha = compute_sha256(meta["path"])
        match = (disk_sha == meta["sha256"])
        print(f"  {fname:42s} | SHA: {disk_sha[:12]}... | Match: {match}")
        if not match:
            all_verified = False

    assert all_verified, "Cryptographic hash mismatch detected!"
    print("\nALL STAGE 5 FILES CRYPTOGRAPHICALLY VERIFIED AGAINST DISK.")

if __name__ == "__main__":
    main()
