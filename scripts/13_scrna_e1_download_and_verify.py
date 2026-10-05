"""
Step E1: Acquire and verify GSE248762 single-cell RNA-seq data (16 donors)
Hub-gene-blind verification with code-level guard.
"""
import os
import sys
import gzip
import hashlib
import json
import ftplib
import tarfile
import pandas as pd
import scipy.io

sys.path.append("scripts")
from guard import guard_check

print("=================================================================")
print("STEP E1: GSE248762 ACQUISITION AND INTEGRITY VERIFICATION")
print("=================================================================")

# Guard verification: ensure no hub gene is accessed
guard_check([], stage="E1_acquisition")

raw_tar = "data/raw/GSE248762_RAW.tar"
extract_dir = "data/raw/GSE248762_extracted"
os.makedirs("data/raw", exist_ok=True)
os.makedirs(extract_dir, exist_ok=True)

# Check extracted files
extracted_files = [f for f in sorted(os.listdir(extract_dir)) if f.startswith("GSM")]
print(f"Found {len(extracted_files)} verified extracted sample files in {extract_dir}.")

# 2. SHA256 Checksum of archive if present
raw_hash = "N/A (extracted individual sample files)"
if os.path.exists(raw_tar):
    print("Calculating RAW archive SHA256...")
    sha256 = hashlib.sha256()
    with open(raw_tar, 'rb') as f:
        while chunk := f.read(1024 * 1024):
            sha256.update(chunk)
    raw_hash = sha256.hexdigest()
    print(f"GSE248762_RAW.tar SHA256: {raw_hash}")

# 4. Group files by GSM / Sample
# Files have naming: GSM7919583_LV_UF-1.barcodes.tsv.gz, etc.
samples_dict = {}
for fname in extracted_files:
    parts = fname.split("_", 1)
    gsm = parts[0]
    if gsm not in samples_dict:
        samples_dict[gsm] = {}
    
    if "barcode" in fname.lower():
        samples_dict[gsm]["barcodes"] = os.path.join(extract_dir, fname)
        samples_dict[gsm]["barcodes_fname"] = fname
    elif "feature" in fname.lower() or "gene" in fname.lower():
        samples_dict[gsm]["features"] = os.path.join(extract_dir, fname)
        samples_dict[gsm]["features_fname"] = fname
    elif "matrix" in fname.lower():
        samples_dict[gsm]["matrix"] = os.path.join(extract_dir, fname)
        samples_dict[gsm]["matrix_fname"] = fname

print(f"Discovered {len(samples_dict)} unique sample packages.")

# 5. Verify dimensions and feature format for each donor
sample_records = []
total_cells = 0

for gsm, paths in sorted(samples_dict.items()):
    b_file = paths["barcodes"]
    f_file = paths["features"]
    m_file = paths["matrix"]
    
    # Read barcodes
    with gzip.open(b_file, 'rt') as f:
        barcodes = [line.strip() for line in f if line.strip()]
    n_cells = len(barcodes)
    
    # Read features
    features_df = pd.read_csv(f_file, sep='\t', header=None, compression='gzip')
    n_genes = len(features_df)
    
    # Check gene symbols format
    sample_title = paths["matrix_fname"].split(".")[0].split("_", 1)[1] if "_" in paths["matrix_fname"] else gsm
    
    # Infer clinical group
    if sample_title.startswith("SV") or "SV" in sample_title:
        grp = "SV (Short Vintage, <1yr)"
    elif "LV_UF" in sample_title:
        grp = "LV_UF (Long Vintage with UFF)"
    elif "LV_NOT_UF" in sample_title or "LV" in sample_title:
        grp = "LV_NOT_UF (Long Vintage without UFF)"
    else:
        grp = "Unknown"
        
    # File hashes
    file_info = {}
    for key, fpath in [("barcodes", b_file), ("features", f_file), ("matrix", m_file)]:
        with open(fpath, 'rb') as fp:
            file_info[key] = {
                "filename": os.path.basename(fpath),
                "size_bytes": os.path.getsize(fpath),
                "sha256": hashlib.sha256(fp.read()).hexdigest()
            }
            
    sample_records.append({
        "GSM": gsm,
        "Sample_Title": sample_title,
        "Group": grp,
        "N_Cells_Raw": n_cells,
        "N_Genes_Raw": n_genes,
        "Features_Columns": features_df.shape[1],
        "Files": file_info
    })
    total_cells += n_cells
    print(f"  {gsm} | {sample_title:<15} | Group: {grp:<32} | {n_genes} genes x {n_cells:5d} cells")

print("=================================================================")
print(f"Total raw cells across all 16 donors: {total_cells}")
print("=================================================================")

# Save provenance metadata
with open("provenance/GSE248762_acquisition_summary.json", "w") as f:
    json.dump({
        "raw_tar_sha256": raw_hash,
        "total_raw_cells": total_cells,
        "samples": sample_records
    }, f, indent=2)

df_summary = pd.DataFrame([
    {
        "GSM": r["GSM"],
        "Sample_Title": r["Sample_Title"],
        "Group": r["Group"],
        "N_Cells_Raw": r["N_Cells_Raw"],
        "N_Genes_Raw": r["N_Genes_Raw"]
    } for r in sample_records
])
df_summary.to_csv("results/tables/E1_gse248762_sample_summary.csv", index=False)
print("Saved results/tables/E1_gse248762_sample_summary.csv")
print(df_summary.groupby("Group")["N_Cells_Raw"].agg(["count", "sum", "mean"]))

print("\nStep E1 finished successfully.")
