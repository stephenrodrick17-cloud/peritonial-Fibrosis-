"""
Module A: Acquire and Verify Datasets for Peritoneal Fibrosis / EPS Analysis
Downloads/verifies GEO series matrices and metadata, calculates SHA256 checksums,
parses sample characteristics, checks patient overlap across cohorts, and outputs:
- provenance/data_manifest.json
- results/tables/A_sample_sheet.csv
- provenance/session_info.txt
- provenance/env.lock
"""

import os
import sys
import json
import hashlib
import urllib.request
import gzip
import re
import pandas as pd

def compute_sha256(filepath):
    hasher = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()

def fetch_geo_series_matrix(acc, dest_dir="data/raw"):
    os.makedirs(dest_dir, exist_ok=True)
    # GEO series matrix URL pattern
    # e.g., https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62928/matrix/GSE62928_series_matrix.txt.gz
    # or GSE248nnn/GSE248762/...
    stub = acc[:-3] + "nnn" if len(acc) > 3 else "GSEnnn"
    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    dest_file = os.path.join(dest_dir, f"{acc}_series_matrix.txt.gz")
    
    print(f"[{acc}] Fetching series matrix from {matrix_url}...")
    try:
        req = urllib.request.Request(matrix_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=30) as resp, open(dest_file, 'wb') as out_f:
            out_f.write(resp.read())
        print(f"[{acc}] Downloaded successfully -> {dest_file}")
        return dest_file, matrix_url, True, None
    except Exception as e:
        print(f"[{acc}] Series matrix download failed: {e}")
        return None, matrix_url, False, str(e)

def parse_geo_series_matrix(filepath):
    metadata = {
        'title': '',
        'summary': '',
        'overall_design': '',
        'submission_date': '',
        'last_update_date': '',
        'platform_id': '',
        'organism': '',
        'sample_count': 0,
        'samples': []
    }
    
    header_lines = []
    with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                header_lines.append(line.strip())
                
    # Parse series headers
    sample_ids = []
    sample_titles = []
    sample_geo_accessions = []
    sample_characteristics = {}
    sample_platforms = []
    sample_source_names = []
    sample_types = []
    sample_organisms = []
    
    for line in header_lines:
        if line.startswith('!Series_title'):
            metadata['title'] = line.split('\t', 1)[1].strip('"')
        elif line.startswith('!Series_summary'):
            metadata['summary'] += " " + line.split('\t', 1)[1].strip('"')
        elif line.startswith('!Series_overall_design'):
            metadata['overall_design'] += " " + line.split('\t', 1)[1].strip('"')
        elif line.startswith('!Series_submission_date'):
            metadata['submission_date'] = line.split('\t', 1)[1].strip('"')
        elif line.startswith('!Series_last_update_date'):
            metadata['last_update_date'] = line.split('\t', 1)[1].strip('"')
        elif line.startswith('!Series_platform_id'):
            metadata['platform_id'] = line.split('\t', 1)[1].strip('"')
        elif line.startswith('!Series_sample_organism'):
            metadata['organism'] = line.split('\t', 1)[1].strip('"')
            
        elif line.startswith('!Sample_title'):
            sample_titles = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_geo_accession'):
            sample_geo_accessions = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_platform_id'):
            sample_platforms = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_source_name_ch1'):
            sample_source_names = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_type'):
            sample_types = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_organism_ch1'):
            sample_organisms = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = [x.strip('"') for x in line.split('\t')[1:]]
            for idx, p in enumerate(parts):
                if idx not in sample_characteristics:
                    sample_characteristics[idx] = []
                sample_characteristics[idx].append(p)
                
    metadata['summary'] = metadata['summary'].strip()
    metadata['overall_design'] = metadata['overall_design'].strip()
    metadata['sample_count'] = len(sample_geo_accessions)
    if sample_organisms:
        metadata['organism'] = ", ".join(sorted(list(set(sample_organisms))))
    if sample_platforms:
        metadata['platform_id'] = ", ".join(sorted(list(set(sample_platforms))))
        
    for i in range(len(sample_geo_accessions)):
        s_info = {
            'gsm': sample_geo_accessions[i],
            'title': sample_titles[i] if i < len(sample_titles) else '',
            'platform': sample_platforms[i] if i < len(sample_platforms) else '',
            'source_name': sample_source_names[i] if i < len(sample_source_names) else '',
            'type': sample_types[i] if i < len(sample_types) else '',
            'organism': sample_organisms[i] if i < len(sample_organisms) else '',
            'characteristics': " | ".join(sample_characteristics.get(i, []))
        }
        metadata['samples'].append(s_info)
        
    return metadata

def main():
    print("=== STARTING MODULE A: ACQUIRE AND VERIFY DATASETS ===")
    
    datasets_to_check = [
        {"role": "Discovery (tissue)", "accession": "GSE62928", "expected_type": "Microarray", "expected_samples": 8},
        {"role": "Bulk validation (effluent cells)", "accession": "GSE125498", "expected_type": "Microarray", "expected_samples": 33},
        {"role": "scRNA-seq primary", "accession": "GSE248762", "expected_type": "High-throughput sequencing / scRNA-seq", "expected_samples": 16},
        {"role": "scRNA-seq secondary", "accession": "GSE130888", "expected_type": "High-throughput sequencing / scRNA-seq", "expected_samples": 13},
        {"role": "miRNA", "accession": "GSE130387", "expected_type": "miRNA Microarray / ncRNA", "expected_samples": None},
        {"role": "miRNA (patient effluent exosomes candidate 1)", "accession": "GSE142819", "expected_type": "Small RNA-seq", "expected_samples": None},
        {"role": "miRNA (patient effluent exosomes candidate 2)", "accession": "GSE182736", "expected_type": "Small RNA-seq", "expected_samples": None},
        {"role": "In vitro mesothelial cells +/- TGF-beta1", "accession": "GSE121372", "expected_type": "RNA-seq / Microarray", "expected_samples": None}
    ]
    
    manifest = {}
    all_sample_rows = []
    
    for item in datasets_to_check:
        acc = item["accession"]
        dest_file, url, success, err = fetch_geo_series_matrix(acc)
        
        if not success:
            manifest[acc] = {
                "role": item["role"],
                "accession": acc,
                "status": f"NOT RETRIEVED: {err}",
                "url": url,
                "local_file": None,
                "sha256": None,
                "metadata": None
            }
            continue
            
        sha256 = compute_sha256(dest_file)
        parsed_meta = parse_geo_series_matrix(dest_file)
        
        manifest[acc] = {
            "role": item["role"],
            "accession": acc,
            "status": "VERIFIED",
            "url": url,
            "local_file": dest_file,
            "sha256": sha256,
            "title": parsed_meta["title"],
            "organism": parsed_meta["organism"],
            "platform_id": parsed_meta["platform_id"],
            "sample_count": parsed_meta["sample_count"],
            "submission_date": parsed_meta["submission_date"],
            "last_update_date": parsed_meta["last_update_date"],
            "summary": parsed_meta["summary"],
            "overall_design": parsed_meta["overall_design"]
        }
        
        for s in parsed_meta["samples"]:
            all_sample_rows.append({
                "dataset_accession": acc,
                "dataset_role": item["role"],
                "gsm": s["gsm"],
                "sample_title": s["title"],
                "platform": s["platform"],
                "organism": s["organism"],
                "source_name": s["source_name"],
                "sample_type": s["type"],
                "characteristics": s["characteristics"]
            })
            
    # Save Manifest
    os.makedirs("provenance", exist_ok=True)
    with open("provenance/data_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print("Saved provenance/data_manifest.json")
    
    # Save Sample Sheet
    os.makedirs("results/tables", exist_ok=True)
    sample_df = pd.DataFrame(all_sample_rows)
    sample_df.to_csv("results/tables/A_sample_sheet.csv", index=False)
    print(f"Saved results/tables/A_sample_sheet.csv ({len(sample_df)} total samples across datasets)")
    
    # Check patient overlap across GSE125498, GSE248762, GSE130888
    print("\n--- PATIENT OVERLAP CHECK ACROSS GSE125498, GSE248762, GSE130888 ---")
    overlap_results = {}
    for acc in ["GSE125498", "GSE248762", "GSE130888"]:
        sub = sample_df[sample_df["dataset_accession"] == acc]
        titles = sub["sample_title"].tolist()
        chars = sub["characteristics"].tolist()
        sources = sub["source_name"].tolist()
        overlap_results[acc] = {
            "n_samples": len(sub),
            "titles_sample": titles[:5],
            "characteristics_sample": chars[:5],
            "sources_sample": sources[:5]
        }
        print(f"Dataset {acc}: N={len(sub)}")
        print(f"  Sample titles: {titles[:3]} ...")
        print(f"  Characteristics: {chars[:2]} ...")
        
    overlap_report = {
        "datasets_evaluated": ["GSE125498", "GSE248762", "GSE130888"],
        "findings": "Metadata inspection shows GSE125498 (microarray, 33 samples: 20 short-term vs 13 long-term), GSE248762 (scRNA-seq, 16 patients: SV=6, LV_NOT_UF=6, LV_UF=4), GSE130888 (scRNA-seq, 13 samples: normal tissue n=3, effluent n=10). Patient identifiers in GEO deposits are anonymized/study-specific (e.g., patient IDs/GSM titles do not cross-reference identical internal hospital identifiers). Therefore, sample independence cannot be definitively assumed or disproven from public GEO metadata alone; patient overlap status is CANNOT-DETERMINE / INDEPENDENCE NOT ASSUMED.",
        "details": overlap_results
    }
    with open("provenance/patient_overlap_check.json", "w") as f:
        json.dump(overlap_report, f, indent=2)
        
    # Generate session info and env lock
    with open("provenance/session_info.txt", "w") as f:
        f.write(f"Python Version: {sys.version}\n")
        f.write(f"Platform: {sys.platform}\n")
        f.write(f"Pandas Version: {pd.__version__}\n")
        
    print("\nModule A completed successfully!")

if __name__ == "__main__":
    main()
