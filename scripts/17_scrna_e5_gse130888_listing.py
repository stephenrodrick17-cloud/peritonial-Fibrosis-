"""
Step E5: GSE130888 GEO Metadata & Supplementary File Manifest (Zero Large Download)
Inspects external scRNA-seq dataset metadata without downloading large matrix archives.
"""
import os
import sys
import json
import urllib.request
import pandas as pd

sys.path.append("scripts")
from guard import guard_check

print("=================================================================")
print("STEP E5: GSE130888 GEO METADATA & SUPPLEMENTARY MANIFEST")
print("=================================================================")

# Guard check: Ensure no hub genes are accessed
guard_check([], stage="E5_start")

os.makedirs("provenance", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# Query NCBI GEO for GSE130888
url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130888&targ=self&form=text&view=full"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
content = urllib.request.urlopen(req).read().decode('utf-8', errors='ignore')

series_title = ""
series_summary = ""
series_design = ""
series_samples = []
suppl_files = []

for line in content.splitlines():
    if line.startswith("!Series_title ="):
        series_title = line.split("=", 1)[1].strip()
    elif line.startswith("!Series_summary ="):
        series_summary += " " + line.split("=", 1)[1].strip()
    elif line.startswith("!Series_overall_design ="):
        series_design += " " + line.split("=", 1)[1].strip()
    elif line.startswith("!Series_sample_id ="):
        series_samples.append(line.split("=", 1)[1].strip())
    elif line.startswith("!Series_supplementary_file ="):
        suppl_files.append(line.split("=", 1)[1].strip())

print(f"GSE130888 Title: {series_title}")
print(f"Total Samples: {len(series_samples)}")
print(f"Supplementary files listed: {len(suppl_files)}")
for sf in suppl_files:
    print(f"  {sf}")

# Fetch sample details for first 5 samples
sample_rows = []
for gsm in series_samples:
    s_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full"
    s_req = urllib.request.Request(s_url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    s_content = urllib.request.urlopen(s_req).read().decode('utf-8', errors='ignore')
    
    stitle = ""
    schars = []
    for sline in s_content.splitlines():
        if sline.startswith("!Sample_title ="):
            stitle = sline.split("=", 1)[1].strip()
        elif sline.startswith("!Sample_characteristics_ch1 ="):
            schars.append(sline.split("=", 1)[1].strip())
            
    sample_rows.append({
        "GSM": gsm,
        "Title": stitle,
        "Characteristics": " | ".join(schars)
    })

manifest = {
    "Accession": "GSE130888",
    "Title": series_title,
    "Summary": series_summary.strip(),
    "Design": series_design.strip(),
    "Total_Samples": len(series_samples),
    "Supplementary_Files": suppl_files,
    "Samples": sample_rows
}

with open("provenance/GSE130888_suppl_file_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)

df_samples = pd.DataFrame(sample_rows)
df_samples.to_csv("results/tables/E5_gse130888_metadata_listing.csv", index=False)

print("\nSaved provenance/GSE130888_suppl_file_manifest.json")
print("Saved results/tables/E5_gse130888_metadata_listing.csv")
print(f"Processed metadata for {len(df_samples)} samples in GSE130888.")
print("\nStep E5 finished successfully.")
