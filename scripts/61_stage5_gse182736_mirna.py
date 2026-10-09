"""
Download GSE182736 via NCBI GEO text API (acc.cgi?form=text) and extract
supplementary file URLs, then download the miRNA expression data.
"""
import requests, os, gzip, re, hashlib, datetime, io, time
import pandas as pd
import numpy as np

RAWDIR = r"d:\Peritoneal Project\data\raw"
OUTDIR = r"d:\Peritoneal Project\results\tables"

# Step 1: Fetch the GEO soft text for GSE182736
geo_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182736&targ=self&form=text&view=full"
print("Fetching GEO soft record ...")
r = requests.get(geo_url, timeout=120)
r.raise_for_status()
soft_text = r.text
print(f"  Got {len(soft_text):,} chars")

# Show key metadata lines
print("\n--- Key metadata ---")
for line in soft_text.split("\n")[:100]:
    if any(k in line for k in ["!Series_title", "!Series_summary", "!Series_overall_design",
                                 "!Series_supplementary_file", "!Sample_geo_accession",
                                 "!Sample_title", "!Sample_supplementary_file"]):
        print(line.strip())

# Step 2: Extract supplementary file URLs
supp_urls = re.findall(r'!Series_supplementary_file\s*=\s*(.*)', soft_text)
supp_urls += re.findall(r'!Sample_supplementary_file\s*=\s*(.*)', soft_text)
supp_urls = [u.strip() for u in supp_urls if u.strip()]
print(f"\nSupplementary file URLs ({len(supp_urls)}):")
for u in supp_urls:
    print(f"  {u}")

# Step 3: Try to download supplementary files
for surl in supp_urls:
    fname = surl.split("/")[-1].strip()
    local = os.path.join(RAWDIR, fname)
    if os.path.exists(local) and os.path.getsize(local) > 1000:
        print(f"  Already downloaded: {fname}")
        continue
    print(f"\nDownloading {fname} ...")
    try:
        dr = requests.get(surl.strip(), stream=True, timeout=300)
        dr.raise_for_status()
        with open(local, "wb") as f:
            for chunk in dr.iter_content(65536):
                f.write(chunk)
        print(f"  Saved {local} ({os.path.getsize(local):,} bytes)")
    except Exception as e:
        print(f"  FAILED: {e}")

# Step 4: Parse sample metadata from the SOFT text
print("\n--- Parsing sample metadata ---")
samples = {}
cur_acc = None
for line in soft_text.split("\n"):
    line = line.strip()
    m_acc = re.match(r'\^SAMPLE\s*=\s*(\S+)', line)
    if m_acc:
        cur_acc = m_acc.group(1)
        samples[cur_acc] = {}
    elif cur_acc:
        for key in ["!Sample_title", "!Sample_characteristics_ch1", "!Sample_source_name_ch1"]:
            if line.startswith(key):
                val = line.split("=", 1)[1].strip()
                samples[cur_acc].setdefault(key, []).append(val)

print(f"Samples found: {len(samples)}")
for gsm, info in samples.items():
    title = info.get("!Sample_title", [""])[0]
    chars = info.get("!Sample_characteristics_ch1", [])
    print(f"  {gsm}: {title} | chars: {chars}")

# Build sample dataframe
rows = []
for gsm, info in samples.items():
    title = info.get("!Sample_title", [""])[0]
    chars = " | ".join(info.get("!Sample_characteristics_ch1", []))
    source = info.get("!Sample_source_name_ch1", [""])[0]
    # Infer group
    t_low = (title + " " + chars + " " + source).lower()
    if re.search(r'uf.fail|ultrafilt.*fail|failure', t_low):
        group = "UF_failure"
    elif re.search(r'non.uf|non_uf|control|ctrl|normal', t_low):
        group = "non_UF"
    else:
        group = "unknown"
    rows.append({"GSM": gsm, "Title": title, "Characteristics": chars,
                 "Source": source, "Group": group})

sample_df = pd.DataFrame(rows)
print(f"\nSample table:\n{sample_df.to_string(index=False)}")
sample_df.to_csv(os.path.join(OUTDIR, "F_GSE182736_sample_metadata.csv"), index=False)

# Step 5: Try to parse any downloaded supplementary files
supp_files = [os.path.join(RAWDIR, surl.split("/")[-1].strip())
              for surl in supp_urls]
supp_files = [f for f in supp_files if os.path.exists(f) and os.path.getsize(f) > 1000]

print(f"\nAvailable supplementary files: {supp_files}")

expr_all = None
for sf in supp_files:
    fname = os.path.basename(sf)
    print(f"\nTrying to parse: {fname}")
    try:
        if fname.endswith(".gz"):
            with gzip.open(sf, "rt", errors="replace") as fh:
                content = fh.read()
        else:
            with open(sf, "r", errors="replace") as fh:
                content = fh.read()
        # Try as tab-delimited
        df = pd.read_csv(io.StringIO(content), sep="\t", index_col=0)
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {df.columns.tolist()}")
        print(df.head(5).to_string())
        expr_all = df
    except Exception as e:
        print(f"  Parse error: {e}")
        # Try comma-delimited
        try:
            df = pd.read_csv(sf, index_col=0)
            print(f"  CSV shape: {df.shape}")
            expr_all = df
        except Exception as e2:
            print(f"  CSV parse also failed: {e2}")

# Step 6: Provenance regardless
prov = {
    "accession": "GSE182736",
    "geo_text_url": geo_url,
    "fetch_date": datetime.datetime.utcnow().isoformat() + "Z",
    "supp_urls_found": len(supp_urls),
    "supp_files_downloaded": len(supp_files),
    "n_samples_metadata": len(samples),
    "groups_found": sample_df["Group"].value_counts().to_dict() if not sample_df.empty else {},
    "note": "n=3 vs 3; descriptive only"
}
prov_df = pd.DataFrame([prov])
prov_df.to_csv(os.path.join(OUTDIR, "F_GSE182736_provenance.csv"), index=False)
print(f"\nProvenance: {prov}")
print("\nDone.")
