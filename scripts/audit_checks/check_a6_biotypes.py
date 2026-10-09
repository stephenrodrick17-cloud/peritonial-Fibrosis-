# scripts/audit_checks/check_a6_biotypes.py
import gzip
import re
import os
import urllib.request
import json
import pandas as pd
import anndata as ad

print("=== CHECK A6: 699 ENSG Features Biotype Breakdown ===")

adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad", backed="r")
var_names = pd.Series(adata_qc.var_names)
ensg_vars = var_names[var_names.str.startswith("ENSG")].tolist()
print(f"Total var_names starting with ENSG: {len(ensg_vars)}")

gtf_path = "data/raw/gencode.v32.long_noncoding_RNAs.gtf.gz"
gencode_lnc_ids = set()
with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.strip().split("\t")
        if len(p) > 8 and p[2] == "gene":
            m = re.search(r'gene_id "([^"]+)"', p[8])
            if m:
                gid = m.group(1).split(".")[0]
                gencode_lnc_ids.add(gid)

clean_ensg_vars = [x.split(".")[0] for x in ensg_vars]
clean_ensg_set = set(clean_ensg_vars)
in_gencode = clean_ensg_set.intersection(gencode_lnc_ids)
not_in_gencode = sorted(list(clean_ensg_set - gencode_lnc_ids))

print(f"ENSG var_names in GENCODE lncRNA: {len(in_gencode)}")
print(f"ENSG var_names NOT in GENCODE lncRNA: {len(not_in_gencode)}")

# Now let's determine their biotypes using Ensembl REST API lookup in batches of 100
print("\nQuerying Ensembl REST API for biotypes...")
biotypes = {}
batch_size = 100
for i in range(0, len(not_in_gencode), batch_size):
    batch = not_in_gencode[i:i+batch_size]
    data = json.dumps({"ids": batch}).encode("utf-8")
    req = urllib.request.Request(
        "https://rest.ensembl.org/lookup/id",
        data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"}
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            res = json.loads(resp.read().decode("utf-8"))
            for gid, info in res.items():
                if info is not None and "biotype" in info:
                    biotypes[gid] = info["biotype"]
                else:
                    biotypes[gid] = "unknown_or_not_found"
    except Exception as e:
        print(f"Batch {i} query error: {e}")
        for gid in batch:
            if gid not in biotypes:
                biotypes[gid] = "query_failed"

df_biotypes = pd.DataFrame({
    "gene_id": not_in_gencode,
    "biotype": [biotypes.get(gid, "unresolved") for gid in not_in_gencode]
})

print("\n=== Biotype breakdown of the 699 ENSG features ===")
counts = df_biotypes["biotype"].value_counts()
print(counts)
print(f"\nTotal resolved: {len(df_biotypes)} features")

df_biotypes.to_csv("results/tables/check_a6_699_ensg_biotypes.csv", index=False)
print("Saved results to results/tables/check_a6_699_ensg_biotypes.csv")
