import gzip, re

lnc_genes = {}
with gzip.open(r"data/raw/gencode.v32.long_noncoding_RNAs.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) > 8 and fields[2] == "gene":
            gid_m = re.search(r'gene_id "([^"]+)"', fields[8])
            gn_m = re.search(r'gene_name "([^"]+)"', fields[8])
            gt_m = re.search(r'gene_type "([^"]+)"', fields[8])
            if gid_m and gn_m:
                gid = gid_m.group(1).split(".")[0]
                lnc_genes[gid] = (gn_m.group(1), gt_m.group(1) if gt_m else "lncRNA")

print(f"Total lncRNA genes in GENCODE v32: {len(lnc_genes)}")

features_gse = []
with gzip.open(r"data/raw/GSE248762_extracted/GSM7919583_LV_UF-1.features.tsv.gz", "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        features_gse.append((parts[0], parts[1]))

print(f"Total features in GSE248762: {len(features_gse)}")
matches = [f for f in features_gse if f[0] in lnc_genes]
print(f"lncRNA features matching GSE248762: {len(matches)}")
print("First 10 matching lncRNAs:")
for m in matches[:10]:
    print(f"  {m[0]} -> {m[1]} ({lnc_genes[m[0]][1]})")
