import pandas as pd, requests, io
from collections import Counter

R = "results/"
deg = pd.read_csv(R + "tables/reconciled_matrisome_DEGs.csv")
genes = [g for g in deg["Gene"] if g not in {"NAT6", "MBP"}]
print("Genes:", len(genes))
base = "https://string-db.org/api"
common = {"identifiers": "\r".join(genes), "species": 9606,
          "required_score": 400, "caller_identity": "matrisome_pd_analysis"}

# Version
print("Version:", requests.get(f"{base}/json/version").text)

# Check each symbol maps to itself
ids = requests.post(f"{base}/json/get_string_ids", data={**common, "limit": 1}).json()
mapped = {x["queryItem"]: x["preferredName"] for x in ids}
print("Mapped:", len(mapped), "| unmapped:", sorted(set(genes) - set(mapped)))
print("Mismatched names:", {k: v for k, v in mapped.items() if k != v})

# PPI enrichment
r = requests.post(f"{base}/json/ppi_enrichment", data=common)
print("PPI enrichment:", r.status_code, r.text)

# Edges (functional network, score >= 0.4)
t = requests.post(f"{base}/tsv/network", data={**common, "network_type": "functional"})
e = pd.read_csv(io.StringIO(t.text), sep="\t")
pairs = {frozenset((a, b)) for a, b in zip(e["preferredName_A"], e["preferredName_B"]) if a != b}
deg_c = Counter()
for p in pairs:
    for g in p:
        deg_c[g] += 1
print("Edges (unique):", len(pairs), "| connected nodes:", len(deg_c),
      "| isolated:", len(genes) - len(deg_c))
print("Degree sum = 2 x edges?", sum(deg_c.values()) == 2 * len(pairs))
print("Top hubs:", deg_c.most_common(6))
print("Isolated genes:", sorted(set(genes) - set(deg_c)))

e.to_csv(R + "revision/matrisome_string_ppi/string_44_edges.tsv", sep="\t", index=False)
