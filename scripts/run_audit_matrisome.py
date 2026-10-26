import pandas as pd, numpy as np, requests
from scipy.stats import hypergeom

R = "results/"
deg = pd.read_csv(R+"tables/reconciled_matrisome_DEGs.csv")
mod = pd.read_csv(R+"revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv")
expr = pd.read_csv(R+"tables/GSE125498_full_expression_matrix.csv")
bg = pd.read_csv(R+"revision/matrisome_volcano/GSE125498_matrisome_volcano_source_data.csv")

# 1. MBP / NAT6 in expression matrix?
print("MBP in matrix:", "MBP" in set(expr["Gene"]), "| NAT6:", "NAT6" in set(expr["Gene"]))

# 2. STRING edge check
e = pd.read_csv(R+"revision/matrisome_string_ppi/functional/string_functional_interactions.tsv", sep="\t")
print("columns:", list(e.columns)); print("edge rows:", len(e))
a, b = [c for c in e.columns if "preferredName" in c][:2]
deg_count = pd.concat([e[a], e[b]]).value_counts()
print("sum of degrees:", deg_count.sum(), "= 2 x edges?", deg_count.sum() == 2*len(e))
print("connected nodes:", len(deg_count), "| isolated of 46:", 46-len(deg_count))
print(deg_count.to_string())

# 3. STRING version + enrichment (identifiers joined by carriage return = %0d)
genes = deg["Gene"].tolist()
print(requests.get("https://string-db.org/api/json/version").text)
r = requests.post("https://string-db.org/api/json/ppi_enrichment",
    data={"identifiers": "\r".join(genes), "species": 9606, "required_score": 400})
print(r.status_code, r.text)

# 4. Cleaned matrisome list (remove NAT6, MBP)
drop = {"NAT6", "MBP"}
M = 11741
K = len(set(bg["Gene"]) - drop)
d = deg[~deg["Gene"].isin(drop)]
for name, n, k in [("all", 813, len(d)),
                   ("up", 496, (d["DEG_Status"]=="Upregulated").sum()),
                   ("down", 317, (d["DEG_Status"]=="Downregulated").sum())]:
    exp = n*K/M
    print(f"{name}: K={K} n={n} k={k} expected={exp:.2f} fold={k/exp:.3f} "
          f"P={hypergeom.sf(k-1, M, K, n):.4g}")
print(d["Matrisome_Division"].value_counts())

# 5. Cleaned module counts
m = mod[~mod["Gene"].isin(drop)] if "Gene" in mod.columns else mod
col = "Module" if "Module" in m.columns else m.columns[1]
print("in input:", m[col].notna().sum(),
      "| non-grey:", (m[col].notna() & (m[col] != "grey")).sum(),
      "| purple:", (m[col] == "purple").sum())
