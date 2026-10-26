import pandas as pd, requests, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "results/"
out = R + "revision/wgcna_paper_figures/"
drop = {"NAT6", "MBP"}

# ---------- Fig 3D: matrisome DEGs per module (cleaned) ----------
mod = pd.read_csv(R + "revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv")
gcol = "Gene.symbol" if "Gene.symbol" in mod.columns else "Annotated Gene"
m = mod[~mod[gcol].isin(drop)].copy()
m["Module"] = m["Module"].fillna("not in network input")
tab = (m.groupby(["Module", "DEG_Status"]).size().unstack(fill_value=0)
         .reindex(columns=["Upregulated", "Downregulated"], fill_value=0))
tab["total"] = tab.sum(axis=1)
tab = tab.sort_values("total", ascending=False)
print(tab)
print("Total:", int(tab["total"].sum()))

fig, ax = plt.subplots(figsize=(8, 4.5))
x = range(len(tab))
ax.bar(x, tab["Upregulated"], color="#0072B2", label="Up in LPD")
ax.bar(x, tab["Downregulated"], bottom=tab["Upregulated"], color="#D55E00", label="Down in LPD")
ax.set_xticks(list(x))
ax.set_xticklabels(tab.index, rotation=45, ha="right")
ax.set_ylabel("Number of matrisome DEGs")
ax.set_title("Matrisome DEGs by WGCNA module (n = 44)")
ax.legend(frameon=False)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
plt.tight_layout()
for ext in ("png", "svg", "pdf"):
    plt.savefig(out + f"Fig3D_matrisome_module_intersection_cleaned.{ext}", dpi=300)

# ---------- Fig 3F: STRING network of the cleaned 44 genes ----------
deg = pd.read_csv(R + "tables/reconciled_matrisome_DEGs.csv")
genes = [g for g in deg["Gene"] if g not in drop]
print("STRING genes:", len(genes))
p = {"identifiers": "\r".join(genes), "species": 9606, "required_score": 400,
     "network_flavor": "confidence", "caller_identity": "matrisome_pd_analysis"}
sdir = R + "revision/matrisome_string_ppi/"
for kind, ext in (("highres_image", "png"), ("svg", "svg")):
    r = requests.post(f"https://string-db.org/api/{kind}/network", data=p)
    print(kind, r.status_code)
    if r.ok:
        open(sdir + f"Fig3F_STRING_44_cleaned.{ext}", "wb").write(r.content)
