import pandas as pd
R = "results/"
mod = pd.read_csv(R+"revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv")
print(mod.columns.tolist())
gcol = [c for c in mod.columns if "gene" in c.lower()][0]
m = mod[~mod[gcol].isin({"NAT6","MBP"})]
print("rows:", len(m), "| in input:", m["Module"].notna().sum(), "| non-grey:", (m["Module"].notna() & (m["Module"]!="grey")).sum(), "| purple:", (m["Module"]=="purple").sum())
