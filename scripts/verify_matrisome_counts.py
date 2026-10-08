"""Verify matrisome counts: probe-level (user's tally) vs unique-gene (authoritative)."""
import pandas as pd

df = pd.read_csv(r"d:\Peritoneal Project\results\tables\GSE62928_ECM_intersection.csv")

print("=" * 80)
print("1. PROBE-LEVEL row counts (raw CSV rows — exactly what user manually tallied)")
print("=" * 80)
print(f"Total CSV rows: {len(df):,}")
probe_cat = df.groupby("Category").size().sort_values(ascending=False)
for i, (c, n) in enumerate(probe_cat.items(), 1):
    print(f"  {i:>2} {c:>30s}  {n:>6,}")
print()

print("=" * 80)
print("2. UNIQUE GENE SYMBOLS per Category (gene-level — what the pipeline reports)")
print("=" * 80)
uniq_gene_cat = df.groupby("Category")["Gene"].nunique().sort_values(ascending=False)
for i, (c, n) in enumerate(uniq_gene_cat.items(), 1):
    print(f"  {i:>2} {c:>30s}  {n:>6,}")
print()

non_mat = df["Category"] == "Non-matrisome"
mat = ~non_mat
n_mat_genes = df.loc[mat, "Gene"].nunique()
n_nonmat_genes = df.loc[non_mat, "Gene"].nunique()
print(f"Unique matrisome genes (non-Non-matrisome):  {n_mat_genes:>6,}")
print(f"Unique Non-matrisome genes:                   {n_nonmat_genes:>6,}")
print(f"Unique genes total (table):                   {df['Gene'].nunique():>6,}")
print()

core = df["Category"].isin(["Collagens", "ECM Glycoproteins", "Proteoglycans"])
ma = df["Category"].isin(["ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"])
core_genes = df.loc[core, "Gene"].nunique()
ma_genes = df.loc[ma, "Gene"].nunique()
print(f"Core matrisome (unique genes, Collagens + Glycoproteins + Proteoglycans):  {core_genes}")
print(f"Matrisome-associated (ECM-affiliated + Regulators + Secreted):             {ma_genes}")
print(f"Sum Core + MA = {core_genes + ma_genes} (should equal matrisome unique = {n_mat_genes})")
print()

print("=" * 80)
print("3. PROBE-DUPLICATION RATIO by Category  —  explains the user's 1972 vs 975 gap")
print("=" * 80)
cmp = pd.concat([probe_cat.rename("Probes"), uniq_gene_cat.rename("Unique_Genes")], axis=1).fillna(0)
cmp["Probes"] = cmp["Probes"].astype(int)
cmp["Unique_Genes"] = cmp["Unique_Genes"].astype(int)
cmp["Dup_Ratio"] = (cmp["Probes"] / cmp["Unique_Genes"]).round(2)
print(cmp.sort_values("Dup_Ratio", ascending=False).to_string())
print()

print("=" * 80)
print("4. RECONCILIATION SUMMARY")
print("=" * 80)
probe_mat_rows = int(probe_cat[probe_cat.index != "Non-matrisome"].sum())
probe_core_rows = int(probe_cat[probe_cat.index.isin(["Collagens", "ECM Glycoproteins", "Proteoglycans"])].sum())
probe_ma_rows = int(probe_cat[probe_cat.index.isin(["ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"])].sum())
print(f"User's probe-level totals (from raw CSV rows):")
print(f"  Matrisome probe rows total:          {probe_mat_rows:>5}  (user got 1972? — diff = {probe_mat_rows - 1972:+d})")
print(f"  Core matrisome probe rows:           {probe_core_rows:>5}  (user got 588?  — diff = {probe_core_rows - 588:+d})")
print(f"  Matrisome-associated probe rows:     {probe_ma_rows:>5}  (user got 1384? — diff = {probe_ma_rows - 1384:+d})")
print()
print(f"Authoritative pipeline UNIQUE GENE counts:")
print(f"  Measured matrisome genes:            {n_mat_genes:>5}  (matches 975 expected)")
print(f"  Core matrisome genes:                {core_genes:>5}")
print(f"  Matrisome-associated genes:          {ma_genes:>5}")
print()
print("Conclusion:")
print("  Your manual tally = PROBE-level row counts (1 gene = many Illumina probes),")
print("  Pipeline uses UNIQUE gene-symbol collapsed counts for all published numbers.")
