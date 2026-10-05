"""
Stage 5G - Step 1: Detect lncRNAs in GSE248762 by GENCODE biotype
Record version: GENCODE v32 (GRCh38.p13, Ensembl 98, 2019-09-18)
Matching 10x Genomics Cell Ranger GRCh38 reference used for GSE248762.

Outputs:
  - G_lncRNA_GSE248762_detection_summary.csv
  - G_lncRNA_GSE248762_by_celltype.csv
  - G_lncRNA_detected_list.csv
  - G_lncRNA_provenance.csv
"""

import os, gzip, re, datetime
import pandas as pd
import numpy as np

DATA_DIR = r"d:\Peritoneal Project\data"
OUTDIR = r"d:\Peritoneal Project\results\tables"
os.makedirs(OUTDIR, exist_ok=True)

# 1. Parse GENCODE v32 lncRNA GTF
gtf_path = os.path.join(DATA_DIR, "raw", "gencode.v32.long_noncoding_RNAs.gtf.gz")
print(f"Loading GENCODE v32 lncRNA GTF from: {gtf_path}")

lnc_dict = {}  # Ensembl_ID_clean -> (gene_name, gene_type, full_id)
with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) > 8 and fields[2] == "gene":
            gid_m = re.search(r'gene_id "([^"]+)"', fields[8])
            gn_m = re.search(r'gene_name "([^"]+)"', fields[8])
            gt_m = re.search(r'gene_type "([^"]+)"', fields[8])
            if gid_m and gn_m:
                full_gid = gid_m.group(1)
                clean_gid = full_gid.split(".")[0]
                gene_name = gn_m.group(1)
                gene_type = gt_m.group(1) if gt_m else "lncRNA"
                lnc_dict[clean_gid] = {
                    "gene_id": full_gid,
                    "clean_id": clean_gid,
                    "gene_name": gene_name,
                    "gene_type": gene_type
                }

print(f"Total lncRNA genes annotated in GENCODE v32: {len(lnc_dict)}")

# 2. Parse GSE248762 features
features_path = os.path.join(DATA_DIR, "raw", "GSE248762_extracted", "GSM7919583_LV_UF-1.features.tsv.gz")
gse_features = []
with gzip.open(features_path, "rt", encoding="utf-8") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            gse_features.append({"ensembl_id": parts[0], "symbol": parts[1], "type": parts[2] if len(parts) > 2 else ""})

df_feat = pd.DataFrame(gse_features)
print(f"Total features in GSE248762 reference: {len(df_feat)}")

# Match features against GENCODE v32 lncRNAs
df_feat["is_lncRNA"] = df_feat["ensembl_id"].isin(lnc_dict)
df_lnc_feat = df_feat[df_feat["is_lncRNA"]].copy()
df_lnc_feat["gencode_name"] = df_lnc_feat["ensembl_id"].apply(lambda x: lnc_dict[x]["gene_name"])
df_lnc_feat["gencode_type"] = df_lnc_feat["ensembl_id"].apply(lambda x: lnc_dict[x]["gene_type"])
print(f"Features annotated as lncRNA in GSE248762: {len(df_lnc_feat)}")

# 3. Read pseudobulk count matrices across cell types
cell_types = {
    "Stromal": "pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv",
    "Monocyte_macrophage": "pb_primary_Monocyte_macrophage_counts.csv",
    "T_cell": "pb_primary_T_cell_counts.csv",
    "NK_cell": "pb_primary_NK_cell_counts.csv",
    "cDC": "pb_primary_cDC_counts.csv"
}

pb_dir = os.path.join(DATA_DIR, "processed", "pseudobulk")
celltype_dfs = {}
for ct_name, fname in cell_types.items():
    fpath = os.path.join(pb_dir, fname)
    if os.path.exists(fpath):
        celltype_dfs[ct_name] = pd.read_csv(fpath, index_col=0)
        print(f"  Loaded {ct_name}: shape={celltype_dfs[ct_name].shape}")

# In pseudobulk matrices, columns are symbols (df_feat["symbol"])
# Map symbols to Ensembl ID and check lncRNAs
# If duplicate symbols exist in features, match by symbol in df_lnc_feat
lnc_symbols = set(df_lnc_feat["symbol"].unique())
print(f"Unique lncRNA symbols in GSE248762: {len(lnc_symbols)}")

# 4. Measure expression and detection per cell type
ct_summary = []
detected_lncrnas_per_ct = {}

for ct_name, pb_df in celltype_dfs.items():
    # Find columns present in lnc_symbols
    common_cols = [c for c in pb_df.columns if c in lnc_symbols]
    sub_counts = pb_df[common_cols]
    
    # Counts sum per gene across libraries
    gene_sums = sub_counts.sum(axis=0)
    # Number of libraries where count > 0
    gene_n_libs = (sub_counts > 0).sum(axis=0)
    gene_means = sub_counts.mean(axis=0)
    
    n_detected_gt0 = (gene_sums > 0).sum()
    n_detected_ge2 = (gene_n_libs >= 2).sum()
    n_detected_ge5 = (gene_n_libs >= 5).sum()
    
    ct_summary.append({
        "Cell_Type": ct_name,
        "Total_lncRNAs_assayed": len(common_cols),
        "Detected_gt0_libraries": int(n_detected_gt0),
        "Detected_ge2_libraries": int(n_detected_ge2),
        "Detected_ge5_libraries": int(n_detected_ge5),
        "Total_counts_all_lncRNAs": float(gene_sums.sum()),
        "Mean_counts_per_detected": float(gene_sums[gene_sums > 0].mean())
    })
    
    detected_lncrnas_per_ct[ct_name] = set(gene_sums[gene_sums > 0].index)

df_ct_summary = pd.DataFrame(ct_summary)
print("\nlmcRNA Detection by Cell Type Summary:")
print(df_ct_summary.to_string(index=False))

# 5. Across all cell types combined (whole GSE248762 dataset)
all_detected_gt0 = set()
all_detected_ge2 = set()
for s in detected_lncrnas_per_ct.values():
    all_detected_gt0 |= s

print(f"\nTotal unique lncRNAs detected (counts > 0 in >= 1 cell type): {len(all_detected_gt0)}")

# Build detailed per-lncRNA expression table
records = []
for sym in common_cols:
    row = {
        "symbol": sym,
        "is_detected_any": sym in all_detected_gt0
    }
    total_count = 0
    total_libs = 0
    for ct_name, pb_df in celltype_dfs.items():
        if sym in pb_df.columns:
            cnt_sum = float(pb_df[sym].sum())
            n_libs = int((pb_df[sym] > 0).sum())
            row[f"{ct_name}_total_counts"] = cnt_sum
            row[f"{ct_name}_n_libs"] = n_libs
            total_count += cnt_sum
            total_libs += n_libs
        else:
            row[f"{ct_name}_total_counts"] = 0.0
            row[f"{ct_name}_n_libs"] = 0
    
    row["Overall_total_counts"] = total_count
    row["Overall_n_libs"] = total_libs
    records.append(row)

df_all_lnc = pd.DataFrame(records)
# Merge annotation info
feat_sub = df_lnc_feat.drop_duplicates(subset=["symbol"])[["symbol", "ensembl_id", "gencode_name", "gencode_type"]]
df_all_lnc = df_all_lnc.merge(feat_sub, on="symbol", how="left")

# Reorder columns
cols_order = ["symbol", "ensembl_id", "gencode_name", "gencode_type", "is_detected_any", 
              "Overall_total_counts", "Overall_n_libs",
              "Stromal_total_counts", "Stromal_n_libs",
              "Monocyte_macrophage_total_counts", "Monocyte_macrophage_n_libs",
              "T_cell_total_counts", "T_cell_n_libs",
              "NK_cell_total_counts", "NK_cell_n_libs",
              "cDC_total_counts", "cDC_n_libs"]
df_all_lnc = df_all_lnc[[c for c in cols_order if c in df_all_lnc.columns]]

# Filter to detected lncRNAs
df_detected = df_all_lnc[df_all_lnc["is_detected_any"]].sort_values(by="Overall_total_counts", ascending=False)
print(f"Top 10 highest expressed lncRNAs overall:")
print(df_detected[["symbol", "ensembl_id", "Overall_total_counts", "Stromal_total_counts", "Monocyte_macrophage_total_counts"]].head(10).to_string(index=False))

# 6. Save tables
summary_path = os.path.join(OUTDIR, "G_lncRNA_GSE248762_detection_summary.csv")
df_ct_summary.to_csv(summary_path, index=False)
print(f"\nSaved {summary_path}")

celltype_path = os.path.join(OUTDIR, "G_lncRNA_GSE248762_by_celltype.csv")
df_all_lnc.to_csv(celltype_path, index=False)
print(f"Saved {celltype_path} ({len(df_all_lnc)} lncRNAs assayed)")

detected_path = os.path.join(OUTDIR, "G_lncRNA_detected_list.csv")
df_detected.to_csv(detected_path, index=False)
print(f"Saved {detected_path} ({len(df_detected)} lncRNAs detected)")

# Provenance file
prov = pd.DataFrame([{
    "Dataset": "GSE248762",
    "Annotation_Source": "GENCODE",
    "GENCODE_Version": "v32",
    "Genome_Build": "GRCh38.p13",
    "Ensembl_Version": "Ensembl 98",
    "Reference_Release_Date": "2019-09-18",
    "CellRanger_Reference": "10x Genomics GRCh38-2020-A (built on GENCODE v32)",
    "Total_Annotated_lncRNA_in_GENCODE_v32": len(lnc_dict),
    "lncRNA_in_GSE248762_Features": len(df_lnc_feat),
    "lncRNA_Detected_gt0": len(df_detected),
    "Analysis_Date": datetime.date.today().isoformat(),
    "Notes": "Detection evaluated on pseudobulk matrices across 5 cell types (Stromal, Monocyte/macrophage, T cell, NK cell, cDC)"
}])
prov_path = os.path.join(OUTDIR, "G_lncRNA_provenance.csv")
prov.to_csv(prov_path, index=False)
print(f"Saved {prov_path}")

print("\n=== Stage 5G Step 1 complete ===")
