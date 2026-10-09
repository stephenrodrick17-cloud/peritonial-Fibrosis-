"""
Stage 5H: GSE121372 Analysis
Human peritoneal mesothelial cells (HPMC) stimulated with TGF-b1 (1 ng/mL)
Timepoints: 6 hours and 24 hours.
Design: n=1 per condition (unreplicated).
Status: DESCRIPTIVE ONLY; no statistical testing or P-values.

Platform: GPL6255 (Illumina HumanRef-8 v2.0 Expression BeadChip)
Target: 11 hub genes (ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX)

Outputs:
  - H_GSE121372_hub_fold_changes.csv
  - H_GSE121372_sample_metadata.csv
  - H_GSE121372_provenance.csv
"""

import os, gzip, datetime
import pandas as pd
import numpy as np

DATA_DIR = r"d:\Peritoneal Project\data\raw"
OUTDIR = r"d:\Peritoneal Project\results\tables"
os.makedirs(OUTDIR, exist_ok=True)

# 1. Platform probe mapping (GPL6255)
annot_path = os.path.join(DATA_DIR, "GPL6255.annot.gz")
print(f"Reading GPL6255 annotation from: {annot_path}")

hub_genes = ["ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"]

probe_to_symbol = {}
symbol_to_probes = {h: [] for h in hub_genes}

with gzip.open(annot_path, "rt", encoding="utf-8", errors="ignore") as f:
    in_table = False
    for line in f:
        if line.startswith("!platform_table_begin"):
            in_table = True
            header = f.readline().strip().split("\t")
            id_idx = header.index("ID")
            sym_idx = header.index("Gene symbol")
            continue
        if not in_table or line.startswith("!"):
            continue
        parts = line.strip().split("\t")
        if len(parts) > max(id_idx, sym_idx):
            pid = parts[id_idx].strip()
            sym = parts[sym_idx].strip()
            if sym in hub_genes:
                symbol_to_probes[sym].append(pid)
                probe_to_symbol[pid] = sym

print("\nHub gene probe mapping on GPL6255:")
for h in hub_genes:
    print(f"  {h:10s}: {symbol_to_probes[h]}")

# 2. Load GSE121372 normalized data
norm_path = os.path.join(DATA_DIR, "GSE121372_normalized.txt.gz")
print(f"\nLoading normalized data from: {norm_path}")
df_norm = pd.read_csv(norm_path, sep="\t", compression="gzip")
print(f"  Matrix shape: {df_norm.shape}")
print(f"  Columns: {df_norm.columns.tolist()}")

# Normalized expression columns
# HPMC_control_6hr, HPMC_TGF-b1_6hr, HPMC_control_24hr, HPMC_TGF-b1_24hr
# Detection Pval columns exist as quality checks

# 3. Calculate fold changes for all 11 hub genes
results = []
pseudocount = 1.0  # standard for linear microarray intensity

for h in hub_genes:
    probes = symbol_to_probes[h]
    if len(probes) == 0:
        # Gene not assayed on platform (e.g. ISM1)
        results.append({
            "Gene_Symbol": h,
            "Probe_ID": "NOT_ON_PLATFORM",
            "Platform": "GPL6255",
            "control_6hr": np.nan,
            "TGFb1_6hr": np.nan,
            "log2FC_6hr": np.nan,
            "direction_6hr": "NOT_ASSAYED",
            "control_24hr": np.nan,
            "TGFb1_24hr": np.nan,
            "log2FC_24hr": np.nan,
            "direction_24hr": "NOT_ASSAYED",
            "det_pval_ctrl_6h": np.nan,
            "det_pval_tgf_6h": np.nan,
            "det_pval_ctrl_24h": np.nan,
            "det_pval_tgf_24h": np.nan,
            "replicate_design": "n=1 per condition",
            "analysis_note": "NOT_ASSAYED_ON_PLATFORM; gene absent from Illumina HumanRef-8 v2.0"
        })
        continue
    
    for pid in probes:
        sub = df_norm[df_norm["ID_REF"] == pid]
        if len(sub) == 0:
            continue
        row = sub.iloc[0]
        
        c6 = float(row["HPMC_control_6hr"])
        t6 = float(row["HPMC_TGF-b1_6hr"])
        c24 = float(row["HPMC_control_24hr"])
        t24 = float(row["HPMC_TGF-b1_24hr"])
        
        dp_c6 = float(row["Detection Pval"])
        dp_t6 = float(row["Detection Pval.1"])
        dp_c24 = float(row["Detection Pval.2"])
        dp_t24 = float(row["Detection Pval.3"])
        
        # log2 fold change: log2( (TGF + pseudo) / (Ctrl + pseudo) )
        # Values can be negative in background-subtracted data, so clip at 0 before pseudocount
        c6_val = max(c6, 0.0) + pseudocount
        t6_val = max(t6, 0.0) + pseudocount
        c24_val = max(c24, 0.0) + pseudocount
        t24_val = max(t24, 0.0) + pseudocount
        
        log2fc_6h = round(np.log2(t6_val / c6_val), 4)
        log2fc_24h = round(np.log2(t24_val / c24_val), 4)
        
        dir_6h = "UP" if log2fc_6h > 0.5 else ("DOWN" if log2fc_6h < -0.5 else "UNCHANGED")
        dir_24h = "UP" if log2fc_24h > 0.5 else ("DOWN" if log2fc_24h < -0.5 else "UNCHANGED")
        
        results.append({
            "Gene_Symbol": h,
            "Probe_ID": pid,
            "Platform": "GPL6255",
            "control_6hr": round(c6, 2),
            "TGFb1_6hr": round(t6, 2),
            "log2FC_6hr": log2fc_6h,
            "direction_6hr": dir_6h,
            "control_24hr": round(c24, 2),
            "TGFb1_24hr": round(t24, 2),
            "log2FC_24hr": log2fc_24h,
            "direction_24hr": dir_24h,
            "det_pval_ctrl_6h": round(dp_c6, 6),
            "det_pval_tgf_6h": round(dp_t6, 6),
            "det_pval_ctrl_24h": round(dp_c24, 6),
            "det_pval_tgf_24h": round(dp_t24, 6),
            "replicate_design": "n=1 per condition",
            "analysis_note": "DESCRIPTIVE_ONLY_n1; no statistical testing or P-values possible"
        })

df_res = pd.DataFrame(results)
print("\nGSE121372 Hub Gene Fold Changes:")
print(df_res[["Gene_Symbol", "Probe_ID", "control_6hr", "TGFb1_6hr", "log2FC_6hr", "control_24hr", "TGFb1_24hr", "log2FC_24hr"]].to_string(index=False))

# 4. Save results table
out_csv = os.path.join(OUTDIR, "H_GSE121372_hub_fold_changes.csv")
df_res.to_csv(out_csv, index=False)
print(f"\nSaved {out_csv} ({len(df_res)} rows)")

# 5. Sample metadata table
meta = pd.DataFrame([
    {
        "Sample_ID": "GSM3433436",
        "Sample_Title": "HPMC_control_6hr",
        "Cell_Type": "Human peritoneal mesothelial cells (HPMC)",
        "Treatment": "Control (none)",
        "Timepoint": "6 hours",
        "Concentration": "0 ng/mL",
        "Assay": "Illumina HumanRef-8 v2.0 Expression BeadChip",
        "Platform": "GPL6255",
        "Replicate": "n=1"
    },
    {
        "Sample_ID": "GSM3433437",
        "Sample_Title": "HPMC_TGF-b1_6hr",
        "Cell_Type": "Human peritoneal mesothelial cells (HPMC)",
        "Treatment": "TGF-b1",
        "Timepoint": "6 hours",
        "Concentration": "1 ng/mL",
        "Assay": "Illumina HumanRef-8 v2.0 Expression BeadChip",
        "Platform": "GPL6255",
        "Replicate": "n=1"
    },
    {
        "Sample_ID": "GSM3433438",
        "Sample_Title": "HPMC_control_24hr",
        "Cell_Type": "Human peritoneal mesothelial cells (HPMC)",
        "Treatment": "Control (none)",
        "Timepoint": "24 hours",
        "Concentration": "0 ng/mL",
        "Assay": "Illumina HumanRef-8 v2.0 Expression BeadChip",
        "Platform": "GPL6255",
        "Replicate": "n=1"
    },
    {
        "Sample_ID": "GSM3433439",
        "Sample_Title": "HPMC_TGF-b1_24hr",
        "Cell_Type": "Human peritoneal mesothelial cells (HPMC)",
        "Treatment": "TGF-b1",
        "Timepoint": "24 hours",
        "Concentration": "1 ng/mL",
        "Assay": "Illumina HumanRef-8 v2.0 Expression BeadChip",
        "Platform": "GPL6255",
        "Replicate": "n=1"
    }
])
meta_csv = os.path.join(OUTDIR, "H_GSE121372_sample_metadata.csv")
meta.to_csv(meta_csv, index=False)
print(f"Saved {meta_csv}")

# 6. Provenance table
prov = pd.DataFrame([{
    "Accession": "GSE121372",
    "Title": "Gene expression profile of human peritoneal mesothelial cells (HPMC) treated with TGF-b1",
    "Organism": "Homo sapiens",
    "Platform_ID": "GPL6255",
    "Platform_Title": "Illumina humanRef-8 v2.0 expression beadchip",
    "PubMed_ID": "30728376",
    "Total_Samples": 4,
    "Replicates_per_Condition": 1,
    "Normalized_Data_Source": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121372/suppl/GSE121372_normalized.txt.gz",
    "Analysis_Date": datetime.date.today().isoformat(),
    "Hub_Genes_Assayed": "10 of 11 (ISM1 absent from platform GPL6255)",
    "Caveat": "DESCRIPTIVE ONLY; unreplicated (n=1); no statistical testing or P-values possible; external validation of in vitro TGF-b1 mesothelial response"
}])
prov_csv = os.path.join(OUTDIR, "H_GSE121372_provenance.csv")
prov.to_csv(prov_csv, index=False)
print(f"Saved {prov_csv}")

print("\n=== Stage 5H analysis complete ===")
