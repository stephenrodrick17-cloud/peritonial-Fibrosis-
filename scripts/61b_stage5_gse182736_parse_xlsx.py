"""
Stage 5F – Parse GSE182736_All_Expressed_miRNA.xlsx with correct header row (row 51)
Groups: disease_*(norm) = UF_failure, control_*(norm) = non_UF (stable PD function)
Descriptive only; n=3 vs 3; no FDR claims.
"""
import os, hashlib, datetime
import pandas as pd
import numpy as np

RAWDIR = r"d:\Peritoneal Project\data\raw"
OUTDIR = r"d:\Peritoneal Project\results\tables"

XLSX = os.path.join(RAWDIR, "GSE182736_All_Expressed_miRNA.xlsx")
sha256 = hashlib.sha256(open(XLSX, "rb").read()).hexdigest()
print(f"SHA256: {sha256}")

# Read with header at row 51 (0-indexed)
df = pd.read_excel(XLSX, sheet_name="Sheet1", header=51)
print(f"Shape after correct header: {df.shape}")
print(f"Columns: {df.columns.tolist()}")
print(df.head(5).to_string())

# miRNA name column
mirna_col = "miR_name"

# Normalized expression columns
ctrl_cols = [c for c in df.columns if "control" in str(c).lower() and "norm" in str(c).lower()]
dis_cols  = [c for c in df.columns if "disease" in str(c).lower() and "norm" in str(c).lower()]

print(f"\nControl (non-UF) norm cols: {ctrl_cols}")
print(f"Disease (UF_failure) norm cols: {dis_cols}")
print(f"n_control={len(ctrl_cols)}, n_disease={len(dis_cols)}")

# Filter to hsa-* miRNAs only (exclude other species)
hsa_mask = df[mirna_col].astype(str).str.startswith("hsa-")
df_hsa = df[hsa_mask].copy()
print(f"\nhsa-* miRNAs: {len(df_hsa)}")

# Compute log2 fold change (disease/UF_failure vs control/non-UF)
expr_ctrl = df_hsa[ctrl_cols].apply(pd.to_numeric, errors="coerce")
expr_dis  = df_hsa[dis_cols].apply(pd.to_numeric, errors="coerce")

# Check value range for log-scaling decision
all_vals = pd.concat([expr_ctrl, expr_dis], axis=1).values.flatten()
all_vals = all_vals[~np.isnan(all_vals) & (all_vals > 0)]
print(f"\nNorm expression range: {all_vals.min():.2f} – {all_vals.max():.2f}")
print(f"Are values already log2? {all_vals.max() < 25}")

# Use log2(norm+1) for FC computation
mean_dis  = np.log2(expr_dis.add(1).mean(axis=1))
mean_ctrl = np.log2(expr_ctrl.add(1).mean(axis=1))
log2fc = mean_dis - mean_ctrl

out_df = pd.DataFrame({
    "miRNA":               df_hsa[mirna_col].astype(str).values,
    "mean_norm_UF_failure": expr_dis.mean(axis=1).values,
    "mean_norm_nonUF":      expr_ctrl.mean(axis=1).values,
    "log2FC_UF_vs_nonUF":  log2fc.values,
    "dataset":             "GSE182736",
    "analysis_note":       "DESCRIPTIVE_ONLY_n3vs3",
    "species":             "Homo sapiens",
})

# Add original expression columns
for c in ctrl_cols + dis_cols:
    out_df[c] = df_hsa[c].values

out_df = out_df.dropna(subset=["log2FC_UF_vs_nonUF"])
out_df = out_df.sort_values("log2FC_UF_vs_nonUF", ascending=False).reset_index(drop=True)

print(f"\nFinal table: {len(out_df)} hsa-miRNAs with expression data")
print(f"\nTop 20 UP in UF_failure:")
print(out_df.head(20)[["miRNA", "mean_norm_UF_failure", "mean_norm_nonUF", "log2FC_UF_vs_nonUF"]].to_string(index=False))
print(f"\nTop 20 DOWN in UF_failure:")
print(out_df.tail(20)[["miRNA", "mean_norm_UF_failure", "mean_norm_nonUF", "log2FC_UF_vs_nonUF"]].to_string(index=False))

out_df.to_csv(os.path.join(OUTDIR, "F_GSE182736_descriptive_logFC.csv"), index=False)
print(f"\nSaved F_GSE182736_descriptive_logFC.csv ({len(out_df)} miRNAs)")

# Sample metadata
smeta = pd.DataFrame(
    [{"Column": c, "Group": "non_UF", "Dataset": "GSE182736", "Type": "normalized"} for c in ctrl_cols] +
    [{"Column": c, "Group": "UF_failure", "Dataset": "GSE182736", "Type": "normalized"} for c in dis_cols]
)
smeta.to_csv(os.path.join(OUTDIR, "F_GSE182736_sample_metadata.csv"), index=False)
print("Saved F_GSE182736_sample_metadata.csv")

# Provenance
prov = pd.DataFrame([{
    "accession":        "GSE182736",
    "file":             "GSE182736_All_Expressed_miRNA.xlsx",
    "sha256":           sha256,
    "download_date":    datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "header_row":       51,
    "n_mirnas_total":   len(df),
    "n_hsa_mirnas":     len(df_hsa),
    "n_with_logFC":     len(out_df),
    "ctrl_cols":        ";".join(ctrl_cols),
    "dis_cols":         ";".join(dis_cols),
    "analysis_type":    "DESCRIPTIVE_ONLY_n3vs3",
    "groups":           "control=non_UF(stable PD function); disease=UF_failure(type I)",
    "note":             "n=3 vs n=3 human peritoneal effluent exosomes; no statistical inference; log2FC descriptive only"
}])
prov.to_csv(os.path.join(OUTDIR, "F_GSE182736_provenance.csv"), index=False)
print("Saved F_GSE182736_provenance.csv")
print("\nDone.")
