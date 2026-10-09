"""
Stage 5F – Parse GSE130387_Normalized_Intensity_File.xlsx (mouse Affymetrix)
Correct parsing: skip row 0 (metadata), use row 1 as header (ID_REF, Name, PD_1..Saline_3)
n=3 PD vs n=3 Saline mice → descriptive only, CROSS-SPECIES
mmu-miR-* → hsa-miR-* by stem homology
"""
import os, hashlib, datetime, re
import pandas as pd
import numpy as np

RAWDIR = r"d:\Peritoneal Project\data\raw"
OUTDIR = r"d:\Peritoneal Project\results\tables"

XLSX = os.path.join(RAWDIR, "GSE130387_Normalized_Intensity_File.xlsx")
sha256 = hashlib.sha256(open(XLSX, "rb").read()).hexdigest()
print(f"GSE130387 Excel SHA256: {sha256}")

# Read with header at row 1 (skip the normalization comment in row 0)
df = pd.read_excel(XLSX, sheet_name="Matrix", header=1)
print(f"Shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")
print(df.head(8).to_string())

# Column names: ID_REF, Name, PD_1, PD_2, PD_3, Saline_1, Saline_2, Saline_3
mirna_col  = "Name"
pd_cols    = ["PD_1", "PD_2", "PD_3"]
saline_cols = ["Saline_1", "Saline_2", "Saline_3"]

# Verify columns exist
for c in [mirna_col] + pd_cols + saline_cols:
    assert c in df.columns, f"Column '{c}' not found; available: {df.columns.tolist()}"
print(f"\nPD cols: {pd_cols}")
print(f"Saline cols: {saline_cols}")

# Drop rows with no Name (NaN in Name column)
df_clean = df[df[mirna_col].notna()].copy()
print(f"Rows with Name: {len(df_clean)} (of {len(df)} total)")

# The Affymetrix miRNA-4.0 array data here uses rno- (Rattus norvegicus) probe names
# despite the study description saying "mouse peritoneal tissue".
# GSE130387 platform GPL17107 contains probes from multiple rodent species.
# We map rno-* probes (rat) to hsa-* by stem homology (same CROSS-SPECIES caveat applies).
rno_mask  = df_clean[mirna_col].astype(str).str.startswith("rno-")
mmu_mask  = df_clean[mirna_col].astype(str).str.startswith("mmu-")
n_total   = len(df_clean)
n_mmu     = mmu_mask.sum()
n_rno     = rno_mask.sum()
print(f"\nTotal named probes: {n_total}")
print(f"mmu-* probes:       {n_mmu}")
print(f"rno-* probes:       {n_rno}")
print(f"NOTE: This dataset uses rno (rat) probes; study description says 'mouse tissue'.")
print(f"      Both are CROSS-SPECIES with respect to human.")

# Use rno probes (dominant species in this data)
df_mmu = df_clean[rno_mask].copy()
n_mmu  = n_rno  # rename for consistent downstream variable use

# rno → hsa name homology mapping
def mmu_to_hsa(name):
    """rno-miR-XXX → hsa-miR-XXX by prefix substitution (same as mmu-* logic)"""
    s = str(name)
    if s.startswith("rno-"):
        return "hsa-" + s[4:]
    if s.startswith("mmu-"):
        return "hsa-" + s[4:]
    return None

df_mmu["hsa_candidate"] = df_mmu[mirna_col].apply(mmu_to_hsa)
df_mmu["species_label"] = "CROSS-SPECIES"
df_mmu["mapping_method"] = "rno_to_hsa_stem_homology"

# Compute log2FC (already normalized; check scale)
expr_pd  = df_mmu[pd_cols].apply(pd.to_numeric, errors="coerce")
expr_sal = df_mmu[saline_cols].apply(pd.to_numeric, errors="coerce")

all_vals = pd.concat([expr_pd, expr_sal], axis=1).values.flatten()
all_vals = all_vals[~np.isnan(all_vals)]
print(f"\nExpression range: {all_vals.min():.4f} – {all_vals.max():.4f}")
print(f"Likely already log2 or normalized (max < 25): {all_vals.max() < 25}")

# log2(x+1) transformation before computing mean difference
log2_pd  = np.log2(expr_pd + 1)
log2_sal = np.log2(expr_sal + 1)

mean_pd  = log2_pd.mean(axis=1)
mean_sal = log2_sal.mean(axis=1)
log2fc   = mean_pd - mean_sal

out_df = pd.DataFrame({
    "Probe_ID_mmu":          df_mmu["ID_REF"].astype(str).values,
    "miRNA_mmu":             df_mmu[mirna_col].astype(str).values,
    "hsa_candidate":         df_mmu["hsa_candidate"].values,
    "species_label":         "CROSS-SPECIES",
    "mapping_method":        "mmu_to_hsa_stem_homology",
    "mean_log2_PD":          mean_pd.values,
    "mean_log2_Saline":      mean_sal.values,
    "log2FC_PDF_vs_saline":  log2fc.values,
    "dataset":               "GSE130387",
    "analysis_note":         "DESCRIPTIVE_ONLY_CROSS_SPECIES_n3vs3"
})
out_df = out_df.dropna(subset=["log2FC_PDF_vs_saline"])
out_df = out_df.sort_values("log2FC_PDF_vs_saline", ascending=False).reset_index(drop=True)

print(f"\nFinal: {len(out_df)} mmu probes with log2FC")
print(f"\nTop 20 UP in PD (CROSS-SPECIES):")
print(out_df.head(20)[["miRNA_mmu", "hsa_candidate", "mean_log2_PD", "mean_log2_Saline",
                         "log2FC_PDF_vs_saline"]].to_string(index=False))
print(f"\nTop 20 DOWN in PD (CROSS-SPECIES):")
print(out_df.tail(20)[["miRNA_mmu", "hsa_candidate", "mean_log2_PD", "mean_log2_Saline",
                         "log2FC_PDF_vs_saline"]].to_string(index=False))

out_df.to_csv(os.path.join(OUTDIR, "F_GSE130387_cross_species_logFC.csv"), index=False)
print(f"\nSaved F_GSE130387_cross_species_logFC.csv ({len(out_df)} probes)")

# miRBase mapping note
note = pd.DataFrame([{
    "mapping_method": "rno_to_hsa_stem_homology",
    "source_species": "Rattus norvegicus (rno); study description says mouse but probes are rno-*",
    "target_species": "Homo sapiens (hsa)",
    "platform": "GPL17107 (Affymetrix miRNA-4.0 / GeneChip miRNA 4.0)",
    "mirbase_version_note": "Probe names follow miRBase v21 convention; rno probes dominant in this dataset",
    "validation_status": "CROSS-SPECIES_NOT_VALIDATED",
    "n_rno_probes": n_mmu,
    "n_mapped_hsa": int(out_df["hsa_candidate"].notna().sum()),
    "note": "Functional equivalence of rno and hsa miRNAs NOT assumed; all rows labeled CROSS-SPECIES; species discrepancy noted"
}])
note.to_csv(os.path.join(OUTDIR, "F_GSE130387_mirbase_mapping_note.csv"), index=False)

# Sample metadata
smeta = pd.DataFrame(
    [{"Column": c, "Group": "PDF_mouse", "Species": "Mus musculus",
      "Species_label": "CROSS-SPECIES", "Dataset": "GSE130387"} for c in pd_cols] +
    [{"Column": c, "Group": "saline_mouse", "Species": "Mus musculus",
      "Species_label": "CROSS-SPECIES", "Dataset": "GSE130387"} for c in saline_cols]
)
smeta.to_csv(os.path.join(OUTDIR, "F_GSE130387_sample_metadata.csv"), index=False)

# Provenance
prov = pd.DataFrame([{
    "accession":      "GSE130387",
    "file":           "GSE130387_Normalized_Intensity_File.xlsx",
    "sha256":         sha256,
    "download_date":  datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "platform":       "GPL17107 Affymetrix miRNA-4.0",
    "n_probes_total": n_total,
    "n_rno_probes":   n_mmu,
    "n_with_logFC":   len(out_df),
    "pd_cols":        ";".join(pd_cols),
    "sal_cols":       ";".join(saline_cols),
    "analysis_type":  "DESCRIPTIVE_ONLY_CROSS_SPECIES_n3vs3",
    "species_note":   "Probe IDs are rno-* (rat); study title says mouse; both rodent species CROSS-SPECIES re human",
    "note":           "n=3 rodent PD vs saline; rno→hsa by stem name; all rows CROSS-SPECIES; no FDR"
}])
prov.to_csv(os.path.join(OUTDIR, "F_GSE130387_provenance.csv"), index=False)
print("Saved F_GSE130387_provenance.csv")
print("\nDone.")
