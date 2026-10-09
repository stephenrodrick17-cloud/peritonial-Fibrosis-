"""
Stage 5F – Step 4: Intersection analysis
Cross-reference:
  (A) Validated miRNAs targeting hub genes (from multiMiR)
  (B) DE miRNAs from GSE182736 (human effluent exosomes; all expressed miRNAs reported)
  (C) DE miRNAs from GSE130387 (rodent PD vs saline; CROSS-SPECIES; all expressed miRNAs)

For each dataset:
  - Universe = all miRNAs tested (all probes in that dataset)
  - Set A = validated miRNAs targeting ≥1 hub gene
  - Set B/C = DE miRNAs in that dataset (|log2FC| > threshold)
  - Overlap = A ∩ B or A ∩ C
  - Hypergeometric test for enrichment
  - Report universe sizes, overlap, P, and full lists

NOTE: For GSE182736 descriptive only (n=3), we use |log2FC| > 1 as "notable" threshold,
      not an FDR threshold. Similarly for GSE130387 (CROSS-SPECIES).
      Results are DESCRIPTIVE and HYPOTHESIS-GENERATING only.
"""

import os, datetime
import pandas as pd
import numpy as np
from scipy.stats import hypergeom

OUTDIR = r"d:\Peritoneal Project\results\tables"

# ------------------------------------------------------------------ #
# 1. Load validated hub miRNA universe (from multiMiR)
# ------------------------------------------------------------------ #
val_hub = pd.read_csv(os.path.join(OUTDIR, "F_hub_to_mirna_validated.csv"))
val_hub_mt = pd.read_csv(os.path.join(OUTDIR, "F_hub_to_mirna_mirtarbase.csv"))

# All unique validated miRNAs targeting ANY hub gene (validated all databases)
validated_mirnas_all = set(val_hub["miRNA"].dropna().unique())
# miRTarBase only
validated_mirnas_mtb = set(val_hub_mt["miRNA"].dropna().unique())

print(f"Validated miRNAs targeting >= 1 hub (all databases): {len(validated_mirnas_all)}")
print(f"Validated miRNAs targeting >= 1 hub (miRTarBase only): {len(validated_mirnas_mtb)}")

# ------------------------------------------------------------------ #
# 2. Load GSE182736 (human, n=3 vs 3)
# ------------------------------------------------------------------ #
gse182736 = pd.read_csv(os.path.join(OUTDIR, "F_GSE182736_descriptive_logFC.csv"))
print(f"\nGSE182736: {len(gse182736)} miRNAs total")
print(f"  Columns: {gse182736.columns.tolist()}")

# Clean miRNA names (strip isoform suffixes like _R-1, _L+1, etc.)
def clean_mirna_name(name):
    """Extract canonical miRNA name (strip sequence variant suffixes)"""
    import re
    s = str(name).strip()
    # Remove variant notations: _R-1, _L+1R-2, _1ss5TC, etc.
    s = re.sub(r'_[LR][+\-]\d+.*$', '', s)
    s = re.sub(r'_\d+ss.*$', '', s)
    return s

gse182736["miRNA_clean"] = gse182736["miRNA"].apply(clean_mirna_name)

# Universe: all miRNAs tested in GSE182736
universe_182736 = set(gse182736["miRNA_clean"].unique())
print(f"  Universe (unique clean names): {len(universe_182736)}")

# DE miRNAs at various thresholds (descriptive — no FDR available)
for thresh in [1.0, 0.5]:
    de_set = set(gse182736.loc[gse182736["log2FC_UF_vs_nonUF"].abs() > thresh,
                               "miRNA_clean"].unique())
    print(f"  |log2FC| > {thresh}: {len(de_set)} miRNAs")

# Use |log2FC| > 1 as "notable" threshold
de_182736 = set(gse182736.loc[gse182736["log2FC_UF_vs_nonUF"].abs() > 1.0,
                               "miRNA_clean"].unique())
print(f"  Selected threshold |log2FC| > 1.0: {len(de_182736)} miRNAs")

# Also: all expressed (any value, even borderline) as most inclusive universe
universe_182736_all = set(gse182736["miRNA_clean"].unique())

# ------------------------------------------------------------------ #
# 3. Load GSE130387 (rodent, CROSS-SPECIES)
# ------------------------------------------------------------------ #
gse130387 = pd.read_csv(os.path.join(OUTDIR, "F_GSE130387_cross_species_logFC.csv"))
print(f"\nGSE130387: {len(gse130387)} probes total")

# Use hsa_candidate as the miRNA ID for intersection
gse130387_named = gse130387.dropna(subset=["hsa_candidate"]).copy()
gse130387_named["hsa_clean"] = gse130387_named["hsa_candidate"].apply(clean_mirna_name)

universe_130387 = set(gse130387_named["hsa_clean"].unique())
print(f"  Universe (hsa candidates with names): {len(universe_130387)}")

de_130387 = set(gse130387_named.loc[
    gse130387_named["log2FC_PDF_vs_saline"].abs() > 0.5, "hsa_clean"].unique())
print(f"  |log2FC| > 0.5 (rodent; CROSS-SPECIES): {len(de_130387)} miRNAs")

# ------------------------------------------------------------------ #
# 4. Hypergeometric test utility
# ------------------------------------------------------------------ #
def hypergeom_test(N, K, n, k):
    """
    N = universe size
    K = # successes in universe (validated hub-targeting miRNAs in universe)
    n = # draws (DE miRNAs)
    k = # observed successes (overlap)
    Returns P(X >= k)
    """
    if k == 0:
        return 1.0
    # P(X >= k) = 1 - P(X < k) = 1 - CDF(k-1)
    p = hypergeom.sf(k - 1, N, K, n)
    return float(p)

# ------------------------------------------------------------------ #
# 5. Intersection: Validated hub miRNAs intersect GSE182736
# ------------------------------------------------------------------ #
print("\n" + "="*60)
print("INTERSECTION: validated hub miRNAs x GSE182736 (human)")
print("="*60)
print("Note: n=3 vs 3; |log2FC| > 1 threshold is descriptive only")
print("      No statistical inference on DE miRNAs themselves")

# Restrict validated miRNAs to those present in the universe
val_in_universe_182736 = validated_mirnas_all & universe_182736_all
print(f"\nValidated hub miRNAs also in GSE182736 universe: {len(val_in_universe_182736)}")
print(f"DE miRNAs (|log2FC|>1) in GSE182736:            {len(de_182736)}")

overlap_182736 = val_in_universe_182736 & de_182736
print(f"Overlap (validated hub x DE):                    {len(overlap_182736)}")

N_182736 = len(universe_182736_all)
K_182736 = len(val_in_universe_182736)
n_182736 = len(de_182736)
k_182736 = len(overlap_182736)
p_182736 = hypergeom_test(N_182736, K_182736, n_182736, k_182736)

print(f"\nHypergeometric test:")
print(f"  Universe N = {N_182736}")
print(f"  Validated in universe K = {K_182736}")
print(f"  DE miRNAs drawn n = {n_182736}")
print(f"  Overlap k = {k_182736}")
print(f"  P(X >= {k_182736}) = {p_182736:.6f}")
print(f"\nOverlapping miRNAs:")
if overlap_182736:
    for m in sorted(overlap_182736):
        # Get direction
        fc = gse182736.loc[gse182736["miRNA_clean"] == m, "log2FC_UF_vs_nonUF"].mean()
        hubs = val_hub.loc[val_hub["miRNA"] == m, "Hub_Gene"].tolist()
        print(f"  {m:30s}  log2FC={fc:+.3f}  Hubs: {hubs}")
else:
    print("  (empty set)")

# ------------------------------------------------------------------ #
# 6. Intersection: Validated hub miRNAs intersect GSE130387 (CROSS-SPECIES)
# ------------------------------------------------------------------ #
print("\n" + "="*60)
print("INTERSECTION: validated hub miRNAs x GSE130387 (CROSS-SPECIES)")
print("="*60)
print("Note: rno -> hsa by stem homology; all inferences CROSS-SPECIES")

val_in_universe_130387 = validated_mirnas_all & universe_130387
print(f"\nValidated hub miRNAs in GSE130387 universe: {len(val_in_universe_130387)}")
print(f"DE miRNAs (|log2FC|>0.5 rodent):           {len(de_130387)}")

overlap_130387 = val_in_universe_130387 & de_130387
print(f"Overlap (validated hub x DE):              {len(overlap_130387)}")

N_130387 = len(universe_130387)
K_130387 = len(val_in_universe_130387)
n_130387 = len(de_130387)
k_130387 = len(overlap_130387)
p_130387 = hypergeom_test(N_130387, K_130387, n_130387, k_130387)

print(f"\nHypergeometric test (CROSS-SPECIES):")
print(f"  Universe N = {N_130387}")
print(f"  Validated in universe K = {K_130387}")
print(f"  DE miRNAs drawn n = {n_130387}")
print(f"  Overlap k = {k_130387}")
print(f"  P(X >= {k_130387}) = {p_130387:.6f}")
print(f"\nOverlapping miRNAs (CROSS-SPECIES):")
if overlap_130387:
    for m in sorted(overlap_130387):
        fc = gse130387_named.loc[gse130387_named["hsa_clean"] == m, "log2FC_PDF_vs_saline"].mean()
        hubs = val_hub.loc[val_hub["miRNA"] == m, "Hub_Gene"].tolist()
        print(f"  {m:30s}  log2FC={fc:+.3f}  Hubs: {hubs}  [CROSS-SPECIES]")
else:
    print("  (empty set)")

# ------------------------------------------------------------------ #
# 7. Save intersection results
# ------------------------------------------------------------------ #
rows = []
# GSE182736
for m in sorted(overlap_182736):
    fc = gse182736.loc[gse182736["miRNA_clean"] == m, "log2FC_UF_vs_nonUF"].mean()
    hubs_val = val_hub.loc[val_hub["miRNA"] == m]
    for _, hr in hubs_val.iterrows():
        rows.append({
            "miRNA": m,
            "Hub_Gene": hr["Hub_Gene"],
            "Source": hr["Source"],
            "Dataset": "GSE182736",
            "species": "Homo sapiens",
            "species_label": "HUMAN",
            "log2FC": round(fc, 4),
            "DE_threshold": "|log2FC| > 1.0",
            "analysis_note": "DESCRIPTIVE_ONLY_n3vs3"
        })

# GSE130387
for m in sorted(overlap_130387):
    fc = gse130387_named.loc[gse130387_named["hsa_clean"] == m, "log2FC_PDF_vs_saline"].mean()
    hubs_val = val_hub.loc[val_hub["miRNA"] == m]
    for _, hr in hubs_val.iterrows():
        rows.append({
            "miRNA": m,
            "Hub_Gene": hr["Hub_Gene"],
            "Source": hr["Source"],
            "Dataset": "GSE130387",
            "species": "Rattus norvegicus -> Homo sapiens (hsa candidate)",
            "species_label": "CROSS-SPECIES",
            "log2FC": round(fc, 4),
            "DE_threshold": "|log2FC| > 0.5",
            "analysis_note": "DESCRIPTIVE_ONLY_CROSS_SPECIES_n3vs3"
        })

overlap_df = pd.DataFrame(rows)
if len(overlap_df) > 0:
    overlap_df.to_csv(os.path.join(OUTDIR, "F_mirna_hub_intersection.csv"), index=False)
    print(f"\nSaved F_mirna_hub_intersection.csv ({len(overlap_df)} rows)")
    print(overlap_df.to_string(index=False))
else:
    # Save empty with headers
    overlap_df = pd.DataFrame(columns=["miRNA", "Hub_Gene", "Source", "Dataset",
                                         "species", "species_label", "log2FC",
                                         "DE_threshold", "analysis_note"])
    overlap_df.to_csv(os.path.join(OUTDIR, "F_mirna_hub_intersection.csv"), index=False)
    print("\nSaved F_mirna_hub_intersection.csv (0 rows - empty intersection)")

# Summary statistics table
summary = pd.DataFrame([
    {
        "Dataset":           "GSE182736",
        "Species":           "Homo sapiens",
        "Species_label":     "HUMAN",
        "N_universe":        N_182736,
        "K_validated_in_universe": K_182736,
        "n_DE":              n_182736,
        "k_overlap":         k_182736,
        "hypergeom_P":       round(p_182736, 6),
        "DE_threshold":      "|log2FC| > 1.0",
        "analysis_note":     "DESCRIPTIVE_ONLY_n3vs3"
    },
    {
        "Dataset":           "GSE130387",
        "Species":           "Rattus norvegicus (rno->hsa)",
        "Species_label":     "CROSS-SPECIES",
        "N_universe":        N_130387,
        "K_validated_in_universe": K_130387,
        "n_DE":              n_130387,
        "k_overlap":         k_130387,
        "hypergeom_P":       round(p_130387, 6),
        "DE_threshold":      "|log2FC| > 0.5",
        "analysis_note":     "DESCRIPTIVE_ONLY_CROSS_SPECIES_n3vs3"
    }
])
summary.to_csv(os.path.join(OUTDIR, "F_mirna_intersection_summary.csv"), index=False)
print("\nSaved F_mirna_intersection_summary.csv")
print(summary.to_string(index=False))

print("\n=== Stage 5F intersection analysis complete ===")
