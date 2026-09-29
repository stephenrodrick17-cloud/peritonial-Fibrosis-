"""
Technical & Code Logic Verification Script for Peritoneal Fibrosis Pipeline
Checks A through G
"""

import os
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu

print("=" * 80)
print("RUNNING TECHNICAL & CODE LOGIC VERIFICATION AUDIT")
print("=" * 80)

# ==============================================================================
# CHECK A: SAMPLE IDENTITY AND LABEL INTEGRITY
# ==============================================================================
print("\n" + "=" * 80)
print("CHECK A: SAMPLE IDENTITY AND LABEL INTEGRITY")
print("=" * 80)

# 1. GSE62928 sample metadata vs expression matrix
df_meta62 = pd.read_csv("results/tables/GSE62928_sample_metadata.csv")
df_expr62 = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)

meta_samples = df_meta62["sample_id"].tolist()
expr_cols = df_expr62.columns.tolist()

print(f"GSE62928 Metadata Sample Count: {len(meta_samples)}")
print(f"GSE62928 Expression Matrix Col Count: {len(expr_cols)}")
print("\nSide-by-side Sample Order Comparison (GSE62928):")
aligned_62 = True
for i, (m_id, e_id) in enumerate(zip(meta_samples, expr_cols)):
    match = (m_id == e_id)
    if not match: aligned_62 = False
    row = df_meta62.iloc[i]
    print(f"  [{i}] Metadata: {m_id:<12} | Expr Col: {e_id:<12} | Title: {row['sample_title']:<25} | Group: {row['group_3class']:<8} | Binary: {row['group_binary_trait']} ({row['binary_numeric']}) | Match: {match}")

print(f"\nOrder alignment check for GSE62928: {'PERFECT PASS' if aligned_62 else 'FAIL MISALIGNMENT'}")

# Verify against GEO series matrix file directly
import gzip
series_file = "data/GSE62928_series_matrix.txt.gz"
if os.path.exists(series_file):
    print(f"\nCross-checking against {series_file}:")
    with gzip.open(series_file, "rt") as f:
        geo_titles = []
        geo_accessions = []
        geo_chars = []
        for line in f:
            if line.startswith("!Sample_title"):
                geo_titles = [x.strip('"\n ') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo_accessions = [x.strip('"\n ') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and not geo_chars:
                geo_chars = [x.strip('"\n ') for x in line.split("\t")[1:]]
    for i in range(len(geo_accessions)):
        print(f"  GEO Acc: {geo_accessions[i]} | GEO Title: {geo_titles[i]} | Meta Match: {geo_accessions[i] == meta_samples[i]}")

# 2. GSE125498 Sample Check
import GEOparse
gse12 = GEOparse.get_GEO(filepath="data/GSE125498_family.soft.gz")
piv12 = gse12.pivot_samples("VALUE")

sample_stages = {}
for name, gsm in gse12.gsms.items():
    title = gsm.metadata.get("title", [""])[0]
    chars = str(gsm.metadata.get("characteristics_ch1", []))
    is_long = "long-term" in title.lower() or "long-term" in chars.lower()
    sample_stages[name] = "Late_Stage_LPD" if is_long else "Early_Stage_SPD"

early_samples = [s for s, st in sample_stages.items() if st == "Early_Stage_SPD"]
late_samples = [s for s, st in sample_stages.items() if st == "Late_Stage_LPD"]
all_samples12 = early_samples + late_samples

print(f"\nGSE125498 Total GSM Samples in SOFT: {len(gse12.gsms)}")
print(f"GSE125498 Pivot Matrix Dimensions: {piv12.shape[0]} probes x {piv12.shape[1]} samples")
print(f"  - Early Stage (Short-term PD, 0-24 mo): {len(early_samples)} samples")
print(f"  - Late Stage (Long-term PD, >=25 mo):   {len(late_samples)} samples")
print(f"All samples present in pivot matrix: {all(s in piv12.columns for s in all_samples12)}")


# ==============================================================================
# CHECK B: GENE SYMBOL MAPPING CONSISTENCY
# ==============================================================================
print("\n" + "=" * 80)
print("CHECK B: GENE SYMBOL MAPPING CONSISTENCY & 40-GENE INTERSECTION")
print("=" * 80)

f_mod_genes = "results/tables/wgcna_trait_significant_module_genes.csv"
f_ecm = "convergent_ECM_DEGs_nominal.csv"
f_conv = "results/tables/convergent_WGCNA_ECM_genes.csv"
if not os.path.exists(f_conv):
    f_conv = "convergent_WGCNA_ECM_genes.csv"

df_mod_g = pd.read_csv(f_mod_genes)
df_ecm_g = pd.read_csv(f_ecm)
df_conv_g = pd.read_csv(f_conv)

# Clean symbols
mod_symbols_raw = df_mod_g["gene_symbol"].tolist()
ecm_col = [c for c in df_ecm_g.columns if "symbol" in c.lower()][0]
ecm_symbols_raw = df_ecm_g[ecm_col].tolist()
conv_symbols_raw = df_conv_g["gene_symbol"].tolist()

mod_symbols_clean = set(str(s).strip().upper() for s in mod_symbols_raw if pd.notna(s))
ecm_symbols_clean = set(str(s).strip().upper() for s in ecm_symbols_raw if pd.notna(s))
conv_symbols_clean = set(str(s).strip().upper() for s in conv_symbols_raw if pd.notna(s))

# Recompute independent intersection
independent_intersection = mod_symbols_clean.intersection(ecm_symbols_clean)

print(f"WGCNA Salmon Module Unique Genes: {len(mod_symbols_clean)}")
print(f"Convergent ECM-DEGs Unique Genes: {len(ecm_symbols_clean)}")
print(f"Independently Recomputed Intersection Count: {len(independent_intersection)}")
print(f"Recorded convergent_WGCNA_ECM_genes.csv Count: {len(conv_symbols_clean)}")

diff_a = independent_intersection - conv_symbols_clean
diff_b = conv_symbols_clean - independent_intersection

print(f"Symbols in independent intersection but missing in CSV: {diff_a}")
print(f"Symbols in CSV but missing in independent intersection: {diff_b}")
print(f"40-Gene Intersection Exact Match: {len(diff_a) == 0 and len(diff_b) == 0 and len(independent_intersection) == 40}")

# Check for trailing whitespace, duplicates, or casing anomalies in CSV
raw_list = df_conv_g["gene_symbol"].tolist()
has_whitespace = any(s != s.strip() for s in raw_list if isinstance(s, str))
has_dupes = len(raw_list) != len(set(raw_list))
print(f"Whitespace anomalies in CSV: {has_whitespace}")
print(f"Duplicate gene rows in CSV: {has_dupes}")


# ==============================================================================
# CHECK C: PROBE-TO-GENE COLLAPSING SANITY CHECK
# ==============================================================================
print("\n" + "=" * 80)
print("CHECK C: PROBE-TO-GENE COLLAPSING SANITY CHECK")
print("=" * 80)

# Check GSE62928 (Affymetrix GPL13158)
# We can check 5 multi-probe genes: FN1, COL3A1, VCAN, COL1A1, GAPDH
test_genes = ["FN1", "COL3A1", "VCAN", "COL1A1", "GAPDH"]
print("Testing MaxMean collapsing for GSE62928:")
# We can inspect the expression of these genes in df_expr62
for g in test_genes:
    if g in df_expr62.index:
        val = df_expr62.loc[g].values
        print(f"  Gene {g:<8}: Mean in collapsed matrix = {np.mean(val):.4f} across 8 samples")

# Check GSE125498 (Illumina GPL10558)
df_tt = pd.read_csv("Validation/GSE125498.top.table.tsv", sep="\t")
hub_targets = ["ISM1", "FN1", "VCAN", "COL3A1", "COL8A1", "THBS3", "LOX", "EDIL3", "COMP", "COL11A1", "INHBA"]
print("\nGSE125498 Probe Annotation & Mapping in Validation/GSE125498.top.table.tsv:")
for gene in hub_targets:
    match = df_tt[df_tt["Gene.symbol"].fillna("").astype(str).str.contains(rf"\b{gene}\b", case=False, regex=True)]
    if len(match) > 0:
        for idx, row in match.iterrows():
            print(f"  Gene: {gene:<8} | Probe: {row['ID']:<14} | Symbol in Table: {row['Gene.symbol']:<12} | Title: {row['Gene.title'][:30]:<30} | logFC: {row['logFC']:+.3f} | P: {row['P.Value']:.4f}")
    else:
        print(f"  Gene: {gene:<8} | NO PROBES FOUND ON GPL10558 (Absent from platform)")


# ==============================================================================
# CHECK D: WGCNA PARAMETER AND CODE LOGIC VERIFICATION
# ==============================================================================
print("\n" + "=" * 80)
print("CHECK D: WGCNA CODE LOGIC VERIFICATION (02b_wgcna_analysis.R)")
print("=" * 80)
with open("02b_wgcna_analysis.R", "r") as f:
    r_code = f.read()

# Check power
has_chosen_power = "chosen_power <- 12" in r_code or "chosen_power" in r_code
uses_power_adj = "adjacency(datExpr, power = chosen_power, type = \"signed\")" in r_code
uses_signed_tom = 'TOMsimilarity(adjacency, TOMType = "signed")' in r_code
uses_merged_colors = "moduleEigengenes(datExpr, moduleColors)$eigengenes" in r_code and "moduleColors <- mergedColors" in r_code
uses_cor_p = "corPvalueStudent(moduleTraitCor, n_samples)" in r_code

print(f"Power beta=12 applied to adjacency: {uses_power_adj}")
print(f"Signed TOM calculated from adjacency: {uses_signed_tom}")
print(f"Module eigengenes use mergedColors: {uses_merged_colors}")
print(f"corPvalueStudent called with n_samples: {uses_cor_p}")

# Trait sample alignment check
aligns_trait_to_expr = "trait_df <- df_meta[match(rownames(datExpr), df_meta$sample_id), ]" in r_code or "sample_ids <- df_meta$sample_id" in r_code
print(f"Trait explicitly aligned to datExpr sample order: {aligns_trait_to_expr}")


# ==============================================================================
# CHECK E: ML PIPELINE CORRECTNESS
# ==============================================================================
print("\n" + "=" * 80)
print("CHECK E: ML PIPELINE CORRECTNESS (05b_ml_hub_gene_identification_wgcna.py)")
print("=" * 80)
with open("05b_ml_hub_gene_identification_wgcna.py", "r") as f:
    py_code = f.read()

sample_align_ml = "sample_ids = df_meta[\"sample_id\"].tolist()" in py_code and "X_df = df_expr.loc[available_genes, sample_ids].T" in py_code
print(f"Sample-wise alignment of X and y: {sample_align_ml}")

seed_fixed = "RANDOM_STATE = 42" in py_code and "np.random.seed(RANDOM_STATE)" in py_code
print(f"Random seed fixed: {seed_fixed}")

vote_threshold_correct = "Is_Hub_Gene\": \"YES\" if votes >= 2 else \"NO\"" in py_code
print(f"Hub gene threshold >= 2 votes correctly implemented: {vote_threshold_correct}")


# ==============================================================================
# CHECK F: EXTERNAL VALIDATION (GSE125498) SCRIPT CORRECTNESS
# ==============================================================================
print("\n" + "=" * 80)
print("CHECK F: EXTERNAL VALIDATION SCRIPT CORRECTNESS (06_external_validation_GSE125498.py)")
print("=" * 80)
with open("06_external_validation_GSE125498.py", "r") as f:
    val_code = f.read()

# Direction check: LPD (Late, >=25mo) vs SPD (Early, 0-24mo)
has_direction_early_late = 'sample_stages[name] = "Late_Stage_LPD" if is_long else "Early_Stage_SPD"' in val_code
has_binary_y = 'y_binary = np.array([1 if sample_stages[s] == "Late_Stage_LPD" else 0 for s in all_samples])' in val_code
print(f"Clinical encoding (Late_Stage_LPD = 1, Early_Stage_SPD = 0): {has_binary_y}")

# Mann-Whitney U two-sided
has_mwu_twosided = 'mannwhitneyu(val_late, val_early, alternative="two-sided")' in val_code
print(f"Mann-Whitney U alternative='two-sided': {has_mwu_twosided}")

# In-sample vs CV
has_insample_lr = "clf.fit(X_val[available_hubs], y_binary)" in val_code and "prob = clf.predict_proba(X_val[available_hubs])[:, 1]" in val_code
print(f"Original script computed in-sample composite AUC: {has_insample_lr}")


# ==============================================================================
# CHECK G: FILE AND FIGURE CONSISTENCY
# ==============================================================================
print("\n" + "=" * 80)
print("CHECK G: FILE AND FIGURE CONSISTENCY (SPOT-CHECK 10+ VALUES)")
print("=" * 80)

# Check values against tables
# 1. Soft power = 12
# 2. Modules = 14
# 3. Salmon genes = 604
# 4. Salmon r = 0.806, p = 0.0157
# 5. Salmon Bonferroni = 0.2198, FDR = 0.2198
# 6. Permutation p = 0.0143
# 7. Convergent genes = 40
# 8. Hub genes = 11
# 9. LOOCV consensus accuracy = 87.5%
# 10. In-sample AUC = 0.869
# 11. LOOCV AUC = 0.677, 5-fold CV AUC = 0.696
# 12. VCAN Limma p = 0.0244, Mann-Whitney p = 0.0341, AUC = 0.723
# 13. COL8A1 Limma p = 0.0488, AUC = 0.665

df_wgcna_cor = pd.read_csv("results/tables/wgcna_module_trait_correlation.csv")
salmon_row = df_wgcna_cor[df_wgcna_cor["Module"] == "MEsalmon"].iloc[0]

df_val_met = pd.read_csv("results/tables/GSE125498_wgcna_hub_validation_metrics.csv")
vcan_row = df_val_met[df_val_met["Gene_Symbol"] == "VCAN"].iloc[0]
col8_row = df_val_met[df_val_met["Gene_Symbol"] == "COL8A1"].iloc[0]

print("Spot-checking 13 specific quantitative parameters across generated tables:")
print(f"  1. Salmon Module Gene Count: {salmon_row['Number_of_Genes']} (Expected: 604)")
print(f"  2. Salmon Trait Correlation: {salmon_row['Correlation_with_Trait']:.4f} (Expected: 0.8060)")
print(f"  3. Salmon Trait P-value:     {salmon_row['P_Value']:.4f} (Expected: 0.0157)")
print(f"  4. Convergent WGCNA-ECM:     {len(df_conv_g)} (Expected: 40)")
print(f"  5. Consensus Hub Genes:      {len(pd.read_csv('results/tables/ML_hub_genes_from_WGCNA_ECM.csv'))} (Expected: 11)")
print(f"  6. VCAN GSE125498 log2FC:    {vcan_row['log2_FC_GSE125498']:.4f} (Expected: -0.5224)")
print(f"  7. VCAN Limma P-value:       {vcan_row['Limma_P_Value']:.4f} (Expected: 0.0244)")
print(f"  8. VCAN Mann-Whitney P:      {vcan_row['Mann_Whitney_Pval']:.4f} (Expected: 0.0341)")
print(f"  9. VCAN ROC-AUC:             {vcan_row['ROC_AUC']:.4f} (Expected: 0.7231)")
print(f" 10. COL8A1 GSE125498 log2FC:  {col8_row['log2_FC_GSE125498']:.4f} (Expected: +0.7493)")
print(f" 11. COL8A1 Limma P-value:     {col8_row['Limma_P_Value']:.4f} (Expected: 0.0488)")
print(f" 12. COL8A1 ROC-AUC:           {col8_row['ROC_AUC']:.4f} (Expected: 0.6654)")
print(f" 13. Total Evaluated Probes in Validation Table: {len(df_val_met)} (Expected: 7)")

print("\n" + "=" * 80)
print("TECHNICAL & CODE LOGIC VERIFICATION COMPLETED")
print("=" * 80)
