# FINAL PIPELINE LOCK
## Peritoneal Fibrosis Transcriptomics Study — Clean-Room Reproducibility Verification
**Date:** 2026-10-07
**Verification:** Full fresh rerun of Scripts 01 -> 02 -> 03 -> 04 -> downstream
**Auditor note:** All numbers below were reproduced from code/data — zero values were assumed or copied from prior audit documents.

---

## §1. Primary Pipeline Verification Table

| Analysis | Expected / Prior Claim | Fresh Rerun Result | Match? | Source Script / File |
|---|---:|---:|:---:|---|
| Gene universe (Limma) | 22,049 | **22,049** | YES | 01_load_qc_preprocess.R |
| Upregulated DEGs (log2FC >= 0.80, P < 0.05) | 367 | **367** | YES | 02_differential_expression.R |
| Downregulated DEGs (log2FC <= -0.80, P < 0.05) | 1,263 | **1,263** | YES | 02_differential_expression.R |
| Total DEGs | 1,630 | **1,630** | YES | 02_differential_expression.R |
| Genome-wide FDR < 0.05 genes | 78 | **78** | YES | results/tables/GSE62928_all_results.csv |
| Curated Matrisome genes (Naba et al.) | 1,027 | **1,027** | YES | 03_matrisome_filtering.R |
| Matrisome genes in expression universe | **975** | **975** | YES | 03_matrisome_filtering.R |
| ECM-DEGs total (up + down) | 132 | **132** | YES | 03_matrisome_filtering.R |
| Pro-fibrotic ECM-DEGs (upregulated) | 71 | **71** | YES | 03_matrisome_filtering.R |
| Fold enrichment (upregulated DEGs) | 4.38x | **4.375x** | YES | Computed: 71/(367*975/22049) |
| Hypergeometric P-value (upregulated) | < 1e-25 | **2.593e-26** | YES | scipy.stats.hypergeom |
| GO BP enriched terms | 537 | **537** | YES | 04_functional_enrichment.R |
| KEGG enriched pathways | 20 | **20** | YES | 04_functional_enrichment.R |

---

## §2. WGCNA Verification

| Analysis | Expected | Fresh Rerun | Match? | Source |
|---|---:|---:|:---:|---|
| WGCNA input genes | 5,038 | **5,038** | YES | 02b_wgcna_analysis.R (top-5000 union 71 ECM-DEGs) |
| Soft-thresholding power beta | 12 | **12** | YES | 02b_wgcna_analysis.R |
| Scale-free R^2 | 0.809 | **0.809** | YES | 02b_wgcna_analysis.R |
| Number of modules | 14 | **14** | YES | results/tables/wgcna_module_trait_correlation.csv |
| Salmon module size | 604 | **604** | YES | results/tables/wgcna_trait_significant_module_genes.csv |
| Salmon module-trait correlation r | 0.806 | **0.806** | YES | results/tables/wgcna_module_trait_correlation.csv |
| Salmon module P-value | 0.0157 | **0.0157** | YES | results/tables/wgcna_module_trait_correlation.csv |
| Exact one-sided permutation P (1/70) | 0.0143 | **0.0143** | YES | Combinatorial: 8 choose 4 = 70 |
| Family-wise permutation P (6/70) | 0.0857 | **0.0857** | YES | 14 modules x exact test |
| Bonferroni-adjusted P (14 modules) | 0.220 | **0.220** | YES | 14 x 0.0157 |

---

## §3. Multi-Tier Convergence & Machine Learning

| Analysis | Expected | Fresh Rerun | Match? | Source |
|---|---:|---:|:---:|---|
| 3-way convergent genes | 40 | **40** | YES | results/tables/convergent_WGCNA_ECM_genes.csv |
| Intersection logic | 367 up-DEGs x 975 measured matrisome x 604 Salmon | Verified identical | YES | Python cross-check |
| LASSO features | 2 (ISM1, FN1) | **2** | YES | 05b_ml_hub_gene_identification_wgcna.py |
| SVM-RFE features | 14 | **14** | YES | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |
| Random Forest features | 17 | **17** | YES | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |
| XGBoost features | 1 (EDIL3) | **1** | YES | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |
| Consensus threshold | >= 2/4 votes | **>= 2/4 votes** | YES | 05b_ml_hub_gene_identification_wgcna.py |
| Consensus hub genes | 11 | **11** | YES | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |

---

## §4. External Validation (GSE125498, N=33)

| Analysis | Expected | Fresh Rerun | Match? | Source |
|---|---:|---:|:---:|---|
| Total samples | 33 | **33** | YES | 06_external_validation_GSE125498.py |
| Cases (long-term PD) | 13 | **13** | YES | results/tables/GSE125498_sample_metadata.csv |
| Controls (short-term PD) | 20 | **20** | YES | results/tables/GSE125498_sample_metadata.csv |
| Mapped hub genes (on GPL10558) | 7 | **7** | YES | results/tables/GSE125498_wgcna_hub_validation_metrics.csv |
| Unmapped hub genes | 4 | **4** | YES | COL11A1, COMP, EDIL3, INHBA |
| 7-gene panel apparent AUC | 0.869 | **0.869** | YES | results/tables/roc_auc_detailed_metrics.csv |
| 50x5-fold repeated CV AUC (primary) | 0.678 | **0.691 pooled / 0.678 per-repeat mean** | YES | results/tables/roc_auc_detailed_metrics.csv |
| LOOCV AUC | 0.658 | **0.658** | YES | results/tables/roc_auc_detailed_metrics.csv |
| Permutation test P (1000 permutations) | 0.0879 | **0.0879 (non-significant)** | YES | results/tables/roc_auc_detailed_metrics.csv |

---

## §5. GSEA Verification

| Analysis | Expected | Fresh Rerun | Match? | Source |
|---|---:|---:|:---:|---|
| EMT NES | 3.203 | **3.203** | YES | results/tables/gsea_preranked_hallmark_results.csv |
| EMT Nominal P | < 0.0001 | **0.0000** | YES | Same file |
| EMT FDR q | < 0.001 | **0.000** | YES | Same file |
| EMT FWER P | 0.000 | **0.000** | YES | Same file |

---

## §6. Stale Value Audit

| Value | Occurrences Found | Classification | Action |
|---|---|---|---|
| 502 | 0 in production scripts | — | None needed |
| 950 | 0 in production scripts | — | None needed |
| 3.26 | 0 in production scripts | — | None needed |
| 1.08e-18 | 0 in production scripts | — | None needed |
| 974 | 2 in results/MANUSCRIPT_NUMBERS.md | E — stale/incorrect | FIXED -> 975 |
| 0.585 | In 02_differential_expression.R lines 49-65 | B — historical comparison block (labelled) | Retained |
| 0.585 | In 03_matrisome_filtering.R volcano vline | E — stale cosmetic | FIXED -> 0.80 |
| 0.585 | In AUDIT_REPORT.md | B — audit record | Retained |

---

## §7. Resolved Discrepancy: 974 vs 975

Root cause: Two different filtering code paths were used in earlier script versions.
- An older version ran distinct() before trimws(), producing 974 when one gene had a trailing-space near-duplicate.
- The production 03_matrisome_filtering.R applies trimws() before distinct(), correctly yielding 975.

Resolution: Fresh rerun of 03_matrisome_filtering.R confirms 975 is authoritative.
All production files (MANUSCRIPT_NUMBERS.md, FINAL_VERIFIED_NUMBERS.md) now state 975.

---

## §8. Downstream Dependency Check (After 367-DEG Correction)

| Downstream Analysis | Regenerated? | Status |
|---|:---:|---|
| ECM-DEGs (71 up, 61 down, 132 total) | YES | Fresh rerun confirms |
| GO/KEGG enrichment | YES | Re-run via 04_functional_enrichment.R |
| WGCNA (5,038 input genes) | YES | Salmon module uses ECM-DEG union |
| 604 Salmon module genes | YES | Consistent |
| 40 convergent tripartite genes | YES | Python cross-check: perfect 40/40 match |
| 11 ML consensus hub genes | YES | Tables consistent with corrected pipeline |
| GSE125498 cross-cohort validation | YES | 7 mapped genes all from confirmed 11 hubs |
| ROC analysis | YES | Metrics derived from confirmed 11 hubs |
| Nomogram (5 genes) | YES | All 5 genes in confirmed 11 hubs |
| GSEA (preranked on all DEGs) | YES | Uses all 1,630 DEGs (not just 367 up) |
| STRING network | YES | All 11 hubs confirmed |
| Immune infiltration | YES | Independent of DEG threshold |

---

## §9. FINAL VERDICT

### READY FOR MANUSCRIPT WRITING

All six readiness criteria are met:

1. 974/975 discrepancy resolved: The correct authoritative count is 975 (confirmed by fresh 03_matrisome_filtering.R rerun). All production files corrected.
2. 367 DEGs reproduced from actual production pipeline (02_differential_expression.R).
3. 71 ECM-DEGs reproduced (03_matrisome_filtering.R fresh rerun).
4. All downstream analyses affected by the DEG correction have been regenerated (see §8).
5. All major manuscript numbers match the fresh rerun (see §1-§5 — all rows pass).
6. No stale production outputs remain.

---

## §10. Mandatory Manuscript Caveats

These must appear in the manuscript and are not negotiable:

- Discovery cohort N = 8 (4 EPS vs. 4 controls). All analyses are exploratory and hypothesis-generating.
- The 367 DEG threshold (log2FC >= 0.80, P < 0.05 nominal, unadjusted) is intentional for a small-N discovery cohort. Only 78 genes pass genome-wide FDR < 0.05.
- Only 5 of the 71 upregulated ECM-DEGs reach FDR < 0.05 (BGN, CLEC11A, COL1A1, MXRA5, SERPINE2).
- The 50x5-fold CV AUC of 0.678 (NOT the apparent 0.869) is the primary reported discrimination metric for GSE125498.
- The 1,000-permutation panel test P = 0.0879 (non-significant) must be reported alongside the apparent AUC.
- WGCNA, ML, nomogram, and cross-cohort analyses are hypothesis-generating given the small N.
- The Bonferroni-corrected WGCNA module P = 0.220 (exploratory, not significant after correction for 14 modules).
