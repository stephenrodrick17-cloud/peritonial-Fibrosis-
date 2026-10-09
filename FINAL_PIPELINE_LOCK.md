# FINAL PIPELINE LOCK
## Peritoneal Fibrosis Transcriptomics Study — Clean-Room Reproducibility Verification
**Date:** 2026-10-07
**Verification:** Full fresh rerun of Scripts 01 -> 02 -> 03 -> convergence scripts -> downstream
**Primary FC Cutoff:** |log2FC| >= 0.585 (1.5-fold raw change) + nominal P < 0.05.
**Auditor note (2026-10-07 correction):** This document is based on freshly rerun outputs from `02_differential_expression.R`, `03_matrisome_filtering.R`, `generate_gse62928_ecm_venn.py`, and `analyze_convergence_wgcna.py`. All numbers below reflect the actual authoritative outputs, not prior stale claims.

---

## §1. Primary Pipeline Verification Table

| Analysis | Fresh Rerun Result (Authoritative) | Source Script / File |
|---|---:|---|
| Gene universe (Limma) | **22,049** | 01_load_qc_preprocess.R |
| Upregulated DEGs (log2FC >= 0.585 [1.5-fold], P < 0.05) | **534** | 02_differential_expression.R |
| Downregulated DEGs (log2FC <= -0.585 [1.5-fold], P < 0.05) | **1,510** | 02_differential_expression.R |
| Total Primary DEGs | **2,044** (534 + 1,510) | 02_differential_expression.R |
| Genome-wide FDR < 0.05 genes (any direction) | **78** | results/tables/GSE62928_all_results.csv |
| Upregulated genes with FDR < 0.05 | **8** | 02_differential_expression.R console output |
| Curated Matrisome genes (Naba et al. reference) | **1,027** | 03_matrisome_filtering.R |
| Matrisome genes in expression universe (22,049) | **975** | 03_matrisome_filtering.R |
| Total ECM-DEGs (both directions, primary cutoff) | **148** (81 up + 67 down) | 03_matrisome_filtering.R |
| Pro-fibrotic ECM-DEGs (upregulated) | **81** | 03_matrisome_filtering.R + CSV count (Upregulated=81) |
| Fold enrichment (upregulated ECM-DEGs) | **3.430x** | 81 / (534 * 975 / 22049) = 81 / 23.614 = 3.430 |
| Hypergeometric P-value (upregulated, one-sided) | **1.321 × 10⁻²²** | scipy.stats.hypergeom (fresh 2026-10-07) |
| GO BP enriched terms | 537 | 04_functional_enrichment.R |
| KEGG enriched pathways | 20 | 04_functional_enrichment.R |

### Historical / Stringent Comparison Cutoff (0.80 log2FC, reference only)
For comparison only — NOT the primary manuscript threshold.
* Up DEGs: 367 / Down DEGs: 1,263 / Total: 1,630
* Up ECM-DEGs: 71 / Fold enrichment ≈ 4.375× / Hypergeom P ≈ 2.593e-26

---

## §2. WGCNA Verification

| Analysis | Fresh Rerun | Source |
|---|---:|---|
| WGCNA input genes | **5,038** | 02b_wgcna_analysis.R (top-5000 most variable + union with the 81 convergent ECM-DEGs) |
| Soft-thresholding power beta | **12** | 02b_wgcna_analysis.R |
| Scale-free R^2 | **0.809** | 02b_wgcna_analysis.R |
| Number of modules | **14** | results/tables/wgcna_module_trait_correlation.csv |
| Salmon module size | **604** | results/tables/wgcna_trait_significant_module_genes.csv |
| Salmon module-trait correlation r | **0.806** | results/tables/wgcna_module_trait_correlation.csv |
| Salmon module P-value (nominal) | **0.0157** | results/tables/wgcna_module_trait_correlation.csv |
| Exact one-sided permutation P (1/70) | **0.0143** | Combinatorial: choose(8,4) = 70 label permutations |
| Family-wise permutation P (6/70, 14 modules) | **0.0857** | 14 modules x exact label-permutation test |
| Bonferroni-adjusted P (14 modules) | **0.220** (exploratory) | 14 x 0.0157 |

---

## §3. Multi-Tier Convergence & Machine Learning

| Analysis | Fresh Rerun Result (Authoritative) | Source |
|---|---:|---|
| 3-way convergent genes (PRIMARY, 0.585 FC) | **44 genes** | results/tables/convergent_WGCNA_ECM_genes.csv (confirmed 44 rows) |
| Intersection logic (PRIMARY) | 534 up-DEGs × 975 measured matrisome × 604 Salmon = **44** | Strict 3-way set cross-check 2026-10-07 |
| Historical 3-way (0.80 FC, reference only) | 40 genes | 367 × 975 × 604 (older pipeline version, not primary) |
| LASSO features | **2** (ISM1, FN1) | 05b_ml_hub_gene_identification_wgcna.py |
| SVM-RFE features | **14** | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |
| Random Forest features | **17** | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |
| XGBoost features | **1** (EDIL3) | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |
| Consensus threshold | >= 2/4 votes | 05b_ml_hub_gene_identification_wgcna.py |
| Consensus hub genes (11) — subset of 44 pool | **11** | results/tables/ML_hub_genes_from_WGCNA_ECM.csv |

All 11 consensus hubs are confirmed present in the 44-gene primary tripartite pool and the 40-gene historical pool.

---

## §4. External Validation (GSE125498, N=33)

| Analysis | Fresh Rerun Result | Source |
|---|---:|---|
| Total samples | **33** | 06_external_validation_GSE125498.py |
| Cases (long-term PD, vintage >= 25 mo) | **13** | results/tables/GSE125498_sample_metadata.csv |
| Controls (short-term PD, vintage <= 24 mo) | **20** | results/tables/GSE125498_sample_metadata.csv |
| Mapped hub genes (on GPL10558) | **7** (FN1, COL3A1, COL8A1, VCAN, THBS3, LOX, ISM1) | results/tables/GSE125498_wgcna_hub_validation_metrics.csv |
| Unmapped hub genes (no GPL10558 probe) | **4** (COL11A1, COMP, EDIL3, INHBA) | Manifest (documented) |
| 7-gene panel apparent (in-sample) AUC | **0.869** | results/tables/roc_auc_detailed_metrics.csv |
| 50x5-fold repeated CV AUC (PRIMARY metric) | **0.678 (per-repeat mean ± 0.058)** | results/tables/roc_auc_detailed_metrics.csv |
| LOOCV AUC | **0.658** | results/tables/roc_auc_detailed_metrics.csv |
| 1,000-iteration label permutation P | **0.0879 (non-significant)** | results/tables/roc_auc_detailed_metrics.csv |

---

## §5. GSEA Verification

| Analysis | Fresh Rerun Result | Source |
|---|---:|---|
| EMT Hallmark NES | **+3.203** | results/tables/gsea_preranked_hallmark_results.csv |
| EMT Nominal P | < 0.0001 | Same file |
| EMT FDR q | < 0.001 | Same file |
| EMT FWER P | 0.000 | Same file |
| EMT Label-Permutation P (N=70 splits) | **0.0143** (1/70, exact minimum) | Pipeline value |

---

## §6. Cutover Audit (0.80 -> 0.585 Primary Threshold)

| Value | Occurrences in production scripts/docs | Classification (2026-10-07) | Action |
|---|---|---|---|
| log2FC = 0.585 in 02_differential_expression.R primary filter | Correct 1 occurrence | A — PRIMARY manuscript cutoff | Retained (confirmed) |
| log2FC = 0.585 in 03_matrisome_filtering.R volcano vline | Correct 1 occurrence | A — PRIMARY manuscript cutoff | Retained (confirmed) |
| log2FC = 0.80 in 02 / 03 comparison blocks | 2 occurrences total | B — historical comparison (clearly labelled) | Retained as reference only |
| "Primary Cutoff log2FC >= 0.80" in audit scripts | Fixed in run_limma_permutation_test.R / permutation_test_ecm.py / run_permutation_suite.py | E — stale label | Fixed: relabelled 0.80 as historical / comparison; 0.585 labelled primary |
| 367 up DEGs claimed at 0.585 in FINAL docs | Corrected in all FINAL_* md files | E — stale/wrong claim | Fixed: 534 at 0.585; 367 only at 0.80 historical |
| 1,263 down DEGs claimed at 0.585 in FINAL docs | Corrected in all FINAL_* md files | E — stale/wrong claim | Fixed: 1,510 at 0.585; 1,263 only at 0.80 historical |
| 71 up ECM-DEGs claimed at 0.585 in FINAL docs | Corrected in all FINAL_* md files | E — stale/wrong claim | Fixed: 81 at 0.585; 71 only at 0.80 historical |
| 132 total ECM-DEGs (old) vs 148 | Corrected in FINAL_PIPELINE_LOCK | E — stale (pre-0.585 recompute) | Fixed: 148 total ECM-DEGs (81 + 67) |
| Tripartite = 40 claimed at 0.585 | Corrected in FINAL_* / AUDIT | E — stale/wrong claim | Fixed: 44 at 0.585; 40 only at 0.80 historical |
| Fold enrichment = 4.38x labelled primary | Corrected in FINAL_* | E — stale 0.80 value | Fixed: 3.430x primary; 4.375x historical |
| Hypergeom P = 2.593e-26 labelled primary | Corrected in FINAL_* | E — stale 0.80 value | Fixed: 1.321e-22 primary; 2.593e-26 historical |

---

## §7. Resolved Discrepancy: 974 vs 975 Matrisome

Root cause: Two different filtering code paths were used in earlier script versions.
* An older version ran `distinct()` before `trimws()`, producing 974 when one gene had a trailing-space near-duplicate.
* The production `03_matrisome_filtering.R` applies `trimws()` before `distinct()`, correctly yielding 975.

Resolution (2026-10-07 rerun): Fresh 03_matrisome_filtering.R confirms 975 is authoritative for "measured matrisome genes in the 22,049-gene expression universe."

---

## §8. Downstream Dependency Check (Post-0.585 Primary Cutover)

| Downstream Analysis | Regenerated with 0.585 primary? | Status Note |
|---|:---:|---|
| DEG tables: GSE62928_DEGs_filtered.csv (534 up, 1,510 down) | YES | Fresh rerun 2026-10-07 |
| Volcano plot (DEG_01_volcano_plot.png) — 0.585 vlines | YES | Regenerated by 02_differential_expression.R |
| MA plot (DEG_02_MA_plot.png) | YES | Regenerated by 02_differential_expression.R |
| ECM intersection tables (975 genes, 148 ECM-DEGs, 81 up) | YES | Fresh rerun 2026-10-07 |
| ECM volcano (ECM_01_volcano_highlight.png) + top ECM bar (ECM_02) | YES | Regenerated by 03_matrisome_filtering.R |
| Venn: 534 DEG ∩ 1027 Matrisome (81 overlap) | YES | `results/figures/venn_gse62928_pro_fibrotic_ecm_81.png` regenerated |
| 3-way Venn: DEG x ECM x Salmon (44 overlap) | YES | `results/figures/venn_wgcna_convergence.png` regenerated (44 genes) |
| Convergence table: convergent_WGCNA_ECM_genes.csv (44 rows) | YES | 44 rows confirmed 2026-10-07 |
| WGCNA (5,038 input genes) — 81-ECM union | YES | Salmon module = 604 genes (unchanged, ECM union at construction) |
| 11 ML consensus hub genes (from >= 2/4 votes) | YES | All 11 hubs are in the new 44-gene tripartite pool; ML outputs stable |
| GSE125498 cross-cohort validation (7 mapped, 4 unmapped) | YES | Independent of primary GSE62928 DEG threshold |
| ROC / AUC metrics (0.869 apparent / 0.678 CV) | YES | Derived from 11 hubs; unchanged |
| Nomogram (5 genes: VCAN, COL8A1, FN1, ISM1, COL3A1) | YES | All 5 in 11 hubs; unchanged |
| GSEA preranked (21,597 ranked universe) | YES | Independent of DEG threshold (uses ranked universe) |
| STRING network (11 hubs: 21 functional / 3 physical edges) | YES | Independent of GSE62928 DEG threshold |
| Immune / Stromal deconvolution | YES | Independent of DEG threshold |

---

## §9. FINAL VERDICT — 2026-10-07 Threshold Cutover

### READY FOR MANUSCRIPT WRITING (0.585 / 1.5-fold PRIMARY THRESHOLD)

All seven readiness criteria are met:

1. Primary threshold consistent: 0.585 (1.5-fold) labelled as PRIMARY across all scripts/docs; 0.80 clearly labelled historical/comparison.
2. Authoritative differential counts reproduced: 534 up, 1,510 down, 2,044 total DEGs (02_differential_expression.R fresh rerun).
3. Authoritative ECM counts reproduced: 81 up, 67 down, 148 total ECM-DEGs; hypergeom P=1.321e-22, 3.430x enrichment (03 + Python cross-check).
4. Tripartite convergence: 44 genes (534 × 975 × 604). All 11 ML consensus hubs confirmed within 44.
5. All downstream analyses verified regenerated (§8). Generated figures reflect correct 0.585 threshold annotations.
6. All FINAL_* markdown source-of-truth documents and README now synchronized with the 0.585 primary numbers.
7. Mandatory caveats updated for 534 / 1,510 / 81 and reported alongside the historical 367 / 1,263 / 71 comparison where appropriate.

---

## §10. Mandatory Manuscript Caveats (Updated 2026-10-07)

These caveats MUST appear in the manuscript and are not negotiable:

* Discovery cohort N = 8 (4 EPS surgical cases vs. 4 controls: 2 first-catheter PD + 2 non-PD uremic). All downstream analyses are explicitly exploratory / hypothesis-generating.
* Primary DEG threshold is |log2FC| >= 0.585 (1.5-fold raw change) + nominal (unadjusted) P < 0.05, chosen for discovery-cohort sensitivity and yielding **534 upregulated and 1,510 downregulated DEGs out of 22,049 tested genes (2,044 total)**. Only 78 genes pass genome-wide Benjamini-Hochberg FDR < 0.05; only 8 of those are upregulated.
* A more stringent historical comparison at log2FC >= 0.80 yields 367 upregulated / 1,263 downregulated / 1,630 total DEGs (retain only as a sensitivity / comparison reference, not as primary result).
* ECM over-representation at primary cutoff: **81 of the 534 upregulated DEGs overlap the 975 measured matrisome genes (15.17%, 3.430-fold enriched, hypergeometric P = 1.321e-22)**. Only 5 of these 81 up ECM-DEGs individually reach genome-wide FDR < 0.05 (BGN, CLEC11A, COL1A1, MXRA5, SERPINE2).
* For GSE125498 (N=33 PD effluent), the PRIMARY reported generalizability metric is the 50-replicate 5-fold stratified cross-validation AUC = 0.678 ± 0.058. The apparent in-sample AUC = 0.869 is explicitly exploratory / optimistic and must not be presented as primary. The 1,000-iteration label-permutation P = 0.0879 (non-significant) MUST be reported alongside any AUC claim.
* WGCNA salmon module (r = 0.806, P = 0.0157, permutation P = 0.0143) is exploratory only; after Bonferroni correction across 14 modules the adjusted P ≈ 0.220.
* ML 11-hub consensus, nomogram (5-gene, apparent C = 0.819, CV AUC = 0.550), cross-cohort classifier, and STRING co-expression claims are all framed as hypothesis-generating derivations for future validation, not validated clinical claims.
