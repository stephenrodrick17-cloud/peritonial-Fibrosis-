# FINAL MANUSCRIPT SOURCE OF TRUTH

## Peritoneal Fibrosis Transcriptomics Study

**Date:** 2026-10-07  
**Status:** locked & verified via Clean-Room Pipeline Rerun (Scripts 01 -> 09)

---

## 1. Primary Discovery Differential Expression (GSE62928)

* **Expression Universe:** 22,049 genes
  * *Source Script:* 01_data_preprocessing.R / 02_differential_expression.R
  * *Output File:* results/tables/GSE62928_all_results.csv
  * *Method:* limma / empirical Bayes moderated t-test
  * *Status:* Freshly reproduced (Discovery)
* **Upregulated DEGs (Pro-fibrotic):** 367 genes
  * *Threshold:* log2FC >= 0.80 AND nominal P < 0.05
  * *Source Script:* 02_differential_expression.R
  * *Output File:* results/tables/GSE62928_DEGs_filtered.csv
  * *Status:* Freshly reproduced (Discovery)
* **Downregulated DEGs:** 1,263 genes
  * *Threshold:* log2FC <= -0.80 AND nominal P < 0.05
  * *Source Script:* 02_differential_expression.R
  * *Output File:* results/tables/GSE62928_DEGs_filtered.csv
  * *Status:* Freshly reproduced (Discovery). (Note: 1,112 was an obsolete historical artifact from an early draft).
* **Total Primary DEGs:** 1,630 genes (367 + 1,263)
* **Genome-wide FDR < 0.05 Count:** 0 genes (under strict genome-wide FDR control)

---

## 2. Matrisome ECM Enrichment Analysis

* **Curated Human Matrisome Total:** 1,027 genes
  * *Source Script:* 03_matrisome_filtering.R
* **Matrisome Genes Represented in Expression Universe:** 975 genes
  * *Source Script:* 03_matrisome_filtering.R
  * *Output File:* results/tables/GSE62928_ECM_intersection.csv
  * *Status:* Freshly reproduced
* **Upregulated ECM-DEGs (Pro-fibrotic ECM):** 71 genes
  * *Definition:* 367 Upregulated DEGs intersect 975 Measured Matrisome genes
  * *Source Script:* 03_matrisome_filtering.R
  * *Output File:* results/tables/GSE62928_ECM_DEGs_candidate_list.csv
* **Downregulated ECM-DEGs:** 16 genes (1,263 Down-DEGs intersect 975 Matrisome)
* **Enrichment Ratio:** 4.375x (observed 71 vs expected 16.23)
  * *Source Script:* 03_matrisome_filtering.R
* **Hypergeometric P-value:** 2.5929384e-26
  * *Source Script:* 03_matrisome_filtering.R
* **Exact Permutation P-value (1-sided):** P = 1/70 = 0.0142857 (P = 0.0143)

---

## 3. Weighted Gene Co-expression Network Analysis (WGCNA)

* **Input Genes:** 5,038 genes (Top 5,000 MAD + 38 Matrisome additions)
  * *Source Script:* 02b_wgcna_analysis.R
* **Module Count:** 14 co-expression modules
* **Parameters:** Soft-thresholding power beta = 12, signed network, Pearson correlation, minModuleSize = 30, mergeCutHeight = 0.25
* **Salmon Module Size:** 604 genes
* **Salmon-EPS Trait Correlation:** r = 0.806 (r = 0.805997)
* **Salmon Nominal P-value:** P = 0.0157
* **Salmon Exact Label-Permutation P-value:** P = 1/70 = 0.0143
* **Family-wise Permutation P-value:** P = 0.0857 (6/70)
* **Bonferroni-adjusted P-value:** P ~ 0.220 (14 x 0.0157)
* **Module Status:** Exploratory co-expression module (not statistically robust after family-wise/Bonferroni correction).

---

## 4. Tripartite Convergence & Machine Learning Prioritization

* **Tripartite Convergence (40 Genes):**
  * *Definition:* 367 Up-DEGs intersect 975 Measured Matrisome intersect 604 Salmon Module = 40 genes
  * *Source Script:* 04_hub_gene_selection.R
* **Multi-Model ML Selection (11 Consensus Candidates):**
  * *Models:* LASSO, SVM-RFE, Random Forest, XGBoost
  * *Threshold:* Selection by >= 2 of 4 ML feature selection algorithms
  * *Consensus Genes:* ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX
  * *Source Script:* 05_machine_learning_selection.py
  * *Status:* Exploratory prioritized candidates

---

## 5. External Cohort Cross-Evaluation (GSE125498)

* **Sample Count:** 33 peritoneal dialysate samples (20 short-term PD, 13 long-term PD)
* **Gene Mapping:** 7 of 11 consensus genes mapped to platform (FN1, COL3A1, COL8A1, VCAN, THBS3, LOX, ISM1)
* **Performance Metrics:**
  * *Apparent / In-Sample Composite 7-Gene AUC:* 0.869 (95% CI: 0.742-0.996) -- Exploratory/Optimistic
  * *Repeated 5-fold CV AUC (Primary Generalization Estimate):* 0.678 (SD 0.058)
  * *LOOCV AUC:* 0.658
  * *Single 5-fold CV AUC:* 0.592
  * *Label Permutation P-value:* P = 0.0879
* **Score Inversion Note:** VCAN and THBS3 score directions were inverted based on observed external cohort direction (Post-hoc / Data-informed limitation).

---

## 6. Exploratory Multivariable Model (Nomogram)

* **Input Genes:** 5 genes (VCAN, COL8A1, FN1, ISM1, COL3A1)
* **Apparent C-index:** 0.819
* **CV AUC:** 0.550
* **Brier Score:** 0.1542
* **Hosmer-Lemeshow P-value:** P = 0.244
* **Status:** Exploratory multivariable point-scoring model (not a validated clinical nomogram).

---

## 7. Pathway & Pathway Enrichment (GSEA)

* **Ranked Universe:** 21,597 genes (ranked by sign(log2FC) * (-log10 P))
* **Collection:** MSigDB Hallmark (50 gene sets)
* **EMT NES:** +3.203 (P < 0.001, FDR < 0.001, Permutation P = 0.0143)

---

## 8. Transcriptomic Correlation vs STRING PPI

* **Hub Pair Co-expression:** 42 of 55 unique pairs (11 x 10 / 2) show |r| >= 0.85 (42/55 = 76.36%).
* **STRING PPI v12.5 (2026-10-01):**
  * Score threshold >= 0.400
  * 21 functional edges, 3 physical edges (FN1-LOX, COL11A1-COL3A1, COMP-FN1).
