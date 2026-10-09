# FINAL MANUSCRIPT SOURCE OF TRUTH

## Peritoneal Fibrosis Transcriptomics Study

**Date:** 2026-10-07
**Status:** Locked & verified via clean-room pipeline rerun (Scripts 01 -> 02 -> 03 -> convergence scripts)
**Primary log2FC Threshold:** |log2FC| >= 0.585 (1.5-fold raw change), nominal P < 0.05. Historical 0.80 FC threshold retained for comparison only.

---

## 1. Primary Discovery Differential Expression (GSE62928)

* **Expression Universe:** 22,049 genes
  * *Source Script:* 01_data_preprocessing.R / 02_differential_expression.R
  * *Output File:* results/tables/GSE62928_all_results.csv
  * *Method:* limma / empirical Bayes moderated t-test
  * *Status:* Freshly reproduced (Discovery)
* **Upregulated DEGs (Pro-fibrotic):** **534 genes**
  * *Threshold:* log2FC >= 0.585 (1.5-fold raw change) AND nominal P < 0.05
  * *Source Script:* 02_differential_expression.R
  * *Output File:* results/tables/GSE62928_DEGs_filtered.csv
  * *Status:* Freshly reproduced (Discovery)
* **Downregulated DEGs:** **1,510 genes**
  * *Threshold:* log2FC <= -0.585 (1.5-fold raw change) AND nominal P < 0.05
  * *Source Script:* 02_differential_expression.R
  * *Output File:* results/tables/GSE62928_DEGs_filtered.csv
  * *Status:* Freshly reproduced (Discovery)
* **Total Primary DEGs:** **2,044 genes** (534 + 1,510)
* **Genome-wide FDR < 0.05 Count (any direction):** **78 genes**
* **Upregulated with FDR < 0.05:** 8 genes

### Historical / Comparison Threshold (reference only)
* **log2FC >= 0.80, P < 0.05:** Up DEGs = 367; Down DEGs = 1,263; Total = 1,630. Used in earlier pipeline versions as a stricter filter.

---

## 2. Matrisome ECM Enrichment Analysis

* **Curated Human Matrisome Total (Naba et al.):** 1,027 genes
  * *Source Script:* 03_matrisome_filtering.R
* **Matrisome Genes Represented in 22,049 Expression Universe:** 975 genes
  * *Source Script:* 03_matrisome_filtering.R
  * *Output File:* results/tables/GSE62928_ECM_intersection.csv
  * *Status:* Freshly reproduced
* **Upregulated (Pro-Fibrotic) ECM-DEGs:** **81 genes**
  * *Definition:* 534 Upregulated DEGs intersect 975 Measured Matrisome genes
  * *Source Script:* 03_matrisome_filtering.R + Python hypergeometric cross-check
  * *Output File:* results/tables/GSE62928_ECM_DEGs_candidate_list.csv (DEG_Status='Upregulated')
* **Downregulated ECM-DEGs:** 67 genes
* **Total ECM-DEGs (both directions):** 148 genes
* **Up ECM-DEG Enrichment Ratio:** 3.430x (observed 81 vs expected 23.61 for 534 × 975 / 22,049)
  * *Source Script:* 03_matrisome_filtering.R + scipy cross-check
* **Up ECM-DEG Hypergeometric P-value (one-sided):** 1.321e-22
  * *Source Script:* scipy.stats.hypergeom (fresh cross-check 2026-10-07)
* **Exact Permutation P-value (1-sided, N=70 splits, Primary cutoff):** P = 1/70 = 0.0142857 (P = 0.0143)

### Historical Threshold (0.80 FC, reference only)
* Upregulated ECM-DEGs: 71 genes, 4.375x enrichment, P(hypergeom) = 2.593e-26.

---

## 3. Weighted Gene Co-expression Network Analysis (WGCNA)

* **Input Genes:** 5,038 genes (Top 5,000 MAD + 38 Matrisome additions; union with the 81 convergent ECM-DEGs)
  * *Source Script:* 02b_wgcna_analysis.R
* **Module Count:** 14 co-expression modules
* **Parameters:** Soft-thresholding power beta = 12, signed hybrid network, Pearson correlation, minModuleSize = 30, mergeCutHeight = 0.25
* **Salmon Module Size:** 604 genes
* **Salmon-EPS Trait Correlation:** r = 0.806 (r = 0.805997)
* **Salmon Nominal P-value:** P = 0.0157
* **Salmon Exact Label-Permutation P-value (1-sided):** P = 1/70 = 0.0143
* **Family-wise Permutation P-value (14 modules):** P = 0.0857 (6/70)
* **Bonferroni-adjusted P-value (14 modules):** P ~ 0.220 (14 x 0.0157)
* **Module Status:** Exploratory co-expression module (not statistically robust after family-wise/Bonferroni correction, discovery cohort framing).

---

## 4. Tripartite Convergence & Machine Learning Prioritization

* **Tripartite Convergence (Primary 0.585 FC threshold):** **44 genes**
  * *Definition:* 534 Up-DEGs intersect 975 Measured Matrisome intersect 604 Salmon Module
  * *44 Gene List:* ADAMTS1, ADAMTS2, ADAMTS6, BGN, CLCF1, CLEC11A, COL11A1, COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, COL8A1, COMP, EDIL3, FBLN7, FGF14, FN1, FRAS1, INHBA, ISM1, LOX, MFAP2, MFGE8, MMP14, MUC1, MXRA5, PCOLCE, PDGFA, POSTN, SDC1, SERPINA1, SERPINA3, SERPINE1, SERPINE2, SNED1, SULF2, TGFB3, THBS1, THBS2, THBS3, THSD4, TNFAIP6, VCAN
  * *Source Script:* analyze_convergence_wgcna.py + strict 3-way cross-check
  * *Output File:* results/tables/convergent_WGCNA_ECM_genes.csv (44 rows)
* **Historical Tripartite (0.80 FC, reference only):** 40 genes (367 Up-DEGs × 975 Matrisome × 604 Salmon)
* **Multi-Model ML Selection (11 Consensus Candidates):**
  * *Candidate Pool:* 44 tripartite convergent genes (primary 0.585); historically computed over 40 (0.80 FC). All 11 consensus hubs are present in both pools.
  * *Models:* LASSO, SVM-RFE, Random Forest, XGBoost
  * *Threshold:* Selection by >= 2 of 4 ML feature selection algorithms
  * *Consensus Genes (11):* ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX
  * *Source Script:* 05_machine_learning_selection.py / 05b_ml_hub_gene_identification_wgcna.py
  * *Status:* Exploratory prioritized candidates (all 11 confirmed in 44-gene pool)

---

## 5. External Cohort Cross-Evaluation (GSE125498)

* **Sample Count:** 33 peritoneal dialysate samples (20 short-term PD, 13 long-term PD)
* **Gene Mapping:** 7 of 11 consensus genes mapped to platform (FN1, COL3A1, COL8A1, VCAN, THBS3, LOX, ISM1)
* **Performance Metrics:**
  * *Apparent / In-Sample Composite 7-Gene AUC:* 0.869 (95% CI: 0.742-0.996) -- Exploratory/Optimistic
  * *Repeated 5-fold CV AUC (Primary Generalization Estimate):* 0.678 (SD 0.058)
  * *LOOCV AUC:* 0.658
  * *Single 5-fold CV AUC:* 0.592
  * *1,000-Iteration Label Permutation P-value:* P = 0.0879 (non-significant)
* **Score Inversion Note:** VCAN and THBS3 score directions were inverted based on observed external cohort direction (Post-hoc / Data-informed limitation).

---

## 6. Exploratory Multivariable Model (Nomogram)

* **Input Genes:** 5 genes (VCAN, COL8A1, FN1, ISM1, COL3A1)
* **Apparent C-index:** 0.819 (95% CI: 0.623-0.968)
* **CV AUC (Generalization):** 0.550 +/- 0.178
* **Brier Score:** 0.1542
* **Hosmer-Lemeshow P-value:** P = 0.244
* **Likelihood-Ratio Test vs. FN1 alone:** LR chi2 = 10.421, df = 4, P = 0.0339
* **Status:** Exploratory multivariable point-scoring model (not a validated clinical nomogram; framing as 'derivation for future validation' only).

---

## 7. Pathway & Gene Set Enrichment Analysis (GSEA)

* **Ranked Universe:** 21,597 genes (ranked by sign(log2FC) * (-log10 P))
* **Collection:** MSigDB Hallmark (50 gene sets)
* **Top Upregulated Pathway — EMT:** NES = +3.203 (Nominal P < 0.0001, FDR q < 0.001, Permutation P = 1/70 = 0.0143)
  * Label: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

---

## 8. Transcriptomic Correlation vs STRING PPI

* **11 Hub Pair Co-expression:** 42 of 55 unique hub-gene pairs show |Pearson r| >= 0.85 (42/55 = 76.36%).
* **STRING PPI v12.5 (retrieved 2026-10-01, score >= 0.400):**
  * 21 functional (combined_score) edges across the 11 hubs.
  * 3 physical (binding/action) edges: FN1-LOX (0.848), COL11A1-COL3A1 (0.720), COMP-FN1 (0.595).
  * ISM1: degree = 0 (isolated node, no STRING edges at threshold 0.400).

---

## 9. FDR-Significant Individual ECM-DEGs

* **5 of 81 Up ECM-DEGs pass genome-wide FDR < 0.05:** BGN, CLEC11A, COL1A1, MXRA5, SERPINE2
* **Plus 2 Downregulated ECM genes:** TIMP4, LEP
* The remaining 76 up ECM-DEGs are only nominally significant (nominal P < 0.05), consistent with small-N discovery framing.
