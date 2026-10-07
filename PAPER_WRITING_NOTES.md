
# PAPER WRITING NOTES — LOCKED MANUSCRIPT NUMBERS & INTERPRETATIONS

## Section 3.1: Discovery Differential Expression (GSE62928)
* **Exact Final Numbers:** Gene universe = 22,049; Upregulated DEGs (pro-fibrotic) = 367; Downregulated DEGs = 1,263 (total DEGs = 1,630). Primary thresholds: log2FC >= 0.80 and nominal P < 0.05. Genome-wide FDR < 0.05 count = 0.
* **Interpretation:** Moderate transcriptomic perturbation in EPS dialysate vs controls, driven by robust pro-fibrotic upregulation.
* **Limitations:** N=8 sample size limit (4 EPS vs 4 Control dialysate samples); strict genome-wide FDR control yields zero genes, necessitating unadjusted P-value filtering with explicit FC thresholds.
* **Figure/Table Reference:** Figure 1 (Volcano plot & DEG heatmap); Table S1.
* **Wording Restrictions:** Use 'pro-fibrotic candidates' or 'upregulated transcriptomic signal'. DO NOT use 'causal drivers' or 'diagnostic biomarkers'.

## Section 3.2: Matrisome ECM Enrichment Analysis
* **Exact Final Numbers:** Curated Matrisome = 1,027; Expression universe Matrisome = 975; Upregulated ECM-DEGs = 71; Downregulated ECM-DEGs = 16. Enrichment ratio = 4.375x (observed 71 vs expected 16.23); Hypergeometric P = 2.5929384e-26; Exact 1-sided permutation P = 0.0143.
* **Interpretation:** Extremely strong selective enrichment of extracellular matrix components within the pro-fibrotic DEG signal.
* **Limitations:** Database-dependent definition of human matrisome.
* **Figure/Table Reference:** Figure 2 (ECM Volcano highlight & enrichment bar chart); Table S2.
* **Wording Restrictions:** Clearly distinguish '71 upregulated ECM-associated genes' from all DEGs.

## Section 3.3: Exploratory WGCNA
* **Exact Final Numbers:** WGCNA input = 5,038 genes; 14 co-expression modules; Salmon module = 604 genes; Trait correlation r = 0.806 (P = 0.0157). Label-permutation P = 0.0143; Family-wise permutation P = 0.0857; Bonferroni P ~ 0.220.
* **Interpretation:** Salmon module represents an exploratory co-expression cluster associated with EPS phenotype.
* **Limitations:** Not statistically robust after multiple-testing correction across 14 modules.
* **Figure/Table Reference:** Figure 3 (WGCNA module-trait heatmap & dendrogram); Table S3.
* **Wording Restrictions:** MUST use 'exploratory co-expression module'. DO NOT call Salmon a 'validated' or 'statistically robust' module.

## Section 3.4: Tripartite Candidate Convergence
* **Exact Final Numbers:** 367 Up-DEGs cap 975 Measured Matrisome cap 604 Salmon Module = 40 convergent genes.
* **Interpretation:** Strict convergence identifies core extracellular matrix genes aligned with co-expression networks.
* **Limitations:** Dependent on WGCNA module boundaries.
* **Figure/Table Reference:** Figure 4 (Venn diagram & 40-gene heatmap); Table S4.
* **Wording Restrictions:** State '975 measured Matrisome genes'.

## Section 3.5: Multi-Model Machine Learning Prioritization
* **Exact Final Numbers:** 11 consensus genes selected by >=2 of 4 ML algorithms (LASSO, SVM-RFE, RF, XGBoost): ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX.
* **Interpretation:** Robust multi-algorithm feature selection prioritizes key ECM remodeling genes.
* **Limitations:** Small discovery N=8 means ML is purely exploratory for prioritization, not predictive modeling.
* **Figure/Table Reference:** Figure 5 (ML selection consensus plot & feature importance); Table 1.
* **Wording Restrictions:** Use 'consensus candidate genes' or 'prioritized candidate genes'. DO NOT use 'validated biomarkers' or 'clinical predictors'.

## Section 3.6: Protein-Protein Interaction & Co-expression Architecture
* **Exact Final Numbers:** 42 of 55 hub gene pairs (76.36%) exhibit co-expression |r| >= 0.85. STRING v12.5 (retrieved 2026-10-01, score >= 0.400) reveals 21 functional edges and 3 physical binding pairs (FN1-LOX, COL11A1-COL3A1, COMP-FN1).
* **Interpretation:** Dense physical and functional connectivity among prioritized ECM candidates.
* **Limitations:** Correlation does not prove direct physical binding for all pairs.
* **Figure/Table Reference:** Figure 6 (PPI network & correlation matrix); Table S5.
* **Wording Restrictions:** Cite exact STRING version and date. Report 42/55 pairs = 76.36%.

## Section 3.7: Pathway & Gene Set Enrichment Analysis (GSEA)
* **Exact Final Numbers:** 21,597 ranked genes. Hallmark EMT NES = +3.203 (FDR < 0.001, permutation P = 0.0143).
* **Interpretation:** Unbiased genome-wide ranking confirms Epithelial-Mesenchymal Transition as the dominant biological pathway.
* **Limitations:** Gene set databases represent general annotations.
* **Figure/Table Reference:** Figure 7 (GSEA EMT plot & pathway summary); Table S6.
* **Wording Restrictions:** Specify ranking metric exactly: sign(log2FC) * (-log10 P).

## Section 3.8: Microenvironmental Cellular Deconvolution
* **Exact Final Numbers:** Myofibroblasts Mann-Whitney U = 16, nominal P = 0.0286, exact permutation P = 0.0286, FDR = 0.1714 (across cell types).
* **Interpretation:** Signal indicates elevated myofibroblast activity in EPS dialysate.
* **Limitations:** Not statistically significant after multiple testing correction across cell types (FDR = 0.1714).
* **Figure/Table Reference:** Figure 8 (Deconvolution boxplots); Table S7.
* **Wording Restrictions:** MUST use 'exploratory myofibroblast-associated signal'. DO NOT claim 'significant expansion'.

## Section 3.9: Cross-Cohort Evaluation (GSE125498)
* **Exact Final Numbers:** N=33 dialysate samples (20 short-term vs 13 long-term PD). 7 of 11 genes mapped (FN1, COL3A1, COL8A1, VCAN, THBS3, LOX, ISM1). Apparent AUC = 0.869; Primary Repeated 5-fold CV AUC = 0.678 (SD 0.058); LOOCV AUC = 0.658; Permutation P = 0.0879.
* **Interpretation:** Moderate cross-cohort biological transferability of the 7-gene signature.
* **Limitations:** Score direction for VCAN and THBS3 was inverted post-hoc based on observed external cohort direction; modest generalization performance (CV AUC = 0.678).
* **Figure/Table Reference:** Figure 9 (External validation ROC curves & expression boxplots); Table 2.
* **Wording Restrictions:** Report 0.678 as primary generalization AUC. MUST use 'cross-cohort evaluation' or 'biological transferability'. DO NOT use 'clinical validation' or 'direct EPS validation'.

## Section 3.10: Exploratory Multivariable Point-Scoring Model
* **Exact Final Numbers:** 5 genes (VCAN, COL8A1, FN1, ISM1, COL3A1). Apparent C-index = 0.819; CV AUC = 0.550; Brier score = 0.1542; Hosmer-Lemeshow P = 0.244.
* **Interpretation:** Demonstrates in-sample feasibility but low cross-validated predictive performance.
* **Limitations:** High optimism gap (CV AUC drops to 0.550).
* **Figure/Table Reference:** Figure 10 (Nomogram & calibration plot); Table S8.
* **Wording Restrictions:** MUST use 'exploratory multivariable point-scoring model'. DO NOT call it a 'clinical nomogram' or 'diagnostic tool'.

## Section 3.11: Single-Cell Cellular Contextualization
* **Status:** SUPPLEMENTARY ONLY (Cellular contextualization, not statistical validation).
* **Wording Restrictions:** DO NOT claim validation of hub genes.

## Section 3.12 & 3.13: Regulatory Network Analyses (miRNA / lncRNA)
* **Status:** SUPPLEMENTARY ONLY (Putative regulatory network based on public databases).
* **Wording Restrictions:** DO NOT claim established ceRNA mechanism. Use 'putative regulatory network'.
