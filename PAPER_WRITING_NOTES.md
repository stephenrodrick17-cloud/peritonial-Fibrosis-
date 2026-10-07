# Peritoneal Membrane Fibrogenesis: Verified Paper Writing Notes

> **Notice to Manuscript Authors:** This document contains strictly verified facts, methodologies, and quantitative results derived from the code and data in this repository. All numerical values reflect `results/manuscript_numbers.json` and `results/MANUSCRIPT_NUMBERS.md`. Cautious, evidence-based scientific language is used throughout.

---

## 1. Title Options

1. **Integrated Matrisome Transcriptomics and Consensus Machine Learning Identify Candidate Extracellular Matrix Biomarkers Associated with Peritoneal Dialysis-Induced Membrane Remodeling** *(Recommended)*
2. **Prioritizing Candidate Matrisome Drivers of Peritoneal Membrane Fibrogenesis Through Co-Expression Networks and Supervised Machine Learning**
3. **Cross-Cohort Evaluation of an Extracellular Matrix Gene Signature in Peritoneal Dialysis-Associated Membrane Fibrosis and Dialysis Vintage**

---

## 2. Structured Abstract

### Background
Encapsulating Peritoneal Sclerosis (EPS) and peritoneal membrane fibrosis represent severe, life-threatening complications of long-term peritoneal dialysis (PD), characterized by extensive mesothelial denudation, submesothelial fibrocollagenous deposition, and neoangiogenesis. Identifying candidate extracellular matrix (ECM) regulators and non-invasive biomarker candidates is critical for understanding disease pathogenesis and tracking membrane injury.

### Methods
We established a multi-tiered bioinformatics and consensus machine learning pipeline. Discovery transcriptomic profiling was performed on human parietal peritoneal biopsy samples (**GSE62928**, $N = 8$: 4 surgical EPS cases vs. 4 non-EPS uremic/PD controls; 22,049 probe-level genes, 20,940 MaxMean-collapsed genes). Pro-fibrotic differentially expressed genes ($\text{DEGs}$, nominal $P < 0.05, \log_2\text{FC} \ge 0.80$) were intersected with the curated Human Matrisome Database (1,027 genes). Exploratory Weighted Gene Co-Expression Network Analysis (WGCNA) was performed on the top 5,000 variable genes unioned with the 71 ECM-DEGs (5,038 features) to identify pro-fibrotic modules. Tripartite convergent candidates were prioritized using four supervised machine learning feature-selection algorithms (LASSO, SVM-RFE, Random Forest, XGBoost). Candidate biomarkers were evaluated for cross-cohort biological transferability in human peritoneal dialysis effluent cell pellets (**GSE125498**, $N = 33$: 20 short-term [$\le 24$ mo] vs. 13 long-term [$\ge 25$ mo] PD vintage). Combinatorial label permutations ($\binom{8}{4} = 70$ splits) and STRING v12.5 functional protein association networks provided empirical statistical baselines and structural connectivity.

### Results
Discovery analysis isolated 367 pro-fibrotic up-regulated DEGs, demonstrating marked over-representation of extracellular matrix genes (**71 ECM-DEGs**, $19.35\%$ vs. $4.42\%$ background; fold enrichment $= 4.38$; hypergeometric $P = 2.59 \times 10^{-26}$; exact one-sided label permutation $P_{\text{perm}} = 1/70 = 0.0143$). Exploratory WGCNA identified the **Salmon module** (604 genes, Pearson $r = +0.806, P = 0.0157$; one-sided $P_{\text{perm}} = 1/70 = 0.0143$, two-sided $P_{\text{perm}} = 2/70 = 0.0286$, Bonferroni $P = 0.220$). Intersecting Salmon module genes with ECM-DEGs yielded **40 tripartite convergent candidates**. Consensus machine learning ($\ge 2/4$ algorithm votes) prioritized **11 consensus candidate hub genes**: `ISM1`, `FN1`, `EDIL3`, `VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, and `LOX`. In external effluent cells ($N = 33$, 7 mapped hubs), a composite 7-gene logistic regression model showed an apparent in-sample $\text{AUC} = 0.869$ ($95\%\text{ CI: } [0.710, 0.992]$), while 50-repeat 5-fold cross-validation yielded an out-of-fold generalization $\text{AUC}_{\text{CV}} = 0.678$ ($\text{SD } 0.058$; Leave-One-Out $\text{AUC} = 0.658$; 1,000-permutation empirical null $P = 0.0879$, reflecting non-significance under stringent cross-validation). Single-cell deconvolution and GSEA confirmed primary enrichment in Epithelial-Mesenchymal Transition ($\text{NES} = +3.203, \text{FDR} < 0.001$) and stromal myofibroblast expansion ($P = 0.0286$).

### Conclusion
Multi-tiered transcriptomic and machine learning analysis prioritizes 11 candidate matrisome genes associated with peritoneal fibrogenesis. While tissue-derived signatures exhibit strong discovery associations, cross-cohort effluent discrimination is moderate ($\text{AUC}_{\text{CV}} \approx 0.68$), reflecting differences between whole-tissue fibrotic biopsies and shed effluent cellular transcriptomes. These candidate genes provide exploratory targets for future mechanistic and prospective biofluid studies.

---

## 3. Introduction

### Clinical Problem
Peritoneal dialysis (PD) is an established renal replacement modality for end-stage kidney disease. However, continuous exposure of the peritoneal membrane to bioincompatible dialysates, glucose degradation products, recurrent peritonitis, and uremic inflammation leads to chronic membrane damage. This manifests as progressive submesothelial fibrosis, vascular hyalinization, neoangiogenesis, and in severe instances, Encapsulating Peritoneal Sclerosis (EPS). EPS is characterized by dense fibrocollagenous cocooning of the intestines, bowel obstruction, and high clinical mortality.

### Knowledge Gap & Rationale
Although transforming growth factor-$\beta$ (TGF-$\beta$) signaling and mesothelial-to-mesenchymal transition (MMT) are recognized drivers of peritoneal injury, the dynamic remodeling of the peritoneal extracellular matrix (the "matrisome") remains incompletely characterized. Understanding which specific structural collagens, matricellular glycoproteins, and remodeling enzymes are co-regulated in fibrotic tissue can reveal candidate pathways involved in membrane failure.

### Research Objective
This study was designed to:
1. Objectively map the pro-fibrotic matrisome in human peritoneal biopsy transcriptomics using differential expression, exact combinatorial permutations, and exploratory co-expression network analysis (WGCNA).
2. Prioritize robust candidate biomarkers using four complementary machine learning algorithms.
3. Assess the cross-cohort biological transferability and discrimination of the candidate signature in human peritoneal effluent cells across dialysis vintage.

---

## 4. Methods

### 4.1 Datasets and Cohort Descriptions
- **Discovery Tissue Dataset (GSE62928):** Parietal peritoneal biopsy samples from 8 human subjects profiled on the Affymetrix Human Genome U133 Plus 2.0 Array (GPL13158/GPL570). The cohort comprises 4 patients with surgical Encapsulating Peritoneal Sclerosis (EPS) undergoing enterolysis and 4 non-EPS controls (2 non-dialysis uremic controls undergoing abdominal surgery and 2 incident PD patients at the time of initial catheter insertion).
- **External Dialysis Effluent Dataset (GSE125498):** Peritoneal effluent cell pellet transcriptomes from 33 human PD patients profiled on the Illumina HumanHT-12 v4.0 Expression BeadChip (GPL10558). Samples are grouped by dialysis vintage: 20 short-term PD patients ($\le 24$ months, SPD) vs. 13 long-term PD patients ($\ge 25$ months, LPD), serving as an independent proxy for cumulative peritoneal membrane injury.

### 4.2 Preprocessing and Gene Universes
- **Probe Collapse & Normalization:** Robust Multichip Average (RMA) normalized log2 intensity data were processed. In probe-to-gene mapping, probe collapse for differential expression was performed using the lowest nominal P-value approach (22,049 unique gene symbols). For WGCNA and network modeling, probe collapse used the `MaxMean` average intensity method (20,940 unique gene symbols).
- **Human Matrisome Reference:** The curated Human Matrisome Database (Naba et al., 1,027 genes) was utilized as the benchmark reference for core matrisome (collagens, glycoproteins, proteoglycans) and matrisome-associated (remodeling enzymes, ECM regulators, secreted factors) categories.

### 4.3 Differential Expression and Over-Representation Testing
- Linear modeling and empirical Bayes moderation was performed using `limma`.
- Pro-fibrotic DEGs were filtered using predefined criteria: nominal unadjusted $P < 0.05$ and $\log_2\text{FC} \ge 0.80$ (upregulated in EPS vs. controls), identifying 367 genes.
- Over-representation of matrisome genes among DEGs was evaluated via Hypergeometric Fisher's exact test ($M = 22,049$, $k = 975$, $N = 367$, $x = 71$) and verified through exact combinatorial label permutations.
- *FDR Note:* Under Benjamini-Hochberg FDR $< 0.05$, 5 of the 71 upregulated ECM-DEGs (`COL1A1`, `BGN`, `CLEC11A`, `SERPINE2`, `MXRA5`) and 2 downregulated ECM genes (`TIMP4`, `LEP`) reach genome-wide significance, confirming the hypothesis-generating nature of the full 71-gene set.

### 4.4 Exploratory WGCNA Network Construction
- Because discovery sample size ($N = 8$) is below standard WGCNA recommendations ($N \ge 15-20$), network analysis was conducted as an exploratory co-expression screening tool.
- Input features: Top 5,000 most variable genes unioned with the 71 convergent ECM-DEGs (5,038 unique features from the 20,940 MaxMean matrix).
- Signed hybrid adjacency was computed using soft-thresholding power $\beta = 12$ (truncated scale-free $R^2 = 0.809$, linear $R^2 = 0.820$).
- Topological Overlap Matrix (TOM) dynamic tree cutting identified 14 modules. Pearson correlation between module eigengenes and binary clinical EPS status prioritized the pro-fibrotic Salmon module (604 genes, $r = +0.806, P = 0.0157$).

### 4.5 Consensus Machine Learning Feature Selection
The 40 tripartite convergent candidates (Salmon module $\cap$ 71 ECM-DEGs) were screened across 4 algorithms (fixed `random_state = 42`):
1. **LASSO Logistic Regression:** `LogisticRegressionCV(Cs=20, cv=4, penalty='l1', solver='liblinear')` with inverse regularization parameter $C = 4.28$ (2 features selected: `ISM1`, `FN1`).
2. **SVM-RFE:** Linear kernel Support Vector Machine with recursive feature elimination, selecting top 35% ranked features (14 features selected).
3. **Random Forest:** 500 estimators with Mean Decrease Gini thresholded at the empirical mean (17 features selected).
4. **XGBoost:** Gradient-boosted decision trees with depth constraints (1 feature selected: `EDIL3`).
- **Consensus Rule:** Features selected by $\ge 2$ algorithms were retained as the 11 consensus candidate hub genes.

### 4.6 Cross-Cohort Effluent Analysis & Validation Metrics
- Hub genes were mapped to Illumina probe IDs in GSE125498. Seven genes were profiled (`FN1`, `COL3A1`, `COL8A1`, `VCAN`, `THBS3`, `LOX`, `ISM1`), while four lacked mapped probes (`COL11A1`, `COMP`, `EDIL3`, `INHBA`).
- Classifiers were trained using L2-penalized logistic regression ($C = 1.0$) inside a cross-validation pipeline with `StandardScaler` to prevent data leakage.
- Metrics evaluated:
  - Apparent in-sample AUC
  - 50-repeat 5-fold stratified cross-validation AUC (primary generalization metric)
  - Leave-One-Out cross-validation (LOOCV) AUC
  - 1,000-iteration label permutation test on the CV pipeline

### 4.7 Exact Combinatorial Permutation Framework
For discovery cohort GSE62928, the balanced assignment of 4 EPS cases vs. 4 controls permits exactly $\binom{8}{4} = 70$ distinct sample label configurations. Exact one-sided ($P = k/70$) and two-sided ($P = 2k/70$) empirical significance values were computed across all 70 splits to establish exact non-parametric empirical bounds.

---

## 5. Programmatically Verified Results

### 5.1 Discovery Transcriptomics & Matrisome Convergence
- **Gene Universes:** 22,049 probe-level genes; 20,940 MaxMean genes; 21,597 GSEA genes.
- **Pro-Fibrotic DEGs:** 367 genes ($P < 0.05, \log_2\text{FC} \ge 0.80$).
- **Matrisome Over-representation:** 71 ECM-DEGs ($19.35\%$ vs. $4.42\%$ background; fold enrichment $= 4.38$; hypergeometric $P = 2.59 \times 10^{-26}$).
- **Exact Permutation Rank:** True clinical split achieved rank 1/70 ($P_{\text{perm}} = 1/70 = 0.0143$, minimum attainable one-sided $P$ at $N=8$).

### 5.2 Co-Expression Modularity & Tripartite Intersection
- **Salmon Module:** 604 genes; trait correlation $r = +0.806, P = 0.0157$; exact one-sided $P_{\text{perm}} = 1/70 = 0.0143$; two-sided $P_{\text{perm}} = 2/70 = 0.0286$; Bonferroni-adjusted FWER across 14 modules $P_{\text{bonf}} = 0.220$.
- **Tripartite Convergence:** 40 genes located at the intersection of 367 DEGs, 1,027 Matrisome, and 604 Salmon module genes.

### 5.3 Consensus Machine Learning Hub Panel
- **Candidate Hubs ($\ge 2/4$ votes):** 11 genes
  - 3 Votes (3 genes): `ISM1`, `FN1` (LASSO + SVM-RFE + RF), `EDIL3` (SVM-RFE + RF + XGBoost)
  - 2 Votes (8 genes, selected by SVM-RFE + RF): `VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`
  - 1 Vote: 9 genes; 0 Votes: 20 genes
- **STRING PPI Connectivity (STRING v12.5, score $\ge 0.400$):** 21 functional association edges across 11 hubs; 3 physical binding edges (`FN1 - LOX` 0.848, `COL11A1 - COL3A1` 0.720, `COMP - FN1` 0.595); `ISM1` has degree 0.
- **Discovery Co-expression ($|r| \ge 0.85$):** 42 / 55 (76.36%) hub pairs exceed threshold (vs. 1.68% random baseline, 5.47% DEG baseline, 9.27% Salmon module baseline).

### 5.4 External Dialysis Effluent Validation (GSE125498, N=33)
- **Cohort:** 20 Short-term PD vs. 13 Long-term PD effluent pellet transcriptomes.
- **Discrimination Metrics (7-gene panel):**
  - In-sample Logistic Regression: $\text{AUC} = 0.869$ ($95\%\text{ CI: } [0.710, 0.992]$)
  - 50-repeat 5-fold Stratified CV (Pipeline Scaled, Primary): $\text{AUC}_{\text{CV}} = 0.678 \pm 0.058$ (per-fold mean $\text{AUC} = 0.701 \pm 0.203$)
  - Leave-One-Out CV: $\text{AUC}_{\text{LOOCV}} = 0.658$
  - Single-Split Seed-42 5-Fold CV: $\text{AUC} = 0.592$
  - 1,000-Permutation Empirical Test: Null mean $\text{AUC} = 0.4883$, $P = 88/1001 = 0.0879$ (non-significant)
- **5-Gene Nomogram:** In-sample C-index $= 0.819$, Brier score $= 0.1542$, Hosmer-Lemeshow $P = 0.2437$; 5-fold cross-validated $\text{AUC} = 0.550 \pm 0.178$.
- **Individual Gene Trends:**
  - `VCAN`: $\text{AUC} = 0.723$ ($P = 0.034$, decreased in late effluent cells, inverted ROC)
  - `COL8A1`: $\text{AUC} = 0.665$ (nominal lowest-P probe `ILMN_2402392` $P = 0.0488$, Bonferroni $P = 0.0976$; second probe `ILMN_1685433` $P = 0.392$)
  - `FN1`: $\text{AUC} = 0.612$ ($P = 0.294$)
  - `THBS3`: $\text{AUC} = 0.596$ ($P = 0.367$, inverted ROC)
  - `COL3A1`: $\text{AUC} = 0.535$, `ISM1`: $\text{AUC} = 0.527$, `LOX`: $\text{AUC} = 0.496$.

### 5.5 Pathway & Microenvironment Signatures
- **Preranked GSEA Top Pathways:**
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`: $\text{NES} = +3.2033, \text{NOM } P < 0.0001, \text{FDR } q < 0.001$ ($P_{\text{perm}} = 1/70 = 0.0143$)
  - `HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB`: $\text{NES} = +2.0820, \text{FDR } q = 0.0010$
  - `HALLMARK_INFLAMMATORY_RESPONSE`: $\text{NES} = +1.5865, \text{FDR } q = 0.0303$
  - `HALLMARK_ADIPOGENESIS`: $\text{NES} = -3.1463, \text{FDR } q < 0.001$
- **Microenvironment Deconvolution:** Peritoneal Myofibroblasts Mann-Whitney $U = 16.0, P = 0.0286$, exact two-sided $P_{\text{perm}} = 2/70 = 0.0286, \text{FDR} = 0.1714$ (exploratory marker set).

---

## 6. Discussion

### Biological Interpretation
The 11 prioritized candidate hub genes cluster functionally around fibrillar matrix organization (`COL3A1`, `COL11A1`, `COMP`), matricellular and basement membrane remodeling (`COL8A1`, `VCAN`, `THBS3`, `EDIL3`), matrix stabilization/crosslinking (`LOX`, `FN1`), and pro-fibrotic signaling (`INHBA`, `ISM1`). In whole peritoneal biopsy tissue, these genes undergo coordinated, high-magnitude upregulation in severe EPS, consistent with profound stromal myofibroblast activation and mesothelial-to-mesenchymal transition (MMT).

### Biological Discordance: Biopsy Tissue vs. Effluent Cells
When evaluating the signature in shed peritoneal effluent cells, marked differences emerge:
1. **Compartmental Differences:** Effluent pellets comprise shed mesothelial cells, inflammatory leukocytes (macrophages, neutrophils, lymphocytes), and minor sloughed stromal cells, whereas biopsy tissue captures the full-thickness peritoneal membrane.
2. **Matrix Anchorage:** Secreted structural proteins (`LOX`, `COMP`) and cell-adherent proteoglycans (`VCAN`) are primarily deposited within the insoluble tissue matrix, explaining why their cellular transcript abundance in effluent does not mirror tissue elevation.
3. **Clinical Endpoint Distinction:** GSE62928 reflects advanced clinical EPS requiring surgery, whereas GSE125498 compares dialysis vintage ($\le 24$ vs $\ge 25$ months) in stable PD patients without diagnosed EPS. Therefore, GSE125498 tests biological transferability to cumulative dialysis stress rather than diagnostic validation of EPS.

### Strengths of the Study
- Complete transparency and programmatic reproducibility.
- Deterministic analysis without synthetic data or post-hoc over-fitting.
- Exact combinatorial permutation baselines ($\binom{8}{4} = 70$) acknowledging sample size boundaries.
- Clear reporting of apparent vs. cross-validated discrimination metrics.

### Study Limitations
- **Small Discovery Sample Size ($N = 8$):** Discovery transcriptomics is constrained to 8 subjects, limiting statistical power and rendering WGCNA network partitioning exploratory.
- **Unadjusted DEG Filtering:** Only 5 pro-fibrotic ECM-DEGs pass genome-wide FDR $< 0.05$; the broader 71-gene set relies on nominal $P < 0.05$ filtering.
- **Probe Mapping Gaps:** 4 of the 11 hub genes could not be evaluated in the external cohort due to microarray platform coverage differences.
- **Cross-Validation Significance:** The 7-gene panel cross-validated AUC ($0.678$) yields an empirical permutation $P = 0.0879$, demonstrating that external discrimination is moderate and does not achieve statistical significance under rigorous out-of-fold permutation.
- **Post-Hoc ROC Direction:** Directions for `VCAN` and `THBS3` were inverted based on observed effluent downregulation, representing an exploratory post-hoc assessment.

---

## 7. Conclusions

By integrating high-dimensional transcriptomic profiling, curated matrisome filtering, exploratory co-expression networks, and consensus machine learning, we identified 11 candidate matrisome genes associated with peritoneal membrane fibrogenesis. While tissue-level associations are strong and biologically coherent, cross-cohort evaluation in peritoneal effluent cells demonstrates moderate generalization ($\text{AUC}_{\text{CV}} \approx 0.68$). These findings provide prioritized, hypothesis-generating candidate biomarkers that warrant targeted prospective validation in standardized clinical cohorts.
