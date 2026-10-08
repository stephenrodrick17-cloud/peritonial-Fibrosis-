
# PAPER WRITING NOTES — LOCKED MANUSCRIPT NUMBERS & INTERPRETATIONS

> **Effective Threshold:** PRIMARY = |log₂FC| ≥ 0.585 (1.5-fold) AND nominal P < 0.05. The 0.80 cutoff is retained ONLY for explicitly-labelled Historical / Stringent Comparison subsections.

---

## Section 3.1: Discovery Differential Expression (GSE62928)
* **Exact Final Numbers (Primary 0.585):** Gene universe = 22,049; Upregulated DEGs (pro-fibrotic) = 534; Downregulated DEGs = 1,510 (total DEGs = 2,044). Primary thresholds: |log₂FC| ≥ 0.585 (1.5-fold) and nominal P < 0.05. Genome-wide FDR < 0.05 count = 78 genes total (8 up / 70 down). **Historical 0.80 comparison (for reference):** 367 up, 1,263 down, 1,630 total DEGs.
* **Interpretation:** Moderate transcriptomic perturbation in EPS dialysate vs controls, driven by robust pro-fibrotic upregulation.
* **Limitations:** N=8 sample size limit (4 EPS vs 4 Control dialysate samples); strict genome-wide FDR control yields only 78 genes (8 up), necessitating unadjusted P-value filtering with explicit FC thresholds as a discovery-phase strategy. Only 5 of the 81 up ECM-DEGs individually pass genome-wide FDR < 0.05 (BGN, CLEC11A, COL1A1, MXRA5, SERPINE2).
* **Figure/Table Reference:** Figure 1 (Volcano plot & DEG heatmap); Table S1.
* **Wording Restrictions:** Use 'pro-fibrotic candidates' or 'upregulated transcriptomic signal'. DO NOT use 'causal drivers' or 'diagnostic biomarkers'. MUST clearly state the threshold as |log₂FC| ≥ 0.585 (1.5-fold) with nominal P < 0.05.

---

## Section 3.2: Matrisome ECM Enrichment Analysis
* **Exact Final Numbers (Primary 0.585):** Curated Matrisome = 1,027; Expression universe Matrisome = 975; Upregulated ECM-DEGs = 81; Downregulated ECM-DEGs = 67 (total ECM-DEGs = 148). Enrichment ratio = 3.430× (observed 81 vs expected ~23.62); Hypergeometric P = 1.321×10⁻²²; Exact 1-sided permutation P = 0.0143 (rank 1/70 splits). **Historical 0.80 comparison (for reference):** 71 up / 16 down ECM-DEGs, 4.375× enrichment, P ≈ 2.593×10⁻²⁶.
* **Interpretation:** Extremely strong selective enrichment of extracellular matrix components within the pro-fibrotic DEG signal.
* **Limitations:** Database-dependent definition of human matrisome. Of the 81 up ECM-DEGs, only 5 individually survive genome-wide BH-FDR < 0.05. Methodological sensitivity: independent third-party MatrisomeAnalyzeR cross-check on the same raw data returns 125 up ECM-DEGs (Δ=+44) and 998 measured matrisome genes (Δ=+23) due to differing probe-collapse logic, bundled matrisome curation, and internal limma re-fit — the qualitative ECM enrichment signal is robust across implementations, but absolute counts are method-dependent. The repository's locked 81/975/22,049 numbers remain primary for all results.
* **Figure/Table Reference:** Figure 2 (ECM Volcano highlight & enrichment bar chart); Table S2.
* **Wording Restrictions:** Clearly distinguish '81 upregulated ECM-associated genes' from all 534 pro-fibrotic DEGs. State the measured matrisome size as 975 (not 1,027) when referring to the expressed subset.

---

## Section 3.3: Exploratory WGCNA
* **Exact Final Numbers:** WGCNA input = 5,038 genes; 14 co-expression modules; Salmon module = 604 genes; Trait correlation r = 0.806 (P = 0.0157). Label-permutation P = 0.0143 (one-sided, rank 1/70); Family-wise permutation P = 0.0857; Bonferroni P ~ 0.220 across 14 modules.
* **Interpretation:** Salmon module represents an exploratory co-expression cluster associated with EPS phenotype.
* **Limitations:** Not statistically robust after multiple-testing correction across 14 modules (Bonferroni P ≈ 0.220); significance is exploratory.
* **Figure/Table Reference:** Figure 3 (WGCNA module-trait heatmap & dendrogram); Table S3.
* **Wording Restrictions:** MUST use 'exploratory co-expression module'. DO NOT call Salmon a 'validated' or 'statistically robust' module.

---

## Section 3.4: Tripartite Candidate Convergence
* **Exact Final Numbers (Primary 0.585):** 534 Up-DEGs ∩ 975 Measured Matrisome ∩ 604 Salmon Module = 44 convergent genes. **Historical 0.80 comparison (for reference):** 367 ∩ 975 ∩ 604 = 40 convergent genes. All 11 consensus ML hub genes are confirmed present in BOTH the 44-gene (0.585) and 40-gene (0.80) pools.
* **Interpretation:** Strict 3-way convergence identifies core extracellular matrix genes aligned with co-expression networks.
* **Limitations:** Dependent on WGCNA module boundaries and the discovery-threshold choices.
* **Figure/Table Reference:** Figure 4 (Venn diagram & convergent-candidate heatmap); Table S4.
* **Wording Restrictions:** State '975 measured Matrisome genes' explicitly. Report 44 as the primary convergence count; 40 only as historical comparison.

---

## Section 3.5: Multi-Model Machine Learning Prioritization
* **Exact Final Numbers:** 11 consensus genes selected by ≥2 of 4 ML algorithms (LASSO C=4.28, SVM-RFE top 35%, RandomForest Gini≥mean, XGBoost): ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX. Vote-count was computed over the historical 0.80 40-candidate pool, but all 11 hubs are confirmed present in the primary 0.585 44-candidate pool as well.
* **Interpretation:** Robust multi-algorithm feature selection prioritizes key ECM remodeling genes.
* **Limitations:** Small discovery N=8 means ML is purely exploratory for prioritization, not predictive modeling. High co-expression rate among hubs (76.36% of pairs |r|≥0.85) is partly by construction from the tight WGCNA Salmon module.
* **Figure/Table Reference:** Figure 5 (ML selection consensus plot & feature importance); Table 1.
* **Wording Restrictions:** Use 'consensus candidate genes' or 'prioritized candidate genes'. DO NOT use 'validated biomarkers' or 'clinical predictors'.

---

## Section 3.6: Protein-Protein Interaction & Co-expression Architecture
* **Exact Final Numbers:** 42 of 55 hub gene pairs (76.36%) exhibit co-expression |r| ≥ 0.85. STRING v12.5 (retrieved 2026-10-01, score ≥ 0.400) reveals 21 functional edges and 3 physical binding pairs (FN1-LOX, COL11A1-COL3A1, COMP-FN1).
* **Interpretation:** Dense physical and functional connectivity among prioritized ECM candidates.
* **Limitations:** Correlation does not prove direct physical binding for all 42 co-expression pairs; ISM1 has zero STRING edges among the 11 hubs.
* **Figure/Table Reference:** Figure 6 (PPI network & correlation matrix); Table S5.
* **Wording Restrictions:** Cite exact STRING version and date. Report 42/55 pairs = 76.36%.

---

## Section 3.7: Pathway & Gene Set Enrichment Analysis (GSEA)
* **Exact Final Numbers:** 21,597 ranked genes. Top pathways: Hallmark EMT NES = +3.2033 (FDR < 0.001); TNF-α signalling via NF-κB NES = +2.0820 (FDR q = 0.0010); Inflammatory Response NES = +1.5865 (FDR q = 0.0303); Angiogenesis NES = +1.4383 (FDR q = 0.0855).
* **Interpretation:** Unbiased genome-wide ranking confirms Epithelial-Mesenchymal Transition as the dominant biological pathway, coupled with pro-inflammatory and angiogenic programmes.
* **Limitations:** Gene set databases represent general annotations. Note the distinction between GSEA gene-label permutations (1,000) and the exact sample-label 1/70 permutation test.
* **Figure/Table Reference:** Figure 7 (GSEA EMT plot & pathway summary); Table S6.
* **Wording Restrictions:** Specify ranking metric exactly: sign(log₂FC) × (−log₁₀ P).

---

## Section 3.8: Microenvironmental Cellular Deconvolution
* **Exact Final Numbers:** Myofibroblasts Mann-Whitney U = 16, nominal P = 0.0286, exact permutation P = 0.0286, FDR = 0.1714 (across cell types).
* **Interpretation:** Exploratory signal indicating elevated myofibroblast-associated activity in EPS dialysate.
* **Limitations:** Not statistically significant after multiple testing correction across cell types (FDR = 0.1714).
* **Figure/Table Reference:** Figure 8 (Deconvolution boxplots); Table S7.
* **Wording Restrictions:** MUST use 'exploratory myofibroblast-associated signal'. DO NOT claim 'significant expansion' or 'validated cellular enrichment'.

---

## Section 3.9: Cross-Cohort Evaluation (GSE125498)
* **Exact Final Numbers:** N=33 dialysate samples (20 short-term SPD vs 13 long-term PD vintage). 7 of 11 genes mapped (FN1, COL3A1, COL8A1, VCAN, THBS3, LOX, ISM1); 4 unmapped (COL11A1, COMP, EDIL3, INHBA) on GPL10558. Apparent AUC = 0.869; **Primary Repeated 5×5-fold CV AUC = 0.678 (SD 0.058)**; LOOCV AUC = 0.658; 1,000-permutation P = 0.0879 (non-significant). Optimism gap = 0.191.
* **Interpretation:** Modest cross-cohort biological transferability of the 7-gene signature; below conventional significance thresholds for external validation.
* **Limitations:** Score direction for VCAN and THBS3 was inverted post-hoc based on observed external cohort direction; modest generalization performance (CV AUC = 0.678); 1,000-perm P = 0.0879 non-significant. Cohort compares SPD vs LPD effluent vintage, NOT direct EPS surgical biopsy validation.
* **Figure/Table Reference:** Figure 9 (External validation ROC curves & expression boxplots); Table 2.
* **Wording Restrictions:** Report 0.678 as PRIMARY generalization AUC. MUST use 'cross-cohort evaluation' or 'biological transferability'. DO NOT use 'clinical validation' or 'direct EPS validation'.

---

## Section 3.10: Exploratory Multivariable Point-Scoring Model
* **Exact Final Numbers:** 5 genes (VCAN, COL8A1, FN1, ISM1, COL3A1). Apparent C-index = 0.819; **CV AUC = 0.550** (near chance); Brier score = 0.1542; Hosmer-Lemeshow P = 0.2437.
* **Interpretation:** Demonstrates in-sample feasibility but low cross-validated predictive performance — extremely high optimism.
* **Limitations:** Very high optimism gap (CV AUC drops to 0.550 from 0.819 in-sample). NOT suitable for clinical deployment.
* **Figure/Table Reference:** Figure 10 (Nomogram & calibration plot); Table S8.
* **Wording Restrictions:** MUST use 'exploratory multivariable point-scoring model'. DO NOT call it a 'clinical nomogram' or 'diagnostic tool'.

---

## Section 3.11: Single-Cell Cellular Contextualization
* **Status:** SUPPLEMENTARY ONLY (Cellular contextualization, not statistical validation).
* **Wording Restrictions:** DO NOT claim validation of hub genes. Frame as 'putative cellular expression contextualization based on public single-cell atlas data'.

---

## Section 3.12 & 3.13: Regulatory Network Analyses (miRNA / lncRNA)
* **Status:** SUPPLEMENTARY ONLY (Putative regulatory network based on public databases).
* **Wording Restrictions:** DO NOT claim established ceRNA mechanism. Use 'putative regulatory network' throughout.

---

## Appendix A: Mandatory Caveats Checklist (must appear in Discussion / Limitations section)
1. **Small N:** Discovery GSE62928 has only 8 samples (4 EPS / 4 Control) — all significance must be framed as exploratory discovery.
2. **FDR Caveat:** Only 5 of 81 (6.2%) up ECM-DEGs individually pass genome-wide BH-FDR < 0.05. The remaining 76/81 rely on nominal-P-with-FC discovery filtering.
3. **Bonferroni WGCNA:** Salmon module Bonferroni across 14 modules = ~0.220. Not statistically robust after multiple-testing correction.
4. **External Validation:** GSE125498 primary 50×5 CV AUC = 0.678 ± 0.058 with 1,000-perm P = 0.0879 (non-significant). NOT statistically validated externally.
5. **Nomogram CV:** Nomogram 5-fold CV AUC = 0.550 — near chance; exploratory illustration only.
6. **Post-Hoc Direction Inversion:** VCAN and THBS3 score directions were inverted in GSE125498 based on observed external-cohort direction.
7. **Biofluid vs Tissue:** GSE62928 and GSE125498 both profile peritoneal effluent (dialysate) cells, not peritoneum surgical biopsies. Tissue-level extrapolation is inferential.
8. **Matrisome/Probe-Collapse Sensitivity:** Absolute matrisome and ECM-DEG counts are method-dependent: independent MatrisomeAnalyzeR cross-check on the same raw dataset yields 998 measured matrisome (vs 975) and 125 up ECM-DEGs (vs 81), with the qualitative ECM enrichment signal preserved across implementations. Report only the repository's locked 975/81 numbers as primary; cite the cross-check delta in Limitations.
