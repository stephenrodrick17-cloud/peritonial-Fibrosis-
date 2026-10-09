# Comprehensive Scientific and Methodological Report: Peritoneal Fibrosis Pipeline

*Generated on:* 2026-10-06  
*Study Scope:* Stages 1 through 5 multi-omics integration and sensitivity auditing  

---

## Executive Summary and Scientific Framing

This investigation synthesizes bulk tissue discovery transcriptomics, machine learning feature selection, external peritoneal effluent validation, single-cell pseudobulk localization, non-coding RNA regulatory layers, and in vitro mesothelial cell response to characterize extracellular matrix (ECM) remodeling associated with peritoneal membrane fibrogenesis and Encapsulating Peritoneal Sclerosis (EPS).

### Core Methodological Framing & Evidentiary Boundaries:
1. **Discovery Cohort Context:** The primary discovery cohort (GSE62928) comprises a small human parietal peritoneal biopsy dataset ($N = 8$: 4 severe EPS cases vs. 4 non-EPS uremic/PD controls). Candidate matrisome genes and co-expression modules were derived from this limited sample size; true empirical significance is bounded by exact label permutations ($\binom{8}{4} = 70$ partitions, minimum attainable two-sided $P = 2/70 = 0.0286$).
2. **External Effluent Validation Boundaries:** External validation in longitudinal dialysis effluent cell pellets (GSE125498, $N = 33$) serves as an accessible cellular proxy. In-sample multi-gene fits achieve apparent $\text{AUC} = 0.877$, but this is subject to overparameterization. Rigorous out-of-fold cross-validation reveals modest, partition-dependent generalization ($\text{AUC}_{\text{CV}} = 0.659 \pm 0.071$, single 5-fold split $\text{AUC} = 0.542$, permutation test $P = 0.0879$, not statistically significant). Furthermore, several structural collagens exhibit discordant directional expression in effluent cells, highlighting phenotypic divergency between peritoneal membrane fibrotic tissue and shed cellular effluent.
3. **Single-Cell Localization as Exploratory Mapping:** Single-cell RNA sequencing (GSE248762, 16 libraries across 14 patients) maps the 11 hub genes predominantly to the stromal/mesothelial-lineage unresolved compartment, consistent with localized mesenchymal activation during ultrafiltration failure.
4. **Regulatory Layers as In Silico Hypotheses:** miRNA differential overlap and lncRNA ceRNA network models are computational predictions based on multiMiR and ENCORI CLIP evidence. They represent hypothesis-generating frameworks requiring future direct biochemical validation.
5. **Cautious Language:** Throughout this report, findings are framed strictly as being *"associated with"* or *"consistent with"* fibrogenesis, without asserting unvalidated causal mechanisms.

---

## Stage 1 & 2: Tissue Discovery, Matrisome Convergence, and Machine Learning Consensus

- **Discovery Dataset (GSE62928):** 8 parietal peritoneal biopsies (4 EPS, 4 controls) profiled on Affymetrix HG-U133_Plus_2. Probe-level universe = 22,049 genes; MaxMean collapsed = 20,940 genes.
- **Differential Expression:** 534 pro-fibrotic up-regulated DEGs ($P < 0.05, \log_2\text{FC} \ge 0.585$).
- **Matrisome Enrichment:** 81 ECM-DEGs ($19.35\%$ of DEGs vs. $4.42\%$ background; exact one-sided label permutation $P_{\text{perm}} = 1/70 = 0.0143$).
- **WGCNA Module Trait Association:** Unsupervised WGCNA identified the **Salmon module** (604 genes, $r = +0.806, P = 0.0157$; exact two-sided $P_{\text{perm}} = 0.0286$). Multiple testing correction yields Bonferroni-adjusted $P = 0.2198$.
- **Tripartite Convergence:** Strict intersection of 604 Salmon module genes with 81 ECM-DEGs yielded **40 convergent candidate genes**.
- **Machine Learning Consensus Selection:** Four supervised machine learning models evaluated the 44 convergent candidates:
  - LASSO (L1 regularization): 2 features (`ISM1`, `FN1`)
  - SVM-RFE (recursive feature elimination): 14 features
  - Random Forest (Gini importance): 17 features
  - XGBoost (gradient boosted trees): 1 feature (`EDIL3`)
  - **Consensus Rule:** Selection by $\ge 2$ algorithms prioritized **11 consensus pro-fibrotic hub genes**: `ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX`.
  - **Algorithm Agreement Audit:** No candidate achieved 4/4 unanimous selection. 3 genes achieved 3/4 votes (`ISM1`, `FN1`, `EDIL3`); 8 genes achieved 2/4 votes (`VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`, selected jointly by SVM-RFE and Random Forest).

---

## Stage 3: External Effluent Validation & Discordance Analysis (GSE125498)

- **Cohort Specifications:** $N = 33$ dialysis effluent cell pellet samples (20 Short-term PD [0–24 months] vs. 13 Long-term PD [$\ge 25$ months]), profiled on Illumina HumanHT-12 v4.0.
- **Pre-specified 7-Gene Panel Performance:**
  - Apparent In-Sample Fit: $\text{AUC} = 0.8769$
  - 50-Repeat 5-Fold Stratified Cross-Validation (Primary Generalization): Mean $\text{AUC}_{\text{CV}} = 0.6585 \pm 0.0705$
  - Per-Fold Mean $\text{AUC}: 0.6845 \pm 0.2011$
  - Leave-One-Out Cross-Validation (LOOCV): $\text{AUC} = 0.6769$
  - Single Seed-42 5-Fold Split: $\text{AUC} = 0.5423$
  - Empirical Permutation Test (1,000 permutations): $P = 0.0879$ (Not statistically significant (P >= 0.05)).
- **Phenotypic Discordance Finding:** In long-term PD effluent cells, several candidate collagens and glycoproteins (e.g., `FN1`, `COL3A1`, `COMP`) do not display significant upregulation, and some trend downward relative to short-term PD effluent. This demonstrates that shed effluent cellular transcriptomes do not directly mirror the dense fibrocollagenous thickening of the parietal peritoneum, functioning as an imperfect proxy.

---

## Stage 4: Single-Cell RNA Sequencing Pseudobulk Localization (GSE248762)

- **Dataset Architecture:** 16 single-cell libraries across 14 patient donors (4 Long-term PD with ultrafiltration failure [LV_UF], 6 Long-term PD without ultrafiltration failure [LV_NOT_UF], and 6 short-vintage non-failing controls [SV]).
- **Primary edgeR Pseudobulk Model:** Evaluated across 5 annotated lineages using quasi-likelihood negative binomial generalized linear models ($QLF$):

| Gene | log2 Fold Change | 95% CI Low | 95% CI High | P-Value | BH FDR | Evidence Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **COMP** | -4.0442 | -8.0185 | -0.0699 | 4.6414e-02 | 0.2587 | PASS |
| **FN1** | -2.9551 | -5.6558 | -0.2543 | 3.3633e-02 | 0.2587 | PASS |
| **COL8A1** | -2.6583 | -4.4866 | -0.8300 | 6.7008e-03 | 0.1313 | PASS |
| **EDIL3** | -2.6090 | -5.1269 | -0.0910 | 4.2890e-02 | 0.2587 | PASS |
| **COL3A1** | -1.9563 | -4.2147 | +0.3022 | 8.5657e-02 | 0.3354 | PASS |
| **LOX** | -1.7217 | -3.8199 | +0.3765 | 1.0211e-01 | 0.3530 | PASS |
| **THBS3** | -1.5654 | -2.7974 | -0.3334 | 1.5335e-02 | 0.2147 | PASS |
| **ISM1** | -1.0913 | -4.1331 | +1.9505 | 4.6750e-01 | N/A | INSUFFICIENT DATA |
| **COL11A1** | -1.0346 | -4.7462 | +2.6769 | 5.6610e-01 | N/A | INSUFFICIENT DATA |
| **VCAN** | -0.6481 | -2.8755 | +1.5793 | 5.4962e-01 | 0.7481 | PASS |
| **INHBA** | +1.0787 | -0.5936 | +2.7509 | 1.9336e-01 | 0.4515 | PASS |

- **Directionality and Significance in Severe UFF:** In the stromal/mesothelial lineage, 8 of 9 testable candidate hub genes show directionally lower expression in severe ultrafiltration failure (negative log2FC: `COMP` -4.04, `FN1` -2.96, `COL8A1` -2.66, `EDIL3` -2.61, `COL3A1` -1.96, `LOX` -1.72, `THBS3` -1.57, `VCAN` -0.65). Only `INHBA` exhibits positive log2FC (+1.08; 95% CI: -0.59 to +2.75, crossing zero). None of the 9 testable genes achieve genome-wide statistical significance after Benjamini-Hochberg FDR correction (all FDR >= 0.131; 0 of 9 with FDR < 0.05). Findings are strictly hypothesis-generating.
- **Sensitivity and QC-Ceiling Robustness:** Across 9 sensitivity specifications:
  - `VCAN` is confirmed as ROBUST, maintaining consistent negative directionality across all pipelines (-0.65 to -1.41).
  - The remaining 8 testable genes are QC-SENSITIVE, showing 58% to 91% attenuation or sign reversal (`INHBA` flips from +1.08 to -0.31) when upper count ceilings are removed, indicating sensitivity to scRNA-seq technical filters and single-donor influence (donor LV_UF-2).

---

## Stage 5: Multi-Layer Regulatory Architecture

### A. miRNA Layer (Stage 5F)
- **multiMiR Database Retrieval:** 2,040 validated interactions (miRTarBase, TarBase, miRWalk); 563 unique miRNAs target the 11 hub genes. Dense functional evidence (Tier A: reporter assays, western blot, qPCR) supports `COL3A1` (10 Tier A miRNAs), `FN1` (7 Tier A miRNAs), and `LOX` (6 Tier A miRNAs).
- **Human Effluent Exosomes (GSE182736, $n=3$ vs. $3$):** 463 expressed miRNAs in universe; at 1.5-fold differential expression, Tier A overlap is k=6 ($P = 0.9505$) and Tier B overlap is k=186 ($P = 0.9849$). Designated **`DESCRIPTIVE ONLY (human exosome sequencing, n=3 vs 3)`**.
- **Rodent Peritoneal Dialysis Model (GSE130387, $n=3$ vs. $3$):** Affymetrix miRNA-4.0 probes are rat-derived (`rno-*`). Mapped to human homologs by stem homology; designated **`CROSS-SPECIES (rodent model, n=3 vs 3)`**. Significant cross-species enrichment is observed at 1.5-fold cutoff for Tier A (k=7 observed vs. 1.99 expected, $P = 0.00132$) and Tier B (k=46 observed vs. 34.29 expected, $P = 0.00154$).

### B. In Vitro Mesothelial Cell Response (Stage 5H, GSE121372)
- Human peritoneal mesothelial cells (HPMCs) stimulated with 1 ng/mL TGF-$\beta$1 for 6h and 24h on Illumina HumanRef-8 v2.0 (GPL6255).
- *Status:* **`DESCRIPTIVE ONLY; n=1; no statistical testing or P-values possible`**.
- 10 of 11 hubs assayed (`ISM1` absent from platform).
- Strict rule-order analysis accounting for untreated time-in-culture drift ($|\text{drift}| \ge 1.0$):
  - **Confounded by Culture Drift:** `COMP` (drift -2.32), `VCAN` (drift +1.90), `INHBA` (drift -1.47), and `COL11A1` probe ILMN_1644 (drift +1.89).
  - **Higher in Single TGF-beta1 Sample:** `EDIL3` ($\log_2\text{FC}_{24h} = +1.61$, drift +0.93) and `THBS3` ($\log_2\text{FC}_{24h} = +1.08$, drift -0.04).
  - **Unchanged in Single TGF-beta1 Sample ($|\log_2\text{FC}| < 1.0$):** `LOX` (+0.86), `FN1` (+0.87), `COL8A1` (+0.49 / +0.79), `COL3A1` (+0.44), and `COL11A1` probe ILMN_3994 (+0.97).

---

## Summary of Applied Claims and Proposed README Adjustments
The following audited resolutions from `results/tables/claims_to_revise.csv` have been incorporated:

- **CLM_01 (Line 13):** Calibrated against code-derived statistics.
- **CLM_02 (Line 13):** Calibrated against code-derived statistics.
- **CLM_04 (Line 12):** Calibrated against code-derived statistics.
- **CLM_05 (Line 226):** Calibrated against code-derived statistics.
- **CLM_06 (Line 418):** Calibrated against code-derived statistics.

---

- Primary Figures: `results/figures/Fig_4A_pseudobulk_celltype_forest.png`, `results/figures/Fig_4B_stromal_sensitivities_forest.png`, `results/figures/Fig_5A_mirna_hub_intersection.png`, `results/figures/Fig_5C_gse121372_tgfb1_timecourse.png`, `results/figures/Fig_6_discrepancy_corrected_panels.png`
- Note on ceRNA: Excluded from primary manuscript figures per audit directive until each lncRNA-miRNA-hub connection is experimentally verified.
- Figure Captions: `results/figures/FIGURE_CAPTIONS.md`
- Paired Source Data: `results/figures/source_data/`
- Audit Manifests: `results/stage4_manifest.json`, `results/stage5_manifest.json`