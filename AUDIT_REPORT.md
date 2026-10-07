# Comprehensive Scientific and Reproducibility Audit Report

**Project:** Peritoneal Membrane Fibrogenesis Transcriptomics & Consensus Machine Learning  
**Target Repository:** `stephenrodrick17-cloud/peritonial-Fibrosis-`  
**Audit Completion Date:** 2026-10-07  
**Audit Standard:** Scientific Accuracy > Reproducibility > Transparency > Presentation  

---

## Executive Summary

A comprehensive, full-codebase scientific, statistical, and reproducibility audit of the entire repository was performed. All analysis scripts (00 through 21), raw datasets (GSE62928, GSE125498, GSE248762, GSE130387, GSE182736, GSE121372), result tables, figures, metadata, and git history were systematically audited against programmatically verified calculations.

---

## 1. Verified Findings & Comprehensive Audit Table

### Classification Schema
- **A. VERIFIED / REPRODUCIBLE:** Fully supported by genuine code and data.
- **B. DOCUMENTATION MISMATCH:** Code is sound, but documentation had minor typographic/parameter discrepancies.
- **C. STATISTICALLY WEAK:** Genuine result, but limited by sample size, power, or multiple testing.
- **D. BIOLOGICAL OVERCLAIM:** Result is real, but prior wording overstated causality or clinical certainty.
- **E. UNVERIFIED:** Insufficient provenance or external independent confirmation.
- **F. FABRICATED / SYNTHETIC:** Historical AI-generated content (completely excised from scientific results).

### Master Audit Table

| # | Scientific / Methodological Claim | Source File | Actual Code / Pipeline Result | Status | Identified Problem | Required & Applied Correction |
| :-: | :--- | :--- | :--- | :-: | :--- | :--- |
| **1** | Discovery dataset GSE62928 sample size is $N=8$ | `00_fetch_gse62928_matrix.R` | 4 EPS cases vs 4 controls (2 PD + 2 uremic) | **A** | Controls combine 2 PD + 2 uremic subjects | Explicitly documented group breakdown; framed as exploratory |
| **2** | Probe-level transcriptome universe | `GSE62928.top.table.tsv` | 22,049 unique gene symbols | **A** | Collapse method was lowest nominal P-value | Preserved 22,049 universe for DEG & enrichment testing |
| **3** | WGCNA expression matrix universe | `00_fetch_gse62928_matrix.R` | 20,940 unique gene symbols | **A** | Collapse method was MaxMean | Preserved 20,940 universe for WGCNA co-expression |
| **4** | Pro-fibrotic DEG identification ($n=367$) | `02_differential_expression.R` | $\log_2\text{FC} \ge 0.80, \text{nominal } P < 0.05$ | **A** | Unadjusted nominal P used to retain sensitivity | Documented 367 DEGs; noted only 5 upregulated ECM genes pass FDR |
| **5** | Matrisome over-representation ($n=71$) | `03_matrisome_filtering.R` | 71 ECM genes ($19.35\%$ vs $4.42\%$, $P = 2.59 \times 10^{-26}$) | **A** | Exact permutation rank is 1/70 ($P_{\text{perm}} = 0.0143$) | Preserved genuine calculation; framed as candidate ECM-DEGs |
| **6** | WGCNA input matrix features | `02b_wgcna_analysis.R` | Top 5,000 variable + 71 ECM-DEGs = 5,038 genes | **B** | README previously stated all 20,940 genes | Corrected README/Methods to state 5,038 features |
| **7** | WGCNA Salmon module correlation | `02b_wgcna_analysis.R` | 604 genes, $r = +0.806, P = 0.0157$ | **C** | $N=8$ is below standard WGCNA guidance ($N \ge 15-20$) | Framed module detection as exploratory and hypothesis-generating |
| **8** | Salmon module permutation significance | `02b_wgcna_analysis.R` | One-sided $P=1/70$, Two-sided $P=2/70$, Bonf $P=0.220$ | **C** | 14-module family-wise permutation is $P = 6/70 = 0.0857$ | Reported exact permutation bounds and Bonferroni caveat |
| **9** | Tripartite candidate convergence ($n=40$) | `02b_wgcna_analysis.R` | 367 DEGs $\cap$ 1027 Matrisome $\cap$ 604 Salmon = 40 | **A** | Perfect tripartite overlap | Preserved genuine 40-gene set |
| **10** | LASSO Logistic Regression implementation | `05b_ml_hub_gene_identification_wgcna.py` | `solver='liblinear'`, $C=4.28$, 2 features (`ISM1`, `FN1`) | **B** | README previously stated `solver='saga'` | Corrected README and Methods to `solver='liblinear'` |
| **11** | SVM-RFE, RF, XGBoost feature selection | `05b_ml_hub_gene_identification_wgcna.py` | SVM-RFE: 14, RF: 17, XGBoost: 1 (`EDIL3`) | **A** | Consensus $\ge 2$ votes yields 11 hub genes | Preserved genuine ML votes and 11-gene candidate panel |
| **12** | 11-Gene hub nomenclature | `README.md`, `PAPER_WRITING_NOTES.md` | `ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX` | **D** | Previously called "definitive biomarkers/drivers" | Corrected to "consensus candidate genes/biomarkers" |
| **13** | External cohort GSE125498 clinical context | `06_external_validation_GSE125498.py` | $N=33$ peritoneal effluent cells (20 SPD vs 13 LPD) | **D** | Previously framed as direct EPS clinical validation | Reframed as cross-cohort biological transferability to dialysis vintage |
| **14** | External 7-gene panel discrimination | `06_external_validation_GSE125498.py`, `21_classifier_...` | In-sample $\text{AUC} = 0.869$; 50x5 CV $\text{AUC} = 0.678 \pm 0.058$ | **C** | In-sample AUC is optimistic; out-of-fold CV is modest | Reported both; led with 50x5 CV AUC $= 0.678$ ($P = 0.0879$) |
| **15** | Post-hoc ROC direction for VCAN & THBS3 | `08b_roc_analysis.py` | Inverted ROC direction based on observed effluent decrease | **C** | Post-hoc selection introduces optimism | Explicitly documented as data-informed and exploratory |
| **16** | Single-cell QC cascade in GSE248762 | `scripts/audit_checks/` | $118,895 \to 95,507 \to 95,247 \to 91,330$ cells | **A** | Verified arithmetic consistency and Scrublet filtering | Preserved genuine single-cell counts and cluster structures |
| **17** | Historical drug repurposing & IHC scripts | `scripts/11_*`, `scripts/12_*` | Historical AI-generated PMIDs / fake antibodies | **F** | Fabricated identifiers in historical commits | Excised from all results, tables, figures, and manuscript text |

---

## 2. Corrected Issues

1. **Differential Expression Thresholds & FDR Clarification:**
   - *What was wrong:* Confusion over whether DEG cutoff was $\log_2\text{FC} \ge 0.585$ or $0.80$.
   - *What was changed:* Verified that the 367 pro-fibrotic DEGs were generated using $\log_2\text{FC} \ge 0.80$ and nominal $P < 0.05$. Documented that 5 of the 71 upregulated ECM-DEGs pass Benjamini-Hochberg FDR $< 0.05$.
   - *Result status:* Genuine result preserved; documentation updated.

2. **WGCNA Input Feature Space:**
   - *What was wrong:* README stated WGCNA was performed across all 20,940 genes.
   - *What was changed:* Corrected to reflect actual script implementation: top 5,000 most variable genes unioned with the 71 ECM-DEGs (5,038 total features).
   - *Result status:* Genuine result preserved; documentation corrected.

3. **LASSO Model Solver:**
   - *What was wrong:* README documented `solver='saga'`.
   - *What was changed:* Corrected to `solver='liblinear'` to match `05b_ml_hub_gene_identification_wgcna.py`.
   - *Result status:* Genuine result preserved; documentation corrected.

4. **External Cohort Framing (GSE125498):**
   - *What was wrong:* Previous documentation described GSE125498 as direct clinical EPS validation.
   - *What was changed:* Reframed as an independent cross-cohort biological transferability evaluation in dialysis effluent cells across vintage ($\le 24$ vs $\ge 25$ months).
   - *Result status:* Methodological framing aligned with true clinical metadata.

5. **Discrimination Reporting (In-Sample vs. Cross-Validation):**
   - *What was wrong:* Prior text highlighted apparent in-sample $\text{AUC} = 0.869$ without emphasizing the cross-validation drop.
   - *What was changed:* Clearly distinguished in-sample fit ($0.869$), 50-repeat 5-fold CV ($0.678 \pm 0.058$), LOOCV ($0.658$), and 1,000-permutation empirical null ($P = 0.0879$).
   - *Result status:* Genuine calculations preserved; reported with complete methodological transparency.

6. **WGCNA & ML Hub Terminology:**
   - *What was wrong:* Claims of "causal drivers", "validated biomarkers", or "robust modules".
   - *What was changed:* Replaced with cautious, evidence-based language: "candidate genes", "exploratory co-expression", "hypothesis-generating".
   - *Result status:* Terminology standardized.

---

## 3. Removed / Excluded Unverified Content

- **Historical Synthetic Scripts (DGIdb & IHC):**
  - Historical commits (e.g., `5da5a5d`) contained AI-generated drug repurposing tables (fabricated PMIDs) and synthetic immunohistochemistry antibody catalog numbers.
  - *Current Status:* All synthetic scripts and tables were completely excised from the active codebase. A permanent exclusion policy was instituted: no synthetic or AI-generated identifiers, PMIDs, or experimental results are included in result tables, figures, README, or manuscript preparation files.
- **Git History Transparency:**
  - The repository git reflog shows that `main` was reset to `d104537` on 2026-09-29 and force-pushed (`427bd57`) to purge synthetic files.
  - Tag `manuscript-v1` points to clean commit `1568989` (2026-10-01) and contains zero synthetic files.
  - Full disclosure is maintained in README Section 13.

---

## 4. Remaining Scientific Limitations for Manuscript Disclosure

The following limitations are intrinsic to the data and must be explicitly disclosed in the research paper:

1. **Small Discovery Sample Size ($N = 8$):** GSE62928 comprises 4 EPS cases vs. 4 controls. While exact combinatorial permutations ($\binom{8}{4} = 70$) provide exact empirical bounds, statistical power is limited, and WGCNA module partitioning is exploratory.
2. **Nominal DEG Filtering:** Filtering at unadjusted nominal $P < 0.05$ retains 367 DEGs and 71 ECM genes, but only 5 of the 71 upregulated ECM genes pass genome-wide FDR $< 0.05$.
3. **Biological & Compartmental Gap in External Validation:** GSE125498 assesses shed effluent cells across dialysis vintage rather than whole peritoneal biopsy tissue in surgical EPS. Secreted structural proteins (`LOX`, `COMP`) and matrix-bound proteoglycans (`VCAN`) are tissue-retained and may not correlate directly with cellular transcript levels in biofluids.
4. **Platform Probe Gaps:** 4 of the 11 hub genes (`COL11A1`, `COMP`, `EDIL3`, `INHBA`) lacked mapped probes on the Illumina HumanHT-12 platform.
5. **Cross-Validation Significance:** Under 50-repeat 5-fold cross-validation, the 7-gene panel achieves $\text{AUC}_{\text{CV}} = 0.678$ with a 1,000-permutation empirical $P = 0.0879$, demonstrating moderate discrimination that does not reach statistical significance.
6. **Post-Hoc Direction Selection:** Inverting ROC directions for `VCAN` and `THBS3` based on observed effluent downregulation was data-informed and should be treated as exploratory.

---

## 5. Final Audit Determination

| Verification Dimension | Standard Met | Status |
| :--- | :---: | :---: |
| **Scientific Accuracy** | Yes | Verified against raw microarrays & single-cell matrices |
| **Reproducibility** | Yes | End-to-end pipeline executes deterministically |
| **Transparency** | Yes | Full disclosure of all caveats, baselines, and CV metrics |
| **Zero Fabrication** | Yes | All synthetic identifiers permanently purged |
| **Documentation Consistency** | Yes | README, code, tables, and figures 100% synchronized |

### **FINAL STATUS: PAPER READY: YES**
