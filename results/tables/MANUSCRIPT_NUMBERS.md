# Master Manuscript Numbers Sheet (Reference Summary)

> **Compilation & Provenance Note:**
> This reference summary was compiled by hand from the listed primary pipeline output CSVs, sensitivity scripts, and literature sources, and formatted via `scripts/audit_checks/generate_manuscript_numbers.py`. It is **not** an independently parsed computational audit output. While individual numbers reflect the cited pipeline outputs, the descriptions and tables were manually structured and maintained. Readers and reviewers should directly consult the primary source CSVs and scripts cited under each section.
Each section cites the exact primary source file and execution command.

## 1. Clinical Cohort Characteristics (GSE248762, N = 16 Donors)
**Source Files:** `data/raw/GSE248762_series_matrix.txt.gz`, `results/tables/pb_primary_metadata.csv`
**Command:** `python -c "import pandas as pd; ..."`

| Group Name | Definition / Phenotype | Donors ($N$) | Vintage Median (Months) | Vintage Range | Peritoneal UF Status | Primary Role |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **`LV_UF`** | Long-vintage with ultrafiltration failure | 4 | 132 | 96 – 168 | $< 400$ mL / 4h | Target pathological group |
| **`LV_NOT_UF`** | Long-vintage with preserved ultrafiltration | 6 | 137 | 96 – 180 | $> 400$ mL / 4h | Matched dialyzed control |
| **`SV`** | Short-vintage dialyzed baseline | 6 | 5 | 1 – 6 | Variable | Vintage-baseline comparator |

*Note: SV contrasts evaluate dialysis vintage; only `LV_UF vs LV_NOT_UF` directly evaluates ultrafiltration failure at matched vintage.*

## 2. Single-Cell Barcode Filtering Hierarchy & Doublet Benchmark
**Source Files:** `data/processed/GSE248762_hubblind_allcells_qc.h5ad`, `data/processed/GSE248762_harmony_annotated_obs.h5ad`, `scripts/29_scrna_doublets_qc_rerun.py`, `scripts/26_scrna_doublet_investigation.py`
**Command:** `python scripts/29_scrna_doublets_qc_rerun.py`

| Method / Stage | Cutoff / Threshold | Flagged (All 118,895 Barcodes) | Flagged (QC-Pass, N=95,507) | Doublet Calls in Final Object (N=95,247) | Used for Primary Exclusion |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Scrublet (Primary Global)** | Cutoff = 0.439 (pooled simulated KDE valley) | 320 (0.27%) | 260 (0.27%) | 0 (0.00%) | **YES** (removed 260 cells from 95,507 to produce 95,247 cells) |
| **scDblFinder** | Cluster-aware dynamic classifier | 10,057 (8.46%) | 8,355 (8.75%) | 8,112 (8.52%) | **NO** (evaluated in sensitivity `sens_c_no_scDblFinder`)* |
| **Legacy Heuristic (Script 14)** | Count/gene percentiles | 4,672 (3.93%) | 977 (1.02%) | 977 (1.03%) | **NO** (superseded heuristic) |

*Barcode Step Accounting:*
- All sequenced barcodes: **118,895**
- Barcodes passing initial library/mito QC filters: **95,507**
- Scrublet predicted doublets among QC-pass cells: **260** (of which 243 were also flagged by scDblFinder)
- Final analysis object: **95,247** ($95,507 - 260 = 95,247$). Exactly 260 Scrublet doublets were removed.
- Exactly **8,112** scDblFinder doublet calls remain in the 95,247-cell object ($8,355 - 243 = 8,112$) because scDblFinder was not used for primary cell exclusion.
- An scDblFinder-removed sensitivity model (`sens_c_no_scDblFinder` in `results/tables/stage4_sensitivity_edger_pseudobulk.csv`) confirmed similar fold-change estimates.

## 3. Cell-Type Demographics & Contrast-Specific Testability
**Source File:** `data/processed/GSE248762_harmony_annotated_obs.h5ad`
**Command:** `python scripts/audit_checks/check_a2_stromal_counts.py`

| Cell Type Annotation | Total Cells | % of Object | Donors >= 20 Cells (LV_UF) | Donors >= 20 Cells (LV_NOT_UF) | Donors >= 20 Cells (SV) | Testable in LV_UF vs LV_NOT_UF? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **T cell** | 44,245 | 46.45% | 4 / 4 | 6 / 6 | 6 / 6 | **YES** (4 vs 6 donors) |
| **Monocyte / macrophage** | 22,925 | 24.07% | 3 / 4 | 6 / 6 | 6 / 6 | **YES** (3 vs 6 donors) |
| **cDC** | 13,672 | 14.35% | 3 / 4 | 6 / 6 | 6 / 6 | **YES** (3 vs 6 donors) |
| **NK cell** | 4,011 | 4.21% | 4 / 4 | 6 / 6 | 6 / 6 | **YES** (4 vs 6 donors) |
| **Stromal / mesothelial-lineage** | 3,670 | 3.85% | 4 / 4 | 5 / 6 | 5 / 6 | **YES** (4 vs 5 donors; 4 vs 6 in edgeR) |
| **Neutrophil** | 3,647 | 3.83% | 2 / 4 | 3 / 6 | 2 / 6 | **NO** (< 3 donors in LV_UF) |
| **B cell** | 1,840 | 1.93% | 2 / 4 | 3 / 6 | 4 / 6 | **NO** (< 3 donors in LV_UF) |
| **Mixed: pDC / Plasma cell** | 1,237 | 1.30% | 3 / 4 | 6 / 6 | 6 / 6 | **YES** (3 vs 6 donors) |

*Rules applied: Pipeline execution rule: sum(counts) >= 10 in >= 2 donors per group. Manuscript-level testability criterion: >= 20 cells per donor in >= 3 donors per group.*
*Pipeline primary pseudobulk tested 5 lineages (TEST_CTS: T cell, Monocyte/macrophage, cDC, NK cell, Stromal; n = 88,523 cells, 92.94% of object). Stromal n = 3,670 cells. Total cells = 95,247 across 8 distinct lineages (no separate Mast cell).*

## 4. Stromal Heterogeneity & Lineage Scoring
**Source Files:** `scripts/audit_checks/check_a4_stromal_heterogeneity.py`, `data/processed/GSE248762_harmony_annotated_obs.h5ad`
**Command:** `python scripts/audit_checks/check_a4_stromal_heterogeneity.py`

| Group | Donors ($N$) | Mean Mesothelial Fraction | Mean UMI Count per Donor | Statistical Test (Between Groups) |
| :--- | :---: | :---: | :---: | :--- |
| **`LV_UF`** | 4 | 60.10% ± 25.10% | 13,008.2 ± 3,901.8 | Kruskal-Wallis: Mesothelial P = 0.1134, UMI P = 0.2115 |
| **`LV_NOT_UF`** | 6 | 82.55% ± 12.35% | 10,751.2 ± 3,463.7 | Kruskal-Wallis: Mesothelial P = 0.1134, UMI P = 0.2115 |
| **`SV`** | 6 | 78.96% ± 9.94% | 17,006.9 ± 10,397.9 | Kruskal-Wallis: Mesothelial P = 0.1134, UMI P = 0.2115 |

*Assessment: No detectable difference; underpowered (4 vs 6 donors); not evidence of equivalence. This lack of difference is not used to justify pooling.*

## 5. Donor-Dominated Cluster Composition
**Source File:** `data/processed/GSE248762_harmony_annotated_obs.h5ad`
**Command:** `python -c "import anndata as ad; ..."`

| Cluster | Dominant Donor | Donor Contribution (%) | Other Donors Contribution (%) | Analytical Action |
| :--- | :--- | :---: | :---: | :--- |
| **Cluster 4** ($n = 2,647$) | `LV_UF-3` | **82.62%** | 17.38% (15 donors) | Excluded from primary pseudobulk |
| **Cluster 12** ($n = 1,270$) | `LV_UF-3` | **100.00%** | 0.00% (0 donors) | Excluded from primary pseudobulk |

## 6. Primary Stromal Pseudobulk Results (LV_UF vs LV_NOT_UF, m = 98 Family)
**Source Files:** `results/tables/stage4_primary_edger_pseudobulk.csv`, `results/tables/stage4_primary_98_bh_family_sorted.csv`
**Command:** `python scripts/audit_checks/check_directive2.py`

| Gene | log2FC (LV_UF vs LV_NOT_UF) | SE | P-Value | BH FDR ($m=98$) | Descriptive Evidence Status | Concordance vs Biopsy Discovery |
| :--- | :---: | :---: | :---: | :---: | :--- | :--- |
| **`VCAN`** | -0.6481 | 1.0637 | 0.5496 | 0.7481 | PASS | Opposite (Negative in Stromal) |
| **`EDIL3`** | -2.6090 | 1.2163 | 0.0429 | 0.2587 | PASS | Opposite (Negative in Stromal) |
| **`COL8A1`** | -2.6583 | 0.8732 | 0.0067 | 0.1313 | PASS | Opposite (Negative in Stromal) |
| **`THBS3`** | -1.5654 | 0.5909 | 0.0153 | 0.2147 | PASS | Opposite (Negative in Stromal) |
| **`FN1`** | -2.9551 | 1.2900 | 0.0336 | 0.2587 | PASS | Opposite (Negative in Stromal) |
| **`COMP`** | -4.0442 | 1.9223 | 0.0464 | 0.2587 | PASS | Opposite (Negative in Stromal) |
| **`COL3A1`** | -1.9563 | 1.0789 | 0.0857 | 0.3354 | PASS | Opposite (Negative in Stromal) |
| **`LOX`** | -1.7217 | 1.0020 | 0.1021 | 0.3530 | PASS | Opposite (Negative in Stromal) |
| **`INHBA`** | 1.0787 | 0.8011 | 0.1934 | 0.4515 | PASS | Concordant (Positive in Stromal) |
| **`COL11A1`** | -1.0346 | 1.7711 | 0.5661 | NA | INSUFFICIENT DATA | Opposite (Negative in Stromal) |
| **`ISM1`** | -1.0913 | 1.4801 | 0.4675 | NA | INSUFFICIENT DATA | Opposite (Negative in Stromal) |

**Key Directional Finding:** In the stromal/mesothelial lineage, **10 of 11 hub genes have negative log2FC** in `LV_UF vs LV_NOT_UF` (8 of 9 testable: VCAN -0.65, EDIL3 -2.61, COL8A1 -2.66, THBS3 -1.57, FN1 -2.96, COMP -4.04, COL3A1 -1.96, LOX -1.72; untestable COL11A1 -1.03, ISM1 -1.09). Only `INHBA` is positive (+1.08). This direction is **opposite in direction to the biopsy contrast (different compartment and groups; not directly comparable)**.

*Note: 'PASS' denotes descriptive testability under the analysis plan (self-documented; committed after primary results) criterion (>= 2 expressing donors per group), NOT statistical significance. None of the 9 testable genes achieve family-wise FDR < 0.05 (all stromal FDR >= 0.131). Findings are exploratory.*

## 7. Stromal Robustness & Sensitivity Suite Across Specifications
**Source File:** `results/tables/stage4_sensitivity_edger_pseudobulk.csv`
**Command:** `python scripts/audit_checks/check_directive2.py`

### A. Selected Sensitivity Models
| Model Specification | LV_UF Donors | LV_NOT_UF Donors | VCAN log2FC (P) | EDIL3 log2FC (P) | FN1 log2FC (P) | THBS3 log2FC (P) | COL8A1 log2FC (P) | Interpretation Role |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Primary** (all clean donors) | 4 | 6 | -0.6481 (0.5496) | -2.6090 (0.0429) | -2.9551 (0.0336) | -1.5654 (0.0153) | -2.6583 (0.0067) | Primary inference |
| **Sens (b):** Drop LV_UF-3 | 3 | 6 | -2.3043 (0.0869) | -3.6685 (0.0154) | -3.4811 (0.0475) | -1.8370 (0.0098) | -2.5286 (0.0240) | Outlier donor check |
| **Sens (e):** >= 50 cells/donor | 3 | 2 | -1.7412 | -3.7840 | -4.1877 | -2.2457 | -4.5844 | **3 vs 2 donors; direction only** (P-value omitted) |
| **Sens (c):** Remove scDblFinder | 4 | 6 | -0.4759 (0.6877) | -2.7029 (0.0387) | -2.8668 (0.0432) | -1.5593 (0.0149) | -2.6933 (0.0059) | Doublet sensitivity |

### B. Stromal Composition edgeR Re-runs (Directive 4; Post Hoc, Exploratory)
**Source Scripts:** `scripts/audit_checks/run_stromal_composition_tests.py`, `scripts/audit_checks/run_stromal_composition_edger.R`, `scripts/audit_checks/inspect_stromal_run_a.R`

| Gene | Run (a): Meso Fraction Covariate (N=16) log2FC | Run (a) edgeR P-value | Run (b1): Meso-Scored Subset log2FC | Run (b1) edgeR P-value | Run (b2): Fibro-Scored Subset log2FC | Run (b2) edgeR P-value | Direction Assessment |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **`COL8A1`** | -2.90 | 7.76e-05 | -4.00 | 0.00124 | -2.72 | 0.00021 | Negative in all 3 runs |
| **`FN1`** | -2.19 | 0.0508 | +0.17 | 0.8406 | -4.52 | 0.00022 | Negative in runs a & b2; positive in b1 |
| **`COL3A1`** | -1.88 | 0.0404 | -3.75 | 0.00156 | -2.27 | 0.00976 | Negative in all 3 runs |
| **`THBS3`** | -1.35 | 0.0220 | -1.67 | 0.0139 | -1.81 | 0.00689 | Negative in all 3 runs |
| **`LOX`** | -1.73 | 0.0386 | -5.03 | 0.00005 | -1.99 | 0.00954 | Negative in all 3 runs |
| **`EDIL3`** | -1.49 | 0.0863 | -2.94 | 0.0512 | -1.86 | 0.1075 | Negative in all 3 runs |
| **`VCAN`** | -0.56 | 0.5266 | +2.01 | 0.1146 | -2.05 | 0.0599 | Negative in runs a & b2; positive in b1 |
| **`COMP`** | -1.45 | 0.4630 | -3.02 | 0.6563 | -2.46 | 0.2284 | Negative in all 3 runs |
| **`COL11A1`** | -0.87 | 0.6676 | -4.58 | 0.1344 | +2.21 | 0.4284 | Negative in runs a & b1; positive in b2 |
| **`INHBA`** | +1.26 | 0.0999 | -2.38 | 0.00147 | +2.80 | 0.0819 | Positive in runs a & b2; negative in b1 |
| **`ISM1`** | +0.26 | 0.8583 | -4.10 | 0.2933 | +1.10 | 0.5254 | Positive in runs a & b2; negative in b1 |

*Directionality and Robustness Assessment:*
Direction is mostly retained in small subsets (3 vs 4 and 3 vs 2 donors); INHBA, ISM1, FN1, VCAN and COL11A1 flip sign in at least one run (5 of 11). Exploratory, no FDR.


## 8. External Validation Cohort Probes & Covariate Analysis (GSE125498, N = 33: 13 Late-Stage LPD vs 20 Early-Stage SPD)
**Source Files:** `data/raw/GPL10558.annot.gz`, `results/tables/C_gse125498_concordance.csv`, `results/tables/roc_auc_detailed_metrics.csv`, `Validation/GSE125498.top.table.tsv`, `scripts/audit_checks/reconcile_p_values.py`
**Command:** `python scripts/audit_checks/reconcile_p_values.py`

| Gene | Probe ID Used | Probe Role | log2FC (LPD vs SPD) | Welch t-test P | Limma Moderated P (Primary) | OLS Unadjusted P | Notes / Interpretation |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **`VCAN`** | `ILMN_1687301` | Single probe | -0.5224 | 0.0315 | **0.0244** | 0.0282 | Significant raw decrease in late PD effluent |
| **`FN1`** | `ILMN_2366463` | Highest-mean | +0.3171 | 0.3278 | **0.3190** | 0.3409 | Non-significant upward trend |
| **`FN1`** | `ILMN_1778237` | Lowest-P | +0.4068 | 0.2961 | **0.2690** | 0.2939 | Non-significant upward trend |
| **`COL8A1`** | `ILMN_1685433` | Highest-mean | +0.3689 | 0.4185 | **0.3920** | 0.4183 | The other probe showed log2FC +0.37, P ≈ 0.4 |
| **`COL8A1`** | `ILMN_2402392` | Lowest-P | +0.7493 | 0.0708 | **0.0488** | 0.0615 | Nominal with lowest-P probe under limma; Bonferroni P = 0.0976 |
| **`COL3A1`** | `ILMN_1773079` | Single probe | +0.1869 | 0.6210 | **0.6160** | 0.6200 | Flat / non-significant |
| **`THBS3`** | `ILMN_1804663` | Single probe | -0.2008 | 0.4053 | **0.4180** | 0.4036 | Non-significant decrease |
| **`LOX`** | `ILMN_1695880` | Single probe | +0.0165 | 0.9781 | **0.9760** | 0.9781 | Inactive in shed cells |
| **`ISM1`** | `ILMN_3239288` | Single probe | +0.0282 | 0.9102 | **0.9040** | 0.9101 | Flat |
| **`COL11A1, COMP, EDIL3, INHBA`** | *None* | Unmapped | — | — | — | — | Absent from Illumina HT-12 v4 array |

### B. Leukocyte Marker Covariate Adjustment Suite for VCAN (Post Hoc, Exploratory)
| Model Specification | CD14 Probe | FCGR3B Probe | CD3E Probe | VCAN log2FC | Limma Moderated P (Primary) | OLS t-statistic | OLS P-Value | Assessment |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **Unadjusted** | — | — | — | -0.5224 | **0.0244** | -2.3510 | 0.0282 | Significant raw decrease in late PD |
| **Adjust CD14 only (Set A)** | `ILMN_2396444` | — | — | -0.2971 | **0.1937** | -1.3258 | 0.1937 | Attenuated |
| **Adjust CD14 only (Set B)** | `ILMN_1740015` | — | — | -0.2845 | **0.1429** | -1.4992 | 0.1429 | Attenuated |
| **Adjust CD3E only** | — | — | `ILMN_1739794` | -0.3241 | **0.1955** | -1.3201 | 0.1955 | Attenuated |
| **Adjust FCGR3B only (Set A)** | — | `ILMN_2134453` | — | -0.6728 | **0.0067** | -2.8871 | 0.0067 | Significant (P = 0.0067) |
| **Adjust FCGR3B only (Set B)** | — | `ILMN_1728639` | — | -0.7043 | **0.0026** | -3.2403 | 0.0026 | Significant (P = 0.0026) |
| **Joint Model (Set A)** | `ILMN_2396444` | `ILMN_2134453` | `ILMN_1739794` | -0.3009 | **0.1944** | -1.3144 | 0.1994 | Attenuated (limma P = 0.1944, OLS P = 0.1994) |
| **Joint Model (Set B)** | `ILMN_1740015` | `ILMN_1728639` | `ILMN_1739794` | -0.0902 | **0.6716** | -0.4360 | 0.6662 | Abolished (limma P = 0.6716, OLS P = 0.6662) |

*Biological Interpretation & Caveats:*
- VCAN's decrease in long-term PD effluent cells is fragile and depends on probe set and adjustment method: attenuated or abolished under CD14/CD3E/FCGR3B marker adjustment, preserved under 5-lineage-score adjustment. A composition effect is one hypothesis for this fragility.
- Probe Cross-Reactivity Note: GPL10558 annotation cannot resolve cross-reactivity between FCGR3A and FCGR3B; the two paralogs share ~97% nucleotide identity. Pairwise local alignment (generated by scripts/audit_checks/run_fcgr3_probe_alignment.py and saved in results/audit/fcgr3_alignment_raw.txt) reveals that FCGR3A probe ILMN_2112580 is the clearest cross-reactive probe, having 0 mismatches to FCGR3A (NM_000569.8) and only 1 mismatch against off-target FCGR3B (NM_000570.5) over 50 bp (score 50 vs 48). Set A's FCGR3B probe ILMN_2134453 fits FCGR3B (NM_000570.5 at 1928–1978) with 1 gap (0 mismatches across 50 bp), but also aligns to FCGR3A at 43/45 bp with 1 mismatch and 1 gap (score 40 vs 48), so it cannot be assumed gene-specific. FCGR3B probe ILMN_1728639 has 2 mismatches to FCGR3B and 2 mismatches + 1 gap to FCGR3A. For ILMN_1703679, no valid alignment found (score <= 15.0); excluded from cross-reactivity claims.
- Set B's CD14 probe (ILMN_1740015) has wider dynamic range and stronger correlation with VCAN (r = 0.66 vs 0.51 for ILMN_2396444). In the leave-one-out covariate analysis for Set B, the attenuation is not driven by CD14 alone but depends on CD14 and FCGR3B together: dropping CD14 alone restores nominal significance (log2FC = -0.4930, OLS P = 0.0429, Limma P = 0.0395), whereas dropping FCGR3B alone completely eliminates the effect (log2FC = +0.0178, OLS P = 0.9263 / Limma P = 0.9276) and dropping CD3E yields OLS P = 0.0745 / Limma P = 0.0736. VCAN remains fragile across covariate adjustments.
- Over-adjustment caveat: Because shed leukocyte subsets covary with clinical dialysis duration, adjusting for multiple cell markers in an N=33 cohort carries risk of collinearity and potential over-adjustment.
- Note on AUC Direction: In the external cohort ROC analysis, AUC direction for VCAN and THBS3 was inverted post hoc after observing that they were downregulated in late-stage PD effluent; reported discrimination estimates are optimistic.

### C. Effluent Lineage Signatures from GSE248762 & Adjusted Regressions (Post Hoc, Exploratory)
**Source Script:** `scripts/audit_checks/check_directive5_lineage_covariates.py`

#### 1. Correlation of Lineage Scores with Group (LPD vs SPD)
| Lineage Signature (Top 20 Markers on GPL10558) | Pearson Correlation ($r$) | P-Value | Correlation Direction |
| :--- | :---: | :---: | :--- |
| **T cell** | +0.457 | 0.0075 | Significant positive correlation with late PD |
| **NK cell** | +0.509 | 0.0025 | Significant positive correlation with late PD |
| **Monocyte / macrophage** | -0.281 | 0.1133 | Negative trend in late PD |
| **cDC** | -0.013 | 0.9442 | No correlation |
| **Neutrophil** | +0.129 | 0.4735 | Weak positive trend |

#### 2. VCAN, FN1, and COL8A1 Regressions Adjusted for All 5 Lineage Scores (OLS)
*Method: Ordinary Least Squares (OLS) multivariable linear regression modeling probe expression ~ group + 5 lineage z-scores.*

| Target Gene / Probe | Unadjusted OLS log2FC | Unadjusted OLS P | Adjusted OLS log2FC | Adjusted OLS t-statistic | Adjusted OLS P-Value | Impact of Lineage Adjustment |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **`VCAN`** (`ILMN_1687301`) | -0.5224 | 0.0282 | -0.4872 | -2.3717 | 0.0254 | Preserved negative direction and nominal significance |
| **`FN1 (MaxMean)`** (`ILMN_2366463`) | +0.3171 | 0.3409 | +0.3735 | +1.0183 | 0.3179 | Preserved positive trend; non-significant |
| **`FN1 (Lowest-P)`** (`ILMN_1778237`) | +0.4068 | 0.2939 | +0.7076 | +1.5874 | 0.1245 | Upregulation magnitude increased; non-significant |
| **`COL8A1 (MaxMean)`** (`ILMN_1685433`) | +0.3689 | 0.4183 | +0.1143 | +0.1997 | 0.8433 | Substantially attenuated |
| **`COL8A1 (Lowest-P)`** (`ILMN_2402392`) | +0.7493 | 0.0615 | +0.5093 | +0.9797 | 0.3362 | limma P = 0.0488 (lowest-P probe); OLS unadjusted 0.0615; adjusted 0.336 |

## 9. Breakdown of 699 Unannotated ENSG Features
**Source Files:** `results/tables/check_a6_699_ensg_biotypes.csv`, `scripts/audit_checks/check_a6_biotypes.py`
**Command:** `python scripts/audit_checks/check_a6_biotypes.py`

| Biotype (Current Ensembl Release) | Feature Count | Percentage | Description / Category |
| :--- | :---: | :---: | :--- |
| **Protein-Coding** | 578 | 82.69% | Genuine protein-coding genes without approved HGNC symbols in reference |
| **lncRNA** | 77 | 11.02% | Long non-coding RNAs not captured in GENCODE v32 lncRNA catalog |
| **Unresolved / Unknown** | 30 | 4.29% | Retired identifiers or novel assemblies |
| **Pseudogene** | 7 | 1.00% | 3 processed, 3 IG_V, 1 transcribed processed pseudogenes |
| **Artifact** | 5 | 0.72% | Genomic sequencing or mapping artifacts |
| **IG_V_gene** | 2 | 0.29% | Immunoglobulin variable segment genes |
| **Total** | **699** | **100.0%** | Biotypes derived from the current Ensembl release, not GENCODE v32 |

## 10. ceRNA Regulatory Layer & Non-Coding RNA Architecture (Predicted, Hypothesis-Generating)
**Source Files:** `results/tables/G_ceRNA_network_predicted.csv`, `results/tables/G_ceRNA_network_summary.csv`, `results/tables/F_hub_to_mirna_tiered.csv`, `results/tables/G_lncRNA_threshold_sensitivity.csv`

- **Inference Constraints:** All ceRNA networks and axes are **predicted, hypothesis-generating** models without prospective wet-lab validation in this system.
- **Network Dimensions:** The predicted network comprises **3,397 candidate tripartite axes** spanning **71 unique lncRNAs**, **10 intermediate miRNAs**, and **11 candidate hub targets**.
- **Literature Selection of miRNAs:** The initial network construction queried **10 literature-selected miRNAs** (`let-7a-5p`, `miR-29a-3p`, `miR-29b-3p`, `miR-29c-3p`, `miR-30a-5p`, `miR-21-5p`, `let-7b-5p`, `miR-34a-5p`, `miR-17-5p`, `miR-93-5p`), selected based on published fibrosis and EMT literature. An objective ranking based strictly on Tier A functional targets is cataloged in `G_ceRNA_alternative_top10_tierA.csv`.
- **Target Evidence Multiplicity:** In multiMiR experimental records, **only 4 of the 11 hub genes have Tier A functional evidence** (reporter assay, western blot, qPCR: `COL3A1` [10 miRNAs], `FN1` [7], `LOX` [6], `VCAN` [1]). The remaining 7 hubs (`ISM1`, `EDIL3`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`) possess **0 Tier A functional miRNAs** (Tier B high-throughput CLIP/sequencing or bioinformatic predictions only).
- **Statistical Rigor & Lack of FDR Support:** In single-cell stromal pseudobulk differential expression across donors, none of the hub genes achieve FDR significance (all stromal hub DE FDR $\ge 0.131$). At the axis level, donor co-expression is nominal-only without FDR support (among robust candidates: 75 axes with nominal log2-CPM Spearman $\rho > 0, P < 0.05$; 40 axes with nominal partial correlation $r > 0, P < 0.05$ controlling for donor cell counts and library depth).
- **Evidence Classification Breakdown:** Across the 3,397 predicted axes, 2,143 are classified as `NOMINAL_COEXPRESSION_AND_MIRNA_DETECTED`, 916 as `LOW_CREDIBILITY_PROMISCUOUS_OR_HOUSEKEEPING` (ubiquitous nuclear transcripts like MALAT1, NEAT1, SNHG*), 195 as `EXCLUDED_DONOR_SEX_BIAS` (XIST), and 143 as `PREDICTED_ONLY` (no co-expression or miRNA detection support).
- **lncRNA Credibility Groups:** 57 unique lncRNAs (2,286 axes) are classified as `ROBUST_CANDIDATE`, 13 unique lncRNAs (916 axes) as `LOW_CREDIBILITY_CANDIDATE`, and 1 unique lncRNA (195 axes) as `EXCLUDED_FROM_INTERPRETATION_DONOR_SEX` (XIST).
- **lncRNA Reference Reconciliation:** Reference GENCODE v32 contains **17,944 annotated lncRNA genes**, of which **16,762** are assayed in the 10x Genomics Cell Ranger GRCh38 feature catalog for GSE248762 (1,182-feature difference not investigated). In single-cell stromal profiling, **867 lncRNAs** pass the primary detection filter (>= 1.0% stromal cells and >= 2 donors: 834 robust candidates, 32 low-credibility transcripts like MALAT1/NEAT1/SNHG*, and 1 XIST excluded for donor sex bias).

## 11. External Reference DOI Resolution Audit
**Source File:** `scripts/audit_checks/verify_dois.py`
**Command:** `python scripts/audit_checks/verify_dois.py`

| Citation Key | Registered DOI | First Author | Container / Journal | Publication Year | Official Title | Resolution Status |
| :--- | :--- | :--- | :--- | :---: | :--- | :---: |
| **Charoentong 2017** | `10.1016/j.celrep.2016.12.019` | Charoentong, Pornpimol | Cell Reports | 2017 | Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade | **SUCCESS (HTTP 200)** |
| **Bindea 2013** | `10.1016/j.immuni.2013.10.003` | Bindea, Gabriela | Immunity | 2013 | Spatiotemporal Dynamics of Intratumoral Immune Cells Reveal the Immune Landscape in Human Cancer | **SUCCESS (HTTP 200)** |
| **Naba 2012** | `10.1074/mcp.M111.014647` | Naba, Alexandra | Molecular & Cellular Proteomics | 2012 | The Matrisome: In Silico Definition and In Vivo Characterization by Proteomics of Normal and Tumor Extracellular Matrices | **SUCCESS (HTTP 200)** |
| **Naba 2016** | `10.1016/j.matbio.2015.06.003` | Naba, Alexandra | Matrix Biology | 2016 | The extracellular matrix: Tools and insights for the “omics” era | **SUCCESS (HTTP 200)** |
| **Ma 2025** | `10.3389/fphar.2025.1646948` | Ma, Xin | Frontiers in Pharmacology | 2025 | Role of biomarker SOCS1 in peritoneal dialysis-associated peritoneal fibrosis and immune infiltration based on machine learning screening | **SUCCESS (HTTP 200)** |
| **Kuppe 2020** | `10.1038/s41586-020-2941-1` | Kuppe, Christoph | Nature | 2020 | Decoding myofibroblast origins in human kidney fibrosis | **SUCCESS (HTTP 200)** |

## 12. Figure File Routing Specification
| Figure / Context | Target Image File | Caption / Description Rule |
| :--- | :--- | :--- |
| **STRING PPI Network Only** | [`results/figures/D_string_ppi_network.png`](results/figures/D_string_ppi_network.png) | Dedicated STRING-only functional protein association network (21 edges) |
| **STRING + Tissue Co-Expression** | [`results/figures/Hub_01_ppi_gene_interaction_network.png`](results/figures/Hub_01_ppi_gene_interaction_network.png) | Dual network displaying live STRING edges and empirical tissue co-expression (|r| >= 0.85, 42 edges) |
| **QC Violins (Per-Sample)** | [`results/figures/E2_qc_violins_per_sample.png`](results/figures/E2_qc_violins_per_sample.png) | Legend: 'Deposited barcodes (>= ~500 UMI)' |
| **QC Violins (Before/After)** | [`results/figures/E2_qc_violins_before_after.png`](results/figures/E2_qc_violins_before_after.png) | Legend: 'Deposited barcodes (>= ~500 UMI)' |

## 13. Synthetic Identifier Excision Disclosure
During an internal quality control audit in September–October 2026, two downstream exploratory scripts (`11_drug_repurposing_dgidb.py` / `11_candidate_drug_repurposing.py` and `12_ihc_protein_validation.py`) were identified as containing AI-generated synthetic identifiers (fabricated antibody catalogue numbers and PMIDs) that do not correspond to authentic database entries or physical reagents.
In accordance with strict scientific integrity standards:
1. Both downstream exploratory scripts—`11_candidate_drug_repurposing.py` (which replaced earlier `11_drug_repurposing_dgidb.py`) and `12_ihc_protein_validation.py`—were excised from the active working tree via `git rm`. However, both files remain present in the public commit history on `origin/main` (and in the local backup branch `backup-main-before-scrub`). A permanent purge from repository history would require a git history rewrite (e.g. `git-filter-repo` or BFG) and forced pushing remote branches and tags (`manuscript-v1`), which has not been performed; history is retained and we rely on this written disclosure.
2. All active code and documentation reflect only verified, reproducible transcriptomics, STRING v12.5 PPI networks, and permutation tests. All outputs associated with the excised scripts (`candidate_drugs_*.csv`, `hub_genes_ihc_*.csv`, `Hub_05_*`, `Hub_06_*`) are superseded and excluded.