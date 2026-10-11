# Changelog: Methodological Audit and Scientific Corrections

All updated estimates below were regenerated with seed 42 from the checked-in discovery tables and GEO series matrices. The scripts and machine-readable outputs are under `scripts/task*.py` and `results/revision/`. The current files also contain unresolved source inconsistencies; claims affected by those inconsistencies have been withdrawn rather than inferred.

## Quantitative changes

| Analysis | Previously reported | Audited / revised | Rationale |
| :--- | :--- | :--- | :--- |
| GSE62928 validation AUC | Apparent AUC = 1.0000; LOOCV AUC = 1.0000 | Fixed seven-hub score AUC = 0.4375 | Legacy pipeline fit on a synthetic Gaussian matrix with a +1.5 case shift; this is not validation on the biopsy matrix. |
| GSE62928 exact permutation P | P = 1.0000, presented as evidence alongside perfect AUC | Withdrawn as uninterpretable; rank-based tests on real expression: one-sided Mann-Whitney P = 0.6571 and two-sided P = 0.8857 | The historical code fit an unpenalized 11-feature model in-sample on eight samples. All 70 four-versus-four label combinations had apparent AUC 1.0000. |
| Validation direction concordance | Not reported | 2/7 genes (28.6%); two-sided binomial sign-test P = 0.4531 | Recomputed per-gene signs in GSE62928. |
| Random 7-gene-set validation null | Not reported; supplied draft said mean AUC = 0.4804 | Mean AUC = 0.4798; empirical P = 0.6040 (604/1,000 null AUCs >= hub AUC) | Re-run of the 1,000-set null with seed 42. The supplied 0.4804 differs from the reproduced 0.4798. |
| Discovery apparent AUC | 1.0000 | Nested repeated 5-fold mean AUC = 0.7700; percentile-bootstrap 95% interval [0.7408, 0.7946], SD = 0.0453 | DEG, Matrisome, LASSO, SVM-RFE, and RF selection moved inside folds. WGCNA was not recomputed inside folds. The interval bootstraps the 10 overlapping repeat-level estimates. |
| Nomogram / selected-model C-index | Apparent = 0.9808; preselected-feature CV = 0.8731 | Nested logistic-score C-index = 0.7700; percentile-bootstrap 95% interval [0.7408, 0.7946] | Binary classification C-index equals ROC AUC here. This is the nested fold-trained logistic score, not a separately regenerated graphical nomogram; the nested pipeline also does not recompute WGCNA. |
| ML algorithms and hub count | Four algorithms including XGBoost; 11 hubs in the legacy summary | Three-model file (LASSO, SVM-RFE, RF): 2 strict 3/3 votes and 7 majority 2/3 votes; separate 31-Up-ECM analysis: 2 strict and 8 majority hubs | The current three-model output is retained as the recorded vote result. Its WGCNA provenance is not verified (see below). |
| FDR sensitivity at 0.05 | Not reported in the primary summary | 135/813 DEGs (16.6%); 6/44 matrisomal DEGs (13.6%); 4/7 recorded hubs (57.1%) | Benjamini-Hochberg results are supplementary; the primary DEG rule remains absolute log2FC >= 0.585 and nominal P < 0.05 without FDR filtering. |
| FDR sensitivity at 0.10 | Not reported | 263/813 DEGs (32.3%); 10/44 matrisomal DEGs (22.7%); 5/7 recorded hubs (71.4%) | Supplementary sensitivity tier only. |
| Nominal false-positive expectation | Not stated prominently | About 587 expected at alpha = 0.05 over 11,741 tested genes, under a global-null calculation | This is an expectation under the null, not a count of actual false discoveries. |
| Leukocyte-composition sensitivity | Not evaluated | Marker-based score higher in LPD (P = 0.0010); 4/7 stage associations remain nominal after score adjustment; EBI3 is borderline | The computed score is a mean of per-gene z-scores, not a validated deconvolution method. Associations do not establish cell-intrinsic expression or causality. |
| Matrisome enrichment | Described as significant without emphasizing magnitude | 1.7145-fold enrichment | Reframed as modest enrichment. |

## Real GSE62928 per-gene effects

Effects are validation minus control log2 expression. Every 95% confidence interval includes zero; none of these per-gene comparisons is nominally significant.

| Hub | Validation log2FC | 95% CI | Welch P | Direction |
| :--- | ---: | :--- | ---: | :--- |
| ADAM19 | +0.4562 | [-0.1883, +1.1008] | 0.1308 | Concordant |
| LTB | +0.1493 | [-0.4418, +0.7404] | 0.5551 | Concordant |
| EBI3 | -0.0252 | [-0.5964, +0.5461] | 0.9171 | Discordant |
| SERPINA10 | -0.0455 | [-0.8770, +0.7859] | 0.8974 | Discordant |
| TNFSF15 | -0.0519 | [-0.4747, +0.3710] | 0.7739 | Discordant |
| FLT3LG | -0.3493 | [-0.9510, +0.2525] | 0.1791 | Discordant |
| CXCL14 | -0.9292 | [-2.8847, +1.0264] | 0.2803 | Discordant |

The composite-score AUC and gene-set null are summarized in the quantitative changes table; the complete source values are in [`results/revision/gse62928_per_gene_concordance.csv`](./results/revision/gse62928_per_gene_concordance.csv) and [`results/revision/gse62928_robustness_metrics.json`](./results/revision/gse62928_robustness_metrics.json).

## WGCNA and hub-provenance corrections

The supplied draft described MEblack as the biological hub cluster because it contains 25 of 31 upregulated ECM DEGs. The audited module-trait results do not support that interpretation:

* MEmagenta: stage correlation r = +0.4804, P = 0.00466; 37 genes and zero upregulated ECM DEGs.
* MEblack: 2,398 genes and 25/31 upregulated ECM DEGs, but stage correlation r = +0.0285, P = 0.8748.
* All seven majority hubs are assigned to black, while their calculated kME-to-black correlations range from -0.2996 to +0.0393 and are nonsignificant. Their stage gene-significance correlations range from +0.4616 to +0.6502.
* The 25-row three-model feature table records the seven majority hubs, but the separate WGCNA convergence table gives every one of them `in_Top_WGCNA_Module=0` and `Convergent=0`. The assertion that these seven hubs derive from a 25-gene WGCNA-convergent pool is therefore not verified by the available files. The README reports this conflict instead of choosing a convenient source.
* The alternate 31-gene Up-ECM analysis remains separately identified and has 8 majority hubs: TNFSF15, EBI3, FLT3LG, ADAM19, CLEC4F, CXCL14, PI3, and TGM3.

Per-hub kME and gene-significance values are in [`results/revision/task1_audit_report.md`](./results/revision/task1_audit_report.md) and `.json`.

## Leukocyte score and adjusted hub associations

The signature score is based on mean per-gene z-scores from canonical leukocyte markers; it is not xCell, CIBERSORT, MCP-counter, or a validated cell-fraction estimate. No Rscript executable or xCell/CIBERSORT/decoupler implementation was available in the checked environment, so the requested true deconvolution was not completed. Mean total-leukocyte marker score was 0.3011 in LPD and -0.1957 in SPD (P = 0.0010); mean lymphocyte score was 0.4154 and -0.2700, respectively (P = 0.0005).

| Hub | Pearson r with total leukocyte score | P-value | Pearson r with lymphocyte score |
| :--- | ---: | ---: | ---: |
| SERPINA10 | +0.7690 | <0.0001 | +0.7972 |
| FLT3LG | +0.7498 | <0.0001 | +0.7978 |
| ADAM19 | +0.6770 | <0.0001 | +0.6407 |
| LTB | +0.6379 | <0.0001 | +0.6201 |
| TNFSF15 | +0.5514 | 0.0009 | +0.5314 |
| EBI3 | +0.4271 | 0.0132 | +0.4559 |
| CXCL14 | +0.2034 | 0.2564 | +0.1561 |

OLS stage associations after adjustment for the total-leukocyte marker score were:

| Hub | Adjusted stage beta | P-value | Interpretation at nominal 0.05 |
| :--- | ---: | ---: | :--- |
| TNFSF15 | +0.7454 | 0.0036 | Retained |
| CXCL14 | +0.6929 | 0.0130 | Retained |
| FLT3LG | +0.7425 | 0.0167 | Retained |
| ADAM19 | +0.6540 | 0.0400 | Retained |
| EBI3 | +0.6049 | 0.0532 | Borderline |
| SERPINA10 | +0.3701 | 0.2162 | Not significant |
| LTB | +0.3285 | 0.3371 | Not significant |

The tests are small-sample sensitivity models, not proof that the retained associations are independent of cell composition.

## Corrections to copied DEG tables and gene claims

The earlier README's 44-gene transcription did not match the current `FINAL_DEG_Matrisomal_Genes.csv`: 20 upregulated rows had discrepant values, and all nine downregulated entries were absent from the current 44-gene source table (some also failed the current primary DEG criteria). The table was replaced with links to the checked-in result CSVs; a row-by-row comparison is preserved in [`results/revision/readme_deg_table_reconciliation.csv`](./results/revision/readme_deg_table_reconciliation.csv).

Under the requested borderline screen (absolute log2FC < 0.62 or nominal P > 0.04), the current source flags LAMB1 (log2FC = -0.6054, P = 0.0431), CCBE1 (+0.6050, P = 0.00594), and COL8A1 (+0.7493, P = 0.0488).

The requested statement that TGFB2 and COL6A2 are downregulated cannot be retained: TGFB2 is absent from current gene-result tables; COL6A2 is present with log2FC = +0.6468 and P = 0.192, making it positive and nonsignificant in the current output. This discrepancy requires corrected upstream source data before it can be resolved.

## Other verified method and metadata details

* GSE125498's GEO processing record states Illumina BeadStudio v3 summarization, quantile normalization in R via `lumi`, and filtering for detectable signal in at least five samples.
* The public GSE125498 matrix includes patient ID, PD duration group, and cell type but not patient age, sex, race, or microarray batch dates. Covariate-adjusted limma analyses and DEG-overlap estimates for those unavailable fields were not possible.
* No experimental validation was performed.
* The repository's existing MIT license was retained. Direct Python dependencies are pinned in `requirements.txt`.

## Revision outputs

* `results/revision/task1_audit_report.md` / `.json`
* `results/revision/nested_cv_results.csv`, `nested_cv_summary.json`, and `nested_cv_fold_log.csv`
* `results/revision/gse62928_per_gene_concordance.csv`, `gse62928_robustness_metrics.json`, and validation figures
* `results/revision/leukocyte_scores_per_sample.csv`, `leukocyte_confounding_hub_analysis.csv`, `task4_leukocyte_summary.json`, and confounding figures
* `results/revision/fdr_sensitivity_table.csv`, `hubs_fdr_breakdown.csv`, and `task5_fdr_summary.json`
* `results/revision/readme_deg_table_reconciliation.csv`

The legacy full pipeline still contains an XGBoost-based ML phase and synthetic Phase 6 validation code. The revision scripts and README explicitly avoid treating those outputs as the audited three-model and real-biopsy results; raw-to-report regeneration remains incomplete until this source path and the conflicting WGCNA provenance are reconciled.

## 2026-10-11 — Panel detectability audit and exploratory panel v2

Added `scripts/task6_serpina10_detectability_audit.py` and the new
`results/revision/panel_v2/` outputs to document detectability of the fixed
eight-gene panel across GSE62928, GSE248762, and GSE130888; preserve the
discovery evidence and rule-by-rule panel decision; and quantify clearly
post-hoc, exploratory panel-v2 scores on existing cohorts. The rule-derived
five-gene panel is FLT3LG, EBI3, LTB, ADAM19, and CST7. This addition addresses
the assay-detectability limitation, including why SERPINA10 was not measurable
under cohort-specific filters; it does not change or replace any result.

All preregistered outputs, thresholds, and verdicts remain unchanged. Existing
cohort score calculations for panel v2 are explicitly labeled post hoc,
exploratory, and not validation claims.
