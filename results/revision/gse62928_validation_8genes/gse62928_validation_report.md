# GSE62928 Preregistered Eight-Gene Validation Report

## Registration and data provenance

The fixed panel, expected direction, tests, probe rules, and decision criteria were preregistered in [PREREGISTRATION.md](PREREGISTRATION.md) before GSE62928 expression values were read. The validation matrix was downloaded from GEO to [GSE62928_series_matrix.txt.gz](GSE62928_series_matrix.txt.gz). It identifies platform GPL13158, Affymetrix HG-U133 Plus 2.0. GEO states that values were processed with Robust Multi-array Average (RMA), including background adjustment, quantile normalization, and median-polish summarization. The deposited RMA matrix is treated as log2-scale; no further normalization was applied.

### GEO sample-to-group mapping

| Sample ID | GEO title | Analysis group | Group-label basis |
|---|---|---|---|
| GSM1536406 | EPS1 | EPS | GEO identifies EPS; severe-EPS status follows supplied cohort description |
| GSM1536407 | EPS2 | EPS | GEO identifies EPS; severe-EPS status follows supplied cohort description |
| GSM1536408 | EPS3 | EPS | GEO identifies EPS; severe-EPS status follows supplied cohort description |
| GSM1536409 | EPS4 | EPS | GEO identifies EPS; severe-EPS status follows supplied cohort description |
| GSM1536410 | PD1 | Early PD | GEO title/treatment says PD; early-PD status follows supplied cohort description |
| GSM1536411 | PD2 | Early PD | GEO title/treatment says PD; early-PD status follows supplied cohort description |
| GSM1536412 | UREMIC1 | Uremic | GEO title/treatment metadata |
| GSM1536413 | UREMIC2 | Uremic | GEO title/treatment metadata |

GEO titles and treatment characteristics identify four EPS biopsies, two PD samples, and two uremic controls. The analysis treats the EPS and PD samples as severe EPS and early PD, respectively, following the supplied cohort description; GEO's sample-level labels themselves do not specify severity or the early-PD stage. The primary comparison combines the two control subgroups (4 EPS versus 4 non-EPS). No samples were excluded.

## Probe selection and QC

For each gene, the primary probe was selected by highest mean expression across all samples, without using group labels. The IQR-selected probe and all mapped probes are recorded in [probe_selection_and_detection.csv](probe_selection_and_detection.csv), [all_panel_probe_inventory.csv](all_panel_probe_inventory.csv), and [sensitivity_all_probes_results.csv](sensitivity_all_probes_results.csv). All 8 genes mapped to GPL13158. The operational detection threshold was mean log2 expression > 5.0; 0 of 8 panel genes exceeded it. The genes below this reporting threshold were TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, SERPINA10, CST7. This threshold is an operational cutoff, not a manufacturer negative-control background estimate or a biological detection limit.

Per-sample distribution and PCA plots are [qc_expression_boxplot.pdf](qc_expression_boxplot.pdf) and [qc_pca_by_subgroup.pdf](qc_pca_by_subgroup.pdf); PCA uses the top 1,000 variable probes without labels. None flagged by the preregistered descriptive 2-SD PC1/PC2 screen. This descriptive screen did not remove samples. PCA coordinates and flags are in [qc_pca_scores.csv](qc_pca_scores.csv).

## Per-gene primary-probe results

Discovery log2FC and nominal P are taken directly from `results/tables/GSE125498_all_results.csv`. Validation log2FC, confidence intervals, and moderated P-values are from limma. Mann-Whitney P-values are exact two-sided permutation/rank probabilities for 4 versus 4. Cohen's d is pooled-SD standardized EPS-minus-control effect; its CI is a percentile bootstrap from 10,000 stratified resamples with seed 42. BH adjustment is across the eight primary-probe limma tests.

| Gene | Primary probe | Discovery log2FC (P) | Validation log2FC (95% limma CI) | limma P | Exact Mann-Whitney P | Cohen's d (bootstrap 95% CI) | BH P (8 genes) | Concordant |
|---|---|---:|---:|---:|---:|---:|---:|---|
| TNFSF15 | 221085_PM_at | 0.9743 (2.73e-05) | -0.0519 [-0.5165, 0.4128] | 0.8066 | 0.8857 | -0.2125 [-2.0569, 1.5228] | 0.9191 | No |
| FLT3LG | 206980_PM_s_at | 1.4377 (2.35e-05) | -0.2663 [-0.9494, 0.4167] | 0.4015 | 0.4857 | -0.5816 [-4.3847, 1.1967] | 0.6423 | No |
| EBI3 | 219424_PM_at | 0.8034 (0.00244) | -0.0252 [-0.5690, 0.5186] | 0.9191 | 1 | -0.0768 [-2.7674, 1.4304] | 0.9191 | No |
| LTB | 207339_PM_s_at | 1.6747 (3.33e-06) | 0.5294 [-0.0125, 1.0712] | 0.05445 | 0.1143 | 1.6251 [0.6433, 5.6554] | 0.1452 | Yes |
| ADAM19 | 209765_PM_at | 1.2131 (0.000133) | 1.0484 [0.2457, 1.8511] | 0.01601 | 0.02857 | 1.8645 [1.3158, 4.9367] | 0.1281 | Yes |
| ECM1 | 209365_PM_s_at | 1.0341 (0.00701) | 1.1892 [-0.3772, 2.7556] | 0.1203 | 0.3429 | 1.0052 [-0.2069, 5.3100] | 0.2406 | Yes |
| SERPINA10 | 220626_PM_at | 1.1659 (0.000653) | -0.0455 [-0.7513, 0.6602] | 0.8874 | 0.8857 | -0.0952 [-2.0704, 1.2033] | 0.9191 | No |
| CST7 | 210140_PM_at | 1.0123 (0.0005) | 1.4287 [0.0694, 2.7881] | 0.04136 | 0.1143 | 1.4030 [0.2475, 4.6575] | 0.1452 | Yes |

Direction concordance was 4/8. The one-sided exact sign-test P-value against a 0.5 concordance probability was 0.6367. The minimum attainable exact two-sided Mann-Whitney P for 4 versus 4 without ties is 2/70 = 0.0286.

The highest-IQR probe sensitivity results are in [sensitivity_highest_iqr_results.csv](sensitivity_highest_iqr_results.csv). Under that probe rule, 3/8 validation effects were positive. All-probe results are reported separately; genes with multiple probes can show probe-dependent effect directions, so the primary outcome-blind probe rule is retained for the prespecified verdict.

## Gene-set tests

| Set | Concordant | Score difference (EPS - non-EPS) | Exact Mann-Whitney P | AUC (bootstrap 95% CI) | Null mean AUC | Matched-null AUC empirical P | Mean gene log2FC | Matched-null log2FC empirical P |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Eight-gene panel | 4/8 | 0.4698 | 0.1143 | 0.875 [0.500, 1.000] | 0.5262 | 0.0150 | 0.4758 | 0.1079 |
| Strict genes: TNFSF15, FLT3LG | 0/2 | -0.4117 | 0.6857 | 0.375 [0.000, 0.812] | 0.5179 | 0.7982 | -0.1591 | 0.7373 |

For each fixed set, 1000 random sets were sampled from genes above the operational detection threshold, matched to the panel's mean-expression decile using the detected-gene universe. Null sampling was without replacement within each set, with seed 42; empirical P-values use the plus-one correction. All eight fixed panel genes fall below the operational detection threshold and below the detected-gene expression range; they are therefore assigned to the lowest decile for this coarse matching. This materially limits how closely the detected-gene null can match their absolute abundance and is reported as a limitation, not hidden.

The score is the unweighted mean of the per-gene z-scores, where each gene was standardized across all eight samples. ROC and matched-null distributions are shown in [composite_roc_and_null.pdf](composite_roc_and_null.pdf); sample-level scores are in [composite_scores_by_sample.csv](composite_scores_by_sample.csv).

## Marker-score proxy, subgroup description, and power

The leukocyte marker-score proxy used 17 mapped genes from the discovery Task 4 `LEUKOCYTE_MARKERS` list, excluding the fixed panel. It is a marker-expression score, not xCell/CIBERSORT or cell-fraction deconvolution. The EPS-minus-control score difference was 0.5193 (exact two-sided Mann-Whitney P = 0.0571); its Pearson correlation with the eight-gene composite was r = 0.6863 (P = 0.0602). Scores are in [leukocyte_marker_score_proxy.csv](leukocyte_marker_score_proxy.csv), and genes/probes used in [leukocyte_marker_genes_used.csv](leukocyte_marker_genes_used.csv).

Expression for EPS, early PD, and uremic subgroups is displayed descriptively (no subgroup hypothesis tests; each control subgroup has n = 2) in [gene_expression_subgroups_descriptive.pdf](gene_expression_subgroups_descriptive.pdf). The minimum detectable standardized effect at 80% power is Cohen's d = 2.3808 for a two-sample equal-variance independent t-test with n = 4 per group, two-sided alpha = 0.05, and 6 degrees of freedom.

## Preregistered verdict

**Not supported.** The fixed criteria require at least 6/8 direction-concordant genes, composite AUC >= 0.80, and matched-null AUC empirical P < 0.05 for “Supportive.” The observed panel meets the AUC and AUC-null conditions but has only 4/8 concordant genes, so it does not meet the conjunction. It also does not meet “Partially supportive” (which requires 5/8 concordant, or AUC from 0.65 to 0.80 without null significance). The secondary strict two-gene set is not supportive.

## Limitations

- The validation comparison has only 4 EPS and 4 controls. The power calculation shows that only very large standardized effects are detectable with the stated test; wide confidence intervals indicate uncertainty, not proof of no effect.
- Discovery used dialysate effluent from long- versus short-term PD; validation uses peritoneal biopsy tissue from severe EPS versus early-PD and uremic controls. Tissue, phenotype, and comparator differences limit direct transportability.
- No age, sex, batch, or other covariate adjustment was possible in this eight-sample analysis.
- The fixed panel genes are all below the prespecified operational mean-expression threshold of 5.0. The threshold is not a platform negative-control limit, but the low means and bottom-decile-only random-gene matching constrain interpretation.
- The continuous composite AUC remains calculable from the deposited expression values, but its 0.875 estimate with a 0.500–1.000 bootstrap interval is highly uncertain. Being below the operational cutoff alone neither establishes noise-dominated signal nor provides evidence of reliable discrimination.
- The eight-gene composite has a broad bootstrap interval, and probe sensitivities show that selecting a different annotated probe can change an individual gene's effect direction.
- This is an exploratory cohort evaluation, not evidence of a validated diagnostic test or causal mechanism.

## Reproducibility

Run `Rscript scripts/revision_gse62928_validation_8genes.R <series_matrix.gz> <GPL13158_annotation> results/revision/gse62928_validation_8genes` from the repository root, then run `python scripts/revision_gse62928_validation_8genes_report.py`. Analysis code: [revision_gse62928_validation_8genes.R](../../../scripts/revision_gse62928_validation_8genes.R). Software versions and machine-readable results are in [validation_summary.json](validation_summary.json).
