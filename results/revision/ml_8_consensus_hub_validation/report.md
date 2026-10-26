# External validation of the revised eight consensus hub genes

**Status:** exploratory external validation  
**Panel:** TNFSF15, EBI3, PI3, FLT3LG, ADAM19, CLEC4F, TGM3, CXCL14  
**Datasets:** GSE62928, GSE248762, and GSE130888

## Executive summary

The 31-gene upregulated ECM candidate table was evaluated by three machine
learning models—LASSO, SVM-RFE, and Random Forest. A gene was retained as a
consensus hub if at least two models voted for it. This yielded **8 consensus
hubs**, of which **2 received all three votes**. The 8-gene panel is listed in
[the panel CSV](../../final_analysis/ML_3_models_8_consensus_hub_genes.csv).
The historical 11-gene result came from an older four-model analysis that
included XGBoost; it is not an intermediate step in this three-model
31-to-8 selection.

The panel was assessed in a peritoneal-tissue dataset (GSE62928) and two
peritoneal-dialysis effluent single-cell datasets (GSE248762, GSE130888).
Results were mixed:

- **GSE62928:** 2/8 validation effects were in the discovery direction.
  The composite score AUC was 0.375 (95% bootstrap interval 0.000–0.750);
  the matched expression-null AUC P-value was 0.887. This did not provide
  convincing panel-level validation.
- **GSE248762:** 7/8 genes met the raw-count detection rule and 6/7 detected
  genes were direction-concordant. The eighth gene was not detected, so the
  fixed eight-gene score and AUC were not calculated.
- **GSE130888:** 5/8 genes were detected and 4/5 detected genes were
  direction-concordant. Three were undetected, so the fixed panel score and
  AUC were not calculated.

**Overall interpretation:** GSE248762 is the stronger *exploratory* effluent
check by sample count, gene detection, and direction concordance, but it is
not a complete validation. No dataset tested the full panel successfully
across all eight genes. These analyses do not justify claiming that all eight
hub genes have been validated.

## 1. Panel selection

The tested panel comes from
[`ML_3_models_hub_genes_all31_up_ecm.csv`](../../final_analysis/ML_3_models_hub_genes_all31_up_ecm.csv),
which contains 31 candidate rows. The three model vote columns are LASSO,
SVM-RFE, and Random Forest. At least two positive votes defined a consensus
hub. The resulting eight genes and their votes are:

| Gene | LASSO | SVM-RFE | Random Forest | Total |
|---|---:|---:|---:|---:|
| TNFSF15 | 1 | 1 | 1 | 3 |
| EBI3 | 1 | 1 | 1 | 3 |
| PI3 | 1 | 1 | 0 | 2 |
| FLT3LG | 1 | 0 | 1 | 2 |
| ADAM19 | 1 | 0 | 1 | 2 |
| CLEC4F | 1 | 1 | 0 | 2 |
| TGM3 | 1 | 1 | 0 | 2 |
| CXCL14 | 1 | 1 | 0 | 2 |

Thus the panel comprises two strict 3/3 genes and six 2/3 genes. The separate
`ML_3_models_hub_genes.csv` is a different candidate pool and should not be
used as the source of this particular eight-gene list.

## 2. Analysis methods

### GSE62928: peritoneal tissue

The GEO GPL13158 Affymetrix HG-U133 Plus 2.0 deposited RMA log2 matrix was
used without additional normalization. The comparison was four EPS biopsies
versus four non-EPS samples (two PD and two uremic controls). For each gene,
the probe with the highest mean expression across all eight samples was
selected without using group labels. limma estimated the EPS-minus-non-EPS
log2 fold change and moderated P-value; Benjamini-Hochberg correction was
applied across the eight genes.

The panel score is the unweighted average of per-gene z-scores across the
eight samples. Its AUC was compared with 1,000 random gene sets matched to
expression deciles. The score AUC interval used 10,000 stratified bootstrap
resamples. The expression cutoff (mean RMA log2 expression >5) is an
operational threshold, not a manufacturer-defined detection limit.

### GSE248762 and GSE130888: effluent single-cell datasets

Previously generated, QC-filtered patient-level raw-count pseudobulks were
used; cells were not treated as independent replicates. edgeR quasi-likelihood
models estimated log2 fold changes, P-values, and approximate 95% confidence
intervals. The contrasts were:

- GSE248762: LV_UF plus LV_NOT_UF (10 patient pseudobulks) versus SV (6).
- GSE130888: long-term PD (4) versus short-term PD (6); three normal controls
  were excluded from the inferential contrast.

An effluent gene was called detected only if raw pseudobulk counts were at
least 10 in at least half of patient samples **in each** contrast group.
Undetected genes were not counted as either concordant or discordant and
were not assigned effect estimates or P-values. Direction concordance means
the validation log2FC has the same sign as the discovery log2FC. BH-adjusted
P-values in the panel tables are calculated across detected panel genes for
that cohort.

For the revised-panel edgeR re-fit, the feature universe included genes
selected by `filterByExpr`, genes meeting the raw-count detection rule in both
groups, and all eight panel genes (forced in so low-count panel genes could
be explicitly assessed for detection). This preserves the existing cohort
pipeline's filtering approach while centering the test on the revised panel.

This raw-count detection criterion is not the same as the operational
mean-expression cutoff used for the microarray. Results across platforms
should not be interpreted as directly comparable effect-size measurements.

## 3. Results

### Dataset-level comparison

| Dataset | Contrast | Sample sizes | Detected | Concordant among detected | Full 8-gene composite |
|---|---|---:|---:|---:|---|
| GSE62928 | EPS vs non-EPS | 4 vs 4 | 1/8 above operational cutoff | 2/8 effects positive | Computed; AUC 0.375 |
| GSE248762 | LV vs SV | 10 vs 6 | 7/8 | 6/7 (85.7%; one-sided sign P=0.0625) | Not estimable |
| GSE130888 | Long-term vs short-term PD | 4 vs 6 | 5/8 | 4/5 (80.0%; one-sided sign P=0.1875) | Not estimable |

The two “concordance” summaries for GSE62928 and effluent cohorts have
different denominators: the GSE62928 value counts positive validation
directions among all eight genes, while effluent concordance excludes genes
that failed the raw-count detection rule.

### Per-gene effects

Effect is validation log2FC in the stated direction (EPS minus non-EPS for
GSE62928; LV minus SV for the effluent cohorts). `q` is the BH-adjusted
P-value. `ND` means the dataset-specific detection rule defined for this
analysis failed,
so no validation estimate/P-value is reported.

| Gene | GSE62928 log2FC [95% CI]; P; q | GSE248762 log2FC [95% CI]; P; q | GSE130888 log2FC [95% CI]; P; q |
|---|---|---|---|
| TNFSF15 | -0.052 [-0.517, 0.413]; 0.8066; 0.9540 | 1.324 [-0.623, 3.271]; 0.1069; 0.1496 | ND |
| EBI3 | -0.025 [-0.569, 0.519]; 0.9191; 0.9540 | 1.371 [0.236, 2.506]; 0.0075; 0.0261 | 0.838 [-1.064, 2.741]; 0.1915; 0.2613 |
| PI3 | 0.024 [-0.881, 0.929]; 0.9540; 0.9540 | 8.749 [4.388, 13.111]; 0.0024; 0.0168 | 1.355 [-1.709, 4.418]; 0.2090; 0.2613 |
| FLT3LG | -0.266 [-0.949, 0.417]; 0.4015; 0.9540 | 1.105 [-0.477, 2.687]; 0.0925; 0.1496 | 0.825 [-1.682, 3.333]; 0.3052; 0.3052 |
| ADAM19 | 1.048 [0.246, 1.851]; 0.0160; 0.1281 | 1.068 [-0.058, 2.194]; 0.0268; 0.0626 | 0.923 [-1.208, 3.054]; 0.1844; 0.2613 |
| CLEC4F | -0.082 [-0.602, 0.437]; 0.7287; 0.9540 | -0.109 [-1.833, 1.615]; 0.8685; 0.8685 | -1.168 [-3.399, 1.062]; 0.1465; 0.2613 |
| TGM3 | -0.133 [-0.583, 0.318]; 0.5229; 0.9540 | ND | ND |
| CXCL14 | -1.714 [-3.819, 0.392]; 0.0989; 0.3954 | 0.274 [-1.397, 1.945]; 0.6811; 0.7946 | ND |

In GSE248762, EBI3 and PI3 had BH-adjusted q<0.05 among the detected panel
genes. PI3's estimated effect is unusually large (log2FC 8.75) and should be
treated cautiously; the approximate 95% interval is wide (4.39 to 13.11).
No GSE130888 detected gene had nominal P<0.05. In GSE62928, ADAM19 had a
nominal P=0.016 but was not significant after correction across eight genes
(q=0.128).

### GSE62928 panel-level statistics

| Statistic | Result |
|---|---:|
| Direction-concordant genes | 2/8 (PI3, ADAM19) |
| Composite mean score difference (EPS minus non-EPS) | -0.1703 |
| Exact two-sided Mann-Whitney P for composite | 0.6857 |
| Composite AUC | 0.375 |
| Stratified-bootstrap 95% AUC interval | 0.000–0.750 |
| Mean per-gene validation log2FC | -0.1500 |
| Matched-null mean AUC | 0.5381 |
| Matched-null AUC empirical P | 0.8871 |
| Genes above operational mean-RMA cutoff | 1/8 |

The revised panel's AUC was not unusual relative to the expression-matched
null. Given the very small cohort and low measured expression for seven genes,
this result does not establish discrimination or reliable per-gene effects.

### Detection counts in the effluent cohorts

These columns count patient pseudobulks with raw counts >=10. A gene passed
only if at least half of samples in both groups met that threshold.

| Gene | GSE248762 LV / SV | GSE130888 LV / SV |
|---|---:|---:|
| TNFSF15 | 9/10, 5/6 | 1/4, 6/6 |
| EBI3 | 10/10, 6/6 | 4/4, 6/6 |
| PI3 | 10/10, 4/6 | 2/4, 5/6 |
| FLT3LG | 10/10, 6/6 | 4/4, 6/6 |
| ADAM19 | 10/10, 6/6 | 4/4, 6/6 |
| CLEC4F | 10/10, 6/6 | 4/4, 6/6 |
| TGM3 | 1/10, 0/6 | 0/4, 0/6 |
| CXCL14 | 5/10, 3/6 | 1/4, 2/6 |

## 4. Interpretation and limitations

1. **No full effluent panel score:** TGM3 was undetected in both effluent
   datasets, and additional genes were undetected in GSE130888. The full
   fixed eight-gene score, AUC, and matched-null composite test therefore
   cannot be calculated there. Scoring only detected genes would change the
   panel and was not substituted.
2. **Small cohorts:** GSE62928 has four samples per arm, GSE248762 has 10
   versus 6 patient pseudobulks, and GSE130888 has 4 versus 6. Confidence
   intervals are consequently broad and power is limited.
3. **Different sample types:** GSE62928 is peritoneal tissue; GSE248762 and
   GSE130888 are effluent cells. The effluent data do not directly reproduce
   the tissue biopsy context.
4. **Composition and cell-type effects:** The RNA signal can reflect cell
   composition as well as expression within a cell type. The pseudobulk
   comparisons do not establish cell-intrinsic expression or causality.
5. **PI3 effect caution:** The large GSE248762 PI3 estimate and interval
   warrant checking raw counts, patient influence, and composition before
   drawing a biological conclusion.
6. **Exploratory panel reuse:** This exact revised 8-gene list was not the
   preregistered panel for the existing GSE62928 and effluent analyses.
   These results are exploratory, not a preregistered replication.
7. **No definitive “best validation”:** GSE248762 has more LV samples and
   higher observed detection/concordance than GSE130888, but neither permits
   full-panel scoring. It is the stronger exploratory follow-up by coverage,
   not a confirmed validation cohort.

## 5. Reproducibility and output map

Start with the [study README](README.md) for exact PowerShell commands and
input requirements.

- [GSE62928 per-gene results](GSE62928/per_gene_validation.csv)
- [GSE62928 panel summary](GSE62928/panel_validation_summary.csv)
- [GSE62928 probe choices](GSE62928/probe_selection.csv)
- [GSE62928 sample scores](GSE62928/sample_composite_scores.csv)
- [GSE62928 matched random-gene null](GSE62928/matched_gene_set_null.csv)
- [GSE248762 per-gene results](GSE248762/per_gene_validation.csv)
- [GSE248762 panel summary](GSE248762/panel_validation_summary.csv)
- [GSE248762 sample mapping](GSE248762/sample_group_mapping.csv)
- [GSE130888 per-gene results](GSE130888/per_gene_validation.csv)
- [GSE130888 panel summary](GSE130888/panel_validation_summary.csv)
- [GSE130888 sample mapping](GSE130888/sample_group_mapping.csv)
- [Effluent comparison memo](effluent_cohort_comparison.md)
- [Validation scripts and guide](../../../scripts/validation/README.md)

The tested R environment was R 4.4.2, limma 3.62.2, and edgeR 4.4.2.
Each validation run writes `sessionInfo.txt` beside its result tables for
environment provenance.
