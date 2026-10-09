# GSE62928 Eight-Gene Validation Preregistration

## Registration status

This analysis plan was written before opening or analyzing GSE62928 expression
values. The cohort mapping will first be checked from GEO sample metadata. The
analysis will stop for clarification if the deposited groups do not match four
severe EPS biopsies, two early-PD controls, and two uremic controls.

## Fixed gene panel and expected direction

The eight discovery-selected genes are fixed and will not be changed after
inspection of validation expression:

`TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, SERPINA10, CST7`

The expected direction for every gene is higher expression in severe EPS
biopsies than in the combined non-EPS controls. This follows the positive
GSE125498 LPD-versus-SPD discovery direction. `TNFSF15` and `FLT3LG` are the
pre-specified strict 3/3 ML genes and will also be analyzed as a secondary
two-gene set.

Discovery log2 fold changes and nominal P-values will be read from
`results/tables/GSE125498_all_results.csv` at execution time.

## Data and sample handling

- Use the real GSE62928 GPL13158 series matrix and GPL13158 probe annotation.
- Confirm each sample ID and group from deposited GEO metadata before reading
  the expression table. No samples will be removed.
- Use deposited matrix values if the series-matrix processing description and
  value distribution support already log2/RMA-normalized expression. Otherwise,
  normalize the raw expression data using RMA and document the implementation.
- Report a per-sample expression boxplot and a PCA plot colored by the
  metadata-derived group.
- EPS versus non-EPS is the primary contrast. Early-PD and uremic controls will
  remain combined for the primary analysis and shown separately only in
  descriptive plots.

## Outcome-blind probe selection and detection

- Map probes to symbols using the GPL13158 annotation, retaining the original
  probe IDs and all mappings for each panel gene.
- Primary probe: among probes mapping to each gene, select the one with the
  highest mean expression across all eight samples. This selection does not
  use group labels.
- Sensitivity 1: select the probe with the highest across-sample IQR.
- Sensitivity 2: report each mapped probe separately.
- For an operational detection flag, the selected probe is detected if its
  mean log2 expression over the eight samples is greater than 5.0. This is an
  explicit reporting threshold, not a manufacturer-provided negative-control
  background estimate; do not interpret it as a biological detection limit.
- Stop and report the missing gene/probe if any fixed gene has no GPL13158
  probe mapping.

## Pre-specified statistical tests

### Per gene, primary probe

- Estimate EPS-minus-non-EPS log2 fold change, 95% confidence interval, and
  moderated t-test P-value using limma.
- Also report a two-sided Mann-Whitney U P-value.
- Calculate pooled-standard-deviation Cohen's d (EPS minus non-EPS) and a
  percentile bootstrap 95% CI using 10,000 stratified resamples (resample four
  EPS and four non-EPS observations independently with replacement; seed 42).
- Call direction concordant when the validation EPS-minus-control estimate is
  positive, matching the positive discovery direction.
- Apply Benjamini-Hochberg adjustment across the eight primary-probe limma
  P-values.

### Gene-set level

- Exact one-sided binomial sign test for the number of direction-concordant
  genes among eight, under null concordance probability 0.5.
- Composite score: z-score each primary-probe gene across all eight samples,
  using the sample standard deviation, then take the unweighted mean of the
  eight z-scores per sample. Report EPS-minus-control mean score difference,
  two-sided Mann-Whitney U P-value, and AUC with a stratified percentile
  bootstrap 95% CI (10,000 resamples, four samples per group, seed 42). AUC is
  oriented so higher scores predict EPS.
- Random-gene-set null: 1,000 sets of eight distinct detected genes, matched
  to the hub genes' decile bins of mean expression over all eight samples.
  Draw one gene for each hub from its same expression decile, without
  replacement within a set. Report null distributions for composite AUC and
  mean gene-level EPS-minus-non-EPS log2 fold change. Use the upper-tail
  empirical P-value `(1 + count(null >= observed)) / 1001`, with seed 42.
  If exact-decile candidates are insufficient to generate 1,000 unique
  matched sets, stop and report the limitation rather than silently changing
  the matching rule.
- Repeat the sign test, composite score, bootstrap, and matched null analysis
  for the two strict genes `TNFSF15` and `FLT3LG`, as secondary results.

### Marker-score proxy, descriptive subgroup plots, and power

- Use the pan-leukocyte marker genes in `scripts/task4_marker_sets.py`, removing
  all eight fixed hub genes before scoring. Average per-gene z-scores across
  the eight validation samples to obtain a marker-score proxy, not cell-fraction
  deconvolution. Report its EPS-versus-combined-control difference and
  two-sided Mann-Whitney P-value, and its correlation with the eight-gene
  composite score.
- Show descriptive expression plots for severe EPS, early-PD, and uremic
  groups; do not conduct subgroup tests because each control subgroup has only
  two samples.
- Calculate the minimum detectable standardized effect for two independent
  groups of four at 80% power and two-sided alpha 0.05 using the two-sample
  t-test power model; state the model and value in the report.

## Fixed interpretation criteria

- **Supportive:** at least 6 of 8 genes concordant in direction **and**
  composite-score AUC >= 0.80 **and** empirical null P < 0.05.
- **Partially supportive:** 5 of 8 concordant **or** composite AUC 0.65-0.80
  without significance.
- **Not supported:** otherwise.

Apply these criteria exactly; do not alter them after examining results.
Report all limitations, including n = 4 versus 4, discovery effluent versus
validation biopsy tissue, phenotype mismatch (LPD versus severe EPS), absence
of covariate adjustment, and the possibility that wide intervals are
inconclusive rather than evidence against an effect.

## Reproducibility and outputs

Use seed 42 for all bootstrap and random gene-set sampling. Use only measured
GSE62928 expression and deposited annotation; do not simulate expression or
hard-code numerical results. Save the analysis script, source tables, figures,
and report under this directory. The README update will use only values from
the resulting saved outputs.
