# Revised Eight-Consensus-Hub Validation in GSE62928

## Panel and analysis

This analysis tests the eight genes in
`results/final_analysis/ML_3_models_8_consensus_hub_genes.csv`: TNFSF15,
EBI3, PI3, FLT3LG, ADAM19, CLEC4F, TGM3, and CXCL14. It is separate from the
previous GSE62928 validation, which used a different eight-gene panel.

GSE62928 contains four EPS biopsies and four non-EPS samples (two PD and two
uremic controls). The deposited GPL13158 RMA log2 matrix was used without
additional normalization. For each gene, the probe with the highest mean
expression across all eight samples was selected without using group labels.
Gene-level EPS-versus-non-EPS effects and moderated P-values were estimated
with limma; P-values were adjusted across the eight genes with Benjamini-Hochberg.

The panel score is the mean of per-gene z-scores across the eight samples.
Discrimination is summarized by AUC, a stratified bootstrap interval (10,000
resamples), an exact two-sided Mann-Whitney test, and a matched random-gene-set
null (1,000 sets; seed 42). Null genes were drawn from genes above the
operational mean-expression threshold (>5), matched to expression deciles.

## Results

- Direction concordance: **2/8** (PI3 and ADAM19).
- Composite AUC: **0.375** (bootstrap 95% interval **0.000–0.750**).
- Composite exact Mann-Whitney P: **0.6857**.
- Expression-matched null AUC empirical P: **0.8871**.
- Mean panel-gene validation log2 fold change: **-0.1500**.
- Only **1/8** genes exceeded the operational mean-expression cutoff.
- ADAM19 had nominal limma P = **0.0160**, but its eight-gene-adjusted
  Benjamini-Hochberg P = **0.1281**.

The composite score did not distinguish EPS from non-EPS in this small cohort,
and its AUC was not unusual relative to the expression-matched null. These
results do **not** provide convincing external confirmation of the revised
eight-gene panel in GSE62928. With only four samples per group and seven genes
below the operational expression cutoff, estimates are highly uncertain.

This revised gene list was not preregistered for GSE62928; these results should
be treated as exploratory and not as a preregistered replication.

## Output files

- [Per-gene validation results](per_gene_validation.csv)
- [Panel-level validation summary](panel_validation_summary.csv)
- [Selected probes and expression cutoff](probe_selection.csv)
- [Sample-level composite scores](sample_composite_scores.csv)
- [Matched-null draws](matched_gene_set_null.csv)
- Analysis script: [`validate_ml_8_consensus_hubs_gse62928.R`](../../../../scripts/validation/validate_ml_8_consensus_hubs_gse62928.R)
