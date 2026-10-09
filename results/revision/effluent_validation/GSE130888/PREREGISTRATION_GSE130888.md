# GSE130888 long-term versus short-term validation

**Status:** fixed before rerunning the accession-specific pipeline.
**Cohort:** real GSE130888 data only; seed 42.

## Contrast and fixed panel

Compare the 4 GEO-labeled **Long-term PD** patient samples with the 6
**Short-term PD** patient samples. The three GEO-labeled normal controls are
descriptive only and are excluded from this contrast.

Test the unchanged eight-gene panel:
`TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, SERPINA10, CST7`.
The strict 3/3 subset remains `TNFSF15, FLT3LG`. No sample, gene, contrast, or
threshold will be changed after inspecting results.

## Pipeline and detection

Run the existing GSE130888 Scanpy/edgeR pipeline with its existing group-blind
QC, fixed cell annotation, patient-level raw-count pseudobulk, edgeR
quasi-likelihood model, and preprocessing. Use a new output directory under
`results/revision/fibrosis_context/GSE130888/`; do not overwrite
`analysis_20261009`.

A panel gene is detected only when its raw pseudobulk count is at least 10 in
at least half of independent patient pseudobulk samples in **each** comparison
group. An undetected gene is neither concordant nor discordant and is excluded
from the concordance denominator. Report the per-gene log2 fold change,
95% CI, nominal P, BH-adjusted P, detection status, and direction.

## Frozen success criteria

Apply the already preregistered GSE effluent validation tiers without
modification:

- **Supportive:** at least 6 of 8 detected genes concordant, the full fixed
  eight-gene composite is estimable and has AUC at least 0.80, its matched-null
  empirical P is below 0.05, and both strict genes are concordant.
- **Partially supportive:** exactly 5 detected genes are concordant, or an
  estimable composite AUC is 0.65 through 0.80 without matched-null
  significance.
- **Not supported:** otherwise.

If any panel gene is undetected, the fixed eight-gene composite and its
matched-null AUC test remain non-estimable; do not replace missing genes with a
detected-gene score for the preregistered verdict. The existing pipeline may
report other fixed diagnostics, but they do not alter these tiers.

No simulated data, alternate gene sets, new cutoffs, sample exclusions, or
post-result rule changes are permitted. Stop and report any input-format,
sample-mapping, or resource problem rather than adapting the analysis.
