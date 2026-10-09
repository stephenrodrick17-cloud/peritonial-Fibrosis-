# Fibrosis-context analysis preregistration

**Cohort:** GSE248762 only. Use observed data and random seed 42. The existing
cohort inclusion and QC decisions remain unchanged.

## Fixed feature sets

- **Candidate genes (8):** TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1,
  SERPINA10, CST7
- **Fibrosis / remodeling markers (8):** ACTA2, COL1A1, COL1A2, FN1, TGFB1,
  CCN2, SERPINE1, VEGFA
- **Composition covariates (2):** per-patient monocyte/macrophage fraction and
  T/NK fraction

## Analyses

1. Aggregate raw counts by patient for all QC-passing cells and separately by
   each available major cell type. Use patient as the independent observation.
   For each candidate-marker pair, report n, Spearman correlation, 95% CI,
   two-sided P, and BH-adjusted q. The composition-adjusted estimate is the
   Pearson correlation of rank-transformed candidate and marker expression
   residuals after separate linear adjustment for the two preregistered
   fractions. Its P value uses the partial-correlation t approximation with
   residual degrees of freedom `n - 2 - 2`. Obtain percentile 95% CIs from
   2,000 patient-level bootstrap resamples (seed 42); resamples with a constant
   residual are non-estimable. BH correction is applied across the 64
   candidate-marker tests within each aggregation stratum and adjustment type.
   Report unavailable or non-estimable tests rather than imputing values.
2. Compare LV_UF with LV_NOT_UF for the eight candidate genes using the
   existing GSE248762 edgeR pseudobulk pipeline and its previously fixed sample
   inclusion, filtering, normalization, and model settings. Do not add SV
   samples to this contrast.
3. Rank genes separately in macrophage/T-cell and fibroblast/mesothelial
   compartments using the same LV_UF versus LV_NOT_UF direction. Test the
   fixed Hallmark EMT, TGF-beta, ECM-organization, and inflammatory-response
   gene sets; report NES, nominal P, FDR, and which of the eight candidate genes
   occur in each leading edge. Do not substitute other gene sets or cohorts.
4. Estimate ligand-receptor communication with LIANA or CellChat from observed
   cell expression only. If resource limits require subsampling, use at most
   20,000 cells with seed 42 and record the sampling. Report TNFSF15-TNFRSF25,
   ADAM19/Notch, and IL-27 signaling only as **inferred** communication, with
   the selected tool, version, database, and scored interactions. Do not
   interpret predictions as measured signaling.
5. Calculate the fixed MMT score in mesothelial and fibroblast cells by PD
   group using the marker set above. State the scoring formula, sample sizes,
   and uncertainty; do not interpret group differences causally.

## Frozen success rule

A candidate gene is **supported** only if both conditions hold:

1. After adjustment for monocyte/macrophage and T/NK fractions, its correlation
   with at least 3 of the 8 fixed markers has a positive Spearman rho > 0.4;
   and
2. Its direction is concordant in the LV_UF versus LV_NOT_UF contrast.

Evaluate the correlation-count condition in the **all-cell patient pseudobulk**
stratum; the cell-type-stratified correlations are contextual/exploratory and
cannot independently confer support. For this directional rule, concordance
means a positive edgeR log2 fold change for LV_UF relative to LV_NOT_UF
(higher candidate expression in LV_UF).

Report both conditions per gene. No gene, marker, covariate, contrast, threshold,
or success criterion may be changed after results are inspected. All analyses
are correlative and cannot establish causality.

## Fixed MMT-score operationalization

Because the request supplied no separate canonical MMT signature or weights, the
MMT measure is defined here, before calculating scores, as a **study-specific
fibrosis/MMT proxy**: within each patient and cell type (mesothelial or
fibroblast), calculate the mean of the eight marker genes' per-cell log1p
CPM-like normalized expression after gene-wise z-scoring across the analyzed
cells in that cell type. A score is summarized per patient and compared by the
reported PD group; it is not represented as a validated canonical MMT score.
