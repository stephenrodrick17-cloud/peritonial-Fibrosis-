# GSE125498 immune/stromal marker-score analysis

Run date: 2026-10-13. All new analysis outputs are written to this directory;
existing project tables and figures are read-only inputs.

Software used: Python 3.13.1, pandas 2.3.3, NumPy 2.4.6, SciPy 1.18.1,
statsmodels 0.15.0, seaborn 0.13.2, and Matplotlib 3.11.2.

## Reproduction

From the project root, use the configured project Python environment:

```powershell
& ".venv\Scripts\python.exe" "results\revision\microenvironment_deconvolution\01_build_signatures.py"
& ".venv\Scripts\python.exe" "results\revision\microenvironment_deconvolution\02_score_signatures.py"
& ".venv\Scripts\python.exe" "results\revision\microenvironment_deconvolution\03_compare_and_permute.py"
& ".venv\Scripts\python.exe" "results\revision\microenvironment_deconvolution\04_make_figures.py"
```

The scripts use the existing 11,741-symbol normalized gene-level expression
matrix and sample metadata for GSE125498 (13 LPD and 20 SPD).

## Signatures and scoring

The source is the human **CellMarker_2024** curated marker collection served
through Enrichr (https://maayanlab.cloud/Enrichr/), retrieved on the run date.
The underlying CellMarker database is described by Hu et al., *CellMarker 2.0:
an updated database of manually curated cell markers in human/mouse and web
tools based on scRNA-seq data*, Nucleic Acids Research 51(D1), D870-D876
(2023), doi:10.1093/nar/gkac947. Population signatures were taken from the
exact CellMarker_2024 terms listed in the source manifest; signatures from
other tissues are database reference markers, not peritoneal-effluent-specific
validated classifiers.
Each source signature was selected by its exact database term; the selected
terms and detectable marker counts are in `signature_source_manifest.csv` and
`signature_marker_counts.csv`. The downloaded library snapshot is retained as
`cellmarker_2024_library_snapshot.txt`. Population-specific markers and source
terms are exported to `final_signature_table_both_variants.csv`.

Thirteen available populations are scored: T cells, B cells, NK cells,
monocytes, M1 macrophages, M2 macrophages, neutrophils, mast cells, immature
dendritic cells, mature dendritic cells, fibroblasts, myofibroblasts, and
mesothelial cells. Populations with fewer than three detectable markers are
excluded and reported in the marker-count table.

For every population, the score is the mean of its markers' row-wise z-scores
across the 33 samples (sample standard deviation, `ddof=1`). These are relative
marker-expression scores, not estimated cell proportions.

The signature audit found no detectable marker overlap with the five-gene
nomogram panel and one overlap with the seven-gene consensus panel: **LTB**
in the T-cell signature. Following user approval, outputs include both the
full signatures and a sensitivity version excluding LTB. The source markers
are not otherwise modified.

## Comparisons and permutation check

Each population is compared between LPD and SPD by a two-sided Mann–Whitney U
test. The effect size is rank-biserial correlation,
`2*U/(n_LPD*n_SPD)-1`, where positive values indicate higher scores in LPD.
BH correction is performed across populations separately within each
signature variant.

The label-permutation check uses 10,000 random permutations with fixed group
sizes (13 LPD and 20 SPD) and seed 20261015. For population-specific empirical
P values, the absolute rank-biserial statistic is compared with its observed
value. Family-wise permutation P is the proportion of random permutations
with at least one population having nominal Mann–Whitney P < 0.05. These are
random permutations, not exhaustive enumeration. The full
`sample_label_permutation_statistics.csv` export records one row per
population per permutation.

## Outputs and interpretation

- `signature_marker_counts.csv`, `signature_marker_counts_by_variant.csv`,
  `signature_panel_overlaps.csv`, and `final_signature_table_both_variants.csv`
  describe signature provenance, detectable counts, and overlap decisions.
- `sample_population_scores_long.csv` and the two
  `sample_population_score_matrix_*.csv` files contain the scores.
- `group_comparisons.csv`, `nominal_P_lt_0_05.csv`, `FDR_lt_0_05.csv`,
  `sample_label_permutation_statistics.csv`, and
  `familywise_permutation_summary.csv` contain test and permutation results.
- Boxplots and sample heatmaps are exported in PNG and PDF for both signature
  variants.

The cohort is small and the marker scores are indirect estimates from
peritoneal effluent. Stromal and mesothelial scores are indirect signals, not
direct measures of tissue composition. Results are exploratory and
hypothesis-generating; they do not establish cell abundance, mechanism, or
causality.
