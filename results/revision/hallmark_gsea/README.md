# GSE125498 MSigDB Hallmark preranked GSEA

This analysis uses the existing discovery-cohort limma result table and creates
new outputs only in this directory.

## Reproduction

Run from the project root with the installed R 4.4.2 executable:

```powershell
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" "results\revision\hallmark_gsea\01_rank_limma.R"
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" "results\revision\hallmark_gsea\02_fgsea_hallmark.R"
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" "results\revision\hallmark_gsea\03_sample_label_permutations.R"
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" "results\revision\hallmark_gsea\04_hallmark_figures.R"
```

## Data and analysis settings

- Cohort: GSE125498 discovery data, 33 samples (13 LPD, 20 SPD).
- Differential expression: reused `Validation/GSE125498.top.table.tsv`, the
  existing probe-level limma result. No covariates were included in the saved
  analysis; the contrast is LPD versus SPD, with LPD as the numerator.
- Ranking: moderated t-statistic in descending order. Duplicate gene symbols
  are collapsed by retaining the probe with the largest absolute moderated
  t-statistic. The resulting CSV includes all retained genes and probe IDs.
- Gene sets: the project-cached human `MSigDB_Hallmark_2020` collection,
  retrieved from Enrichr on 2026-10-11. This is the exact collection identifier
  available locally (50 Hallmark sets), rather than a newly downloaded
  current MSigDB release.
- GSEA: `fgseaSimple`, 10,000 gene-set permutations for observed results,
  seed 20261013, minimum/maximum set size 15/500, scoring exponent 1, BH FDR
  cutoff 0.05. A large but non-significant NES is defined descriptively as
  absolute NES at least 1.5 with FDR at least 0.05.
- Sample-label check: 1,000 random label permutations with group sizes held
  fixed (13 LPD, 20 SPD), seed 20261014. Limma `lmFit` and `eBayes` are rerun
  per assignment. GSEA is rerun for each assignment with 500 gene-set
  permutations to estimate NES; empirical P values compare absolute permuted
  and observed NES for each set. These are random permutations, not exhaustive
  enumeration.
- Software versions and exact run settings are recorded in
  `hallmark_gsea_metadata.csv` and `hallmark_sample_label_permutation_metadata.csv`.
  The analysis was run on 2026-10-13 with R 4.4.2, limma 3.62.2, and
  fgsea 1.32.4.

## Outputs

- `hallmark_ranked_genes.csv`: unique gene ranking from existing limma table.
- `hallmark_gsea_all.csv`: all tested Hallmark sets, including NES, P, FDR,
  set size, and leading-edge genes.
- `hallmark_gsea_upregulated.csv` and `hallmark_gsea_downregulated.csv`:
  significant sets separated by NES direction. Positive NES indicates
  enrichment in LPD; negative NES indicates enrichment in SPD.
- `hallmark_gsea_large_NES_non_significant.csv`: sets meeting the prespecified
  descriptive large-NES threshold but not the FDR significance cutoff.
- `hallmark_sample_label_permutation_NES.csv`: pathway NES values for every
  random label permutation.
- `hallmark_sample_label_empirical_P.csv` and
  `hallmark_sample_label_top_set_rank.csv`: empirical pathway P values and
  rank of the observed top pathway relative to permuted top sets.
- Figures: significant-set NES bar plot and top significant up/down
  enrichment plots, each in PNG and PDF format.

## Limitations

Results are exploratory and hypothesis-generating because the cohort is small.
The locally available Hallmark collection is identified by Enrichr as
`MSigDB_Hallmark_2020`; no later MSigDB release was downloaded. For sample-label
permutations, the saved gene-level normalized expression matrix has one
representative row per symbol, whereas the observed ranking reuses the
probe-level limma table and collapses probes by absolute moderated t-statistic.
Accordingly, the label-permutation check is a sensitivity analysis of the
limma/GSEA pipeline using the available gene-level matrix, not an exact
probe-level re-fit of the original analysis. fgsea also reports a small
fraction of tied moderated-t statistics (0.01%); ties can make their relative
ordering arbitrary.
