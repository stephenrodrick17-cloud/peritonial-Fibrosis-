# Revised 8-gene consensus-hub validation

This folder is the entry point for the discovery panel and its external
checks. Start with [the research report](report.md).

## Folder guide

- `GSE62928/`: per-gene limma results, panel score and matched-null results,
  probe choices, sample scores, and dataset-specific report.
- `GSE248762/`: per-gene edgeR results and sample-level group mapping.
- `GSE130888/`: per-gene edgeR results and sample-level group mapping.
- `effluent_cohort_comparison.md`: concise comparison of the two effluent
  cohorts.
- [Eight-gene panel CSV](../../final_analysis/ML_3_models_8_consensus_hub_genes.csv):
  the exact gene list and model votes tested here.
- [Validation scripts and guide](../../../scripts/validation/README.md):
  R code used to generate the GSE62928 and effluent cohort results.

Inputs for the analyses are the saved GEO matrix/annotation and the
patient-level pseudobulk files from `results/revision/effluent_validation/`.
The commands below run from the repository root in PowerShell. Set `$Rscript`
to the Rscript executable for your R installation.

```powershell
$Rscript = "C:\Program Files\R\R-4.4.2\bin\Rscript.exe"

& $Rscript --vanilla `
  "scripts\validation\validate_ml_8_consensus_hubs_gse62928.R" `
  "results\revision\gse62928_validation_8genes\GSE62928_series_matrix.txt.gz" `
  "data\GPL13158.soft.gz" `
  "results\final_analysis\ML_3_models_8_consensus_hub_genes.csv" `
  "results\revision\ml_8_consensus_hub_validation\GSE62928"

& $Rscript --vanilla `
  "scripts\validation\validate_ml_8_hubs_effluent_cohorts.R" `
  "GSE248762" `
  "results\revision\effluent_validation\GSE248762\analysis_20261009" `
  "results\final_analysis\ML_3_models_8_consensus_hub_genes.csv" `
  "results\revision\ml_8_consensus_hub_validation\GSE248762"

& $Rscript --vanilla `
  "scripts\validation\validate_ml_8_hubs_effluent_cohorts.R" `
  "GSE130888" `
  "results\revision\effluent_validation\GSE130888\analysis_20261009" `
  "results\final_analysis\ML_3_models_8_consensus_hub_genes.csv" `
  "results\revision\ml_8_consensus_hub_validation\GSE130888"
```

Required R packages: `limma` for GSE62928 and `edgeR` for the effluent
cohorts. The effluent commands consume already-generated patient pseudobulk
counts; they do not redo cell QC or pseudobulk construction.

These analyses are exploratory for this revised panel. Neither effluent
cohort detected every panel gene, so neither supports a full eight-gene
composite validation.
