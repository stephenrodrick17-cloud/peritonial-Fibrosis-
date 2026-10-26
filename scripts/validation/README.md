# Eight-hub validation scripts

The scripts in this folder validate the 8-gene consensus panel defined in
`results/final_analysis/ML_3_models_8_consensus_hub_genes.csv`.

- `validate_ml_8_consensus_hubs_gse62928.R` uses limma on the deposited
  GPL13158 RMA series matrix.
- `validate_ml_8_hubs_effluent_cohorts.R` uses edgeR on existing patient-level
  pseudobulk count, metadata, and gene-symbol mapping tables for GSE248762 or
  GSE130888. It does not redo single-cell QC or pseudobulk aggregation.

See `results/revision/ml_8_consensus_hub_validation/README.md` for runnable
PowerShell commands, inputs, and output locations.

The output directory also receives `sessionInfo.txt` to record the R and
package versions used for that run.
