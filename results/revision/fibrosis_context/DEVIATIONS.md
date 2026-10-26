# Fibrosis-context analysis deviations and post-run corrections

This log records corrections to `run_fibrosis_context.py` and
`compartment_edger_gsea_input.R` after the initial GSE248762 fibrosis-context
execution. The existing raw and pseudobulk inputs, cohort membership, gene and
marker panels, contrasts, composition covariates, and frozen success rule were
not changed. No new fibrosis-context analysis was run while preparing this log.

| Script / edit | Reason | Effect on reported numbers |
|---|---|---|
| `run_fibrosis_context.py`: replaced 2,000 serial bootstrap-resample correlation calls with a batched/vectorized patient bootstrap. | The first run spent an impractical amount of time in the per-resample Python/SciPy loop. The vectorized calculation applies the same patient-resampling and rank/residual-correlation procedure. | No intended change to the estimand or point estimates. Bootstrap interval values may differ only if the resampling sequence/numerical evaluation differs; the final saved intervals are those in `gene_marker_correlations.csv`. The partial first run did not produce a complete comparison set. |
| `run_fibrosis_context.py`: used a Moore-Penrose pseudoinverse for batched bootstrap least-squares residualization. | Some bootstrap resamples make the rank-transformed covariate design singular (for example, a covariate becomes constant); a direct matrix solve raised `LinAlgError` and stopped execution. | No change to the defined residual-correlation estimate; singular resamples are evaluated with the least-squares pseudoinverse and constant residual results remain non-estimable. It enabled valid resamples to finish. |
| `run_fibrosis_context.py`: added progress messages for each aggregation/cell-type stratum. | The lengthy patient-level correlation calculation otherwise provided no indication that it was progressing. | None. |
| `run_fibrosis_context.py`: explicitly labeled constant-input unadjusted tests and constant-residual adjusted tests as non-estimable. | Scipy correctly returns an undefined correlation when a gene has no variation, but earlier status handling did not distinguish those rows from ordinary estimates. | The 72 rows now have explicit non-estimable statuses; no finite correlation estimate was replaced. All 72 are associated with zero raw candidate-gene counts in that patient-by-cell-type stratum, hence zero variance. There were 48 constant-input rows and 24 constant-residual rows. None was due to fewer than four patients; the contributing strata had n=9, 10, or 16. |
| `run_fibrosis_context.py`: recorded Hallmark/GO gene-set membership in `gsea_gene_set_membership.csv`. | Makes the exact gene-set contents used by preranked enrichment reviewable and reproducible. | None to the rankings or enrichment. |
| `run_fibrosis_context.py`: represented zero-permutation-exceedance GSEA P/FDR values at the empirical resolution floor, 1/1,001, rather than exact zero. | The preranked tests use 1,000 permutations; zero sampled exceedances do not establish a probability of exactly zero. | Yes: any reported zero P or FDR from the 1,000-permutation calculation is displayed at the minimum resolvable value (approximately 0.000999). NES, leading edges, and permutation procedure are unchanged. |
| `run_fibrosis_context.py`: added input pseudobulk, metadata, composition, and gene-map copies with SHA-256 digests in the manifest. | Preserve the exact patient-level analysis inputs alongside the outputs without copying raw cell matrices. | None. |
| `run_fibrosis_context.py`: clarified LIANA route summaries, best-ranked source/target examples, and the absence of a direct ADAM19-Notch interaction in the selected resource. | The initial report could be read as affirmative evidence of the specifically requested signaling routes. In fact, the consensus resource lacks a direct ADAM19-Notch pair, and immune-to-stromal TNFSF15-TNFRSF25 and IL-27 candidates had uninformative rank/P outputs. | No LIANA scores changed. The report interpretation was narrowed to what the saved inferred scores and selected database support. |
| `run_fibrosis_context.py`: displayed the saved MMT group-summary table, including its Kruskal-Wallis P value, rather than recomputing a separate report table. | Keep the report synchronized with the saved MMT results. | No MMT scores or tests changed. |
| `compartment_edger_gsea_input.R`: no post-run edits. | The cell-type pseudobulk edgeR design, contrast, filtering, TMM normalization, and robust quasi-likelihood test were unchanged after its successful run. | None. |

The 1,024 correlation rows comprise 8 candidate genes × 8 markers × 8
aggregation strata × 2 adjustment types. Exactly 72 were non-estimable:
48 unadjusted tests had zero-variance candidate-gene input, and 24
composition-adjusted tests had a zero-variance residual for that same
no-count candidate gene. These are detection/zero-variance cases, not
insufficient-patient cases. In the contributing cell-type strata, candidate
gene raw pseudobulk counts were zero across the included patients; the all-cell
analyses and other cell types remain in the results. No additional filtering
or imputation was applied.
