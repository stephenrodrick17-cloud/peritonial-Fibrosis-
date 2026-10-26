# GSE130888 preregistered analysis (Steps 2–4)

**Preregistered verdict for long-term vs short-term PD: Partially supportive.**

- Primary-comparison pseudobulks: long-term PD 4, short-term PD 6; normal controls 3 (descriptive only).
- Detected genes direction-concordant: 5/6.
- TNFSF15 and SERPINA10 were not detected by the frozen raw-count rule; neither is concordant nor discordant, and both are excluded from the concordance denominator.
- Full-panel composite: not estimable under the frozen rule because one or more of the eight genes were not detected.
- Full-panel matched-null empirical P: not estimable because one or more fixed panel genes were not detected; the frozen full-panel composite was not computed.
- No full-panel composite ROC or matched-null histogram was drawn because the fixed eight-gene composite is not estimable.
- Strict hubs: 1/1 detected genes concordant.
- Strict-hub composite/null: not estimable under the preregistered detection or decile requirements.

## Per-gene results

Discovery confidence intervals are not present in the authoritative discovery results CSV, so discovery log2FC points are shown without intervals in `discovery_vs_gse130888_forest.png`; GSE130888 intervals are shown from edgeR.

| gene      |   discovery_log2FC |   discovery_P | detection_status   | log2FC             | ci_lower          | ci_upper         | P                 | BH_q_across_tested_panel_genes   | direction_concordant   |
|:----------|-------------------:|--------------:|:-------------------|:-------------------|:------------------|:-----------------|:------------------|:---------------------------------|:-----------------------|
| TNFSF15   |           0.97431  |      2.73e-05 | not detected       | not estimated      | not estimated     | not estimated    | not estimated     | not estimated                    | not estimated          |
| FLT3LG    |           1.43771  |      2.35e-05 | detected           | 0.825211703493845  | -1.70111331209599 | 3.35153671908368 | 0.309631476085139 | 0.5516030674105455               | True                   |
| EBI3      |           0.803423 |      0.00244  | detected           | 0.838215750135198  | -1.07155942645569 | 2.74799092672609 | 0.194049883116702 | 0.5516030674105455               | True                   |
| LTB       |           1.67472  |      3.33e-06 | detected           | 0.196322794440222  | -2.7174331973353  | 3.11007878621574 | 0.831029482474644 | 0.8310294824746439               | True                   |
| ADAM19    |           1.21309  |      0.000133 | detected           | 0.923061238458685  | -1.22313726520747 | 3.06925974212484 | 0.188312230056285 | 0.5516030674105455               | True                   |
| ECM1      |           1.03405  |      0.00701  | detected           | -0.630148371837476 | -2.73247230997387 | 1.47217556629892 | 0.367735378273697 | 0.5516030674105455               | False                  |
| SERPINA10 |           1.16585  |      0.000653 | not detected       | not estimated      | not estimated     | not estimated    | not estimated     | not estimated                    | not estimated          |
| CST7      |           1.01227  |      0.0005   | detected           | 0.158321516991908  | -2.03750032234611 | 2.35414335632993 | 0.819583783502453 | 0.8310294824746439               | True                   |

## Normal controls (descriptive only)

Normal-control rows show raw sample-level pseudobulk counts and TMM-normalized log-CPM; no inferential test uses the 3 normal controls.

| gene      | gene_id         | sample_id   |   raw_pseudobulk_count |    logCPM |
|:----------|:----------------|:------------|-----------------------:|----------:|
| TNFSF15   | ENSG00000181634 | GSM3755693  |                    104 |  0.904623 |
| TNFSF15   | ENSG00000181634 | GSM3755694  |                     49 |  0.87958  |
| TNFSF15   | ENSG00000181634 | GSM3755695  |                     55 |  0.867646 |
| FLT3LG    | ENSG00000090554 | GSM3755693  |                    185 |  1.71638  |
| FLT3LG    | ENSG00000090554 | GSM3755694  |                    160 |  2.55573  |
| FLT3LG    | ENSG00000090554 | GSM3755695  |                    169 |  2.45672  |
| EBI3      | ENSG00000105246 | GSM3755693  |                     12 | -1.91061  |
| EBI3      | ENSG00000105246 | GSM3755694  |                     46 |  0.791318 |
| EBI3      | ENSG00000105246 | GSM3755695  |                      0 | -4.14579  |
| LTB       | ENSG00000227507 | GSM3755693  |                     73 |  0.412379 |
| LTB       | ENSG00000227507 | GSM3755694  |                    175 |  2.68383  |
| LTB       | ENSG00000227507 | GSM3755695  |                    504 |  4.02321  |
| ADAM19    | ENSG00000135074 | GSM3755693  |                    283 |  2.32104  |
| ADAM19    | ENSG00000135074 | GSM3755694  |                    161 |  2.56464  |
| ADAM19    | ENSG00000135074 | GSM3755695  |                    138 |  2.16769  |
| ECM1      | ENSG00000143369 | GSM3755693  |                    784 |  3.78063  |
| ECM1      | ENSG00000143369 | GSM3755694  |                    513 |  4.22705  |
| ECM1      | ENSG00000143369 | GSM3755695  |                   1003 |  5.01355  |
| SERPINA10 | ENSG00000140093 | GSM3755693  |                      0 | -4.14579  |
| SERPINA10 | ENSG00000140093 | GSM3755694  |                      0 | -4.14579  |
| SERPINA10 | ENSG00000140093 | GSM3755695  |                      0 | -4.14579  |
| CST7      | ENSG00000077984 | GSM3755693  |                     26 | -0.970347 |
| CST7      | ENSG00000077984 | GSM3755694  |                    570 |  4.37862  |
| CST7      | ENSG00000077984 | GSM3755695  |                     93 |  1.60708  |

## Expression by cell type

The top cell type is selected by mean normalized expression across all QC-retained cells, without group labels. The full gene-by-cell-type profile is saved in `gene_expression_by_cell_type.csv`.

| gene      | gene_top_cell_type   |   cells |   fraction_cells_nonzero |   mean_log1p_normalized_expression |
|:----------|:---------------------|--------:|-------------------------:|-----------------------------------:|
| TNFSF15   | fibroblast           |    5199 |               0.00788613 |                         0.00910368 |
| FLT3LG    | T/NK                 |   24596 |               0.219751   |                         0.353369   |
| EBI3      | B                    |    9557 |               0.0169509  |                         0.0220527  |
| LTB       | T/NK                 |   24596 |               0.617865   |                         1.4723     |
| ADAM19    | B                    |    9557 |               0.15915    |                         0.21799    |
| ECM1      | fibroblast           |    5199 |               0.21658    |                         0.322578   |
| SERPINA10 | T/NK                 |   24596 |               0          |                         0          |
| CST7      | T/NK                 |   24596 |               0.513498   |                         1.07305    |

Top-cell-type assignments are descriptive marker-based annotations; a relative maximum is not evidence of meaningful expression.

## QC and analysis notes

Cells were retained using the preregistered group-blind thresholds of at least 200 detected genes and at most 20% mitochondrial counts. The GEO groups were short-term PD n=6, long-term PD n=4, and normal controls n=3. Samples passing the blind QC were aggregated per GEO sample; each GSM is used as an anonymous patient proxy because participant IDs are not exposed in the GEO metadata. A panel gene was detected only if raw counts were at least 10 in at least half of patient pseudobulks in each group of the primary contrast. Raw counts were aggregated per sample; cells were not treated as independent replicates. The one-GSM-per-person assumption cannot be confirmed from the available GEO identifiers.

Cell types use the fixed preregistered marker-score rules because cell-level author annotations were not present in the downloaded 10x archive. All resampling and permutation-style tests use seed 42. The composite-AUC random-gene null is used for the supportive decision rule; the preregistered mean-log2FC null is also reported.

edgeR P values are quasi-likelihood F-test results. Approximate 95% intervals use the quasi-likelihood variance-scaled GLM coefficient covariance and residual-adjusted degrees of freedom, expressed on the log2 scale.

## Composition

| cell_type           |   n_LV |   n_SV |   mean_fraction_LV |   mean_fraction_SV |   wilcoxon_p |     bh_q |
|:--------------------|-------:|-------:|-------------------:|-------------------:|-------------:|---------:|
| T/NK                |      4 |      6 |         0.42057    |         0.215542   |    0.257143  | 0.6      |
| B                   |      4 |      6 |         0.120779   |         0.083425   |    0.257143  | 0.6      |
| monocyte/macrophage |      4 |      6 |         0.231291   |         0.543981   |    0.0666667 | 0.466667 |
| neutrophil          |      4 |      6 |         0.0575042  |         0.0618739  |    0.914286  | 0.914286 |
| mesothelial         |      4 |      6 |         0.151492   |         0.0923816  |    0.914286  | 0.914286 |
| fibroblast          |      4 |      6 |         0.00911075 |         0.00101241 |    0.609524  | 0.853333 |
| other               |      4 |      6 |         0.00925391 |         0.00178457 |    0.47619   | 0.833333 |

## Composition-adjusted per-gene results

When a gene is not detected under the preregistered raw-count rule, its adjusted estimate and P value are suppressed as not detected. Any detected gene whose adjusted coefficient is no longer positive is described as reflecting immune-cell abundance.

| gene      | gene_id         | detected   |   discovery_log2FC | contrast      | gene_symbol   | logFC         | ci_lower      | ci_upper      | PValue        | detection_status   | FDR_all_genes   | F             | direction_concordant   | BH_q_across_tested_panel_genes   |
|:----------|:----------------|:-----------|-------------------:|:--------------|:--------------|:--------------|:--------------|:--------------|:--------------|:-------------------|:----------------|:--------------|:-----------------------|:---------------------------------|
| TNFSF15   | ENSG00000181634 | False      |           0.97431  | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | not detected       | not estimated   | not estimated | not estimated          | not estimated                    |
| FLT3LG    | ENSG00000090554 | True       |           1.43771  | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | detected           | not estimated   | not estimated | not estimated          | not estimated                    |
| EBI3      | ENSG00000105246 | True       |           0.803423 | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | detected           | not estimated   | not estimated | not estimated          | not estimated                    |
| LTB       | ENSG00000227507 | True       |           1.67472  | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | detected           | not estimated   | not estimated | not estimated          | not estimated                    |
| ADAM19    | ENSG00000135074 | True       |           1.21309  | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | detected           | not estimated   | not estimated | not estimated          | not estimated                    |
| ECM1      | ENSG00000143369 | True       |           1.03405  | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | detected           | not estimated   | not estimated | not estimated          | not estimated                    |
| SERPINA10 | ENSG00000140093 | False      |           1.16585  | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | not detected       | not estimated   | not estimated | not estimated          | not estimated                    |
| CST7      | ENSG00000077984 | True       |           1.01227  | not estimated | not estimated | not estimated | not estimated | not estimated | not estimated | detected           | not estimated   | not estimated | not estimated          | not estimated                    |

Composition adjustment was not feasible under the preregistered minimum sample-size rule; direction changes after adjustment cannot be assessed.

## Top-cell-type per-gene results

| gene      | gene_id         | top_cell_type   | status                        |   n_LV |   n_SV |   cells_min_count_LV |   cells_min_count_SV |      logFC |   ci_lower |    ci_upper |      PValue |   BH_q_across_panel_genes |
|:----------|:----------------|:----------------|:------------------------------|-------:|-------:|---------------------:|---------------------:|-----------:|-----------:|------------:|------------:|--------------------------:|
| TNFSF15   | ENSG00000181634 | fibroblast      | not detected in top cell type |      2 |      0 |                    0 |                    0 | nan        | nan        | nan         | nan         |               nan         |
| FLT3LG    | ENSG00000090554 | T/NK            | tested                        |      4 |      6 |                    4 |                    6 |  -0.159257 |  -0.908544 |   0.590029  |   0.603186  |                 0.603186  |
| EBI3      | ENSG00000105246 | B               | tested                        |      4 |      6 |                    4 |                    3 |   1.90449  |   0.382456 |   3.42652   |   0.0138036 |                 0.0527364 |
| LTB       | ENSG00000227507 | T/NK            | tested                        |      4 |      6 |                    4 |                    6 |  -0.969224 |  -1.87808  |  -0.0603693 |   0.0210946 |                 0.0527364 |
| ADAM19    | ENSG00000135074 | B               | tested                        |      4 |      6 |                    4 |                    6 |   0.568586 |  -0.363753 |   1.50093   |   0.17747   |                 0.221838  |
| ECM1      | ENSG00000143369 | fibroblast      | not estimable by edgeR        |      2 |      0 |                    2 |                    0 | nan        | nan        | nan         | nan         |               nan         |
| SERPINA10 | ENSG00000140093 | T/NK            | not detected in top cell type |      4 |      6 |                    0 |                    0 | nan        | nan        | nan         | nan         |               nan         |
| CST7      | ENSG00000077984 | T/NK            | tested                        |      4 |      6 |                    4 |                    6 |  -0.498293 |  -1.17936  |   0.182774  |   0.0907849 |                 0.151308  |

## Software versions

- Python 3.13.1
- Scanpy 1.12.4; AnnData 0.13.4.
- R and edgeR versions are recorded in `edger_analysis_config.csv`.

This is the GSE130888 effluent-cell cohort. The preregistered independence assessment labels its patient overlap with GSE248762 unknown; it is also not formally confirmed independent of discovery GSE125498. See `../independence_assessment.md`. The small primary groups, single-cell dropout, lack of participant identifiers, and wide intervals limit interpretation; wide intervals are inconclusive, not proof of no effect. Composition adjustment is not feasible with four long-term and six short-term samples under the preregistered minimum of five samples per group (long-term n=4, short-term n=6). No GSE92455 analysis was performed.
