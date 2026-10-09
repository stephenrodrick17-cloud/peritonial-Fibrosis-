# GSE248762 preregistered analysis (Steps 2–4)

**Preregistered verdict for LV vs SV: Not supported.**

- Pseudobulk patients: LV 10, SV 6.
- Detected genes direction-concordant: 7/7.
- SERPINA10 was not detected by the frozen raw-count rule; it is neither concordant nor discordant and is excluded from that denominator.
- Full-panel composite: not estimable under the frozen rule because one or more of the eight genes were not detected.
- Verdict interpretation: 7/7 detected genes are concordant, so the supportive direction-count component (at least 6) is met. The fixed composite and its null are not estimable because SERPINA10 is not detected. The partial tier's exact-five branch does not apply, and its AUC branch cannot be evaluated. Thus the literal fallback is Not supported; this non-monotonic label does not mean the directional result is weaker than a 5/6 result.
- Full-panel matched-null empirical P: not estimable because SERPINA10 was not detected; the frozen full-panel composite was not computed.
- Strict hubs: 2/2 detected genes concordant.
- Strict-hub composite (secondary): AUC 0.7167 (95% bootstrap CI 0.4167–1.0000); Mann–Whitney P 0.180569; matched-null empirical P for AUC 0.1968 and mean gene log2FC 0.0470.

## Per-gene results

| gene      |   discovery_log2FC |   discovery_P | detection_status   | log2FC            | ci_lower            | ci_upper         | P                  | BH_q_across_tested_panel_genes   | direction_concordant   |
|:----------|-------------------:|--------------:|:-------------------|:------------------|:--------------------|:-----------------|:-------------------|:---------------------------------|:-----------------------|
| TNFSF15   |           0.97431  |      2.73e-05 | detected           | 1.3234540510612   | -0.624306466535089  | 3.27121456865749 | 0.106868408807977  | 0.1496157723311678               | True                   |
| FLT3LG    |           1.43771  |      2.35e-05 | detected           | 1.10484599970192  | -0.478582915383247  | 2.68827491478709 | 0.0925255419588824 | 0.1496157723311678               | True                   |
| EBI3      |           0.803423 |      0.00244  | detected           | 1.37070669067421  | 0.235600737105744   | 2.50581264424269 | 0.0074523990422494 | 0.026083396647872903             | True                   |
| LTB       |           1.67472  |      3.33e-06 | detected           | 0.820745042342153 | -0.705167451838138  | 2.34665753652244 | 0.181709779778873  | 0.18382667134536                 | True                   |
| ADAM19    |           1.21309  |      0.000133 | detected           | 1.06774162883591  | -0.0593807925721233 | 2.19486405024394 | 0.0268239173527543 | 0.06258914048976004              | True                   |
| ECM1      |           1.03405  |      0.00701  | detected           | 5.88930480921725  | 2.29769471704987    | 9.48091490138464 | 0.0048560694739974 | 0.026083396647872903             | True                   |
| SERPINA10 |           1.16585  |      0.000653 | not detected       | not estimated     | not estimated       | not estimated    | not estimated      | not estimated                    | not estimated          |
| CST7      |           1.01227  |      0.0005   | detected           | 0.452212511249147 | -0.405105544159116  | 1.30953056665741 | 0.18382667134536   | 0.18382667134536                 | True                   |

## Per-patient count and influence audit

Counts are raw pseudobulk sums over QC-retained cells; expressing-cell counts use raw count > 0. Fibroblast count and ECM1-expressing fibroblasts use the preregistered marker annotation.

| sample_id   | group     |   raw_pseudobulk_ECM1 |   fibroblast_cells |   fibroblasts_expressing_ECM1 |   all_cells_expressing_ECM1 |   raw_pseudobulk_TNFSF15 |   cells_expressing_TNFSF15 |   raw_pseudobulk_SERPINA10 |   cells_expressing_SERPINA10 |
|:------------|:----------|----------------------:|-------------------:|------------------------------:|----------------------------:|-------------------------:|---------------------------:|---------------------------:|-----------------------------:|
| GSM7919583  | LV_UF     |                   142 |                 30 |                             8 |                          65 |                       43 |                         25 |                          0 |                            0 |
| GSM7919584  | LV_UF     |                115870 |                  5 |                             3 |                        1331 |                       23 |                         15 |                          1 |                            1 |
| GSM7919585  | LV_UF     |                 72808 |                 39 |                            18 |                         731 |                      214 |                         72 |                          1 |                            1 |
| GSM7919586  | LV_UF     |                   183 |                  3 |                             1 |                          74 |                        0 |                          0 |                          0 |                            0 |
| GSM7919587  | LV_NOT_UF |                 20334 |                832 |                           710 |                        1491 |                       73 |                         60 |                          0 |                            0 |
| GSM7919588  | LV_NOT_UF |                   635 |                 91 |                            32 |                         255 |                       33 |                         26 |                          0 |                            0 |
| GSM7919589  | LV_NOT_UF |                   961 |                  6 |                             2 |                         489 |                       48 |                         35 |                          0 |                            0 |
| GSM7919590  | LV_NOT_UF |                   199 |                  1 |                             0 |                          98 |                       43 |                         21 |                          0 |                            0 |
| GSM7919591  | LV_NOT_UF |                   590 |                  4 |                             2 |                          83 |                       87 |                         18 |                          0 |                            0 |
| GSM7919592  | LV_NOT_UF |                   468 |                 19 |                             2 |                         223 |                       21 |                         16 |                          0 |                            0 |
| GSM7919593  | SV        |                   469 |                 29 |                             6 |                         305 |                       58 |                         42 |                          0 |                            0 |
| GSM7919594  | SV        |                   346 |                 27 |                             3 |                         222 |                       97 |                         73 |                          0 |                            0 |
| GSM7919595  | SV        |                   552 |                 86 |                            35 |                         238 |                       30 |                         24 |                          0 |                            0 |
| GSM7919596  | SV        |                   749 |                 26 |                             7 |                         510 |                       29 |                         23 |                          0 |                            0 |
| GSM7919597  | SV        |                  1306 |                 77 |                            34 |                         475 |                       23 |                         21 |                          0 |                            0 |
| GSM7919598  | SV        |                   262 |                  3 |                             2 |                         198 |                        7 |                          7 |                          0 |                            0 |

ECM1 full-model log2FC was 5.8893. The most influential single patient was GSM7919584 (LV_UF); excluding that patient gave log2FC 5.0152. All leave-one-out estimates remained positive (5.0152 to 6.0993).

## Secondary contrasts (descriptive)

Detection is evaluated separately for each contrast using the same raw-count threshold; undetected genes are shown without estimates.

| contrast           | gene      | detection_status   | log2FC              | ci_lower           | ci_upper         | P                  | BH_q_across_tested_panel_genes   | direction_concordant   |
|:-------------------|:----------|:-------------------|:--------------------|:-------------------|:-----------------|:-------------------|:---------------------------------|:-----------------------|
| LV_UF_vs_LV_NOT_UF | TNFSF15   | detected           | 1.47741747508545    | -0.928258714966521 | 3.88309366513742 | 0.0991960447927488 | 0.3471861567746208               | True                   |
| LV_UF_vs_LV_NOT_UF | FLT3LG    | detected           | -0.548253732814737  | -2.73554058937173  | 1.63903312374226 | 0.482697050208908  | 0.99278217352652                 | False                  |
| LV_UF_vs_LV_NOT_UF | EBI3      | detected           | 0.0801388278227417  | -1.50275206917344  | 1.66302972481892 | 0.885338534645547  | 0.99278217352652                 | True                   |
| LV_UF_vs_LV_NOT_UF | LTB       | detected           | 0.0067530288120259  | -2.13689130230387  | 2.15039735992792 | 0.99278217352652   | 0.99278217352652                 | True                   |
| LV_UF_vs_LV_NOT_UF | ADAM19    | detected           | -0.302884258130328  | -1.87201018085686  | 1.2662416645962  | 0.584597230260521  | 0.99278217352652                 | False                  |
| LV_UF_vs_LV_NOT_UF | ECM1      | detected           | 5.03357864479517    | 1.2601164937471    | 8.80704079584324 | 0.0026330623909931 | 0.0184314367369517               | True                   |
| LV_UF_vs_LV_NOT_UF | SERPINA10 | not detected       | not estimated       | not estimated      | not estimated    | not estimated      | not estimated                    | not estimated          |
| LV_UF_vs_LV_NOT_UF | CST7      | detected           | -0.0863598976070065 | -1.29211859920621  | 1.1193988039922  | 0.837520370863493  | 0.99278217352652                 | False                  |
| LV_UF_vs_SV        | TNFSF15   | detected           | 2.0312272898759     | -0.388127667057196 | 4.450582246809   | 0.0297242820874419 | 0.0693566582040311               | True                   |
| LV_UF_vs_SV        | FLT3LG    | detected           | 0.746887873720079   | -1.44214612941491  | 2.93592187685507 | 0.330222600033483  | 0.34966029902459                 | True                   |
| LV_UF_vs_SV        | EBI3      | detected           | 1.41350448743709    | -0.175295530632273 | 3.00230450550646 | 0.0195263105393396 | 0.0683420868876886               | True                   |
| LV_UF_vs_SV        | LTB       | detected           | 0.820702974653001   | -1.32299934480378  | 2.96440529410978 | 0.275003121359     | 0.34966029902459                 | True                   |
| LV_UF_vs_SV        | ADAM19    | detected           | 0.874547829037191   | -0.695195250760228 | 2.44429090883461 | 0.119910697164565  | 0.20984372003798873              | True                   |
| LV_UF_vs_SV        | ECM1      | detected           | 7.1359645356895     | 3.35766798995097   | 10.914261081428  | 0.0002874780905162 | 0.0020123466336134               | True                   |
| LV_UF_vs_SV        | SERPINA10 | not detected       | not estimated       | not estimated      | not estimated    | not estimated      | not estimated                    | not estimated          |
| LV_UF_vs_SV        | CST7      | detected           | 0.395369217668576   | -0.810409726588908 | 1.60114816192606 | 0.34966029902459   | 0.34966029902459                 | True                   |

## Expression by cell type

The top cell type is selected by mean normalized expression across all QC-retained cells, without group labels. The full gene-by-cell-type profile is saved in `gene_expression_by_cell_type.csv`.

| gene      | gene_top_cell_type   |   cells |   fraction_cells_nonzero |   mean_log1p_normalized_expression |
|:----------|:---------------------|--------:|-------------------------:|-----------------------------------:|
| TNFSF15   | neutrophil           |    9912 |              0.00696126  |                        0.00870882  |
| FLT3LG    | T/NK                 |   49401 |              0.192668    |                        0.222707    |
| EBI3      | B                    |   12402 |              0.0420094   |                        0.0385921   |
| LTB       | T/NK                 |   49401 |              0.706889    |                        1.50479     |
| ADAM19    | B                    |   12402 |              0.287212    |                        0.257234    |
| ECM1      | fibroblast           |    1278 |              0.676839    |                        1.04015     |
| SERPINA10 | mesothelial          |    6339 |              0.000157754 |                        3.71905e-05 |
| CST7      | T/NK                 |   49401 |              0.657477    |                        1.24704     |

SERPINA10's listed top cell type is a relative maximum, not evidence of meaningful expression: only 0.016% of mesothelial cells had a nonzero count, and it failed the patient-level pseudobulk detection rule.

## QC and analysis notes

Cells were retained using the preregistered group-blind thresholds of at least 200 detected genes and at most 20% mitochondrial counts. The GEO groups were SV n=6, LV_NOT_UF n=6, and LV_UF n=4; all 16 samples passed the 20-cell patient pseudobulk inclusion rule. A panel gene was detected only if raw counts were at least 10 in at least half of patient pseudobulks in each group of the primary contrast. Raw counts were aggregated by sample/patient; cells were not treated as independent replicates. Each GEO sample is used as one anonymous patient pseudobulk because the linked study reports 16 patients and one GEO sample per participant; GEO does not expose participant IDs.

Cell types use the fixed preregistered marker-score rules because cell-level author annotations were not present in the downloaded 10x archive. All resampling and permutation-style tests use seed 42. The composite-AUC random-gene null is used for the supportive decision rule; the preregistered mean-log2FC null is also reported.

edgeR P values are quasi-likelihood F-test results. Approximate 95% intervals use the quasi-likelihood variance-scaled GLM coefficient covariance and residual-adjusted degrees of freedom, expressed on the log2 scale.

## Composition

| cell_type           |   n_LV |   n_SV |   mean_fraction_LV |   mean_fraction_SV |   wilcoxon_p |     bh_q |
|:--------------------|-------:|-------:|-------------------:|-------------------:|-------------:|---------:|
| T/NK                |     10 |      6 |         0.476956   |         0.353567   |    0.427822  | 0.598951 |
| B                   |     10 |      6 |         0.106363   |         0.129352   |    0.367632  | 0.598951 |
| monocyte/macrophage |     10 |      6 |         0.202078   |         0.383591   |    0.0726773 | 0.508741 |
| neutrophil          |     10 |      6 |         0.138968   |         0.0538612  |    0.313187  | 0.598951 |
| mesothelial         |     10 |      6 |         0.0542515  |         0.0708549  |    0.313187  | 0.598951 |
| fibroblast          |     10 |      6 |         0.0165854  |         0.00486141 |    0.792458  | 0.792458 |
| other               |     10 |      6 |         0.00479824 |         0.00391295 |    0.562188  | 0.655886 |

## Composition-adjusted per-gene results

When a gene is not detected under the preregistered raw-count rule, its adjusted estimate and P value are suppressed as not detected. Any detected gene whose adjusted coefficient is no longer positive is described as reflecting immune-cell abundance.

| gene      | gene_id         | detected   |   discovery_log2FC | contrast                      | gene_symbol   | logFC               | ci_lower           | ci_upper          | PValue            | FDR_all_genes     | F                  |   mean_logCPM | detection_status   | direction_concordant   | BH_q_across_tested_panel_genes   |
|:----------|:----------------|:-----------|-------------------:|:------------------------------|:--------------|:--------------------|:-------------------|:------------------|:------------------|:------------------|:-------------------|--------------:|:-------------------|:-----------------------|:---------------------------------|
| TNFSF15   | ENSG00000181634 | True       |           0.97431  | LV_vs_SV_composition_adjusted | TNFSF15       | 0.40228668259488    | -1.57055091884899  | 2.37512428403875  | 0.602971450206611 | 0.908284861304184 | 0.2825402305866    |     -0.934737 | detected           | True                   | 0.8441600302892555               |
| FLT3LG    | ENSG00000090554 | True       |           1.43771  | LV_vs_SV_composition_adjusted | FLT3LG        | 0.262341075806309   | -0.782583975075288 | 1.30726612668791  | 0.520218479674831 | 0.880318991827642 | 0.433870469605829  |      3.27242  | detected           | True                   | 0.8441600302892555               |
| EBI3      | ENSG00000105246 | True       |           0.803423 | LV_vs_SV_composition_adjusted | EBI3          | 1.047814118903      | -0.438075247592933 | 2.53370348539893  | 0.079043588663928 | 0.524565536706056 | 3.55821567379744   |      1.49008  | detected           | True                   | 0.515055752886437                |
| LTB       | ENSG00000227507 | True       |           1.67472  | LV_vs_SV_composition_adjusted | LTB           | -0.0895400828156575 | -1.23525969328444  | 1.05617952765313  | 0.829776268375977 | 0.970972252066781 | 0.0478847380173384 |      7.41802  | detected           | False                  | 0.910849710975844                |
| ADAM19    | ENSG00000135074 | True       |           1.21309  | LV_vs_SV_composition_adjusted | ADAM19        | 0.44648624267357    | -0.566684986592654 | 1.45965747193979  | 0.223408939129333 | 0.711506371034306 | 1.61503112078911   |      4.40965  | detected           | True                   | 0.5212875246351103               |
| ECM1      | ENSG00000143369 | True       |           1.03405  | LV_vs_SV_composition_adjusted | ECM1          | 2.9878617115625     | -1.5641779362542   | 7.53990135937921  | 0.147158786538982 | 0.636427508032712 | 2.34040069338931   |      4.00055  | detected           | True                   | 0.515055752886437                |
| SERPINA10 | ENSG00000140093 | False      |           1.16585  | LV_vs_SV_composition_adjusted | SERPINA10     | not estimated       | not estimated      | not estimated     | not estimated     | not estimated     | not estimated      |     -5.07091  | not detected       | not estimated          | not estimated                    |
| CST7      | ENSG00000077984 | True       |           1.01227  | LV_vs_SV_composition_adjusted | CST7          | 0.0279949807242618  | -0.670149853277969 | 0.726139814726493 | 0.910849710975844 | 0.986850046772784 | 0.0129728874620982 |      7.26322  | detected           | True                   | 0.910849710975844                |

LTB was positive before adjustment but negative after including monocyte/macrophage and T/NK proportions; under the preregistered rule its unadjusted signal is described as reflecting immune-cell abundance. No other detected panel gene reversed direction.

## Top-cell-type per-gene results

| gene      | gene_id         | top_cell_type   | status                        |   n_LV |   n_SV |   cells_min_count_LV |   cells_min_count_SV |       logFC |    ci_lower |   ci_upper |       PValue |   BH_q_across_panel_genes |
|:----------|:----------------|:----------------|:------------------------------|-------:|-------:|---------------------:|---------------------:|------------:|------------:|-----------:|-------------:|--------------------------:|
| TNFSF15   | ENSG00000181634 | neutrophil      | not detected in top cell type |     10 |      6 |                    3 |                    0 | nan         | nan         | nan        | nan          |              nan          |
| FLT3LG    | ENSG00000090554 | T/NK            | tested                        |     10 |      6 |                   10 |                    6 |   0.0215937 |  -0.721088  |   0.764275 |   0.934351   |                0.934351   |
| EBI3      | ENSG00000105246 | B               | tested                        |     10 |      6 |                   10 |                    6 |   1.98189   |   0.85307   |   3.11071  |   0.00070729 |                0.00424374 |
| LTB       | ENSG00000227507 | T/NK            | tested                        |     10 |      6 |                   10 |                    6 |  -0.293212  |  -0.831214  |   0.244789 |   0.123612   |                0.185418   |
| ADAM19    | ENSG00000135074 | B               | tested                        |     10 |      6 |                   10 |                    6 |   1.03158   |  -0.0172247 |   2.08038  |   0.0282535  |                0.0565069  |
| ECM1      | ENSG00000143369 | fibroblast      | tested                        |      4 |      5 |                    4 |                    4 |   0.575138  |  -0.95528   |   2.10556  |   0.389183   |                0.467019   |
| SERPINA10 | ENSG00000140093 | mesothelial     | not detected in top cell type |     10 |      6 |                    0 |                    0 | nan         | nan         | nan        | nan          |              nan          |
| CST7      | ENSG00000077984 | T/NK            | tested                        |     10 |      6 |                   10 |                    6 |  -0.418004  |  -0.936467  |   0.100459 |   0.0276628  |                0.0565069  |

## Software versions

- Python 3.13.1
- Scanpy 1.12.4; AnnData 0.13.4.
- R and edgeR versions are recorded in `edger_analysis_config.csv`.

This analysis is limited to the GSE248762 cohort and does not establish independence from GSE130888; the preregistered pairwise assessment remains unknown. GSE248762 is from Guangzhou, China, whereas discovery GSE125498 is associated with IKEM in Prague; their patient overlap is unlikely but not formally confirmed. The cohort-pair assessment is in `../independence_assessment.md`. No expression analysis of GSE130888 or GSE92455 was performed. The small patient counts and wide intervals make the estimates uncertain, not evidence that an effect is absent; models did not include covariates beyond the preregistered cell-type proportions.
