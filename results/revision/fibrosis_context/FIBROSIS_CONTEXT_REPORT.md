# GSE248762 fibrosis-context analysis

## Scope and status

This report covers GSE248762 only. It uses the fixed patient/sample inclusion and Scanpy QC already applied in the effluent analysis. Patients are the independent units; the all-cell pseudobulk analysis has n=16, so uncertainty is wide. The frozen support rule is in [PREREGISTRATION_FIBROSIS.md](PREREGISTRATION_FIBROSIS.md).

All findings are **correlative, not causal**. Cell-cell communication results are **inferred**, not measured. The MMT quantity is explicitly a study-specific marker-panel proxy, not a validated canonical MMT score.

## Frozen success rule

Support is assessed only from all-cell pseudobulk composition-adjusted correlations and requires at least three markers with rho > 0.4 plus positive LV_UF-versus-LV_NOT_UF edgeR log2FC.

| gene      |   n_markers_with_adjusted_rho_gt_0_4 | correlation_condition_met   |   LV_UF_vs_LV_NOT_UF_log2FC |   log2FC_CI_lower |   log2FC_CI_upper |   edgeR_P |   edgeR_FDR | contrast_concordant_positive   | supported   |
|:----------|-------------------------------------:|:----------------------------|----------------------------:|------------------:|------------------:|----------:|------------:|:-------------------------------|:------------|
| TNFSF15   |                                    1 | False                       |                    1.477    |           -0.9283 |             3.883 |  0.0992   |     0.4407  | True                           | False       |
| FLT3LG    |                                    2 | False                       |                   -0.5483   |           -2.736  |             1.639 |  0.4827   |     0.7904  | False                          | False       |
| EBI3      |                                    2 | False                       |                    0.08014  |           -1.503  |             1.663 |  0.8853   |     0.966   | True                           | False       |
| LTB       |                                    2 | False                       |                    0.006753 |           -2.137  |             2.15  |  0.9928   |     0.9979  | True                           | False       |
| ADAM19    |                                    2 | False                       |                   -0.3029   |           -1.872  |             1.266 |  0.5846   |     0.8445  | False                          | False       |
| ECM1      |                                    0 | False                       |                    5.034    |            1.26   |             8.807 |  0.002633 |     0.06749 | True                           | False       |
| SERPINA10 |                                    0 | False                       |                    2.894    |          nan      |           nan     |  0.07373  |     0.3941  | True                           | False       |
| CST7      |                                    2 | False                       |                   -0.08636  |           -1.292  |             1.119 |  0.8375   |     0.949   | False                          | False       |

No changes were made to the candidate genes, markers, covariates, contrast, thresholds, or success rule after analysis.

## Patient-level gene-marker correlations

The complete all-cell and per-major-cell-type results (n, rho, patient-bootstrap 95% CI, two-sided P, and BH q) are in [gene_marker_correlations.csv](gene_marker_correlations.csv). The patient-level normalized expression used is in [patient_gene_marker_expression.csv](patient_gene_marker_expression.csv). Composition adjustment residualizes rank-transformed expression on rank-transformed monocyte/macrophage and T/NK fractions. Bootstrap intervals use 2,000 patient resamples (seed 42).

All-cell composition-adjusted correlations with rho > 0.4:

| gene    | marker   |   n |    rho |   ci_lower |   ci_upper |        p |       q |
|:--------|:---------|----:|-------:|-----------:|-----------:|---------:|--------:|
| TNFSF15 | VEGFA    |  16 | 0.6901 |    0.1145  |     0.9483 | 0.006306 | 0.08072 |
| FLT3LG  | TGFB1    |  16 | 0.7759 |    0.3209  |     0.9527 | 0.001106 | 0.0354  |
| FLT3LG  | VEGFA    |  16 | 0.4812 |   -0.2947  |     0.8199 | 0.08151  | 0.3478  |
| EBI3    | TGFB1    |  16 | 0.5625 |   -0.04139 |     0.8767 | 0.03627  | 0.211   |
| EBI3    | VEGFA    |  16 | 0.5337 |   -0.2341  |     0.8343 | 0.04936  | 0.2632  |
| LTB     | TGFB1    |  16 | 0.4987 |   -0.1427  |     0.9052 | 0.0695   | 0.3177  |
| LTB     | VEGFA    |  16 | 0.5931 |   -0.07602 |     0.8705 | 0.02539  | 0.1788  |
| ADAM19  | TGFB1    |  16 | 0.6435 |    0.1437  |     0.8866 | 0.01303  | 0.139   |
| ADAM19  | VEGFA    |  16 | 0.5851 |    0.08813 |     0.8519 | 0.02794  | 0.1788  |
| CST7    | TGFB1    |  16 | 0.7572 |    0.2152  |     0.9346 | 0.001714 | 0.03656 |
| CST7    | VEGFA    |  16 | 0.5055 |   -0.2735  |     0.8076 | 0.06519  | 0.3177  |

## LV_UF versus LV_NOT_UF edgeR contrast

The existing edgeR pseudobulk pipeline was rerun in a temporary copy of its input files; its earlier outputs were not overwritten. The two LV groups were compared without SV samples. Package/model details: R_version=R version 4.4.2 (2024-10-31 ucrt); edgeR_version=4.4.2; composition_adjustment=feasible: both groups have at least five patients; design full rank; primary_patients_LV=10; primary_patients_SV=6

| gene_symbol   |     logFC |   ci_lower |   ci_upper |   PValue |   FDR_all_genes |
|:--------------|----------:|-----------:|-----------:|---------:|----------------:|
| ECM1          |  5.034    |     1.26   |      8.807 | 0.002633 |         0.06749 |
| ADAM19        | -0.3029   |    -1.872  |      1.266 | 0.5846   |         0.8445  |
| LTB           |  0.006753 |    -2.137  |      2.15  | 0.9928   |         0.9979  |
| TNFSF15       |  1.477    |    -0.9283 |      3.883 | 0.0992   |         0.4407  |
| SERPINA10     |  2.894    |   nan      |    nan     | 0.07373  |         0.3941  |
| CST7          | -0.08636  |    -1.292  |      1.119 | 0.8375   |         0.949   |
| EBI3          |  0.08014  |    -1.503  |      1.663 | 0.8853   |         0.966   |
| FLT3LG        | -0.5483   |    -2.736  |      1.639 | 0.4827   |         0.7904  |

## Compartment-level preranked GSEA

Gene rankings use edgeR robust quasi-likelihood statistics for LV_UF relative to LV_NOT_UF, separately for monocyte/macrophage, T/NK, fibroblast, and mesothelial pseudobulks. Gene sets are the Hallmark 2020 EMT, TGF-beta, and inflammatory-response sets and GO BP 2023 Extracellular Matrix Organization (GO:0030198). Full NES, P, FDR, and panel leading-edge membership are in [compartment_gsea_results.csv](compartment_gsea_results.csv).

Gene-set P/FDR estimates with zero permutation exceedances are reported at their empirical resolution floor (1/1,001), not as exact zero.

| compartment         | gene_set                                   |    NES |        P |       FDR |   leading_edge_size | panel_genes_in_leading_edge   |
|:--------------------|:-------------------------------------------|-------:|---------:|----------:|--------------------:|:------------------------------|
| monocyte/macrophage | HALLMARK_INFLAMMATORY_RESPONSE             |  1.882 | 0.000999 | 0.000999  |                  43 | TNFSF15                       |
| monocyte/macrophage | HALLMARK_TGF_BETA_SIGNALING                |  1.638 | 0.005128 | 0.00452   |                  20 | nan                           |
| monocyte/macrophage | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION |  1.451 | 0.0191   | 0.02195   |                  32 | ECM1                          |
| monocyte/macrophage | GO_EXTRACELLULAR_MATRIX_ORGANIZATION       | -1.402 | 0.04054  | 0.02922   |                  17 | nan                           |
| T/NK                | GO_EXTRACELLULAR_MATRIX_ORGANIZATION       | -2.209 | 0.000999 | 0.000999  |                  34 | nan                           |
| T/NK                | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | -1.986 | 0.000999 | 0.000999  |                  52 | nan                           |
| T/NK                | HALLMARK_INFLAMMATORY_RESPONSE             | -1.829 | 0.001458 | 0.0005156 |                  42 | nan                           |
| T/NK                | HALLMARK_TGF_BETA_SIGNALING                | -1.405 | 0.03865  | 0.029     |                  12 | nan                           |
| fibroblast          | GO_EXTRACELLULAR_MATRIX_ORGANIZATION       | -2.168 | 0.000999 | 0.000999  |                  39 | nan                           |
| fibroblast          | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | -2.117 | 0.000999 | 0.000999  |                  56 | ECM1                          |
| fibroblast          | HALLMARK_INFLAMMATORY_RESPONSE             |  1.907 | 0.000999 | 0.000999  |                  25 | nan                           |
| fibroblast          | HALLMARK_TGF_BETA_SIGNALING                |  1.198 | 0.1772   | 0.1411    |                  11 | nan                           |
| mesothelial         | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | -2.289 | 0.000999 | 0.000999  |                  88 | nan                           |
| mesothelial         | GO_EXTRACELLULAR_MATRIX_ORGANIZATION       | -2.215 | 0.000999 | 0.000999  |                  60 | nan                           |
| mesothelial         | HALLMARK_INFLAMMATORY_RESPONSE             |  1.186 | 0.1535   | 0.3538    |                  36 | nan                           |
| mesothelial         | HALLMARK_TGF_BETA_SIGNALING                |  1.183 | 0.2272   | 0.1787    |                  16 | nan                           |

## Inferred cell-cell communication

LIANA 1.10.0 consensus resource analyzed 20,000 observed cells subsampled with seed 42 (target at most 20,000). Predictions connecting immune and stromal cell types are in [liana_inferred_targeted_interactions.csv](liana_inferred_targeted_interactions.csv).

Predicted ligand-receptor scores are not evidence of signaling activity. The consensus resource includes TNFSF15-TNFRSF25 and EBI3/IL-27 receptor pairs. It does not contain a direct ADAM19-Notch interaction; JAG/DLL-Notch rankings are separately inferred and do not attribute receptor processing to ADAM19.

The route-level summary below aggregates the rank of the best resource-supported prediction within each direction; lower rank is better, and these rankings are not hypothesis-test P values.

| route_class      | direction         |   resource_pairs_evaluated |   best_cellphone_p |   best_magnitude_rank |   best_specificity_rank | best_source   | best_target         | best_ligand   | best_receptor   |
|:-----------------|:------------------|---------------------------:|-------------------:|----------------------:|------------------------:|:--------------|:--------------------|:--------------|:----------------|
| IL-27            | immune-to-stromal |                          8 |                  1 |             1         |               1         | B             | fibroblast          | EBI3_IL27     | IL27RA_IL6ST    |
| IL-27            | stromal-to-immune |                          8 |                  1 |             1         |               1         | fibroblast    | B                   | EBI3_IL27     | IL27RA_IL6ST    |
| JAG/DLL-Notch    | immune-to-stromal |                        160 |                  1 |             1         |               1         | B             | fibroblast          | JAG2          | NOTCH3          |
| JAG/DLL-Notch    | stromal-to-immune |                        160 |                  0 |             0.0001183 |               6.354e-07 | fibroblast    | monocyte/macrophage | JAG1          | NOTCH2          |
| TNFSF15-TNFRSF25 | immune-to-stromal |                          8 |                  1 |             1         |               1         | B             | fibroblast          | TNFSF15       | TNFRSF25        |
| TNFSF15-TNFRSF25 | stromal-to-immune |                          8 |                  1 |             1         |               1         | fibroblast    | B                   | TNFSF15       | TNFRSF25        |

No immune-to-stromal TNFSF15-TNFRSF25 or IL-27 pair had an informative LIANA score in this analysis: all such CellPhoneDB p-values and magnitude/specificity ranks were 1.0. Any reverse-direction or Notch ranking remains inferred and does not establish signaling.

Full targeted output has 352 immune/stromal directed predictions, including non-significant resource-supported pairs.

All LIANA results are inferred ligand-receptor predictions, not measured signaling. ADAM19 was not a direct interaction in the selected consensus resource; Notch-route predictions do not establish ADAM19-mediated Notch signaling.

## MMT proxy by PD group

The patient-level score is the mean of the eight marker-gene z-scores within each cell type after Scanpy total-count normalization to 10,000 and log1p transformation. Patient is the unit of comparison; group summaries and Kruskal-Wallis P values are in [mmt_proxy_group_summary.csv](mmt_proxy_group_summary.csv), and individual patient scores are in [mmt_proxy_patient_scores.csv](mmt_proxy_patient_scores.csv).

| cell_type   | group     |   n_patients |   median |        q1 |      q3 |   omnibus_Kruskal_P |
|:------------|:----------|-------------:|---------:|----------:|--------:|--------------------:|
| fibroblast  | LV_UF     |            4 | -0.6248  | -0.6543   | -0.54   |              0.1353 |
| fibroblast  | LV_NOT_UF |            6 | -0.0738  | -0.4301   |  0.1645 |              0.1353 |
| fibroblast  | SV        |            6 | -0.3717  | -0.5769   | -0.2444 |              0.1353 |
| mesothelial | LV_UF     |            4 | -0.2671  | -0.3716   | -0.0912 |              0.1353 |
| mesothelial | LV_NOT_UF |            6 |  0.01425 | -0.02097  |  0.1676 |              0.1353 |
| mesothelial | SV        |            6 |  0.1287  |  0.005915 |  0.2053 |              0.1353 |

## Interpretation

No candidate gene met the frozen support rule. With only 16 independent
patients overall (and 4 versus 6 in the LV_UF versus LV_NOT_UF contrast), power
is low and uncertainty is substantial. These patient-level associations,
expression contrasts, pathway enrichments, inferred communications, and the
study-specific MMT proxy are correlative and inconclusive; they do not establish
causation or validate a mechanism.

## Reproducibility files

The analysis driver is [run_fibrosis_context.py](run_fibrosis_context.py); compartment-level edgeR ranking code is [compartment_edger_gsea_input.R](compartment_edger_gsea_input.R). The pseudobulk input matrices and metadata used for this analysis are included in this folder as [all_cell_pseudobulk_counts.csv](all_cell_pseudobulk_counts.csv), [celltype_pseudobulk_counts.csv](celltype_pseudobulk_counts.csv), [all_cell_pseudobulk_metadata.csv](all_cell_pseudobulk_metadata.csv), [celltype_pseudobulk_metadata.csv](celltype_pseudobulk_metadata.csv), and [patient_celltype_composition.csv](patient_celltype_composition.csv). Exact GSEA gene-set memberships are in [gsea_gene_set_membership.csv](gsea_gene_set_membership.csv). The original edgeR pipeline and source data remain in their existing locations and were not edited by this analysis.
