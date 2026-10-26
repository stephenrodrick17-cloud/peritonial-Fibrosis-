# GSE125498 panel context analysis

**Label:** exploratory, post hoc, same data as discovery, not validation.

All analyses in this report reuse the discovery cohort that informed gene selection. They are descriptive/exploratory and are not validation.

## Cohort and scoring

- Samples: 33 total; 13 LPD and 20 SPD.
- Immune cell-type scores: GSEApy ssGSEA using human CellMarker_2024 signatures retrieved from the Enrichr library. For each broad cell type, matching human tissue/context marker sets were unioned before scoring; signatures are expression enrichment estimates, not cell counts or validated deconvolution.
- Dominant leukocyte covariate: mean per-gene z-score of the 17 repository pan-leukocyte markers (17 observed in the matrix: PTPRC, CD3D, CD3E, CD3G, CD4, CD8A, CD8B, CD19, MS4A1, CD79A, CD79B, NCAM1, KLRD1, CD14, CD68, ITGAM, FCGR3A); this is a proxy, not an independently measured cell fraction.
- Composite scores: average of gene-wise z-scores across these 33 samples. Original panel has eight genes; panel v2 has five.
- Group differences use two-sided Mann–Whitney tests; BH q-values are adjusted across the 13 immune cell-type comparisons.

## Immune cell-type differences

| Cell type | LPD mean ssGSEA | SPD mean ssGSEA | LPD − SPD | Mann–Whitney P | BH q |
|---|---:|---:|---:|---:|---:|
| NK cell | 0.3372 | 0.2401 | +0.0971 | 0.002675 | 0.01739 |
| B cell | 0.3157 | 0.2870 | +0.0287 | 0.001438 | 0.01739 |
| T cell (other) | 0.4250 | 0.3580 | +0.0670 | 0.005407 | 0.01757 |
| CD8 T cell | 0.4422 | 0.3938 | +0.0484 | 0.005407 | 0.01757 |
| NKT cell | 0.2638 | 0.0665 | +0.1973 | 0.006769 | 0.0176 |
| CD4 T cell | 0.4560 | 0.4095 | +0.0465 | 0.01288 | 0.02792 |
| Dendritic cell | 0.4700 | 0.4539 | +0.0161 | 0.03737 | 0.06941 |
| Plasma cell | 0.6137 | 0.5926 | +0.0211 | 0.04871 | 0.07915 |
| Neutrophil | 0.4780 | 0.4410 | +0.0370 | 0.2037 | 0.2942 |
| Macrophage | 0.6966 | 0.6805 | +0.0161 | 0.311 | 0.4043 |
| Eosinophil | 0.4566 | 0.4326 | +0.0240 | 0.3475 | 0.4106 |
| Monocyte | 0.7414 | 0.7435 | -0.0021 | 0.8974 | 0.9722 |
| Mast cell | 0.1948 | 0.1947 | +0.0001 | 0.9853 | 0.9853 |

## Correlation and leukocyte adjustment

The table gives Pearson correlations for each v2 gene and the v2 composite against each immune cell-type ssGSEA score. Spearman results and original-8 correlations are in `panel_celltype_correlations.csv`; Pearson BH q is adjusted within each score. Group-adjusted model details are in `composite_leukocyte_adjustment.csv`.

| V2 score | Cell type | Pearson r | P | BH q within score |
|---|---|---:|---:|---:|
| panel_v2_FLT3LG | CD4 T cell | 0.7540 | 4.054e-07 | 1.317e-06 |
| panel_v2_FLT3LG | CD8 T cell | 0.7852 | 6.263e-08 | 2.714e-07 |
| panel_v2_FLT3LG | NKT cell | 0.8541 | 2.596e-10 | 3.375e-09 |
| panel_v2_FLT3LG | NK cell | 0.8116 | 9.98e-09 | 6.487e-08 |
| panel_v2_FLT3LG | T cell (other) | 0.7319 | 1.299e-06 | 3.378e-06 |
| panel_v2_FLT3LG | B cell | 0.3860 | 0.0265 | 0.04922 |
| panel_v2_FLT3LG | Plasma cell | 0.3988 | 0.02152 | 0.04662 |
| panel_v2_FLT3LG | Monocyte | 0.0787 | 0.6632 | 0.7185 |
| panel_v2_FLT3LG | Macrophage | 0.2480 | 0.164 | 0.2132 |
| panel_v2_FLT3LG | Dendritic cell | 0.3394 | 0.05335 | 0.0867 |
| panel_v2_FLT3LG | Neutrophil | 0.3070 | 0.08229 | 0.1189 |
| panel_v2_FLT3LG | Eosinophil | 0.2253 | 0.2075 | 0.2452 |
| panel_v2_FLT3LG | Mast cell | -0.0552 | 0.7602 | 0.7602 |
| panel_v2_EBI3 | CD4 T cell | 0.3588 | 0.04033 | 0.0937 |
| panel_v2_EBI3 | CD8 T cell | 0.3914 | 0.02428 | 0.0937 |
| panel_v2_EBI3 | NKT cell | 0.3595 | 0.03987 | 0.0937 |
| panel_v2_EBI3 | NK cell | 0.4040 | 0.01973 | 0.0937 |
| panel_v2_EBI3 | T cell (other) | 0.3540 | 0.04325 | 0.0937 |
| panel_v2_EBI3 | B cell | 0.2967 | 0.09355 | 0.1737 |
| panel_v2_EBI3 | Plasma cell | 0.1905 | 0.2882 | 0.3482 |
| panel_v2_EBI3 | Monocyte | 0.1781 | 0.3214 | 0.3482 |
| panel_v2_EBI3 | Macrophage | 0.2575 | 0.1479 | 0.2403 |
| panel_v2_EBI3 | Dendritic cell | 0.4399 | 0.01041 | 0.0937 |
| panel_v2_EBI3 | Neutrophil | 0.2145 | 0.2306 | 0.333 |
| panel_v2_EBI3 | Eosinophil | 0.1811 | 0.3131 | 0.3482 |
| panel_v2_EBI3 | Mast cell | 0.1671 | 0.3528 | 0.3528 |
| panel_v2_LTB | CD4 T cell | 0.6389 | 6.29e-05 | 0.0001635 |
| panel_v2_LTB | CD8 T cell | 0.5750 | 0.0004648 | 0.000964 |
| panel_v2_LTB | NKT cell | 0.3955 | 0.0227 | 0.02952 |
| panel_v2_LTB | NK cell | 0.5497 | 0.0009204 | 0.001496 |
| panel_v2_LTB | T cell (other) | 0.6974 | 6.477e-06 | 2.807e-05 |
| panel_v2_LTB | B cell | 0.7732 | 1.33e-07 | 1.251e-06 |
| panel_v2_LTB | Plasma cell | 0.7670 | 1.925e-07 | 1.251e-06 |
| panel_v2_LTB | Monocyte | 0.1924 | 0.2833 | 0.2833 |
| panel_v2_LTB | Macrophage | 0.3737 | 0.03215 | 0.038 |
| panel_v2_LTB | Dendritic cell | 0.6468 | 4.761e-05 | 0.0001547 |
| panel_v2_LTB | Neutrophil | 0.4037 | 0.01983 | 0.02864 |
| panel_v2_LTB | Eosinophil | 0.5711 | 0.0005191 | 0.000964 |
| panel_v2_LTB | Mast cell | 0.2565 | 0.1495 | 0.162 |
| panel_v2_ADAM19 | CD4 T cell | 0.5543 | 0.0008161 | 0.001326 |
| panel_v2_ADAM19 | CD8 T cell | 0.6088 | 0.0001703 | 0.0005616 |
| panel_v2_ADAM19 | NKT cell | 0.6083 | 0.0001728 | 0.0005616 |
| panel_v2_ADAM19 | NK cell | 0.6779 | 1.463e-05 | 0.0001902 |
| panel_v2_ADAM19 | T cell (other) | 0.6437 | 5.323e-05 | 0.000346 |
| panel_v2_ADAM19 | B cell | 0.5610 | 0.0006824 | 0.001267 |
| panel_v2_ADAM19 | Plasma cell | 0.5758 | 0.000455 | 0.001183 |
| panel_v2_ADAM19 | Monocyte | 0.1183 | 0.512 | 0.512 |
| panel_v2_ADAM19 | Macrophage | 0.3210 | 0.06852 | 0.08098 |
| panel_v2_ADAM19 | Dendritic cell | 0.4945 | 0.00344 | 0.004968 |
| panel_v2_ADAM19 | Neutrophil | 0.3278 | 0.06254 | 0.08098 |
| panel_v2_ADAM19 | Eosinophil | 0.5669 | 0.0005813 | 0.00126 |
| panel_v2_ADAM19 | Mast cell | 0.1797 | 0.3171 | 0.3435 |
| panel_v2_CST7 | CD4 T cell | 0.7283 | 1.551e-06 | 4.032e-06 |
| panel_v2_CST7 | CD8 T cell | 0.7650 | 2.164e-07 | 9.378e-07 |
| panel_v2_CST7 | NKT cell | 0.8590 | 1.591e-10 | 1.034e-09 |
| panel_v2_CST7 | NK cell | 0.8683 | 5.889e-11 | 7.655e-10 |
| panel_v2_CST7 | T cell (other) | 0.7549 | 3.859e-07 | 1.254e-06 |
| panel_v2_CST7 | B cell | 0.4584 | 0.007296 | 0.01355 |
| panel_v2_CST7 | Plasma cell | 0.4955 | 0.003366 | 0.007292 |
| panel_v2_CST7 | Monocyte | -0.0192 | 0.9154 | 0.9154 |
| panel_v2_CST7 | Macrophage | 0.2038 | 0.2554 | 0.332 |
| panel_v2_CST7 | Dendritic cell | 0.3102 | 0.07893 | 0.1283 |
| panel_v2_CST7 | Neutrophil | 0.1829 | 0.3082 | 0.3642 |
| panel_v2_CST7 | Eosinophil | 0.2365 | 0.1851 | 0.2674 |
| panel_v2_CST7 | Mast cell | -0.0407 | 0.822 | 0.8905 |
| panel_v2_composite | CD4 T cell | 0.7623 | 2.542e-07 | 6.609e-07 |
| panel_v2_composite | CD8 T cell | 0.7852 | 6.286e-08 | 2.724e-07 |
| panel_v2_composite | NKT cell | 0.7729 | 1.358e-07 | 4.415e-07 |
| panel_v2_composite | NK cell | 0.8319 | 1.987e-09 | 2.583e-08 |
| panel_v2_composite | T cell (other) | 0.7993 | 2.43e-08 | 1.58e-07 |
| panel_v2_composite | B cell | 0.6219 | 0.0001119 | 0.0002425 |
| panel_v2_composite | Plasma cell | 0.6099 | 0.0001647 | 0.0003058 |
| panel_v2_composite | Monocyte | 0.1377 | 0.4446 | 0.4797 |
| panel_v2_composite | Macrophage | 0.3527 | 0.04408 | 0.05209 |
| panel_v2_composite | Dendritic cell | 0.5604 | 0.0006942 | 0.001128 |
| panel_v2_composite | Neutrophil | 0.3607 | 0.03918 | 0.05093 |
| panel_v2_composite | Eosinophil | 0.4474 | 0.00904 | 0.01306 |
| panel_v2_composite | Mast cell | 0.1274 | 0.4797 | 0.4797 |

- **panel_v2_composite:** unadjusted group coefficient 1.1072 (P=8.87e-06); after adding the dominant leukocyte marker proxy, group coefficient 0.5392 (P=0.003078), proxy coefficient 1.1433 (P=7.829e-07), and group-coefficient attenuation 51.3%.
- **original_8_composite:** unadjusted group coefficient 1.1026 (P=4.289e-07); after adding the dominant leukocyte marker proxy, group coefficient 0.6395 (P=9.422e-05), proxy coefficient 0.9323 (P=1.514e-06), and group-coefficient attenuation 42.0%.

**Interpretation:** The panel-v2 group coefficient is attenuated by 51.3% after adding the leukocyte-marker proxy, but remains nominally significant (P=0.003078). This is consistent with a substantial leukocyte-composition contribution, but the proxy does not fully explain the association; neither the attenuation nor residual association establishes causality or cell-intrinsic expression.

## Preranked Hallmark GSEA

The complete 11,741-gene limma t-statistic ranking was analyzed with GSEApy prerank (1,000 permutations, seed 125498) using `MSigDB_Hallmark_2020`, the same library source used elsewhere in this repository. The table below includes BH adjustment of nominal pathway P values across the 50 Hallmark sets; GSEApy's FDR q-values are retained separately in the CSV.

| Pathway | NES | Nominal P | BH q | GSEA FDR q | Panel-v2 genes in leading edge |
|---|---:|---:|---:|---:|---|
| Myc Targets V1 | -2.2084 | <0.001 | <0.001 | <0.001 | — |
| IL-2/STAT5 Signaling | 2.1427 | <0.001 | <0.001 | <0.001 | FLT3LG;LTB;ADAM19;CST7 |
| Interferon Alpha Response | 2.1338 | <0.001 | <0.001 | <0.001 | — |
| Allograft Rejection | 2.0376 | <0.001 | <0.001 | <0.001 | LTB |
| Interferon Gamma Response | 1.9789 | <0.001 | <0.001 | <0.001 | — |
| Oxidative Phosphorylation | -1.8333 | <0.001 | <0.001 | 0.01103 | — |
| Protein Secretion | -1.8287 | <0.001 | <0.001 | 0.00776 | — |
| Adipogenesis | -1.7730 | <0.001 | <0.001 | 0.006371 | — |
| TNF-alpha Signaling via NF-kB | 1.7280 | <0.001 | <0.001 | 0.008554 | — |
| mTORC1 Signaling | -1.6572 | <0.001 | <0.001 | 0.01327 | — |
| E2F Targets | -1.5824 | <0.001 | <0.001 | 0.0224 | — |
| Myc Targets V2 | -1.7975 | 0.002203 | 0.009178 | 0.006739 | — |
| Fatty Acid Metabolism | -1.4865 | 0.002445 | 0.009404 | 0.04717 | — |
| KRAS Signaling Up | 1.5833 | 0.003442 | 0.01229 | 0.02618 | — |
| p53 Pathway | 1.5078 | 0.006734 | 0.02196 | 0.03794 | — |

All 50 sets, full leading-edge strings, and panel-gene membership are in `hallmark_preranked_gsea.csv`; the retrieved gene sets are saved in `hallmark_gene_sets.json`. Across all 50 pathways, 18 have BH q < 0.05. The panel-v2 leading-edge memberships are: IL-2/STAT5 Signaling (FLT3LG, LTB, ADAM19, CST7; BH q < 0.001), Allograft Rejection (LTB; BH q < 0.001), Inflammatory Response (EBI3; BH q=0.0251), and IL-6/JAK/STAT3 Signaling (EBI3, LTB; BH q=0.3712).

Nominal permutation P values of 0 are displayed as <0.001; with 1,000 permutations they are below the reporting resolution, not literal zero.

## Descriptive ROC

**circular: genes were selected on this data.** ROC estimates below are in-sample descriptions, not performance validation. Bootstrap CIs resample within LPD/SPD groups (5,000 replicates; seed 125498).

| Score | AUC | Bootstrap 95% CI | Mann–Whitney P |
|---|---:|---:|---:|
| panel_v2_FLT3LG | 0.9038 | 0.7845–0.9885 | 0.000118 |
| panel_v2_EBI3 | 0.7923 | 0.6192–0.9346 | 0.005407 |
| panel_v2_LTB | 0.8038 | 0.6231–0.9615 | 0.003825 |
| panel_v2_ADAM19 | 0.8769 | 0.7346–0.9846 | 0.0003278 |
| panel_v2_CST7 | 0.8308 | 0.6808–0.9538 | 0.001632 |
| panel_v2_composite | 0.9308 | 0.8346–0.9923 | 3.99e-05 |
| original_8_TNFSF15 | 0.9038 | 0.7885–0.9846 | 0.0001179 |
| original_8_FLT3LG | 0.9038 | 0.7885–0.9885 | 0.000118 |
| original_8_EBI3 | 0.7923 | 0.6231–0.9346 | 0.005407 |
| original_8_LTB | 0.8038 | 0.6154–0.9615 | 0.003825 |
| original_8_ADAM19 | 0.8769 | 0.7308–0.9846 | 0.0003278 |
| original_8_ECM1 | 0.7654 | 0.5731–0.9231 | 0.01161 |
| original_8_SERPINA10 | 0.8231 | 0.6615–0.9538 | 0.002095 |
| original_8_CST7 | 0.8308 | 0.6731–0.9500 | 0.001632 |
| original_8_composite | 0.9500 | 0.8654–1.0000 | 1.768e-05 |

Figure: `discovery_circular_roc.png`.

## DGIdb interactions (retrieved live)

Queried DGIdb GraphQL on 2026-10-11. No drug ranking was performed. Empty interaction rows mean the endpoint returned no interaction for that queried gene at retrieval time; this is not proof that no interaction is known.

| Gene | Drug | Interaction type | Directionality | Source(s) | Publication(s) |
|---|---|---|---|---|---|
| FLT3LG | — | No interaction returned | — | — | — |
| EBI3 | — | No interaction returned | — | — | — |
| LTB | BAMINERCEPT | inhibitor | INHIBITORY | ChEMBL | not reported |
| LTB | ABACAVIR | not reported | not reported | PharmGKB | 16998491 |
| LTB | BAMINERCEPT ALFA | not reported | not reported | TdgClinicalTrial | not reported |
| ADAM19 | — | No interaction returned | — | — | — |
| CST7 | — | No interaction returned | — | — | — |

The parsed result is in `dgidb_panel_v2_interactions.csv`; the raw response, query, and retrieval timestamp are in `dgidb_raw_response.json`.

## Caveats

- CellMarker tissue-context markers were unioned to obtain broad immune-cell signatures; tissue effects and overlapping marker biology may influence scores.
- The cohort is small (13 LPD, 20 SPD); nominal associations are unstable and exploratory.
- ROC and discovery GSEA are computed on the same data used in feature selection.
- GSEApy reported a tie affecting 0.01% of ranked genes; the tied order may be arbitrary and could affect results at that tie boundary.
- DGIdb is a database snapshot queried live; database coverage and annotation can change.
