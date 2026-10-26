# Figure 7: exploratory LPD subgroups

**Exploratory; n=13 LPD samples; not validated.**

## Gene set and samples

Current candidate panel: FLT3LG, EBI3, LTB, ADAM19, CST7 (panel-v2 rule-derived set).
Each gene was z-scored across these 13 samples using the sample SD (ddof=1).

GSM3575609, GSM3575613, GSM3575614, GSM3575615, GSM3575617, GSM3575620, GSM3575621, GSM3575622, GSM3575623, GSM3575624, GSM3575629, GSM3575635, GSM3575639

## Consensus clustering

Resampling: 1000 subsamples, pItem=0.8, seed=20261016; each subsample uses Ward agglomerative clustering with Euclidean distance. For each k, consensus is the fraction of co-sampled resamples in which a sample pair co-clustered. Final consensus partitions use complete linkage on distance 1-consensus. k was selected by the largest Euclidean silhouette on the full z-scored expression matrix (ties to smaller k), not by visual inspection of the figure.

CDF curves for all k and the selected-k consensus matrix are exported to `consensus_cdf.csv` and `consensus_matrix_selected_k.csv`.

| k | CDF area | Delta area | Euclidean silhouette | Consensus-distance silhouette | PAC (0.1, 0.9) | Cluster sizes |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 0.2539 | NA | 0.3930 | 0.7894 | 0.4744 | 1:12;2:1 |
| 3 | 0.4785 | 0.2246 | 0.3216 | 0.5706 | 0.6154 | 1:11;2:1;3:1 |
| 4 | 0.6741 | 0.1956 | 0.2031 | 0.4835 | 0.4744 | 1:8;2:1;3:3;4:1 |
| 5 | 0.7944 | 0.1203 | 0.1805 | 0.4501 | 0.3590 | 1:6;2:1;3:3;4:2;5:1 |

Selected k=2; cluster sizes: C1 n=12, C2 n=1.

Selected-k sample membership:

| Cluster | n | Sample IDs |
| --- | --- | --- |
| 1 | 12 | GSM3575609;GSM3575613;GSM3575614;GSM3575615;GSM3575617;GSM3575620;GSM3575621;GSM3575623;GSM3575624;GSM3575629;GSM3575635;GSM3575639 |
| 2 | 1 | GSM3575622 |

## Random-gene-set null controls

Each null set contains 5 genes, sampled without replacement and excluding the observed genes. There are 1,000 sets from the 43-gene pool and 1,000 from all 11741 expressed genes. Each set uses the same 1,000 sample subsamples and is re-z-scored across the 13 LPD samples. All selected-k null metrics are in `null_distributions.csv`; empirical two-sided p-values use the plus-one correction.

| Universe | Metric | Observed | Null median (95% range) | Observed percentile | Empirical two-sided p |
| --- | --- | --- | --- | --- | --- |
| 43_gene_pool | cdf_area | 0.2539 | 0.4685 (0.2783-0.5296) | 1.2% | 0.0260 |
| 43_gene_pool | mean_pairwise_consensus | 0.7469 | 0.5324 (0.4710-0.7225) | 98.8% | 0.0260 |
| 43_gene_pool | pac_0_1 | 0.4744 | 0.6282 (0.1410-0.8333) | 22.8% | 0.4575 |
| 43_gene_pool | silhouette_euclidean_z | 0.3930 | 0.2362 (0.1499-0.3707) | 98.4% | 0.0340 |
| 43_gene_pool | silhouette_consensus_distance | 0.7894 | 0.7481 (0.4730-0.9637) | 62.3% | 0.7552 |
| all_expressed_genes | cdf_area | 0.2539 | 0.4431 (0.2456-0.5231) | 3.5% | 0.0719 |
| all_expressed_genes | mean_pairwise_consensus | 0.7469 | 0.5578 (0.4776-0.7549) | 96.5% | 0.0719 |
| all_expressed_genes | pac_0_1 | 0.4744 | 0.6154 (0.1410-0.8333) | 24.2% | 0.4855 |
| all_expressed_genes | silhouette_euclidean_z | 0.3930 | 0.2547 (0.1586-0.4554) | 91.3% | 0.1758 |
| all_expressed_genes | silhouette_consensus_distance | 0.7894 | 0.7557 (0.4824-0.9602) | 60.3% | 0.7952 |

## PCA and descriptive gene panel

The PCA is calculated on the 5 gene-wise z-scored panel and colored by selected consensus cluster (`Fig7_exploratory.png`). The first two components explain 44.6% and 28.6% of variance, respectively. The cluster-wise panel plot is descriptive only: clusters were constructed from these same genes (`Fig7_gene_panel_descriptive.png`).

## Immune/stromal ssGSEA scores

Thirteen exact CellMarker_2024 human signatures were scored by ssGSEA. The leukocyte proxy is the mean of 17 canonical leukocyte-marker z-scores, each standardized across all 33 discovery samples. For k=2, the reported unadjusted tests are two-sided Mann-Whitney U tests with BH q across 13 signatures. For k>2, all pairwise two-sided Mann-Whitney U tests are BH-adjusted together; the `ssgsea_tests` file also reports an omnibus cluster-term F test after leukocyte-proxy adjustment, BH-adjusted across the 13 signatures. Scores are descriptive marker-enrichment scores, not cell fractions.

| Signature | Unadjusted test p | MW BH q (or minimum pairwise q) | Leukocyte-adjusted cluster p | Adjusted BH q | Survives q<0.05? |
| --- | --- | --- | --- | --- | --- |
| T cells | 0.1538 | 0.8 | 0.2513 | 0.8817 | no |
| B cells | 0.9231 | 0.9231 | 0.4682 | 0.8817 | no |
| NK cells | 0.3077 | 0.8 | 0.956 | 0.956 | no |
| Monocytes | 0.3077 | 0.8 | 0.6217 | 0.8981 | no |
| M1 macrophages | 0.9231 | 0.9231 | 0.6981 | 0.9075 | no |
| M2 macrophages | 0.1538 | 0.8 | 0.19 | 0.8817 | no |
| Neutrophils | 0.6154 | 0.9231 | 0.5426 | 0.8817 | no |
| Mast cells | 0.4615 | 0.8571 | 0.4941 | 0.8817 | no |
| Immature dendritic cells | 0.4615 | 0.8571 | 0.3598 | 0.8817 | no |
| Mature dendritic cells | 0.7692 | 0.9231 | 0.4686 | 0.8817 | no |
| Fibroblasts | 0.3077 | 0.8 | 0.2608 | 0.8817 | no |
| Myofibroblasts | 0.9231 | 0.9231 | 0.9339 | 0.956 | no |
| Mesothelial cells | 0.9231 | 0.9231 | 0.9339 | 0.956 | no |

Pairwise Mann-Whitney values and BH q-values are in `ssgsea_pairwise_mann_whitney_selected_k.csv`; sample scores are in `ssgsea_scores_selected_k.csv`.

| Signature | Cluster A | Cluster B | Median A | Median B | MW p | BH q |
| --- | --- | --- | --- | --- | --- | --- |
| T cells | 1 | 2 | 0.1761 | 0.3164 | 0.1538 | 0.8 |
| B cells | 1 | 2 | 0.2017 | 0.1842 | 0.9231 | 0.9231 |
| NK cells | 1 | 2 | -0.2078 | -0.0243 | 0.3077 | 0.8 |
| Monocytes | 1 | 2 | 0.4777 | 0.5021 | 0.3077 | 0.8 |
| M1 macrophages | 1 | 2 | 0.3753 | 0.3565 | 0.9231 | 0.9231 |
| M2 macrophages | 1 | 2 | 0.3092 | 0.4411 | 0.1538 | 0.8 |
| Neutrophils | 1 | 2 | 0.3160 | 0.3414 | 0.6154 | 0.9231 |
| Mast cells | 1 | 2 | 0.2007 | 0.2195 | 0.4615 | 0.8571 |
| Immature dendritic cells | 1 | 2 | 0.2843 | 0.3083 | 0.4615 | 0.8571 |
| Mature dendritic cells | 1 | 2 | 0.2304 | 0.1499 | 0.7692 | 0.9231 |
| Fibroblasts | 1 | 2 | -0.2022 | -0.3319 | 0.3077 | 0.8 |
| Myofibroblasts | 1 | 2 | -0.0933 | -0.0846 | 0.9231 | 0.9231 |
| Mesothelial cells | 1 | 2 | -0.0933 | -0.0846 | 0.9231 | 0.9231 |

## Clinical covariates and stability

The GEO series matrix reports patient ID, short/long-term group, and peritoneal-cell type only; it does not provide per-patient dialysis duration, age, sex, or ultrafiltration-failure status. Its `PRAGUE_1` through `PRAGUE_33` sample descriptions do not define a varying batch factor. Therefore, tests for these clinical variables and batch are not estimable; see `clinical_covariate_availability.csv`.

Across 20 additional resampling seeds, selected k counts were k=2: 20, k=3: 0, k=4: 0, k=5: 0. The selected 12/1 split includes a cluster smaller than four; all tested k values also contain at least one cluster smaller than four. Therefore, the subgroup definition is unstable despite no selected-k changes across these seeds. See `seed_stability.csv` and `cluster_size_lt_4` in `run_metadata.json`.

## Figure legend

**Figure 7. Exploratory molecular subgroups within long-term PD (LPD; n=13).** Consensus clustering used the current five-gene candidate panel (FLT3LG, EBI3, LTB, ADAM19, CST7), gene-wise z-scored across the 13 LPD samples, with 1000 subsamples (pItem=0.8; seed 20261016), Ward linkage with Euclidean distance, and final complete-linkage clustering of 1-consensus. The silhouette-based rule selected k=2 with cluster sizes 12 and 1. CDF/delta-area and silhouette statistics are reported for k=2-5. The PCA is colored by selected cluster. The gene panel is descriptive because the same genes defined the clusters. ssGSEA signatures are exploratory marker scores. A selected cluster has n<4. This analysis is exploratory (n=13), unvalidated, and not suitable for clinical subgroup inference.

## Results text

Within the 13 LPD samples, exploratory consensus clustering of the current five-gene panel selected k=2 by the pre-specified maximum Euclidean silhouette rule (silhouette=0.3930); the resulting cluster sizes were 12 and 1. Results for k=2-5, including consensus CDF area, delta-area, silhouette, and memberships, are shown in the metrics and membership tables. The selected-k Euclidean silhouette was 0.3930; 98.4% of the 1,000 five-gene null sets from the 43-gene pool were below the observed score (empirical two-sided p=0.0340). Against 11,741 expressed-gene null sets, 91.3% were below the observed score (p=0.1758). No immune/stromal signature difference survived BH adjustment, including after leukocyte-proxy adjustment. Duration, age, sex, UF failure, and batch could not be tested because patient-level covariates were unavailable. This analysis is exploratory, has n=13, includes a singleton cluster, and is not validated.

## Code line map

- Input validation and current-panel/LPD checks: [run_fig7_lpd.py](./run_fig7_lpd.py#L84).
- Fixed resampling and Ward/Euclidean consensus clustering: [run_fig7_lpd.py](./run_fig7_lpd.py#L118) and [run_fig7_lpd.py](./run_fig7_lpd.py#L129).
- Fast Ward implementation, verified against SciPy partitions: [run_fig7_lpd.py](./run_fig7_lpd.py#L186).
- Silhouette-based k selection: [run_fig7_lpd.py](./run_fig7_lpd.py#L273).
- Frozen 13-signature definitions and ssGSEA: [run_fig7_lpd.py](./run_fig7_lpd.py#L278) and [run_fig7_lpd.py](./run_fig7_lpd.py#L313).
- Matched random-gene-set null clustering: [run_fig7_lpd.py](./run_fig7_lpd.py#L352) and [run_fig7_lpd.py](./run_fig7_lpd.py#L594).
- ssGSEA tests and covariate availability: [run_fig7_lpd.py](./run_fig7_lpd.py#L708) and [run_fig7_lpd.py](./run_fig7_lpd.py#L793).
