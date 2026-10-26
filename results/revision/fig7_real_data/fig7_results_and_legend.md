# Figure 7: exploratory peritoneal effluent panel-score analysis

Analysis run: 2026-10-16 (UTC).

## Data and methods

Frozen panel genes: FLT3LG, EBI3, LTB, ADAM19, CST7.
Discovery was restricted to the 13 long-term PD (LPD) peritoneal-effluent samples in GSE125498. The 13 GEO sample accessions were: GSM3575609, GSM3575613, GSM3575614, GSM3575615, GSM3575617, GSM3575620, GSM3575621, GSM3575622, GSM3575623, GSM3575624, GSM3575629, GSM3575635, GSM3575639. The 20 short-term samples were not included in the discovery score or correlations. For each cohort separately, gene-wise z-scores were calculated across that cohort's analyzed samples and the panel score is their mean.
GSE125498 raw non-normalized Illumina intensities were log2 transformed, quantile normalized across all 33 deposited arrays, mapped using GPL10558 official annotation, and duplicate gene symbols were collapsed by retaining the probe with highest mean normalized expression. The 13 LPD samples were then selected. scRNA-seq raw gene-by-cell matrices were aggregated to one pseudobulk count vector per GEO sample; counts were library-size normalized to CPM and log2(CPM+1) transformed. Cohorts were processed and scored separately without merging or batch correction. For GSE57070, GEO's processed log2 series-matrix values were used as deposited, without a second quantile normalization; its raw CEL archive was downloaded and checksum-recorded but was not the matrix used for scoring.
The frozen 13 human CellMarker_2024 signatures were scored with GSEApy ssGSEA (rank normalization, minimum detectable signature size 3). Scores are marker-enrichment scores, not cell fractions. For cohorts with n>=6, all 13 Spearman tests have within-cohort Benjamini-Hochberg q-values; 95% CIs for rho are percentile bootstrap intervals (2,000 resamples, fixed seed). For n<6, only descriptive rho is reported; p-values, q-values, and bootstrap CIs are withheld because asymptotic and bootstrap inference is not credible at that sample size.
Panel/signature overlap was audited: T cells: LTB. Correlations with signatures sharing panel genes should be interpreted with this dependence in mind.
The LTB circularity sensitivity was run by removing LTB from the T-cell marker set while keeping the five-gene panel unchanged; all other signature definitions and scoring settings were unchanged.
No consensus clustering was run: the requested primary inference is continuous, and with small cohorts unsupervised grouping would not improve the prespecified continuous-score test. Consequently, cluster sizes, cluster null distributions, and clustering figures are not applicable.

## Dataset audit and inclusion

Audited 13 relevant GEO/ArrayExpress candidates. Source series matrices and supplementary files are staged in `data/raw/fig7_real_data/` and checksum-verified. The GSE130888 and GSE248762 scRNA-seq archives, and other pre-existing source files, have unknown original download dates. GSE57070_RAW.tar was downloaded on the run date, but scoring used the GEO processed log2 matrix. See `download_manifest.csv` for URLs, byte sizes, SHA-256 hashes, and acquisition-date notes.
SRA search (human peritoneal dialysis/effluent/fibrosis terms) returned 95 run records representing nine BioProjects. Three duplicate GEO projects were GSE130888, GSE248762, and GSE297829. Other projects included noneligible biopsy small-RNA, effluent exosomal miRNA, heart tissue, sepsis blood TCR, an animal model, and E-MTAB-5462. Its SDRF/IDF are public, but no processed gene-level matrix is listed, and reprocessing its FASTQ files would require a separate alignment/quantification workflow, so it is recorded but not analyzed.
GSE130888 contributed 10 effluent samples (6 short-term, 4 long-term); the 3 normal tissue samples were excluded. GSE248762 contributed 16 overnight-effluent samples (6 SV, 6 LV_NOT_UF, 4 LV_UF). GSE57070 contributed six paired treatment/control arrays from three patients; the paired arrays were averaged within patient before scoring (n=3). It is not treated or cited as validation and is omitted from the Figure 7 validation/concordance panels because its rank correlations are necessarily limited to 1, 0.5, or -0.5 and are not informative at this n. GSE297829 is a single patient after three years of PD and is not a valid association validation cohort. E-MTAB-5439 has 10 effluent-derived sorted-cell arrays (five CD14+ macrophage and five CD1c+ dendritic-cell samples from five donors), but it is not bulk effluent and BioStudies lists raw CEL files and a final gene list rather than a processed expression matrix. E-MTAB-5461 is a microRNA array and cannot score this mRNA panel. The biopsy and in-vitro/cell-selected studies in the audit were not treated as bulk effluent validation.

## Correlations

| Dataset | n | Signature | rho (95% bootstrap CI) | p | BH q |
|---|---:|---|---:|---:|---:|
| GSE125498 | 13 | T cells | 0.857 [0.521, 0.977] | 0.0001809 | 0.002352 |
| GSE125498 | 13 | B cells | 0.484 [-0.138, 0.842] | 0.09413 | 0.4079 |
| GSE125498 | 13 | NK cells | 0.390 [-0.302, 0.865] | 0.1876 | 0.6096 |
| GSE125498 | 13 | Monocytes | 0.522 [-0.020, 0.832] | 0.06729 | 0.4079 |
| GSE125498 | 13 | M1 macrophages | -0.225 [-0.659, 0.419] | 0.4593 | 0.8998 |
| GSE125498 | 13 | M2 macrophages | 0.099 [-0.531, 0.667] | 0.7479 | 0.8998 |
| GSE125498 | 13 | Neutrophils | -0.115 [-0.608, 0.413] | 0.7074 | 0.8998 |
| GSE125498 | 13 | Mast cells | -0.066 [-0.801, 0.613] | 0.8305 | 0.8998 |
| GSE125498 | 13 | Immature dendritic cells | 0.093 [-0.459, 0.703] | 0.7615 | 0.8998 |
| GSE125498 | 13 | Mature dendritic cells | -0.022 [-0.678, 0.569] | 0.9432 | 0.9432 |
| GSE125498 | 13 | Fibroblasts | 0.225 [-0.503, 0.905] | 0.4593 | 0.8998 |
| GSE125498 | 13 | Myofibroblasts | 0.099 [-0.555, 0.821] | 0.7479 | 0.8998 |
| GSE125498 | 13 | Mesothelial cells | 0.099 [-0.560, 0.803] | 0.7479 | 0.8998 |
| GSE130888 | 10 | T cells | 0.939 [0.623, 1.000] | 5.484e-05 | 0.0007129 |
| GSE130888 | 10 | B cells | 0.842 [0.371, 1.000] | 0.00222 | 0.00962 |
| GSE130888 | 10 | NK cells | 0.709 [0.143, 0.960] | 0.02167 | 0.07041 |
| GSE130888 | 10 | Monocytes | 0.503 [-0.208, 0.875] | 0.1383 | 0.2569 |
| GSE130888 | 10 | M1 macrophages | 0.539 [-0.203, 0.961] | 0.1076 | 0.2331 |
| GSE130888 | 10 | M2 macrophages | -0.139 [-0.873, 0.593] | 0.7009 | 0.9074 |
| GSE130888 | 10 | Neutrophils | 0.370 [-0.438, 0.846] | 0.2931 | 0.4762 |
| GSE130888 | 10 | Mast cells | 0.042 [-0.709, 0.723] | 0.9074 | 0.9074 |
| GSE130888 | 10 | Immature dendritic cells | -0.248 [-0.827, 0.559] | 0.4888 | 0.706 |
| GSE130888 | 10 | Mature dendritic cells | 0.879 [0.444, 1.000] | 0.0008139 | 0.00529 |
| GSE130888 | 10 | Fibroblasts | -0.673 [-0.988, 0.033] | 0.03304 | 0.08591 |
| GSE130888 | 10 | Myofibroblasts | -0.067 [-0.720, 0.650] | 0.8548 | 0.9074 |
| GSE130888 | 10 | Mesothelial cells | -0.067 [-0.701, 0.635] | 0.8548 | 0.9074 |
| GSE248762 | 16 | T cells | 0.956 [0.818, 0.997] | 7.76e-09 | 1.009e-07 |
| GSE248762 | 16 | B cells | 0.715 [0.315, 0.895] | 0.001863 | 0.006054 |
| GSE248762 | 16 | NK cells | 0.891 [0.646, 0.982] | 3.619e-06 | 2.352e-05 |
| GSE248762 | 16 | Monocytes | 0.341 [-0.244, 0.845] | 0.1959 | 0.3184 |
| GSE248762 | 16 | M1 macrophages | 0.194 [-0.389, 0.739] | 0.4713 | 0.557 |
| GSE248762 | 16 | M2 macrophages | -0.341 [-0.735, 0.221] | 0.1959 | 0.3184 |
| GSE248762 | 16 | Neutrophils | -0.203 [-0.671, 0.376] | 0.451 | 0.557 |
| GSE248762 | 16 | Mast cells | -0.306 [-0.713, 0.231] | 0.2493 | 0.36 |
| GSE248762 | 16 | Immature dendritic cells | -0.459 [-0.805, 0.070] | 0.07384 | 0.16 |
| GSE248762 | 16 | Mature dendritic cells | 0.694 [0.234, 0.937] | 0.002854 | 0.00742 |
| GSE248762 | 16 | Fibroblasts | -0.747 [-0.925, -0.349] | 0.0008821 | 0.003822 |
| GSE248762 | 16 | Myofibroblasts | -0.074 [-0.583, 0.489] | 0.7867 | 0.7867 |
| GSE248762 | 16 | Mesothelial cells | -0.074 [-0.611, 0.533] | 0.7867 | 0.7867 |

### T-cell result after removing LTB

The sensitivity analysis below rescored the T-cell signature without LTB. The five-gene panel remains unchanged. GSE57070 is omitted because n=3 supports no useful inferential check.

| Dataset | n | T-cell rho (95% bootstrap CI) | p | BH q (13-signature family) |
|---|---:|---:|---:|---:|
| GSE125498 | 13 | 0.857 [0.573, 0.977] | 0.0001809 | 0.002352 |
| GSE130888 | 10 | 0.939 [0.623, 1.000] | 5.484e-05 | 0.0007129 |
| GSE248762 | 16 | 0.956 [0.826, 0.997] | 7.76e-09 | 1.009e-07 |

In discovery, removing LTB left rho at 0.857 (p=0.0001809; BH q=0.002352). Thus the observed T-cell association is not driven solely by the shared LTB marker, although that does not rule out broader immune-content confounding.

### Leukocyte-content caution

The five-gene panel contains several immune-associated genes. A panel-score association with an immune-marker proxy is therefore compatible with the panel reflecting leukocyte content rather than a distinct fibrosis program; ssGSEA scores are enrichment scores, not cell fractions. The following exploratory check correlates the panel score with the mean of the nine non-T-cell immune ssGSEA scores (the T-cell signature is left out to avoid the direct LTB overlap):

| Dataset | n | Spearman rho | 95% bootstrap CI | p | BH q across cohorts |
|---|---:|---:|---:|---:|---:|
| GSE125498 | 13 | 0.330 | [-0.270, 0.784] | 0.2713 | 0.2713 |
| GSE130888 | 10 | 0.612 | [-0.038, 0.974] | 0.05997 | 0.1799 |
| GSE248762 | 16 | 0.318 | [-0.262, 0.834] | 0.2306 | 0.2713 |

None of these three panel-score/proxy correlations survived BH correction (q range 0.180–0.271); this does not establish that the panel is independent of immune content, given the small samples and marker-score proxy. Conversely, these proxy results do not demonstrate that leukocyte abundance explains the panel associations. The negative fibroblast association in GSE248762 is compatible with a compositional dilution interpretation but does not establish fibroblast depletion or a causal mechanism.

The full machine-readable table includes all 13 signatures, including non-significant associations and descriptive GSE57070 values; that n=3 cohort is not validation and is excluded from the report table and figures. Directional agreement in GSE130888 and GSE248762 is exploratory cross-cohort concordance, not identical-method replication: both validation cohorts are single-cell-derived pseudobulk, whereas discovery is bulk microarray. Small n (13, 10, 16) and wide intervals limit precision.

## Covariates and limitations

Clinical association models use panel score as the outcome and include the leukocyte proxy (mean of the 10 immune ssGSEA scores) as an adjustment term. Reported coefficients are unstandardized OLS effects with 95% confidence intervals, p-values, and within-dataset BH q-values. None of the five fitted clinical-covariate associations survived BH correction (all q >= 0.0655). These models test clinical covariate associations conditional on the proxy; they do not directly test whether the panel is independent of leukocyte content. Where metadata do not provide a patient-level value, no value was inferred.

| Dataset | Covariate | n | Adjusted effect (95% CI) | p | BH q |
|---|---|---:|---:|---:|---:|
| GSE130888 | long_vintage_vs_short | 10 | 0.665 [-0.056, 1.386] | 0.06548 | 0.06548 |
| GSE248762 | age_years | 16 | 0.001 [-0.029, 0.031] | 0.9381 | 0.9381 |
| GSE248762 | sex_M_vs_F | 16 | -0.063 [-1.109, 0.983] | 0.8983 | 0.9381 |
| GSE248762 | long_vintage_vs_short | 16 | 0.644 [-0.249, 1.536] | 0.1432 | 0.5729 |
| GSE248762 | UF_failure_vs_no_UF_in_long_vintage | 10 | 0.514 [-1.206, 2.235] | 0.5024 | 0.9381 |

Exact dialysis vintage and batch were not available as per-sample values in the analyzed cohorts. Previous peritonitis was available only in GSE57070. In GSE248762, long/short vintage and UF status are encoded by design group; age and sex were assessed only if present as sample-level series-matrix fields. In the LPD-only discovery subset, a between-vintage-group test is impossible because the subset contains only long-term patients. These results remain underpowered and exploratory; neither the discovery correlations nor the validation cohorts constitute clinical validation. GSE57070 has only three independent patients; previous peritonitis is in its sample metadata, but no clinical covariate model was fitted because n=3. Age, sex, vintage, UF status, transport status, and batch were not available for that cohort.

## Figure legend

**Figure 7. Exploratory association of the frozen five-gene panel score with immune/stromal marker-enrichment scores in peritoneal dialysis effluent.** (A) Spearman rho heatmap for the 13 frozen CellMarker_2024 signatures in GSE125498 (n=13), GSE130888 (n=10), and GSE248762 (n=16); these are separate cohorts, not identical-method validation. GSE57070 (n=3) is excluded from these panels and is not cited as validation. (B) Three strongest absolute discovery associations, shown descriptively; these signatures are not independent of the score-construction process. (C) Cohort-specific Spearman rho and bootstrap 95% confidence intervals across the discovery and two separate single-cell pseudobulk cohorts. The T-cell sensitivity analysis excludes LTB from that signature. ssGSEA values represent marker enrichment, not cell fractions. All analyses are exploratory, small-sample, and not validated for clinical use.
