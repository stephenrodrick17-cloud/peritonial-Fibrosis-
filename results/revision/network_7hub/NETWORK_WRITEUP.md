# lncRNA-miRNA-mRNA regulatory network and single-cell co-expression

## Methods

The analysis was restricted to the 7 validated hub genes read from `results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`. The frozen 5-gene panel was read from `results/revision/panel_v2/panel_v2.json`; no gene symbols are repeated in this prose.

Positive validated miRNA-mRNA evidence was queried through multiMiR for miRTarBase, miRecords, and TarBase, and predicted evidence was kept separate. Negative records were excluded from positive support. Tier A required a positive validated record whose experiment field explicitly named a reporter/luciferase assay, western blot, qPCR, or qRT-PCR. Tier B required at least two ENCORI CLIP-seq experiments for the hub mRNA. Tier C was restricted to multiMiR prediction-only pairs with no positive validated record or Tier A/B support. Validated rows lacking a qualifying functional assay were counted separately as unclassified, not upgraded to Tier A.

ENCORI interactions were retrieved against hg38. The hub-mRNA API query requested records with at least one CLIP experiment and no target-prediction program requirement; Tier B was then assigned using `clipExpNum >= 2`. The ENCORI miRNA-lncRNA archive was retrieved as a raw archive and matched to miRNAs linked to the validated hubs. High-confidence axes required a Tier A or B hub-mRNA edge and an lncRNA edge with `clipExpNum >= 2`. All remaining database-backed candidate axes were low confidence. A recorded target site with zero CLIP support is a prediction-level lncRNA candidate, not experimentally established binding.

The ENCORI ceRNA endpoint was also queried separately for each hub (`miRNAnum >= 2`, `pval <= 0.01`, `fdr <= 0.01`, `pancancerNum = 0`). These shared-miRNA network results were not used to define the lncRNA-miRNA-hub axes because they are a separate, cancer-network-derived evidence product.

The relevant external miRNA-seq dataset was GSE182736 / PRJNA757643. ENA run metadata were saved, but sample-level miRNA abundance was not quantified. Thus miRNA detection is NOT VERIFIED and was not used to filter candidate network axes.

For GSE248762, the saved QC-retained cell annotations and sample-lineage pseudobulk count matrices were reused. An lncRNA was eligible for co-expression only if it was detected in at least 1.0% of QC-retained cells and had nonzero counts in at least 3 distinct samples. Each eligible lncRNA-hub pair was tested per lineage using Spearman correlation on log2(CPM + 1) values from sample-level pseudobulk profiles; samples with fewer than 20 cells in a lineage were excluded and at least 5 samples were required. The adjusted analysis is partial Spearman correlation controlling for the rank of log(1 + lineage cell count). BH correction was applied over all testable pairs across lineages separately for unadjusted and adjusted P values. A cell-level correlation was calculated only as a sensitivity analysis; cells are not independent biological replicates.

## Results

The primary network result is 1260 high-confidence lncRNA-miRNA-hub axes across 7 hubs, 172 miRNAs, and 14950 lncRNAs. The 322201 low-confidence axes are reported only as a candidate background, not as supported interactions; many include prediction-level lncRNA sites. The complete table contains 323461 axes in total, reconciled in `network_axes.csv` and `network_reconciliation.csv`.

MultiMiR/ENCORI tier counts (best-class, unique miRNAs summed over hubs) were Tier A = 5, Tier B = 155, and Tier C = 389; Tier C may be incomplete because the predicted query returned exactly 1,000 rows and truncation by a service cap could not be verified. Per-hub Tier A counts are shown below; full hub-specific counts, unclassified validation, and source detail are in `hub_miRNA_tier_counts.csv` and `miRNA_mRNA_evidence.csv`. These are database evidence counts, not evidence of regulation in peritoneal samples.

| Hub gene | Tier A miRNAs |
|---|---:|
| TNFSF15 | 0 |
| FLT3LG | 0 |
| EBI3 | 0 |
| LTB | 0 |
| ADAM19 | 3 |
| ECM1 | 2 |
| CST7 | 0 |

Of the 14 multiMiR sources listed in the saved source metadata, 3 were not queried: miR2Disease, Pharmaco-miR (Verified Sets), PhenomiR. Current-release miRTarBase was not retrieved directly; the queried multiMiR miRTarBase snapshot is older and must not be represented as the current release.


The separate ENCORI ceRNA queries returned parsable results for 3 of 7 queried hubs; 4 queries returned API errors. The successful responses contained 1247 records representing 1247 unique gene-pair links, including 4 rows with an lncRNA partner. These are not counted as validated regulation and are tabulated separately in `ENCORI_ceRNA_pairs.csv` and `ENCORI_ceRNA_query_status.csv`.

Of 14950 lncRNAs present in the candidate axes, 13949 matched an exact Ensembl gene ID in the GSE248762 count matrix and 1271 met the predeclared co-expression detection threshold. The pseudobulk analysis produced 38446 testable lncRNA-hub-lineage pairs; 4072 had nominal unadjusted P < 0.05 and 175 had global BH q < 0.05. After adjusting for log cell count, 3108 had nominal P < 0.05 and 14 had global BH q < 0.05. All nominal and adjusted results, including nonsignificant pairs, are in `lncRNA_mRNA_correlation_results.csv`. The adjusted results are the most stringent co-expression evidence available in this analysis, but remain exploratory because lineage-level pseudobulk sample counts are small. The nominal unadjusted fraction is 4072/38446 (10.59%); this exceeds the nominal 5% reference, but is descriptive rather than a calibrated chance estimate because tests are dependent and share genes and samples.

The separate cell-level sensitivity table contains 39272 estimable lncRNA-hub-lineage correlations. It is reported descriptively because individual cells are nested within donors and are not independent replicates.

The external database query status, including the unavailable direct current miRTarBase release and the older source versions available through multiMiR, is documented in `database_query_status.csv`. Current-release results must not be inferred from the older multiMiR backend.

## Limitations

This is a database-supported candidate network, not proof of post-transcriptional regulation in peritoneal dialysis. Only the high-confidence class is treated as the primary network result; low-confidence axes are a background set and many rely on prediction-level lncRNA sites. The 1,000-row predicted multiMiR result may be truncated, so Tier C counts may be incomplete. The direct current miRTarBase release was not retrieved; the multiMiR backend contains older source snapshots, and three of the fourteen listed sources were not queried. ENCORI CLIP support does not establish functional repression, and prediction-only edges and ceRNA-network associations remain hypotheses. 4 ceRNA API queries returned errors and cannot be interpreted as zero biological interactions. The relevant public effluent miRNA-seq study has six archived runs, but no sample-level miRNA abundance table was verified, so network miRNA expression and co-expression-dependent axes remain unverified. The single-cell data are sparse for many lncRNAs, use a small number of donors, and represent effluent pseudobulk rather than tissue. Pseudobulk correlations are exploratory; the cell-level sensitivity is especially vulnerable to pseudoreplication and donor composition.

## Figure legends

**Network figure.** Candidate lncRNA-miRNA-hub axes are shown by confidence class. Node colors indicate RNA class; edge styles distinguish the best hub-mRNA evidence class. Predicted and CLIP-supported interactions are candidates, not validated regulation.

**Tier bar plot.** Counts of unique miRNAs by hub and best evidence tier. Tier definitions and evidence source versions are given in Methods.

**Detection plot.** Distribution of the fraction of QC-retained GSE248762 cells with nonzero lncRNA counts. The vertical line indicates the cell-level detection component of the stated eligibility threshold.

**Correlation plot.** Unadjusted sample-level Spearman coefficients are compared with partial Spearman coefficients adjusted for log lineage cell count. Each point is an eligible lncRNA-hub-lineage pair. The layered network plot foregrounds high-confidence axes and displays at most 30 low-confidence background axes per hub for readability; the full axis set is in `network_axes.csv`.

**Current miRTarBase access.** NOT RETRIEVED is recorded in the query-status table; the retrieved multiMiR source snapshot is not a substitute for the current release.
