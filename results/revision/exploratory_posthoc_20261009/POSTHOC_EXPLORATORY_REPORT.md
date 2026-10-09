# Post hoc / exploratory audit outputs

**Date:** 2026-10-09  
**Scope:** GSE248762 ECM1 raw-count localization and sensitivity fits; detected-gene
composites for GSE248762 and GSE130888; GEO metadata-based donor overlap check.

These analyses were run after the preregistered outputs were finalized. They do
not change the frozen gene panel, thresholds, sample inclusion, or verdicts.
All outputs and scripts are contained in this directory; source cohorts were
read from their existing saved GEO data.

## Detected-gene composites

The score is the mean of per-gene standardized edgeR log-CPM values, with each
gene's mean and standard deviation estimated from LV and SV samples in the
primary contrast. GSE130888 NORMAL samples are retained in the patient-score
file for display, transformed using those LV/SV reference statistics, and
excluded from the primary AUC and matched-null AUC. Nulls consist of 1,000
detected gene sets matched to the score genes by mean-expression decile. Seed
42 was used for bootstrap resampling and null-set selection.

| Cohort and contrast | Genes in exploratory score | AUC (95% bootstrap CI) | Mann–Whitney P | Matched-null AUC empirical P |
|---|---:|---:|---:|---:|
| GSE248762, LV vs SV | 7 | 0.7833 (0.5000–1.0000) | 0.0727 | 0.0460 |
| GSE130888, long-term vs short-term PD | 6 | 0.7083 (0.2917–1.0000) | 0.3524 | 0.0519 |

For completeness, matched-null empirical P for the mean detected-gene log2FC
was 0.0010 in GSE248762 and 0.0609 in GSE130888. The GSE130888 AUC is 17/24
correctly ordered LV–SV pairs, equal to 0.7083; its no-tie AUC grid increment
is 1/24 (and the general half-tie increment is 1/48). The previous 0.8056
value came from treating all samples not labeled LV as controls, which
incorrectly included the three NORMAL samples. The reported sample sizes had
shown 4 LV and 6 SV, but the score helper had compared 4 LV with 9 non-LV.
The corrected implementation restricts score standardization, observed AUC,
and matched-null AUC calculations to LV and SV. The GSE130888 score CSV
retains all NORMAL rows for transparency and marks them
`auc_included=False`.

GSE130888 per-sample composite scores (NORMAL samples are retained here for
descriptive transparency but excluded from the primary LV-vs-SV AUC) are:

| Sample | Group | Composite | In AUC |
|---|---|---:|---|
| GSM3755687 | SV | 0.550069 | Yes |
| GSM3755688 | SV | 0.570549 | Yes |
| GSM3755689 | SV | -0.204885 | Yes |
| GSM3755690 | SV | -1.236510 | Yes |
| GSM3755691 | SV | 0.418011 | Yes |
| GSM3755692 | SV | -0.961697 | Yes |
| GSM3755693 | NORMAL | -2.671807 | No |
| GSM3755694 | NORMAL | -0.999649 | No |
| GSM3755695 | NORMAL | -1.989649 | No |
| GSM3755696 | LV | -0.037438 | Yes |
| GSM3755697 | LV | 0.740187 | Yes |
| GSM3755698 | LV | -0.659753 | Yes |
| GSM3755699 | LV | 0.821467 | Yes |

The GSE130888 bootstrap interval is 0.2917–1.0000. The GSE248762 contrast
contains only LV and SV samples and was unchanged; its 47/60 correctly ordered
LV–SV pairs give AUC 0.7833. These wide intervals reflect the small cohorts.
These are explicitly **post hoc detected-gene composites, not the
preregistered test**. They do not change GSE248762's preregistered
**Not supported** verdict or GSE130888's **Partially supportive** verdict.

Files: `detected_gene_composites_summary.json`,
`detected_gene_composites_summary.csv`,
`detected_gene_composites_patient_scores.csv`, and
`detected_gene_composites_matched_null.csv`. The metadata-only sample
identifiability findings are saved in `cohort_independence_audit.json` and
`cohort_independence_audit.csv`.

## GSE248762 ECM1 count localization

The real GEO 10x count archive was re-read for all 16 samples. The original
group-blind QC thresholds and marker-score cell labels were reused without
modification. QC-retained cell counts matched the saved primary run for every
sample; no sample was removed.

In GSM7919584 (LV_UF-2), 115,870 raw ECM1 counts were detected across 5,732
QC-retained cells; 1,331 cells (23.2%) had nonzero ECM1 counts. The highest
1% of cells carried 30.5% of total ECM1 counts and the highest 10% carried
97.7%; the maximum was 1,909 counts in one cell. This is not uniform
low-level expression across all retained cells, but it is also distributed
across many cells rather than being confined to one cell.

| Marker-based cell type | QC cells | ECM1 raw counts | ECM1-positive cells |
|---|---:|---:|---:|
| T/NK | 4,136 | 883 | 361 |
| B | 162 | 18 | 17 |
| Monocyte/macrophage | 241 | 117 | 57 |
| Neutrophil | 383 | 38,545 | 269 |
| Mesothelial | 753 | 76,200 | 614 |
| Fibroblast | 5 | 87 | 3 |
| Other | 52 | 20 | 10 |

Thus 114,745 of 115,870 ECM1 counts in this sample were assigned by the
marker-based labels to mesothelial or neutrophil cells. The fibroblast label
contains only five cells and is not a reliable basis for cell-type inference.

The same sample's labeled mesothelial cells included 5,060 MSLN counts
(557 positive cells) and 2,483 UPK3B counts (315 positive cells); among
ECM1-positive mesothelial cells, 497 also expressed MSLN and 280 UPK3B.
Neutrophil-labeled cells included 48,689 S100A8 counts (328 positive cells),
280,153 S100A9 counts (380 positive cells), and 128,157 LCN2 counts
(303 positive cells); among ECM1-positive neutrophil-labeled cells, 249
expressed S100A8, 267 S100A9, and 256 LCN2. MSLN/UPK3B signal also appeared
in neutrophil-labeled cells, and S100A8/S100A9/LCN2 signal appeared in
mesothelial-labeled cells. This mixed marker pattern makes the labels
descriptive rather than definitive in this sample.

Other LV samples also had high ECM1 totals: GSM7919585 had 72,808 counts in
731 positive cells, and GSM7919587 had 20,334 counts in 1,491 positive cells.
GSM7919584 is therefore not the only sample with elevated totals. The
sample-by-cell-type table and all-sample concentration metrics are in
`ecm1_posthoc_by_sample_celltype.csv` and
`ecm1_posthoc_cell_concentration.csv`; marker counts and ECM1 co-expression
are in `ecm1_posthoc_marker_comparison.csv`.

The filtered-cell matrices do not contain empty-droplet profiles or an
ambient-RNA estimate. Neither these count distributions nor the marker overlap
proves ambient RNA, doublets, annotation error, or a cell-intrinsic biological
effect. They show that ECM1 counts are concentrated in cells assigned to
mesothelial and neutrophil categories, with substantial mixed marker
expression; more source-level validation would be needed to explain why.

## ECM1 edgeR sensitivity fits

Both exploratory models use patient-level raw-count pseudobulks and edgeR
quasi-likelihood tests, retaining all samples except where the sensitivity
explicitly omits GSM7919584. Filtering applies `filterByExpr` and then forces
ECM1 to remain in the fit, even if the generic filter would exclude it. The
approximate 95% intervals use one-degree-of-freedom quasi-likelihood
F-statistic inversion (`SE = |log2FC| / sqrt(F)` and a t critical value with
the adjusted residual degrees of freedom); they are not the primary script's
fitted-covariance intervals.

The mesothelial-only patient-level pseudobulk counts underlying the fit were:

| Group | Patients | Raw ECM1 pseudobulk counts |
|---|---:|---|
| SV | 6 | 78, 115, 80, 68, 331, 52 |
| LV | 10 | 34, 76,200, 17,694, 126, 905, 112, 224, 53, 276, 111 |

The low and zero SV counts are visible in these raw sums. The cell-type
selection and patient restriction are exploratory and do not change the
preregistered patient inclusion or verdict.

| Exploratory analysis | n SV | n LV | log2FC (LV−SV) | Approx. 95% CI | P |
|---|---:|---:|---:|---:|---:|
| All-cell pseudobulk, excluding GSM7919584 | 6 | 9 | 5.0152 | [1.4166, 8.6138] | 0.00839 |
| Mesothelial-only pseudobulk | 6 | 10 | 7.2997 | [3.1948, 11.4047] | 0.00139 |

These exploratory estimates remain positive after the specified sensitivity
changes. They do not identify the cause of the ECM1 signal or establish
cell-intrinsic expression. The preregistered all-sample analysis still includes
GSM7919584; no primary sample was removed.

The saved fit table is `ecm1_posthoc_edgeR_sensitivities.csv`, and the
reproducible R script is `ecm1_sensitivity_edgeR.R`.

## Cohort independence check

The GEO sample metadata for GSE248762 (16 samples) and GSE130888 (13 samples)
contains group labels, GSM accessions, and, for GSE248762, demographic
characteristics, but no reported patient/donor/participant crosswalk.
Sample titles use within-study group labels and replicate numbers only. The
sample-level `supplementary_files` column in the derived metadata CSV is empty,
but the downloaded GEO family SOFT files do list raw 10x supplementary
matrices (three files per sample); their names are keyed by GSM accession and
within-cohort group labels, not a shared donor identifier. Consequently,
patient overlap cannot be checked from these GEO records or downloaded
supplementary files, and independence remains **unknown**.

Resolving this requires a de-identified participant crosswalk shared by the
study teams, explicit overlap disclosure from the authors, or another
reliable linking key present in both cohorts. No overlap is inferred from
matching demographic values or disease-group labels.
