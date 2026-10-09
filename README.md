# Peritoneal Dialysis Effluent and Peritoneal Fibrosis: Transcriptomic Analysis

## Publication-facing summary

This repository documents an exploratory transcriptomic study of peritoneal
dialysis (PD) membrane remodeling, with a fixed candidate-gene panel evaluated
in independent-accession effluent and biopsy cohorts. The analyses do **not**
establish a diagnostic test, causal mechanism, or clinically validated
biomarker panel. Validation cohorts are small, and patient overlap between the
two effluent cohorts cannot be excluded using available GEO metadata.

The analyses summarized here use real GEO expression data. Synthetic shifted
expression matrices used by an earlier validation workflow were withdrawn and
are not evidence. Numerical results below are from saved program outputs.

## Study questions and cohorts

| Accession | Material and comparison | Role | Sample sizes |
|---|---|---|---:|
| [GSE125498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125498) | PD dialysate effluent; long-term PD (LPD) vs short-term PD (SPD) | Discovery | 13 LPD, 20 SPD |
| [GSE62928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62928) | Peritoneal biopsies; severe EPS vs combined non-EPS controls | External biopsy evaluation | 4 EPS, 4 controls |
| [GSE248762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248762) | PD effluent cells; long vintage (LV_NOT_UF + LV_UF) vs short vintage (SV) | Primary effluent cohort | 10 LV, 6 SV |
| [GSE130888](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130888) | Peritoneal cells; long-term vs short-term PD | Secondary effluent cohort | 4 long-term, 6 short-term; 3 normal controls descriptive only |
| [GSE92455](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92455) | Effluent-derived mesothelial cells; epithelioid vs non-epithelioid | Not analyzed | No results reported |

GSE92455 uses a different phenotype contrast and must not be pooled with
long-versus-short-vintage PD analyses or described as validation of that
contrast.

## Discovery and candidate selection

GSE125498 contains 33 effluent samples measured on the Illumina HumanHT-12 v4
platform. The deposited BeadStudio v3 summarized intensities were quantile
normalized with `lumi`; detectable probes were retained according to the GEO
processing description. Differential expression compared 13 LPD with 20 SPD
samples using an absolute log2 fold-change threshold of 0.585 and nominal
P < 0.05. This yielded 813 DEGs (496 higher and 317 lower in LPD). Of these,
44 intersected the measured Naba matrisome genes (31 higher, 13 lower).

Three feature-selection methods—LASSO, linear SVM-RFE, and random forest—were
used for exploratory prioritization. The recorded majority-vote set was
FLT3LG, TNFSF15, LTB, ADAM19, SERPINA10, EBI3, and CXCL14. The strict three-of-
three candidates were TNFSF15 and FLT3LG. The fixed eight-gene panel evaluated
in the subsequent cohorts was **TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1,
SERPINA10, and CST7**. This panel was not altered after examining validation
expression.

The discovery DEG screen uses nominal P-values and is exploratory. In the
supplementary BH-FDR sensitivity analysis, 135/813 total DEGs, 6/44
matrisome DEGs, and 4/7 recorded majority-vote hubs had q < 0.05 (263/813,
10/44, and 5/7, respectively, at q < 0.10).

The nested repeated stratified five-fold analysis used 10 repeats (50 outer
folds; seed 42), with DEG filtering, matrisome intersection, and LASSO,
SVM-RFE, and random-forest feature selection conducted inside the training
folds. Mean out-of-fold AUC/C-index was **0.7700** (repeat-level bootstrap
95% CI **0.7408–0.7946**, SD **0.0453**). WGCNA was not rerun within folds, so
some upstream selection optimism may remain.

## Preregistered external-cohort results

The eight candidate genes and expected direction (higher in long-term PD)
were fixed before the effluent-cohort expression analyses. For scRNA-seq,
Scanpy 1.12.4 and AnnData 0.13.4 were used with group-blind QC thresholds of
at least 200 detected genes and at most 20% mitochondrial counts. Marker-score
rules assigned broad cell types because author cell annotations were not
available in the downloaded matrices. Raw counts were aggregated by GEO
sample, treating each sample as an anonymous patient proxy; cells were not
treated as independent replicates. edgeR 4.4.2 was used for pseudobulk
differential expression. The available GEO metadata does not confirm that
each sample represents a unique participant.

| Cohort and primary contrast | Per-gene direction results | Fixed-panel composite | Frozen cohort verdict |
|---|---|---|---|
| GSE62928: biopsy, severe EPS vs combined non-EPS controls (4 vs 4) | 4/8 directions concordant with discovery | AUC 0.875 (bootstrap 95% CI 0.500–1.000); exact two-sided Mann–Whitney P = 0.1143 | **Not supported.** Biopsy/EPS contrast differs from PD vintage; all eight genes were below the operational expression threshold, which is not a biological detection limit. |
| GSE248762: LV vs SV (10 vs 6) | 7/7 detected genes concordant; SERPINA10 not detected | Not estimable under the frozen rule because the full panel was not detected | **Not supported.** The directional criterion was met, but the full-panel composite and its null were unavailable. This literal label does not mean direction agreement was absent. |
| GSE130888: long-term vs short-term PD (4 vs 6) | 5/6 detected genes concordant; TNFSF15 and SERPINA10 not detected | Not estimable under the frozen rule because the full panel was not detected | **Partially supportive** under the preregistered exact-five-concordant branch. The three normal controls were descriptive only. |

### Effluent per-gene estimates

Validation values are log2 fold-change (95% CI), followed by nominal P and
BH-adjusted q across tested panel genes. “Not detected” means the frozen raw-
count rule was not met; no effect or direction was estimated for that gene.

| Gene | Discovery log2FC (P) | GSE248762 LV vs SV | GSE130888 long vs short |
|---|---:|---|---|
| TNFSF15 | +0.9743 (2.73e-05) | +1.3235 (-0.6243, 3.2712); P=0.1069, q=0.1496 | Not detected |
| FLT3LG | +1.4377 (2.35e-05) | +1.1048 (-0.4786, 2.6883); P=0.0925, q=0.1496 | +0.8252 (-1.7011, 3.3515); P=0.3096, q=0.5516 |
| EBI3 | +0.8034 (0.00244) | +1.3707 (0.2356, 2.5058); P=0.00745, q=0.0261 | +0.8382 (-1.0716, 2.7480); P=0.1940, q=0.5516 |
| LTB | +1.6747 (3.33e-06) | +0.8207 (-0.7052, 2.3467); P=0.1817, q=0.1838 | +0.1963 (-2.7174, 3.1101); P=0.8310, q=0.8310 |
| ADAM19 | +1.2131 (0.000133) | +1.0677 (-0.0594, 2.1949); P=0.0268, q=0.0626 | +0.9231 (-1.2231, 3.0693); P=0.1883, q=0.5516 |
| ECM1 | +1.0341 (0.00701) | +5.8893 (2.2977, 9.4809); P=0.00486, q=0.0261 | -0.6301 (-2.7325, 1.4722); P=0.3677, q=0.5516 |
| SERPINA10 | +1.1659 (0.000653) | Not detected | Not detected |
| CST7 | +1.0123 (0.000500) | +0.4522 (-0.4051, 1.3095); P=0.1838, q=0.1838 | +0.1583 (-2.0375, 2.3541); P=0.8196, q=0.8310 |

These cohort-level verdicts are reported as preregistered. In particular, the
GSE248762 rule is non-monotonic: 7/7 detected genes meet the direction-count
component, but the full-panel score cannot be tested; the exact-five partial
branch does not apply. The verdict has not been changed post hoc.

## Post hoc / exploratory results

These analyses were performed after the preregistered outputs and do not
change the gene panel, thresholds, sample inclusion, or verdicts.

### Detected-gene-only composite scores

The following AUCs use only detected genes and the intended LV-versus-SV
samples. They are **not** the preregistered full-panel tests.

| Cohort | Genes in score | AUC (bootstrap 95% CI) | Mann–Whitney P | Matched-null AUC empirical P |
|---|---:|---:|---:|---:|
| GSE248762, LV vs SV | 7 | 0.7833 (0.5000–1.0000) | 0.0727 | 0.0460 |
| GSE130888, long-term vs short-term | 6 | 0.7083 (0.2917–1.0000) | 0.3524 | 0.0519 |

The earlier GSE130888 exploratory AUC of 0.8056 incorrectly included its
three normal controls as non-LV controls. The corrected estimate uses the
registered four long-term and six short-term PD samples only. No detected-
gene-only score replaces the fixed eight-gene test.

### ECM1 sensitivity and cell-type distribution

In the preregistered GSE248762 all-cell fit, ECM1 had log2FC +5.8893 (95% CI
2.2977–9.4809; P=0.00486). Exploratory edgeR sensitivity fits estimated
log2FC +5.0152 (95% CI 1.4166–8.6138; P=0.00839) after excluding
GSM7919584, and +7.2997 (95% CI 3.1948–11.4047; P=0.00139) in
mesothelial-only pseudobulk.

For GSM7919584, 115,870 raw ECM1 counts were observed in 5,732 QC-retained
cells. Marker-labeled mesothelial cells accounted for 76,200 counts and
neutrophil-labeled cells for 38,545; the remaining marker categories
accounted for 1,125. These count distributions and marker overlaps do not
establish ambient contamination, misannotation, or a cell-intrinsic
biological effect. No sample was removed from the preregistered primary
analysis.

## Independence, interpretation, and limitations

The available GEO metadata, sample titles, and GSM-keyed supplementary
filenames do not provide a participant crosswalk between GSE248762 and
GSE130888. Their independence is **unknown**, not confirmed. Author
confirmation or a de-identified participant crosswalk would be needed to
resolve overlap. The overall evidence therefore does not establish supportive
independent validation of the full panel.

Interpret results in light of the small patient counts, wide confidence
intervals, single-cell dropout, different tissue/phenotype contrasts, and
possible donor overlap. Composition sensitivity analyses cannot establish
cell-intrinsic expression and used no covariates beyond cell-type proportions.
The nominal-P discovery filter and unresolved WGCNA/feature-pool provenance
are additional limitations. WGCNA was not repeated inside the nested CV.
Wide intervals indicate uncertainty, not proof of no effect. No experimental
follow-up or clinical-performance study is reported here.

All counts and numerical summaries above are analysis outputs; rounded values
are shown for readability. The data accessions link to GEO for the deposited
cohort records and supplementary data.
