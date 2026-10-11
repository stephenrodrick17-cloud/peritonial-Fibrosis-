# Corrected quantitative summary and audit status

**Updated:** 2026-10-09. This file supersedes the legacy summary that reported
synthetic GSE62928 validation, four-model/XGBoost selection, and unverified
network/enrichment results. Those legacy claims are withdrawn. The numbers
below refer only to saved analyses cited here; this is not a claim that every
biological uncertainty has been resolved.

## Discovery: GSE125498

| Measure | Audited result |
|---|---:|
| Samples | 33 (13 long-term PD, 20 short-term PD) |
| Unique-gene analysis universe | 11,741 |
| Primary nominal DEGs | 813 (496 up, 317 down) |
| Up/down ECM DEGs | 31 / 13 (44 total) |
| ECM genes measured in the universe | 428 |
| Up-ECM fold enrichment | 1.7145x |
| Hypergeometric over-representation P | 0.0024028 |
| 10,000-label-permutation Monte Carlo empirical P | 0.0306969 |

The hypergeometric and permutation P-values test different null models. The
hypergeometric test samples 496 up-DEG labels from the 11,741-gene universe,
with 428 measured Matrisome genes; its null expected overlap is
`496 * 428 / 11741 = 18.09`, versus 31 observed. The permutation test shuffles
sample labels and re-runs the expression/DEG thresholding on the measured
Matrisome genes. The latter is Monte Carlo, not an exact permutation
enumeration. The P-values are not expected to match.

The FDR sensitivity table reports 135/813 nominal primary DEGs, 6/44 ECM
DEGs, and 4/7 recorded hubs at BH q<0.05. The separate 143 figure is the
number of BH-significant genes across the full 11,741-gene tested universe;
these are different denominators, not contradictory counts. The global-null
expectation `11,741 * 0.05 = 587` applies to nominal P<0.05 tests before the
fold-change filter and is not an estimate of false discoveries among the 813
selected DEGs.

## WGCNA and machine-learning feature selection

Dynamic WGCNA assigned 4,549 of 5,001 input genes to non-grey modules. The
top trait-associated module was MEmagenta (37 genes; r=0.4804, P=0.00466),
with no ECM genes. MEblack contained 25 of the 31 upregulated ECM DEGs but
was not significantly associated with stage (r=0.0285, P=0.8748). WGCNA was
not recomputed inside the repeated-CV folds; the reported nested estimate can
therefore retain some upstream selection optimism.

The saved three-model summaries verify:

| Feature pool | LASSO | SVM-RFE | RF | >=2/3 consensus | 3/3 |
|---|---:|---:|---:|---:|---:|
| 624 DEG ∩ non-grey WGCNA genes | 24 | 10 | 10 | 8 | 0 |
| 38 ECM-DEG ∩ non-grey WGCNA genes | 11 | 10 | 10 | 8 | 2 |

The two unanimous genes in the 38-gene analysis are TNFSF15 and FLT3LG.
These full-cohort feature-selection outputs are exploratory and are not
out-of-sample performance estimates.

## External-cohort results

- **GSE62928 biopsy cohort (4 EPS, 4 combined controls):** 4/8 fixed-panel
  genes were direction-concordant; all eight selected probes were below the
  preregistered operational mean-expression cutoff of 5.0. The continuous
  composite AUC was 0.875 (bootstrap 95% CI 0.500–1.000), but this score is
  based on below-cutoff probes and is highly uncertain; the operational cutoff
  does not establish that the score is noise-dominated. The fixed
  preregistered verdict was Not supported.
- **GSE248762, primary effluent contrast (LV n=10, SV n=6):** 7/7 detected
  genes were direction-concordant; SERPINA10 was not detected. EBI3 and ECM1
  passed BH q<0.05; ADAM19 was nominally significant (P=0.0268) but not after
  BH adjustment. The fixed eight-gene composite and its null were not
  estimable. The literal preregistered verdict remains Not supported, not
  because the direction-count criterion failed, but because the composite
  evidence was unavailable and neither partial-support branch applied.
- **GSE130888, secondary effluent contrast (long-term n=4, short-term n=6):**
  6/8 genes were detected and 5/6 detected genes were concordant. TNFSF15 and
  SERPINA10 were not detected; the full-panel composite/null were not
  estimable. Its frozen cohort-level label is Partially supportive because
  exactly five detected genes were concordant. Patient independence from
  GSE248762 remains unknown, so this is not confirmed independent replication.
- **Post hoc / exploratory detected-gene scores (not the preregistered tests):**
  the 7-gene GSE248762 score had AUC 0.7833 (bootstrap 95% CI 0.5000–1.0000),
  Mann–Whitney P=0.0727, and matched-null AUC empirical P=0.0460. The
  corrected 6-gene GSE130888 score had AUC 0.7083 (bootstrap 95% CI
  0.2917–1.0000), Mann–Whitney P=0.3524, and matched-null AUC empirical
  P=0.0519. The earlier AUC 0.8056 (29/36 pairwise orderings), P=0.1063, and
  matched-null P=0.0330 were erroneous: the calculation treated the three
  normal controls as non-LV cases. The corrected LV-versus-SV AUC is 17/24
  pairwise orderings and excludes all three normal controls from both the
  observed comparison and matched-null AUC. These detected-gene scores are
  post hoc and do not change either preregistered verdict.
- **GSE92455:** expression analysis has not been performed and it is not
  included in any validation verdict.

The preregistered tiers are intentionally preserved; no threshold or verdict
was changed after seeing expression results. The GSE248762 tier behavior is
non-monotonic: 7/7 concordant with an unestimable composite is Not supported,
whereas exactly 5 concordant is Partially supportive. This should be made
more coherent in a future preregistration, not retroactively revised here.

## ECM1, CV, and remaining limitations

In GSE248762 the ECM1 all-cell log2FC was 5.8893; all leave-one-patient-out
estimates remained positive (5.0152–6.0993). This rules out dependence on one
patient alone, but not composition or technical effects. In the most
influential sample, GSM7919584, total pseudobulk ECM1 count was 115,870 among
5,732 QC cells, while only five cells were classified as fibroblasts (three
ECM1-positive). Of the marker-annotated cell-type pseudobulk, 76,200 counts
were mesothelial and 38,545 neutrophil; fibroblasts were excluded from this
cell-type pseudobulk because there were fewer than 20. The composition-adjusted
all-cell estimate was 2.9879 (95% CI -1.5642 to 7.5399, P=0.1472), and the
fibroblast-restricted estimate was 0.5751 (95% CI -0.9553 to 2.1056,
P=0.3892). These results do not establish a fibroblast-intrinsic effect or
prove ambient RNA contamination; the unusually concentrated counts require
raw-data/annotation QC before biological interpretation. No sample was
removed.

Repeated 5-fold CV (10 repeats) gave mean AUC/C-index 0.7700 (bootstrap
95% interval 0.7408–0.7946). The folds repeated DEG, Matrisome, and
three-model selection but not WGCNA; this does not eliminate all selection
optimism. Small sample sizes, unknown cohort overlap, single-cell dropout,
marker-based cell annotations, absent covariates, and wide intervals remain
material limitations. See the referenced reports and saved result tables
under `results/revision/` for full methods and results.

## Panel detectability and panel v2 sensitivity (added 2026-10-11)

This new section reports an assay-detectability audit and rule-derived
panel-v2 sensitivity analysis. It does not revise the preregistered analyses,
thresholds, outputs, or verdicts. The v2 panel was defined using existing
cohort data before any new test data and is not preregistered; all v2 scores
below are post hoc, exploratory, and not validation claims.

### Detectability

The operational GSE62928 rule was a primary-probe mean log2 RMA expression
greater than 5.0; **0/8** panel genes passed. This is an operational cutoff,
not a biological limit of detection. In the effluent cohorts, the rule was a
raw pseudobulk count of at least 10 in at least half (ceiling) of samples in
each of the LV and SV groups: **7/8** genes passed in GSE248762 and **6/8** in
GSE130888. SERPINA10 was annotated in all three cohorts. It had two raw counts
total across 16 GSE248762 samples (mean 0.125; 0/16 samples at count >=10),
and zero raw counts across all 10 GSE130888 contrast samples. Thus it was
below the effluent expression rule in GSE248762 and had zero counts in
GSE130888, rather than being absent from annotation.

### Discovery evidence and mechanical panel selection

`results/revision/panel_v2/panel_discovery_evidence.csv` records each panel
gene's discovery logFC, nominal P, BH-adjusted P across the saved
11,741-gene universe, probe ID, and BH q<0.05 flag. The rule table applies
(a) detected in at least 2/3 validation cohorts and (b) direction concordant
with discovery in every detected cohort; (c) discovery BH q<0.05 is recorded
as a flag only, not a filter. Rule (b) is marked not applicable when no
cohort detects the gene. Applying (a) AND (b), without manual selection,
produced **5 genes: FLT3LG, EBI3, LTB, ADAM19, CST7**.

### Original verdicts versus post-hoc panel-v2 scores

| Cohort | Original 8-gene verdict (unchanged) | Original AUC (95% CI), MW P, matched-null P | V2 exploratory AUC (95% CI), MW P, matched-null P |
|---|---|---|---|
| GSE62928 | Not supported | 0.875 (0.500–1.000), 0.1143, 0.015 | 1.000 (1.000–1.000), 0.02857, not computed |
| GSE248762 | Not supported | Full-panel composite not estimable | 0.6667 (0.3833–0.9333), 0.3132, 0.1209 |
| GSE130888 | Partially supportive | Full-panel composite not estimable | 0.6667 (0.2500–1.000), 0.4762, 0.06194 |

The GSE62928 v2 score is especially limited: all five v2 genes are below the
operational expression cutoff, and its matched-null P was not computed. Its
perfect point AUC and degenerate interval must not be interpreted as robust
validation. Effluent v2 intervals use 10,000 stratified bootstrap resamples
(seed 42); effluent matched-null P values use 1,000 detected-gene,
expression-decile-matched sets. GSE130888 comparisons exclude normal controls
and use LV versus SV only. The existing detected-gene exploratory scores are
reported separately in `panel_v2_score_sensitivity.csv`; they are not the
original fixed-panel test.

See `panel_detectability_audit.csv`, `panel_detectability_audit.md`,
`panel_v2_rule_results.csv`, `panel_v2.json`, `panel_v2_report.md`, and
`panel_v2_sample_scores.csv` in `results/revision/panel_v2/` for the complete
saved tables, methods, sample-level scores, and output provenance.

### First post-freeze GSE92455 panel-v2 test (added 2026-10-11)

The panel-v2 definition was frozen in commit
`26da2cfc081b64be841686f8ee2cf4e6e1014a82` before GSE92455 expression was
analyzed. The supplementary epithelioid-versus-non-epithelioid HPMC contrast
included 9 versus 8 patient labels. The five-gene panel score was estimable:
**AUC 0.3056** (10,000-patient-bootstrap 95% CI **0.0833–0.5972**) and exact
two-sided Mann–Whitney **P=0.1996**, with higher score fixed as epithelioid;
that direction was not predeclared and the AUC was not flipped. This is not a
long-vintage versus short-vintage test and does not count as a validation
verdict. The preregistered GSE92455 analysis remains per-gene only.

On the selected probe per gene, limma results (Epithelioid minus
Non-epithelioid, submitted GEO value scale; BH across the five panel genes)
were significant at q<0.05 for EBI3 (estimate -0.3068, q=0.02444) and ADAM19
(estimate -0.9785, q=0.02444). No matched-gene null was computed for this
distinct phenotype contrast. One consistently oriented array per patient
was analyzed; reciprocal technical arrays and pooled-omentum reference
samples were excluded. See
`results/revision/panel_v2/GSE92455_20261011/` for the saved report, per-gene
results, probe selection, sample scores, and analysis manifest.

### GSE248762–GSE130888 patient-overlap check (added 2026-10-11)

The saved GEO metadata contain 16 GSE248762 sample records and 13 GSE130888
sample records; none of their GSM sample IDs are shared. Neither metadata
table provides patient identifiers or a cross-cohort donor crosswalk.
Therefore, zero shared GSM IDs does **not** establish zero patient overlap:
patient-level independence remains **unknown** and requires a crosswalk or
author confirmation.
