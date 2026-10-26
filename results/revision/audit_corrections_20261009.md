# Audit corrections and remaining limitations

**Audit date:** 2026-10-09  
**Scope:** User-raised numerical, statistical, interpretation, and file-date
concerns. No expression data were changed; no sample was removed from the
preregistered analyses; and no threshold or preregistered verdict was changed.

| Concern | Audit result and correction |
|---|---|
| Effluent verdicts appear non-monotonic | Confirmed in the frozen tier wording. GSE248762 has 7/7 detected genes concordant, so it meets the direction-count component for Supportive, but SERPINA10 is undetected and therefore the required fixed eight-gene composite and matched-null AUC are not estimable. Its partial tier requires exactly 5 concordant genes or an estimable AUC of 0.65–0.80 without null significance; neither branch applies. The literal verdict remains **Not supported**. GSE130888 is **Partially supportive** because exactly five detected genes are concordant. The GSE248762 script previously short-circuited to Not supported whenever the composite was unavailable, incorrectly skipping the exact-five partial branch; that implementation defect is fixed. It does not change either saved cohort verdict. The non-monotonic decision rule remains a limitation to address in a future preregistration, not retroactively. |
| Primary/secondary cohort language and overall verdict | Corrected summary wording: GSE248762 is primary; GSE130888 is secondary. Their cohort labels are reported separately. Under the frozen overall rule neither contrast is Supportive, so the overall label is Not supported. GSE130888's partial cohort label is not replaced by the overall label. The two datasets' independence remains unknown, so no confirmed independent-replication claim is made. |
| ML selection counts could not be verified | Verified from both saved JSON summaries and their matching gene tables. For the 624-gene pool: LASSO 24, SVM-RFE 10, RF 10, eight >=2/3 consensus, zero 3/3. For the 38-gene ECM/WGCNA pool: 11, 10, 10, eight >=2/3 consensus, two 3/3 (TNFSF15 and FLT3LG). The old four-model/XGBoost and 11-hub summary is superseded, not a result of these requested three-model analyses. |
| Hypergeometric and permutation P-values differ | No arithmetic discrepancy found. The hypergeometric test uses M=11,741 measured genes, K=428 measured Matrisome genes, n=496 up-DEGs, and k=31 observed up-ECM DEGs. Expected overlap is 496*428/11,741=18.09; fold enrichment is 1.7145 and the saved one-sided P is 0.0024028. The 10,000-draw permutation instead shuffles sample labels and re-runs the DEG threshold on measured Matrisome genes, giving Monte Carlo empirical P=0.0306969. These are different null models. Wording calling the latter “exact” has been corrected to “Monte Carlo empirical P.” |
| ECM1 log2FC=5.89 may reflect low counts, composition, or ambient RNA | The saved real-count audits confirm log2FC=5.8893 and that all leave-one-patient-out estimates remain positive (5.0152–6.0993), so the result is not attributable to a single patient alone. However, the cell-level distribution raises a legitimate unresolved QC concern: GSM7919584 has 115,870 ECM1 raw pseudobulk counts among 5,732 QC cells, but only five marker-classified fibroblasts (three expressing ECM1). Its marker-classified cell-type pseudobulk has 76,200 mesothelial and 38,545 neutrophil ECM1 counts; fibroblasts were omitted from that table because n=5 (<20). After adjustment for monocyte/macrophage and T/NK proportions, log2FC=2.9879 (95% CI -1.5642 to 7.5399, P=0.1472); the fibroblast-restricted estimate is 0.5751 (95% CI -0.9553 to 2.1056, P=0.3892). These observations do not prove ambient RNA or establish a fibroblast-intrinsic effect. They are now documented as a material caveat; no sample was removed and no causal interpretation is made. |
| WGCNA leakage, FDR, immune genes, and biopsy AUC | Preserved explicitly: WGCNA was not rerun inside CV folds, so 0.7700 nested AUC may retain upstream selection optimism. The 135/813 q<0.05 count is within nominal DEGs; 143 is the full-universe q<0.05 count across 11,741 tested genes. The global-null expectation of 587 nominal P<0.05 tests is not an estimate of false discoveries among the 813 DEGs and precedes the fold-change filter. Several panel genes are immune-associated, and composition remains a confounder. In GSE62928 all eight primary probes are below the operational expression cutoff, but that cutoff alone neither establishes a noise-dominated AUC=0.875 nor provides evidence of reliable discrimination; the 95% bootstrap interval [0.500, 1.000] indicates substantial uncertainty. |
| GSE130888 analysis folder date | Verified against the system date and file modification times: outputs were written on 2026-10-09, and the output directory name was corrected to `analysis_20261009`. The script's default output path was updated. Saved cohort contents were not re-run or altered. |

## Updated summary locations

- [README.md](../../README.md) now separates the two effluent-cohort labels,
  explains the composite limitation, distinguishes the enrichment nulls, and
  corrects the cross-validation and below-cutoff interpretation.
- [FINAL_VERIFIED_NUMBERS.md](../../FINAL_VERIFIED_NUMBERS.md) replaces the
  legacy “final” summary that contained the synthetic validation and other
  withdrawn claims.
- The historical `nested_cv_summary.json` now carries the same percentile
  bootstrap interval [0.7408, 0.7946] as the authoritative `_taskE_v2`
  output; the README points to the authoritative result.
- The legacy `manuscript_numbers.json` is marked non-authoritative. The
  uncommitted legacy-runner replacement has unresolved provenance and is
  excluded from the proposed commits pending separate review.
- [GSE248762 analysis report](effluent_validation/GSE248762/analysis_20261009/gse248762_analysis_report.md)
  and [deviation log](effluent_validation/GSE248762/DEVIATIONS.md) explain the
  literal verdict and implementation correction.
- [GSE130888 report](effluent_validation/GSE130888/analysis_20261009/gse130888_analysis_report.md)
  now resides under the correctly dated folder.
- [GSE62928 report](gse62928_validation_8genes/gse62928_validation_report.md)
  distinguishes an operational cutoff from proof of noise.

## Still unresolved

The ECM1 raw-count concentration warrants verification against the original
10x files, feature annotation, and cell labels; current saved summaries cannot
determine whether it reflects biology, annotation, ambient RNA, or a technical
artifact. Patient overlap between GSE248762 and GSE130888 is unknown. GSE92455
has not been analyzed. The nested CV does not include WGCNA in-fold. Therefore
the repository should not be described as error-free or as establishing
independent validation. No commit was made during this audit.

## Continued audit: verification and post hoc / exploratory analyses

### Verification of the prior working-tree edits

`git diff` was inspected for the requested tracked files. The README and
`FINAL_VERIFIED_NUMBERS.md` contain the corrected validation wording; the
GSE248762 Scanpy diff contains the corrected verdict branches and rationale;
the GSE248762 deviations log records the post-run fixes. The GSE130888 Scanpy
script and deviations log, and this audit report, are untracked files, so Git
has no HEAD version against which to show their contents. The current files
were read directly instead. The tracked legacy-runner diff is unexpectedly
broad (about 1,082 added and 72 removed lines), far beyond the described
output-path/wording adjustments; it should be reviewed separately before it
is staged. This audit did not rewrite or revert that diff.

The stale-string search found no residual outdated dual-primary-cohort
wording, the former folder-date label, claims of exactness for the Monte
Carlo test, or unsupported characterization of the biopsy score. The only
remaining literal match for the searched primary-cohort wording is a sentence
in the frozen preregistration at line 152 describing its two statistical
contrasts; it does not call both datasets primary. The previous interval values in `CHANGELOG_REVISION.md` lines 13–14,
`nested_cv_summary.json` lines 14 and 20, and `nested_cv_results.csv` lines
4–5 were corrected to the authoritative [0.7408, 0.7946]. Remaining
substring matches are unrelated data or score values, not CV confidence
intervals.

The initial phrase audit also located low-expression/noise wording in
`README.md` line 92, `FINAL_VERIFIED_NUMBERS.md` line 64,
`scripts/revision_gse62928_validation_8genes_report.py` line 139,
`gse62928_validation_report.md` line 74, and this report's first table at
line 15. Each was revised to state the inferential limitation without
overstating what the operational expression cutoff means. The Monte Carlo
permutation description in the legacy runner was likewise corrected. The prior output-folder date mention at this report's original
line 16 was rewritten as a dated correction rather than a stale active path.

### Post hoc / exploratory results

The separate output folder
`exploratory_posthoc_20261009/` contains scripts, configuration, tables, and
this report. The detected-gene composite uses seven genes in GSE248762 and six
in GSE130888. A hand-calculation audit found that the earlier GSE130888
helper classified every non-LV sample as a control, including the three
NORMAL samples. Thus its reported AUC used 4 LV versus 9 non-LV samples
(29/36 = 0.8056), despite the saved summary reporting n=4 vs n=6. The script
now restricts the observed and matched-null AUCs to LV and SV; the summary
lists the three excluded NORMAL samples. Pairwise recalculation gives 17
correct orderings among 24 LV–SV pairs, AUC 17/24 = 0.7083 (increment 1/24
without ties; half-tie increment 1/48). The bootstrap 95% CI is
0.2917–1.0000, Mann–Whitney P=0.3524, and matched-null AUC empirical
P=0.0519.

GSE248762 has only LV and SV samples: 47/60 correct pairwise orderings give
AUC 0.7833, unchanged; its bootstrap 95% CI is 0.5000–1.0000,
Mann–Whitney P=0.0727, and matched-null AUC empirical P=0.0460. Both are
labelled **post hoc detected-gene composite, not the preregistered test** and
do not change GSE248762's **Not supported** or GSE130888's **Partially
supportive** verdict.

The full ECM1 raw-cell audit re-read all 16 GSE248762 sample matrices using
the same QC and annotation rules and confirmed each saved QC cell total.
GSM7919584 had 115,870 ECM1 counts in 5,732 cells (1,331 ECM1-positive);
the top 1% of cells contained 30.5% of its counts and the top 10% contained
97.7%. The largest assigned cell-type totals were 76,200 counts in
mesothelial-labeled cells (614 positive) and 38,545 in neutrophil-labeled
cells (269 positive). The same sample showed mixed expression of epithelial
and neutrophil markers across these labels. Since empty-droplet/background
profiles were unavailable, these patterns neither prove nor rule out ambient
RNA, and they do not establish a cell-intrinsic biological effect.

In post hoc edgeR sensitivities, ECM1 remained positive after excluding
GSM7919584 (log2FC 5.0152, approximate 95% CI [1.4166, 8.6138], P=0.00839)
and in mesothelial-only pseudobulk (log2FC 7.2997, approximate 95% CI
[3.1948, 11.4047], P=0.00139). The sensitivity script uses `filterByExpr`
then forces ECM1 into the fit; it uses unadjusted `~ 0 + groups`, TMM
normalization, robust edgeR quasi-likelihood fitting, and approximate CIs
from inversion of the one-d.f. QL F statistic. A fresh R process reproduced
the saved CSV byte-for-byte (SHA-256
`42D55D19106136E4BBD149C107134016468095049715C727E0F259517E1C5353`).
Mesothelial pseudobulk raw ECM1 counts were SV 78, 115, 80, 68, 331, 52 and
LV 34, 76200, 17694, 126, 905, 112, 224, 53, 276, 111. These are raw
patient-level sums before TMM normalization. GSM7919584's marker comparison
and the all-sample cell-concentration table are included in the exploratory
outputs. These are exploratory contrasts; the preregistered full-sample
analysis and inclusion remain unchanged. No ambient contamination or
biological cause is inferred from these data.

The GEO sample metadata has no patient/donor crosswalk for either effluent
cohort; titles contain only within-cohort group and replicate labels. The
family SOFT files list raw 10x supplementary matrices, but those filenames are
GSM/group keyed and contain no shared donor identifier. The available
downloaded data therefore cannot resolve cross-cohort overlap; independence
remains unknown.
This metadata-only audit is saved in
`exploratory_posthoc_20261009/cohort_independence_audit.csv` and `.json`.

The existing README and this report now point to the separate exploratory
outputs. History identifies the tracked runner at HEAD as a short Stage 4/5
version and shows only three historical committed updates (initial runner,
BH-sync addition, and portable-manifest addition); it does not explain the
current 1,082-line worktree replacement. The replacement is an eight-phase
legacy pipeline, but Git history cannot establish who made it or whether it
was intentional. A separate untracked near-copy also exists. Both current
copies have a LEGACY warning and guarded output destination, but these runner
changes are excluded from the proposed commits pending provenance review. The
README Phase 6 description now matches the inspected code: it loads saved
real GSE62928 outputs and rejects a hub-panel mismatch. No GSE92455 analysis,
staging, or commit was performed.
