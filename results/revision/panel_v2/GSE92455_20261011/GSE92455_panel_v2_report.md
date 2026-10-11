# GSE92455: first panel-v2 test

**Analysis date:** 2026-10-11
**Panel freeze commit:** `26da2cfc081b64be841686f8ee2cf4e6e1014a82`.
**Label:** post-freeze, exploratory panel-v2 test; not preregistered.

## Cohort and processing

The tested contrast is the supplementary phenotype comparison of epithelioid versus non-epithelioid effluent-derived HPMCs. It is not a vintage contrast and cannot establish long-vintage/short-vintage concordance. GPL6480 ex-vivo data were used; GPL6848 in-vitro treatment arrays and pooled-omentum reference samples were excluded.

GEO metadata identified 17 distinct patient labels (9 epithelioid and 8 non-epithelioid). One array per patient was used, selecting the consistently oriented layout with the effluent sample in Cy3/channel 1 and pooled reference in Cy5/channel 2; reciprocal technical arrays were excluded so technical replicates were not treated as independent patients.

The series-matrix processing metadata states background subtraction and LOESS/LOWESS plus quantile normalization. The GEO record does not explicitly identify the logarithm base for the submitted matrix values. No additional normalization or channel sign reversal was applied; all analyzed arrays have the same channel layout. Limma estimates are reported on the submitted GEO value scale.

## Fixed panel results

The panel was read from the committed `panel_v2.json` before expression analysis. The fixed five-gene composite was estimable.

AUC (higher score = epithelioid): 0.3056 (bootstrap 95% CI 0.0833–0.5972); exact two-sided Wilcoxon rank-sum (no ties), P=0.1996. Bootstrap resampling was stratified by patient group (10,000 resamples; seed 42). The phenotype direction was not predeclared, so the AUC was not flipped to maximize performance.

| Gene | Selected probe | Limma estimate (Epithelioid − Non-epithelioid) | 95% CI | P | BH q (panel) |
|---|---|---:|---:|---:|---:|
| FLT3LG | A_23_P78742 | 0.0758 | -0.1810–0.3327 | 0.5444 | 0.6805 |
| EBI3 | A_23_P119478 | -0.3068 | -0.5306–-0.0829 | 0.0098 | 0.0244 |
| LTB | A_23_P93348 | 0.1339 | -0.2116–0.4794 | 0.4277 | 0.6805 |
| ADAM19 | A_23_P374082 | -0.9785 | -1.6310–-0.3261 | 0.0054 | 0.0244 |
| CST7 | A_23_P68601 | 0.0364 | -0.1874–0.2602 | 0.7375 | 0.7375 |

No random-gene matched-null P was calculated. This distinct phenotype contrast is not covered by the frozen vintage-specific matched-null rule. The original preregistered GSE92455 analysis remains per-gene only; this additional composite is a new exploratory panel-v2 analysis and does not replace or alter it.

## Reproducibility

Saved outputs in this directory include the panel-v2 sample scores, patient-level metadata, outcome-blind probe selection, per-gene limma results, machine-readable score summary, and input/package manifest. The script refuses to overwrite existing non-empty outputs.
