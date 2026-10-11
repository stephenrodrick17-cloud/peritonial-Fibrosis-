# Panel v2 rule and score report

**Panel date:** 2026-10-11. **Status:** defined before any new test data.

Panel inclusion is determined mechanically by rule (a) AND rule (b). BH q<0.05 is only a flag (rule c), not a filter. All v2 scores below reuse existing cohorts and are labeled post hoc, exploratory, and not preregistered; they are not new validation claims.

## Rule-derived genes

| Gene | Detected in ≥2/3 | Direction concordant wherever detected | Discovery BH q<0.05 flag | Included |
|---|---|---|---|---|
| TNFSF15 | FAIL | PASS | Yes | No |
| FLT3LG | PASS | PASS | Yes | Yes |
| EBI3 | PASS | PASS | No | Yes |
| LTB | PASS | PASS | Yes | Yes |
| ADAM19 | PASS | PASS | Yes | Yes |
| ECM1 | PASS | FAIL | No | No |
| SERPINA10 | FAIL | NOT APPLICABLE (no detected cohort) | No | No |
| CST7 | PASS | PASS | Yes | Yes |

**Panel v2 (5 genes):** FLT3LG, EBI3, LTB, ADAM19, CST7.

## Original verdicts and v2 exploratory scores

| Cohort | Original 8-gene verdict (unchanged) | Original 8-gene AUC (95% CI) | Original MW P | Original matched-null P | V2 exploratory AUC (95% CI) | V2 MW P | V2 matched-null P |
|---|---|---:|---:|---:|---:|---:|---:|
| GSE62928 | Not supported | 0.875 (0.5–1) | 0.1143 | 0.015 | 1 (1–1) | 0.02857 | not estimable/not computed |
| GSE248762 | Not supported | not estimable | not estimable/not computed | not estimable/not computed | 0.6667 (0.3833–0.9333) | 0.3132 | 0.1209 |
| GSE130888 | Partially supportive | not estimable | not estimable/not computed | not estimable/not computed | 0.6667 (0.25–1) | 0.4762 | 0.06194 |

Original effluent 8-gene composites were not estimable under their frozen detection rule. Existing detected-gene exploratory scores are separately retained in `panel_v2_score_sensitivity.csv`; they are not the original fixed-panel test. The historical GSE130888 detected-gene AUC of 0.8056 included normal controls; this audit uses the corrected saved result of 0.7083 for LV vs SV only. The biopsy v2 genes remain below the operational expression threshold, so its v2 matched-null P is not calculated. For effluent v2 score nulls, 1,000 detected-gene expression-decile matched sets and the plus-one empirical-P formula were used.

The v2 score uses the mean of per-gene z-scored expression across the primary contrast samples; AUC confidence intervals use 10,000 stratified bootstrap resamples (seed 42). GSE130888 normal controls are excluded from the primary LV-vs-SV score and test.

## Integrity

No preregistered result, verdict, threshold, sample inclusion, or file under the previous cohort analysis folders was modified. The panel is rule-derived; no gene was selected by manual judgment.
