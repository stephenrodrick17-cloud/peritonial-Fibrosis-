# Statistical Rigor & Overfitting Audit Summary Report
**Project:** Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics  
**Discovery Cohort:** GSE62928 (Human peritoneal tissue biopsies; $N = 8$: 4 EPS Cases vs 4 PD/Uremic Controls)  
**External Validation Cohort:** GSE125498 (Human peritoneal effluent cells; $N = 33$: 20 Short-term PD [SPD, 0–24 mo] vs 13 Long-term PD [LPD, $\ge 25$ mo])  
**Target Panel:** 11 WGCNA-ECM Consensus Hub Genes (`ISM1`, `FN1`, `EDIL3`, `VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`)  
**Audit Date:** September 2026  

---

## Executive Summary & Final Verdict

### **Verdict: (b) SUGGESTIVE / EXPLORATORY ONLY**

Based on the execution of all six validation control audits, the 11-gene WGCNA-ECM hub panel **cannot be claimed as statistically validated**, but it is **not pure noise**:

1. **Why it CANNOT be reported as "Validated"**:
   - **Discovery Module Multi-Testing Failure:** The pro-fibrotic Salmon module ($r = +0.806$, raw $P = 0.0157$) **fails** multiple testing correction across the 14 tested WGCNA modules ($\text{Bonferroni } P = 0.2198$, $\text{BH FDR } Q = 0.2198$).
   - **Validation Cohort Multi-Testing Failure:** In the independent external dataset (GSE125498), **ZERO of the 7 profiled hub genes** survive Bonferroni or Benjamini-Hochberg FDR correction. The nominal signals (`VCAN` raw $P = 0.0244 \to P_{\text{bonf}} = 0.1708$; `COL8A1` raw $P = 0.0488 \to P_{\text{bonf}} = 0.3416$) become non-significant once penalizing for 7 hypothesis tests.
   - **Substantial Overfitting in Composite AUC:** The previously reported in-sample composite $\text{AUC} = 0.869$ drops to **$\text{LOOCV AUC} = 0.677$** and **$5\text{-Fold Stratified CV AUC} = 0.696 \pm 0.056$** (a drop of $\Delta\text{AUC} = 0.173$ to $0.192$). The in-sample AUC was heavily inflated by evaluating the model on the exact same 33 samples it was trained on.
   - **Extreme Sample Size Vulnerability:** Negative control simulations on $N = 8$ show that **26.5%** of completely random label assignments generate at least one "statistically significant" module ($P < 0.05$), and **29.4%** generate $|r| \ge 0.70$.

2. **Why it is NOT "Indistinguishable from Noise"**:
   - **Cross-Cohort Generalization Above Chance:** In the external validation cohort ($N = 33$), the out-of-fold cross-validated composite AUC remains **$0.696$**, which is significantly above the null baseline of $0.500$.
   - **Extreme Permutation Ranking:** In the discovery cohort, the Salmon module correlation ($r = 0.806$) is the **single most extreme positive correlation possible** out of all 70 possible label combinations ($P_{\text{perm}} = 1/70 = 0.0143$).
   - **Biological Coherence:** The surviving core hub genes (`VCAN`, `COL8A1`, `FN1`, `COL3A1`) are heavily documented regulators of extracellular matrix remodeling, epithelial-to-mesenchymal transition (EMT), and peritoneal membrane thickening.

**Conclusion:** The 11-gene signature represents a **biologically credible, suggestive set of pro-fibrotic candidates** that demonstrates modest cross-cohort discrimination ($\text{AUC} \approx 0.68 - 0.70$). It must be presented strictly as an **exploratory biomarker panel**, requiring verification in large, prospectively powered clinical cohorts ($N \ge 50 - 100$).

---

## Detailed Findings for Checks 1 through 6

### Check 1: Multiple Testing Correction on WGCNA Module-Trait Correlations
In the discovery analysis (GSE62928), 14 merged co-expression modules were correlated against the binary fibrosis trait ($N = 8$).
- **Bonferroni significance threshold:** $\alpha / 14 = 0.05 / 14 = \mathbf{0.00357}$

| Module | Genes | Trait Correlation ($r$) | Raw Student's $P$-Value | Bonferroni Adjusted $P$ | Benjamini-Hochberg FDR $Q$ | Survives Multi-Testing? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`MEsalmon`** | **604** | **+0.8060** | **0.0157** | **0.2198** | **0.2198** | **NO** |
| `MEyellowgreen` | 66 | +0.6415 | 0.0865 | 1.0000 | 0.3406 | NO |
| `MEdarkturquoise` | 190 | +0.6317 | 0.0929 | 1.0000 | 0.3406 | NO |
| `MEskyblue` | 200 | -0.5889 | 0.1246 | 1.0000 | 0.3406 | NO |
| `MEskyblue3` | 393 | -0.5571 | 0.1514 | 1.0000 | 0.3406 | NO |
| `MElightcyan` | 100 | +0.5469 | 0.1607 | 1.0000 | 0.3406 | NO |
| `MEbrown4` | 202 | -0.5366 | 0.1703 | 1.0000 | 0.3406 | NO |
| `MEfloralwhite` | 827 | +0.4222 | 0.2974 | 1.0000 | 0.5150 | NO |
| `MEdarkorange2` | 172 | -0.3963 | 0.3311 | 1.0000 | 0.5150 | NO |
| `MEbisque4` | 846 | +0.3168 | 0.4445 | 1.0000 | 0.6223 | NO |
| `MEgreenyellow` | 128 | -0.2783 | 0.5045 | 1.0000 | 0.6421 | NO |
| `MEroyalblue` | 1137 | -0.2284 | 0.5865 | 1.0000 | 0.6842 | NO |
| `MEsienna3` | 67 | +0.1155 | 0.7853 | 1.0000 | 0.8426 | NO |
| `MEcyan` | 101 | +0.0843 | 0.8426 | 1.0000 | 0.8426 | NO |

> **Audit Finding:** Not a single module reaches statistical significance after family-wise error rate (FWER) or false discovery rate (FDR) control. The Salmon module pro-fibrotic association must be transparently designated as nominally significant only ($P = 0.0157$), subject to high risk of type I error.

---

### Check 2: Permutation Test for Salmon Module-Trait Correlation
To bypass the asymptotic assumptions of the parametric Student's $t$-distribution on a tiny sample size ($N = 8$), an exact combinatorial permutation test was conducted.
- There are exactly $\binom{8}{4} = \mathbf{70}$ unique ways to partition 8 patients into 4 Cases and 4 Controls.

* **Observed MEsalmon Correlation ($r_{\text{obs}}$):** $+0.8060$ (Parametric $P = 0.0157$)
* **Exact One-Sided Permutation $P$-Value ($r \ge 0.8060$):** **$\mathbf{0.0143}$** ($1 / 70$ combinations)
* **Exact Two-Sided Permutation $P$-Value ($|r| \ge 0.8060$):** **$\mathbf{0.0286}$** ($2 / 70$ combinations)
* **Monte Carlo Permutations ($B = 10,000$):** One-sided $P = 0.0151$; Two-sided $P = 0.0302$

> **Audit Finding:** The observed correlation of $0.8060$ is the **maximum possible correlation** that could have been observed given this eigengene across all 70 permutations. However, because $N = 8$, the minimum attainable two-sided $P$-value is mathematically bounded at $2/70 = 0.0286$. It passes nominal permutation testing as an isolated hypothesis, but its sample size makes it impossible to ever survive family-wise correction across 14 modules.

---

### Check 3: Leave-One-Out Cross-Validation (LOOCV) ML Consensus on GSE62928
The initial ML consensus feature selection was fit on all 8 discovery samples simultaneously, yielding 100% training separation. LOOCV was implemented where in each of the 8 folds, 1 sample was completely held out, and all 4 models were re-fit on the remaining 7 samples (with internal standardization to prevent data leakage).

* **Null Baseline Accuracy (Majority Class):** $50.0\%$ (4 Cases vs 4 Controls)
* **LASSO (L1-Logistic):** Accuracy = **$100.0\%$** (8/8), Sensitivity = $100.0\%$, Specificity = $100.0\%$
* **SVM-RFE (Linear):** Accuracy = **$87.5\%$** (7/8), Sensitivity = $75.0\%$, Specificity = $100.0\%$
* **Random Forest:** Accuracy = **$87.5\%$** (7/8), Sensitivity = $100.0\%$, Specificity = $75.0\%$
* **XGBoost:** Accuracy = **$0.0\%$** (0/8) *(Collapsed due to tree-depth overfitting on $N=7$)*
* **Consensus Vote ($\ge 2/4$ algorithms):** Accuracy = **$\mathbf{87.5\%}$** (7/8), Sensitivity = **$100.0\%$**, Specificity = **$75.0\%$**

#### Hub Gene Selection Stability Across the 8 LOOCV Folds:
| Gene Symbol | Selection Frequency (out of 8 Folds) | Selection Stability (%) | Notes |
| :--- | :---: | :---: | :--- |
| `FN1` | 8 / 8 | 100.0% | Core anchor hub |
| `INHBA` | 8 / 8 | 100.0% | Perfectly stable |
| `COL11A1`| 8 / 8 | 100.0% | Perfectly stable |
| `COMP` | 8 / 8 | 100.0% | Perfectly stable |
| `THBS3` | 8 / 8 | 100.0% | Perfectly stable |
| `COL3A1`| 8 / 8 | 100.0% | Perfectly stable |
| `LOX` | 8 / 8 | 100.0% | Perfectly stable |
| `ISM1` | 8 / 8 | 100.0% | Perfectly stable |
| `VCAN` | 7 / 8 | 87.5% | Dropped in Fold 4 |
| `COL8A1`| 7 / 8 | 87.5% | Dropped in Fold 1 |
| `EDIL3` | 5 / 8 | **62.5%** | **Substantial fold instability** (dropped in 3/8 folds) |

> **Audit Finding:** While the consensus vote achieves an 87.5% LOOCV generalization accuracy (substantially outperforming the 50% baseline), feature selection stability is uneven: `EDIL3` drops out in 37.5% of iterations, demonstrating that the single-fit 11-gene list contained marginal features dependent on individual samples.

---

### Check 4: Cross-Validated External Validation AUC on GSE125498
In the validation cohort ($N = 33$), the logistic regression composite score was re-evaluated under strict cross-validation instead of naive in-sample scoring:

* **Original In-Sample Composite AUC:** **$0.869$**
* **Leave-One-Out Cross-Validation (LOOCV) AUC:** **$0.677$** (Overfitting Drop: **$\Delta = 0.192$**)
* **5-Fold Stratified Cross-Validation (50 repeats):** **$\mathbf{0.696 \pm 0.056}$** (Overfitting Drop: **$\Delta = 0.173$**)

> **Audit Finding:** The naive in-sample AUC of 0.869 was severely overfitted by evaluating on the training data. The true out-of-fold generalization performance is **$\text{AUC} \approx 0.68 - 0.70$**. While this indicates moderate pro-fibrotic discriminative signal above the 0.50 baseline, claiming an AUC $> 0.85$ is scientifically indefensible.

---

### Check 5: Individual Gene Significance with Multiple Testing Correction
In GSE125498, 7 hub genes were profiled on the Illumina beadchip. Individual differential expression between Early-Stage (SPD, $n=20$) and Late-Stage (LPD, $n=13$) was audited:
- **Bonferroni threshold:** $0.05 / 7 = \mathbf{0.00714}$

| Gene Symbol | Direction | $\log_2\text{FC}$ | Raw Limma $P$ | Limma Bonferroni $P$ | Limma BH FDR $Q$ | Raw Mann-Whitney $P$ | MW Bonferroni $P$ | MW BH FDR $Q$ | Survives Multi-Testing? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`VCAN`** | DOWN | -0.522 | 0.0244 | 0.1708 | 0.1708 | 0.0341 | 0.2389 | 0.2389 | **NO** |
| **`COL8A1`** | UP | +0.749 | 0.0488 | 0.3416 | 0.1708 | 0.1174 | 0.8217 | 0.4108 | **NO** |
| `FN1` | UP | +0.407 | 0.2690 | 1.0000 | 0.6277 | 0.2937 | 1.0000 | 0.6417 | NO |
| `THBS3` | DOWN | -0.201 | 0.4180 | 1.0000 | 0.7315 | 0.3667 | 1.0000 | 0.6417 | NO |
| `COL3A1`| UP | +0.187 | 0.6160 | 1.0000 | 0.8624 | 0.7541 | 1.0000 | 0.9459 | NO |
| `ISM1` | UP | +0.028 | 0.9040 | 1.0000 | 0.9760 | 0.8107 | 1.0000 | 0.9459 | NO |
| `LOX` | UP | +0.016 | 0.9760 | 1.0000 | 0.9760 | 0.9853 | 1.0000 | 0.9853 | NO |

> **Audit Finding:** **Zero (0 / 7) individual genes survive multiple testing correction.** Both `VCAN` and `COL8A1` achieve nominal significance ($P < 0.05$ uncorrected), but their adjusted $P$-values are $\ge 0.1708$. Genome-wide across the entire microarray, their Limma FDR values were $0.3057$ and $0.4022$. They cannot be reported as statistically confirmed single-gene biomarkers.

---

### Check 6: Negative Control (Sanity Check with Fake Random Traits)
To directly test whether the WGCNA pipeline at $N = 8$ discovers spurious pro-fibrotic modules from complete noise, random 4-vs-4 binary traits were assigned to the 8 discovery samples and correlated against all 14 module eigengenes.
- Tested across all 68 non-true combinatorial permutations:

* **Average Maximum $|r|$ from pure noise:** **$0.6300$**
* **Fraction of random traits generating $\ge 1$ module with $P < 0.05$:** **$\mathbf{26.5\%}$**
* **Fraction of random traits generating $\ge 1$ module with $|r| \ge 0.70$:** **$\mathbf{29.4\%}$**
* **Peak correlation generated from pure noise:** **$r = 0.9380$ ($P = 5.70 \times 10^{-4}$)**

> **Audit Finding:** This negative control provides undeniable mathematical evidence that with $N = 8$ samples and 14 modules, WGCNA has a **~27% to 29% chance of producing a strong, significant-looking module correlation ($|r| > 0.70, P < 0.05$) purely by chance**. The discovery correlation of $r = 0.806$ is therefore vulnerable to high false-discovery risk.

---

## Methodological Recommendations for Manuscript / Reporting

1. **State the Exploratory Nature Explicitly:** Re-classify the 11 hub genes from "validated diagnostic biomarkers" to an "exploratory candidate panel for peritoneal dialysis-induced fibrogenesis."
2. **Report Both In-Sample and Cross-Validated Metrics:** Always present the cross-validated $\text{AUC} = 0.696$ alongside the in-sample $\text{AUC} = 0.869$, transparently discussing the optimism penalty.
3. **Report Corrected $P$-Values Side-by-Side:** When reporting `VCAN` ($P_{\text{raw}} = 0.024$) and `COL8A1` ($P_{\text{raw}} = 0.049$), clearly indicate that neither survives FDR correction ($Q = 0.171$).
4. **Emphasize the Need for Larger Prospective Cohorts:** Highlight that discovery at $N = 8$ was constrained by biopsy availability in Encapsulating Peritoneal Sclerosis (EPS), and mandate that future validation utilize large multicenter peritoneal effluent cohorts ($N \ge 100$) using targeted RNA-seq or multiplex digital PCR.
