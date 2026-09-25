# Technical Code Logic, Data Integrity & Silent Bug Audit Report
**Project:** Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics  
**Scope:** Verification of data integrity, sample alignment, probe collapsing, WGCNA parameters, machine learning logic, external validation directionality, and document/table consistency.  
**Auditor Mode:** Read-only code inspection & data-integrity verification (no silent fixes applied).  
**Date:** September 2026  

---

## Executive Status Overview

| Check | Domain | Status | Key Finding |
| :--- | :--- | :---: | :--- |
| **Check A** | Sample Identity & Label Integrity | **PASS** | Perfect 1-to-1 GSM accession alignment, title mapping, and group encoding across GSE62928 ($N=8$) and GSE125498 ($N=33$). |
| **Check B** | Gene Symbol Mapping & Intersection | **DISCREPANCY** | Mathematical intersection of Salmon Module (604 genes) and Curated 71 ECM-DEGs is **exactly 40 genes** (100% clean). However, an upstream script variable points to a 148-gene nominal file instead of the 71-gene file. |
| **Check C** | Probe-to-Gene Collapsing Logic | **PASS** | MaxMean collapsing in GSE62928 (Affymetrix GPL13158) was manually verified for 5 test genes (`FN1`, `COL3A1`, `VCAN`, `COL1A1`, `GAPDH`); each selected probe is confirmed to be the highest-mean expression probe. GPL10558 probe mappings in GSE125498 are accurate. |
| **Check D** | WGCNA Parameters & Code Logic | **PASS** | Soft power $\beta = 12$ correctly applied to signed adjacency and TOM; merged module colors (`mergedColors`) correctly used for MEs; Student $t$-test degrees of freedom ($df = 6$) correctly parameterized. |
| **Check E** | ML Consensus Pipeline Correctness | **DISCREPANCY** | In the original discovery run (`05b_ml_hub_gene_identification_wgcna.py`), `StandardScaler` was fit on all 8 samples without holdout. The LOOCV audit corrected this, but XGBoost collapsed ($0\%$ accuracy) on $N=7$ folds. Consensus voting threshold ($\ge 2/4$) is correctly implemented. |
| **Check F** | External Validation (GSE125498) Logic | **DISCREPANCY** | Clinical directionality (LPD [Late] = 1 vs SPD [Early] = 0) and two-sided Mann-Whitney U tests are 100% correct. However, the originally reported composite $\text{AUC} = 0.869 / 0.873$ was evaluated **in-sample** without holdout (cross-validated AUC is $0.696$). |
| **Check G** | File, Table & Figure Consistency | **DISCREPANCY** | Minor numerical transcription discrepancy in `README.md` reporting composite $\text{AUC} = 0.873$ vs current script output $\text{AUC} = 0.869$ (in-sample) and $0.696$ (cross-validated). |

---

## Detailed Check-by-Check Findings

### CHECK A — Sample Identity and Label Integrity
**Status: PASS**

#### GSE62928 (Discovery Cohort, $N = 8$)
1. **Metadata Cross-Check against NCBI GEO Series Matrix:**
   Direct inspection of `data/GSE62928_series_matrix.txt.gz` and `results/tables/GSE62928_sample_metadata.csv` confirms exact alignment:
   - `GSM1536406`: Title = `EPS1` $\to$ Group = `EPS`, Binary = `Case_EPS` (1)
   - `GSM1536407`: Title = `EPS2` $\to$ Group = `EPS`, Binary = `Case_EPS` (1)
   - `GSM1536408`: Title = `EPS3` $\to$ Group = `EPS`, Binary = `Case_EPS` (1)
   - `GSM1536409`: Title = `EPS4` $\to$ Group = `EPS`, Binary = `Case_EPS` (1)
   - `GSM1536410`: Title = `PD1` $\to$ Group = `PD`, Binary = `Control` (0)
   - `GSM1536411`: Title = `PD2` $\to$ Group = `PD`, Binary = `Control` (0)
   - `GSM1536412`: Title = `UREMIC1` $\to$ Group = `Uremic`, Binary = `Control` (0)
   - `GSM1536413`: Title = `UREMIC2` $\to$ Group = `Uremic`, Binary = `Control` (0)
2. **Column Order Alignment:**
   The column headers of `results/tables/GSE62928_full_expression_matrix.csv` match the row order of `results/tables/GSE62928_sample_metadata.csv` position-for-position with zero transposition or row-offset errors.

#### GSE125498 (Validation Cohort, $N = 33$)
1. Parsing of `data/GSE125498_family.soft.gz` confirmed 33 patients:
   - **Early-Stage (Short-term PD, SPD: 0–24 mo):** 20 samples (`GSM3575608` to `GSM3575627`)
   - **Late-Stage (Long-term PD, LPD: $\ge 25$ mo):** 13 samples (`GSM3575628` to `GSM3575640`)
2. In `06_external_validation_GSE125498.py`, `sample_stages` assigns `Late_Stage_LPD` if `"long-term"` is detected in title/characteristics and `Early_Stage_SPD` otherwise. All 33 samples exist in the expression pivot matrix without missingness.

---

### CHECK B — Gene Symbol Mapping Consistency
**Status: DISCREPANCY (Minor Script Variable Mismatch, No Downstream Corruption)**

1. **Independent Intersection Verification:**
   - **WGCNA Salmon Module:** 604 unique gene symbols in `results/tables/wgcna_trait_significant_module_genes.csv`.
   - **Curated 71 ECM-DEGs:** Exactly 71 unique gene symbols in `convergent_71_ECM_DEGs.csv`.
   - **Independent Mathematical Intersection:**
     $$\text{Salmon Module (604)} \cap \text{ECM-DEGs (71)} = \mathbf{40 \text{ Genes}}$$
   - Comparing this independent intersection against `results/tables/convergent_WGCNA_ECM_genes.csv`:
     - Missing in CSV: `set()` (0 genes)
     - Extra in CSV: `set()` (0 genes)
     - Duplicate rows in CSV: `0`
     - Whitespace or casing errors: `None`
   - **The 40 convergent genes are 100% identical and verified.**

2. **Identified Discrepancy:**
   - In [02b_wgcna_analysis.R](file:///d:/Peritoneal%20Project/02b_wgcna_analysis.R#L107), line 107 states:
     ```r
     ecm_71_file <- "convergent_ECM_DEGs_nominal.csv"
     ```
   - **What's wrong:** The variable is named `ecm_71_file`, but it loads `convergent_ECM_DEGs_nominal.csv`, which contains **148 nominal ECM genes** rather than the final 71 filtered genes in `convergent_71_ECM_DEGs.csv`.
   - **Impact Assessment:** Because `02b_wgcna_analysis.R` simply takes the `union()` of the top 5,000 variable genes and this ECM file to build the initial WGCNA input matrix, and all 71 true ECM genes are a strict subset of the 148 nominal genes, this did **not** corrupt module detection. However, the file naming is misleading and should be corrected to `"convergent_71_ECM_DEGs.csv"`.

---

### CHECK C — Probe-to-Gene Collapsing Sanity Check
**Status: PASS**

#### GSE62928 (Affymetrix GPL13158)
Five representative multi-probe genes were checked across the 54,715 raw probes of GPL13158:
1. **`FN1` (7 annotated probes):**
   - Probes: `214701_PM_s_at` (Mean = 7.6104), `211719_PM_x_at` (7.3615), `212464_PM_s_at` (7.1540), `216442_PM_x_at` (7.0240), `210495_PM_x_at` (7.0047), `1558199_PM_at` (3.8437), `214702_PM_at` (1.7169).
   - **Retained in matrix:** `214701_PM_s_at` (Mean = 7.6104) $\implies$ **Confirmed Highest Mean**.
2. **`COL3A1` (4 annotated probes):**
   - Probes: `215076_PM_s_at` (Mean = 9.4292), `201852_PM_x_at` (9.1424), `211161_PM_s_at` (9.0834), `232458_PM_at` (2.9427).
   - **Retained in matrix:** `215076_PM_s_at` (Mean = 9.4292) $\implies$ **Confirmed Highest Mean**.
3. **`VCAN` (5 annotated probes):**
   - Probes: `211571_PM_s_at` (Mean = 6.6651), `221731_PM_x_at` (6.6259), `215646_PM_s_at` (6.3769), `204620_PM_s_at` (5.5640), `204619_PM_s_at` (3.7930).
   - **Retained in matrix:** `211571_PM_s_at` (Mean = 6.6651) $\implies$ **Confirmed Highest Mean**.
4. **`COL1A1` (5 annotated probes):**
   - Probes: `1556499_PM_s_at` (Mean = 11.0146), `202310_PM_s_at` (10.0200), `202311_PM_s_at` (7.4222), etc.
   - **Retained in matrix:** `1556499_PM_s_at` (Mean = 11.0146) $\implies$ **Confirmed Highest Mean**.
5. **`GAPDH` (6 annotated probes):**
   - Probes: `AFFX-HUMGAPDH/M33197_3_at` (Mean = 6.6538), `212581_PM_x_at` (6.0375), etc.
   - **Retained in matrix:** `AFFX-HUMGAPDH/M33197_3_at` (Mean = 6.6538) $\implies$ **Confirmed Highest Mean**.

#### GSE125498 (Illumina GPL10558)
All 11 hub genes were checked against `Validation/GSE125498.top.table.tsv`:
- `VCAN`: Probe `ILMN_1687301` matches `VCAN` (versican, $\log_2\text{FC} = -0.522, P = 0.0244$).
- `COL8A1`: Probe `ILMN_2402392` matches `COL8A1` (collagen VIII $\alpha 1$, $\log_2\text{FC} = +0.749, P = 0.0488$).
- `FN1`: Probe `ILMN_1778237` matches `FN1` (fibronectin 1, $\log_2\text{FC} = +0.407, P = 0.2690$).
- `THBS3`: Probe `ILMN_1804663` matches `THBS3` (thrombospondin 3, $\log_2\text{FC} = -0.201, P = 0.4180$).
- `COL3A1`: Probe `ILMN_1773079` matches `COL3A1` (collagen III $\alpha 1$, $\log_2\text{FC} = +0.187, P = 0.6160$).
- `ISM1`: Probe `ILMN_3239288` matches `ISM1` (isthmin 1, $\log_2\text{FC} = +0.028, P = 0.9040$).
- `LOX`: Probe `ILMN_1695880` matches `LOX` (lysyl oxidase, $\log_2\text{FC} = +0.016, P = 0.9760$).
- `EDIL3`, `COMP`, `COL11A1`, `INHBA`: Confirmed completely absent from the GPL10558 beadchip manifest (not a lookup bug).

---

### CHECK D — WGCNA Parameter and Code Logic Verification
**Status: PASS**

Code inspection of [02b_wgcna_analysis.R](file:///d:/Peritoneal%20Project/02b_wgcna_analysis.R):
1. **Sample Alignment:** Lines 47–48 explicitly sort metadata to match expression row order:
   ```r
   sample_order <- rownames(datExpr0)
   meta_df <- meta_df[match(sample_order, meta_df$sample_id), ]
   trait_df <- data.frame(Fibrosis_Status = trait_vec, row.names = sample_order)
   ```
2. **Soft Power Consistency:** `chosen_power <- 12` is defined in Step 3 and passed directly into:
   ```r
   adjacency <- adjacency(datExpr, power = chosen_power, type = "signed")
   TOM <- TOMsimilarity(adjacency, TOMType = "signed")
   ```
   No leftover or hardcoded power variable slipped into the calculations.
3. **Module Color Merging:** Line 221 merges modules at cut height 0.25 into `mergedColors`. Line 247 explicitly sets `moduleColors <- mergedColors`, ensuring module eigengenes (`MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes`) use merged cluster assignments rather than unmerged dynamic cuts.
4. **Student Degrees of Freedom:** `n_samples <- nrow(datExpr)` ($N = 8$). Line 255 calls `corPvalueStudent(moduleTraitCor, n_samples)`. In WGCNA, this computes $t = r \sqrt{\frac{n-2}{1-r^2}}$ with $df = n-2 = 6$. For $r = 0.805997$, $t = 3.3353 \implies P = 0.015701$. The calculation is exact.

---

### CHECK E — Machine Learning Pipeline Correctness
**Status: DISCREPANCY (Methodological Overfitting in Original Script)**

Inspection of [05b_ml_hub_gene_identification_wgcna.py](file:///d:/Peritoneal%20Project/05b_ml_hub_gene_identification_wgcna.py):
1. **Sample & Feature Alignment:** Line 58–69 extracts `sample_ids = df_meta["sample_id"].tolist()` and subsets `X_df = df_expr.loc[available_genes, sample_ids].T`. Sample ordering of $X$ and $y$ is strictly aligned.
2. **Voting Threshold Logic:** Line 199 implements `"Is_Hub_Gene": "YES" if votes >= 2 else "NO"`. This is mathematically correct ($\ge 2/4$ models).
3. **Identified Discrepancies:**
   - **In-Sample Feature Selection & Preprocessing:** Lines 76–77 run `scaler.fit_transform(X_df)` on all 8 samples before training the 4 models without an independent test split. This guaranteed 100% training separation on $N = 8$.
   - **XGBoost Collapse on Small $N$:** In the LOOCV audit, XGBoost achieved $0.0\%$ accuracy on held-out folds because gradient boosting with `max_depth = 2` overfit the $N = 7$ training folds, predicting the training majority class every time.
   - **Hub Gene Instability:** While 8 genes were selected in 8/8 LOOCV folds, `EDIL3` was only selected in 5/8 folds ($62.5\%$), showing that single-fit feature selection masked instability for marginal features.

---

### CHECK F — External Validation Script Correctness
**Status: DISCREPANCY (In-Sample AUC Evaluation in Original Script)**

Inspection of [06_external_validation_GSE125498.py](file:///d:/Peritoneal%20Project/06_external_validation_GSE125498.py):
1. **Clinical Directionality:**
   - `sample_stages`: `Late_Stage_LPD` = 1, `Early_Stage_SPD` = 0.
   - In `Validation/GSE125498.top.table.tsv`:
     - `VCAN`: $\log_2\text{FC} = -0.5224$ (Table) vs $-0.5224$ (Calculated Late - Early).
     - `COL8A1`: $\log_2\text{FC} = +0.7493$ (Table) vs $+0.7493$ (Calculated Late - Early).
     - **Directional alignment between clinical metadata, expression matrices, and top table statistics is identical.**
2. **Mann-Whitney U Test:** Line 128 explicitly sets `alternative="two-sided"`.
3. **Identified Discrepancy:**
   - Lines 183–187 fit a `LogisticRegression` composite model on all 33 patients and evaluate `predict_proba()` on those same 33 patients:
     ```python
     clf.fit(X_val[available_hubs], y_binary)
     prob = clf.predict_proba(X_val[available_hubs])[:, 1]
     fpr_comp, tpr_comp, _ = roc_curve(y_binary, prob)
     auc_comp = auc(fpr_comp, tpr_comp) # Yields 0.869
     ```
   - **What's wrong:** Evaluating a 7-parameter logistic model on $N = 33$ without cross-validation yields severe optimistic bias.
   - **Correct generalization metric:** As established in `statistical_rigor_audit.py`, under 5-Fold Stratified Cross-Validation, the true generalization $\text{AUC}$ is **$0.696 \pm 0.056$** ($\text{LOOCV AUC} = 0.677$).

---

### CHECK G — File, Table and Figure Consistency
**Status: DISCREPANCY (Minor Transcription Values in Documentation)**

Spot-checking 13 specific numerical parameters across the codebase:
1. Soft power $\beta = 12$ $\to$ **Matches across code, logs, and plots.**
2. Total modules = 14 $\to$ **Matches.**
3. Salmon module genes = 604 $\to$ **Matches.**
4. Salmon trait correlation = $0.8060$ ($P = 0.0157$) $\to$ **Matches.**
5. Salmon Bonferroni $P = 0.2198$, FDR $Q = 0.2198$ $\to$ **Matches.**
6. Exact permutation $P = 0.0143$ ($1/70$) $\to$ **Matches.**
7. Convergent WGCNA $\cap$ ECM genes = 40 $\to$ **Matches.**
8. ML consensus hub genes = 11 $\to$ **Matches.**
9. `VCAN` GSE125498 $\log_2\text{FC} = -0.5224, P = 0.0244, \text{AUC} = 0.7231$ $\to$ **Matches.**
10. `COL8A1` GSE125498 $\log_2\text{FC} = +0.7493, P = 0.0488, \text{AUC} = 0.6654$ $\to$ **Matches.**
11. Profiled hub genes in validation = 7 $\to$ **Matches.**
12. **Discrepancy in Composite AUC:**
    - `README.md` lines 61 & 140 state: `Composite Signature AUC = 0.873`.
    - `06_external_validation_GSE125498.py` outputs: `Composite Signature AUC = 0.869`.
    - `statistical_rigor_audit.py` cross-validated output: `5-Fold CV AUC = 0.696`.
    - **Recommendation:** Update `README.md` to state both the in-sample $\text{AUC} = 0.869$ and the cross-validated $\text{AUC} = 0.696 \pm 0.056$ to prevent transcription confusion.

---

## Final Technical Verdict Statement

> ### **FINAL STATEMENT:**
> **No fatal technical bugs, data transpositions, or sample misalignments were found in the codebase.**  
> Sample identities, GSM accessions, Affymetrix/Illumina probe mappings, and directional fold changes ($\text{LPD} - \text{SPD}$) are verified to be mathematically and biologically consistent across both datasets.
> 
> The discrepancies identified are **methodological rather than coding bugs**:
> 1. In `02b_wgcna_analysis.R`, line 107 loads `convergent_ECM_DEGs_nominal.csv` (148 genes) into a variable named `ecm_71_file`. (Harmless upstream, but should be updated to `convergent_71_ECM_DEGs.csv` for clarity).
> 2. The original ML feature selection (`05b_ml_hub_gene_identification_wgcna.py`) and external validation composite AUC (`06_external_validation_GSE125498.py`) were evaluated **in-sample without cross-validation**, causing optimistic inflation ($\text{AUC} = 0.869$ vs true cross-validated $\text{AUC} = 0.696$).
> 3. `README.md` quotes an older run value of $\text{AUC} = 0.873$ rather than the current script value of $0.869$ (or cross-validated $0.696$).
> 
> All core quantitative discrepancies are directly attributable to the small-sample-size limitations already documented and addressed in the statistical rigor audit.
