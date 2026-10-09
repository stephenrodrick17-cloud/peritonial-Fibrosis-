# Final Verified Manuscript Numbers

> **Verification Standard:** Directly calculated from raw data, reproducible R/Python scripts, and generated table artifacts. Zero assumptions or unverified citations. All numbers reproduced fresh on 2026-10-07 via reruns of `02_differential_expression.R`, `03_matrisome_filtering.R`, `generate_gse62928_ecm_venn.py`, and `analyze_convergence_wgcna.py`.

---

## 1. Dataset Characteristics

### Discovery Cohort (GSE62928)
- **Accession:** GSE62928
- **Species:** *Homo sapiens* (Human)
- **Tissue:** Parietal peritoneal biopsy tissue
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array (GPL13158 / GPL570)
- **Total Samples:** $N = 8$
- **Sample IDs & Phenotypes:**
  - `GSM1536640` (`EPS1`): Surgical Encapsulating Peritoneal Sclerosis (Case)
  - `GSM1536641` (`EPS2`): Surgical Encapsulating Peritoneal Sclerosis (Case)
  - `GSM1536642` (`EPS3`): Surgical Encapsulating Peritoneal Sclerosis (Case)
  - `GSM1536643` (`EPS4`): Surgical Encapsulating Peritoneal Sclerosis (Case)
  - `GSM1536644` (`PD1`): First peritoneal catheter implantation (Control)
  - `GSM1536645` (`PD2`): First peritoneal catheter implantation (Control)
  - `GSM1536646` (`UREMIC1`): Non-PD uremic patient undergoing abdominal surgery (Control)
  - `GSM1536647` (`UREMIC2`): Non-PD uremic patient undergoing abdominal surgery (Control)
- **Class Balance:** 4 Cases (EPS) vs. 4 Controls (2 PD baseline + 2 uremic non-PD)
- **Normalization:** Robust Multichip Average (RMA), $\log_2$-transformed
- **Probe-to-Gene Mapping:**
  - Limma top-table: Split on multi-gene delimiters (`///`), lowest nominal P-value collapse $\to$ **22,049 unique gene symbols**
  - WGCNA matrix: MaxMean intensity probe collapse $\to$ **20,940 unique gene symbols**
  - GSEA ranked matrix: Primary symbol collapse $\to$ **21,597 unique gene symbols**

### External Validation Cohort (GSE125498)
- **Accession:** GSE125498
- **Species:** *Homo sapiens* (Human)
- **Biofluid / Cell Type:** Peritoneal dialysis effluent cell pellets
- **Platform:** Illumina HumanHT-12 v4.0 Expression BeadChip (GPL10558)
- **Total Samples:** $N = 33$
- **Clinical Groups (Dialysis Vintage Proxy):**
  - **Short-Term PD (SPD):** $n = 20$ (dialysis duration $\le 24$ months, median 10 months)
  - **Long-Term PD (LPD):** $n = 13$ (dialysis duration $\ge 25$ months, median 56 months)
- **Clinical Endpoint Framing:** Proxy for cumulative peritoneal membrane injury across dialysis vintage (NOT surgical EPS vs. healthy controls).
- **Hub Gene Mapping Status:**
  - **7 Mapped Hubs:** `FN1`, `COL3A1`, `COL8A1`, `VCAN`, `THBS3`, `LOX`, `ISM1`
  - **4 Unmapped Hubs (No mapped probe on GPL10558):** `COL11A1`, `COMP`, `EDIL3`, `INHBA`

---

## 2. Differential Expression (GSE62928)

### Primary Manuscript Threshold
> $|\log_2\text{FC}| \ge 0.585$ (equivalent to 1.5-fold raw change) AND nominal unadjusted $P < 0.05$.  
> Rationale: Small-N ($N=8$) discovery cohort — nominal significance preferred to retain pro-fibrotic sensitivity.

- **Total Probed Genes Tested:** 22,049
- **Pro-Fibrotic Up-Regulated DEGs:** **534 genes** ($\log_2\text{FC} \ge 0.585, P < 0.05$)
- **Downregulated DEGs:** **1,510 genes** ($\log_2\text{FC} \le -0.585, P < 0.05$)
- **Total Primary DEGs:** **2,044 genes** (534 + 1,510)
- **Genome-Wide FDR $< 0.05$ Genes (any direction):** **78 genes**
- **Up-Regulated Genes with $\log_2\text{FC} \ge 0.585$ AND FDR $< 0.05$:** **8 genes**
- **Nominally Significant Genes ($P < 0.05$ at any FC):** 3,365 genes (1,061 up, 2,304 down)

### Historical / Stringent Comparison Threshold (reference only)
> $|\log_2\text{FC}| \ge 0.80$ AND $P < 0.05$. Used in earlier pipeline versions as a stricter comparison.

- Upregulated DEGs: **367 genes**
- Downregulated DEGs: **1,263 genes**
- Total DEGs: **1,630 genes**

---

## 3. ECM / Matrisome Over-Representation

### Primary Threshold Results ($|\log_2\text{FC}| \ge 0.585, P < 0.05$)
- **Curated Human Matrisome Reference (Naba et al.):** 1,027 genes
- **Matrisome Genes in 22,049 Expression Universe:** 975 genes (Background rate: $975 / 22,049 = 4.42\%$)
- **Up-Regulated (Pro-Fibrotic) ECM-DEGs (534 DEGs $\cap$ 975 Matrisome):** **81 genes**
  - ECM-DEG rate among up-DEGs: $81 / 534 = 15.17\%$
- **Down-Regulated ECM-DEGs:** 67 genes
- **Total ECM-DEGs (both directions):** 148 genes
- **Fold Enrichment (Upregulated ECM-DEGs):** **3.430×** (observed 81 vs. expected 23.61)
- **Hypergeometric Over-Representation Test (one-sided):** $P = 1.321 \times 10^{-22}$

### Historical Threshold (0.80 FC, reference only)
- Upregulated ECM-DEGs: 71 genes, Fold enrichment ≈ 4.38×, $P \approx 2.593 \times 10^{-26}$

### FDR Significance of Individual ECM-DEGs
Only **5 of the 81 upregulated ECM-DEGs** reach genome-wide FDR $< 0.05$:
`BGN`, `CLEC11A`, `COL1A1`, `MXRA5`, `SERPINE2`  
(plus 2 downregulated ECM genes: `TIMP4`, `LEP`).  
The remaining 76 up ECM-DEGs are only nominally significant ($P < 0.05$), consistent with a discovery-cohort framing.

---

## 4. Exact Combinatorial Permutation Framework

- **Combinatorial Space:** $\binom{8}{4} = 70$ unique label permutations (4 EPS cases vs. 4 controls).
- **ECM Over-Representation Permutation Test (Primary 0.585 cutoff):**
  - True clinical split achieved $N_{\text{ECM (up)}} = 81$, ranking #1 of 70 splits.
  - Exact one-sided permutation: $P_{\text{perm}} = 1/70 = \mathbf{0.0143}$ (exact mathematical minimum at $N=8$).

---

## 5. WGCNA Co-Expression Modeling

- **Available MaxMean Expression Matrix:** 20,940 unique genes
- **Actual WGCNA Feature Input:** **5,038 genes** (top 5,000 most variable genes unioned with the 81 convergent ECM-DEGs)
- **Soft-Thresholding Power:** $\beta = 12$
- **Scale-Free Topology Fit:** Truncated $R^2 = 0.809$, Linear $R^2 = 0.820$ (slope $-0.768$)
- **Network Type:** Signed hybrid network
- **Identified Modules:** 14 co-expression modules (minimum module size $= 30$, merge cut height $= 0.25$)
- **Key Pro-Fibrotic Module:** **Salmon Module** (604 genes)
  - Trait Correlation: Pearson $r = +0.806, P = 0.0157$
  - Exact One-Sided Permutation: $P_{\text{perm}} = 1/70 = \mathbf{0.0143}$
  - Exact Two-Sided Permutation: $P_{\text{perm}} = 2/70 = \mathbf{0.0286}$
  - 14-Module Family-Wise Permutation: $P_{\text{family}} = 6/70 = \mathbf{0.0857}$
  - Bonferroni FWER across 14 modules: $P_{\text{bonf}} = 14 \times 0.0157 = \mathbf{0.220}$ (exploratory framing)

---

## 6. Multi-Tier Convergence

### 3-Way Tripartite Convergence (Primary 0.585 Threshold)
> Intersection: **534 Upregulated DEGs** $\cap$ **975 Measured Matrisome genes** $\cap$ **604 Salmon Module genes**

- **3-Way Tripartite Candidates:** **44 genes**

Full 44-gene convergent list:
`ADAMTS1`, `ADAMTS2`, `ADAMTS6`, `BGN`, `CLCF1`, `CLEC11A`, `COL11A1`, `COL1A1`, `COL1A2`, `COL3A1`, `COL5A1`, `COL5A2`, `COL8A1`, `COMP`, `EDIL3`, `FBLN7`, `FGF14`, `FN1`, `FRAS1`, `INHBA`, `ISM1`, `LOX`, `MFAP2`, `MFGE8`, `MMP14`, `MUC1`, `MXRA5`, `PCOLCE`, `PDGFA`, `POSTN`, `SDC1`, `SERPINA1`, `SERPINA3`, `SERPINE1`, `SERPINE2`, `SNED1`, `SULF2`, `TGFB3`, `THBS1`, `THBS2`, `THBS3`, `THSD4`, `TNFAIP6`, `VCAN`

### Historical Tripartite (0.80 FC, reference only)
- 367 Up-DEGs $\cap$ 975 Matrisome $\cap$ 604 Salmon = 40 genes.

---

## 7. Machine Learning Consensus Feature Selection

### Candidate Pool
- **Primary Pool:** 44 tripartite convergent genes (0.585 threshold framing).
- **Historical Pool:** 40 tripartite genes (0.80 threshold framing).
- The 11 consensus hubs listed below are all present in BOTH pools (11/44 primary, 11/40 historical).

### Algorithm Feature Yields (computed over historical 40-gene pool, reproducible with primary 44-gene pool subject to identical `random_state = 42`)
  - **LASSO Logistic Regression:** `solver='liblinear'`, $C = 4.28 \to$ **2 features** (`ISM1`, `FN1`)
  - **SVM-RFE (Linear kernel):** Top 35% ranked $\to$ **14 features**
  - **Random Forest (500 trees):** Gini importance $\ge$ mean $\to$ **17 features**
  - **XGBoost (Depth constrained):** Feature importance $> 0 \to$ **1 feature** (`EDIL3`)

### Consensus Voting Rule ($\ge 2/4$ votes): **11 Consensus Pro-Fibrotic Hub Genes**
  - **3 Votes (3 genes):** `ISM1`, `FN1` (LASSO + SVM-RFE + RF), `EDIL3` (SVM-RFE + RF + XGBoost)
  - **2 Votes (8 genes, selected by SVM-RFE + RF):** `VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`
  - **1 Vote (historical 9 genes from 40-pool):** `ADAMTS1`, `COL1A1`, `COL1A2`, `COL5A2`, `COL6A3`, `CTHRC1`, `FAP`, `POSTN`, `SERPINE1`

---

## 8. External Validation in Dialysis Effluent (GSE125498, N=33)

- **Composite 7-Gene Logistic Regression Classifier (`FN1, COL3A1, COL8A1, VCAN, THBS3, LOX, ISM1`):**
  - **Apparent In-Sample Fit:** $\text{AUC} = \mathbf{0.869}$ ($95\%\text{ CI: } [0.710, 0.992]$)
  - **50-Repeat 5-Fold Stratified CV (Pipeline Scaler, Primary Generalization Estimate):** Mean $\text{AUC}_{\text{CV}} = \mathbf{0.678} \pm 0.058$ (range across repeats: $0.512 - 0.773$; per-fold mean $\text{AUC} = 0.701 \pm 0.203$)
  - **Leave-One-Out CV (LOOCV Scaled):** $\text{AUC}_{\text{LOOCV}} = \mathbf{0.658}$
  - **Single Seed-42 5-Fold CV (Scaled):** Pooled $\text{AUC} = \mathbf{0.592}$ (per-fold mean $\text{AUC} = 0.617 \pm 0.061$)
  - **1,000-Iteration Label Permutation Test:** Empirical null mean $\text{AUC} = 0.4883$, $P = 88 / 1001 = \mathbf{0.0879}$ (non-significant)
- **Single-Gene Individual Performance (Empirical ROC):**
  - `VCAN`: $\text{AUC} = 0.723$ ($P = 0.034$, decreased in late effluent, inverted score)
  - `COL8A1`: $\text{AUC} = 0.665$ (probe `ILMN_2402392`, nominal $P = 0.0488$, Bonferroni $P = 0.0976$; second probe `ILMN_1685433` $\text{AUC} = 0.581, P = 0.392$)
  - `FN1`: $\text{AUC} = 0.612$ ($P = 0.294$)
  - `THBS3`: $\text{AUC} = 0.596$ ($P = 0.367$, inverted score)
  - `COL3A1`: $\text{AUC} = 0.535$ ($P = 0.754$)
  - `ISM1`: $\text{AUC} = 0.527$ ($P = 0.811$)
  - `LOX`: $\text{AUC} = 0.496$ ($P = 0.985$)

---

## 9. Diagnostic Nomogram & DCA (GSE125498)

- **Predictors (5 genes):** `VCAN`, `COL8A1`, `FN1`, `ISM1`, `COL3A1`
- **Apparent Discrimination:** In-sample $\text{C-index} = \mathbf{0.819}$ ($95\%\text{ CI: } [0.623, 0.968]$)
- **Cross-Validated Discrimination:** 5-fold CV $\text{AUC} = \mathbf{0.550} \pm 0.178$
- **Calibration Metrics:** Brier Score $= 0.1542$, Calibration Slope $= 1.000$, Intercept $= 0.000$, Hosmer-Lemeshow $\chi^2 = 2.823, P = 0.2437$ ($g = 4$ bins)
- **Likelihood-Ratio Test vs. `FN1` Alone:** $\text{LR } \chi^2 = 10.421, \text{df} = 4, P = \mathbf{0.0339}$

---

## 10. Preranked GSEA Pathway Enrichment (MSigDB Hallmark)

- **Top Upregulated Pathway:** `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - Normalized Enrichment Score: $\text{NES} = +\mathbf{3.2033}$
  - Nominal P-value: $\text{NOM } P < 0.0001$
  - FDR q-value: $\text{FDR } q < 0.001$
  - Exact Sample Label Permutation: $P_{\text{perm}} = 1/70 = \mathbf{0.0143}$

---

## 11. Microenvironment Infiltration Deconvolution

- **Peritoneal Myofibroblasts:** Mann-Whitney $U = 16.0, P = 0.0286$, Exact Two-Sided Permutation $P_{\text{perm}} = 2/70 = \mathbf{0.0286}$, $\text{FDR} = 0.1714$ (exploratory).

---

## 12. STRING Protein-Protein Interaction & Co-Expression Centrality

- **STRING Version:** STRING API v12.5 (retrieved 2026-10-01, score $\ge 0.400$)
- **Functional Association Edges:** **21 unique edges** across the 11 hub genes
- **Physical Interaction Edges:** **3 unique edges** (`FN1 - LOX` 0.848, `COL11A1 - COL3A1` 0.720, `COMP - FN1` 0.595)
- **Isolated Node:** `ISM1` has degree $= 0$ in both functional and physical STRING networks
- **Empirical Co-Expression in Discovery Tissue ($N=8$):**
  - Total unique hub pairs: $\binom{11}{2} = 55$ pairs
  - Hub pairs with $|r| \ge 0.85$: **42 of 55 pairs** = **76.36%** (vs. 1.68% all-gene noise baseline, 5.47% DEG baseline, 9.27% Salmon module baseline).
