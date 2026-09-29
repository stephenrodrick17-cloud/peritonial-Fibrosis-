# High-Dimensional Matrisome Convergence & Machine Learning Transcriptomics in Peritoneal Fibrosis

[![Pipeline Status](https://img.shields.io/badge/Pipeline-Locked%20%26%20Reproducible-success.svg)](#-pipeline-execution-guide)
[![Discovery Dataset](https://img.shields.io/badge/Discovery-GSE62928_%28N%3D8%29-blue.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62928)
[![Validation Dataset](https://img.shields.io/badge/Validation-GSE125498_%28N%3D33%29-indigo.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125498)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## 🔬 Study Overview & Graphical Abstract

Encapsulating Peritoneal Sclerosis (EPS) is the most severe and life-threatening complication of long-term peritoneal dialysis (PD), characterized by extensive fibrocollagenous thickening, neoangiogenesis, and progressive encapsulation of the peritoneal membrane.

This repository provides the complete, deterministic bioinformatics workflow for identifying extracellular matrix (ECM) drivers of peritoneal fibrogenesis. By integrating discovery transcriptomics of human peritoneal tissue (**GSE62928**, $N = 8$: 4 EPS cases vs 4 uremic/PD controls), curated Matrisome masterlists (Naba et al., 1,027 ECM genes), weighted gene co-expression network analysis (WGCNA), four-way consensus machine learning (LASSO, SVM-RFE, Random Forest, XGBoost), external cohort evaluation (**GSE125498**, $N = 33$ peritoneal effluent cell profiles), live protein-protein interaction networking (STRING v12.0), permutation-controlled gene set enrichment analysis (GSEA), and microenvironmental immune deconvolution.

![Study Graphical Abstract](results/figures/graphical_abstract.jpg)
* **Figure 0: Comprehensive Study Workflow.** Discovery in human peritoneal biopsy transcriptomics ($N = 8$), intersection with the human Matrisome, WGCNA co-expression modeling, machine-learning consensus selection, external validation in dialysis effluent ($N = 33$), live STRING protein interaction networks, and permutation-controlled pathway/immune deconvolution.

---

## 📊 Summary of Key Findings

| Domain / Finding | Metric / Statistic | Validation / Permutation Control | Scientific Interpretation |
| :--- | :--- | :--- | :--- |
| **Primary Finding: ECM Over-Representation** | **71 ECM-DEGs** among **367 up-regulated DEGs** ($19.35\%$ vs $4.42\%$ background) | Exact Combinatorial Label Permutation ($N=70$ splits): **$P_{\text{perm}} = 1/70 = \mathbf{0.0143}$** (Rank 1/70) | **Robust Primary Result:** True EPS tissue split produces more ECM-DEGs than any other label split (minimum attainable $P$-value at $N=8$). Hypergeometric $P = 2.59 \times 10^{-26}$ is a secondary descriptive statistic (assumes gene independence). |
| **Secondary Finding: WGCNA Salmon Module** | Module size: 604 genes; Trait correlation: **$r = +0.806$** (Student's $P = 0.0157$) | Permutation: $P_{\text{perm}} = 1/70 = \mathbf{0.0143}$; FWER across 14 modules: **$P_{\text{bonf}} = \mathbf{0.2198}$** | **Suggestive / Exploratory:** The observed $r=+0.806$ is the maximum possible correlation across all 70 splits, but fails family-wise error control. Random label splits yield significant modules $26.5\%$ of the time. |
| **Convergent Candidate Layer** | **40 Pro-Fibrotic ECM Genes** | Strict intersection: 604 Salmon Module $\cap$ 71 ECM-DEGs | Core candidate pool capturing interstitial matrix, basement membrane, and fibrogenic regulators. |
| **Consensus ML Hub Panel** | **11 Hub Genes** (`ISM1`, `FN1`, `EDIL3`, `VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`) | $\ge 2/4$ ML model consensus votes (LASSO, SVM-RFE, RF, XGBoost) | Cross-validated ML feature stability: 8/11 genes selected in 100% of LOOCV folds; `EDIL3` selected in 62.5% of folds. |
| **External Cohort Discrimination (GSE125498)** | Primary 7-Gene Model (all profiled hubs) | In-Sample: **$\text{AUC} = \mathbf{0.869}$**; 50-Repeat 5-Fold Stratified CV: **$\text{AUC} = \mathbf{0.696} \pm 0.056$** | **Moderate Cross-Cohort Generalization ($\text{AUC} \approx 0.70$):** Clear signal above null ($0.50$), but in-sample scores $>0.85$ reflect overfitting ($\Delta\text{AUC} = 0.173$). 0/7 individual genes survive FDR. |
| **Pathway Activation (GSEA)** | Epithelial-Mesenchymal Transition (EMT) | Preranked GSEA: $\text{NES} = +\mathbf{3.115}, \text{FDR} < 10^{-4}$; Permutation: **$P_{\text{perm}} = 1/70 = \mathbf{0.0143}$** | **Focal Pathway Robust:** EMT is the #1 enriched hallmark across all 70 permutations. Global pathway count ($\ge 35$ nominal Hallmarks) is caveated ($P_{\text{perm}} = 0.1286$). |
| **Immune Microenvironment** | Peritoneal Myofibroblast Expansion | Mann-Whitney $U = 16.0, P_{\text{raw}} = 0.0286$; Permutation: **$P_{\text{perm}} = 2/70 = \mathbf{0.0286}$** | **Suggestive Myofibroblast Signal:** Global deconvolution ($31.4\%$ baseline chance of $\ge 1$ significant cell type, family-wise $P=0.3143$) requires prospective single-cell validation. |

---

## 📈 Phase-by-Phase Analytical Pipeline

### Phase 1: Quality Control, Normalization & Discovery Differential Expression
Discovery transcriptomics was profiled from human peritoneal biopsies (**GSE62928**, Affymetrix Human Genome U133 Plus 2.0 Array, $N = 8$: 4 severe Encapsulating Peritoneal Sclerosis cases vs 4 non-EPS uremic/peritoneal dialysis controls). Probe summarization with MaxMean collapsing yielded 20,940 unique genes.

```
Total Probed Genes: 20,940
Discovery DEG Rule: Nominal P < 0.05 and log2FC >= 0.80 (Up-regulated only)
Total Up-regulated DEGs: 367
```

![DEG Volcano & MA Plots](results/figures/DEG_01_volcano_plot.png)
* **Figure 1: Discovery Differential Expression Landscape.** Volcano plot of 20,940 genes highlighting the 367 up-regulated pro-fibrotic DEGs ($P < 0.05, \log_2\text{FC} \ge 0.80$) identified in human peritoneal tissue.

---

### Phase 2: Curated Matrisome Filtering & ECM Enrichment Permutation Test
The 367 up-regulated DEGs were intersected with the Human Matrisome Project reference database (Naba et al., 1,027 curated ECM structural and regulatory genes).

* **Matrisome Overlap:** **71 ECM-DEGs** ($19.35\%$ of DEGs vs $4.42\%$ genome-wide background rate; 4.38-fold enrichment).
* **Primary Evidence (Exact Label Permutation Test):** All $\binom{8}{4} = 70$ combinatorial sample label partitions were evaluated. The true clinical split ranked **#1 out of 70** ($P_{\text{perm}} = 1/70 = \mathbf{0.0143}$, the minimum achievable $p$-value at $N=8$).
* **Secondary Descriptive Metric:** Hypergeometric test $P = 2.5929 \times 10^{-26}$ (reported as supportive descriptive context, assuming gene independence).

![Matrisome Venn Diagram](results/figures/venn_gse62928_pro_fibrotic_ecm_71.png)
* **Figure 2: Human Matrisome Overlap.** Venn diagram displaying the intersection of 367 up-regulated DEGs with the 1,027-gene curated Matrisome, isolating the 71 pro-fibrotic ECM-DEGs.

---

### Phase 3: Weighted Gene Co-Expression Network Analysis (WGCNA)
A signed co-expression network was constructed across all 20,940 genes ($\beta = 12, R^2 = 0.82$, minimum module size $= 30$, merge cut height $= 0.25$), resolving 14 co-expression modules.

* **Salmon Module:** 604 genes, positively correlated with the EPS clinical trait ($r = +0.806, P = 0.0157$).
* **Permutation Significance:** $P_{\text{perm}} = 1/70 = \mathbf{0.0143}$ (Rank 1/70 across permutations).
* **Multiple Testing Caveat:** Bonferroni adjusted $P = \mathbf{0.2198}$ (fails family-wise control across 14 modules). Random label permutations produce a significant module ($P < 0.05$) in **$26.5\%$** of trials.
* **Convergence Candidates:** Intersecting the 604 Salmon module genes with the 71 ECM-DEGs identified **40 Convergent WGCNA-ECM Pro-Fibrotic Candidates**.

![WGCNA Module Heatmap](results/figures/WGCNA_03_module_trait_heatmap.png)
* **Figure 3: Module-Trait Association Heatmap.** Pearson correlation between 14 WGCNA module eigengenes and binary peritoneal fibrosis status ($N = 8$).

![WGCNA Convergence Venn](results/figures/venn_wgcna_convergence.png)
* **Figure 4: 3-Way Convergence Venn Diagram.** Tripartite overlap between GSE62928 up-regulated DEGs ($n=367$), Curated Matrisome ($n=1,027$), and WGCNA Salmon Module ($n=604$), defining the 40 convergent candidate genes.

---

### Phase 4: Consensus Machine Learning Feature Selection (11 Hub Genes)
Four supervised machine learning algorithms were trained on the 40 convergent candidate genes across the discovery cohort:
1. **LASSO (L1-Penalized Logistic Regression):** $\alpha = 0.042$, 10 genes selected.
2. **SVM-RFE (Support Vector Machine Recursive Feature Elimination):** Linear kernel, 11 genes selected.
3. **Random Forest (Mean Decrease Gini):** 500 trees, top 11 features.
4. **XGBoost (Extreme Gradient Boosting):** Depth-constrained gradient boosted trees, top 11 features.

**Consensus Rule:** Candidates selected by $\ge 2$ algorithms were designated consensus hub genes, yielding **11 Consensus Pro-Fibrotic Hub Genes**:
$$\text{Hub Panel: } \mathbf{ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX}$$

```
Consensus ML Votes:
- 4/4 Models: ISM1, FN1, EDIL3, VCAN
- 3/4 Models: COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA
- 2/4 Models: LOX
```

![ML Consensus Barchart](results/figures/WGCNA_ML_01_consensus_votes_barchart.png)
* **Figure 5: Machine Learning Consensus Selection.** Vote aggregation across LASSO, SVM-RFE, Random Forest, and XGBoost identifying the 11 hub genes.

![Hub Gene Heatmap](results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png)
* **Figure 6: Discovery Expression Profiles.** Clustered heatmap of the 11 consensus hub genes in discovery tissue ($N = 8$: 4 EPS vs 4 Control).

---

### Phase 5: External Cohort Cross-Validation (GSE125498, N = 33)
The hub genes were evaluated in an independent clinical dataset (**GSE125498**, Illumina HumanHT-12 v4.0 Expression BeadChip, $N = 33$ human peritoneal dialysis effluent cell samples: 20 short-term PD [SPD, 0–24 months] vs 13 long-term PD [LPD, $\ge 25$ months]).

* **Probe Coverage:** 7 of 11 hub genes were represented on the Illumina array (`FN1`, `COL3A1`, `COL8A1`, `VCAN`, `THBS3`, `LOX`, `ISM1`). Four genes (`COL11A1`, `COMP`, `EDIL3`, `INHBA`) lacked mapped probes.
* **Primary Composite 7-Gene Model:** 
  - In-Sample (training fit): $\text{AUC} = \mathbf{0.869}$ ($95\%\text{ CI: } [0.710, 0.992]$)
  - 50-Repeat 5-Fold Stratified Cross-Validation: $\text{AUC} = \mathbf{0.696} \pm 0.056$ (True Out-of-Fold Generalization)
  - Leave-One-Out Cross-Validation (LOOCV): $\text{AUC} = \mathbf{0.677}$
* **Secondary 5-Gene Sub-Model (Nomogram):** Excluding `THBS3` (downregulated in effluent) and `LOX` (flat in effluent) yields an in-sample $\text{AUC} = \mathbf{0.800 - 0.819}$ and 5-fold cross-validated $\text{AUC} = \mathbf{0.550}$.

![External Validation Boxplots](results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png)
* **Figure 7: External Cohort Expression Boxplots.** Expression of profiled hub genes in peritoneal effluent cells across Early (SPD) vs Late (LPD) stages ($N = 33$).

---

### Master Reconciliation of Discrimination Metrics

| Cohort / Model | Feature(s) / Predictors | $N$ | Validation Method | AUC / C-index | 95% Confidence Interval | Methodological Assessment |
| :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| **GSE125498 (Effluent)** | `VCAN` (Single Gene) | 33 | Empirical ROC | **0.723** ($P=0.034$) | [0.540, 0.880] | Nominal $P<0.05$; fails multi-testing FDR ($Q=0.171$) |
| **GSE125498 (Effluent)** | `COL8A1` (Single Gene) | 33 | Empirical ROC | **0.665** ($P=0.117$) | [0.470, 0.850] | Limma $P=0.049$; fails multi-testing FDR ($Q=0.171$) |
| **GSE125498 (Effluent)** | `FN1` (Single Gene) | 33 | Empirical ROC | **0.612** ($P=0.294$) | [0.410, 0.800] | Non-significant in effluent |
| **GSE125498 (Effluent)** | `THBS3` (Single Gene) | 33 | Empirical ROC | **0.596** ($P=0.367$) | [0.380, 0.790] | Inverted direction in effluent vs biopsy |
| **GSE125498 (Effluent)** | `COL3A1` (Single Gene) | 33 | Empirical ROC | **0.535** ($P=0.754$) | [0.320, 0.740] | Non-significant in effluent |
| **GSE125498 (Effluent)** | `ISM1` (Single Gene) | 33 | Empirical ROC | **0.527** ($P=0.811$) | [0.310, 0.730] | Non-significant in effluent |
| **GSE125498 (Effluent)** | `LOX` (Single Gene) | 33 | Empirical ROC | **0.496** ($P=0.985$) | [0.280, 0.710] | Indistinguishable from chance baseline |
| **GSE125498 (Primary 7-Gene Panel)** | 7 Profiled Hub Genes | 33 | In-Sample Logistic Fit | **0.869** | [0.710, 0.992] | **Optimistic / In-Sample (Subject to overfitting)** |
| **GSE125498 (Primary 7-Gene Panel)** | 7 Profiled Hub Genes | 33 | 50x 5-Fold Stratified CV | **0.696** | [0.640, 0.752] | **Lead Generalization Metric (Realistic performance)** |
| **GSE125498 (Primary 7-Gene Panel)** | 7 Profiled Hub Genes | 33 | Leave-One-Out CV | **0.677** | N/A | Consistent out-of-fold generalization drop ($\Delta=0.192$) |
| **GSE125498 (Secondary 5-Gene Nomogram)** | `VCAN, COL8A1, FN1, ISM1, COL3A1` | 33 | In-Sample Statsmodels `Logit` | **0.819** | [0.623, 0.968] | **In-Sample Nomogram C-index (DCA model fit)** |
| **GSE125498 (Secondary 5-Gene Nomogram)** | `VCAN, COL8A1, FN1, ISM1, COL3A1` | 33 | 5-Fold Stratified CV | **0.550** | [0.322, 0.759] | Out-of-fold generalization drop ($\Delta=0.250$) |

> [!IMPORTANT]
> **Saturated-Model Note:** Multi-gene models evaluated on tiny discovery samples ($p=11 > n=8$) achieve mathematical separation ($\text{AUC} = 1.0$) trivially due to parameter saturation. In accordance with strict statistical standards, discovery multi-gene AUC claims are excluded as uninformative.

![ROC Analysis Multi-Panel](results/figures/Hub_02b_roc_analysis.png)
* **Figure 8: External Validation ROC Curves.** (A) Single-gene ROC curves in GSE125498 ($N=33$). (B) Multi-gene composite ROC comparison demonstrating the gap between optimistic in-sample fits ($\text{AUC} = 0.869$) and cross-validated out-of-fold performance ($\text{AUC} = 0.696$).

---

### Peritoneal Biopsy vs. Effluent Discordance

| Gene Symbol | Discovery Biopsy (GSE62928, $N=8$) | Effluent Cells (GSE125498, $N=33$) | Directional Concordance | Biological & Compartmental Rationale |
| :--- | :--- | :--- | :--- | :--- |
| **`COL8A1`** | $\log_2\text{FC} = +1.11, P = 0.003$ | $\log_2\text{FC} = +0.75, P = 0.049$ | **Concordant (UP)** | Robust vascular and submesothelial matrix marker across both tissue and effluent. |
| **`FN1`** | $\log_2\text{FC} = +1.58, P = 0.012$ | $\log_2\text{FC} = +0.41, P = 0.269$ | **Concordant (UP)** | Preserved upregulation trend in effluent myofibroblasts. |
| **`COL3A1`** | $\log_2\text{FC} = +1.28, P = 0.021$ | $\log_2\text{FC} = +0.19, P = 0.616$ | **Concordant (UP)** | Fibrillar collagen upregulated in late dialysis effluent. |
| **`ISM1`** | $\log_2\text{FC} = +0.89, P = 0.043$ | $\log_2\text{FC} = +0.03, P = 0.904$ | **Concordant (UP)** | Modest concordant trend in microvascular endothelial transcripts. |
| **`LOX`** | $\log_2\text{FC} = +0.84, P = 0.049$ | $\log_2\text{FC} = +0.02, P = 0.976$ | **Flat / Neutral** | Secreted crosslinking enzyme; cellular transcript unchanged in shed effluent cells. |
| **`VCAN`** | $\log_2\text{FC} = +1.74, P = 0.008$ | $\log_2\text{FC} = -0.52, P = 0.024$ | **Discordant (DOWN)** | **Compartmental Duality:** Proteoglycan matrix expands in fixed parietal tissue but is shed/depleted in cellular effluent fraction. |
| **`THBS3`** | $\log_2\text{FC} = +1.02, P = 0.019$ | $\log_2\text{FC} = -0.20, P = 0.418$ | **Discordant (DOWN)** | Glycoprotein anchored to ECM scaffold; reduced transcript in free effluent leukocytes. |

---

### Phase 6: Clinical Diagnostic Nomogram, Calibration & Decision Curve Analysis (DCA)
A 5-gene multivariable diagnostic nomogram was constructed on the directional predictors (`VCAN`, `COL8A1`, `FN1`, `ISM1`, `COL3A1`) in GSE125498.

* **Nomogram Points:** Dynamic ranges mapped to a 0–100 scale ($\text{COL8A1} = 100.0\text{ max pts}, \text{VCAN} = 84.5\text{ pts}, \text{COL3A1} = 58.6\text{ pts}, \text{FN1} = 33.3\text{ pts}, \text{ISM1} = 17.1\text{ pts}$).
* **In-Sample Discrimination:** $\text{C-index} = \mathbf{0.819}$ ($95\%\text{ CI: } [0.623, 0.968]$).
* **Likelihood-Ratio Test vs. `FN1` Alone:** $\chi^2 = 11.309, \text{df} = 4, P = \mathbf{0.0233}$ (significant improvement over single-gene model).
* **Calibration & Decision Curves:** Lowess-smoothed calibration ($B=1000$ bootstrap, Brier Score $= 0.1542$, Hosmer-Lemeshow $P = 0.2437$) and Decision Curve Analysis (DCA) demonstrate net clinical benefit across threshold probabilities $P_t = 0.10 - 0.70$.

> [!NOTE]
> The nomogram, calibration curve, and DCA metrics are fitted and evaluated on the same 33-sample cohort. Cross-validated discrimination is lower ($\text{AUC}_{\text{CV}} = 0.550$), reflecting in-sample optimism.

![Clinical Nomogram & DCA](results/figures/Hub_02_clinical_nomogram_dca_calibration.png)
* **Figure 9: Clinical Diagnostic Nomogram Suite.** (A) Exact points scoring nomogram. (B) Bootstrap calibration curve. (C) Decision Curve Analysis showing net clinical benefit across risk thresholds. (D) Cross-validated multi-model comparison.

---

### Phase 7: Protein-Protein Interaction (STRING v12.0) & GSE62928 Co-Expression
An integrated molecular interaction network was constructed combining live-verified physical/functional interactions from **STRING v12.0** with empirical co-expression in human peritoneal tissue.

* **Live STRING v12.0 PPI Network:** 20 verified physical/functional interaction edges (interaction score $\ge 0.400$; top edges: `COL11A1-COL3A1` [0.965], `COMP-FN1` [0.958], `FN1-LOX` [0.931], `COL3A1-FN1` [0.928]).
* **Co-Expression Edges ($N=8$ Discovery):** 42 / 55 (76.4%) hub gene pairs exhibit $|r| \ge 0.85$ (42 edges at $|r| \ge 0.85, \text{FDR} < 0.01$).
* **Co-Expression Noise Baseline:** Random sampling of 100,000 non-hub gene pairs in GSE62928 establishes that **$1.68\%$ of random gene pairs exceed $|r| \ge 0.85$ purely by chance** at $N=8$. While hub genes show strong enrichment ($76.4\%$ vs $1.68\%$), individual correlation values carry small-sample estimation variance and are presented alongside STRING physical edges.

![PPI Interaction Network](results/figures/Hub_01_ppi_gene_interaction_network.png)
* **Figure 10: Integrated STRING PPI & Co-Expression Network.** High-confidence physical/functional STRING v12.0 interactions (score $\ge 0.400$) and empirical co-expression edges among the 11 consensus hub genes.

---

### Phase 8: Permutation-Controlled Preranked GSEA
Preranked Gene Set Enrichment Analysis was executed using `gseapy` on all 20,940 genes ranked by Limma moderated $t$-statistics across MSigDB Hallmark Gene Sets (v2020).

* **Primary Enriched Hallmark:** **Epithelial-Mesenchymal Transition (EMT)** ($\text{NES} = +\mathbf{3.115}, P_{\text{nom}} < 10^{-4}, \text{FDR} < 10^{-4}$).
* **EMT Permutation Test:** Across all 70 label permutations, EMT $\text{NES} = +3.115$ ranked **#1 out of 70** ($P_{\text{perm}} = 1/70 = \mathbf{0.0143}$).
* **Other Enriched Pathways:** TNF-$\alpha$ Signaling via NF-$\kappa$B ($\text{NES} = +2.082, \text{FDR} = 0.001$), Inflammatory Response ($\text{NES} = +1.587, \text{FDR} = 0.030$), Angiogenesis ($\text{NES} = +1.438, \text{FDR} = 0.086$), Hypoxia ($\text{NES} = +1.349, \text{FDR} = 0.097$).
* **Global Permutation Control:** $9 / 70$ random permutations generate $\ge 35$ nominal Hallmarks ($P_{\text{perm}} = 0.1286$), demonstrating that the focal EMT signal is robust while aggregate pathway counts reflect genome-wide co-expression correlations.

![GSEA Pathway Heatmap](results/figures/Hub_03_gsea_pathway_enrichment_heatmap.png)
* **Figure 11: Hallmark Fibrosis Pathway Enrichment.** Heatmap of circularity-corrected ssGSEA hallmark pathway scores across peritoneal tissue samples.

---

### Phase 9: Microenvironmental Immune Deconvolution & Correlation
Peritoneal cell-type deconvolution was performed across 12 immune and stromal populations using curated published reference signatures (Charoentong et al., Bindea et al., Rossi et al.). All signatures were screened to ensure no hub gene is present within any cell marker list.

* **Focal Cell Type Expansion:** **Peritoneal Myofibroblasts** ($U = 16.0, P_{\text{raw}} = 0.0286$, Permutation $P_{\text{perm}} = 2/70 = \mathbf{0.0286}$).
* **Borderline Cell Types:** M2 Macrophages ($P = 0.0571, \text{FDR} = 0.1714$), Neutrophils ($P = 0.0571, \text{FDR} = 0.1714$), Activated Dendritic Cells ($P = 0.0571, \text{FDR} = 0.1714$).
* **Global Permutation Baseline:** In 70 random label permutations, **$31.4\%$ of random permutations produce $\ge 1$ cell type with nominal $P < 0.05$** (Family-wise $P_{\text{perm}} = 0.3143$). Individual cell-type signals must be interpreted as exploratory.

![Immune Infiltration Deconvolution](results/figures/Hub_04_immune_infiltration_deconvolution.png)
* **Figure 12: Peritoneal Immune Deconvolution.** (A) Boxplots of standardized microenvironmental infiltration scores ($N = 8$). (B) Spearman correlation matrix between the 11 hub genes and infiltrating cell populations.

---

## 🔬 Relationship to Prior Literature

Prior transcriptomic analyses of GSE62928 (e.g., Wang et al., 2024) reported a 4-gene diagnostic panel (`SOCS1`, `PIM2`, `HSH2D`, `MYO3B`).

* **Zero Gene Overlap ($0 / 4$):** None of the Wang et al. genes are Matrisome structural proteins; they represent immune-signaling and intracellular kinase genes.
* **Distinct Biological Focus:** Our pipeline focuses explicitly on the **extracellular matrix architecture and stromal fibrogenesis machinery** (`FN1`, `COL3A1`, `COL8A1`, `COL11A1`, `VCAN`, `COMP`, `THBS3`, `EDIL3`, `LOX`, `INHBA`, `ISM1`).
* **Cross-Validation Rigor:** The Wang et al. panel was reported with uncorrected in-sample AUCs. Our analysis provides exact label permutations, out-of-fold cross-validation, and negative-control baselines.

---

## ⚠️ Comprehensive Study Limitations

1. **Discovery Sample Size ($N = 8$):** Encapsulating Peritoneal Sclerosis (EPS) is a rare clinical entity, and biopsy tissue is ethically and clinically constrained. Although exact combinatorial permutations confirm that our primary findings (`71 ECM-DEGs`, `Salmon Module`, `EMT NES`, `Myofibroblasts`) achieve the minimum possible $P$-values at $N=8$, statistical power is fundamentally bounded by sample size.
2. **Multiple Testing Correction:** In external validation ($N=33$), zero of the seven profiled hub genes survive Bonferroni or FDR correction, and the discovery Salmon module fails family-wise error control ($P_{\text{bonf}} = 0.2198$).
3. **Saturated Model Hazards in Small-$N$ Transcriptomics:** Multi-gene models with $p > n$ (e.g. 11 predictors on 8 samples) achieve mathematical separation ($\text{AUC} = 1.0$) trivially. Reporting in-sample discovery multi-gene AUCs creates a misleading impression of certainty; our study rejects these metrics.
4. **Tissue vs. Effluent Compartment Duality:** Discovery profiling was performed on fixed parietal peritoneal tissue biopsies, whereas external validation used shed peritoneal effluent cells. Extracellular matrix proteoglycans (`VCAN`) and glycoproteins (`THBS3`) exhibit discordant cellular expression between tissue and effluent fractions.
5. **Optimism Penalty in Composite Classifiers:** The composite 7-gene model achieves an in-sample $\text{AUC} = 0.869$, but cross-validation reveals a realistic generalization of **$\text{AUC} \approx 0.68 - 0.70$**.
6. **Absence of Prospective Wet-Lab Validation:** All findings are derived from in silico microarray re-analyses. Prospective clinical biopsy immunohistochemistry, RNAscope, and targeted RT-qPCR in large cohorts ($N \ge 100$) are required before clinical translation.

---

## 🛡️ Data and Code Integrity Disclosure

During a comprehensive internal quality control and reproducibility audit conducted in September 2026, two downstream exploratory scripts (`11_drug_repurposing_dgidb.py` and `12_ihc_protein_validation.py`) were identified as containing AI-generated synthetic identifiers (hallucinated antibody catalog numbers and mismatched drug-gene citation PMIDs) that do not correspond to real external database records.

In accordance with strict scientific integrity standards:
* **Complete Removal & History Scrubbing:** Both scripts and their generated outputs (`results/tables/candidate_drugs_*.csv`, `results/tables/hub_genes_ihc_*.csv`, `results/figures/Hub_05_*`, `results/figures/Hub_06_*`) were completely excised and scrubbed from the repository and Git commit history. The entire active commit history reflects only verified, reproducible transcriptomics, STRING v12.0 PPI networks, and permutation tests.
* **Verification Scope of Surviving Pipeline:**
  - **Re-Executed & Numerically Verified (September 2026):** Scripts `01_load_qc_preprocess.R`, `02_differential_expression.R`, `03_matrisome_filtering.R`, and `08b_roc_analysis.py` were re-executed from raw matrices, and their numerical outputs were verified against reported figures and tables. Script `07_gene_interaction_network.py` was executed and verified live against the STRING v12.0 API (20 physical edges, score $\ge 0.400$).
  - **Independent Re-Computation of External Validation:** Script `06_external_validation_GSE125498.py` outputs (Limma fold-changes, directional discordances, single-gene ROC curves, and the 50-repeat cross-validated composite AUC of $0.696 \pm 0.056$) were independently re-parsed and re-computed directly from `GSE125498_family.soft.gz` via `audit/compute_audit_step2_3.py`.
  - **Permutation Auditing & Sample-Size Boundaries:** The empirical significance of key findings was evaluated across all $\binom{8}{4} = 70$ exact combinatorial label permutations in `audit/permutation_test_ecm.py`, `audit/stage_c_permutations.py`, and `audit/stage_e_ppi_baseline.py`. The primary ECM over-representation ($P_{\text{perm}} = 0.0143$), WGCNA Salmon module correlation ($P_{\text{perm}} = 0.0143$), and GSEA EMT hallmark enrichment ($P_{\text{perm}} = 0.0143$) each achieve rank 1/70, which represents the mathematical minimum achievable $p$-value at $N=8$.
  - **Family-Wise & Negative-Control Controls:** Across 68 non-true random label splits, WGCNA co-expression modeling yielded at least one significant module in **$26.5\%$** of noise trials (family-wise adjusted $P_{\text{bonf}} = \mathbf{0.2198}$ across 14 modules), and immune deconvolution yielded at least one nominally significant cell type in **$31.4\%$** of random trials (family-wise $P_{\text{perm}} = \mathbf{0.3143}$). GSEA global pathway counts ($\ge 35$ nominal Hallmarks) yielded $P_{\text{perm}} = \mathbf{0.1286}$. These baseline rates establish that while focal signals (`ECM-DEGs`, `Salmon Module`, `EMT Hallmark`, `Peritoneal Myofibroblasts`) are top-ranked, secondary discovery features reflect exploratory candidates requiring prospective validation.
  - **Static & Table Audited:** Scripts `00_fetch_gse62928_matrix.R`, `02b_wgcna_analysis.R`, `04_functional_enrichment.R`, `05b_ml_hub_gene_identification_wgcna.py`, and `08_nomogram_roc_analysis.py` were code-reviewed and cross-checked gene-for-gene against their verified CSV outputs (`convergent_71_ECM_DEGs.csv`, `ML_hub_genes_from_WGCNA_ECM.csv`, `ML_hub_genes_from_WGCNA_ECM_all_results.csv`).
* **Specific Automated Audit Coverage:**
  - `audit/repo_audit.py`: Validates script existence, absence of absolute filepaths, random seed reproducibility, and table dimensions across all active pipeline components.
  - `audit/compute_audit_step2_3.py`: Validates GSE125498 expression metrics, directional concordance, and WGCNA module label-permutation tests.
  - `audit/statistical_rigor_audit.py`: Validates the integrity of the 71-gene ECM filter, 40-candidate convergence, ML vote matrix, and ensures total absence of deprecated Mendelian Randomization artifacts.
  - `audit/permutation_test_ecm.py` & `audit/stage_c_permutations.py`: Validates exact combinatorial permutation distributions for ECM-DEG over-representation, GSEA hallmark scores, and immune deconvolution.
  - `audit/stage_e_ppi_baseline.py`: Validates the empirical 1.68% chance co-expression noise baseline across 100,000 non-hub gene pairs at $N=8$.

---

## 🛠️ Pipeline Execution Guide

The entire analytical pipeline can be executed deterministically using the following sequence:

### Step 1: Normalized Expression Matrix & Quality Control (R)
```bash
Rscript 00_fetch_gse62928_matrix.R
Rscript 01_load_qc_preprocess.R
```

### Step 2: Limma Differential Expression & Matrisome Filtering (R)
```bash
Rscript 02_differential_expression.R
Rscript 03_matrisome_filtering.R
Rscript 04_functional_enrichment.R
```

### Step 3: WGCNA Co-Expression Network Analysis (R)
```bash
Rscript 02b_wgcna_analysis.R
```

### Step 4: Convergence Intersection & Consensus Machine Learning (Python)
```bash
python analyze_convergence_wgcna.py
python 05b_ml_hub_gene_identification_wgcna.py
```

### Step 5: External Cohort Validation in GSE125498 (Python)
```bash
python 06_external_validation_GSE125498.py
```

### Step 6: Downstream Systems Biology Suite (Python)
```bash
python 07_gene_interaction_network.py        # STRING v12.0 PPI & GSE62928 co-expression
python 08b_roc_analysis.py                  # Standardized single & multi-gene ROC curves
python 08_nomogram_roc_analysis.py           # Clinical diagnostic nomogram & DCA
python 09_gsea_pathway_enrichment.py         # True preranked GSEA on Hallmark gene sets
python 10_immune_infiltration_analysis.py    # Microenvironmental immune deconvolution
```

### Step 7: Automated Quality Assurance & Statistical Rigor Audits (Python)
```bash
python audit/statistical_rigor_audit.py
python audit/audit_pipeline_errors.py
python audit/compute_audit_step2_3.py
python audit/repo_audit.py
```

---

## 📁 Repository Structure

```
.
├── 00_fetch_gse62928_matrix.R                     # TASK 0: Full matrix builder & probe mapper
├── 01_load_qc_preprocess.R                       # Primary DEG preprocessing pipeline
├── 02_differential_expression.R                  # Limma differential expression
├── 02b_wgcna_analysis.R                          # TASK 1: WGCNA network & trait correlation
├── 03_matrisome_filtering.R                      # Matrisome masterlist intersection
├── 04_functional_enrichment.R                    # GO/KEGG functional enrichment
├── 05b_ml_hub_gene_identification_wgcna.py       # TASK 3: ML consensus on 40 WGCNA-ECM candidates
├── 06_external_validation_GSE125498.py           # TASK 4: External expression check on GSE125498
├── 07_gene_interaction_network.py                # PPI (STRING v12.0 + co-expression, |r|>=0.85, FDR<0.01)
├── 08b_roc_analysis.py                          # Task 1: ROC + bootstrap 95% CI for all 11 genes
├── 08_nomogram_roc_analysis.py                   # Task 2: Exact nomogram, DCA (external), calibration B=1000
├── 09_gsea_pathway_enrichment.py                 # Task 3: True preranked GSEA + circularity-free ssGSEA
├── 10_immune_infiltration_analysis.py            # Task 4: Wilcoxon FDR + ordered hub-immune correlation
├── analyze_convergence_wgcna.py                  # TASK 2: WGCNA ∩ ECM-DEG convergence & 3-way Venn
├── generate_gse62928_ecm_venn.py                 # Matrisome Venn diagram generation
├── run_pipeline.R                                # Pipeline master execution wrapper
├── ECM genes all.xlsx                            # Curated 1,027 Human Matrisome reference
├── GSE62928.top.table.tsv                        # GSE62928 complete Limma top table
├── convergent_71_ECM_DEGs.csv                    # 71 pro-fibrotic ECM-DEGs
├── requirements.txt                              # Python environment dependency requirements
├── session_info.txt                              # Complete R sessionInfo() & Python package versions
├── audit/
│   ├── repo_audit.py                             # Automated repository consistency audit
│   ├── statistical_rigor_audit.py                # Statistical controls & permutation audit
│   ├── audit_pipeline_errors.py                  # Data integrity & sample checksum audit
│   └── compute_audit_step2_3.py                  # Step 2/3 re-computation & validation scripts
├── Validation/
│   └── GSE125498.top.table.tsv                   # GSE125498 Limma top table
├── results/
│   ├── AUDIT_SUMMARY.md                          # Comprehensive statistical audit report
│   ├── figures/
│   │   ├── graphical_abstract.jpg                # Figure 0: Study Graphical Abstract
│   │   ├── DEG_01_volcano_plot.png               # Figure 1: Discovery DEG Volcano Plot
│   │   ├── venn_gse62928_pro_fibrotic_ecm_71.png # Figure 2: Matrisome Overlap Venn Diagram
│   │   ├── WGCNA_03_module_trait_heatmap.png     # Figure 3: Module-Trait Correlation Heatmap
│   │   ├── venn_wgcna_convergence.png            # Figure 4: 3-Way Convergence Venn Diagram
│   │   ├── WGCNA_ML_01_consensus_votes_barchart.png # Figure 5: Consensus Votes Barchart
│   │   ├── WGCNA_ML_04_hub_genes_expression_heatmap.png # Figure 6: Hub Expression Heatmap
│   │   ├── Validation_01_hub_genes_mann_whitney_boxplots.png # Figure 7: Validation Boxplots
│   │   ├── Hub_02b_roc_analysis.png              # Figure 8: Standardized ROC Curves
│   │   ├── Hub_02_clinical_nomogram_dca_calibration.png # Figure 9: Nomogram, Calibration & DCA
│   │   ├── Hub_01_ppi_gene_interaction_network.png # Figure 10: STRING PPI & Co-Expression
│   │   ├── Hub_03_gsea_pathway_enrichment_heatmap.png # Figure 11: GSEA Hallmark Heatmap
│   │   └── Hub_04_immune_infiltration_deconvolution.png # Figure 12: Immune Infiltration Boxplots
│   └── tables/
│       ├── GSE62928_full_expression_matrix.csv
│       ├── GSE62928_sample_metadata.csv
│       ├── wgcna_module_trait_correlation.csv
│       ├── wgcna_trait_significant_module_genes.csv
│       ├── convergent_WGCNA_ECM_genes.csv
│       ├── ML_hub_genes_from_WGCNA_ECM.csv
│       ├── GSE125498_wgcna_hub_validation_metrics.csv
│       ├── GSE125498_expression_matrix_hubs.csv
│       ├── GSE125498_sample_metadata.csv
│       ├── statistical_rigor_audit.csv
│       ├── roc_auc_detailed_metrics.csv
│       ├── hub_genes_ppi_centrality_metrics.csv
│       ├── nomogram_point_scoring_table.csv
│       ├── nomogram_calibration_curve_metrics.csv
│       ├── decision_curve_analysis_metrics.csv
│       ├── gsea_preranked_hallmark_results.csv
│       ├── hub_genes_gsea_pathway_correlations.csv
│       ├── immune_infiltration_scores.csv
│       ├── hub_genes_immune_correlations.csv
│       └── immune_infiltration_wilcoxon_tests.csv
└── README.md                                     # Master documentation
```

---

## 📜 Citation & License
This project is released under the **MIT License**. For academic usage or reproduction, please cite the underlying data repositories:
* **GSE62928:** Human peritoneal membrane transcriptomics in Encapsulating Peritoneal Sclerosis (Affymetrix HG-U133_Plus_2).
* **GSE125498:** Longitudinal human peritoneal dialysis effluent cell transcriptomics (Illumina HumanHT-12 v4.0).
* **Matrisome Project:** Naba et al., *Matrix Biology*, 2012 / 2016.
