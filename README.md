# Peritoneal Dialysis-Associated Peritoneal Fibrosis: Multi-Omics, Causal Discovery & Machine Learning Pipeline

A comprehensive computational biology and genetic epidemiology framework integrating bulk transcriptomics, extracellular matrix (ECM) matrisome biology, consensus machine learning (LASSO, SVM-RFE, Random Forest, XGBoost), transcriptome-wide two-sample Mendelian Randomization (TWMR), and external cohort validation (GSE125498) to identify diagnostic and therapeutic biomarkers for peritoneal dialysis (PD)-induced peritoneal membrane injury and fibrosis.

---

## 📌 Project Overview & Roadmap

Long-term peritoneal dialysis (LPD) leads to chronic mesothelial injury, epithelial-to-mesenchymal transition (EMT), extracellular matrix accumulation, and ultrafiltration failure. This pipeline identifies early-stage detection biomarkers and late-stage progression drivers across multi-omics layers:

```
  ┌────────────────────────────────────────────────────────────────────────────────────────────────┐
  │                                    STUDY ROADMAP & PROGRESS                                    │
  └────────────────────────────────────────────────────────────────────────────────────────────────┘
                                                  │
             ┌────────────────────────────────────┴────────────────────────────────────┐
             ▼                                                                         ▼
   [PHASE 1: GSE62928 & MATRISOME]                                           [PHASE 2: TWMR CAUSAL SCREEN]
   - GSE62928 DEGs: 1,526 genes                                              - eQTLGen whole-blood cis-eQTLs (n=31,684)
   - Curated Matrisome: 1,027 genes                                          - CKDGen eGFR GWAS (n=567,460)
   - Convergence: 71 ECM-DEGs (956 exclusive ECM)                            - Transcriptome-wide causal screen (FDR < 0.05)
             └────────────────────────────────────┬────────────────────────────────────┘
                                                  ▼
                         ┌─────────────────────────────────────────────────┐
                         │   CONVERGENCE: 15 CAUSAL MATRISOME GENES        │
                         │   6 Upregulated: P4HA2, ADAMTS1, TNC, etc.      │
                         │   9 Downregulated: BMP6, IGFBP3, WNT11, etc.    │
                         └────────────────────────┬────────────────────────┘
                                                  │
             ┌────────────────────────────────────┴────────────────────────────────────┐
             ▼                                                                         ▼
   [PANEL A: 4 CAUSAL HUB GENES]                                             [PANEL B: 12 MATRISOME HUB GENES]
   - Source: 15 Causal ECM Genes                                             - Source: 71 Convergent ECM-DEGs
   - 4 ML Models: LASSO, SVM-RFE, RF, XGBoost                                - 4 ML Models: LASSO, SVM-RFE, RF, XGBoost
   - Members (>= 2 votes):                                                   - Members (>= 2 votes):
     P4HA2 (3), ADAMTS1 (3), WNT11 (2), TNC (2)                                ISM1 (3), TGM2 (3), MXRA5 (2), COL3A1 (2),
                                                                               COL5A2 (2), COL11A1 (2), EDIL3 (2), POSTN (2),
                                                                               LOX (2), INHBA (2), COMP (2), THBS3 (2)
             └────────────────────────────────────┬────────────────────────────────────┘
                                                  ▼
                         ┌─────────────────────────────────────────────────┐
                         │    PHASE 4: EXTERNAL VALIDATION (GSE125498)     │
                         │    Early-Stage (SPD, n=20) vs Late (LPD, n=13)  │
                         │    - Mann-Whitney U Tests                       │
                         │    - ROC Discrimination (Composite AUC = 0.731) │
                         │    - Longitudinal Stage Trajectories            │
                         │    - Clustered Patient Heatmaps                 │
                         └─────────────────────────────────────────────────┘
```

---

## 🔬 Gene Panels Identified by Consensus Machine Learning

Two distinct biomarker panels were generated using a rigorous **4-algorithm consensus machine learning ensemble** (LASSO L1-regularization, Support Vector Machine Recursive Feature Elimination [SVM-RFE], Random Forest Gini Importance, and XGBoost Gain Importance):

### 1. Panel A: The 4-Hub Causal Biomarker Panel
* **Origin**: Screened from the **15 Causal ECM Genes** identified by converging GSE62928 DEGs with transcriptome-wide Mendelian Randomization (eQTLGen $\times$ CKDGen eGFR).
* **Consensus Selection Threshold**: Selected by $\ge 2$ out of 4 independent ML models.
* **Panel Members**:
  1. **`P4HA2`** (3 votes: LASSO, SVM-RFE, Random Forest) — Key prolyl 4-hydroxylase essential for collagen triple-helix stabilization; top genetically causal driver ($\text{FDR} = 3.61 \times 10^{-7}$).
  2. **`ADAMTS1`** (3 votes: LASSO, SVM-RFE, Random Forest) — Matrix metalloproteinase regulator of collagen assembly and cell-matrix interactions.
  3. **`WNT11`** (2 votes: SVM-RFE, XGBoost) — Non-canonical Wnt morphogen regulating mesothelial polarity and EMT.
  4. **`TNC`** (2 votes: SVM-RFE, Random Forest) — Tenascin-C, mechanosensitive hexameric matricellular glycoprotein induced during fibrotic stress.
* **Associated Data Table**: [results/tables/ML_hub_genes_final_list.csv](file:///d:/Peritoneal%20Project/results/tables/ML_hub_genes_final_list.csv)
* **Associated Visualizations & Graphs**:
  * **Consensus Votes Bar Chart**: [results/figures/ML_01_consensus_votes_barchart.png](file:///d:/Peritoneal%20Project/results/figures/ML_01_consensus_votes_barchart.png) — Displays vote counts across the 4 ML models for the 15 causal candidates.
  * **Model Agreement Heatmap**: [results/figures/ML_02_model_selection_heatmap.png](file:///d:/Peritoneal%20Project/results/figures/ML_02_model_selection_heatmap.png) — Binary presence/absence matrix showing exact model agreement.
  * **Per-Model Feature Importance (2x2 Panel)**: [results/figures/ML_03_per_model_importance_2x2.png](file:///d:/Peritoneal%20Project/results/figures/ML_03_per_model_importance_2x2.png) — Ranked feature weights for LASSO coefficients, SVM-RFE ranking, RF mean decrease in impurity, and XGBoost gain.
  * **Hub Gene Expression Heatmap**: [results/figures/ML_04_hub_genes_expression_heatmap.png](file:///d:/Peritoneal%20Project/results/figures/ML_04_hub_genes_expression_heatmap.png) — Clustered expression levels across peritoneal samples in GSE62928.
  * **15-Gene Causal Convergence Venn Diagram**: [venn_15_causal_convergence.png](file:///d:/Peritoneal%20Project/venn_15_causal_convergence.png) — Triple Venn diagram showing the convergence between GSE62928 DEGs (1,526), Curated Matrisome (1,027), and TWMR Causal Genes (15).

---

### 2. Panel B: The 12-Hub Matrisome Structural Panel
* **Origin**: Screened across **all 71 Convergent ECM-DEGs** overlapping between GSE62928 ($P < 0.05$) and the complete Human Matrisome reference (`ECM genes all.xlsx`).
* **Consensus Selection Threshold**: Selected by $\ge 2$ out of 4 independent ML models.
* **Panel Members**:
  1. **`ISM1`** (3 votes: LASSO, SVM-RFE, RF) — Isthmin 1, angiogenesis and mesothelial survival factor.
  2. **`TGM2`** (3 votes: LASSO, SVM-RFE, RF) — Transglutaminase 2, primary collagen/fibronectin cross-linking enzyme driving peritoneal stiffening.
  3. **`MXRA5`** (2 votes: SVM-RFE, RF) — Matrix-remodelling associated 5, anti-inflammatory and fibrotic matrix protector.
  4. **`COL3A1`** (2 votes: SVM-RFE, RF) — Collagen type III alpha 1 chain, dominant fibrillar collagen in peritoneal membrane expansion.
  5. **`COL5A2`** (2 votes: SVM-RFE, RF) — Collagen type V alpha 2 chain, regulator of collagen fibrillogenesis.
  6. **`COL11A1`** (2 votes: SVM-RFE, RF) — Collagen type XI alpha 1 chain.
  7. **`EDIL3`** (2 votes: SVM-RFE, RF) — EGF-like repeats and discoidin I-like domains 3, endothelial/mesothelial adhesion.
  8. **`POSTN`** (2 votes: SVM-RFE, RF) — Periostin, mechanosensitive matricellular ligand promoting cell motility and EMT.
  9. **`LOX`** (2 votes: SVM-RFE, RF) — Lysyl oxidase, initiator of covalent collagen cross-linking.
  10. **`INHBA`** (2 votes: SVM-RFE, RF) — Inhibin beta A (Activin A), major upstream driver of peritoneal fibroblast activation.
  11. **`COMP`** (2 votes: SVM-RFE, RF) — Cartilage oligomeric matrix protein.
  12. **`THBS3`** (2 votes: SVM-RFE, XGBoost) — Thrombospondin 3, adhesive glycoprotein regulating cell-matrix interactions.
* **Associated Data Table**: [results/tables/ML_12_hub_genes_from_71_ECM_DEGs.csv](file:///d:/Peritoneal%20Project/results/tables/ML_12_hub_genes_from_71_ECM_DEGs.csv)
* **Associated Visualizations & Graphs**:
  * **71-Gene Consensus Vote Distribution**: [results/figures/ML_01_consensus_votes_71_genes.png](file:///d:/Peritoneal%20Project/results/figures/ML_01_consensus_votes_71_genes.png) — Barchart of all 71 ECM-DEGs ranked by model votes, isolating the top 12 consensus hub genes.
  * **12-Hub Gene Expression Heatmap**: [results/figures/ML_04_hub_genes_from_71_heatmap.png](file:///d:/Peritoneal%20Project/results/figures/ML_04_hub_genes_from_71_heatmap.png) — Clustered expression heatmap of the 12 Matrisome Hub Genes in GSE62928.

---

## 📊 External Cohort Validation: Early- vs Late-Stage Detection (GSE125498)

To evaluate exploratory diagnostic utility for discriminating early peritoneal dialysis exposure from progressive membrane fibrosis, prioritized biomarkers were tested in the independent human effluent cohort **GSE125498**:
* **Platform**: Illumina HumanHT-12 V4.0 (`GPL10558`, 19,164 probes present in dataset).
* **Sample Stratification**:
  * **Early Stage (Short-Term PD, SPD: 0–24 Months)**: $n = 20$ patients (preserved membrane transport).
  * **Late Stage (Long-Term PD, LPD: $\ge 25$ Months)**: $n = 13$ patients (progressive membrane injury, high solute transport, established fibrotic remodeling).

### Complete Validation Metrics (All 19 Evaluated Candidates)

Every hub and causal candidate profiled in GSE125498 is reported below, sorted by discrimination AUC-ROC. No candidates with non-significant or null results are omitted:

| Biomarker | Panel Classification | $\log_2\text{FC}$ (Late vs Early) | Expression Trend | Mann-Whitney $U$ | Mann-Whitney $P$-value | Individual AUC-ROC | Validation Outcome |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **ADAM28** | Causal ECM | +0.55 | Upregulated in Late | 180.0 | 0.068 | 0.692 | Trend toward significance |
| **LTBP4** | Causal ECM | +0.65 | Upregulated in Late | 174.5 | 0.105 | 0.671 | Trend toward significance |
| **IGFBP3** | Causal ECM | +0.51 | Upregulated in Late | 169.0 | 0.156 | 0.650 | Non-significant |
| **TGM2** | 12-Hub Matrisome | +0.40 | Upregulated in Late | 168.0 | 0.167 | 0.646 | Non-significant (Top 12-Hub) |
| **FGL2** | Causal ECM | -0.18 | Downregulated in Late | 99.0 | 0.261 | 0.619 | Non-significant |
| **MXRA5** | 12-Hub Matrisome | +0.35 | Upregulated in Late | 157.0 | 0.329 | 0.604 | Non-significant |
| **COL4A2** | Causal ECM | +0.22 | Upregulated in Late | 156.5 | 0.338 | 0.602 | Non-significant |
| **COL5A2** | 12-Hub Matrisome | +0.40 | Upregulated in Late | 156.0 | 0.347 | 0.600 | Non-significant |
| **BMP6** | Causal ECM | +0.34 | Upregulated in Late | 156.0 | 0.347 | 0.600 | Non-significant |
| **THBS3** | 12-Hub Matrisome | -0.20 | Downregulated in Late | 105.0 | 0.367 | 0.596 | Non-significant |
| **ADAMTS1** | 4-Hub Causal | -0.16 | Downregulated in Late | 112.0 | 0.519 | 0.569 | Non-significant |
| **P4HA2** | 4-Hub Causal | +0.17 | Upregulated in Late | 145.0 | 0.593 | 0.558 | Non-significant |
| **FBLN5** | Causal ECM | +0.02 | Stable | 144.0 | 0.619 | 0.554 | Non-significant |
| **POSTN** | 12-Hub Matrisome | +0.12 | Upregulated in Late | 141.0 | 0.699 | 0.542 | Non-significant |
| **COL3A1** | 12-Hub Matrisome | +0.19 | Upregulated in Late | 139.0 | 0.754 | 0.535 | Non-significant |
| **ISM1** | 12-Hub Matrisome | +0.03 | Stable | 137.0 | 0.811 | 0.527 | Non-significant |
| **LOX** | 12-Hub Matrisome | +0.02 | Stable | 129.0 | 0.985 | 0.496 | Non-significant |
| **CRISPLD2** | Causal ECM | -0.07 | Stable | 130.0 | 1.000 | 0.500 | Null discrimination |
| **TNC** | 4-Hub Causal | +0.04 | Stable | 130.0 | 1.000 | 0.500 | **Validation Failure** |

### Individual vs Multi-Gene Composite Performance

* **Individual-Gene Analysis (Lack of Stand-Alone Significance)**:  
  **No single gene achieved nominal statistical significance ($p < 0.05$)** in the two-sided Mann-Whitney U tests. The strongest individual discriminator in the 12-Hub panel was `TGM2` ($p = 0.167$, $\text{AUC} = 0.646$), while `ADAM28` ($p = 0.068$, $\text{AUC} = 0.692$) and `LTBP4` ($p = 0.105$, $\text{AUC} = 0.671$) in the causal ECM candidate set showed marginal non-significant trends. Stand-alone individual genes therefore lack sufficient diagnostic power in this effluent cohort.
* **TNC Validation Failure**:  
  Tenascin-C (`TNC`), a core member of the 4-Hub Causal Panel and canonical fibrotic matricellular protein in discovery tissue biopsies, showed virtually no differential signal in effluent cells ($\log_2\text{FC} = +0.04$, Mann-Whitney $U = 130.0$, $p = 1.000$, $\text{AUC} = 0.500$). This indicates an explicit **validation failure for TNC** in this cellular effluent dataset, possibly due to detachment-induced dilution or non-shedding of matrix-anchored tenascin into peritoneal effluent.
* **Multi-Gene Composite Signatures (Exploratory Discrimination)**:  
  While individual genes were non-significant, combining them into multi-gene classifiers provided exploratory separation between early and late stages:
  * **12-Hub Matrisome Composite Signature (8 profiled genes)**:  
    $$\mathbf{AUC = 0.731} \quad (95\%\ \text{CI: } 0.54 - 0.89)$$  
    Captures combined matrix remodeling across collagen expansion (`COL3A1`, `COL5A2`) and cross-linking (`TGM2`, `MXRA5`).
  * **4-Hub Causal Composite Signature (3 profiled genes)**:  
    $$\mathbf{AUC = 0.642} \quad (95\%\ \text{CI: } 0.45 - 0.82)$$  
    Anchored by genetically prioritized causal drivers (`P4HA2`, `ADAMTS1`, `TNC`).

### ⚠️ Validation Limitations

1. **Small Sample Size & Limited Power**: The GSE125498 cohort comprises only 33 total patients ($n = 20$ SPD vs $n = 13$ LPD). This small sample size severely restricts statistical power, making these validation results exploratory.
2. **Lack of Individual-Gene Significance**: None of the 19 candidates tested reached $p < 0.05$ individually; diagnostic separation was only detectable when evaluating multi-gene composite models.
3. **Platform & Dataset Coverage Gaps**:
   * **4-Hub Causal Panel**: Only 3 of the 4 genes (`P4HA2`, `ADAMTS1`, `TNC`) could be evaluated. `WNT11` was excluded from validation because no probe exists in the uploaded GSE125498 expression matrix (19,164 probes).
   * **12-Hub Matrisome Panel**: Only 8 of 12 genes (`ISM1`, `TGM2`, `MXRA5`, `COL3A1`, `COL5A2`, `POSTN`, `LOX`, `THBS3`) were profiled; `COL11A1`, `EDIL3`, `INHBA`, and `COMP` lacked probes in the uploaded matrix.
4. **Discordance for Specific Hub Drivers**: As noted above, `TNC` showed complete lack of discrimination ($p = 1.000, \text{AUC} = 0.500$), demonstrating that not all tissue-derived causal targets translate directly into cellular effluent biomarkers.
5. **Requirement for External Confirmation**: These findings are preliminary and hypothesis-generating. They require rigorous prospective confirmation in larger, multi-center cohorts with matched peritoneal biopsy histology and longitudinal ultrafiltration tracking before considering any clinical or diagnostic application.

### Associated Validation Visualizations & Graphs:
* **Mann-Whitney U Boxplots**: [results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png](file:///d:/Peritoneal%20Project/results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png) — Distribution boxplots comparing Early-Stage (SPD) vs Late-Stage (LPD) expression levels with exact non-parametric $p$-values.
* **ROC Curves (Early vs Late Detection)**: [results/figures/Validation_02_roc_curves_early_vs_late.png](file:///d:/Peritoneal%20Project/results/figures/Validation_02_roc_curves_early_vs_late.png) — Sensitivity vs 1-Specificity curves for individual genes and multi-gene composite models ($\text{AUC} = 0.731$ for 12-Hub, $\text{AUC} = 0.642$ for 4-Hub).
* **Stage Progression Trajectories**: [results/figures/Validation_03_stage_progression_trajectories.png](file:///d:/Peritoneal%20Project/results/figures/Validation_03_stage_progression_trajectories.png) — Line charts illustrating continuous biomarker shifts from short-term to long-term dialysis.
* **Clustered Patient Cohort Heatmap**: [results/figures/Validation_04_patient_cohort_heatmap.png](file:///d:/Peritoneal%20Project/results/figures/Validation_04_patient_cohort_heatmap.png) — Unsupervised hierarchical clustering of all 33 patients annotated by clinical stage (top spacing adjusted to prevent title/annotation overlap).
* **Validation Metrics Table**: [results/tables/GSE125498_hub_genes_validation_metrics.csv](file:///d:/Peritoneal%20Project/results/tables/GSE125498_hub_genes_validation_metrics.csv)

---

## 📁 Repository Structure & Artifact Guide

```
d:/Peritoneal Project/
├── data/
│   ├── GSE62928_family.soft.gz                 # GSE62928 Discovery/Validation cohort SOFT
│   ├── GSE125498_family.soft.gz                # GSE125498 Early/Late validation cohort SOFT
│   └── ECM genes all.xlsx                      # Curated 1,027 Human Matrisome master database
│
├── results/
│   ├── figures/
│   │   ├── ML_01_consensus_votes_barchart.png      # Panel A: 4-Hub consensus votes
│   │   ├── ML_02_model_selection_heatmap.png       # Panel A: 4 ML model agreement heatmap
│   │   ├── ML_03_per_model_importance_2x2.png      # Panel A: 2x2 Feature importance
│   │   ├── ML_04_hub_genes_expression_heatmap.png  # Panel A: 4 Hub genes clustered heatmap
│   │   ├── ML_01_consensus_votes_71_genes.png      # Panel B: 12-Hub consensus votes from 71 genes
│   │   ├── ML_04_hub_genes_from_71_heatmap.png     # Panel B: 12 Hub genes clustered heatmap
│   │   ├── Validation_01_hub_genes_mann_whitney_boxplots.png # GSE125498 Mann-Whitney U boxplots
│   │   ├── Validation_02_roc_curves_early_vs_late.png        # GSE125498 ROC discrimination curves
│   │   ├── Validation_03_stage_progression_trajectories.png  # GSE125498 Stage progression trajectories
│   │   ├── Validation_04_patient_cohort_heatmap.png          # GSE125498 Clustered patient heatmap
│   │   ├── venn_15_causal_convergence.png                    # 15 Causal ECM convergence Venn
│   │   └── venn_mr_deg_convergence.png                       # TWMR vs Peritoneal DEG Venn
│   │
│   └── tables/
│       ├── ML_hub_genes_final_list.csv             # Panel A: 4 Causal Hub Genes table
│       ├── ML_12_hub_genes_from_71_ECM_DEGs.csv    # Panel B: 12 Matrisome Hub Genes table
│       ├── convergent_15_causal_genes.csv          # 15 Causal ECM genes complete statistics
│       ├── GSE125498_hub_genes_validation_metrics.csv # GSE125498 Early vs Late validation metrics
│       └── convergent_ECM_DEGs_nominal.csv         # 71 Convergent ECM-DEGs from GSE62928
│
├── 05_ml_hub_gene_identification.py            # 4-algorithm ML consensus hub selection script
├── 06_external_validation_GSE125498.py         # GSE125498 Mann-Whitney U, ROC, and trajectory validation
├── analyze_convergence.py                      # 71-gene and 15-gene mathematical convergence script
├── generate_12_gene_venn.py                    # Venn diagram generation suite
├── transcriptome_wide_mr.py                     # Transcriptome-wide Two-Sample MR screen
└── README.md                                   # Master project documentation
```

---

## 🚀 Execution Instructions

### 1. Extract Machine Learning Hub Genes
```bash
python 05_ml_hub_gene_identification.py
```
*Generates both the 4 Causal Hub Genes and the 12 Matrisome Hub Genes panels with complete importance plots and consensus tables.*

### 2. Run GSE125498 External Validation
```bash
python 06_external_validation_GSE125498.py
```
*Calculates Mann-Whitney U test statistics, ROC-AUC metrics, stage trajectories, and clinical heatmaps.*

### 3. Generate Convergence Venn Diagrams
```bash
python generate_12_gene_venn.py
```
