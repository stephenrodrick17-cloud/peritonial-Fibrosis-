# Peritoneal Dialysis-Associated Peritoneal Fibrosis: Multi-Omics, WGCNA Network Biology & Consensus Machine Learning Pipeline

A comprehensive computational biology and systems medicine framework integrating bulk transcriptomics (GSE62928), extracellular matrix (ECM) matrisome biology, **Weighted Gene Co-expression Network Analysis (WGCNA)**, consensus machine learning (LASSO, SVM-RFE, Random Forest, XGBoost), and external clinical validation (GSE125498) to identify diagnostic and therapeutic biomarkers for peritoneal dialysis (PD)-induced peritoneal membrane injury and fibrosis.

---

## 📌 Executive Summary & Methodological Evolution

### Methodological Shift: From Mendelian Randomization to WGCNA Co-expression Networks
Originally, this pipeline employed transcriptome-wide Two-Sample Mendelian Randomization (TWMR) using systemic whole-blood cis-eQTLs (eQTLGen) and renal function GWAS (CKDGen eGFR) as a genetic causal filtering step. While genetically informative, systemic blood eQTLs often fail to capture localized peritoneal tissue-specific regulatory architectures and coordinated extracellular matrix co-expression modules.

**To overcome this limitation, the MR causal filtering step has been replaced by Weighted Gene Co-expression Network Analysis (WGCNA).**  
This represents a deliberate methodological shift from germline causal inference to **localized biological network co-expression and clinical trait correlation**. WGCNA directly identifies clusters (modules) of highly co-regulated genes within peritoneal tissue that correlate with the clinical disease phenotype (Peritoneal Fibrosis / Encapsulating Peritoneal Sclerosis vs Controls).

```
  ┌────────────────────────────────────────────────────────────────────────────────────────────────┐
  │                               STUDY PIPELINE ARCHITECTURE (WGCNA)                              │
  └────────────────────────────────────────────────────────────────────────────────────────────────┘
                                                  │
             ┌────────────────────────────────────┴────────────────────────────────────┐
             ▼                                                                         ▼
   [DISCOVERY COHORT (GSE62928)]                                             [CURATED HUMAN MATRISOME]
   - Affymetrix Human Gene 1.0 ST (GPL13158)                                 - Naba et al. Extracellular Matrix Database
   - 20,940 unique genes collapsed by MaxMean                                - 1,027 Curated ECM Glycoproteins, Collagens,
   - 1,526 DEGs (Nominal P < 0.05)                                             Proteoglycans & Regulators
             └────────────────────────────────────┬────────────────────────────────────┘
                                                  ▼
                         ┌─────────────────────────────────────────────────┐
                         │      71 CONVERGENT ECM-DEGs (Nominal P < 0.05)   │
                         └────────────────────────┬────────────────────────┘
                                                  │
             ┌────────────────────────────────────┴────────────────────────────────────┐
             ▼                                                                         ▼
   [TASK 1: WGCNA NETWORK BIOLOGY]                                           [TASK 2: CONVERGENCE SCREEN]
   - Signed co-expression network                                            - Trait-Significant Module (Salmon, n=604)
   - Soft-threshold power β = 12 (R² = 0.809)                                  ∩ 71 Convergent ECM-DEGs
   - 14 biologically discrete merged modules                                 - Result: 40 Convergent WGCNA-ECM Candidates
   - Salmon Module: r = +0.806, P = 0.0157 (n = 604)                         - 3-Way Venn: GSE62928 ∩ Matrisome ∩ WGCNA
             └────────────────────────────────────┬────────────────────────────────────┘
                                                  ▼
                         ┌─────────────────────────────────────────────────┐
                         │   TASK 3: CONSENSUS MACHINE LEARNING (4 MODELS) │
                         │   - LASSO (L1 Regularization)                   │
                         │   - SVM-RFE (Recursive Feature Elimination)     │
                         │   - Random Forest (Gini Importance)             │
                         │   - XGBoost (Gain Importance)                   │
                         │   - Voting Threshold: ≥ 2 / 4 Models            │
                         └────────────────────────┬────────────────────────┘
                                                  ▼
                         ┌─────────────────────────────────────────────────┐
                         │    11 WGCNA-ECM CONSENSUS HUB BIOMARKERS        │
                         │    ISM1, FN1, EDIL3, VCAN, COL3A1, COMP,        │
                         │    COL8A1, THBS3, COL11A1, INHBA, LOX           │
                         └────────────────────────┬────────────────────────┘
                                                  ▼
                         ┌─────────────────────────────────────────────────┐
                         │   TASK 4: EXTERNAL CLINICAL VALIDATION          │
                         │   Independent Cohort GSE125498 (N = 33 Patients)│
                         │   - Early Stage (SPD, n=20) vs Late (LPD, n=13) │
                         │   - Mann-Whitney U Tests & ROC Discrimination   │
                         │   ★ Composite Signature AUC = 0.873             │
                         └─────────────────────────────────────────────────┘
```

---

## 🔬 Core Discoveries & Biomarker Panels

### 1. The WGCNA-ECM Consensus Hub Panel (Primary Recommended Panel)
* **Screening Origin**: Selected from the **40 Convergent WGCNA $\cap$ ECM-DEGs** identified by intersecting the trait-significant **Salmon module** ($r = +0.806, p = 0.0157$) with the 71 verified ECM-DEGs.
* **Selection Criterion**: Selected by $\ge 2$ out of 4 independent machine learning algorithms.
* **The 11 Consensus Hub Genes**:
  1. **`ISM1`** (3 votes: LASSO, SVM-RFE, RF) — Isthmin 1; Secreted factor; WGCNA $MM = 0.704, GS = 0.903, \log_2\text{FC} = +1.90$.
  2. **`FN1`** (3 votes: LASSO, SVM-RFE, RF) — Fibronectin 1; Core ECM glycoprotein; WGCNA $MM = 0.920, GS = 0.952, \log_2\text{FC} = +1.93$.
  3. **`EDIL3`** (3 votes: SVM-RFE, RF, XGBoost) — EGF-like repeats and discoidin I-like domains 3; ECM Glycoprotein; WGCNA $MM = 0.987, GS = 0.826, \log_2\text{FC} = +2.51$.
  4. **`VCAN`** (2 votes: SVM-RFE, RF) — Versican; Proteoglycan; WGCNA $MM = 0.877, GS = 0.849, \log_2\text{FC} = +2.75$.
  5. **`COL3A1`** (2 votes: SVM-RFE, RF) — Collagen Type III Alpha 1; Fibrillar Collagen; WGCNA $MM = 0.813, GS = 0.888, \log_2\text{FC} = +2.84$.
  6. **`COMP`** (2 votes: SVM-RFE, RF) — Cartilage Oligomeric Matrix Protein; ECM Glycoprotein; WGCNA $MM = 0.753, GS = 0.864, \log_2\text{FC} = +4.08$.
  7. **`COL8A1`** (2 votes: SVM-RFE, RF) — Collagen Type VIII Alpha 1; Short-chain Collagen; WGCNA $MM = 0.885, GS = 0.871, \log_2\text{FC} = +2.68$.
  8. **`THBS3`** (2 votes: SVM-RFE, RF) — Thrombospondin 3; Adhesive ECM Glycoprotein; WGCNA $MM = 0.763, GS = 0.834, \log_2\text{FC} = +1.16$.
  9. **`COL11A1`** (2 votes: SVM-RFE, RF) — Collagen Type XI Alpha 1; Minor Fibrillar Collagen; WGCNA $MM = 0.901, GS = 0.867, \log_2\text{FC} = +3.79$.
  10. **`INHBA`** (2 votes: SVM-RFE, RF) — Inhibin Subunit Beta A (Activin A); Secreted Growth Factor; WGCNA $MM = 0.906, GS = 0.878, \log_2\text{FC} = +2.52$.
  11. **`LOX`** (2 votes: SVM-RFE, RF) — Lysyl Oxidase; Covalent Matrix Cross-Linking Enzyme; WGCNA $MM = 0.736, GS = 0.823, \log_2\text{FC} = +2.16$.

* **Associated Data Table**: [results/tables/ML_hub_genes_from_WGCNA_ECM.csv](file:///d:/Peritoneal%20Project/results/tables/ML_hub_genes_from_WGCNA_ECM.csv)
* **Associated Visualizations**:
  * **WGCNA ML Consensus Barchart**: [results/figures/WGCNA_ML_01_consensus_votes_barchart.png](file:///d:/Peritoneal%20Project/results/figures/WGCNA_ML_01_consensus_votes_barchart.png)
  * **WGCNA ML Model Selection Heatmap**: [results/figures/WGCNA_ML_02_model_selection_heatmap.png](file:///d:/Peritoneal%20Project/results/figures/WGCNA_ML_02_model_selection_heatmap.png)
  * **Per-Model Importance (2x2)**: [results/figures/WGCNA_ML_03_per_model_importance_2x2.png](file:///d:/Peritoneal%20Project/results/figures/WGCNA_ML_03_per_model_importance_2x2.png)
  * **Hub Genes Expression Heatmap**: [results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png](file:///d:/Peritoneal%20Project/results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png)
  * **WGCNA 3-Way Convergence Venn Diagram**: [results/figures/venn_wgcna_convergence.png](file:///d:/Peritoneal%20Project/results/figures/venn_wgcna_convergence.png)

---

## 📈 External Cohort Validation in GSE125498

To evaluate exploratory diagnostic utility for discriminating early peritoneal dialysis exposure from progressive membrane fibrosis, prioritized biomarkers were tested in human effluent-derived peritoneal cells from **GSE125498** ($N = 33$ patients):
* **Early Stage (Short-Term PD, SPD: 0–24 Months)**: $n = 20$ patients (preserved membrane transport).
* **Late Stage (Long-Term PD, LPD: $\ge 25$ Months)**: $n = 13$ patients (high solute transport, established fibrotic remodeling).

### Validation Performance Comparison Across Panels

| Biomarker Panel | Screening Methodology | Available Profiled Genes in GSE125498 | Composite ROC-AUC | Validation Performance Summary |
| :--- | :--- | :--- | :---: | :--- |
| **WGCNA-ECM Hub Panel** | **WGCNA + Consensus ML** | **`ISM1`, `FN1`, `VCAN`, `COL3A1`, `COL8A1`, `THBS3`, `LOX`** (7 genes) | **0.873** | **Superior discrimination; `VCAN` significant ($p=0.034, \text{AUC}=0.723$)** |
| 12-Hub Matrisome Panel | DEG $\cap$ Matrisome + ML | `ISM1`, `TGM2`, `MXRA5`, `COL3A1`, `COL5A2`, `POSTN`, `LOX`, `THBS3` (8 genes) | 0.731 | Moderate composite separation |
| 4-Hub Causal Panel | TWMR + Consensus ML | `P4HA2`, `ADAMTS1`, `TNC` (3 genes) | 0.642 | Weak separation; `TNC` validation failure ($\text{AUC}=0.500$) |

### Detailed Performance of Evaluated WGCNA Hub Biomarkers in GSE125498

| Biomarker | Probe ID | Direction in GSE125498 | $\log_2\text{FC}$ (Late vs Early) | Mann-Whitney $U$ | Mann-Whitney $P$-value | Individual AUC-ROC | Clinical Interpretation |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **`VCAN`** | `ILMN_1687301` | DOWN | -0.52 | 72.0 | **0.0341** | **0.723** | **Statistically significant stage discriminator** |
| **`FN1`** | `ILMN_2366463` | UP | +0.32 | 162.5 | 0.2384 | 0.625 | Upregulated in progressive injury |
| **`THBS3`** | `ILMN_1804663` | DOWN | -0.20 | 105.0 | 0.3667 | 0.596 | Moderate negative correlation with dialysis vintage |
| **`COL8A1`** | `ILMN_1685433` | UP | +0.37 | 151.0 | 0.4501 | 0.581 | Progressive collagen accumulation |
| **`COL3A1`** | `ILMN_1773079` | UP | +0.19 | 139.0 | 0.7541 | 0.535 | Progressive interstitial collagen deposition |
| **`ISM1`** | `ILMN_3239288` | UP | +0.03 | 137.0 | 0.8107 | 0.527 | Preserved baseline expression |
| **`LOX`** | `ILMN_1695880` | UP | +0.02 | 129.0 | 0.9853 | 0.496 | Baseline cross-linking signal in cellular fraction |

*(Note: `EDIL3`, `COMP`, `COL11A1`, and `INHBA` lacked corresponding probes on Illumina HumanHT-12 V4.0 / GPL10558).*

---

## 🛠️ Step-by-Step Pipeline Execution Guide

### Task 0: Build Full Normalized GSE62928 Expression Matrix (R)
```bash
Rscript 00_fetch_gse62928_matrix.R
```
* **Function**: Fetches GSE62928 from GEO (`GPL13158`), extracts `exprs()` and `pData()`, maps probes to gene symbols, applies MaxMean probe collapsing, and outputs:
  * `results/tables/GSE62928_full_expression_matrix.csv` (20,940 genes $\times$ 8 samples)
  * `results/tables/GSE62928_sample_metadata.csv` (8 samples, Case/Control labels)

### Task 1: WGCNA Module Detection & Clinical Trait Correlation (R)
```bash
Rscript 02b_wgcna_analysis.R
```
* **Function**: Executes WGCNA signed network analysis, determines soft-threshold power ($\beta = 12, R^2 = 0.809$), performs hierarchical clustering with dynamic tree cut and module merging (`MEDissThres = 0.25`), and correlates module eigengenes with peritoneal fibrosis.
* **Outputs**:
  * `results/tables/wgcna_module_trait_correlation.csv`
  * `results/tables/wgcna_trait_significant_module_genes.csv` (604 genes in Salmon module)
  * `results/figures/WGCNA_00_sample_outlier_dendrogram.png`
  * `results/figures/WGCNA_01_soft_threshold_selection.png`
  * `results/figures/WGCNA_02_gene_dendrogram_modules.png`
  * `results/figures/WGCNA_03_module_trait_heatmap.png`

### Task 2: WGCNA Trait Module $\cap$ 71 ECM-DEGs Convergence (Python)
```bash
python analyze_convergence_wgcna.py
```
* **Function**: Intersects the 604 Salmon module genes with the 71 verified ECM-DEGs, yielding 40 convergent candidates, and generates a 3-way Venn diagram.
* **Outputs**:
  * `convergent_WGCNA_ECM_genes.csv` & `results/tables/convergent_WGCNA_ECM_genes.csv`
  * `venn_wgcna_convergence.png` & `results/figures/venn_wgcna_convergence.png`

### Task 3: Consensus Machine Learning Hub Identification (Python)
```bash
python 05b_ml_hub_gene_identification_wgcna.py
```
* **Function**: Trains LASSO, SVM-RFE, Random Forest, and XGBoost models on GSE62928 normalized expression across the 40 convergent candidates to identify hub genes with $\ge 2/4$ votes.
* **Outputs**:
  * `results/tables/ML_hub_genes_from_WGCNA_ECM.csv`
  * `results/figures/WGCNA_ML_01_consensus_votes_barchart.png`
  * `results/figures/WGCNA_ML_02_model_selection_heatmap.png`
  * `results/figures/WGCNA_ML_03_per_model_importance_2x2.png`
  * `results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png`

### Task 4: External Clinical Validation in GSE125498 (Python)
```bash
python 06_external_validation_GSE125498.py
```
* **Function**: Tests all biomarker panels in the independent GSE125498 patient cohort ($n = 33$), generating Mann-Whitney U test statistics, ROC discrimination curves, stage trajectories, and clustered patient heatmaps.
* **Outputs**:
  * `results/tables/GSE125498_wgcna_hub_validation_metrics.csv`
  * `results/tables/GSE125498_hub_genes_validation_metrics.csv`
  * `results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png`
  * `results/figures/Validation_02_roc_curves_early_vs_late.png`
  * `results/figures/Validation_03_stage_progression_trajectories.png`
  * `results/figures/Validation_04_patient_cohort_heatmap.png`

---

## 📁 Repository Structure & Directory Map

```
d:/Peritoneal Project/
├── 00_fetch_gse62928_matrix.R                     # TASK 0: Full matrix builder & probe mapper
├── 01_load_qc_preprocess.R                       # Primary DEG preprocessing pipeline
├── 02_differential_expression.R                  # Limma differential expression
├── 02b_wgcna_analysis.R                          # TASK 1: WGCNA network & trait correlation
├── 03_matrisome_filtering.R                      # Matrisome masterlist intersection
├── 04_functional_enrichment.R                    # GO/KEGG functional enrichment
├── 05_ml_hub_gene_identification.py              # ML consensus (Original 15 causal & 71 ECM)
├── 05b_ml_hub_gene_identification_wgcna.py       # TASK 3: ML consensus on 40 WGCNA-ECM candidates
├── 06_external_validation_GSE125498.py           # TASK 4: External validation on GSE125498
├── analyze_convergence_wgcna.py                  # TASK 2: WGCNA ∩ ECM-DEG convergence & 3-way Venn
├── data/
│   ├── GSE62928_series_matrix.txt.gz             # GSE62928 series matrix from GEO
│   ├── GSE125498_family.soft.gz                  # GSE125498 SOFT file from GEO
│   └── ECM genes all.xlsx                        # Curated 1,027 Human Matrisome reference
├── results/
│   ├── figures/
│   │   ├── WGCNA_00_sample_outlier_dendrogram.png
│   │   ├── WGCNA_01_soft_threshold_selection.png
│   │   ├── WGCNA_02_gene_dendrogram_modules.png
│   │   ├── WGCNA_03_module_trait_heatmap.png
│   │   ├── venn_wgcna_convergence.png
│   │   ├── WGCNA_ML_01_consensus_votes_barchart.png
│   │   ├── WGCNA_ML_02_model_selection_heatmap.png
│   │   ├── WGCNA_ML_03_per_model_importance_2x2.png
│   │   ├── WGCNA_ML_04_hub_genes_expression_heatmap.png
│   │   ├── Validation_01_hub_genes_mann_whitney_boxplots.png
│   │   ├── Validation_02_roc_curves_early_vs_late.png
│   │   ├── Validation_03_stage_progression_trajectories.png
│   │   └── Validation_04_patient_cohort_heatmap.png
│   └── tables/
│       ├── GSE62928_full_expression_matrix.csv
│       ├── GSE62928_sample_metadata.csv
│       ├── wgcna_module_trait_correlation.csv
│       ├── wgcna_trait_significant_module_genes.csv
│       ├── convergent_WGCNA_ECM_genes.csv
│       ├── ML_hub_genes_from_WGCNA_ECM.csv
│       ├── GSE125498_wgcna_hub_validation_metrics.csv
│       └── GSE125498_hub_genes_validation_metrics.csv
└── README.md                                     # Master documentation
```

---

## ⚖️ Methodological Notes & Limitations

1. **Cohort Sample Sizes**:
   * **Discovery (GSE62928)**: Consists of $N = 8$ human peritoneal tissue samples (4 Encapsulating Peritoneal Sclerosis vs 4 non-fibrotic controls). In small cohort transcriptomics, signed WGCNA with soft power $\beta = 12$ successfully satisfied scale-free topology ($R^2 = 0.809$) and identified robust co-expression modules without sample outliers.
   * **Validation (GSE125498)**: Consists of $N = 33$ human effluent-derived peritoneal cell samples (20 short-term vs 13 long-term PD).
2. **Biological Context**:
   * Discovery samples reflect full-thickness peritoneal membrane tissue biopsies, whereas validation samples reflect cellular components shed into peritoneal dialysis effluent. Despite this biological difference in tissue compartments, the 7-gene WGCNA-ECM composite signature achieved **$\text{AUC} = 0.873$**, underscoring strong translational relevance.
3. **Causal vs Network Shift**:
   * Mendelian Randomization operates under strict instrumental variable assumptions ($Z \to X \to Y$) to deduce lifelong unconfounded genetic causality. In contrast, WGCNA identifies co-regulated gene modules exhibiting strong phenotypic correlation. The WGCNA approach captures active localized tissue pathophysiology, which translates directly into high-accuracy diagnostic biomarker signatures.
