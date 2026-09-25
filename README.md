# Peritoneal Dialysis-Associated Peritoneal Fibrosis: Multi-Omics & Causal Discovery Pipeline

A comprehensive bioinformatics and genetic epidemiology pipeline integrating bulk transcriptomics, extracellular matrix (ECM) matrisome biology, cross-cohort validation, and transcriptome-wide two-sample Mendelian Randomization (TWMR) to identify causal drivers and therapeutic biomarkers of peritoneal dialysis (PD)-induced peritoneal fibrosis and ultrafiltration failure.

---

## 📌 Project Overview & Biological Rationale

Long-term peritoneal dialysis (LPD) is frequently complicated by peritoneal membrane injury, chronic inflammation, epithelial-to-mesenchymal transition (EMT) of mesothelial cells, and progressive peritoneal fibrosis, ultimately culminating in ultrafiltration failure or encapsulating peritoneal sclerosis (EPS). 

This pipeline provides a rigorous, multi-tiered framework to discover, validate, and prioritize causal genes:
1. **Discovery Transcriptomics (GSE125498)**: Bulk expression profiling of human effluent-derived cells comparing Long-term PD (LPD, $\ge 25$ months) vs Short-term PD (SPD, $0-24$ months).
2. **Validation Cohort (GSE62928)**: Independent validation in peritoneal biopsies / effluent mesothelial cells.
3. **ECM Matrisome In Silico Intersection**: Benchmarking against the Naba et al. Human Matrisome database to isolate core structural proteins, glycoproteins, proteoglycans, regulators, and secreted factors.
4. **Transcriptome-Wide Two-Sample Mendelian Randomization (TWMR)**: Harnessing large-scale genetic instruments (eQTLGen, $n=31,684$) and kidney function GWAS (CKDGen eGFR, $n=567,460$) to establish causal directionality and protect against reverse causation and environmental confounding.

---

## 🔬 Pipeline Workflow & Phases Completed

```
   ┌─────────────────────────────────────────────────────────────┐
   │                PHASE 1: DISCOVERY (GSE125498)               │
   │  - 33 effluent cell samples (13 LPD vs 20 SPD, GPL10558)    │
   │  - QC, MaxMean probe collapsing (11,741 unique genes)       │
   │  - Limma Empirical Bayes DEG analysis                       │
   └──────────────────────────────┬──────────────────────────────┘
                                  ▼
   ┌─────────────────────────────────────────────────────────────┐
   │            PHASE 2: VALIDATION & ECM CONVERGENCE            │
   │  - GSE62928 validation cohort                               │
   │  - Naba et al. Human Matrisome Reference (706-1,027 genes)  │
   │  - Identification of 86 convergent ECM-DEGs                │
   └──────────────────────────────┬──────────────────────────────┘
                                  ▼
   ┌─────────────────────────────────────────────────────────────┐
   │       PHASE 3: TRANSCRIPTOME-WIDE MENDELIAN RANDOMIZATION   │
   │  - Exposure: eQTLGen whole-blood cis-eQTLs (n=31,684)       │
   │  - Outcome: CKDGen eGFR GWAS summary statistics (n=567,460) │
   │  - Methods: IVW, Wald Ratio, MR-Egger, Weighted Median      │
   │  - QC: F-stat > 10, Cochran's Q, MR-Egger pleiotropy test   │
   └──────────────────────────────┬──────────────────────────────┘
                                  ▼
   ┌─────────────────────────────────────────────────────────────┐
   │           TRIANGULATED CAUSAL ECM-FIBROSIS TARGETS          │
   │  Top candidates: P4HA2, BMP6, IGFBP3, WNT11, COL4A2,        │
   │  ADAMTS1, LTBP4, TNC, VCAN, ADAM19, SEMA3E                  │
   └─────────────────────────────────────────────────────────────┘
```

---

## 📊 Key Findings & Discovery Metrics

### 1. Discovery Transcriptomics (GSE125498)
- **Cohort**: 13 Long-term PD (LPD, $\ge 25$ mo) vs 20 Short-term PD (SPD, $0-24$ mo).
- **Probes & Genes**: 19,164 probes mapped $\rightarrow$ 11,741 unique HGNC genes via **Maximum Mean Expression (`MaxMean`)** rule.
- **Significant DEGs** ($|\log_2\text{FC}| > 1.0, \text{FDR} < 0.05$): **41 genes** (39 upregulated in LPD, 2 downregulated).
- **Core ECM-DEGs (Strict)**: 
  - **`ADAM19`** ($\log_2\text{FC} = +1.213, \text{FDR} = 0.0277$, ECM Regulator / Metalloproteinase)
  - **`SEMA3E`** ($\log_2\text{FC} = +1.046, \text{FDR} = 0.0327$, Core Matrisome Glycoprotein)
- **Borderline ECM Candidates** ($P < 0.05$): **31 genes**, including canonical peritoneal fibrosis regulators `MMP14` ($P=0.0031$), `TIMP1` ($P=0.0084$), `CTGF`/`CCN2` ($P=0.0120$), `COL4A1` ($P=0.0185$), and `LOXL2` ($P=0.0241$).

### 2. Validation & Matrisome Convergence (GSE62928)
- **Validation DEGs** ($P < 0.05, \log_2\text{FC} \ge 0.585$): 559 genes.
- **Matrisome Overlap**: **86 convergent ECM-DEGs** overlapping between GSE62928 and the curated Human Matrisome master list.

### 3. Transcriptome-Wide Mendelian Randomization (TWMR)
- **Exposure**: eQTLGen Consortium ($n=31,684$ individuals, genome-wide cis-eQTLs $P < 5 \times 10^{-8}$ and relaxed $P < 10^{-5}$).
- **Outcome**: CKDGen eGFR GWAS (Wuttke et al. 2019, $n=567,460$ European ancestry individuals).
- **Harmonization & QC**: Automatic strand alignment, removal of ambiguous palindromic SNPs (MAF $0.42-0.58$), LD clumping ($500\text{ kb}$ window), and instrument strength verification ($F\text{-statistic} > 10$).
- **Triangulated Causal Candidates**:
  - **`P4HA2`** (Prolyl 4-hydroxylase alpha II): $\text{MR Beta} = -0.0176, P = 2.37 \times 10^{-9}, \text{FDR} = 3.61 \times 10^{-7}$ ($\log_2\text{FC} = +1.22$). Essential for collagen triple-helix stabilization and peritoneal collagen deposition.
  - **`BMP6`** (Bone morphogenetic protein 6): $\text{MR Beta} = -0.0102, P = 4.79 \times 10^{-4}, \text{FDR} = 0.0129$ ($\log_2\text{FC} = -1.48$).
  - **`IGFBP3`** (IGF binding protein 3): $\text{MR Beta} = +0.0041, P = 1.08 \times 10^{-3}, \text{FDR} = 0.0234$ ($\log_2\text{FC} = -1.08$).
  - **`WNT11`**, **`COL4A2`**, **`ADAMTS1`**, **`LTBP4`**, **`TNC`**, **`VCAN`**, **`FGL2`**, **`ADAM28`**.

---

## 📁 Repository Structure

```
.
├── 01_load_qc_preprocess.R          # GSE125498 Data loading, metadata QC, PCA, MaxMean collapsing
├── 02_differential_expression.R      # Limma linear modeling, contrast fit (LPD - SPD), Volcano & Heatmap
├── 03_matrisome_filtering.R          # Naba et al. Human Matrisome intersection, Venn & ECM Volcano plots
├── 04_functional_enrichment.R        # clusterProfiler GO (BP, CC, MF) and KEGG pathway enrichment
├── run_pipeline.R                    # Master execution orchestrator for Phase 1
│
├── analyze_convergence.py           # Cross-dataset convergence between GSE62928 & Matrisome
├── data_fetch.py                     # Data loading, eQTLGen & CKDGen GWAS parsing and harmonizer
├── mr_stats.py                       # High-precision MR statistical estimators (IVW, Egger, Median, Q, F-stat)
├── mr_plots.py                       # Publication-grade MR visualization suite (Scatter, Forest, Funnel, LOO)
├── run_mr_pipeline.py                # Pipeline driver for candidate-focused Mendelian Randomization
├── transcriptome_wide_mr.py          # High-throughput genome-wide/transcriptome-wide TWMR screen
│
├── convergent_ECM_DEGs_nominal.csv   # 86 Convergent ECM-DEGs from validation cohort
├── mr_significant_genes.csv          # Transcriptome-wide MR significant causal genes (FDR < 0.05)
├── mr_deg_ecm_intersection.csv       # Multi-omics triangulated candidate table
├── mr_transcriptome_wide_all_results.csv # Complete TWMR statistical results
│
├── venn_diagram_convergence.png      # Venn diagram: GSE62928 DEGs vs Human Matrisome
├── venn_mr_deg_convergence.png       # Venn diagram: TWMR causal genes vs Peritoneal DEGs
│
├── results/                          # Discovery phase outputs
│   ├── figures/                      # High-resolution PDF and PNG plots
│   └── tables/                       # Statistical summary CSV tables
└── README.md                         # Project documentation and progress report
```

---

## 🛠️ Methodological Standards & Justifications

1. **Probe Collapsing Strategy (`MaxMean`)**:
   - For genes targeted by multiple Illumina probes, the probe with the highest mean hybridization signal across samples was retained. This avoids dilution of real biological signal by non-specific or low-affinity probes and maximizes signal-to-noise ratio.
2. **Harmonization & Allele Matching in MR**:
   - Effect alleles and effect directions across eQTLGen and CKDGen GWAS are harmonized. Z-scores are mapped to standardized effect sizes using sample size $N$ and empirical allele frequencies. Palindromic variants with intermediate frequencies are discarded to eliminate strand ambiguity.
3. **Piotropy and Heterogeneity Diagnostics**:
   - Multi-instrument genes are evaluated via Cochran's $Q$ test for heterogeneity and MR-Egger intercept test for directional horizontal pleiotropy. Weak instruments ($F < 10$) are automatically flagged.

---

## 🚀 Getting Started & Execution

### Prerequisites
- **R** ($\ge 4.4.0$) with Bioconductor packages: `GEOquery`, `limma`, `Biobase`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `pheatmap`, `ggplot2`, `ggrepel`, `dplyr`.
- **Python** ($\ge 3.9$) with: `numpy`, `pandas`, `scipy`, `matplotlib`, `matplotlib-venn`, `openpyxl`.

### 1. Run Discovery Bulk Transcriptomics (R)
```bash
Rscript run_pipeline.R
```

### 2. Run Convergence Analysis (Python)
```bash
python analyze_convergence.py
```

### 3. Run Transcriptome-Wide Mendelian Randomization (Python)
```bash
python transcriptome_wide_mr.py
```

---

## 📈 Next Steps

- [ ] **WGCNA Co-expression Network Analysis**: Construct unsigned/signed co-expression modules in GSE125498 to identify ECM-enriched hub genes correlated with PD duration and peritoneal solute transport rate (PSTR).
- [ ] **Machine Learning Feature Selection**: Implement LASSO, Random Forest (Boruta), Support Vector Machine Recursive Feature Elimination (SVM-RFE), and XGBoost to extract core biomarker panels.
- [ ] **Independent Validation**: Validate candidate panel expression in external cohorts (GSE62928, single-cell/spatial datasets) and experimental models of PD fibrosis.

---

## 📜 References
1. **GSE125498**: Bulk transcriptomics of peritoneal effluent cells in long-term vs short-term peritoneal dialysis.
2. **Naba et al.**: The Matrisome: in silico definition and in vivo characterization by proteomics of normal and diseased extracellular matrices. *Matrix Biol* (2012).
3. **eQTLGen Consortium**: Large-scale cis- and trans-eQTL mapping in 31,684 individuals. *Nat Genet* (2021).
4. **CKDGen Consortium (Wuttke et al.)**: A catalog of genetic loci associated with kidney function from analyses of a million individuals. *Nat Genet* (2019).
