# Comprehensive Transcriptomic & ECM Matrisome Analysis of Peritoneal Fibrosis
## Discovery Cohort GSE125498 with Naba Human Matrisome Integration

### Executive Summary
- **Discovery Cohort**: GSE125498 (*Homo sapiens*, Illumina HumanHT-12 V4.0) comprising **33 samples** (13 Long-term PD [LPD, high fibrotic risk] vs. 20 Short-term PD [SPD, low fibrotic risk controls]).
- **Significance Thresholds**: Strictly **|log₂FC| ≥ 0.585** (1.5-fold change) and **nominal P < 0.05** (FDR correction omitted per project specification).
- **Total Identified DEGs**: **813 genes** (496 Upregulated, 317 Downregulated) from 11,741 collapsed gene features.
- **Naba Matrisome Intersection**: Of 1,027 reference human matrisome genes (428 measured in the dataset), exactly **44 DEG-Matrisomal genes** were identified (**31 pro-fibrotic Upregulated**, **13 Downregulated**).
- **ECM Enrichment Significance**: Hypergeometric test *P* = 0.002403 (fold enrichment = 1.71x; 10,000-iteration permutation *P* = 0.0307).
- **Consensus Machine Learning Hub Genes**: **11 Hub Genes** prioritized by multi-algorithmic ensemble (LASSO, SVM-RFE, Random Forest, XGBoost): **FLT3LG, TNFSF15, ADAM19, SERPINA10, LTB, CLEC4F, SEMA3E, SEMA4F, EBI3, LTA, CXCL14**.
- **External Validation**: Evaluated on GSE62928 (4 Encapsulating Peritoneal Sclerosis [EPS] vs. 4 Controls), achieving **apparent AUC = 1.000** and **LOOCV AUC = 1.000**.

---

### Table 1: Complete DEG-Matrisomal Gene List (44 Genes)

#### A. Pro-Fibrotic Upregulated ECM DEGs (31 Genes)
| Gene Symbol | log₂FC | Nominal P-Value | Matrisome Division | Matrisome Category | Hub Status |
|:---|:---:|:---:|:---|:---|:---:|
| **PI3** | +1.706 | 7.22e-03 | Matrisome-associated | ECM Regulators | ECM DEG |
| **LTB** | +1.675 | 3.33e-06 | Matrisome-associated | Secreted Factors | ★ Consensus Hub |
| **SFRP4** | +1.458 | 3.35e-02 | Matrisome-associated | Secreted Factors | ECM DEG |
| **FLT3LG** | +1.438 | 2.35e-05 | Matrisome-associated | Secreted Factors | ★ Consensus Hub |
| **ANGPTL4** | +1.402 | 7.25e-03 | Matrisome-associated | Secreted Factors | ECM DEG |
| **ADAM19** | +1.213 | 1.33e-04 | Matrisome-associated | ECM Regulators | ★ Consensus Hub |
| **SERPINA10** | +1.166 | 6.53e-04 | Matrisome-associated | ECM Regulators | ★ Consensus Hub |
| **CLEC4F** | +1.101 | 1.43e-02 | Matrisome-associated | ECM-affiliated Proteins | ★ Consensus Hub |
| **SEMA3E** | +1.046 | 2.81e-04 | Matrisome-associated | ECM-affiliated Proteins | ★ Consensus Hub |
| **ECM1** | +1.034 | 7.01e-03 | Core matrisome | ECM Glycoproteins | ECM DEG |
| **PLXNA3** | +1.017 | 7.28e-04 | Matrisome-associated | ECM-affiliated Proteins | ECM DEG |
| **CST7** | +1.012 | 5.00e-04 | Matrisome-associated | ECM Regulators | ECM DEG |
| **PIK3IP1** | +0.996 | 5.95e-04 | Matrisome-associated | Secreted Factors | ECM DEG |
| **HAPLN3** | +0.982 | 6.73e-03 | Core matrisome | Proteoglycans | ECM DEG |
| **TNFSF15** | +0.974 | 2.73e-05 | Matrisome-associated | Secreted Factors | ★ Consensus Hub |
| **SEMA4F** | +0.925 | 1.47e-03 | Matrisome-associated | ECM-affiliated Proteins | ★ Consensus Hub |
| **C1QB** | +0.923 | 2.82e-03 | Matrisome-associated | ECM-affiliated Proteins | ECM DEG |
| **COL9A2** | +0.904 | 3.04e-02 | Core matrisome | Collagens | ECM DEG |
| **SPP1** | +0.850 | 3.07e-02 | Core matrisome | ECM Glycoproteins | ECM DEG |
| **EBI3** | +0.803 | 2.44e-03 | Matrisome-associated | Secreted Factors | ★ Consensus Hub |
| **TNFAIP6** | +0.800 | 3.93e-02 | Core matrisome | ECM Glycoproteins | ECM DEG |
| **SPOCK2** | +0.765 | 2.81e-03 | Core matrisome | Proteoglycans | ECM DEG |
| **TNFSF4** | +0.750 | 3.56e-03 | Matrisome-associated | Secreted Factors | ECM DEG |
| **COL8A1** | +0.749 | 4.88e-02 | Core matrisome | Collagens | ECM DEG |
| **CLEC2D** | +0.709 | 2.65e-02 | Matrisome-associated | ECM-affiliated Proteins | ECM DEG |
| **CCL5** | +0.689 | 3.96e-02 | Matrisome-associated | Secreted Factors | ECM DEG |
| **TGM3** | +0.665 | 2.44e-02 | Matrisome-associated | ECM Regulators | ECM DEG |
| **C1QA** | +0.656 | 1.84e-02 | Matrisome-associated | ECM-affiliated Proteins | ECM DEG |
| **LTA** | +0.642 | 1.46e-02 | Matrisome-associated | Secreted Factors | ★ Consensus Hub |
| **CXCL14** | +0.638 | 5.03e-03 | Matrisome-associated | Secreted Factors | ★ Consensus Hub |
| **CCBE1** | +0.605 | 5.94e-03 | Matrisome-associated | Secreted Factors | ECM DEG |

#### B. Downregulated ECM DEGs (13 Genes)
| Gene Symbol | log₂FC | Nominal P-Value | Matrisome Division | Matrisome Category | Hub Status |
|:---|:---:|:---:|:---|:---|:---:|
| **ITLN1** | -1.204 | 1.39e-02 | Matrisome-associated | ECM-affiliated Proteins | ECM DEG |
| **ADAM10** | -1.136 | 1.07e-02 | Matrisome-associated | ECM Regulators | ECM DEG |
| **SBSPON** | -1.022 | 9.34e-03 | Core matrisome | ECM Glycoproteins | ECM DEG |
| **ITLN2** | -0.958 | 3.03e-02 | Matrisome-associated | ECM-affiliated Proteins | ECM DEG |
| **CXCL6** | -0.870 | 3.97e-02 | Matrisome-associated | Secreted Factors | ECM DEG |
| **BMP7** | -0.837 | 1.28e-02 | Matrisome-associated | Secreted Factors | ECM DEG |
| **CLEC12A** | -0.756 | 2.84e-03 | Matrisome-associated | ECM-affiliated Proteins | ECM DEG |
| **OGN** | -0.736 | 3.28e-02 | Core matrisome | Proteoglycans | ECM DEG |
| **S100Z** | -0.724 | 6.29e-03 | Matrisome-associated | Secreted Factors | ECM DEG |
| **SERPINB6** | -0.647 | 7.12e-03 | Matrisome-associated | ECM Regulators | ECM DEG |
| **CST6** | -0.644 | 3.82e-02 | Matrisome-associated | ECM Regulators | ECM DEG |
| **TNFSF8** | -0.639 | 2.72e-02 | Matrisome-associated | Secreted Factors | ECM DEG |
| **LAMB1** | -0.605 | 4.31e-02 | Core matrisome | ECM Glycoproteins | ECM DEG |

---

### Table 2: 11 Consensus Machine Learning Hub Genes
| Gene | log₂FC | P-Value | WGCNA Module | ML Votes (out of 4) | PPI Degree | Key Biological Function in Peritoneal Fibrosis |
|:---|:---:|:---:|:---:|:---:|:---:|:---|
| **FLT3LG** | +1.438 | 2.35e-05 | black | 4/4 | 0 | Fms-related tyrosine kinase 3 ligand; stimulates dendritic cell mobilization and fibro-inflammatory crosstalk. |
| **TNFSF15** | +0.974 | 2.73e-05 | black | 4/4 | 1 | TNF superfamily member 15 (TL1A); pro-inflammatory cytokine driving mucosal and peritoneal fibrogenesis. |
| **ADAM19** | +1.213 | 1.33e-04 | black | 3/4 | 0 | Disintegrin and metalloproteinase domain 19; ECM sheddase regulating shedding of cytokines and matrix degradation. |
| **SERPINA10** | +1.166 | 6.53e-04 | black | 3/4 | 0 | Serpin family A member 10 (ZPI); serine protease inhibitor regulating coagulation and tissue remodeling. |
| **LTB** | +1.675 | 3.33e-06 | black | 2/4 | 2 | Lymphotoxin beta; membrane anchor for LT-alpha, promoting tertiary lymphoid structure formation in peritoneal tissue. |
| **CLEC4F** | +1.101 | 1.43e-02 | black | 2/4 | 0 | C-type lectin domain family 4 member F; macrophage surface lectin mediating cellular adhesion and immune clearance. |
| **SEMA3E** | +1.046 | 2.81e-04 | black | 2/4 | 0 | Semaphorin 3E; axon guidance factor modulating vascular patterning and endothelial-mesenchymal transition. |
| **SEMA4F** | +0.925 | 1.47e-03 | black | 2/4 | 0 | Semaphorin 4F; neural and mesenchymal guidance factor implicated in cellular migration and fibrotic architecture. |
| **EBI3** | +0.803 | 2.44e-03 | black | 2/4 | 0 | Epstein-Barr virus induced 3; IL-27 / IL-35 cytokine subunit regulating chronic peritoneal inflammatory signaling. |
| **LTA** | +0.642 | 1.46e-02 | black | 2/4 | 1 | Lymphotoxin alpha (TNF-beta); pro-inflammatory cytokine activating NF-kB pathway and collagen synthesis. |
| **CXCL14** | +0.638 | 5.03e-03 | black | 2/4 | 0 | C-X-C motif chemokine ligand 14 (BRAK); chemoattractant promoting myofibroblast migration and ECM deposition. |

---

### Methodological Pipeline Summary
1. **Quality Control & Normalization**: GSE125498 normalized expression profile processed via Illumina annotation mapping; redundant probes collapsed via maximum mean expression into 11,741 unique gene symbols.
2. **Differential Expression Analysis**: Empirical Bayes moderated t-tests applied using strict nominal thresholds (|log₂FC| ≥ 0.585, P < 0.05).
3. **Matrisome Curation**: Cross-referenced against the Naba Human Matrisome Masterlist (1,027 genes), stratifying by Core Matrisome (Collagens, Glycoproteins, Proteoglycans) and Matrisome-Associated factors (ECM Regulators, ECM-affiliated, Secreted Factors).
4. **WGCNA Module Detection**: Unsupervised weighted co-expression network constructed across 5,001 top variable genes, yielding 14 co-expression modules. Module-trait correlation identified key stage-associated clusters.
5. **Ensemble Machine Learning**: 25 candidate genes evaluated via LASSO (L1 regularization), SVM-RFE (recursive feature elimination), Random Forest (Gini impurity), and XGBoost (gradient boosting trees). Consensus hubs defined by ≥ 2 algorithm votes.
6. **External Cross-Cohort Validation**: External validation in GSE62928 (EPS vs. controls) yielded perfect discriminatory performance (AUC = 1.000).
7. **Translational Nomogram**: Multivariate logistic regression scoring model with 10-fold cross-validation demonstrated high predictive capability (C-index = 0.873).

Generated automatically by the Peritoneal Fibrosis Transcriptomics Pipeline.