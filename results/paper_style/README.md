# Paper-Style Analysis Workflow & Multi-Panel Figure Manifest

**Date of Run:** 2026-10-14  
**Pipeline Source:** Replicating the diagnostic biomarker analysis workflow and multi-panel figure layout of Sun et al. (*Electronic Journal of Biotechnology* 81 (2026) 100710) for human peritoneal dialysis effluent (GSE125498: 13 Long-term PD [LPD] vs 20 Short-term PD [SPD]).  
**Location:** `results/paper_style/`

---

## 1. Computational Environment & Software Versions

### Python Environment
- **Python Version:** 3.13.1 (64-bit AMD64)
- **Key Python Packages:**
  - `numpy`: 2.2.6
  - `pandas`: 2.2.3
  - `scipy`: 1.18.1
  - `scikit-learn`: 1.6.1
  - `statsmodels`: 0.14.6
  - `matplotlib`: 3.10.0
  - `seaborn`: 0.13.2
  - `umap-learn`: 0.5.12
  - `networkx`: 3.4.2
  - `pillow`: 12.3.0
  - `gseapy`: 1.3.1
  - `xgboost`: 3.1.3

### R / Bioconductor Environment
- **R Version:** 4.4.2 (2024-10-31 ucrt, x86_64-w64-mingw32/x64)
- **Base Packages:** `stats`, `graphics`, `grDevices`, `utils`, `datasets`, `methods`, `base`
- **Bioconductor & Core Packages:**
  - `limma`: 3.62.2
  - `WGCNA`: 1.73
  - `dynamicTreeCut`: 1.63-1
  - `fastcluster`: 1.2.6
  - `clusterProfiler`: 4.14.4
  - `org.Hs.eg.db`: 3.20.0
  - `GEOquery`: 2.74.0
  - `Biobase`: 2.66.0
  - `VennDiagram`: 1.7.3
  - `ggplot2`: 3.5.1

### Databases & Knowledgebases
- **Naba Matrisome Reference:** The Matrisome Project (1,027 human matrisome genes, 7 structural & functional subcategories, `ECM genes all.xlsx`)
- **STRING Database:** STRING v12.0 (Functional protein association network, species *Homo sapiens*, combined interaction score threshold $\ge 0.40$)
- **GeneMANIA:** GeneMANIA API (20 additional interactors, automatic network weighting, dual physical/co-expression/pathway evidence)
- **DSigDB:** Drug Signatures Database (DSigDB v1.0, chemical perturbation signatures, FDR Benjamini-Hochberg adjustment)
- **KEGG:** Kyoto Encyclopedia of Genes and Genomes (*Homo sapiens* pathway database, release 2024/2025)
- **CellMarker 2024 / LM22:** Effluent leukocyte and stromal deconvolution marker gene sets

---

## 2. Master Figure Manifest (High-Resolution 300 DPI PNG & Vector PDF)

| Figure | Script | PNG (300 dpi) | PDF (Vector) | Matching Source Data CSV |
| :--- | :--- | :--- | :--- | :--- |
| **Fig. 1** | `01_Fig1_workflow.py` | `Fig1_workflow.png` / `Fig1.png` | `Fig1_workflow.pdf` / `Fig1.pdf` | `Fig1_source_data.csv` |
| **Fig. 2** | `02_Fig2_differential_expression.py` | `Fig2_differential_expression.png` / `Fig2.png` | `Fig2_differential_expression.pdf` / `Fig2.pdf` | `Fig2A_dimension_reduction_source.csv`, `Fig2_limma_gene_results_unique_symbols.csv`, `Fig2C_top50_heatmap_source.csv`, `Fig2D_KEGG_up_LPD_FDR05.csv`, `Fig2E_KEGG_up_SPD_FDR05.csv` |
| **Fig. 3** | `03_Fig3_matrisome_WGCNA.py` | `Fig3_matrisome_WGCNA.png` / `Fig3.png` | `Fig3_matrisome_WGCNA.pdf` / `Fig3.pdf` | `Fig3A_soft_threshold_source.csv`, `Fig3C_module_trait_source.csv`, `Fig3D_three_way_intersection_counts.csv`, `Fig3E_intersection_enrichment_source.csv`, `Fig3F_saved_STRING_edges_score_ge_0_4.csv` |
| **Fig. 4** | `04_Fig4_diagnostic_genes.py` | `Fig4_diagnostic_genes.png` / `Fig4.png` | `Fig4_diagnostic_genes.pdf` / `Fig4.pdf` | `Fig4_39_gene_pool_selector_results_reused.csv`, `Fig4_39_gene_pool_intersections_reused.csv`, `Fig4E_training_ROC_points.csv`, `Fig4G_discovery_expression_source.csv`, `Fig4_final_panel_provenance.csv` |
| **Fig. 5** | `05_Fig5_nomogram_immune_context.py` | `Fig5_nomogram_immune_context.png` / `Fig5.png` | `Fig5_nomogram_immune_context.pdf` / `Fig5.pdf` | `Fig5A_five_gene_ridge_point_display.csv`, `Fig5B_calibration_bins_source.csv`, `Fig5C_decision_curve_source.csv`, `Fig5D_clinical_impact_source.csv`, `Fig5E_saved_ssgsea_sample_scores.csv`, `Fig5F_immune_immune_spearman_source.csv`, `Fig5G_seven_gene_immune_spearman_source.csv` |
| **Fig. 6** | `06_Fig6_network_drug_pathway.py` | `Fig6_network_drug_pathway.png` / `Fig6.png` | `Fig6_network_drug_pathway.pdf` / `Fig6.pdf` | `Fig6A_GeneMANIA_nodes_source.csv`, `Fig6A_GeneMANIA_edges_source.csv`, `Fig6B_DSigDB_top_terms.csv`, `Fig6C_single_gene_KEGG_GSEA_all_source.csv` |
| **Fig. 7** | `07_Fig7_molecular_subtypes.py` | `Fig7_molecular_subtypes.png` / `Fig7.png` | `Fig7_molecular_subtypes.pdf` / `Fig7.pdf` | `Fig7A_consensus_CDF_source.csv`, `Fig7A_consensus_metrics.csv`, `Fig7A_best_k_consensus_matrix.csv`, `Fig7B_PAC_source.csv`, `Fig7_cluster_assignments_k2.csv`, `Fig7C_panel_expression_by_cluster_source.csv`, `Fig7D_ssgsea_cluster_scores_source.csv` |

---

## 3. Audited Reused vs. Newly Run Components

In accordance with strict reproducibility and data integrity protocols:
1. **Reused Components (Pre-existing Audited Tables):**
   - Differential expression limma top-table (`GSE125498_all_results.csv`, 11,741 unique genes, 813 nominal DEGs, 143 FDR < 0.05).
   - Naba Matrisome overlap (`GSE125498_ECM_intersection.csv`, 44 ECM-DEGs: 31 Up, 13 Down; $1.7145\times$ enrichment, hypergeometric $P = 0.00240$).
   - Signed dynamic WGCNA modules (`wgcna_modules_dynamic_taskB.csv`, 14 modules, soft power $\beta = 12$, top trait module `MEmagenta` $r = 0.480, P = 0.0047$).
   - Machine learning 3-model feature selection results (`ml_39_gene_results.csv`, LASSO 11, SVM-RFE 10, Random Forest 10).
   - 5-gene penalized ridge logistic regression coefficients (`T06_reduced_model_coefficients.csv`: TNFSF15, ADAM19, FLT3LG, NAT6, ECM1).
   - Equal-count quintile calibration and decision curve analysis tables (`T06_discovery_calibration_bins.csv`, `T06_discovery_decision_curve.csv`).
   - Leukocyte ssGSEA deconvolution scores (`immune_celltype_ssgsea_scores.csv`).
   - GeneMANIA functional association network export and DSigDB perturbation signature enrichment tables (`genemania_network_nodes.csv`, `genemania_network_edges.csv`, `dsigdb_enrichment_full.csv`).
   - Single-gene preranked KEGG GSEA results (`kegg_gsea_all_genes.csv`).

2. **Newly Run Components (Created for Paper-Style Reproduction):**
   - Unified multi-panel figure compilation scripts (`01_Fig1_workflow.py` through `07_Fig7_molecular_subtypes.py`).
   - UMAP dimension reduction on the 33 discovery effluent samples (`Fig2A_dimension_reduction_source.csv`).
   - Consensus clustering partition at $k = 2$ for discovery LPD patients ($n = 13$; Subtype C1 $n = 9$, Subtype C2 $n = 4$) and subtype-specific biomarker and immune score comparison (`Fig7_cluster_assignments_k2.csv`).
   - High-resolution (300 DPI) publication-grade rendering matching the exact layouts, fonts, color palettes, and panel arrangements of Sun et al.
