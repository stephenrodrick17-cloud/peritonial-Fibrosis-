# 🔒 FINAL PIPELINE LOCK: Peritoneal Fibrosis Transcriptomics Framework

**Status**: LOCKED & VERIFIED  
**Discovery Cohort**: **GSE125498** (Illumina HumanHT-12 v4, $N=33$ peritoneal effluent cells: 13 Long-term PD cases vs 20 Short-term PD controls)  
**Validation Cohort**: **GSE62928** (Affymetrix HG-U133 Plus 2.0, $N=8$ human parietal peritoneal biopsies: 4 severe EPS cases vs 4 non-EPS controls)  
**Primary Cutoffs**: $|\log_2\text{FC}| \ge 0.585$ (1.5-fold), nominal $P < 0.05$  
**Fixed Seed**: 42 across all feature selectors, cross-validations, and permutations.

---

## 1. Master Pipeline Number Lock

| Phase | Scientific Metric | Exact Value | Verification Table / Source |
| :--- | :--- | :--- | :--- |
| **Phase 1** | Discovery Samples ($N$) | **33** (13 LPD / 20 SPD) | `results/tables/GSE125498_sample_metadata.csv` |
| **Phase 1** | Measured Probe Sets | **19,164** probes | `results/tables/GSE125498_all_results.csv` |
| **Phase 1** | Gene Universe A (Best-Probe-by-$P$) | **11,741** unique genes | `results/tables/GSE125498_all_results.csv` |
| **Phase 1** | Gene Universe B (MaxMean Collapse) | **11,741** unique genes | `results/tables/GSE125498_full_expression_matrix.csv` |
| **Phase 2** | Primary Pro-fibrotic Up-DEGs | **496** genes | `results/tables/GSE125498_DEGs_filtered.csv` |
| **Phase 2** | Primary Down-DEGs | **317** genes | `results/tables/GSE125498_DEGs_filtered.csv` |
| **Phase 2** | Total Primary DEGs | **813** genes (7.0% genome) | `results/tables/GSE125498_DEGs_filtered.csv` |
| **Phase 2** | Genome-wide FDR $< 0.05$ | **143** genes | `results/tables/GSE125498_all_results.csv` |
| **Phase 3** | Naba Reference Matrisome Genes | **1,027** genes (7 subcats) | `ECM genes all.xlsx` |
| **Phase 3** | Measured Matrisome in Discovery | **428** genes | `results/tables/GSE125498_ECM_intersection.csv` |
| **Phase 3** | Pro-fibrotic Up ECM-DEGs | **31** genes | `results/tables/GSE125498_ECM_intersection.csv` |
| **Phase 3** | Down ECM-DEGs | **13** genes | `results/tables/GSE125498_ECM_intersection.csv` |
| **Phase 3** | ECM Fold Enrichment (Up) | **$1.7145\times$** | `results/tables/GSE125498_Matrisome_Overlap_Counts.csv` |
| **Phase 3** | Hypergeometric Over-representation $P$ | **$2.403 \times 10^{-3}$** | `results/tables/phase3_counts.json` |
| **Phase 3** | 10,000 Monte Carlo Permutation $P$ | **$0.0307$** | `results/tables/phase3_counts.json` |
| **Phase 4** | WGCNA Input Matrix | **5,001** genes $\times 33$ samples | `results/tables/wgcna_gene_module_membership.csv` |
| **Phase 4** | Co-expression Modules Cut | **14** modules ($\beta = 12$) | `results/tables/wgcna_module_eigengenes.csv` |
| **Phase 4** | Top Trait-Correlated Module | **`MEmagenta`** ($r=+0.480, P=4.66\times 10^{-3}$) | `results/tables/wgcna_module_trait_correlation.csv` |
| **Phase 4** | Top Module Core Gene Count | **37** genes | `results/tables/wgcna_trait_significant_module_genes.csv` |
| **Phase 5** | ML Candidate Feature Pool | **25** genes | `results/tables/ML_hub_genes_from_WGCNA_ECM_all_results.csv` |
| **Phase 5** | **11 Consensus Hub Genes** ($\ge 2/4$ votes) | **`FLT3LG, TNFSF15, ADAM19, SERPINA10, LTB, CLEC4F, SEMA3E, SEMA4F, EBI3, LTA, CXCL14`** | `results/tables/ML_hub_genes_from_WGCNA_ECM.csv` |
| **Phase 6** | Validation Cohort Size ($N$) | **8** (4 EPS / 4 Non-EPS) | `results/tables/GSE125498_hub_validation_metrics.csv` |
| **Phase 6** | Mapped Hubs in GSE62928 | **11 / 11** (100% mapped) | `results/tables/GSE62928_hub_probe_mapping.csv` |
| **Phase 6** | Apparent In-Sample AUC | **1.0000** | `results/tables/GSE125498_hub_validation_metrics.csv` |
| **Phase 6** | Leave-One-Out CV AUC (LOOCV) | **1.0000** | `results/tables/GSE125498_hub_validation_metrics.csv` |
| **Phase 6** | 50-Replicate 2-Fold CV AUC | **$0.7575 \pm 0.2585$** | `results/tables/GSE125498_hub_validation_metrics.csv` |
| **Phase 7** | STRING PPI Interacting Edges | **2** edges (3.6% connectivity) | `results/tables/STRING_edges_functional.csv` |
| **Phase 7** | 5-Gene Nomogram Apparent C-index | **0.9808** | `results/tables/decision_curve_analysis_metrics.csv` |
| **Phase 7** | 5-Gene Nomogram 5-Fold CV C-index | **0.8731** (Optimism = 0.1077) | `results/tables/nomogram_calibration_curve_metrics.csv` |
| **Phase 7** | Myofibroblast Deconvolution $P$ | **MWU $P = 0.3632$** | `results/tables/immune_infiltration_scores.csv` |

---

## 2. Master Figure Manifest (Stored in `results/figures/`)

1. **`GSE125498_01_volcano_plot.png`** / **`.pdf`**: Discovery volcano plot (496 Up / 317 Down DEGs).
2. **`GSE125498_02_matrisome_enrichment.png`** / **`.pdf`**: Matrisome subcategory distribution ($1.71\times$ enrichment, $P = 0.0024$).
3. **`GSE125498_03_wgcna_module_trait.png`** / **`.pdf`**: Signed WGCNA dendrogram and module-trait correlations.
4. **`GSE125498_04_ml_consensus_hubs.png`** / **`.pdf`**: 4-algorithm ML feature selection consensus barchart and vote heatmap.
5. **`GSE125498_05_external_validation_roc.png`** / **`.pdf`**: External validation ROC curves on biopsy cohort GSE62928.
6. **`GSE125498_06_functional_networks_composite.png`** / **`.pdf`**: Multi-panel composite (STRING PPI, GSEA EMT, Nomogram DCA, Myofibroblast score).
