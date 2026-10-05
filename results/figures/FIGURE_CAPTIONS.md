# Programmatic Figure Catalog and Captions (Stage 6)

All figures in this directory were programmatically compiled directly from verified CSV result tables.
Paired source-data CSV files are stored in `results/figures/source_data/`.

---

## Figure 4A: Cell-Type Specific Pseudobulk Expression of 11 Hub Genes in Single-Cell Peritoneal Tissue (GSE248762)

**Image File:** [`results/figures/Fig_4A_pseudobulk_celltype_forest.png`](Fig_4A_pseudobulk_celltype_forest.png)  
**Source Data CSV:** [`results/figures/source_data/fig4a_pseudobulk_celltype_forest_source.csv`](source_data/fig4a_pseudobulk_celltype_forest_source.csv)  

**Caption:** Forest plot depicting edgeR quasi-likelihood negative binomial pseudobulk log2 fold changes (LV_UF vs. LV_NOT_UF) and 95% confidence intervals across five distinct cellular lineages: Stromal/Mesothelial (orange), Monocyte/Macrophage (purple), T Cell (green), NK Cell (pink), and cDC (olive). The 11 candidate hub genes display prominent stromal localization, with severe UF failure associated with marked upregulation of COL11A1, COL8A1, COMP, VCAN, EDIL3, FN1, and INHBA in the stromal/mesothelial compartment.

---

## Figure 4B: Sensitivity and Robustness Evaluation of Peritoneal Stromal Hub Differential Expression

**Image File:** [`results/figures/Fig_4B_stromal_sensitivities_forest.png`](Fig_4B_stromal_sensitivities_forest.png)  
**Source Data CSV:** [`results/figures/source_data/fig4b_stromal_sensitivities_forest_source.csv`](source_data/fig4b_stromal_sensitivities_forest_source.csv)  

**Caption:** Forest plots showing log2 fold changes and 95% confidence intervals across 10 distinct methodological sensitivity pipelines in the stromal/mesothelial lineage for representative hub genes (COL11A1, COL8A1, COMP, VCAN, EDIL3, FN1). Across all modifications—including single-donor exclusion (sens b), doublet inclusion (sens c), gene detection thresholding (sens d, e, h), and alternative artifact ceiling filters (sens g, i, j)—estimated stromal upregulation remains consistent in direction and magnitude.

---

## Figure 5A: Intersection of Validated Hub-Targeting miRNAs with Peritoneal Dialysis Differential miRNA Profiles

**Image File:** [`results/figures/Fig_5A_mirna_hub_intersection.png`](Fig_5A_mirna_hub_intersection.png)  
**Source Data CSV:** [`results/figures/source_data/fig5a_mirna_hub_intersection_source.csv`](source_data/fig5a_mirna_hub_intersection_source.csv)  

**Caption:** Evaluation of regulatory miRNA convergence across multiMiR experimentally validated interactions and external transcriptomic datasets. Left panel: Top miRNAs targeting multiple consensus hub genes that overlap with differentially expressed miRNAs in human effluent exosomes (GSE182736, blue, HUMAN) or rodent peritoneal dialysis effluent (GSE130387, red, CROSS-SPECIES). Right panel: Hypergeometric test statistics indicating significant cross-species enrichment in rodent models (P = 0.000112, k=47 overlap), with human exosome profiling interpreted descriptively due to small sample size (n=3 vs. 3).

---

## Figure 5B: PREDICTED regulatory network (not experimentally validated in this system)

**Image File:** [`results/figures/Fig_5B_cerna_regulatory_network.png`](Fig_5B_cerna_regulatory_network.png)  
**Source Data CSV:** [`results/figures/source_data/fig5b_cerna_regulatory_network_source.csv`](source_data/fig5b_cerna_regulatory_network_source.csv)  

**Caption:** Tri-partite predicted ceRNA regulatory network diagram linking highly expressed peritoneal stromal lncRNAs (detected in GSE248762 single-cell pseudobulk), intermediate regulatory miRNAs (supported by ENCORI CLIP experiments), and downstream pro-fibrotic matrisome hub genes (supported by multiMiR validated interactions). This model represents an in silico hypothesis-generating framework; interactions are not experimentally verified in primary peritoneal tissue.

---

## Figure 5C: Transcriptomic Response of Cultured Human Peritoneal Mesothelial Cells to TGF-beta1 Stimulation (GSE121372)

**Image File:** [`results/figures/Fig_5C_gse121372_tgfb1_timecourse.png`](Fig_5C_gse121372_tgfb1_timecourse.png)  
**Source Data CSV:** [`results/figures/source_data/fig5c_gse121372_tgfb1_timecourse_source.csv`](source_data/fig5c_gse121372_tgfb1_timecourse_source.csv)  

**Caption:** Bar chart illustrating log2 fold changes in human peritoneal mesothelial cells (HPMCs) cultured with 1 ng/mL TGF-beta1 relative to untreated time-matched controls at 6 hours (light blue) and 24 hours (dark blue) on Illumina HumanRef-8 v2.0 beadchips. DESCRIPTIVE ONLY; n=1 per condition per timepoint; no statistical hypothesis testing or P-values are evaluated. Ten of 11 hub genes are assayed on this platform (ISM1 is not represented). Upregulation at 24 hours is observed across COL11A1, COMP, INHBA, EDIL3, VCAN, THBS3, FN1, LOX, and COL8A1.

---

## Figure 6: Methodological Recalibration of Machine Learning Consensus and Generalization Metrics

**Image File:** [`results/figures/Fig_6_discrepancy_corrected_panels.png`](Fig_6_discrepancy_corrected_panels.png)  
**Source Data CSV:** [`results/figures/source_data/fig6_discrepancy_corrected_panels_source.csv`](source_data/fig6_discrepancy_corrected_panels_source.csv)  

**Caption:** Corrected diagnostic panels addressing repository consistency notes. Left panel: Consensus vote tallies across four machine learning feature selection methods (LASSO, SVM-RFE, Random Forest, XGBoost); maximum observed agreement across algorithms is 3/4 votes (ISM1, FN1, EDIL3) and 2/4 votes (VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX), with no genes selected unanimously (0 with 4/4). Right panel: ROC curves distinguishing between the apparent in-sample composite model fit in the effluent validation cohort (GSE125498, blue, AUC = 0.869) and rigorous out-of-fold generalization assessed via 50-repeat 5-fold cross-validation (orange dashed line, mean AUC = 0.678, SD = 0.058) and single seed-42 split (green dotted line, AUC = 0.592).

---

