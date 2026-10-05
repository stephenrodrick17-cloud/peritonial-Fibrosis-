# Programmatic Figure Catalog and Captions (Stage 6)

All figures in this catalog were programmatically compiled directly from verified CSV result tables.
Every value, statistic, and confidence interval is authentic and synchronized with the underlying analysis tables.
Paired source-data CSV files are stored in `results/figures/source_data/`.

---

## Figure 4A: Cell-Type-Specific Pseudobulk Expression of 11 Hub Genes in Single-Cell Peritoneal Tissue (GSE248762)

**Image File:** [`results/figures/Fig_4A_pseudobulk_celltype_forest.png`](Fig_4A_pseudobulk_celltype_forest.png)  
**Source Data CSV:** [`results/figures/source_data/fig4a_pseudobulk_celltype_forest_source.csv`](source_data/fig4a_pseudobulk_celltype_forest_source.csv)  

**Caption:** Forest plot depicting edgeR quasi-likelihood negative binomial pseudobulk log2 fold changes (LV_UF vs. LV_NOT_UF) and 95% confidence intervals across five distinct cellular lineages: Stromal/Mesothelial (orange), Monocyte/Macrophage (purple), T Cell (green), NK Cell (pink), and cDC (olive). Solid circles indicate genes passing expression threshold filters (PASS); open diamonds indicate insufficient counts for pseudobulk testing. In the peritoneal stromal/mesothelial lineage, 8 of 9 testable candidate hub genes exhibit directionally lower expression in severe ultrafiltration failure (negative log2FC: COMP -4.04, FN1 -2.96, COL8A1 -2.66, EDIL3 -2.61, COL3A1 -1.96, LOX -1.72, THBS3 -1.57, VCAN -0.65). Only INHBA displays a positive point estimate (+1.08; 95% CI: -0.59 to +2.75, crossing zero). None of the 9 testable genes achieve statistical significance after genome-wide Benjamini-Hochberg false discovery rate correction (all stromal FDR >= 0.131; COL8A1 FDR=0.131, THBS3 FDR=0.215, FN1/EDIL3/COMP FDR=0.259, COL3A1 FDR=0.335, LOX FDR=0.353, INHBA FDR=0.451, VCAN FDR=0.748). Findings are exploratory and hypothesis-generating.

---

## Figure 4B: Stromal Sensitivity and Robustness Evaluation Across 9 Quality Control and Analytical Pipelines

**Image File:** [`results/figures/Fig_4B_stromal_sensitivities_forest.png`](Fig_4B_stromal_sensitivities_forest.png)  
**Source Data CSV:** [`results/figures/source_data/fig4b_stromal_sensitivities_forest_source.csv`](source_data/fig4b_stromal_sensitivities_forest_source.csv)  

**Caption:** Multi-panel forest plots displaying edgeR quasi-likelihood pseudobulk log2 fold changes and 95% confidence intervals in the peritoneal stromal/mesothelial lineage for all 9 testable consensus hub genes across the primary pipeline (red squares), standard sensitivity models (blue circles: sens b [drop single donor LV_UF-3], sens c [keep doublets], sens d [>=700 genes], sens e [>=50 cells], sens f [top 2000 HVGs], sens h [strict filter]), and technical ceiling models (amber triangles: sens i [unconstrained upper count ceilings], sens j [pooled sample fixed ceiling]). Only VCAN is classified as ROBUST across specifications, maintaining consistent negative directionality (-0.65 in primary to -1.41 in sens i). The remaining 8 genes are classified as QC-SENSITIVE, exhibiting 58% to 91% attenuation of point estimates upon ceiling removal (FN1 attenuates from -2.96 to -0.27 [90.8%], LOX from -1.72 to -0.17 [90.3%], EDIL3 from -2.61 to -0.73 [71.9%], COMP from -4.04 to -1.14 [71.7%], COL8A1 from -2.66 to -0.82 [69.0%], THBS3 from -1.57 to -0.57 [63.9%], COL3A1 from -1.96 to -0.82 [57.9%]), while INHBA reverses sign (flips from +1.08 to -0.31). This confirms that primary pseudobulk estimates are sensitive to technical count thresholds and donor composition in small single-cell cohorts.

---

## Figure 5A: Evidence-Tiered miRNA-Target Interactions and Cross-Species Peritoneal Intersections

**Image File:** [`results/figures/Fig_5A_mirna_hub_intersection.png`](Fig_5A_mirna_hub_intersection.png)  
**Source Data CSV:** [`results/figures/source_data/fig5a_mirna_hub_intersection_source.csv`](source_data/fig5a_mirna_hub_intersection_source.csv)  

**Caption:** Multi-tiered evaluation of microRNA convergence across multiMiR experimentally validated databases and transcriptomic profiles from peritoneal models. Left panel: Distribution of validated targeting miRNAs across 11 consensus hub genes categorized by evidentiary tier: Tier A (dark blue; functional evidence from luciferase reporter assays, western blot, or qPCR in miRTarBase) vs. Tier B (light blue; high-throughput CLIP-Seq or curated predictions in TarBase/miRWalk). COL3A1 (10 Tier A, 266 Tier B), FN1 (7 Tier A, 492 Tier B), and LOX (6 Tier A, 374 Tier B) demonstrate dense functional regulatory coverage. Right panel: Hypergeometric overlap testing across external peritoneal profiles at a 1.5-fold differential expression threshold (|log2FC| >= 0.58). Rodent peritoneal dialysis effluent (GSE130387; C57BL/6 mouse effluent assayed on Affymetrix miRNA-4.0 rat probes mapped by stem homology to human) displays significant cross-species convergence for both Tier A (k=7 observed vs. 1.99 expected, P = 0.00133) and Tier B (k=46 observed vs. 34.28 expected, P = 0.00154). In human dialysis effluent exosomes (GSE182736; n=3 vs. 3), Tier A overlap is k=6 (expected 8.90, P = 0.9505) and Tier B overlap is k=186 (expected 193.88, P = 0.9849), which are non-significant and appropriately interpreted as descriptive only.

---

## Figure 5C: In Vitro Mesothelial Cell Response to TGF-beta1 and Gating Against Untreated Culture Drift (GSE121372)

**Image File:** [`results/figures/Fig_5C_gse121372_tgfb1_timecourse.png`](Fig_5C_gse121372_tgfb1_timecourse.png)  
**Source Data CSV:** [`results/figures/source_data/fig5c_gse121372_tgfb1_timecourse_source.csv`](source_data/fig5c_gse121372_tgfb1_timecourse_source.csv)  

**Caption:** Bar chart illustrating log2 fold changes in human peritoneal mesothelial cells (HPMCs) exposed to 1 ng/mL TGF-beta1 relative to time-matched untreated controls at 6 hours (light blue) and 24 hours (dark blue), side-by-side with untreated culture drift (hatched amber; untreated control at 24h vs. 6h) on Illumina HumanRef-8 v2.0 beadchips. DESCRIPTIVE ONLY; unreplicated n=1 per condition; no hypothesis testing or P-values are evaluated. Ten of 11 hub genes are assayed on this platform (ISM1 is not represented). Applying the audited rule order demonstrates that baseline culture drift (|drift| >= 1.0) confounds apparent 24h induction for VCAN (drift = +1.90), COMP (drift = -2.32), INHBA (drift = -1.47), and COL11A1 probe ILMN_1644 (drift = +1.89). Only EDIL3 (log2FC_24h = +1.61, drift = +0.93) and THBS3 (log2FC_24h = +1.08, drift = -0.04) meet criteria for being higher in the single TGF-beta1 sample without culture drift confounding. FN1 (+0.87), LOX (+0.86), COL8A1 (+0.49 / +0.79), COL3A1 (+0.44), and COL11A1 probe ILMN_3994 (+0.97) have |log2FC_24h| < 1.0 and are classified as unchanged in the single TGF-beta1 sample.

---

## Figure 6: Machine Learning Consensus Feature Selection and Cross-Validated Performance Calibration

**Image File:** [`results/figures/Fig_6_discrepancy_corrected_panels.png`](Fig_6_discrepancy_corrected_panels.png)  
**Source Data CSV:** [`results/figures/source_data/fig6_discrepancy_corrected_panels_source.csv`](source_data/fig6_discrepancy_corrected_panels_source.csv)  

**Caption:** Methodological calibration of machine learning consensus selection and external validation diagnostics. Left panel: Consensus vote tallies across four distinct feature selection algorithms (LASSO, SVM-RFE, Random Forest, XGBoost) applied to 40 convergent candidates in the discovery cohort (GSE62928). Maximum observed agreement across algorithms is 3/4 votes (ISM1, FN1, EDIL3) and 2/4 votes (VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX); zero genes were selected unanimously (0 of 11 with 4/4 votes). Right panel: Receiver operating characteristic (ROC) curves in the external peritoneal dialysis effluent validation cohort (GSE125498, N=33). While the apparent in-sample composite model achieves an apparent AUC of 0.8769 (blue solid line), rigorous 50-repeat 5-fold cross-validation with within-fold scaling yields a modest mean generalization AUC of 0.6585 ± 0.0705 (orange dashed line), a single seed-42 split yields AUC of 0.5423 (green dotted line), and 1,000 label permutations yield an empirical P = 0.0879 (not statistically significant). The shaded region illustrates the optimism gap (ΔAUC = 0.2184), demonstrating that high apparent in-sample classifier performance reflects overparameterization on small sample sizes rather than robust out-of-sample diagnostic utility.

---

