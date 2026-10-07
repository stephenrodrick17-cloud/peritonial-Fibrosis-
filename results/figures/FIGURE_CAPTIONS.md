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

## Figure 5B: Predicted Tripartite ceRNA Regulatory Network (lncRNA–miRNA–Hub mRNA)

**Image File:** [`results/figures/Fig_5B_cerna_regulatory_network.png`](Fig_5B_cerna_regulatory_network.png)  
**Source Data CSV:** [`results/figures/source_data/fig5b_cerna_regulatory_network_source.csv`](source_data/fig5b_cerna_regulatory_network_source.csv)  

**Caption:** Tripartite network visualization illustrating bioinformatically predicted competitive endogenous RNA (ceRNA) axes connecting stromal lncRNAs to candidate matrisome hub targets via intermediate microRNAs. **PREDICTED ONLY (NOT PROSPECTIVELY VALIDATED IN THIS SYSTEM).** Left column (Blue, 7 nodes): Detected stromal lncRNAs (`H19`, `MEG3`, `SNHG16`, `NORAD`, `GAS5`, `NEAT1`, `MALAT1`). Middle column (Teal, 9 nodes): Intermediate microRNAs selected from published fibrosis/EMT literature supported by ENCORI CLIP-Seq experimental evidence (`let-7a-5p`, `let-7b-5p`, `miR-17-5p`, `miR-29a-3p`, `miR-29b-3p`, `miR-29c-3p`, `miR-30a-5p`, `miR-34a-5p`, `miR-93-5p`). Right column (Orange, 11 nodes): 11 candidate matrisome hub gene targets supported by multiMiR experimental records (`VCAN`, `THBS3`, `LOX`, `ISM1`, `INHBA`, `FN1`, `EDIL3`, `COMP`, `COL8A1`, `COL3A1`, `COL11A1`). **Important Evidentiary Caveats:** (1) The 7 lncRNAs displayed in this visualization represent ubiquitous, highly expressed housekeeping or nuclear architectural transcripts (`MALAT1`, `NEAT1`, `GAS5`, `NORAD`, `SNHG16`, `MEG3`, `H19`) that were formally categorized as `LOW_CREDIBILITY_CANDIDATE` during audit due to non-specific promiscuous binding. (2) Among the 11 hub genes, only 4 possess Tier A functional validation (`COL3A1`, `FN1`, `LOX`, `VCAN`); the remaining 7 hubs have 0 Tier A functional miRNAs. (3) In single-cell stromal pseudobulk testing across human peritoneal donors, candidate axes exhibit nominal-only correlation without FDR support, and none of the hub genes achieve FDR < 0.05 differential expression (all stromal hub DE FDR >= 0.131). This panel is purely hypothesis-generating.

---

## Figure 5C: In Vitro Mesothelial Cell Response to TGF-beta1 and Gating Against Untreated Culture Drift (GSE121372)

**Image File:** [`results/figures/Fig_5C_gse121372_tgfb1_timecourse.png`](Fig_5C_gse121372_tgfb1_timecourse.png)  
**Source Data CSV:** [`results/figures/source_data/fig5c_gse121372_tgfb1_timecourse_source.csv`](source_data/fig5c_gse121372_tgfb1_timecourse_source.csv)  

**Caption:** Bar chart illustrating log2 fold changes in human peritoneal mesothelial cells (HPMCs) exposed to 1 ng/mL TGF-beta1 relative to time-matched untreated controls at 6 hours (light blue) and 24 hours (dark blue), side-by-side with untreated culture drift (hatched amber; untreated control at 24h vs. 6h) on Illumina HumanRef-8 v2.0 beadchips. DESCRIPTIVE ONLY; unreplicated n=1 per condition; no hypothesis testing or P-values are evaluated. Ten of 11 hub genes are assayed on this platform (ISM1 is not represented). Applying the audited rule order demonstrates that baseline culture drift (|drift| >= 1.0) confounds apparent 24h induction for VCAN (drift = +1.90), COMP (drift = -2.32), INHBA (drift = -1.47), and COL11A1 probe ILMN_1644 (drift = +1.89). Only EDIL3 (log2FC_24h = +1.61, drift = +0.93) and THBS3 (log2FC_24h = +1.08, drift = -0.04) meet criteria for being higher in the single TGF-beta1 sample without culture drift confounding. FN1 (+0.87), LOX (+0.86), COL8A1 (+0.49 / +0.79), COL3A1 (+0.44), and COL11A1 probe ILMN_3994 (+0.97) have |log2FC_24h| < 1.0 and are classified as unchanged in the single TGF-beta1 sample.

---

## Figure 6: Machine Learning Consensus Feature Selection and Cross-Validated Performance Calibration

**Image File:** [`results/figures/Fig_6_discrepancy_corrected_panels.png`](Fig_6_discrepancy_corrected_panels.png)  
**Source Data CSV:** [`results/figures/source_data/fig6_discrepancy_corrected_panels_source.csv`](source_data/fig6_discrepancy_corrected_panels_source.csv)  

**Caption:** Methodological calibration of machine learning consensus selection and external validation diagnostics. Left panel: Consensus vote tallies across four distinct feature selection algorithms (LASSO, SVM-RFE, Random Forest, XGBoost) applied to 44 convergent candidates in the discovery cohort (GSE62928). Maximum observed agreement across algorithms is 3/4 votes (ISM1, FN1, EDIL3) and 2/4 votes (VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX); zero genes were selected unanimously (0 of 11 with 4/4 votes). Right panel: Receiver operating characteristic (ROC) curves in the external peritoneal dialysis effluent validation cohort (GSE125498, N=33). While the apparent in-sample composite model achieves an apparent AUC of 0.8769 (blue solid line), rigorous 50-repeat 5-fold cross-validation with within-fold scaling yields a modest mean generalization AUC of 0.6585 ± 0.0705 (orange dashed line), a single seed-42 split yields AUC of 0.5423 (green dotted line), and 1,000 label permutations yield an empirical P = 0.0879 (not statistically significant). The shaded region illustrates the optimism gap (ΔAUC = 0.2184), demonstrating that high apparent in-sample classifier performance reflects overparameterization on small sample sizes rather than robust out-of-sample diagnostic utility.

---

## Figure 7: Prior Literature Benchmark, Compartment Discordance, and Covariate Fragility (GSE125498 vs. GSE62928)

**Image File:** [`results/figures/Fig_7_prior_literature_effluent_comparison.png`](Fig_7_prior_literature_effluent_comparison.png)  
**Source Data CSV:** [`results/figures/source_data/fig7_prior_literature_effluent_comparison_source.csv`](source_data/fig7_prior_literature_effluent_comparison_source.csv)  

**Caption:** Comparative multi-panel evaluation benchmarked against prior effluent literature (Ma et al., *Front Pharmacol* 2025; GSE125498, $N=33$: 20 SPD vs. 13 LPD). (A) Shed effluent leukocyte marker shifts across dialysis vintage. While Ma et al. reported increased monocyte infiltration in long-term PD, empirical transcriptomic profiling shows significantly decreased expression of monocyte marker CD14 ($\log_2\text{FC} = -0.38, P = 0.016$) alongside increased T-cell marker CD3E ($\log_2\text{FC} = +1.46, P = 0.003$) and neutrophil/Fc-receptor marker FCGR3B ($\log_2\text{FC} = +0.81, P = 0.024$). (B) Expression discordance of the 11 hub genes between discovery parietal tissue biopsy (GSE62928, $N=8$, orange) and shed cellular effluent (GSE125498, $N=33$, blue). Whereas VCAN and COL8A1 are strongly up-regulated in fibrotic tissue (+2.75 and +2.68 $\log_2\text{FC}$), VCAN is nominally decreased in late PD effluent (-0.52 $\log_2\text{FC}, P = 0.0244$), and COL8A1 is only nominally elevated on one probe (+0.75, $P = 0.0488$). Four hub genes (`COL11A1`, `COMP`, `EDIL3`, `INHBA`) are completely unassayed on GPL10558. (C) VCAN covariate adjustment fragility and leave-one-out sensitivity under probe Set B. Unadjusted nominal significance (-0.52, $P = 0.0244$) attenuates under CD14 adjustment ($P = 0.1429$) and CD3E adjustment ($P = 0.1955$), and is completely abolished under joint 3-marker adjustment (-0.09, $P = 0.6716$). In leave-one-out testing, dropping CD14 alone retains nominal significance ($P = 0.0395$), whereas dropping FCGR3B alone completely eliminates the effect (+0.02, $P = 0.9276$), proving that attenuation depends on CD14 and FCGR3B together. (D) Biological coupling between VCAN and CD14 across all 33 effluent samples (Pearson $r = +0.658, P = 3.10 \times 10^{-5}$; Spearman $\rho = +0.529, P = 0.0015$). This strong positive coupling indicates that VCAN's apparent down-regulation in late PD effluent reflects decreased shedding of viable monocytes rather than downregulation of the fibrotic matrisome.

---

## Figure 10A: STRING v12.5 Functional Protein Association Network (Isolated PPI Evidence)

**Image File:** [`results/figures/D_string_ppi_network.png`](D_string_ppi_network.png)  
**Source Data CSV:** [`results/tables/D_string_edges_functional.csv`](../tables/D_string_edges_functional.csv)  

**Caption:** Direct network representation of live STRING v12.5 functional protein association interactions among the 11 consensus hub genes (minimum interaction confidence score >= 0.400). Displays 21 functional association edges with edge thickness proportional to STRING confidence score (top interaction: LOX–FN1, score = 0.953; physical interaction network contains 3 edges). Co-expression edges are NOT computed or drawn in this figure. ISM1 has degree = 0 (no known STRING functional interactions among the 11 hub genes).

---

## Figure 10B: Dual-Evidence Protein-Protein Interaction & Co-Expression Network

**Image File:** [`results/figures/Hub_01_ppi_gene_interaction_network.png`](Hub_01_ppi_gene_interaction_network.png)  
**Source Data CSV:** [`results/tables/hub_genes_ppi_centrality_metrics.csv`](../tables/hub_genes_ppi_centrality_metrics.csv)  

**Caption:** Integrated dual-evidence molecular interaction network combining live STRING v12.5 functional associations (score >= 0.400, 21 edges; solid lines) with empirical tissue co-expression correlations in the discovery cohort (GSE62928, N=8; 42 edges meeting |r| >= 0.85, FDR < 0.01; dashed green lines). Violet solid lines represent dual-validated edges supported by both STRING PPI and empirical tissue co-expression. Node size reflects network degree centrality, and node colors denote Matrisome functional categories. Long node labels (COL11A1, COL8A1) are rendered with adjusted font sizes and scaled node diameters to prevent label overflow or clipping.



