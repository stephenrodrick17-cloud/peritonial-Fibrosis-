# Discrepancy Report: Figure 4 Diagnostic Biomarker Analysis Audit

**Date of Audit:** October 15, 2026  
**Project Path:** `d:/Peritoneal Project`  
**Focus:** Resolution of discrepancies between legacy written results text (Version A) and current verified analysis/figure (Version B).

---

## 1. Executive Summary of Discrepancies

The peritoneal dialysis diagnostic biomarker workflow underwent a critical transition from an earlier legacy analysis (Version A, based on 38/39 features and simulated figure curves) to a fully reconciled, empirical analysis (Version B, based on 43 module-assigned matrisome DEGs and 100% real patient data).

| Conflict ID | Feature / Metric | Legacy Written Text (Version A) | Current Figure & Code (Version B) | Root Cause | Resolution |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **DISC-01** | **Input Feature Pool** | 38 or 39 candidate genes (`ml_39_gene_results.csv`) | 43 module-assigned matrisome DEGs (`ml_43_matrisome_results.csv`) | Version A excluded grey-module genes (39 non-grey genes). Version B retained all 43 DEGs assigned to WGCNA modules (including 4 grey-module genes: CXCL14, PI3, SPP1, TGM3). | Adopt the full 43 module-assigned matrisome DEG pool as the unified candidate set. |
| **DISC-02** | **LASSO Selected Count** | 11 genes selected | 17 genes selected | Candidate pool expansion from 39 to 43 genes caused 5-fold CV to select $C = 5.0$ ($\lambda = 0.20$) instead of $C = 0.50$ ($\lambda = 2.00$). At $C = 5.0$, 17 genes have non-zero coefficients. | Report 17 selected genes at optimal regularisation parameter ($\lambda_{\min} = 0.20$), with 10 non-zeros at the continuous $\lambda_{\min} = 3.51$. |
| **DISC-03** | **Three-Way Venn Overlap** | 5 consensus genes (TNFSF15, ADAM19, FLT3LG, NAT6, ECM1) | 0 genes in 3-way intersection (LASSO $\cap$ SVM-RFE $\cap$ RF = 0) | Version A hand-drew a simulated Venn diagram placing 5 arbitrary genes in the center. In the actual 43-gene run, SVM-RFE ($n=10$) and RF ($n=10$) have zero intersection ($S \cap R = 0$). | Report the exact mathematical set overlap: 3-way = 0; consensus is defined by majority vote ($\ge 2/3$ algorithms = 10 hub genes). |
| **DISC-04** | **Consensus Diagnostic Panel** | 5 core genes: TNFSF15, ADAM19, FLT3LG, NAT6, ECM1 | 10 hub genes: TNFSF15, LTB, EBI3, CCBE1, NAT6, BMP7, SFRP4, TGM3, COL8A1, CLEC2D | In the 43-gene run, FLT3LG and ADAM19 were selected only by RF (1 vote); ECM1 was selected by 0 models. Conversely, 8 genes were selected by LASSO+SVM and 2 by LASSO+RF. | Update clinical biomarker text to the 10 verified majority-vote consensus genes ($8 + 2 = 10$). |
| **DISC-05** | **Discovery Individual AUCs** | 0.858–0.946 (reported as "strong discriminatory capacity") | 0.535–0.904 (apparent AUCs across all 10 hub genes) | Version A reported apparent AUCs only for the hand-picked 5 genes. Version B computes empirical AUCs across all 10 hub genes, including moderate performers (COL8A1: 0.581, CLEC2D: 0.535). | Tabulate apparent, LOOCV, and repeated 5-fold CV AUCs with 95% bootstrap CIs for all 10 genes. |
| **DISC-06** | **Discovery Multigene Composite Model** | 5-gene ridge model apparent AUC = 0.9808 (95% CI: 0.942–1.000) | 10-hub ridge model apparent AUC = 1.000 (LOOCV AUC = 0.9692) | Version A trained a 5-gene ridge model. Version B evaluates single diagnostic genes in Figure 4 and audits apparent optimism. | Explicitly report apparent optimism (apparent AUC = 1.000 vs LOOCV AUC = 0.969) in the text. |
| **DISC-07** | **Verification Cohort Nature** | Described as "effluent bulk RNA-seq cohort with LPD/SPD labels" | Primary mesothelial cell single-cell RNA-seq (scRNA-seq) pseudobulk (GSE248762) | Metadata mislabeling in earlier text. GSE248762 consists of 16 participants: Long-Vintage (LV, $n=10$) vs Short-Vintage (SV, $n=6$). | Accurately describe GSE248762 as single-cell pseudobulk comparing Long-Vintage ($\ge 3$–$5$ years) against Short-Vintage ($< 1$ year). |
| **DISC-08** | **NAT6 Verification Status & Curve** | Stated as validated with smooth ROC curve (AUC = 0.700 or 0.500) | Missing from mesothelial pseudobulk under symbol "NAT6" (renamed NAA80/NAA50); smooth curve was an artificial parametric plot | `generate_fig4_v2.py` used a fallback power formula $tpr = fpr^{1/AUC}$ when NAT6 was missing from the count table. In all-cell pseudobulk, NAA80 AUC = 0.500. | Replaced artificial smooth curve with empirical step-function ROC curves; clearly state that NAT6 was absent in mesothelial pseudobulk. |
| **DISC-09** | **CCBE1 Verification Inversion** | Reported as AUC = 0.983 (claimed as strong validation) | Inverted prediction! True discovery-consistent AUC = 0.017 (mesothelial) / 0.067 (all-cell) | In discovery, CCBE1 is upregulated in LPD (log2FC = +0.605, $P = 0.006$). In verification, CCBE1 is markedly elevated in Short-Vintage (log2FC = -3.129, $P = 0.00125$). Inverting the label gave 0.983. | Report the honest, direction-consistent AUC (0.017 / 0.067), documenting the complete directional reversal as a failure of external replication. |
| **DISC-10** | **Statistical Footnotes & Annotations** | Code used Mann-Whitney U test, but footnotes had conflicting star thresholds and overlapping annotations | Mann-Whitney U two-sided test with exact thresholds (`*p<0.05`, `**p<0.01`, `***p<0.001`, `****p<0.0001`, `ns`) | Legend and axis overlapping in Panels A, D, and H. | Cleaned up all figure annotations, eliminated label collisions, and harmonized footnotes with code. |

---

## 2. Evidence Traceability Table

| Identifier | Variable / Value | Source File | Source Script & Line |
| :--- | :--- | :--- | :--- |
| **DISC-01** | 43 candidate genes | `results/revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv` | `scripts/run_ml_43_genes.py:32` |
| **DISC-02** | LASSO $n=17$, $C=5.0$, $\lambda_{\min}=0.20$ | `results/revision/ml_43_matrisome_wgcna/ml_43_matrisome_results.csv` | `scripts/run_ml_43_genes.py:45-56` |
| **DISC-03** | Venn: L=17, S=10, R=10, L∩S=8, L∩R=2, S∩R=0, 3-way=0 | `results/revision/ml_43_matrisome_wgcna/ml_43_matrisome_results.csv` | `scripts/generate_fig4_perfect.py:65-74` |
| **DISC-04** | 10 consensus hub genes | `results/revision/ml_43_matrisome_wgcna/ml_43_matrisome_summary.json` | `scripts/run_ml_43_genes.py:120-138` |
| **DISC-05** | Discovery AUCs (0.535–0.904) | `results/tables/Table_10_hub_genes_discovery_validation_stats.csv` | `scripts/generate_fig4_perfect.py:126-135` |
| **DISC-06** | 10-hub Ridge Model LOOCV AUC = 0.969 | `results/tables/Final_numbers.csv` | `scripts/audit_figure4_complete.py:LOOCV` |
| **DISC-07** | GSE125498 (13 LPD, 20 SPD); GSE248762 (10 LV, 6 SV) | `results/tables/GSE125498_sample_metadata.csv` & `pb_metadata.csv` | `scripts/generate_fig4_perfect.py:77-90` |
| **DISC-08** | NAT6 / NAA80 absent in mesothelial pseudobulk | `data/processed/pseudobulk/pb_counts.csv` | `scripts/generate_fig4_perfect.py:118` |
| **DISC-09** | CCBE1 direction-consistent AUC = 0.017 (meso) / 0.067 (all-cell) | `results/tables/Table_10_hub_genes_discovery_validation_stats.csv` | `scripts/generate_fig4_perfect.py:270` |
| **DISC-10** | Mann-Whitney U test p-values & stars | `results/tables/Table_10_hub_genes_discovery_validation_stats.csv` | `scripts/generate_fig4_perfect.py:320-345` |
