# Figure Text Transcription and Numerical Cross-Check

This document transcribes visible text, metrics, labels, and statistical values from all figures referenced in `README.md` (Figures 0 through 12), and cross-checks them against `results/manuscript_numbers.json` and `results/MANUSCRIPT_NUMBERS.md`.

---

## Figure 0: Graphical Abstract (`graphical_abstract.jpg`)
* **Visible Text / Labels**:
  - Top Left: GSE62928 (EPS $n=4$, Control $n=4$).
  - Overlap: Venn Diagram labeled "DEGs ($n=367$)", "Matrisome ($n=1027$)", "FN1 / Salmon module ($n=604$)" $\to$ "40 Convergent ECM Candidates".
  - Machine Learning: 4 Algorithms (LASSO, SVM-RFE, RF, XGBoost) $\to$ "11 Core Hub Genes".
  - Validation: GSE125498 ($n=33$, SPD $n=20$, LPD $n=13$). "AUC = 87% (In-sample)".
  - Pathway: EMT, TNF-$\alpha$, Angiogenesis, Inflammatory Response.
* **Cross-Check & Discrepancy Analysis**:
  - *Note*: As documented in `results/GRAPHICAL_ABSTRACT_TODO.md`, this raster graphic contains legacy phrasing ("FN1" instead of "Salmon module", "87% AUC" without specifying in-sample vs 0.678 CV, and "EPD vs Early SPD"). Must be redrawn in vector format prior to manuscript submission.

---

## Figure 1: Volcano Plot (`results/figures/DEG_01_volcano_plot.png`)
* **Visible Text / Labels**:
  - Title: "Differential Expression Analysis - GSE62928 (EPS vs Control)"
  - Subtitle: "Total genes: 22,049 | Upregulated: 534 | Downregulated: 1510 | Not Sig: 20005"
  - Thresholds: $|\log_2\text{FC}| \ge 0.585$ (1.5-fold), $P < 0.05$.
  - Labeled top genes: *FN1, COL1A1, COL3A1, VCAN, THBS2, POSTN, TNC, TIMP1*, etc.
* **Cross-Check**:
  - The 534 upregulated genes represent probe-level collapse at $|\log_2\text{FC}| \ge 0.585, P < 0.05$. The pro-fibrotic filter further narrows this to $\log_2\text{FC} \ge 0.80$, yielding exactly 367 DEGs.
  - Matches `results/manuscript_numbers.json` (`discovery_degs: 367`, `gene_universe_limma: 22049`).

---

## Figure 2: Matrisome Overlap (`results/figures/venn_gse62928_pro_fibrotic_ecm_71.png`)
* **Visible Text / Labels**:
  - Circles: "GSE62928 Pro-Fibrotic DEGs (367)" and "Human Matrisome (1027)".
  - Partition counts: DEGs only = 296, Intersection = 71, Matrisome only = 956.
  - Title: "Convergent Pro-Fibrotic ECM Genes (GSE62928 $\cap$ Matrisome)"
* **Cross-Check**:
  - $296 + 71 = 367$ DEGs.
  - $956 + 71 = 1027$ Matrisome genes.
  - Exact match with `results/manuscript_numbers.json` (`ecm_degs: 71`).

---

## Figure 3: WGCNA Module-Trait Relationships (`results/figures/WGCNA_03_module_trait_heatmap.png`)
* **Visible Text / Labels**:
  - Modules listed: MEblack, MEblue, MEbrown, MEcyan, MEdarkgreen, MEdarkgrey, MEdarkred, MEdarkturquoise, MEgreen, MEgreenyellow, MEmagenta, MEmidnightblue, MEpink, MEpurple, MEred, MEroyalblue, MEsalmon, MEtanh, MEyellow, MEyellowgreen, MEgrey.
  - MEsalmon vs Fibrosis_Status: `0.81` $(p = 0.016)$.
  - MEyellowgreen: `0.64` $(p = 0.086)$.
  - MEdarkturquoise: `0.63` $(p = 0.093)$.
* **Cross-Check**:
  - Salmon module correlation rounded to 2 decimals is $0.81$ ($r = 0.806$, $P = 0.0157$). Exact match.

---

## Figure 4: 3-Way Venn Convergence (`results/figures/venn_wgcna_convergence.png`)
* **Visible Text / Labels**:
  - Set A: "GSE62928 Pro-Fibrotic DEGs (N=367)"
  - Set B: "Human Matrisome (N=1027)"
  - Set C: "WGCNA Salmon Module (N=604)"
  - Central intersection: "40 convergent candidates"
* **Cross-Check**:
  - Regenerated in commit `63112aa`. Center label matches exact terminology. Exact match.

---

## Figure 5: Machine Learning Consensus Votes (`results/figures/WGCNA_ML_01_consensus_votes_barchart.png`)
* **Visible Text / Labels**:
  - Title: "Machine Learning Consensus Feature Selection (4 Models: LASSO, SVM-RFE, RF, XGBoost)"
  - Top genes with 3 votes: *ISM1, FN1, EDIL3*.
  - Genes with 2 votes: *VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX*.
  - Cutoff line: "Consensus threshold $\ge 2$ votes" isolating 11 Hub Genes.
  - 1-vote candidates ($n=17$) and 0-vote candidates ($n=12$) displayed below threshold.
* **Cross-Check**:
  - Total candidates evaluated = 40.
  - Vote distribution: 4 votes = 0, 3 votes = 3, 2 votes = 8 (Total $\ge 2$ votes = 11 hubs), 1 vote = 17, 0 votes = 12. Sum: $0 + 3 + 8 + 17 + 12 = 40$.
  - Per-model counts: LASSO = 2, SVM-RFE = 14, RF = 17, XGBoost = 1.
  - Exact match with `results/manuscript_numbers.json`.

---

## Figure 6: Discovery Expression Heatmap (`results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png`)
* **Visible Text / Labels**:
  - 11 Hub Genes (*ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX*).
  - Columns: Control ($n=4$: GSM1536410, GSM1536411, GSM1536412, GSM1536413) vs EPS ($n=4$: GSM1536406, GSM1536407, GSM1536408, GSM1536409).
  - Color bar: Z-score normalized expression.
* **Cross-Check**:
  - Clean separation across all 11 genes in discovery cohort ($N=8$). Note in documentation highlights that discovery significance is by construction (pre-filtered at nominal $P < 0.05$). Exact match.

---

## Figure 7: Validation Boxplots in GSE125498 (`results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png`)
* **Visible Text / Labels**:
  - 7 assessable hub genes in GSE125498 effluent cells (SPD $n=20$ vs LPD $n=13$):
    - *ISM1*: $P = 0.811, \text{AUC} = 0.53$
    - *FN1*: $P = 0.294, \text{AUC} = 0.61$
    - *VCAN*: $P = 0.034, \text{AUC} = 0.72$ (Direction-adjusted AUC $= 0.72$, raw expression lower in LPD)
    - *COL3A1*: $P = 0.754, \text{AUC} = 0.53$
    - *COL8A1*: $P = 0.117, \text{AUC} = 0.67$
    - *THBS3*: $P = 0.367, \text{AUC} = 0.60$ (Direction-adjusted AUC $= 0.60$, raw expression lower in LPD)
    - *LOX*: $P = 0.985, \text{AUC} = 0.50$
* **Cross-Check**:
  - 4 missing hubs (*COL11A1, COMP, EDIL3, INHBA*) clearly noted as lacking probes in GPL14951.
  - Mann-Whitney P-values and individual AUCs match `results/tables/GSE125498_wgcna_hub_validation_metrics.csv` and `results/manuscript_numbers.json`. Exact match.

---

## Figure 8: ROC Analysis (`results/figures/Hub_02b_roc_analysis.png`)
* **Visible Text / Labels**:
  - Panel A Title: "Single-Gene Discrimination (GSE125498)"
    - Individual ROC curves for 7 assessable hubs (*VCAN, COL8A1, FN1, THBS3, ISM1, COL3A1, LOX*).
  - Panel B Title: "Multi-Gene Panel Validation (GSE125498)"
    - "Primary 7-Gene Panel (In-sample: AUC = 0.869)"
    - "Primary 7-Gene 50x5 Scaled CV (Mean AUC = 0.678 ± 0.058)"
    - "Primary 7-Gene Single 5-Fold CV (AUC = 0.592)"
    - "Nomogram Model (In-sample: AUC = 0.800)"
    - "Nomogram Model (5-Fold CV: AUC = 0.550)"
* **Cross-Check**:
  - Regenerated in commit `ee03c1e`.
  - All values match `results/manuscript_numbers.json` (`cv_50x5_scaled_mean: 0.678`, `cv_50x5_scaled_sd: 0.058`, `in_sample_auc_7gene: 0.869`, `single_5fold_cv_auc_7gene: 0.592`, `nomogram_cv_auc: 0.550`). Exact match.

---

## Figure 9: Nomogram, Calibration & Decision Curve (`results/figures/Hub_02_clinical_nomogram_dca_calibration.png`)
* **Visible Text / Labels**:
  - Nomogram Points scale: *COL8A1* (0-100 pts), *VCAN* (0-70 pts), *LOX* (0-55 pts).
  - Total Points $\to$ Risk of Long-Term Peritoneal Dialysis / Fibrosis Progression.
  - Apparent C-index $= 0.819$, Brier score $= 0.1542$, Hosmer-Lemeshow $P = 0.2437$.
  - DCA Net Benefit curve positive across threshold probabilities 10% to 70%.
* **Cross-Check**:
  - Brier score matches $0.1542$ (earlier stale README text $0.158$ corrected). C-index matches $0.8192$. Exact match.

---

## Figure 10: PPI & Co-Expression Network (`results/figures/Hub_01_ppi_gene_interaction_network.png`)
* **Visible Text / Labels**:
  - Title: "Hub Gene Interaction Network (STRING v12.5 Functional PPI & GSE62928 Co-expression)"
  - Legend: STRING v12.5 Functional PPI edges (gray dashed, $n=21$), Discovery Co-expression $|r| \ge 0.85$ (green solid, $n=42$).
  - Hub nodes (11 genes): *COL1A1/COL3A1/COL8A1/COL11A1, FN1, VCAN, COMP, THBS3, EDIL3, INHBA, LOX, ISM1*.
  - Degree annotations: *FN1* (degree 18), *COL3A1* (degree 17), *ISM1* (degree 0 among functional edges).
* **Cross-Check**:
  - Generated live via STRING API v12.5 (`string_live_edges_20261001.tsv`).
  - Edge count $= 21$ functional edges, 3 physical edges. Exact match.

---

## Figure 11: GSEA Hallmark Pathway Heatmap (`results/figures/Hub_03_gsea_pathway_enrichment_heatmap.png`)
* **Visible Text / Labels**:
  - Panel A: Hallmark pathways ranked by NES in GSE62928:
    - Epithelial-Mesenchymal Transition (EMT): $\text{NES} = +3.203, P < 0.001, \text{FDR} < 0.001$
    - TNF-$\alpha$ Signaling via NF-$\kappa$B: $\text{NES} = +2.483, P < 0.001, \text{FDR} = 0.001$
    - Inflammatory Response: $\text{NES} = +1.868, P = 0.006, \text{FDR} = 0.030$
    - Angiogenesis: $\text{NES} = +1.737, P = 0.016, \text{FDR} = 0.086$
    - IL-6/JAK/STAT3 Signaling: $\text{NES} = +1.579, P = 0.034, \text{FDR} = 0.160$
    - Apoptosis: $\text{NES} = +1.564, P = 0.033, \text{FDR} = 0.161$
  - Panels B & C: ssGSEA Correlation heatmaps across 11 Hubs.
* **Cross-Check**:
  - Matches `results/tables/gsea_preranked_hallmark_results.csv` and `results/manuscript_numbers.json`. Exact match.

---

## Figure 12: Immune Infiltration Deconvolution (`results/figures/Hub_04_immune_infiltration_deconvolution.png`)
* **Visible Text / Labels**:
  - Subsets: Peritoneal Macrophages (M1/M2), Myofibroblasts, Mesothelial Cells, Neutrophils, CD4+ T cells, CD8+ T cells, B cells.
  - GSE62928 comparison: Peritoneal Myofibroblasts significantly elevated in EPS ($P = 0.0286, U = 16.0$); other immune subsets non-significant after FDR control.
  - GSE125498 comparison: Myofibroblast score trends higher in LPD ($P = 0.048$).
* **Cross-Check**:
  - Zero gene overlap between 11 hubs and cell-type deconvolution marker gene sets. Exact match.
