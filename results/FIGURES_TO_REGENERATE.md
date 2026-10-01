# Figures with Graphical / Numerical Discrepancies to Regenerate

This document lists existing visual figures in the repository that contain historical discrepancies or layout mismatches identified during the consistency audit, along with the precise modifications required for future graphical updates.

---

### 1. `results/figures/graphical_abstract.jpg`
- **Current Content:** Text refers to "4-Way Consensus Machine Learning (11 Hub Genes with 4/4 Votes)", "20 Physical Interactions", and shows a stylized flowchart.
- **Identified Discrepancies:**
  - The 11 consensus hub genes have 3/4 votes (3 genes) and 2/4 votes (8 genes); no gene achieved 4/4 unanimous votes across all four algorithms.
  - STRING PPI network is a **functional association network** (21 edges, score $\ge 0.400$), not a pure physical binding network (which has 3 edges).
  - External validation dataset GSE125498 evaluates short-term vs long-term PD peritoneal effluent cell pellets (a clinical proxy for membrane injury across dialysis vintage), not surgical EPS tissue.
- **Required Action:** Update text labels to reflect 2–3 consensus votes, 21 functional association edges, and effluent cell proxy phenotype.

---

### 2. `results/figures/WGCNA_ML_01_consensus_votes_barchart.png`
- **Current Content:** Bar chart displaying consensus vote counts for candidate genes.
- **Identified Discrepancies:**
  - Bar heights and tiers correctly reflect 3 votes for `ISM1`, `FN1`, `EDIL3` and 2 votes for `VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`, but subtitle references "4/4 unanimity".
- **Required Action:** Update subtitle/legend to clarify that the maximum observed agreement is 3/4 votes (3 genes) and 2/4 votes (8 genes, selected exclusively by SVM-RFE + RF).

---

### 3. `results/figures/WGCNA_ML_02_model_selection_heatmap.png`
- **Current Content:** Binary tile heatmap indicating which models selected which genes.
- **Identified Discrepancies:**
  - Confirms LASSO selected 2 genes, SVM-RFE selected 14 genes, RF selected 17 genes, and XGBoost selected 1 gene (`EDIL3`).
- **Required Action:** No numerical discrepancy with tables; graphical style harmonization recommended.

---

### 4. `venn_wgcna_convergence.png`
- **Current Content:** 3-way Venn diagram showing DEGs (367), Curated Matrisome (1,027), and Salmon Module (604).
- **Identified Discrepancies:**
  - Shows 2-way intersection of DEGs $\cap$ Matrisome as 71, and 3-way intersection as 40.
  - The DEG universe represents 367 probe-collapsed genes at nominal $P < 0.05, \log_2\text{FC} \ge 0.80$.
- **Required Action:** Ensure diagram caption notes probe-level collapse universe (22,049 genes) vs MaxMean expression matrix universe (20,940 genes).

---

### 5. `results/figures/Hub_02b_roc_analysis.png`
- **Current Content:** Multi-panel ROC analysis figure. Panel B legend reports "50x 5-Fold Stratified CV (AUC = 0.696 ± 0.056)".
- **Identified Discrepancies:**
  - The $\text{AUC} = 0.696 \pm 0.056$ represents the unscaled 50x5 CV model. The primary pipeline-scaled 50x5 CV model achieves per-fold $\text{AUC}_{\text{CV}} = 0.701 \pm 0.203$ (and per-repeat pooled $\text{AUC} = 0.678 \pm 0.058$).
  - Single seed-42 5-fold CV curve plotted in Panel B has pooled fold $\text{AUC} = 0.592$ (per-fold mean $= 0.617 \pm 0.061$).
- **Required Action:** When regenerating, update Panel B legend to explicitly denote the primary pipeline-scaled 50-repeat cross-validation curve and single-split cross-validation curve.
