# Figures with Graphical / Numerical Discrepancies to Regenerate

This document lists existing visual figures in the repository that contain historical discrepancies or layout mismatches identified during the consistency audit, along with the precise modifications required for future graphical updates.

---

### 1. `results/figures/graphical_abstract.jpg`
- **Current Content:** Stylized summary infographic of the study workflow.
- **Identified Discrepancies:**
  - **Venn Labeling:** Contains floating label `"FN1"` and duplicated `"extracellular matrix"` labels in the candidate overlap section.
  - **Cohort Cohort Grouping Typo:** Contains label `"EPD vs. Early SPD Late LPD"` instead of standard clinical terminology (e.g., Short-term vs. Long-term PD peritoneal effluent cell pellets).
  - **Consensus Vote Count:** Mentions "4-Way Consensus Machine Learning (11 Hub Genes with 4/4 Votes)"; true maximum agreement is 3/4 votes (3 genes) and 2/4 votes (8 genes).
  - **In-sample AUC Overstatement:** Highlights composite `"AUC 0.87"` without designating it as an uncorrected in-sample fit subject to overfitting (vs. realistic out-of-fold generalization $\text{AUC} = 0.678$, $\text{SD } 0.058$).
  - **Network Physical Label:** Refers to "20 Physical Interactions"; the true network is a **functional association network** (21 edges, score $\ge 0.400$), whereas physical binding interactions comprise only 3 edges.
  - **Validation Claim:** Uses absolute wording `"Validated ... remodeling"` rather than acknowledging effluent cellular proxy discordances.
- **Required Action:** Redraw graphical abstract with corrected labels (2–3 consensus votes, 21 functional STRING edges, in-sample vs. cross-validated AUC $0.869$ vs. $0.678$, and effluent cell proxy phenotype).

---

### 2. `venn_wgcna_convergence.png`
- **Current Content:** 3-way Venn diagram showing DEGs (367), Curated Matrisome (1,027), and Salmon Module (604).
- **Identified Discrepancies:**
  - Central 3-way intersection is labeled `"40 (Core Hubs)"`.
  - These 40 genes are **convergent candidate genes** (WGCNA Salmon module $\cap$ 81 ECM-DEGs), from which machine learning consensus filtering prioritized the **11 consensus hub biomarkers**.
  - Diagram should annotate the probe-level collapse universe (22,049 genes) vs. the MaxMean matrix universe (20,940 genes).
- **Required Action:** Replace `"40 (Core Hubs)"` with `"44 Convergent Candidates"` and clarify the probe universe context.

---

### 3. `results/figures/Hub_02b_roc_analysis.png`
- **Current Content:** Multi-panel ROC analysis figure.
- **Identified Discrepancies:**
  - **Panel A (Discovery ROC):** Reports multi-gene composite $\text{AUC} = 1.000$ in GSE62928 ($N=8$). This multi-gene fit is mathematically saturated/trivial for $p=11 > n=8$ and does not reflect true predictive performance.
  - **Panel B (Validation Cross-Validation):** Legend reports `"50x 5-Fold Stratified CV (AUC = 0.696 ± 0.056)"` from an unscaled model with sorted sample ordering.
  - The primary pipeline-scaled cross-validation yields per-repeat pooled $\text{AUC} = 0.678$ ($\text{SD } 0.058$ across repeats; per-fold mean $\text{AUC} = 0.701$, $\text{SD } 0.203$), LOOCV $\text{AUC} = 0.658$, and single seed-42 5-fold CV pooled $\text{AUC} = 0.592$ (fold extremes: $0.500 - 0.667$).
- **Required Action:** When regenerating, annotate Panel A as mathematically saturated, and update Panel B legend to explicitly denote the primary pipeline-scaled 50-repeat cross-validation curve ($\text{AUC} = 0.678$, $\text{SD } 0.058$) and single-split cross-validation curve ($\text{AUC} = 0.592$).

---

### 4. `results/figures/Hub_03_gsea_pathway_enrichment_heatmap.png` & `.pdf`
- **Current Content:** Heatmap of ssGSEA hallmark pathway enrichment scores across samples.
- **Identified Discrepancies:**
  - Regenerated with fresh deterministic GSEA run (EMT $\text{NES} = +3.203, \text{FDR} < 0.001$).
  - Ensure pathway row labels strictly match official MSigDB Hallmark terms and confirm absence of legacy Hypoxia mislabeling.
- **Required Action:** Check label formatting and visual alignment against `results/tables/gsea_preranked_hallmark_results.csv`.

---

### 5. `results/figures/WGCNA_ML_01_consensus_votes_barchart.png`
- **Current Content:** Bar chart displaying consensus vote counts for candidate genes.
- **Identified Discrepancies:**
  - Bar heights and tiers correctly reflect 3 votes for `ISM1`, `FN1`, `EDIL3` and 2 votes for `VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`, but subtitle references "4/4 unanimity".
- **Required Action:** Update subtitle/legend to clarify that the maximum observed agreement is 3/4 votes (3 genes) and 2/4 votes (8 genes, selected exclusively by SVM-RFE + RF).
