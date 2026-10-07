# Graphical Abstract Redraw Specification & Correction Checklist

The current raster image `results/figures/graphical_abstract.jpg` was created outside the programmatic Python/R pipeline. Before publication, it must be redrawn in vector illustration software (e.g., Adobe Illustrator, BioRender, Inkscape) with the following mandatory corrections:

---

## 1. Panel 2 (Systems Biology & WGCNA)
* **Venn Diagram Labeling Error:**
  - **Current text:** Right circle of the 44 convergent genes Venn diagram is mislabeled as `"FN1"`.
  - **Correction:** Redraw label as `"WGCNA Salmon Module (604 genes)"`.

---

## 2. Panel 3 (Machine Learning Consensus)
* **Feature Selection & Vote Counts Clarification:**
  - **Current text:** Lists model names (LASSO, SVM-RFE, Random Forest, XGBoost) and a partial hub list with `"..."`.
  - **Correction:** Include exact algorithm selection counts:
    - LASSO: 2 features (`ISM1`, `FN1`)
    - SVM-RFE: 14 features
    - Random Forest: 17 features
    - XGBoost: 1 feature (`EDIL3`)
    - Consensus rule: $\ge 2$ model votes (0 with 4/4, 3 with 3/4, 8 with 2/4).
  - List all 11 consensus hubs: `ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX`.

---

## 3. Panel 4 (External Clinical Assessment in Effluent Cells)
* **Cohort Labeling Typo:**
  - **Current text:** `"EPD vs. Early SPD Late LPD"`.
  - **Correction:** Update to `"Short-term PD (SPD, 0–24 mos, n=20) vs. Long-term PD (LPD, ≥25 mos, n=13)"` and clarify that this is a *dialysis effluent cell proxy cohort*, not surgical EPS biopsy.
* **ROC Curve Optimization & Optimism Reporting:**
  - **Current text:** `"AUC 0.87"` (presented without context as if it were validated generalization).
  - **Correction:** Explicitly label:
    - **In-sample fit:** $\text{AUC} = 0.87$ ($95\%\text{ CI: } [0.71, 0.99]$)
    - **50-Repeat 5-Fold Cross-Validated (Primary):** $\text{AUC}_{\text{CV}} = 0.68$ ($\text{SD } 0.06$)
    - **Single 5-Fold Split:** $\text{AUC} = 0.59$

---

## 4. Bottom Summary Banner / Overclaim Correction
* **Typo & Overclaim:**
  - **Current text:** `"Validated extracellular matrix extracellular matrix remodeling in the peritoneal membrane"` (contains accidental duplicate phrase and overclaims "validation").
  - **Correction:** Update banner to `"High-Dimensional Matrisome Convergence & Consensus Machine Learning Define Extracellular Matrix Drivers of Peritoneal Membrane Fibrogenesis"`.
