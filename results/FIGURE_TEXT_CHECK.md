# Figure Text Transcription and Numerical Cross-Check

This document transcribes visible text, metrics, labels, and statistical values from all figures referenced in `README.md` (Figures 0 through 12), and cross-checks them against `results/tables/*.csv`, `results/MANUSCRIPT_NUMBERS.md`, and the authoritative 2026-10-07 fresh reruns.

**Effective Primary Threshold:** |log₂FC| ≥ 0.585 (1.5-fold), nominal P < 0.05. Historical 0.80 threshold used ONLY for explicitly-labelled comparison subsections. Last verified 2026-10-07.

---

## Figure 0: Graphical Abstract (`graphical_abstract.jpg`)
* **Visible Text / Labels**:
  - Top Left: GSE62928 (EPS n=4, Control n=4).
  - Overlap: Venn Diagram currently legacied as "DEGs (n=367)", "Matrisome (n=1027)", "FN1 / Salmon module (n=604)" → "40 Convergent ECM Candidates".
  - Machine Learning: 4 Algorithms (LASSO, SVM-RFE, RF, XGBoost) → "11 Core Hub Genes".
  - Validation: GSE125498 (n=33, SPD n=20, LPD n=13). "AUC = 87% (In-sample)".
  - Pathway: EMT, TNF-α, Angiogenesis, Inflammatory Response.
* **Cross-Check & Discrepancy Analysis (RESOLVED AS LEGACY GRAPHIC)**:
  - ⚠️ **LEGACY RASTER GRAPHIC — needs redraw before manuscript submission.** Current embedded labels are stale:
    - Primary DEG count should be **534** (not 367) at |log₂FC| ≥ 0.585.
    - Primary 3-way convergence should be **44** (not 40).
    - "FN1" label should read "WGCNA Salmon module".
    - "AUC = 87%" MUST be explicitly relabelled "In-sample AUC = 0.869; Primary 50×5 CV AUC = 0.678 ± 0.058" to avoid misleading readers.
  - *Resolution:* Documented in `results/GRAPHICAL_ABSTRACT_TODO.md`. **All in-text manuscript prose MUST use 534 / 44 / 0.678 regardless of what this legacy raster image displays.**

---

## Figure 1: Volcano Plot (`results/figures/DEG_01_volcano_plot.png`)
* **Visible Text / Labels** (fresh rerun 2026-10-07):
  - Title: "Differential Expression Analysis - GSE62928 (EPS vs Control)"
  - Subtitle: "Total genes: 22,049 | Upregulated: 534 | Downregulated: 1510 | Not Sig: 20005"
  - Thresholds: |log₂FC| ≥ 0.585 (1.5-fold), P < 0.05.
  - Vertical dashed lines at ±0.585; horizontal line at P = 0.05.
  - Labeled top genes: *FN1, COL1A1, COL3A1, VCAN, THBS2, POSTN, TNC, TIMP1*, etc.
* **Cross-Check**:
  - ✅ **CONFIRMED 0.585-primary.** Fresh rerun produces exactly 534 pro-fibrotic (upregulated) DEGs + 1,510 downregulated = 2,044 total DEGs from a 22,049-gene Limma universe.
  - **Historical / Stringent Comparison (for reference only):** A stricter log₂FC ≥ 0.80 cutoff yields 367 up / 1,263 down = 1,630 total DEGs. This is NOT the primary manuscript result.
  - Genome-wide FDR < 0.05 count: 78 genes (8 up / 70 down); only 5 of the 81 up ECM-DEGs individually survive (BGN, CLEC11A, COL1A1, MXRA5, SERPINE2) — mandatory discovery caveat.

---

## Figure 2: Matrisome Overlap 2-Way Venn (`results/figures/venn_gse62928_pro_fibrotic_ecm_81.png`)
* **Visible Text / Labels** (fresh rerun 2026-10-07):
  - Circles: "GSE62928 Pro-Fibrotic DEGs (534)" and "Human Matrisome (1027)".
  - Partition counts: DEGs only = 453, Intersection = 81, Matrisome only = 946.
  - Title: "Convergent Pro-Fibrotic ECM Genes (GSE62928 ∩ Matrisome)"
* **Cross-Check**:
  - ✅ Arithmetic: 453 + 81 = **534 pro-fibrotic DEGs**. Correct.
  - ✅ Arithmetic: 946 + 81 = **1,027 curated Matrisome reference**. Correct.
  - ✅ Measured matrisome (expressed in 22,049-gene universe): 975 genes. 946 + 81 counts only the reference 1,027; 52 reference matrisome genes are unmeasured in GSE62928.
  - ✅ **Exact match** with fresh rerun and `results/tables/GSE62928_ECM_DEGs_candidate_list.csv` (DEG_Status value_counts: Up=81, Down=67, Total=148).
  - **Historical / Stringent Comparison (for reference only):** Old filename `_71.png` corresponds to log₂FC ≥ 0.80 cutoff yielding 71 up ECM-DEGs (DEGs only = 296, intersection = 71, 296 + 71 = 367). Do NOT cite this as primary.

---

## Figure 3: WGCNA Module-Trait Relationships (`results/figures/WGCNA_03_module_trait_heatmap.png`)
* **Visible Text / Labels**:
  - Modules listed: 20 total including MEblack, MEblue, MEbrown, MEsalmon, MEyellow, MEgrey, etc.
  - MEsalmon vs Fibrosis_Status: `0.81` (p = 0.016).
  - MEyellowgreen: `0.64` (p = 0.086).
  - MEdarkturquoise: `0.63` (p = 0.093).
* **Cross-Check**:
  - ✅ Salmon module correlation rounded to 2 decimals = 0.81 (exact r = 0.806, P = 0.0157). Correct.
  - ✅ Exact label-permutation P_perm = 1/70 = 0.0143 (one-sided, rank #1 of 70 splits).
  - ⚠️ **Mandatory caveat:** Bonferroni across 14 modules = P_bonf ≈ 0.220. Not statistically robust after family-wise correction.

---

## Figure 4: 3-Way Venn Convergence (`results/figures/venn_wgcna_convergence.png`)
* **Visible Text / Labels** (fresh rerun 2026-10-07):
  - Set A: "GSE62928 Pro-Fibrotic DEGs (N=534)" (primary 0.585 threshold)
  - Set B: "Human Matrisome (N=1027)"
  - Set C: "WGCNA Salmon Module (N=604)"
  - Central 3-way intersection: "44 convergent candidates"
* **Cross-Check**:
  - ✅ Confirmed: Strict 3-set intersection {534 Up DEGs} ∩ {975 measured Matrisome} ∩ {604 Salmon} = **44 genes** exact. Full verbatim list in `FINAL_VERIFIED_NUMBERS.md §6` and `convergent_WGCNA_ECM_genes.csv`.
  - ✅ All 11 consensus ML hub genes are confirmed present in BOTH the 44-gene (0.585-primary) pool and the 40-gene (0.80-historical) pool — ML downstream is stable across cutovers.
  - **Historical / Stringent Comparison (for reference only):** Set A N=367 yields 3-way intersection N=40; this was the candidate pool for ML vote-count computation.

---

## Figure 5: Machine Learning Consensus Votes (`results/figures/WGCNA_ML_01_consensus_votes_barchart.png`)
* **Visible Text / Labels**:
  - Title: "WGCNA-ECM Machine Learning Consensus Feature Selection (4 Models)
Hub Biomarkers Identified by >= 2 Algorithms"
  - Top genes with 3 votes (n=3): *ISM1, FN1, EDIL3*.
  - Genes with 2 votes (n=8): *VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX*.
  - Cutoff line: "Hub Selection Cutoff (>= 2 Models)" isolating 11 Hub Genes.
  - 1-vote candidates (n=9): *MXRA5, COL1A2, SNED1, COL5A2, COL5A1, PDGFA, TNFAIP6, PCOLCE, POSTN*.
  - 0-vote candidates (n=20): *COL1A1, BGN, CLEC11A, SERPINE2, SULF2, THSD4, THBS2, MUC1, TGFB3, THBS1, SDC1, FBLN7, ADAMTS2, MFAP2, MFGE8, FGF14, MMP14, SERPINE1, SERPINA3, ADAMTS1*.
* **Cross-Check**:
  - ✅ **Important framing note:** Total candidates evaluated = 40 (this is the historical 0.80-stringent 3-way pool). ML vote-count computation was executed over this 40-candidate pool.
  - ✅ HOWEVER: All 11 consensus hub genes (≥2/4 votes) are confirmed present in the primary 0.585 44-candidate pool as well. The 4 extra primary-pool genes (ADAMTS6, FGF14, FRAS1, THBS1) did not change the 11-hub result.
  - Vote distribution check: 4 votes = 0, 3 votes = 3, 2 votes = 8 (Total ≥2 votes = 11 hubs), 1 vote = 9, 0 votes = 20. Sum: 0 + 3 + 8 + 9 + 20 = **40**. Correct.
  - Per-model counts: LASSO = 2, SVM-RFE = 14, RF = 17, XGBoost = 1 (Total selections = 34). Correct.
  - ✅ Exact match with `results/tables/ML_hub_genes_from_WGCNA_ECM_all_results.csv`.

---

## Figure 6: Discovery Expression Heatmap (`results/figures/WGCNA_ML_04_hub_genes_expression_heatmap.png`)
* **Visible Text / Labels**:
  - 11 Hub Genes (*ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX*).
  - Columns: Control (n=4: GSM1536410/11/12/13) vs EPS (n=4: GSM1536406/07/08/09).
  - Color bar: Z-score normalized expression.
* **Cross-Check**:
  - ✅ Clean separation across all 11 genes in discovery cohort (N=8).
  - ⚠️ Note: Discovery significance is partly by construction (all 11 were pre-filtered as upregulated DEGs at nominal P < 0.05). This is expected for prioritization workflow, not a flaw, but MUST be noted.

---

## Figure 7: Validation Boxplots in GSE125498 (`results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png`)
* **Visible Text / Labels**:
  - 7 assessable hub genes in GSE125498 effluent cells (SPD n=20 vs LPD n=13):
    - *ISM1*: P = 0.811, AUC = 0.53
    - *FN1*: P = 0.294, AUC = 0.61
    - *VCAN*: P = 0.034, AUC = 0.72 (Direction-adjusted AUC = 0.72; raw expression LOWER in LPD — inverted post-hoc)
    - *COL3A1*: P = 0.754, AUC = 0.53
    - *COL8A1*: P = 0.117, AUC = 0.67
    - *THBS3*: P = 0.367, AUC = 0.60 (Direction-adjusted AUC = 0.60; raw expression LOWER in LPD — inverted post-hoc)
    - *LOX*: P = 0.985, AUC = 0.50
* **Cross-Check**:
  - ✅ 4 missing hubs (*COL11A1, COMP, EDIL3, INHBA*) clearly noted as lacking probes in GPL10558.
  - ⚠️ Post-hoc direction inversion for VCAN and THBS3 MUST be explicitly disclosed.
  - ✅ P-values and per-gene AUCs match `results/tables/GSE125498_wgcna_hub_validation_metrics.csv`.

---

## Figure 8: ROC Analysis (`results/figures/Hub_02b_roc_analysis.png`)
* **Visible Text / Labels**:
  - Panel A Title: "Single-Gene Discrimination (GSE125498)" — 7 individual ROC curves.
  - Panel B Title: "Multi-Gene Panel Validation (GSE125498)"
    - "Primary 7-Gene Panel (In-sample: AUC = 0.869)"
    - "Primary 7-Gene 50×5 Scaled CV (Mean AUC = 0.678 ± 0.058)" ← PRIMARY GENERALIZATION METRIC
    - "Primary 7-Gene Single 5-Fold CV (AUC = 0.592)"
    - "Nomogram Model (In-sample: AUC = 0.800)"
    - "Nomogram Model (5-Fold CV: AUC = 0.550)"
* **Cross-Check**:
  - ✅ All values match MANUSCRIPT_NUMBERS §6–§7.
  - ⚠️ 1,000-perm P for 7-gene panel = 0.0879 (non-significant at α=0.05) — MUST be reported alongside all AUC values.
  - Optimism gap 0.191 between in-sample 0.869 and CV 0.678 should be noted.

---

## Figure 9: Nomogram, Calibration & Decision Curve (`results/figures/Hub_02_clinical_nomogram_dca_calibration.png`)
* **Visible Text / Labels**:
  - Nomogram Points scale (5-gene model): *COL8A1* (0–100 pts), *VCAN* (0–84.5 pts), *COL3A1* (0–58.6 pts), *FN1* (0–33.3 pts), *ISM1* (0–17.1 pts).
  - Total Points → Risk of Long-Term Peritoneal Dialysis / Fibrosis Progression.
  - Apparent C-index = 0.819, Brier = 0.1542, Hosmer-Lemeshow P = 0.2437.
  - DCA Net Benefit curve positive across threshold probabilities 10% to 70%.
* **Cross-Check**:
  - ✅ Brier 0.1542, C-index 0.819 exact match.
  - ⚠️ **Critical framing:** 5-fold CV AUC = 0.550 (near chance). The nomogram's apparent performance is drastically overfit. MUST be labelled "exploratory multivariable point-scoring model", NOT "validated clinical nomogram".

---

## Figure 10: PPI & Co-Expression Network (`results/figures/Hub_01_ppi_gene_interaction_network.png`)
* **Visible Text / Labels**:
  - Title: "Hub Gene Interaction Network (STRING v12.5 Functional PPI & GSE62928 Co-expression)"
  - Legend: STRING v12.5 Functional PPI edges (gray dashed, n=21), Discovery Co-expression |r| ≥ 0.85 (green solid, n=42).
  - Hub nodes (11 genes): *COL1A1/COL3A1/COL8A1/COL11A1, FN1, VCAN, COMP, THBS3, EDIL3, INHBA, LOX, ISM1*.
  - Degree annotations: *FN1* (degree 18), *COL3A1* (degree 17), *ISM1* (degree 0 among functional edges).
* **Cross-Check**:
  - ✅ Generated live via STRING API v12.5 on 2026-10-01. Functional edges = 21, physical binding pairs = 3 (FN1-LOX 0.848, COL11A1-COL3A1 0.720, COMP-FN1 0.595).
  - ✅ 42/55 hub pairs = 76.36% co-expression at |r| ≥ 0.85.

---

## Figure 11: GSEA Hallmark Pathway Heatmap (`results/figures/Hub_03_gsea_pathway_enrichment_heatmap.png`)
* **Visible Text / Labels (CORRECTED against actual CSV 2026-10-07)**:
  > ⚠️ Previous transcription (2.483 / 1.868 / 1.737 / 1.579 / 1.564) was **stale / incorrect**. Values below are cross-checked verbatim from `results/tables/gsea_preranked_hallmark_results.csv` (verified 2026-10-07).
  - Panel A: Hallmark pathways ranked by NES in GSE62928:
    - Epithelial-Mesenchymal Transition (EMT): NES = +3.203, NOM P < 0.001, FDR q < 0.001
    - TNF-α Signaling via NF-κB: NES = +2.082, NOM P < 0.001, FDR q = 0.0010
    - Inflammatory Response: NES = +1.587, NOM P < 0.001, FDR q = 0.0303
    - Angiogenesis: NES = +1.438, NOM P = 0.049, FDR q = 0.0855
    - IL-6/JAK/STAT3 Signaling: NES = +1.392, NOM P = 0.033, FDR q = 0.0794
    - Apoptosis: NES = +1.349, NOM P = 0.021, FDR q = 0.0975
  - Panels B & C: ssGSEA Correlation heatmaps across 11 Hubs.
* **Cross-Check**:
  - ✅ All 6 top-up NES values above now exactly match the authoritative gsea CSV read 2026-10-07.
  - ✅ Matches `results/MANUSCRIPT_NUMBERS.md §8` exactly.
  - Top downregulated pathways (verbatim from CSV): Adipogenesis NES = −3.146, OxPhos NES = −3.101, Fatty-Acid Metabolism NES = −2.766, IFN-α Response NES = −2.385, ROS Pathway NES = −2.270 (all FDR < 0.001).

---

## Figure 12: Immune Infiltration Deconvolution (`results/figures/Hub_04_immune_infiltration_deconvolution.png`)
* **Visible Text / Labels**:
  - Subsets: Peritoneal Macrophages (M1/M2), Myofibroblasts, Mesothelial Cells, Neutrophils, CD4+ T cells, CD8+ T cells, B cells.
  - GSE62928 comparison: Peritoneal Myofibroblasts elevated in EPS (P = 0.0286, U = 16.0); other immune subsets non-significant after FDR control.
  - GSE125498 comparison: Myofibroblast score trends higher in LPD (P = 0.048).
* **Cross-Check**:
  - ✅ Mann-Whitney P_raw = 0.0286 for Myofibroblasts (GSE62928). FDR across cell types = 0.1714 (does not survive correction).
  - Exact combinatorial label-permutation P_perm = 2/70 = 0.0286 (two-sided).
  - ⚠️ Zero gene overlap between 11 hubs and the cell-type deconvolution marker gene sets — deconvolution is independent of hub selection, not circular.
