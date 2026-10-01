# Changelog: Numerical & Methodological Consistency Fix

**Branch:** `readme-consistency-fix`  
**Date:** 2026-10-01  
**Audit Scope:** Full repository numerical consistency, documentation harmonization, live API provenance, and permutation rigor.

---

## 1. Summary of Changes

### A. STRING Protein-Protein Interaction Network Provenance (Script 07)
- **Live Query Enforcement:** Updated `07_gene_interaction_network.py` to raise an explicit `RuntimeError` on connection failure rather than falling back to an internal hardcoded cache.
- **Raw API Artifact Export:** Script 07 now saves the live raw TSV response from the STRING database directly to `results/tables/string_live_edges_20261001.tsv`.
- **Regenerated PPI Artifacts:** Re-executed script 07 against STRING v12.5 (retrieved 2026-10-01), regenerating `results/tables/hub_genes_ppi_centrality_metrics.csv`, `results/figures/Hub_01_ppi_gene_interaction_network.png`, and `results/figures/Hub_01_ppi_gene_interaction_network.pdf`.
- **Documentation Updated:** `README.md` updated to report **21 verified functional association edges** (score $\ge 0.400$), listing top edges (`FN1 - LOX` 0.953, `COL3A1 - FN1` 0.911, `COL11A1 - COL3A1` 0.877, `COL3A1 - LOX` 0.802), removing "live physical" terminology, and noting physical interaction network results (3 edges).

### B. GSEA Hallmark Pathway Clarification & Hypoxia Correction
- **Hypoxia Mislabling Corrected:** Rectified historical documentation error where Apoptosis statistics ($\text{NES} = +1.3490, \text{FDR} = 0.0975$) had been erroneously attributed to Hypoxia. `HALLMARK_HYPOXIA` is confirmed in `results/tables/gsea_preranked_hallmark_results.csv` as downregulated and non-significant ($\text{NES} = -1.2097, \text{NOM } P = 0.1407, \text{FDR } q = 0.2026$).
- **Downregulated Hallmarks Added:** Explicitly listed top suppressed metabolic and immune hallmarks: Adipogenesis ($\text{NES} = -3.1463, \text{FDR} < 0.001$), Oxidative Phosphorylation ($\text{NES} = -3.1010, \text{FDR} < 0.001$), Fatty Acid Metabolism ($\text{NES} = -2.7658, \text{FDR} < 0.001$), Interferon Alpha Response ($\text{NES} = -2.3846, \text{FDR} < 0.001$).

### C. Permutation Test Sidedness & Minimum Achievable P-Values
- **ECM Over-representation:** Explicitly designated as one-sided count rank ($P_{\text{perm}} = 1/70 = 0.0143$).
- **WGCNA Module Trait:** Reported both one-sided ($r \ge +0.806, P_{\text{perm}} = 1/70 = 0.0143$) and two-sided ($|r| \ge 0.806, P_{\text{perm}} = 2/70 = 0.0286$). Noted that minimum achievable two-sided $P$ at $N=8$ is $2/70 = 0.0286$.
- **GSEA EMT:** One-sided NES rank ($P_{\text{perm}} = 1/70 = 0.0143$).
- **Immune Myofibroblasts:** Two-sided rank ($P_{\text{perm}} = 2/70 = 0.0286$), noting exact equality with the two-sided Wilcoxon Mann-Whitney test ($P = 0.028571$).

### D. Cross-Validation Harmonization for 7-Gene Panel (GSE125498)
- **Primary Generalization Metric:** Designated 50-repeat 5-fold stratified CV with `StandardScaler` inside a `Pipeline` as primary: Per-fold $\text{AUC}_{\text{CV}} = 0.701 \pm 0.203$ ($N=250$), Per-repeat pooled $\text{AUC} = 0.678 \pm 0.058$ ($N=50$, range: $0.512 - 0.773$).
- **Supplementary CV Models Reported:**
  - 50x5 CV without scaling: Per-fold $\text{AUC}_{\text{CV}} = 0.706 \pm 0.209$, Per-repeat pooled $\text{AUC} = 0.684 \pm 0.063$ (range: $0.523 - 0.781$).
  - LOOCV with pipeline scaling: Pooled $\text{AUC} = 0.658$.
  - Single seed-42 5-fold CV with pipeline scaling: Per-fold Mean $\text{AUC} = 0.617 \pm 0.061$, Pooled $\text{AUC} = 0.592$.
- **Figure 8 Caption:** Synchronized Figure 8 caption with primary and single-split metrics.

### E. Co-Expression Noise Baselines
- **Empirical Threshold Rates ($|r| \ge 0.85$ in GSE62928, $N=8$):**
  - All non-hub random pairs ($N=99,993$): $1.68\%$ (`audit/stage_e_ppi_baseline.py` seed=42).
  - Among 367 Pro-Fibrotic DEGs ($N=55,945$ pairs): $5.47\%$ ($14.0\times$ enrichment for hubs vs DEG baseline).
  - Among 604 Salmon-Module Genes ($N=182,106$ pairs): $9.27\%$ ($8.2\times$ enrichment for hubs vs Salmon-module baseline).
  - Among 11 Consensus Hub Genes ($N=55$ pairs): $76.36\%$ ($42/55$ pairs). High hub co-expression is partly by construction because hubs were prioritized from the co-expressed WGCNA Salmon module.

### F. Gene Universes & Probe Collapse
- **Gene Universe Distinction:** Clarified that probe-level Limma DE and Volcano Plot analyze 22,049 genes (best-probe-by-P collapse), whereas WGCNA and full matrix analyze 20,940 genes (MaxMean collapse), and GSEA evaluates 21,597 primary symbols.
- **Probe Collapse Note:** Added explicit methodological disclosure regarding probe collapse sensitivity.

### G. Machine Learning Consensus Harmonization
- **Algorithm Feature Selections:** Clarified LASSO (2 features, $C=4.2813, \alpha=0.2336$), SVM-RFE (14 features), Random Forest (17 features), and XGBoost (1 feature, `EDIL3`).
- **Vote Tiers:** 4/4 votes = 0 genes; 3/4 votes = 3 genes (`ISM1`, `FN1`, `EDIL3`); 2/4 votes = 8 genes (`VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX` — all 8 selected by SVM-RFE + RF).

### H. Discovery Fold-Changes in Discordance Table
- **True Limma Values Transcribed:** Replaced discordance table discovery fold-changes with true values from `GSE62928_gene_level_toptable.csv` (`COL8A1 +2.68, FN1 +1.93, COL3A1 +2.84, ISM1 +1.90, LOX +2.16, VCAN +2.75, THBS3 +1.16, COMP +4.08, COL11A1 +3.79`).

### I. Programmatic Manuscript Numbers
- **Automated Generation:** Created `scratch/generate_manuscript_numbers.py` to programmatically extract and format all pipeline numbers directly from CSV tables into `results/MANUSCRIPT_NUMBERS.md`.
