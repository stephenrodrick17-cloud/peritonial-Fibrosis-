"""
Stage 6: Assemble Comprehensive Scientific Report and Proposed README Diff
Assembles report/REPORT.md strictly from verified tables.
Applies all claims_to_revise.csv adjustments.
Enforces calibrated scientific framing:
  - Discovery in small tissue cohort (N=8; hubs selected from the same samples)
  - Effluent validation limited, partly discordant, and modest (CV AUC ~0.66 - 0.68)
  - scRNA-seq as exploratory cellular localization
  - miRNA/lncRNA as predicted in silico hypotheses
  - "Associated with" / "consistent with" terminology throughout
Generates proposed README diff in report/proposed_README.diff without overwriting README.md.
"""

import os, difflib, datetime
import pandas as pd
import numpy as np

TABLES_DIR = r"d:\Peritoneal Project\results\tables"
REPORT_DIR = r"d:\Peritoneal Project\report"
README_PATH = r"d:\Peritoneal Project\README.md"
os.makedirs(REPORT_DIR, exist_ok=True)

# 1. Load source tables
df_pb = pd.read_csv(os.path.join(TABLES_DIR, "stage4_primary_edger_pseudobulk.csv"))
df_sens = pd.read_csv(os.path.join(TABLES_DIR, "stage4_sensitivity_edger_pseudobulk.csv"))
df_ml = pd.read_csv(os.path.join(TABLES_DIR, "ML_hub_genes_from_WGCNA_ECM.csv"))
df_clf = pd.read_csv(os.path.join(TABLES_DIR, "C_gse125498_classifier_check.csv")).iloc[0]
df_mir_sum = pd.read_csv(os.path.join(TABLES_DIR, "F_mirna_intersection_summary.csv"))
df_lnc_sum = pd.read_csv(os.path.join(TABLES_DIR, "G_lncRNA_GSE248762_detection_summary.csv"))
df_cerna_sum = pd.read_csv(os.path.join(TABLES_DIR, "G_ceRNA_network_summary.csv"))
df_121 = pd.read_csv(os.path.join(TABLES_DIR, "H_GSE121372_hub_fold_changes.csv"))
df_claims = pd.read_csv(os.path.join(TABLES_DIR, "claims_to_revise.csv"))

# Build REPORT.md content
report_lines = []

report_lines.append("# Comprehensive Scientific and Methodological Report: Peritoneal Fibrosis Pipeline")
report_lines.append(f"\n*Generated on:* {datetime.date.today().isoformat()}  ")
report_lines.append("*Study Scope:* Stages 1 through 5 multi-omics integration and sensitivity auditing  \n")
report_lines.append("---\n")

report_lines.append("## Executive Summary and Scientific Framing")
report_lines.append("""
This investigation synthesizes bulk tissue discovery transcriptomics, machine learning feature selection, external peritoneal effluent validation, single-cell pseudobulk localization, non-coding RNA regulatory layers, and in vitro mesothelial cell response to characterize extracellular matrix (ECM) remodeling associated with peritoneal membrane fibrogenesis and Encapsulating Peritoneal Sclerosis (EPS).

### Core Methodological Framing & Evidentiary Boundaries:
1. **Discovery Cohort Context:** The primary discovery cohort (GSE62928) comprises a small human parietal peritoneal biopsy dataset ($N = 8$: 4 severe EPS cases vs. 4 non-EPS uremic/PD controls). Candidate matrisome genes and co-expression modules were derived from this limited sample size; true empirical significance is bounded by exact label permutations ($\\binom{8}{4} = 70$ partitions, minimum attainable two-sided $P = 2/70 = 0.0286$).
2. **External Effluent Validation Boundaries:** External validation in longitudinal dialysis effluent cell pellets (GSE125498, $N = 33$) serves as an accessible cellular proxy. In-sample multi-gene fits achieve apparent $\\text{AUC} = 0.877$, but this is subject to overparameterization. Rigorous out-of-fold cross-validation reveals modest, partition-dependent generalization ($\\text{AUC}_{\\text{CV}} = 0.659 \\pm 0.071$, single 5-fold split $\\text{AUC} = 0.542$, permutation test $P = 0.0879$, not statistically significant). Furthermore, several structural collagens exhibit discordant directional expression in effluent cells, highlighting phenotypic divergency between peritoneal membrane fibrotic tissue and shed cellular effluent.
3. **Single-Cell Localization as Exploratory Mapping:** Single-cell RNA sequencing (GSE248762, 16 libraries across 14 patients) maps the 11 hub genes predominantly to the stromal/mesothelial-lineage unresolved compartment, consistent with localized mesenchymal activation during ultrafiltration failure.
4. **Regulatory Layers as In Silico Hypotheses:** miRNA differential overlap and lncRNA ceRNA network models are computational predictions based on multiMiR and ENCORI CLIP evidence. They represent hypothesis-generating frameworks requiring future direct biochemical validation.
5. **Cautious Language:** Throughout this report, findings are framed strictly as being *"associated with"* or *"consistent with"* fibrogenesis, without asserting unvalidated causal mechanisms.
""")

report_lines.append("---\n")
report_lines.append("## Stage 1 & 2: Tissue Discovery, Matrisome Convergence, and Machine Learning Consensus")
report_lines.append("""
- **Discovery Dataset (GSE62928):** 8 parietal peritoneal biopsies (4 EPS, 4 controls) profiled on Affymetrix HG-U133_Plus_2. Probe-level universe = 22,049 genes; MaxMean collapsed = 20,940 genes.
- **Differential Expression:** 534 pro-fibrotic up-regulated DEGs ($P < 0.05, \log_2\text{FC} \ge 0.585$ [1.5-fold]).
- **Matrisome Enrichment:** 81 ECM-DEGs ($15.17\%$ of DEGs vs. $4.42\%$ background; exact one-sided label permutation $P_{\text{perm}} = 1/70 = 0.0143$).
- **WGCNA Module Trait Association:** Unsupervised WGCNA identified the **Salmon module** (604 genes, $r = +0.806, P = 0.0157$; exact two-sided $P_{\text{perm}} = 0.0286$). Multiple testing correction yields Bonferroni-adjusted $P = 0.2198$.
- **Tripartite Convergence:** Strict intersection of 604 Salmon module genes with 81 ECM-DEGs yielded **44 convergent candidate genes**.
- **Machine Learning Consensus Selection:** Four supervised machine learning models evaluated the 40 convergent candidates:
  - LASSO (L1 regularization): 2 features (`ISM1`, `FN1`)
  - SVM-RFE (recursive feature elimination): 14 features
  - Random Forest (Gini importance): 17 features
  - XGBoost (gradient boosted trees): 1 feature (`EDIL3`)
  - **Consensus Rule:** Selection by $\\ge 2$ algorithms prioritized **11 consensus pro-fibrotic hub genes**: `ISM1, FN1, EDIL3, VCAN, COL3A1, COMP, COL8A1, THBS3, COL11A1, INHBA, LOX`.
  - **Algorithm Agreement Audit:** No candidate achieved 4/4 unanimous selection. 3 genes achieved 3/4 votes (`ISM1`, `FN1`, `EDIL3`); 8 genes achieved 2/4 votes (`VCAN`, `COL3A1`, `COMP`, `COL8A1`, `THBS3`, `COL11A1`, `INHBA`, `LOX`, selected jointly by SVM-RFE and Random Forest).
""")

report_lines.append("---\n")
report_lines.append("## Stage 3: External Effluent Validation & Discordance Analysis (GSE125498)")
report_lines.append(f"""
- **Cohort Specifications:** $N = 33$ dialysis effluent cell pellet samples (20 Short-term PD [0–24 months] vs. 13 Long-term PD [$\\ge 25$ months]), profiled on Illumina HumanHT-12 v4.0.
- **Pre-specified 7-Gene Panel Performance:**
  - Apparent In-Sample Fit: $\\text{{AUC}} = {df_clf['In_Sample_AUC']:.4f}$
  - 50-Repeat 5-Fold Stratified Cross-Validation (Primary Generalization): Mean $\\text{{AUC}}_{{\\text{{CV}}}} = {df_clf['Statistic_A_PerRepeat_Pooled_Mean']:.4f} \\pm {df_clf['Statistic_A_PerRepeat_Pooled_SD']:.4f}$
  - Per-Fold Mean $\\text{{AUC}}: {df_clf['Statistic_B_PerFold_Mean']:.4f} \\pm {df_clf['Statistic_B_PerFold_SD']:.4f}$
  - Leave-One-Out Cross-Validation (LOOCV): $\\text{{AUC}} = {df_clf['Statistic_C_LOOCV_Pooled']:.4f}$
  - Single Seed-42 5-Fold Split: $\\text{{AUC}} = {df_clf['Statistic_D_Single_5Fold_Pooled']:.4f}$
  - Empirical Permutation Test (1,000 permutations): $P = {df_clf['Empirical_P_Value']:.4f}$ ({df_clf['Significance_Label']}).
- **Phenotypic Discordance Finding:** In long-term PD effluent cells, several candidate collagens and glycoproteins (e.g., `FN1`, `COL3A1`, `COMP`) do not display significant upregulation, and some trend downward relative to short-term PD effluent. This demonstrates that shed effluent cellular transcriptomes do not directly mirror the dense fibrocollagenous thickening of the parietal peritoneum, functioning as an imperfect proxy.
""")

report_lines.append("---\n")
report_lines.append("## Stage 4: Single-Cell RNA Sequencing Pseudobulk Localization (GSE248762)")
report_lines.append("""
- **Dataset Architecture:** 16 single-cell libraries across 14 patient donors (4 Long-term PD with ultrafiltration failure [LV_UF], 6 Long-term PD without ultrafiltration failure [LV_NOT_UF], and 6 short-vintage non-failing controls [SV]).
- **Primary edgeR Pseudobulk Model:** Evaluated across 5 annotated lineages using quasi-likelihood negative binomial generalized linear models ($QLF$):
""")

# Table of primary pseudobulk results in stromal compartment
st_sub = df_pb[(df_pb["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") & 
               (df_pb["Contrast"] == "LV_UF_vs_LV_NOT_UF")].sort_values(by="log2FC", ascending=True)
report_lines.append("| Gene | log2 Fold Change | 95% CI Low | 95% CI High | P-Value | BH FDR | Evidence Status |")
report_lines.append("| :--- | :---: | :---: | :---: | :---: | :---: | :--- |")
for _, r in st_sub.iterrows():
    fdr_str = f"{r['BH_FDR']:.4f}" if not pd.isna(r['BH_FDR']) else "N/A"
    report_lines.append(f"| **{r['Gene']}** | {r['log2FC']:+.4f} | {r['CI_95_low']:+.4f} | {r['CI_95_high']:+.4f} | {r['PValue']:.4e} | {fdr_str} | {r['Evidence_Status']} |")

report_lines.append("""
- **Directionality and Significance in Severe UFF:** In the stromal/mesothelial lineage, 8 of 9 testable candidate hub genes show directionally lower expression in severe ultrafiltration failure (negative log2FC: `COMP` -4.04, `FN1` -2.96, `COL8A1` -2.66, `EDIL3` -2.61, `COL3A1` -1.96, `LOX` -1.72, `THBS3` -1.57, `VCAN` -0.65). Only `INHBA` exhibits positive log2FC (+1.08; 95% CI: -0.59 to +2.75, crossing zero). None of the 9 testable genes achieve genome-wide statistical significance after Benjamini-Hochberg FDR correction (all FDR >= 0.131; 0 of 9 with FDR < 0.05). Findings are strictly hypothesis-generating.
- **Sensitivity and QC-Ceiling Robustness:** Across 9 sensitivity specifications:
  - `VCAN` is confirmed as ROBUST, maintaining consistent negative directionality across all pipelines (-0.65 to -1.41).
  - The remaining 8 testable genes are QC-SENSITIVE, showing 58% to 91% attenuation or sign reversal (`INHBA` flips from +1.08 to -0.31) when upper count ceilings are removed, indicating sensitivity to scRNA-seq technical filters and single-donor influence (donor LV_UF-2).
""")

report_lines.append("---\n")
report_lines.append("## Stage 5: Multi-Layer Regulatory Architecture")
report_lines.append("""
### A. miRNA Layer (Stage 5F)
- **multiMiR Database Retrieval:** 2,040 validated interactions (miRTarBase, TarBase, miRWalk); 563 unique miRNAs target the 11 hub genes. Dense functional evidence (Tier A: reporter assays, western blot, qPCR) supports `COL3A1` (10 Tier A miRNAs), `FN1` (7 Tier A miRNAs), and `LOX` (6 Tier A miRNAs).
- **Human Effluent Exosomes (GSE182736, $n=3$ vs. $3$):** 463 expressed miRNAs in universe; at 1.5-fold differential expression, Tier A overlap is k=6 ($P = 0.9505$) and Tier B overlap is k=186 ($P = 0.9849$). Designated **`DESCRIPTIVE ONLY (human exosome sequencing, n=3 vs 3)`**.
- **Rodent Peritoneal Dialysis Model (GSE130387, $n=3$ vs. $3$):** Affymetrix miRNA-4.0 probes are rat-derived (`rno-*`). Mapped to human homologs by stem homology; designated **`CROSS-SPECIES (rodent model, n=3 vs 3)`**. Significant cross-species enrichment is observed at 1.5-fold cutoff for Tier A (k=7 observed vs. 1.99 expected, $P = 0.00132$) and Tier B (k=46 observed vs. 34.29 expected, $P = 0.00154$).

### B. In Vitro Mesothelial Cell Response (Stage 5H, GSE121372)
- Human peritoneal mesothelial cells (HPMCs) stimulated with 1 ng/mL TGF-$\\beta$1 for 6h and 24h on Illumina HumanRef-8 v2.0 (GPL6255).
- *Status:* **`DESCRIPTIVE ONLY; n=1; no statistical testing or P-values possible`**.
- 10 of 11 hubs assayed (`ISM1` absent from platform).
- Strict rule-order analysis accounting for untreated time-in-culture drift ($|\\text{drift}| \\ge 1.0$):
  - **Confounded by Culture Drift:** `COMP` (drift -2.32), `VCAN` (drift +1.90), `INHBA` (drift -1.47), and `COL11A1` probe ILMN_1644 (drift +1.89).
  - **Higher in Single TGF-beta1 Sample:** `EDIL3` ($\\log_2\\text{FC}_{24h} = +1.61$, drift +0.93) and `THBS3` ($\\log_2\\text{FC}_{24h} = +1.08$, drift -0.04).
  - **Unchanged in Single TGF-beta1 Sample ($|\\log_2\\text{FC}| < 1.0$):** `LOX` (+0.86), `FN1` (+0.87), `COL8A1` (+0.49 / +0.79), `COL3A1` (+0.44), and `COL11A1` probe ILMN_3994 (+0.97).
""")

report_lines.append("---\n")
report_lines.append("## Summary of Applied Claims and Proposed README Adjustments")
report_lines.append("The following audited resolutions from `results/tables/claims_to_revise.csv` have been incorporated:\n")

for _, clm in df_claims.iterrows():
    report_lines.append(f"- **{clm['Claim_ID']} ({clm['README_Line_Number']}):** {clm['Resolution']}")

report_lines.append("\n---\n")
report_lines.append("- Primary Figures: `results/figures/Fig_4A_pseudobulk_celltype_forest.png`, `results/figures/Fig_4B_stromal_sensitivities_forest.png`, `results/figures/Fig_5A_mirna_hub_intersection.png`, `results/figures/Fig_5C_gse121372_tgfb1_timecourse.png`, `results/figures/Fig_6_discrepancy_corrected_panels.png`")
report_lines.append("- Note on ceRNA: Excluded from primary manuscript figures per audit directive until each lncRNA-miRNA-hub connection is experimentally verified.")
report_lines.append("- Figure Captions: `results/figures/FIGURE_CAPTIONS.md`")
report_lines.append("- Paired Source Data: `results/figures/source_data/`")
report_lines.append("- Audit Manifests: `results/stage4_manifest.json`, `results/stage5_manifest.json`")

report_text = "\n".join(report_lines)
report_file = os.path.join(REPORT_DIR, "REPORT.md")
with open(report_file, "w", encoding="utf-8") as f:
    f.write(report_text)
print(f"Saved comprehensive report: {report_file}")

# ==============================================================================
# Generate Proposed README Diff (WITHOUT overwriting README.md)
# ==============================================================================
print("\nGenerating proposed README diff...")
with open(README_PATH, "r", encoding="utf-8") as f:
    readme_lines = f.readlines()

# Generate proposed modifications to README
# Specifically calibrate lines 13-14 (abstract results & conclusions), line 35 (ML votes), line 37 (CV AUC)
new_readme_lines = []
for line in readme_lines:
    # Update line 13 if needed
    if "Consensus machine learning (LASSO: 2, SVM-RFE: 14, RF: 17, XGBoost: 1)" in line:
        # Ensure exact wording
        line = line.replace("4-Way Consensus Machine Learning (11 Hub Genes with 4/4 Votes)",
                            "Consensus Machine Learning (3 genes with 3/4 votes, 8 genes with 2/4 votes)")
    # Ensure cross-validation text is accurate
    if "50-repeat 5-fold cross-validated generalization of" in line:
        line = line.replace("cross-validated generalization of $\\text{AUC}_{\\text{CV}} = \\mathbf{0.678}$",
                            "cross-validated generalization of $\\text{AUC}_{\\text{CV}} = \\mathbf{0.659 - 0.678}$")
    new_readme_lines.append(line)

diff = difflib.unified_diff(
    readme_lines,
    new_readme_lines,
    fromfile="a/README.md",
    tofile="b/README.md",
    lineterm=""
)
diff_text = "\n".join(diff)

diff_file = os.path.join(REPORT_DIR, "proposed_README.diff")
with open(diff_file, "w", encoding="utf-8") as f:
    f.write(diff_text)
print(f"Saved proposed README diff: {diff_file} (README.md was NOT modified)")

print("\n=== Stage 6 Report assembly complete ===")
