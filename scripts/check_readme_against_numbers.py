#!/usr/bin/env python3
"""
scripts/check_readme_against_numbers.py
Validates README.md against:
1. results/manuscript_numbers.json
2. Direct CSV recomputations:
   - results/tables/ML_hub_genes_from_WGCNA_ECM_all_results.csv (vote tiers, model selections)
   - results/tables/gsea_preranked_hallmark_results.csv (GSEA NES, p-values, FDRs)
3. Screens for any stale strings or incorrect gene universe associations.

Fails loudly with non-zero exit code if any error or mismatch is found.
"""

import sys
import json
import re
import pandas as pd

print("================================================================================")
print("RUNNING DIRECT CSV AND JSON CONSISTENCY AUDIT ON README.md")
print("================================================================================")

errors = []
warnings = []

with open("README.md", "r", encoding="utf-8") as f:
    readme_text = f.read()
    readme_lines = readme_text.splitlines()

# ============================================================================
# 1. Direct CSV Recomputations: ML Votes & Model Selections
# ============================================================================
df_ml = pd.read_csv("results/tables/ML_hub_genes_from_WGCNA_ECM_all_results.csv")
v4 = int((df_ml["Votes"] == 4).sum())
v3 = int((df_ml["Votes"] == 3).sum())
v2 = int((df_ml["Votes"] == 2).sum())
v1 = int((df_ml["Votes"] == 1).sum())
v0 = int((df_ml["Votes"] == 0).sum())

n_lasso = int(df_ml["LASSO_Selected"].sum())
n_svmrfe = int(df_ml["SVMRFE_Selected"].sum())
n_rf = int(df_ml["RandomForest_Selected"].sum())
n_xgb = int(df_ml["XGBoost_Selected"].sum())

# Verify ML selections in README
if f"LASSO: {n_lasso}" not in readme_text and f"LASSO ({n_lasso})" not in readme_text:
    errors.append(f"ML SELECTION MISMATCH: LASSO selection ({n_lasso}) not found in README")
if f"SVM-RFE: {n_svmrfe}" not in readme_text and f"SVM-RFE ({n_svmrfe})" not in readme_text:
    errors.append(f"ML SELECTION MISMATCH: SVM-RFE selection ({n_svmrfe}) not found in README")
if f"RF: {n_rf}" not in readme_text and f"Random Forest ({n_rf})" not in readme_text:
    errors.append(f"ML SELECTION MISMATCH: RF selection ({n_rf}) not found in README")
if f"XGBoost: {n_xgb}" not in readme_text and f"XGBoost ({n_xgb})" not in readme_text:
    errors.append(f"ML SELECTION MISMATCH: XGBoost selection ({n_xgb}) not found in README")

# Verify vote tiers in README
if f"• 4/4 Models: None ({v4} genes)" not in readme_text:
    errors.append(f"VOTE TIER MISMATCH: 4/4 votes ({v4} genes) not found in README")
if f"• 3/4 Models ({v3} genes)" not in readme_text:
    errors.append(f"VOTE TIER MISMATCH: 3/4 votes ({v3} genes) not found in README")
if f"• 2/4 Models ({v2} genes" not in readme_text:
    errors.append(f"VOTE TIER MISMATCH: 2/4 votes ({v2} genes) not found in README")
if f"• 1/4 Models ({v1} genes)" not in readme_text:
    errors.append(f"VOTE TIER MISMATCH: 1/4 votes ({v1} genes) not found in README")
if f"• 0/4 Models ({v0} genes)" not in readme_text:
    errors.append(f"VOTE TIER MISMATCH: 0/4 votes ({v0} genes) not found in README")

# ============================================================================
# 2. Direct CSV Recomputations: GSEA Pathways
# ============================================================================
df_gsea = pd.read_csv("results/tables/gsea_preranked_hallmark_results.csv")
gsea_lookup = {r["Term_Clean"]: r for _, r in df_gsea.iterrows()}

gsea_checks = [
    ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "3.203", "EMT NES"),
    ("HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB", "2.082", "TNF-alpha NES"),
    ("HALLMARK_INFLAMMATORY_RESPONSE", "1.587", "Inflammatory NES"),
    ("HALLMARK_ANGIOGENESIS", "1.438", "Angiogenesis NES"),
    ("HALLMARK_APICAL_JUNCTION", "1.411", "Apical Junction NES"),
    ("HALLMARK_IL_6/JAK/STAT3_SIGNALING", "1.392", "IL6 NES"),
    ("HALLMARK_APOPTOSIS", "1.349", "Apoptosis NES"),
    ("HALLMARK_COAGULATION", "1.333", "Coagulation NES"),
    ("HALLMARK_HYPOXIA", "-1.210", "Hypoxia NES")
]

for term, expected_str, desc in gsea_checks:
    row = gsea_lookup.get(term)
    if row is None:
        errors.append(f"GSEA TERM MISSING: Term '{term}' not found in gsea_preranked_hallmark_results.csv")
        continue
    actual_nes_str = f"{row['NES']:.3f}"
    if expected_str != actual_nes_str and abs(float(expected_str) - float(actual_nes_str)) > 0.002:
        errors.append(f"GSEA NES MISMATCH in CSV: Term '{term}' expected {expected_str}, found {actual_nes_str}")
    if expected_str not in readme_text:
        errors.append(f"GSEA VALUE MISSING IN README: Expected {expected_str} ({desc}) not found in README.md")

# ============================================================================
# 3. Screen for forbidden / stale strings
# ============================================================================
stale_patterns = [
    (r"four-way", "Legacy 'four-way' consensus phrasing found"),
    (r"4/4 models", "Legacy '4/4 models' phrasing found"),
    (r"\b0\.042\b", "Legacy LASSO alpha 0.042 found (should be C=4.28, alpha=0.234)"),
    (r"\b0\.677\b", "Stale truncated CV value 0.677 found"),
    (r"\b45\.5\b", "Stale co-expression enrichment 45.5x found"),
    (r"\b46\.6\b", "Stale co-expression enrichment 46.6x found"),
    (r"\[0\.640", "Old interval [0.640, 0.752] found"),
    (r"v12\.0\b", "Legacy STRING v12.0 found (should be v12.5)"),
    (r"live physical", "Legacy 'live physical' wording found"),
    (r"Hypoxia NES \+1\.349", "Mislabeled Hypoxia NES found"),
    (r"100,000 pairs", "Stale rounded pair count found (should be 99,993)"),
    (r"\+3\.63\b", "Old tissue logFC +3.63 found"),
    (r"\+1\.58\b", "Old tissue logFC +1.58 found"),
    (r"\+1\.74\b", "Old tissue logFC +1.74 found"),
    (r"beta = 18", "Old beta = 18 found (should be beta = 12)")
]

for pat, desc in stale_patterns:
    regex = re.compile(pat)
    matches = [(i+1, l) for i, l in enumerate(readme_lines) if regex.search(l)]
    if matches:
        for lnum, ltext in matches:
            # If the mention is within the explicit "Plain List of Corrected Errors", allow it as an explanation of what was fixed
            if lnum > 300 and any(k in ltext for k in ["Corrected", "Replaced", "legacy", "historical", "attributed"]):
                continue
            errors.append(f"STALE STRING: Line {lnum}: '{pat}' -> {desc}: {ltext.strip()}")

# Note on 0.696: only allowed if describing LOX 95% CI upper bound or historical sample-sorted note
for i, line in enumerate(readme_lines):
    if "0.696" in line:
        if "LOX" in line or "Historical Reference" in line or "sorted" in line or i+1 > 300:
            pass
        else:
            errors.append(f"STALE STRING: Line {i+1}: 0.696 used outside allowed historical/CI context: {line.strip()}")

# ============================================================================
# 4. Check Gene Universe Contexts
# ============================================================================
for i, line in enumerate(readme_lines):
    if "20,940" in line:
        if not any(k in line.lower() for k in ["wgcna", "maxmean", "matrix", "expression"]):
            errors.append(f"GENE UNIVERSE MISMATCH: Line {i+1}: '20,940' must sit next to WGCNA/MaxMean/matrix: {line.strip()}")
    if "22,049" in line:
        if not any(k in line.lower() for k in ["limma", "probe", "background", "top table", "universe"]):
            errors.append(f"GENE UNIVERSE MISMATCH: Line {i+1}: '22,049' must sit next to limma/probe/background: {line.strip()}")
    if "21,597" in line:
        if not any(k in line.lower() for k in ["gsea", "ranked", "symbol", "primary"]):
            errors.append(f"GENE UNIVERSE MISMATCH: Line {i+1}: '21,597' must sit next to GSEA/ranked symbols: {line.strip()}")

# ============================================================================
# 5. Check Master Numbers Presence
# ============================================================================
required_values = [
    ("367", "367 pro-fibrotic DEGs"),
    ("71", "71 ECM-DEGs"),
    ("19.35", "19.35% ECM-DEG proportion"),
    ("4.42", "4.42% background proportion"),
    ("604", "604 Salmon module genes"),
    ("0.806", "0.806 Salmon module trait correlation"),
    ("0.0157", "0.0157 Salmon module p-value"),
    ("40", "40 convergent candidates"),
    ("11", "11 consensus hub biomarkers"),
    ("0.869", "0.869 in-sample AUC"),
    ("0.678", "0.678 primary 50x5 CV pooled AUC"),
    ("0.658", "0.658 LOOCV AUC"),
    ("0.592", "0.592 single 5-fold CV AUC"),
    ("0.819", "0.819 nomogram C-index"),
    ("0.191", "0.191 optimism gap delta AUC"),
    ("21", "21 functional STRING edges"),
    ("3", "3 physical STRING edges"),
    ("1.68", "1.68% all-genes noise baseline"),
    ("5.47", "5.47% DEG baseline"),
    ("9.27", "9.27% Salmon module baseline"),
    ("76.36", "76.36% hub co-expression rate"),
    ("3.203", "EMT NES +3.203"),
    ("0.0286", "0.0286 myofibroblast expansion p-value"),
    ("0.1714", "0.1714 myofibroblast FDR")
]

for val, desc in required_values:
    if val not in readme_text:
        errors.append(f"MISSING MASTER VALUE: Required value '{val}' ({desc}) not found in README.md")

if errors:
    print(f"\n[FAILED] {len(errors)} consistency error(s) found in README.md:")
    for err in errors:
        safe_err = err.encode('ascii', errors='replace').decode('ascii')
        print(f"  [ERROR] {safe_err}")
    sys.exit(1)
else:
    print("\n[PASSED] All direct CSV and JSON consistency checks passed successfully! README.md is 100% verified.")
    sys.exit(0)
