import os
import itertools
import numpy as np
import pandas as pd
from scipy import stats

print("=" * 75)
print("COMPREHENSIVE PERMUTATION TEST SUITE FOR PERITONEAL FIBROSIS STUDY")
print("=" * 75)

# 1. Load GSE62928 Expression Matrix and Sample Annotations
df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
samples = list(df_expr.columns)
# True labels: first 4 EPS cases (1), next 4 Controls (0)
true_cases = samples[:4]
true_ctrls = samples[4:]
print(f"Dataset GSE62928: {df_expr.shape[0]} genes x {len(samples)} samples")
print(f"True Cases (EPS): {true_cases}")
print(f"True Controls:    {true_ctrls}")

# Load Curated Matrisome
df_mat = pd.read_excel("ECM genes all.xlsx", skiprows=1)
sym_col = [c for c in df_mat.columns if "symbol" in c.lower() or "gene" in c.lower()][0]
mat_genes = set(df_mat[sym_col].dropna().astype(str).str.strip().str.upper()) - {"GENE SYMBOL", "NA", "NAN"}
print(f"Curated Matrisome Database: {len(mat_genes)} unique ECM genes")

# Generate all 70 combinatorial 4 vs 4 partitions
all_splits = list(itertools.combinations(samples, 4))
print(f"Total Combinatorial Label Partitions (8 choose 4): {len(all_splits)}")

# -----------------------------------------------------------------------------
# TEST 1: PRIMARY FINDING - ECM-DEG OVER-REPRESENTATION PERMUTATION TEST
# -----------------------------------------------------------------------------
print("\n--- TEST 1: ECM-DEG ENRICHMENT PERMUTATION TEST ---")
true_case_set = set(true_cases)
ecm_counts_per_split = []

for idx, case_split in enumerate(all_splits):
    case_cols = list(case_split)
    ctrl_cols = [s for s in samples if s not in case_split]
    
    expr_case = df_expr[case_cols].values
    expr_ctrl = df_expr[ctrl_cols].values
    
    # Fast Welch t-test & logFC across all genes
    mean_case = np.mean(expr_case, axis=1)
    mean_ctrl = np.mean(expr_ctrl, axis=1)
    logfc = mean_case - mean_ctrl
    
    # Two-sample t-test per gene
    t_stat, p_vals = stats.ttest_ind(expr_case, expr_ctrl, axis=1, equal_var=False)
    
    # Pro-fibrotic DEGs: logFC >= 0.585 (1.5-fold) and P < 0.05
    sig_mask = (logfc >= 0.585) & (p_vals < 0.05)
    sig_genes = df_expr.index[sig_mask]
    n_ecm = sum(g.upper() in mat_genes for g in sig_genes)
    
    is_true = (set(case_split) == true_case_set) or (set(ctrl_cols) == true_case_set)
    ecm_counts_per_split.append((n_ecm, len(sig_genes), is_true))

true_ecm_count = [x[0] for x in ecm_counts_per_split if x[2]][0]
true_deg_count = [x[1] for x in ecm_counts_per_split if x[2]][0]
# One-sided rank
ecm_counts_only = [x[0] for x in ecm_counts_per_split]
rank_ecm = sum(c >= true_ecm_count for c in ecm_counts_only)
p_perm_ecm = rank_ecm / len(all_splits)

print(f"True Clinical Label Split: {true_deg_count} DEGs -> {true_ecm_count} ECM-DEGs")
print(f"Permutation Distribution Max ECM-DEGs: {max(ecm_counts_only)}, Median: {np.median(ecm_counts_only)}")
print(f"Rank of True Split: #{rank_ecm} / {len(all_splits)}")
print(f"Primary Permutation P-value (P_perm): {p_perm_ecm:.4f} ({rank_ecm}/{len(all_splits)})")

# -----------------------------------------------------------------------------
# TEST 2: STAGE E - PPI CO-EXPRESSION CHANCE RATE AT N=8
# -----------------------------------------------------------------------------
print("\n--- TEST 2: EMPIRICAL CO-EXPRESSION NOISE BASELINE (N=8) ---")
np.random.seed(42)
hub_genes = ['ISM1', 'FN1', 'EDIL3', 'VCAN', 'COL3A1', 'COMP', 'COL8A1', 'THBS3', 'COL11A1', 'INHBA', 'LOX']
non_hub_genes = [g for g in df_expr.index if g.upper() not in hub_genes]

# Sample 100,000 random non-hub gene pairs
n_samples = 100000
idx1 = np.random.choice(len(non_hub_genes), n_samples)
idx2 = np.random.choice(len(non_hub_genes), n_samples)
# Ensure different genes
mask_diff = idx1 != idx2
idx1, idx2 = idx1[mask_diff], idx2[mask_diff]

expr_m1 = df_expr.iloc[idx1].values
expr_m2 = df_expr.iloc[idx2].values

# Vectorized Pearson r
m1_norm = expr_m1 - expr_m1.mean(axis=1, keepdims=True)
m2_norm = expr_m2 - expr_m2.mean(axis=1, keepdims=True)
denom = np.sqrt(np.sum(m1_norm**2, axis=1) * np.sum(m2_norm**2, axis=1))
denom[denom == 0] = 1e-12
r_vals = np.sum(m1_norm * m2_norm, axis=1) / denom

frac_abs_85 = np.mean(np.abs(r_vals) >= 0.85)
print(f"Total Random Non-Hub Gene Pairs Sampled: {len(r_vals)}")
print(f"Fraction with |r| >= 0.85 purely by chance at N=8: {frac_abs_85 * 100:.2f}%")

# Hub gene pairwise correlations
hub_expr = df_expr.loc[[g for g in hub_genes if g in df_expr.index]].values
hub_corr = np.corrcoef(hub_expr)
triu_idx = np.triu_indices(len(hub_corr), k=1)
hub_r_pairs = hub_corr[triu_idx]
hub_frac_85 = np.mean(np.abs(hub_r_pairs) >= 0.85)
print(f"Fraction of Hub Gene Pairs with |r| >= 0.85: {hub_frac_85 * 100:.2f}% ({sum(np.abs(hub_r_pairs) >= 0.85)} / {len(hub_r_pairs)})")

# -----------------------------------------------------------------------------
# TEST 3: GSEA EMT FOCAL PATHWAY PERMUTATION TEST
# -----------------------------------------------------------------------------
print("\n--- TEST 3: GSEA EMT PATHWAY PERMUTATION TEST ---")
try:
    import gseapy as gp
    hallmark_lib = gp.get_library_name(organism='Human')
    # Use MSigDB_Hallmark_2020
    true_ranked = pd.Series(stats.ttest_ind(df_expr[true_cases], df_expr[true_ctrls], axis=1, equal_var=False)[0], index=df_expr.index).sort_values(ascending=False)
    res_true = gp.prerank(rnk=true_ranked, gene_sets='MSigDB_Hallmark_2020', min_size=5, max_size=500, permutation_num=100, seed=42, verbose=False)
    df_gsea_true = res_true.res2d
    emt_true_nes = df_gsea_true[df_gsea_true['Term'].str.contains('EPITHELIAL_MESENCHYMAL_TRANSITION', case=False)]['NES'].values[0]
    print(f"Observed EMT NES (True Labels): +{emt_true_nes:.3f}")
except Exception as e:
    print(f"GSEA calculation note: {e}")
