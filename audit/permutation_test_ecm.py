"""
Permutation Test on ECM Enrichment among Up-Regulated DEGs
Method 1: Gene-Level Collapsed Expression Matrix (results/tables/GSE62928_full_expression_matrix.csv, 20,940 genes)
Method 2: Probe-Level Series Matrix (54,715 probes -> 22,049 unique genes collapsed per permutation)
"""

import os
import gzip
import itertools
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, hypergeom

# 1. Load ECM Masterlist
df_ecm = pd.read_excel('ECM genes all.xlsx', skiprows=1)
sym_col = [c for c in df_ecm.columns if 'symbol' in c.lower() or 'gene' in c.lower()][0]
ecm_genes = set(df_ecm[sym_col].dropna().astype(str).str.strip().str.upper())
ecm_genes = {g for g in ecm_genes if g and g not in ['GENE SYMBOL', 'NA', 'NAN']}
print(f'Total Unique Curated ECM Genes: {len(ecm_genes)}')

# ------------------------------------------------------------------------------
# METHOD 1: On Gene-Level Expression Matrix (20,940 unique genes)
# ------------------------------------------------------------------------------
print('\n' + '=' * 75)
print('METHOD 1: GENE-LEVEL COLLAPSED EXPRESSION MATRIX (20,940 Genes x 8 Samples)')
print('=' * 75)

df_expr = pd.read_csv('results/tables/GSE62928_full_expression_matrix.csv', index_col=0)
print('Gene-level matrix shape:', df_expr.shape)
genes_m1 = np.array([g.strip().upper() for g in df_expr.index])
is_ecm_m1 = np.isin(genes_m1, list(ecm_genes))
mat_m1 = df_expr.values # (20940, 8)

all_combos = list(itertools.combinations(range(8), 4))
true_cases = [0, 1, 2, 3]

def run_m1(logfc_cut=0.80):
    counts_ecm = []
    counts_up = []
    for case_idx in all_combos:
        case_idx = list(case_idx)
        ctrl_idx = [i for i in range(8) if i not in case_idx]
        
        c_mat = mat_m1[:, case_idx]
        k_mat = mat_m1[:, ctrl_idx]
        
        mean_c = np.mean(c_mat, axis=1)
        mean_k = np.mean(k_mat, axis=1)
        lfc = mean_c - mean_k
        
        t_stat, p_val = ttest_ind(c_mat, k_mat, axis=1, equal_var=True)
        
        up_mask = (lfc >= logfc_cut) & (p_val < 0.05)
        n_up = np.sum(up_mask)
        n_ecm = np.sum(up_mask & is_ecm_m1)
        
        counts_up.append(n_up)
        counts_ecm.append(n_ecm)
        
    return np.array(counts_up), np.array(counts_ecm)

# Primary Cutoff (log2FC >= 0.80)
up_80_m1, ecm_80_m1 = run_m1(logfc_cut=0.80)
obs_up_80_m1 = up_80_m1[0]
obs_ecm_80_m1 = ecm_80_m1[0]
rank_80_m1 = np.sum(ecm_80_m1 > obs_ecm_80_m1) + 1
p_perm_80_m1 = np.mean(ecm_80_m1 >= obs_ecm_80_m1)

N_pop_m1 = len(genes_m1) # 20,940
K_ecm_m1 = np.sum(is_ecm_m1)
p_hyper_80_m1 = hypergeom.sf(obs_ecm_80_m1 - 1, N_pop_m1, K_ecm_m1, obs_up_80_m1)

print(f'Primary Cutoff (log2FC >= 0.80, P < 0.05):')
print(f'  - Observed Up-Regulated DEGs: {obs_up_80_m1}')
print(f'  - Observed ECM-DEGs: {obs_ecm_80_m1}')
print(f'  - Rank within 70 Permutations: {rank_80_m1} / 70')
print(f'  - Exact Permutation P-Value: {p_perm_80_m1:.4f} ({np.sum(ecm_80_m1 >= obs_ecm_80_m1)} / 70)')
print(f'  - Hypergeometric P-Value (Background N={N_pop_m1}, K={K_ecm_m1}): {p_hyper_80_m1:.4e}')
print(f'  - Permutation Distribution: Min = {np.min(ecm_80_m1)}, Median = {np.median(ecm_80_m1):.1f}, Mean = {np.mean(ecm_80_m1):.2f}, Max = {np.max(ecm_80_m1)}')
print('  - Full Frequency Table:')
for val, count in pd.Series(ecm_80_m1).value_counts().sort_index().items():
    print(f'      Count {val:2d}: {count:2d} permutations')

# Sensitivity Cutoff (log2FC >= 0.585)
up_585_m1, ecm_585_m1 = run_m1(logfc_cut=0.585)
obs_up_585_m1 = up_585_m1[0]
obs_ecm_585_m1 = ecm_585_m1[0]
rank_585_m1 = np.sum(ecm_585_m1 > obs_ecm_585_m1) + 1
p_perm_585_m1 = np.mean(ecm_585_m1 >= obs_ecm_585_m1)
p_hyper_585_m1 = hypergeom.sf(obs_ecm_585_m1 - 1, N_pop_m1, K_ecm_m1, obs_up_585_m1)

print(f'\nSensitivity Cutoff (log2FC >= 0.585, P < 0.05):')
print(f'  - Observed Up-Regulated DEGs: {obs_up_585_m1}')
print(f'  - Observed ECM-DEGs: {obs_ecm_585_m1}')
print(f'  - Rank within 70 Permutations: {rank_585_m1} / 70')
print(f'  - Exact Permutation P-Value: {p_perm_585_m1:.4f} ({np.sum(ecm_585_m1 >= obs_ecm_585_m1)} / 70)')
print(f'  - Hypergeometric P-Value: {p_hyper_585_m1:.4e}')

# ------------------------------------------------------------------------------
# METHOD 2: Probe-Level with Probe Collapsing (22,049 Unique Genes)
# ------------------------------------------------------------------------------
print('\n' + '=' * 75)
print('METHOD 2: PROBE-LEVEL DYNAMIC COLLAPSE (22,049 Unique Genes x 8 Samples)')
print('=' * 75)

# Load series matrix
with gzip.open('data/GSE62928_series_matrix.txt.gz', 'rt') as f:
    lines = []
    in_table = False
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            in_table = True
            continue
        if line.startswith('!series_matrix_table_end'):
            break
        if in_table:
            lines.append(line)

header = [h.replace('"', '').strip() for h in lines[0].strip().split('\t')]
data = [[val.replace('"', '').strip() for val in l.strip().split('\t')] for l in lines[1:]]
df_mat = pd.DataFrame(data, columns=header).set_index(header[0])
for c in df_mat.columns:
    df_mat[c] = pd.to_numeric(df_mat[c])

# Probe annotation mapping
df_tt = pd.read_csv('GSE62928.top.table.tsv', sep='\t')
df_annot = df_tt.dropna(subset=['Gene.symbol']).copy()
df_annot = df_annot[~df_annot['Gene.symbol'].isin(['', '---'])]

probe_gene_pairs = []
for _, r in df_annot.iterrows():
    pid = str(r['ID']).strip()
    for s in str(r['Gene.symbol']).split('///'):
        sc = s.strip().upper()
        if sc:
            probe_gene_pairs.append((pid, sc))

df_pairs = pd.DataFrame(probe_gene_pairs, columns=['ID', 'Gene'])
df_mat_annot = df_mat.loc[df_mat.index.intersection(df_pairs['ID'].unique())]
df_mat_merged = df_pairs.merge(df_mat_annot, left_on='ID', right_index=True)

samples = list(df_mat.columns)
expr_vals_m2 = df_mat_merged[samples].values
genes_m2 = df_mat_merged['Gene'].values

def run_m2(logfc_cut=0.80):
    counts_ecm = []
    counts_up = []
    for case_idx in all_combos:
        case_idx = list(case_idx)
        ctrl_idx = [i for i in range(8) if i not in case_idx]
        
        c_mat = expr_vals_m2[:, case_idx]
        k_mat = expr_vals_m2[:, ctrl_idx]
        
        mean_c = np.mean(c_mat, axis=1)
        mean_k = np.mean(k_mat, axis=1)
        lfc = mean_c - mean_k
        
        t_stat, p_val = ttest_ind(c_mat, k_mat, axis=1, equal_var=True)
        
        df_temp = pd.DataFrame({'Gene': genes_m2, 'logFC': lfc, 'P.Value': p_val})
        df_best = df_temp.sort_values('P.Value').drop_duplicates('Gene')
        
        up_degs = df_best[(df_best['logFC'] >= logfc_cut) & (df_best['P.Value'] < 0.05)]
        n_up = len(up_degs)
        n_ecm = np.sum(up_degs['Gene'].isin(ecm_genes))
        
        counts_up.append(n_up)
        counts_ecm.append(n_ecm)
        
    return np.array(counts_up), np.array(counts_ecm)

up_80_m2, ecm_80_m2 = run_m2(logfc_cut=0.80)
obs_up_80_m2 = up_80_m2[0]
obs_ecm_80_m2 = ecm_80_m2[0]
rank_80_m2 = np.sum(ecm_80_m2 > obs_ecm_80_m2) + 1
p_perm_80_m2 = np.mean(ecm_80_m2 >= obs_ecm_80_m2)
N_pop_m2 = 22049
K_ecm_m2 = 975
p_hyper_80_m2 = hypergeom.sf(obs_ecm_80_m2 - 1, N_pop_m2, K_ecm_m2, obs_up_80_m2)

print(f'Primary Cutoff (log2FC >= 0.80, P < 0.05):')
print(f'  - Observed Up-Regulated DEGs: {obs_up_80_m2}')
print(f'  - Observed ECM-DEGs: {obs_ecm_80_m2}')
print(f'  - Rank within 70 Permutations: {rank_80_m2} / 70')
print(f'  - Exact Permutation P-Value: {p_perm_80_m2:.4f} ({np.sum(ecm_80_m2 >= obs_ecm_80_m2)} / 70)')
print(f'  - Hypergeometric P-Value (Background N={N_pop_m2}, K={K_ecm_m2}): {p_hyper_80_m2:.4e}')
print(f'  - Permutation Distribution: Min = {np.min(ecm_80_m2)}, Median = {np.median(ecm_80_m2):.1f}, Mean = {np.mean(ecm_80_m2):.2f}, Max = {np.max(ecm_80_m2)}')
print('  - Full Frequency Table:')
for val, count in pd.Series(ecm_80_m2).value_counts().sort_index().items():
    print(f'      Count {val:2d}: {count:2d} permutations')

# Sensitivity Cutoff (log2FC >= 0.585)
up_585_m2, ecm_585_m2 = run_m2(logfc_cut=0.585)
obs_up_585_m2 = up_585_m2[0]
obs_ecm_585_m2 = ecm_585_m2[0]
rank_585_m2 = np.sum(ecm_585_m2 > obs_ecm_585_m2) + 1
p_perm_585_m2 = np.mean(ecm_585_m2 >= obs_ecm_585_m2)
p_hyper_585_m2 = hypergeom.sf(obs_ecm_585_m2 - 1, N_pop_m2, K_ecm_m2, obs_up_585_m2)

print(f'\nSensitivity Cutoff (log2FC >= 0.585, P < 0.05):')
print(f'  - Observed Up-Regulated DEGs: {obs_up_585_m2}')
print(f'  - Observed ECM-DEGs: {obs_ecm_585_m2}')
print(f'  - Rank within 70 Permutations: {rank_585_m2} / 70')
print(f'  - Exact Permutation P-Value: {p_perm_585_m2:.4f} ({np.sum(ecm_585_m2 >= obs_ecm_585_m2)} / 70)')
print(f'  - Hypergeometric P-Value: {p_hyper_585_m2:.4e}')
