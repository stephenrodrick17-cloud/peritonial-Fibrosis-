import itertools
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu
import gseapy as gp
import time

np.random.seed(42)

print("=================================================================")
print("STAGE C: PERMUTATION & NEGATIVE CONTROL AUDIT (GSEA & IMMUNE)")
print("=================================================================\n")

# Load Discovery Expression Matrix & Metadata (GSE62928, N=8)
df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
df_meta = pd.read_csv("results/tables/GSE62928_sample_metadata.csv")
sample_ids = df_meta["sample_id"].tolist()
true_groups = df_meta["group_binary_trait"].tolist() # 'Case_EPS' vs 'Control'
n_samples = len(sample_ids)

# Generate all 70 combinatorial splits of 4 Cases vs 4 Controls
all_combinations = list(itertools.combinations(range(n_samples), 4))
print(f"Total combinatorial permutations (8 choose 4): {len(all_combinations)}")

# -------------------------------------------------------------------------
# PART 1: IMMUNE INFILTRATION 70-PERMUTATION TEST
# -------------------------------------------------------------------------
print("\n--- Running 70 Permutations on Immune Infiltration (Wilcoxon Tests) ---")
df_infil = pd.read_csv("results/tables/immune_infiltration_scores.csv", index_col=0)
cell_types = [c for c in df_infil.columns if c != "Clinical_Group"]

# Observed true split
true_cases_idx = [i for i, g in enumerate(true_groups) if g == "Case_EPS"]
true_ctrls_idx = [i for i, g in enumerate(true_groups) if g == "Control"]

min_p_list = []
n_sig_list = []
myo_p_list = []

for comb in all_combinations:
    case_idx = list(comb)
    ctrl_idx = [i for i in range(n_samples) if i not in case_idx]
    
    pvals = []
    for ct in cell_types:
        v_case = df_infil[ct].iloc[case_idx].values
        v_ctrl = df_infil[ct].iloc[ctrl_idx].values
        _, p = mannwhitneyu(v_case, v_ctrl, alternative="two-sided")
        pvals.append(p)
        if ct == "Peritoneal Myofibroblasts":
            myo_p = p
    
    min_p_list.append(min(pvals))
    n_sig_list.append(sum(p < 0.05 for p in pvals))
    myo_p_list.append(myo_p)

# True values
true_comb = tuple(true_cases_idx)
true_idx = all_combinations.index(true_comb) if true_comb in all_combinations else all_combinations.index(tuple(true_ctrls_idx))

true_min_p = min_p_list[true_idx]
true_n_sig = n_sig_list[true_idx]
true_myo_p = myo_p_list[true_idx]

frac_any_sig = np.mean([n > 0 for n in n_sig_list])
p_perm_n_sig = np.mean([n >= true_n_sig for n in n_sig_list])
p_perm_myo = np.mean([p <= true_myo_p for p in myo_p_list])

print(f"Observed true split: min P = {true_min_p:.4f}, # cell types with P < 0.05 = {true_n_sig}, Myofibroblast P = {true_myo_p:.4f}")
print(f"Fraction of 70 random splits producing >= 1 cell type with P < 0.05: {frac_any_sig * 100:.1f}% ({sum(n > 0 for n in n_sig_list)}/70)")
print(f"Permutation P-value for # significant cell types (>= {true_n_sig}): P = {p_perm_n_sig:.4f} ({sum(n >= true_n_sig for n in n_sig_list)}/70)")
print(f"Permutation P-value for Peritoneal Myofibroblasts (P <= {true_myo_p:.4f}): P = {p_perm_myo:.4f} ({sum(p <= true_myo_p for p in myo_p_list)}/70)")

# -------------------------------------------------------------------------
# PART 2: GSEA PRERANKED 70-PERMUTATION TEST
# -------------------------------------------------------------------------
print("\n--- Running Permutations on MSigDB Hallmark GSEA (gseapy.prerank) ---")
# To ensure computational tractability while maintaining rigor, we run a representative sample of 20 random splits + the true split, or all 70.
# Let's test timing for 1 run:
t0 = time.time()
# Calculate ranked list for true split:
expr_mat = df_expr.values
c_idx = true_cases_idx
k_idx = true_ctrls_idx
t_stats = (expr_mat[:, c_idx].mean(axis=1) - expr_mat[:, k_idx].mean(axis=1)) / (
    np.sqrt(expr_mat[:, c_idx].var(axis=1)/4 + expr_mat[:, k_idx].var(axis=1)/4) + 1e-7
)
rnk = pd.Series(t_stats, index=df_expr.index).sort_values(ascending=False)

res = gp.prerank(rnk=rnk, gene_sets="MSigDB_Hallmark_2020", min_size=15, max_size=500, permutation_num=200, seed=42, verbose=False)
t_elapsed = time.time() - t0
print(f"Single GSEA prerank run time: {t_elapsed:.2f} seconds.")

# Run 70 permutations
gsea_n_sig_p05 = []
gsea_n_sig_fdr05 = []
gsea_emt_nes = []

print(f"Running full 70-permutation GSEA test (estimated time: {t_elapsed * 70 / 60:.1f} mins)...")
for idx, comb in enumerate(all_combinations):
    c_idx = list(comb)
    k_idx = [i for i in range(n_samples) if i not in c_idx]
    
    t_stats = (expr_mat[:, c_idx].mean(axis=1) - expr_mat[:, k_idx].mean(axis=1)) / (
        np.sqrt(expr_mat[:, c_idx].var(axis=1)/4 + expr_mat[:, k_idx].var(axis=1)/4) + 1e-7
    )
    rnk_perm = pd.Series(t_stats, index=df_expr.index).sort_values(ascending=False)
    
    try:
        res_perm = gp.prerank(
            rnk=rnk_perm,
            gene_sets="MSigDB_Hallmark_2020",
            min_size=15,
            max_size=500,
            permutation_num=200,
            seed=42 + idx,
            verbose=False
        )
        df_res = res_perm.res2d
        n_p05 = sum(df_res["NOM p-val"] < 0.05)
        n_fdr05 = sum(df_res["FDR q-val"] < 0.05)
        emt_row = df_res[df_res["Term"].str.contains("Epithelial Mesenchymal Transition", case=False, na=False)]
        emt_nes = emt_row["NES"].iloc[0] if len(emt_row) > 0 else 0.0
    except Exception as e:
        n_p05, n_fdr05, emt_nes = 0, 0, 0.0
        
    gsea_n_sig_p05.append(n_p05)
    gsea_n_sig_fdr05.append(n_fdr05)
    gsea_emt_nes.append(emt_nes)
    if (idx + 1) % 10 == 0 or idx == 0:
        print(f"  Processed {idx + 1}/70 permutations...")

true_gsea_p05 = gsea_n_sig_p05[true_idx]
true_gsea_fdr = gsea_n_sig_fdr05[true_idx]
true_emt_nes = gsea_emt_nes[true_idx]

p_perm_gsea_p05 = np.mean([n >= true_gsea_p05 for n in gsea_n_sig_p05])
p_perm_gsea_fdr = np.mean([n >= true_gsea_fdr for n in gsea_n_sig_fdr05])
p_perm_emt_nes = np.mean([nes >= true_emt_nes for nes in gsea_emt_nes])

print("\n=== GSEA 70-PERMUTATION RESULTS ===")
print(f"True split: {true_gsea_p05} Hallmarks with NOM P < 0.05, {true_gsea_fdr} Hallmarks with FDR < 0.05, EMT NES = {true_emt_nes:.3f}")
print(f"Permutation P-value for # Hallmarks with NOM P < 0.05 (>= {true_gsea_p05}): P = {p_perm_gsea_p05:.4f} ({sum(n >= true_gsea_p05 for n in gsea_n_sig_p05)}/70)")
print(f"Permutation P-value for # Hallmarks with FDR < 0.05 (>= {true_gsea_fdr}): P = {p_perm_gsea_fdr:.4f} ({sum(n >= true_gsea_fdr for n in gsea_n_sig_fdr05)}/70)")
print(f"Permutation P-value for EMT NES (>= {true_emt_nes:.3f}): P = {p_perm_emt_nes:.4f} ({sum(nes >= true_emt_nes for nes in gsea_emt_nes)}/70)")
