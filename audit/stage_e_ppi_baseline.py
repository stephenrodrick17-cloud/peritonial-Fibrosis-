import numpy as np
import pandas as pd
from scipy.stats import pearsonr

print("=================================================================")
print("STAGE E: PPI CO-EXPRESSION CORRELATION BASELINE AT N=8")
print("=================================================================\n")

# Load GSE62928 Expression Matrix (N=8, 20,940 genes)
df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
hub_genes = ["FN1", "COL3A1", "COL11A1", "COL8A1", "VCAN", "COMP", "THBS3", "EDIL3", "LOX", "INHBA", "ISM1"]

# Filter out the 11 hub genes
non_hub_genes = [g for g in df_expr.index if g not in hub_genes]
print(f"Total non-hub genes available in GSE62928: {len(non_hub_genes)}")

# Extract matrix for non-hub genes
mat_non_hub = df_expr.loc[non_hub_genes].values
n_non_hub = len(non_hub_genes)

# Standardize rows with sample std (ddof=1) so dot product / (N - 1) is exact sample Pearson r
mat_std = (mat_non_hub - mat_non_hub.mean(axis=1, keepdims=True)) / np.std(mat_non_hub, axis=1, keepdims=True, ddof=1)

# Randomly sample 100,000 pairs of non-hub genes
np.random.seed(42)
N_PAIRS = 100000
idx_A = np.random.randint(0, n_non_hub, size=N_PAIRS)
idx_B = np.random.randint(0, n_non_hub, size=N_PAIRS)

# Ensure distinct indices
mask_diff = idx_A != idx_B
idx_A = idx_A[mask_diff]
idx_B = idx_B[mask_diff]
n_actual_pairs = len(idx_A)

# Compute exact Pearson r for all pairs: r = sum(z_A * z_B) / (8 - 1)
r_vals = np.sum(mat_std[idx_A] * mat_std[idx_B], axis=1) / 7.0
r_vals = np.clip(r_vals, -1.0, 1.0)

# Calculate baseline statistics
frac_abs_r_85 = np.mean(np.abs(r_vals) >= 0.85)
frac_pos_r_85 = np.mean(r_vals >= 0.85)
frac_abs_r_70 = np.mean(np.abs(r_vals) >= 0.70)
mean_abs_r = np.mean(np.abs(r_vals))

print(f"Computed Pearson r for {n_actual_pairs:,} randomly sampled non-hub gene pairs in GSE62928 (N=8):")
print(f"  * Mean |r| across random pairs: {mean_abs_r:.4f}")
print(f"  * Fraction of random pairs with |r| >= 0.70: {frac_abs_r_70 * 100:.2f}% ({int(frac_abs_r_70 * n_actual_pairs):,}/{n_actual_pairs:,})")
print(f"  * Fraction of random pairs with |r| >= 0.85: {frac_abs_r_85 * 100:.2f}% ({int(frac_abs_r_85 * n_actual_pairs):,}/{n_actual_pairs:,})")
print(f"  * Fraction of random pairs with r >= +0.85:  {frac_pos_r_85 * 100:.2f}% ({int(frac_pos_r_85 * n_actual_pairs):,}/{n_actual_pairs:,})")

# Check hub-gene co-expression pairs among the 11 hub genes (55 total possible pairs)
hub_expr = df_expr.loc[hub_genes].values
hub_corr = np.corrcoef(hub_expr)

hub_pairs_r = []
for i in range(len(hub_genes)):
    for j in range(i + 1, len(hub_genes)):
        hub_pairs_r.append(hub_corr[i, j])

hub_pairs_r = np.array(hub_pairs_r)
hub_frac_abs_85 = np.mean(np.abs(hub_pairs_r) >= 0.85)
hub_n_abs_85 = np.sum(np.abs(hub_pairs_r) >= 0.85)

print("\n--- Hub Gene Co-expression (55 total pairs among 11 hub genes) ---")
print(f"  * Total possible hub-hub pairs: {len(hub_pairs_r)}")
print(f"  * Hub pairs with |r| >= 0.85: {hub_n_abs_85}/55 ({hub_frac_abs_85 * 100:.1f}%)")
print(f"  * Mean |r| among hub genes: {np.mean(np.abs(hub_pairs_r)):.4f}")
