import gzip, os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats

# 1. Load series matrix expression
series_file = "data/raw/GSE125498_series_matrix.txt.gz"
with gzip.open(series_file, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
    header = f.readline().strip().split("\t")
    sample_ids = [s.strip('"') for s in header[1:]]
    expr_dict = {}
    for line in f:
        if line.startswith("!series_matrix_table_end"):
            break
        parts = line.strip().split("\t")
        pid = parts[0].strip('"')
        expr_dict[pid] = [float(x) for x in parts[1:]]

df_expr = pd.DataFrame(expr_dict, index=sample_ids).T # probes x samples

# 2. Metadata
df_meta = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
meta_map = dict(zip(df_meta["sample_id"], df_meta["stage_binary"]))
y = np.array([meta_map[s] for s in df_expr.columns]) # 1=LPD, 0=SPD

top = pd.read_csv('Validation/GSE125498.top.table.tsv', sep='\t', index_col=0)

target_probes = {
    'VCAN': ['ILMN_1687301'],
    'FN1': ['ILMN_2366463', 'ILMN_1778237'],
    'COL8A1': ['ILMN_1685433', 'ILMN_2402392']
}

print(f"Total samples: {len(y)} ({sum(y == 1)} LPD, {sum(y == 0)} SPD)")
print(f"Gene\tProbe\tRole\tlog2FC\tWelch_P\tLimma_P\tOLS_unadj_P")

rows = []
for gene, pids in target_probes.items():
    for pid in pids:
        vals = df_expr.loc[pid].values
        spd_vals = vals[y == 0]
        lpd_vals = vals[y == 1]
        
        fc = np.mean(lpd_vals) - np.mean(spd_vals)
        t_welch, p_welch = stats.ttest_ind(lpd_vals, spd_vals, equal_var=False)
        limma_p = top.loc[pid, 'P.Value'] if pid in top.index else np.nan
        
        X = sm.add_constant(y)
        model = sm.OLS(vals, X).fit()
        ols_p = model.pvalues[1]
        
        role = 'Lowest-P' if pid in ['ILMN_1778237', 'ILMN_2402392'] else ('MaxMean' if pid in ['ILMN_2366463', 'ILMN_1685433'] else 'Single')
        rows.append({
            'Gene': gene,
            'Probe': pid,
            'Role': role,
            'log2FC': fc,
            'Welch_P': p_welch,
            'Limma_P': limma_p,
            'OLS_unadj_P': ols_p
        })
        print(f"{gene}\t{pid}\t{role}\t{fc:+.4f}\t{p_welch:.4f}\t{limma_p:.4f}\t{ols_p:.4f}")

df_res = pd.DataFrame(rows)
df_res.to_csv("results/tables/reconciled_p_values_vcan_fn1_col8a1.csv", index=False)
