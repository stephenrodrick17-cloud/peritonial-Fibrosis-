import gzip
import pandas as pd
import numpy as np
import statsmodels.api as sm

series_file = 'data/raw/GSE125498_series_matrix.txt.gz'
pids = ['ILMN_2396444', 'ILMN_1740015', 'ILMN_2134453', 'ILMN_1728639', 'ILMN_1739794', 'ILMN_1687301']
expr = {}

with gzip.open(series_file, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            header = [s.strip('"') for s in f.readline().strip().split('\t')[1:]]
            break
    for line in f:
        if line.startswith('!series_matrix_table_end'):
            break
        parts = line.strip().split('\t')
        pid = parts[0].strip('"')
        if pid in pids:
            expr[pid] = [float(x) for x in parts[1:]]

meta = pd.read_csv('results/tables/GSE125498_sample_metadata.csv')
meta_map = dict(zip(meta['sample_id'], meta['stage_binary']))
df_p = pd.DataFrame(expr, index=header)
df_p['Group_LPD'] = [meta_map[s] for s in header]

y = df_p['ILMN_1687301']

print("=== LEAVE-ONE-COVARIATE-OUT REGRESSION FOR VCAN ===")

for sname, cov_dict in [
    ('Set A', {'CD14': 'ILMN_2396444', 'FCGR3B': 'ILMN_2134453', 'CD3E': 'ILMN_1739794'}),
    ('Set B', {'CD14': 'ILMN_1740015', 'FCGR3B': 'ILMN_1728639', 'CD3E': 'ILMN_1739794'})
]:
    print(f"\n--- {sname} ---")
    cov_names = list(cov_dict.keys())
    cov_pids = list(cov_dict.values())
    
    # Full joint model
    X_full = sm.add_constant(df_p[['Group_LPD'] + cov_pids])
    res_full = sm.OLS(y, X_full).fit()
    print(f"Full Joint (all 3 covs): Group log2FC = {res_full.params['Group_LPD']:+.4f}, t = {res_full.tvalues['Group_LPD']:+.4f}, OLS P = {res_full.pvalues['Group_LPD']:.4f}")
    
    # Drop each covariate in turn
    for drop_name, drop_pid in cov_dict.items():
        rem_pids = [p for p in cov_pids if p != drop_pid]
        X_sub = sm.add_constant(df_p[['Group_LPD'] + rem_pids])
        res_sub = sm.OLS(y, X_sub).fit()
        print(f"  Drop {drop_name} ({drop_pid}): Group log2FC = {res_sub.params['Group_LPD']:+.4f}, t = {res_sub.tvalues['Group_LPD']:+.4f}, OLS P = {res_sub.pvalues['Group_LPD']:.4f}")
