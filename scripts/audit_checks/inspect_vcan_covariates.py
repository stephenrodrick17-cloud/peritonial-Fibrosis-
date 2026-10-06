import gzip
import pandas as pd
import numpy as np
import statsmodels.api as sm

series_file = 'data/raw/GSE125498_series_matrix.txt.gz'
pids = ['ILMN_2396444', 'ILMN_1740015', 'ILMN_2134453', 'ILMN_1728639', 'ILMN_1739794', 'ILMN_1687301']
expr = {}
with gzip.open(series_file, 'rt', encoding='utf-8', errors='replace') as f:
    in_tab = False
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            in_tab = True
            continue
        if line.startswith('!series_matrix_table_end'):
            break
        if in_tab:
            parts = line.strip().split('\t')
            pid = parts[0].strip('"')
            if pid in pids:
                expr[pid] = [float(x) for x in parts[1:]]

meta = pd.read_csv('results/tables/GSE125498_sample_metadata.csv')
meta_map = dict(zip(meta['sample_id'], meta['stage_binary']))

# header has sample ids
with gzip.open(series_file, 'rt') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
    header = [s.strip('"') for s in f.readline().strip().split('\t')[1:]]

y = np.array([meta_map[s] for s in header])

df_p = pd.DataFrame(expr, index=header)
df_p['Group_LPD'] = y

print('Probe Means and Standard Deviations:')
for col in pids:
    print(f"  {col}: Mean = {df_p[col].mean():.3f}, SD = {df_p[col].std():.3f}")

print('\nCorrelations with VCAN (ILMN_1687301):')
for col in pids:
    r = df_p[col].corr(df_p['ILMN_1687301'])
    print(f"  {col} vs VCAN: r = {r:+.4f}")

print('\nCorrelations with Group (LPD vs SPD):')
for col in pids:
    r = df_p[col].corr(df_p['Group_LPD'])
    print(f"  {col} vs Group: r = {r:+.4f}")

# Compare Set A and Set B regressions
print("\n--- OLS Regression Models for VCAN ---")
X_A = sm.add_constant(df_p[['Group_LPD', 'ILMN_2396444', 'ILMN_2134453', 'ILMN_1739794']])
res_A = sm.OLS(df_p['ILMN_1687301'], X_A).fit()
print("Set A (CD14 ILMN_2396444, FCGR3B ILMN_2134453, CD3E ILMN_1739794):")
print(f"  Group coef (log2FC): {res_A.params['Group_LPD']:+.4f}, t = {res_A.tvalues['Group_LPD']:+.4f}, P = {res_A.pvalues['Group_LPD']:.4f}")
print("  Coefficients:\n", res_A.summary().tables[1])

X_B = sm.add_constant(df_p[['Group_LPD', 'ILMN_1740015', 'ILMN_1728639', 'ILMN_1739794']])
res_B = sm.OLS(df_p['ILMN_1687301'], X_B).fit()
print("\nSet B (CD14 ILMN_1740015, FCGR3B ILMN_1728639, CD3E ILMN_1739794):")
print(f"  Group coef (log2FC): {res_B.params['Group_LPD']:+.4f}, t = {res_B.tvalues['Group_LPD']:+.4f}, P = {res_B.pvalues['Group_LPD']:.4f}")
print("  Coefficients:\n", res_B.summary().tables[1])
