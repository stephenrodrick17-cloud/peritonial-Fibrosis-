import os
import json
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.makedirs("results/revision", exist_ok=True)
os.makedirs("scripts", exist_ok=True)
OUTPUT_SUFFIX = "_taskA_fig5"

# 1. Parse GPL13158 mapping
probe_to_gene = {}
with open('data/raw/GPL13158.soft.gz', 'r', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('!platform_table_begin'):
            break
    header = f.readline().strip().split('\t')
    id_idx = header.index('ID')
    sym_idx = header.index('Gene Symbol')
    for line in f:
        if line.startswith('!platform_table_end'):
            break
        parts = line.strip().split('\t')
        if len(parts) > sym_idx:
            pid = parts[id_idx].strip()
            sym = parts[sym_idx].strip()
            if sym and sym not in ['---', 'NA', '']:
                probe_to_gene[pid] = sym.split('///')[0].strip()

# 2. Parse GSE62928 expression matrix
matrix_gz = 'data/raw/GSE62928_series_matrix.txt.gz'
with gzip.open(matrix_gz, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
    expr_df = pd.read_csv(f, sep='\t')

expr_df = expr_df.rename(columns={expr_df.columns[0]: 'ID_REF'})
expr_df['Gene'] = expr_df['ID_REF'].map(probe_to_gene)
expr_df = expr_df[expr_df['Gene'].notna()].copy()

sample_cols = [c for c in expr_df.columns if c.startswith('GSM')]
expr_gene = expr_df.groupby('Gene')[sample_cols].mean()

eps_cols = ['GSM1536406', 'GSM1536407', 'GSM1536408', 'GSM1536409']
ctrl_cols = ['GSM1536410', 'GSM1536411', 'GSM1536412', 'GSM1536413']
all_cols = eps_cols + ctrl_cols
y_true = np.array([1, 1, 1, 1, 0, 0, 0, 0])

# Load discovery statistics
discovery_deg = pd.read_csv('results/tables/GSE125498_all_results.csv').set_index('Gene')

hubs_7 = ['FLT3LG', 'TNFSF15', 'LTB', 'ADAM19', 'SERPINA10', 'EBI3', 'CXCL14']
per_gene_records = []

for g in hubs_7:
    v_eps = expr_gene.loc[g, eps_cols].values.astype(float)
    v_ctrl = expr_gene.loc[g, ctrl_cols].values.astype(float)
    
    m_eps = float(np.mean(v_eps))
    m_ctrl = float(np.mean(v_ctrl))
    val_logfc = float(m_eps - m_ctrl)
    
    s_eps = float(np.var(v_eps, ddof=1))
    s_ctrl = float(np.var(v_ctrl, ddof=1))
    se_diff = float(np.sqrt(s_eps / 4.0 + s_ctrl / 4.0))
    
    # Welch degrees of freedom
    num = (s_eps / 4.0 + s_ctrl / 4.0) ** 2
    den = ((s_eps / 4.0)**2) / 3.0 + ((s_ctrl / 4.0)**2) / 3.0
    df_welch = float(num / den) if den > 0 else 6.0
    
    t_stat = val_logfc / max(1e-9, se_diff)
    p_val = float(2.0 * stats.t.sf(abs(t_stat), df=df_welch))
    
    ci_t = stats.t.ppf(0.975, df=df_welch)
    ci_lower = float(val_logfc - ci_t * se_diff)
    ci_upper = float(val_logfc + ci_t * se_diff)
    
    disc_logfc = float(discovery_deg.loc[g, 'logFC'])
    disc_p = float(discovery_deg.loc[g, 'P.Value'])
    
    is_concordant = (disc_logfc * val_logfc) > 0
    
    per_gene_records.append({
        'Gene': g,
        'Discovery_log2FC': round(disc_logfc, 4),
        'Discovery_P': disc_p,
        'Validation_log2FC': round(val_logfc, 4),
        'Validation_CI95_Lower': round(ci_lower, 4),
        'Validation_CI95_Upper': round(ci_upper, 4),
        'Validation_SE': round(se_diff, 4),
        'Validation_Welch_P': round(p_val, 4),
        'Direction_Concordant': 'Yes' if is_concordant else 'No'
    })

df_per_gene = pd.DataFrame(per_gene_records)
df_per_gene.to_csv(f'results/revision/gse62928_per_gene_concordance{OUTPUT_SUFFIX}.csv', index=False)

# Sign concordance test
n_concordant = int((df_per_gene['Direction_Concordant'] == 'Yes').sum())
n_total = len(hubs_7)
# Binomial test H0: p=0.5
binom_res = stats.binomtest(n_concordant, n_total, p=0.5, alternative='two-sided')
binom_p = float(binom_res.pvalue)

# Rank-based / Composite score test on N=8
hub_mat = expr_gene.loc[hubs_7, all_cols].values # (7, 8)
hub_z = (hub_mat - hub_mat.mean(axis=1, keepdims=True)) / np.clip(hub_mat.std(axis=1, ddof=1, keepdims=True), 1e-9, None)
hub_score = hub_z.mean(axis=0)
hub_score_source = pd.DataFrame({
    "sample_id": all_cols,
    "outcome": ["EPS"] * len(eps_cols) + ["Control"] * len(ctrl_cols),
    "stage_binary": y_true,
    "fixed_hub_score": hub_score,
})
for gene, z_score in zip(hubs_7, hub_z):
    hub_score_source[f"{gene}_z"] = z_score
hub_score_source.to_csv(
    f"results/revision/gse62928_hub_score_by_sample{OUTPUT_SUFFIX}.csv",
    index=False,
)

hub_auc = float(roc_auc_score(y_true, hub_score))
u_stat, mw_p = stats.mannwhitneyu(hub_score[:4], hub_score[4:], alternative='greater')
u_two_sided, mw_two_p = stats.mannwhitneyu(hub_score[:4], hub_score[4:], alternative='two-sided')

# 1,000 random 7-gene sets null
np.random.seed(42)
all_genes = list(expr_gene.index)
null_aucs = []
for _ in range(1000):
    rand_genes = np.random.choice(all_genes, 7, replace=False)
    rand_mat = expr_gene.loc[rand_genes, all_cols].values
    rand_z = (rand_mat - rand_mat.mean(axis=1, keepdims=True)) / np.clip(rand_mat.std(axis=1, ddof=1, keepdims=True), 1e-9, None)
    rand_score = rand_z.mean(axis=0)
    a = roc_auc_score(y_true, rand_score)
    null_aucs.append(float(a))

null_aucs = np.array(null_aucs)
emp_p = float(np.mean(null_aucs >= hub_auc))
n_null_at_least_hub_auc = int(np.sum(null_aucs >= hub_auc))

robustness_results = {
    "validation_cohort": "GSE62928 (Affymetrix HG-U133 Plus 2.0 / GPL13158)",
    "n_samples": 8,
    "n_cases_eps": 4,
    "n_controls": 4,
    "mapped_genes_universe": len(expr_gene),
    "hub_genes_evaluated": hubs_7,
    "per_gene_concordance": {
        "n_concordant": n_concordant,
        "n_total": n_total,
        "concordance_rate": round(n_concordant / n_total, 4),
        "binomial_sign_test_p": round(binom_p, 4)
    },
    "composite_hub_signature": {
        "hub_signature_auc": round(hub_auc, 4),
        "mann_whitney_u_statistic": float(u_stat),
        "mann_whitney_one_sided_p": round(float(mw_p), 4),
        "mann_whitney_two_sided_p": round(float(mw_two_p), 4)
    },
    "random_gene_set_null_1000": {
        "n_null_sets_at_least_hub_auc": n_null_at_least_hub_auc,
        "null_auc_mean": round(float(np.mean(null_aucs)), 4),
        "null_auc_sd": round(float(np.std(null_aucs)), 4),
        "null_auc_median": round(float(np.median(null_aucs)), 4),
        "empirical_p_value": round(emp_p, 4)
    },
    "scientific_interpretation": (
        f"In GSE62928 biopsy data, {n_concordant}/{n_total} hub genes have concordant direction "
        f"(binomial P = {binom_p:.4f}). The fixed composite score has AUC = {hub_auc:.4f} "
        f"(one-sided Mann-Whitney P = {mw_p:.4f}; two-sided P = {mw_two_p:.4f}); "
        f"the 1,000 random 7-gene-set null has mean AUC = {np.mean(null_aucs):.4f} "
        f"and empirical P = {emp_p:.4f}. The discovery effluent and validation biopsy cohorts "
        "differ in tissue and phenotype, limiting transportability. Historical synthetic "
        "validation outputs from the legacy pipeline are withdrawn."
    )
}

with open(f'results/revision/gse62928_robustness_metrics{OUTPUT_SUFFIX}.json', 'w') as f:
    json.dump(robustness_results, f, indent=2)

# Generate Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#F8FAFC')

# Panel A: Forest plot of log2FC
y_pos = np.arange(len(df_per_gene))
axes[0].errorbar(df_per_gene['Validation_log2FC'], y_pos,
                xerr=[df_per_gene['Validation_log2FC'] - df_per_gene['Validation_CI95_Lower'],
                      df_per_gene['Validation_CI95_Upper'] - df_per_gene['Validation_log2FC']],
                fmt='o', color='#2563EB', ecolor='#64748B', elinewidth=2, capsize=4, markersize=7)
axes[0].axvline(0, color='#DC2626', linestyle='--', linewidth=1.2, label='No effect (log2FC = 0)')
axes[0].set_yticks(y_pos)
axes[0].set_yticklabels(df_per_gene['Gene'], fontsize=10, fontweight='bold')
axes[0].set_xlabel('Validation log2FC (EPS vs Control, 95% CI)', fontsize=10, fontweight='bold')
axes[0].set_title(f'A. Per-Gene Effect Size in GSE62928 (N=8 Biopsies)\nSign Concordance: {n_concordant}/7 ({n_concordant/n_total*100:.1f}%), Binomial P = {binom_p:.3f}',
                  fontsize=11, fontweight='bold')
axes[0].grid(axis='x', linestyle=':', alpha=0.5)
axes[0].legend(loc='lower right', fontsize=9)

# Panel B: Random 7-gene null distribution
axes[1].hist(null_aucs, bins=25, color='#94A3B8', edgecolor='#334155', alpha=0.8, density=True, label='1,000 Random 7-Gene Sets')
axes[1].axvline(hub_auc, color='#DC2626', linewidth=2.5, linestyle='-', label=f'7-Hub Signature AUC = {hub_auc:.4f}')
axes[1].set_xlabel('ROC AUC (4 EPS vs 4 Controls)', fontsize=10, fontweight='bold')
axes[1].set_ylabel('Density', fontsize=10, fontweight='bold')
axes[1].set_title(f'B. Random Gene-Set Null Distribution (1,000 Sets)\nEmpirical P = {emp_p:.4f} (Null Mean AUC = {np.mean(null_aucs):.3f})',
                  fontsize=11, fontweight='bold')
axes[1].legend(loc='upper right', fontsize=9)
axes[1].grid(linestyle=':', alpha=0.5)

plt.tight_layout()
plt.savefig(f'results/revision/Fig_rev_validation_robustness{OUTPUT_SUFFIX}.png', dpi=300)
plt.savefig(f'results/revision/Fig_rev_validation_robustness{OUTPUT_SUFFIX}.pdf')
plt.close()

# Separate ROC output for use by the legacy figure manifest. The score is a
# fixed, per-gene z-score average; no classifier is fitted to the validation set.
fpr, tpr, _ = roc_curve(y_true, hub_score)
fig, ax = plt.subplots(figsize=(7, 6), facecolor='#F8FAFC')
ax.plot(fpr, tpr, color='#2563EB', linewidth=2.5,
        label=f'Fixed 7-gene score (AUC = {hub_auc:.4f})')
ax.plot([0, 1], [0, 1], color='#64748B', linestyle=':', label='Chance reference')
ax.set(
    xlim=(0, 1), ylim=(0, 1),
    xlabel='False positive rate',
    ylabel='True positive rate',
    title='GSE62928 Biopsy Cohort: Fixed Hub-Score ROC\n'
          f'4 EPS cases vs 4 controls; two-sided Mann-Whitney P = {mw_two_p:.4f}'
)
ax.legend(loc='lower right')
ax.grid(linestyle=':', alpha=0.5)
fig.tight_layout()
fig.savefig(f'results/revision/Fig_rev_validation_roc{OUTPUT_SUFFIX}.png', dpi=300)
fig.savefig(f'results/revision/Fig_rev_validation_roc{OUTPUT_SUFFIX}.pdf')
plt.close(fig)

print("Task 3 Validation Robustness completed and saved successfully!")
