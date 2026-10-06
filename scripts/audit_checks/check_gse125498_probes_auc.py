# scripts/audit_checks/check_gse125498_probes_auc.py
import gzip
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score

# 1. Load GPL10558 probe mappings
annot_file = "data/raw/GPL10558.annot.gz"
probe_gene_map = {}
with gzip.open(annot_file, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!platform_table_begin"):
            break
    header = f.readline().strip().split("\t")
    id_idx = header.index("ID")
    symbol_idx = header.index("Gene symbol")
    for line in f:
        if line.startswith("!platform_table_end"):
            break
        parts = line.strip().split("\t")
        if len(parts) > symbol_idx:
            pid = parts[id_idx]
            sym = parts[symbol_idx]
            if sym:
                # Handle multiple symbols separated by ///
                syms = [s.strip() for s in sym.split("///")]
                for s in syms:
                    probe_gene_map.setdefault(s, []).append(pid)

# Hub genes
hub_genes = ["FN1", "COL3A1", "COL11A1", "COL8A1", "VCAN", "COMP", "THBS3", "EDIL3", "LOX", "INHBA", "ISM1"]

# 2. Load GSE125498 series matrix
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

# Load metadata to get group labels (SPD vs LPD)
df_meta = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
# Map sample to group
# stage_binary: 1 = LPD (13), 0 = SPD (20)
meta_map = dict(zip(df_meta["sample_id"], df_meta["stage_binary"]))
y_true = np.array([meta_map[s] for s in df_expr.columns])

print("Total samples:", len(y_true), "LPD:", sum(y_true == 1), "SPD:", sum(y_true == 0))

# 3. Check each hub gene's candidate probes
audit_rows = []
for gene in hub_genes:
    cand_probes = probe_gene_map.get(gene, [])
    cand_in_matrix = [p for p in cand_probes if p in df_expr.index]
    if not cand_in_matrix:
        audit_rows.append({
            "Gene": gene,
            "Probe_ID": "NO PROBE",
            "Mean_Expr": np.nan,
            "log2FC": np.nan,
            "t_stat": np.nan,
            "P_Welch": np.nan,
            "MW_P": np.nan,
            "AUC": np.nan,
            "Role": "Missing on GPL10558"
        })
        continue
    
    for pid in cand_in_matrix:
        vals = df_expr.loc[pid].values
        vals_lpd = vals[y_true == 1]
        vals_spd = vals[y_true == 0]
        mean_all = vals.mean()
        lfc = vals_lpd.mean() - vals_spd.mean()
        tt = stats.ttest_ind(vals_lpd, vals_spd, equal_var=False)
        mw = stats.mannwhitneyu(vals_lpd, vals_spd, alternative='two-sided')
        auc = roc_auc_score(y_true, vals)
        audit_rows.append({
            "Gene": gene,
            "Probe_ID": pid,
            "Mean_Expr": mean_all,
            "log2FC": lfc,
            "t_stat": tt.statistic,
            "P_Welch": tt.pvalue,
            "MW_P": mw.pvalue,
            "AUC": auc,
            "Role": ""
        })

df_res = pd.DataFrame(audit_rows)
print(df_res.to_string(index=False))

# Now check which probes were actually used in README and tables:
print("\n--- Compare highest-mean vs lowest-P probe for genes with multiple probes ---")
for gene in ["FN1", "COL8A1"]:
    sub = df_res[df_res["Gene"] == gene].copy()
    max_mean_row = sub.sort_values("Mean_Expr", ascending=False).iloc[0]
    lowest_p_row = sub.sort_values("P_Welch", ascending=True).iloc[0]
    print(f"\nGene: {gene}")
    print(f"  Highest-mean probe: {max_mean_row['Probe_ID']} (Mean={max_mean_row['Mean_Expr']:.3f}, log2FC={max_mean_row['log2FC']:+.4f}, P={max_mean_row['P_Welch']:.4f}, AUC={max_mean_row['AUC']:.4f})")
    print(f"  Lowest-P probe:     {lowest_p_row['Probe_ID']} (Mean={lowest_p_row['Mean_Expr']:.3f}, log2FC={lowest_p_row['log2FC']:+.4f}, P={lowest_p_row['P_Welch']:.4f}, AUC={lowest_p_row['AUC']:.4f})")
