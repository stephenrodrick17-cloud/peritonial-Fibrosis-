import pandas as pd

print("=== 1. Total Probed Genes (20,940) ===")
df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
print(f"File: results/tables/GSE62928_full_expression_matrix.csv | Index length = {len(df_expr)}")

print("\n=== 2. Pro-Fibrotic Up-Regulated DEGs (534) ===")
df_tt = pd.read_csv("GSE62928.top.table.tsv", sep="\t")
df_clean = df_tt.dropna(subset=["Gene.symbol"]).copy()
df_clean = df_clean[~df_clean["Gene.symbol"].isin(["", "---"])]
rows = []
for idx, r in df_clean.iterrows():
    for sym in str(r["Gene.symbol"]).split("///"):
        s = sym.strip().upper()
        if s:
            rc = r.to_dict()
            rc["Gene_clean"] = s
            rows.append(rc)
df_exp = pd.DataFrame(rows)
df_best_p = df_exp.sort_values("P.Value").drop_duplicates("Gene_clean")
degs_534 = df_best_p[(df_best_p["P.Value"] < 0.05) & (df_best_p["logFC"] >= 0.585)]
print(f"File: GSE62928.top.table.tsv | Unique Gene.symbol with P.Value < 0.05 and logFC >= 0.585 [1.5-fold] = {len(degs_534)}")

print("\n=== 3. Curated Matrisome Reference (1,027) ===")
df_mat = pd.read_excel("ECM genes all.xlsx", skiprows=1)
sym_col = [c for c in df_mat.columns if "symbol" in c.lower() or "gene" in c.lower()][0]
mat_genes = set(df_mat[sym_col].dropna().astype(str).str.strip().str.upper()) - {"GENE SYMBOL", "NA", "NAN"}
print(f"File: ECM genes all.xlsx | Unique genes in col '{sym_col}' = {len(mat_genes)}")

print("\n=== 4. ECM-DEGs (81) ===")
df_71 = pd.read_csv("convergent_81_ECM_DEGs.csv")
print(f"File: convergent_81_ECM_DEGs.csv | Rows = {len(df_71)} | Columns = {list(df_71.columns)}")

print("\n=== 5. WGCNA Salmon Module Genes (604) ===")
df_wgcna = pd.read_csv("results/tables/wgcna_trait_significant_module_genes.csv")
print(f"File: results/tables/wgcna_trait_significant_module_genes.csv | Rows = {len(df_wgcna)} | Columns = {list(df_wgcna.columns)}")

print("\n=== 6. Convergent WGCNA-ECM Candidates (44) ===")
df_conv = pd.read_csv("results/tables/convergent_WGCNA_ECM_genes.csv")
print(f"File: results/tables/convergent_WGCNA_ECM_genes.csv | Rows = {len(df_conv)} | Columns = {list(df_conv.columns)}")

print("\n=== 7. 11 ML Consensus Hub Genes & Votes ===")
df_hubs = pd.read_csv("results/tables/ML_hub_genes_from_WGCNA_ECM.csv")
print(f"File: results/tables/ML_hub_genes_from_WGCNA_ECM.csv | Rows = {len(df_hubs)}")
for _, r in df_hubs.iterrows():
    print(f"  Gene: {r['Gene_Symbol']:10} | Votes: {r['Votes']} | Models: {r['Selecting_Models']}")

print("\n=== 8. All Discrimination Metrics & AUCs from results/tables/roc_auc_detailed_metrics.csv ===")
df_auc = pd.read_csv("results/tables/roc_auc_detailed_metrics.csv")
print(f"Columns: {list(df_auc.columns)}")
for _, r in df_auc.iterrows():
    ci_low = r.get('CI_95_Low', r.get('CI_95_low', 'N/A'))
    ci_high = r.get('CI_95_High', r.get('CI_95_high', 'N/A'))
    auc_val = r.get('AUC', 'N/A')
    ci_str = f"[{ci_low:.3f}, {ci_high:.3f}]" if (pd.notna(ci_low) and isinstance(ci_low, (int, float))) else "N/A"
    auc_str = f"{auc_val:.4f}" if (pd.notna(auc_val) and isinstance(auc_val, (int, float))) else "N/A"
    print(f"  Model: {str(r.get('Model_Type'))[:38]:38} | Feature: {str(r.get('Feature_or_Panel'))[:35]:35} | Status: {str(r.get('Platform_Status'))[:22]:22} | AUC: {auc_str} | 95% CI: {ci_str}")
