import os
import json
import gzip
import pandas as pd
import numpy as np

os.makedirs("results/revision", exist_ok=True)
os.makedirs("scripts", exist_ok=True)

# 1. Load tables
all_res = pd.read_csv('results/tables/GSE125498_all_results.csv').set_index('Gene')
degs_813 = pd.read_csv('results/tables/GSE125498_DEGs_filtered.csv')['Gene'].tolist()
ecm_44 = pd.read_csv('results/final_analysis/FINAL_DEG_Matrisomal_Genes.csv')['Gene'].tolist()
hubs_7 = ['FLT3LG', 'TNFSF15', 'LTB', 'ADAM19', 'SERPINA10', 'EBI3', 'CXCL14']

# 2. FDR Sensitivity Breakdown
def get_fdr_counts(gene_list, label):
    sub = all_res.loc[gene_list]
    n_tot = len(gene_list)
    n_fdr05 = int((sub['adj.P.Val'] < 0.05).sum())
    n_fdr10 = int((sub['adj.P.Val'] < 0.10).sum())
    return {
        "Gene_Set": label,
        "Total_Genes": n_tot,
        "Surviving_FDR_0_05": n_fdr05,
        "Pct_Surviving_FDR_0_05": round(100.0 * n_fdr05 / n_tot, 2),
        "Surviving_FDR_0_10": n_fdr10,
        "Pct_Surviving_FDR_0_10": round(100.0 * n_fdr10 / n_tot, 2),
        "FDR_Non_Significant_p10": n_tot - n_fdr10,
        "Pct_Non_Significant_p10": round(100.0 * (n_tot - n_fdr10) / n_tot, 2)
    }

fdr_summary_rows = [
    get_fdr_counts(degs_813, "All Primary DEGs (N=813)"),
    get_fdr_counts(ecm_44, "DEG-Matrisomal Genes (N=44)"),
    get_fdr_counts(hubs_7, "Consensus Hub Genes (N=7)")
]
df_fdr = pd.DataFrame(fdr_summary_rows)
df_fdr.to_csv("results/revision/fdr_sensitivity_table.csv", index=False)

# Detail on the 7 hubs
hub_detail_rows = []
for g in hubs_7:
    r = all_res.loc[g]
    hub_detail_rows.append({
        "Gene": g,
        "log2FC": round(float(r["logFC"]), 4),
        "Nominal_P": float(r["P.Value"]),
        "BH_adj_P_FDR": float(r["adj.P.Val"]),
        "Survives_FDR_0_05": "Yes" if r["adj.P.Val"] < 0.05 else "No",
        "Survives_FDR_0_10": "Yes" if r["adj.P.Val"] < 0.10 else "No"
    })
df_hub_fdr = pd.DataFrame(hub_detail_rows)
df_hub_fdr.to_csv("results/revision/hubs_fdr_breakdown.csv", index=False)

# 3. GEO Metadata Audit for Age, Sex, Batch
matrix_gz = "data/raw/GSE125498_series_matrix.txt.gz"
char_lines = []
with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_description"):
            char_lines.append(line.strip())

meta_audit = {
    "geo_accession": "GSE125498",
    "metadata_lines_inspected": len(char_lines),
    "sample_characteristics_found": [
        "patient id (1 to 33)",
        "time (peritoneal dialysis): short-term vs long-term",
        "cell type: peritoneal cells"
    ],
    "age_available": False,
    "sex_available": False,
    "batch_id_available": False,
    "statement": (
        "The public GEO deposit for GSE125498 contains only 'patient id', 'time (peritoneal dialysis): short-term / long-term', "
        "and 'cell type: peritoneal cells'. Clinical covariates including patient age, sex, race, diabetes status, residual renal function, "
        "and technical microarray processing batch dates are not provided in the GEO series matrix. "
        "Consequently, multivariable adjustment in limma for age, sex, or batch cannot be performed. "
        "This absence is explicitly documented as a potential source of unmeasured confounding in the study limitations."
    )
}

with open("results/revision/task5_fdr_summary.json", "w") as f:
    json.dump({
        "fdr_breakdown": fdr_summary_rows,
        "hubs_fdr": hub_detail_rows,
        "metadata_covariate_audit": meta_audit
    }, f, indent=2)

print("Task 5 completed and saved successfully!")
