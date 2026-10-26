import os
import json
from pathlib import Path
import pandas as pd
from scipy import stats

os.makedirs("results/revision", exist_ok=True)
os.makedirs("scripts", exist_ok=True)

report = {}

# -------------------------------------------------------------------------
# Item 1: Legacy validation source inspection (no synthetic rerun)
# -------------------------------------------------------------------------
pipeline_source = Path("run_pipeline.py").read_text(encoding="utf-8")
synthetic_validation_source_present = (
    "np.random.RandomState(RANDOM_STATE).randn(8, len(mapped_hubs))" in pipeline_source
    and "X_val_raw[idx_c, :] += 1.5" in pipeline_source
)

report["item1_permutation_p_audit"] = {
    "synthetic_validation_source_present": synthetic_validation_source_present,
    "analysis_method": "Static inspection of historical source; no synthetic data was generated",
    "verdict": "INVALID HISTORICAL VALIDATION IMPLEMENTATION",
    "explanation": (
        "The historical validation phase generated Gaussian features and shifted the first four samples before "
        "calculating apparent, cross-validation, and permutation metrics. These computations were not based on "
        "GSE62928 expression. The historical permutation output is therefore withdrawn; the real-cohort analysis "
        "is reported by scripts/task3_validation_robustness.py."
    )
}

# -------------------------------------------------------------------------
# Item 2: Verify 7 consensus hubs from 25-gene WGCNA pool vs 31-gene Up-ECM
# -------------------------------------------------------------------------
ml_3 = pd.read_csv("results/final_analysis/ML_3_models_hub_genes.csv")
ml_31 = pd.read_csv("results/final_analysis/ML_3_models_hub_genes_all31_up_ecm.csv")
convergence = pd.read_csv("results/tables/convergent_WGCNA_ECM_genes.csv")

hubs_25_strict = ml_3[ml_3["Is_Strict_Hub_3of3"] == 1]["Gene"].tolist()
hubs_25_majority = ml_3[ml_3["Is_Consensus_Hub_2of3"] == 1]["Gene"].tolist()

hubs_31_strict = ml_31[ml_31["Is_Strict_Hub_3of3"] == 1]["Gene"].tolist()
hubs_31_majority = ml_31[ml_31["Is_Consensus_Hub_2of3"] == 1]["Gene"].tolist()
convergent_genes = set(convergence.loc[convergence["Convergent"] == 1, "Gene"])
hub_convergence_flags = convergence[convergence["Gene"].isin(hubs_25_majority)][
    ["Gene", "in_Top_WGCNA_Module", "Convergent"]
].to_dict(orient="records")

report["item2_hub_derivation_audit"] = {
    "three_model_feature_rows_n": int(len(ml_3)),
    "three_model_candidate_genes": ml_3["Gene"].astype(str).tolist(),
    "primary_strict_3of3_hubs": hubs_25_strict,
    "primary_majority_2of3_hubs": hubs_25_majority,
    "alternate_up_ecm_pool_n": 31,
    "alternate_strict_3of3_hubs": hubs_31_strict,
    "alternate_majority_2of3_hubs": hubs_31_majority,
    "primary_hub_convergence_flags": hub_convergence_flags,
    "primary_hubs_in_convergent_pool": sorted(set(hubs_25_majority) & convergent_genes),
    "difference_commentary": (
        "The three-model feature table contains 25 rows and marks 7 genes as majority hubs. However, the separate "
        "convergent_WGCNA_ECM_genes.csv table marks all seven as in_Top_WGCNA_Module=0 and Convergent=0. "
        "Therefore, the available files do not verify that these hubs were selected from a 25-gene WGCNA-convergent pool. "
        "The alternate 31-gene Up-ECM analysis is a distinct analysis and yields 8 majority hubs "
        "(TNFSF15, EBI3, FLT3LG, ADAM19, CLEC4F, CXCL14, PI3, TGM3)."
    )
}

# -------------------------------------------------------------------------
# Item 3: MEblack vs MEmagenta Documentation
# -------------------------------------------------------------------------
mt = pd.read_csv("results/tables/wgcna_module_trait_correlation.csv")
mm = pd.read_csv("results/tables/wgcna_gene_module_membership.csv")
up_ecm = pd.read_csv("results/final_analysis/ECM_Up_DEGs_profibrotic.csv")["Gene"].tolist()
all_results = pd.read_csv("results/tables/GSE125498_all_results.csv")
measured_ecm = pd.read_csv("results/tables/GSE125498_ECM_intersection.csv")
expression = pd.read_csv("results/tables/GSE125498_full_expression_matrix.csv").set_index("Gene")
metadata = pd.read_csv("results/tables/GSE125498_sample_metadata.csv").set_index("sample_id")
eigengenes = pd.read_csv("results/tables/wgcna_module_eigengenes.csv").set_index("sample_id")
sample_ids = [sample for sample in metadata.index if sample in expression.columns and sample in eigengenes.index]
stage = metadata.loc[sample_ids, "stage_binary"].astype(float).to_numpy()

up_degs = set(all_results.loc[all_results["DEG_Status"] == "Upregulated", "Gene"].astype(str))
measured_ecm_genes = set(measured_ecm["Gene"].astype(str))
positive_module_colors = set(mt.loc[mt["Correlation"] > 0, "Module_Color"].astype(str))
positive_module_genes = set(mm.loc[mm["Module"].isin(positive_module_colors), "Gene"].astype(str))
fallback_candidate_pool = sorted(up_degs & measured_ecm_genes & positive_module_genes)
top_trait_module_pool = sorted(
    up_degs
    & measured_ecm_genes
    & set(mm.loc[mm["Module"] == mt.iloc[0]["Module_Color"], "Gene"].astype(str))
)
candidate_table_genes = set(ml_3["Gene"].astype(str))
report["item2_hub_derivation_audit"].update({
    "top_trait_module": str(mt.iloc[0]["Module"]),
    "top_trait_module_candidate_n": len(top_trait_module_pool),
    "positive_module_fallback_candidate_n": len(fallback_candidate_pool),
    "positive_module_fallback_matches_25_row_feature_table": (
        set(fallback_candidate_pool) == candidate_table_genes
    ),
    "positive_module_fallback_genes": fallback_candidate_pool,
    "positive_trait_module_colors": sorted(positive_module_colors)
})

ecm_in_magenta = sorted(list(set(up_ecm) & set(mm[mm["Module"] == "magenta"]["Gene"])))
ecm_in_black = sorted(list(set(up_ecm) & set(mm[mm["Module"] == "black"]["Gene"])))
hub_module_metrics = []
for gene in hubs_25_majority:
    values = expression.loc[gene, sample_ids].astype(float).to_numpy()
    kme = stats.pearsonr(values, eigengenes.loc[sample_ids, "black"].astype(float).to_numpy())
    gene_significance = stats.pearsonr(values, stage)
    hub_module_metrics.append({
        "Gene": gene,
        "Assigned_Module": mm.loc[mm["Gene"] == gene, "Module"].iloc[0],
        "kME_black": float(kme.statistic),
        "kME_black_p": float(kme.pvalue),
        "GS_stage": float(gene_significance.statistic),
        "GS_stage_p": float(gene_significance.pvalue)
    })

magenta_row = mt.loc[mt["Module"] == "MEmagenta"].iloc[0]
black_row = mt.loc[mt["Module"] == "MEblack"].iloc[0]
magenta_gene_n = int((mm["Module"] == "magenta").sum())
black_gene_n = int((mm["Module"] == "black").sum())
black_hub_kme = [row["kME_black"] for row in hub_module_metrics]
report["item3_module_justification"] = {
    "MEmagenta": {
        "correlation_with_stage": float(mt[mt["Module"] == "MEmagenta"]["Correlation"].iloc[0]),
        "p_value": float(mt[mt["Module"] == "MEmagenta"]["P_value"].iloc[0]),
        "total_genes": int((mm["Module"] == "magenta").sum()),
        "up_ecm_genes_contained": len(ecm_in_magenta),
        "ecm_gene_list": ecm_in_magenta
    },
    "MEblack": {
        "correlation_with_stage": float(mt[mt["Module"] == "MEblack"]["Correlation"].iloc[0]),
        "p_value": float(mt[mt["Module"] == "MEblack"]["P_value"].iloc[0]),
        "total_genes": int((mm["Module"] == "black").sum()),
        "up_ecm_genes_contained": len(ecm_in_black),
        "ecm_gene_list": ecm_in_black
    },
    "hub_module_membership_and_gene_significance": hub_module_metrics,
    "scientific_justification": (
        f"MEmagenta is the strongest module-level stage association (r = {magenta_row['Correlation']:+.4f}, "
        f"P = {magenta_row['P_value']:.4g}; {magenta_gene_n} genes) and contains "
        f"{len(ecm_in_magenta)} upregulated ECM DEGs. MEblack contains {black_gene_n} genes and "
        f"{len(ecm_in_black)} upregulated ECM DEGs, but its eigengene is not associated with stage "
        f"(r = {black_row['Correlation']:+.4f}, P = {black_row['P_value']:.4g}). All seven majority hubs are "
        f"assigned to black, yet their calculated kME_black correlations range from {min(black_hub_kme):+.4f} "
        f"to {max(black_hub_kme):+.4f} and are nonsignificant. The files support a distinct MEmagenta "
        "module-trait association and individual hub-stage associations, but do not establish black-module "
        "co-expression or top-module convergence of the hubs."
    )
}

# -------------------------------------------------------------------------
# Item 4: Illumina Normalization Documentation
# -------------------------------------------------------------------------
illum_norm_info = (
    "According to the official GEO series matrix header (!Sample_data_processing) for GSE125498, "
    "bead-level raw data were summarized using Illumina BeadStudio Software v3, imported into R, and normalized "
    "using the quantile normalization method via the Lumi package (Du et al., Bioinformatics 2008). "
    "Probes were filtered to retain only those with detectable signal intensity in at least five samples."
)

report["item4_illumina_normalization"] = {
    "software": "Illumina BeadStudio v3",
    "package": "R lumi",
    "algorithm": "Quantile normalization",
    "filtering_rule": "Detectable signal intensity in >= 5 samples",
    "documentation": illum_norm_info
}

# Save a new audit version without replacing prior revision outputs.
with open("results/revision/task1_audit_report_taskA.json", "w") as f:
    json.dump(report, f, indent=2)

with open("results/revision/task1_audit_report_taskA.md", "w") as f:
    f.write("# Task 1 Audit Report: Verification and Discrepancy Analysis\n\n")
    f.write("## 1. Exact Permutation P in GSE62928 Validation\n")
    f.write(report["item1_permutation_p_audit"]["explanation"] + "\n\n")
    f.write(f"- Legacy synthetic source pattern remains in pipeline: {synthetic_validation_source_present}\n\n")
    f.write("## 2. ML Candidate-Pool and Hub Derivation\n")
    f.write(report["item2_hub_derivation_audit"]["difference_commentary"] + "\n\n")
    f.write(f"- Three-model feature rows: {len(ml_3)}; strict hubs ({hubs_25_strict}); majority hubs ({hubs_25_majority})\n")
    f.write(
        f"- Top-trait-module up-DEG/matrisome intersection: "
        f"{report['item2_hub_derivation_audit']['top_trait_module_candidate_n']}\n"
    )
    f.write(
        f"- Positive-module fallback pool: "
        f"{report['item2_hub_derivation_audit']['positive_module_fallback_candidate_n']}; "
        f"matches three-model candidate table: "
        f"{report['item2_hub_derivation_audit']['positive_module_fallback_matches_25_row_feature_table']}\n"
    )
    f.write(f"- Hubs marked convergent in the WGCNA table: {report['item2_hub_derivation_audit']['primary_hubs_in_convergent_pool']}\n")
    f.write(f"- Alternate Up-ECM pool (N=31): 2 strict hubs ({hubs_31_strict}) and 8 majority hubs ({hubs_31_majority})\n\n")
    f.write("## 3. Module Justification: MEblack vs MEmagenta\n")
    f.write(report["item3_module_justification"]["scientific_justification"] + "\n\n")
    f.write("| Hub | kME to MEblack | kME P | Stage gene significance (r) | GS P |\n|---|---:|---:|---:|---:|\n")
    for row in hub_module_metrics:
        f.write(
            f"| {row['Gene']} | {row['kME_black']:.4f} | {row['kME_black_p']:.4g} | "
            f"{row['GS_stage']:.4f} | {row['GS_stage_p']:.4g} |\n"
        )
    f.write("\n")
    f.write("## 4. Illumina Normalization Protocol\n")
    f.write(report["item4_illumina_normalization"]["documentation"] + "\n")

print("Task 1 Audit completed and saved to new taskA report files under results/revision.")
