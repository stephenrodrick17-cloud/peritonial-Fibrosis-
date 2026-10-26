import json
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm

from task4_marker_sets import HUB_GENES, LEUKOCYTE_MARKERS, LYMPHOCYTE_MARKERS


OUTPUT_DIR = "results/revision"
os.makedirs(OUTPUT_DIR, exist_ok=True)

metadata = pd.read_csv("results/tables/GSE125498_sample_metadata.csv")
expression = pd.read_csv("results/tables/GSE125498_full_expression_matrix.csv").set_index("Gene")
sample_ids = [sample for sample in metadata["sample_id"] if sample in expression.columns]
metadata = metadata.set_index("sample_id").loc[sample_ids]
stage = metadata["stage_binary"].to_numpy(dtype=float)


def get_present_markers(marker_set):
    return [gene for gene in marker_set if gene in expression.index]


def signature_score(marker_genes):
    if not marker_genes:
        raise ValueError("Cannot calculate a marker score: no marker genes are present.")
    values = expression.loc[marker_genes, sample_ids].to_numpy(dtype=float)
    scale = values.std(axis=1, ddof=1, keepdims=True)
    if np.any(scale <= 0):
        raise ValueError("Marker score includes a gene with zero sample variance.")
    return ((values - values.mean(axis=1, keepdims=True)) / scale).mean(axis=0)


present_leukocyte = get_present_markers(LEUKOCYTE_MARKERS)
present_lymphocyte = get_present_markers(LYMPHOCYTE_MARKERS)
hub_overlap_leukocyte = sorted(set(HUB_GENES) & set(present_leukocyte))
hub_overlap_lymphocyte = sorted(set(HUB_GENES) & set(present_lymphocyte))
filtered_leukocyte = [gene for gene in present_leukocyte if gene not in HUB_GENES]
filtered_lymphocyte = [gene for gene in present_lymphocyte if gene not in HUB_GENES]

score_pairs = {
    "total_leukocyte": (
        signature_score(present_leukocyte),
        signature_score(filtered_leukocyte),
    ),
    "lymphocyte": (
        signature_score(present_lymphocyte),
        signature_score(filtered_lymphocyte),
    ),
}

records = []
for gene in HUB_GENES:
    expression_values = expression.loc[gene, sample_ids].to_numpy(dtype=float)
    record = {
        "Gene": gene,
        "Included_in_total_leukocyte_marker_set": gene in present_leukocyte,
        "Included_in_lymphocyte_marker_set": gene in present_lymphocyte,
    }
    for score_name, (original_score, hub_excluded_score) in score_pairs.items():
        for version, score in (
            ("original", original_score),
            ("hub_excluded", hub_excluded_score),
        ):
            design = sm.add_constant(np.column_stack([stage, score]))
            model = sm.OLS(expression_values, design).fit()
            record[f"{score_name}_{version}_stage_beta"] = float(model.params[1])
            record[f"{score_name}_{version}_stage_p"] = float(model.pvalues[1])
            ci = model.conf_int()[1]
            record[f"{score_name}_{version}_stage_ci95_lower"] = float(ci[0])
            record[f"{score_name}_{version}_stage_ci95_upper"] = float(ci[1])
    records.append(record)

result = pd.DataFrame(records)
result.to_csv(
    os.path.join(OUTPUT_DIR, "leukocyte_hub_overlap_sensitivity.csv"),
    index=False,
)

summary = {
    "analysis": "Marker-score covariate sensitivity with consensus hubs excluded",
    "marker_score_method": (
        "Mean of per-gene z-scores across samples; expression-signature proxy, "
        "not cell-fraction deconvolution"
    ),
    "hub_genes": list(HUB_GENES),
    "present_total_leukocyte_markers": present_leukocyte,
    "present_lymphocyte_markers": present_lymphocyte,
    "hub_overlap_total_leukocyte": hub_overlap_leukocyte,
    "hub_overlap_lymphocyte": hub_overlap_lymphocyte,
    "total_leukocyte_score_unchanged_after_exclusion": bool(
        np.array_equal(score_pairs["total_leukocyte"][0], score_pairs["total_leukocyte"][1])
    ),
    "lymphocyte_score_unchanged_after_exclusion": bool(
        np.array_equal(score_pairs["lymphocyte"][0], score_pairs["lymphocyte"][1])
    ),
    "output_table": "leukocyte_hub_overlap_sensitivity.csv",
}
with open(
    os.path.join(OUTPUT_DIR, "leukocyte_hub_overlap_sensitivity_summary.json"),
    "w",
    encoding="utf-8",
) as output_file:
    json.dump(summary, output_file, indent=2)

print(json.dumps(summary, indent=2))
