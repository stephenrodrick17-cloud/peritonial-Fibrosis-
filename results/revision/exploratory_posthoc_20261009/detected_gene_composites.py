#!/usr/bin/env python
"""Post hoc / exploratory detected-gene scores; never replace fixed-panel tests."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    helper_dir = (
        root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE248762"
        / "scripts"
    )
    sys.path.insert(0, str(helper_dir))
    from gse248762_step2_4_scanpy import (
        N_NULL_SETS,
        calculate_set_summary,
        empirical_upper_p,
        make_null_sets,
    )

    base = root / "results" / "revision" / "effluent_validation"
    summaries: list[dict] = []
    null_rows: list[pd.DataFrame] = []
    sample_score_rows: list[pd.DataFrame] = []

    for accession in ("GSE248762", "GSE130888"):
        analysis = base / accession / "analysis_20261009"
        primary = pd.read_csv(analysis / "edger_all_gene_results.csv")
        primary = primary[primary["contrast"] == "LV_vs_SV"].copy()
        if primary.empty:
            raise ValueError(f"{accession}: no saved LV_vs_SV edgeR contrast")
        detected_panel = pd.read_csv(analysis / "primary_per_gene_results.csv")
        detected_panel = detected_panel[detected_panel["detected"].astype(bool)]
        if len(detected_panel) != (7 if accession == "GSE248762" else 6):
            raise ValueError(
                f"{accession}: detected panel size differs from saved result"
            )
        gene_ids = detected_panel["gene_id"].astype(str).tolist()
        genes = detected_panel["gene"].astype(str).tolist()
        if not set(gene_ids).issubset(set(primary["gene_id"].astype(str))):
            raise ValueError(f"{accession}: detected gene absent from edgeR fit")

        metadata = pd.read_csv(
            analysis / "all_cell_pseudobulk_metadata.csv"
        )[["sample_id", "primary_group"]]
        auc_metadata = metadata.loc[
            metadata["primary_group"].isin(["LV", "SV"])
        ].reset_index(drop=True)
        expected_group_sizes = (
            {"LV": 10, "SV": 6}
            if accession == "GSE248762"
            else {"LV": 4, "SV": 6}
        )
        observed_group_sizes = auc_metadata["primary_group"].value_counts().to_dict()
        if observed_group_sizes != expected_group_sizes:
            raise ValueError(
                f"{accession}: expected LV/SV group sizes "
                f"{expected_group_sizes}, found {observed_group_sizes}"
            )
        logcpm = pd.read_csv(analysis / "edger_logCPM.csv", index_col=0)
        detection = pd.read_csv(analysis / "pseudobulk_gene_detection.csv")
        if "detected" not in detection or "gene_id" not in detection:
            raise ValueError(f"{accession}: unexpected detection table schema")
        detection_ids = detection.loc[
            detection["detected"].astype(bool), "gene_id"
        ].astype(str)
        fit_ids = set(primary["gene_id"].astype(str))
        eligible_ids = [
            gene_id
            for gene_id in detection_ids
            if gene_id in fit_ids and gene_id in logcpm.index
        ]
        eligible = pd.DataFrame(
            {
                "gene_id": eligible_ids,
                "mean_logCPM": [
                    float(logcpm.loc[gene_id].mean())
                    for gene_id in eligible_ids
                ],
            }
        )
        observed = calculate_set_summary(
            gene_ids,
            primary,
            logcpm,
            auc_metadata,
            seed=42,
            bootstrap=True,
        )
        if not observed.get("estimable", False):
            raise ValueError(f"{accession}: detected-gene composite not estimable")

        null_sets, null_details = make_null_sets(
            eligible,
            gene_ids,
            N_NULL_SETS,
            seed=42,
        )
        if len(null_sets) != N_NULL_SETS:
            raise ValueError(
                f"{accession}: expected {N_NULL_SETS} matched null sets; "
                f"generated {len(null_sets)}"
            )
        primary_by_id = primary.set_index("gene_id")
        local_null_rows = []
        for index, null_set in enumerate(null_sets, start=1):
            null_summary = calculate_set_summary(
                null_set,
                primary,
                logcpm,
                auc_metadata,
                seed=42,
                bootstrap=False,
            )
            local_null_rows.append(
                {
                    "accession": accession,
                    "set_number": index,
                    "matched_null_mean_gene_log2FC": float(
                        primary_by_id.loc[null_set, "logFC"].mean()
                    ),
                    "matched_null_composite_auc": float(
                        null_summary["auc"]
                    ),
                }
            )
        null_table = pd.DataFrame(local_null_rows)
        null_rows.append(null_table)
        auc_null_p = empirical_upper_p(
            null_table["matched_null_composite_auc"].to_numpy(),
            float(observed["auc"]),
        )
        mean_lfc_null_p = empirical_upper_p(
            null_table["matched_null_mean_gene_log2FC"].to_numpy(),
            float(observed["mean_gene_log2FC"]),
        )
        expression = logcpm.loc[
            gene_ids, metadata["sample_id"]
        ].to_numpy(dtype=float)
        auc_columns = metadata["sample_id"].isin(
            auc_metadata["sample_id"]
        ).to_numpy()
        expression_reference = expression[:, auc_columns]
        z_scores = (
            expression - expression_reference.mean(axis=1, keepdims=True)
        ) / expression_reference.std(axis=1, ddof=0, keepdims=True)
        composite = z_scores.mean(axis=0)
        score_table = metadata.copy()
        score_table["detected_gene_composite"] = composite
        score_table["auc_included"] = score_table["primary_group"].isin(
            ["LV", "SV"]
        )
        score_table["analysis_label"] = (
            "post hoc detected-gene composite, not the preregistered test"
        )
        score_table["accession"] = accession
        sample_score_rows.append(score_table)

        summaries.append(
            {
                "analysis_label": (
                    "post hoc detected-gene composite, not the preregistered test"
                ),
                "accession": accession,
                "primary_contrast": (
                    "LV (LV_NOT_UF + LV_UF) vs SV"
                    if accession == "GSE248762"
                    else "Long-term PD vs short-term PD"
                ),
                "detected_genes_in_score": genes,
                "n_genes": len(genes),
                "n_LV": int((auc_metadata["primary_group"] == "LV").sum()),
                "n_SV": int((auc_metadata["primary_group"] == "SV").sum()),
                "excluded_noncontrast_group_sizes": {
                    str(group): int(count)
                    for group, count in metadata.loc[
                        ~metadata["primary_group"].isin(["LV", "SV"]),
                        "primary_group",
                    ].value_counts().items()
                },
                "excluded_noncontrast_sample_ids": metadata.loc[
                    ~metadata["primary_group"].isin(["LV", "SV"]),
                    "sample_id",
                ].astype(str).tolist(),
                "composite_auc": float(observed["auc"]),
                "auc_bootstrap_ci_95": [
                    float(observed["auc_ci_lower"]),
                    float(observed["auc_ci_upper"]),
                ],
                "mann_whitney_p": float(observed["mann_whitney_p"]),
                "mean_score_LV": float(observed["mean_score_LV"]),
                "mean_score_SV": float(observed["mean_score_SV"]),
                "mean_detected_gene_log2FC": float(
                    observed["mean_gene_log2FC"]
                ),
                "matched_null_sets": len(null_sets),
                "matched_null_auc_empirical_p": auc_null_p,
                "matched_null_mean_log2FC_empirical_p": mean_lfc_null_p,
                "null_draw_details": null_details,
                "seed": 42,
                "preregistered_verdict_changed": False,
            }
        )

    (output / "detected_gene_composites_summary.json").write_text(
        json.dumps(summaries, indent=2), encoding="utf-8"
    )
    pd.DataFrame(summaries).to_csv(
        output / "detected_gene_composites_summary.csv", index=False
    )
    pd.concat(null_rows, ignore_index=True).to_csv(
        output / "detected_gene_composites_matched_null.csv", index=False
    )
    pd.concat(sample_score_rows, ignore_index=True).to_csv(
        output / "detected_gene_composites_patient_scores.csv", index=False
    )
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
