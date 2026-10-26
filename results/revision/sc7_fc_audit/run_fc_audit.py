from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SINGLE_CELL = ROOT / "results" / "revision" / "single_cell_context"
NETWORK = ROOT / "results" / "revision" / "network_7hub"
SEED = 20261019
N_PERMUTATIONS = 200
R_SCRIPT = OUT / "audit_edgeR_sensitivities.R"
EXACT_SCRIPT = OUT / "run_exact_ecm1_permutations_parallel.R"
DECONT_X_SCRIPT = OUT / "run_decontx_ecm1.R"
AMBIENT_SCRIPT = OUT / "run_ambient_ecm1_screen.py"
PRE_QC_SCRIPT = OUT / "run_pre_qc_ecm1_check.py"
RSCRIPT = Path(r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe")
EXACT_LINEAGES = ("myeloid", "mesothelial", "epithelial")
MANUSCRIPT = NETWORK / "MANUSCRIPT_WRITEUP_HUB_CONTEXT_AND_NETWORK.md"
RAW_GSE248762 = ROOT / "data" / "raw" / "GSE248762_extracted"

INPUTS = [
    ROOT / "results" / "revision" / "hub_gene_network_drug_gsea" / "validated7_kegg_gsea_summary.csv",
    NETWORK / "sc7_hub_gene_lineage_results.csv",
    SINGLE_CELL / "pseudobulk_counts" / "counts_myeloid.csv",
    SINGLE_CELL / "pseudobulk_counts" / "counts_mesothelial.csv",
    SINGLE_CELL / "pseudobulk_counts" / "counts_epithelial.csv",
    SINGLE_CELL / "sample_lineage_cell_counts.csv",
    SINGLE_CELL / "sample_group_qc_audit.csv",
    SINGLE_CELL / "cell_type_annotation_table.csv",
    SINGLE_CELL / "hub_gene_expression_per_cell.csv",
    SINGLE_CELL / "04_pseudobulk_edger.R",
    SINGLE_CELL / "edger_gene_lineage_results_raw.csv",
]

R_OUTPUTS = [
    "sc7_nominal_pair_sample_counts.csv",
    "sc7_pair_group_diagnostics.csv",
    "sc7_fold_change_comparison.csv",
    "sc7_leave_one_out.csv",
    "sc7_three_group_descriptive.csv",
    "sc7_confounding_correlations.csv",
    "sc7_sample_lineage_expression.csv",
    "sc7_negative_control_all_genes.csv",
    "sc7_permutation_results.csv",
    "edgeR_software_versions.csv",
]

PYTHON_READ_OUTPUTS = [
    "sc7_nominal_pair_sample_counts.csv",
    "sc7_pair_group_diagnostics.csv",
    "sc7_fold_change_comparison.csv",
    "sc7_leave_one_out.csv",
    "sc7_three_group_descriptive.csv",
    "sc7_confounding_correlations.csv",
    "sc7_negative_control_all_genes.csv",
    "sc7_permutation_results.csv",
    "edgeR_software_versions.csv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, index=False)


def format_value(value: object, digits: int = 15) -> str:
    if pd.isna(value):
        return "NOT VERIFIED"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}g}"
    return str(value)


def format_sig(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return "NOT VERIFIED"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}g}"
    return str(value)


def make_summaries() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample = pd.read_csv(OUT / "sc7_nominal_pair_sample_counts.csv")
    groups = pd.read_csv(OUT / "sc7_pair_group_diagnostics.csv")
    fc = pd.read_csv(OUT / "sc7_fold_change_comparison.csv")
    loo = pd.read_csv(OUT / "sc7_leave_one_out.csv")
    permutations = pd.read_csv(OUT / "sc7_permutation_results.csv")
    nominal = pd.read_csv(NETWORK / "sc7_hub_gene_lineage_results.csv")
    nominal = nominal[
        nominal["status"].eq("testable")
        & nominal["P_value"].lt(0.05)
    ].copy()
    nominal["pair_key"] = nominal["gene"].astype(str) + "|" + nominal["lineage"].astype(str)

    if not permutations.empty:
        permutation_summary = (
            permutations.groupby(["gene", "lineage"], as_index=False)
            .agg(
                n_permutations=("permutation", "nunique"),
                n_abs_permuted_logfc_ge_observed=(
                    "permuted_log2FC",
                    lambda values: int(
                        np.sum(
                            np.abs(values.to_numpy())
                            >= abs(float(permutations.loc[values.index[0], "observed_log2FC"]))
                        )
                    ),
                ),
                observed_log2FC=("observed_log2FC", "first"),
                permuted_log2FC_median=("permuted_log2FC", "median"),
                permuted_log2FC_min=("permuted_log2FC", "min"),
                permuted_log2FC_max=("permuted_log2FC", "max"),
            )
        )
        permutation_summary["empirical_two_sided_p"] = (
            permutation_summary["n_abs_permuted_logfc_ge_observed"] + 1
        ) / (permutation_summary["n_permutations"] + 1)
    else:
        permutation_summary = pd.DataFrame()
    write_csv(permutation_summary, "sc7_permutation_summary.csv")
    exact_summary = pd.concat(
        [
            pd.read_csv(OUT / f"sc7_exact_permutation_summary_{lineage}.csv")
            for lineage in EXACT_LINEAGES
        ],
        ignore_index=True,
    )

    summaries = []
    for pair in nominal.itertuples(index=False):
        key = f"{pair.gene}|{pair.lineage}"
        pair_sample = sample[(sample["gene"] == pair.gene) & (sample["lineage"] == pair.lineage)]
        group_rows = groups[(groups["gene"] == pair.gene) & (groups["lineage"] == pair.lineage)]
        fc_row = fc[(fc["gene"] == pair.gene) & (fc["lineage"] == pair.lineage)].iloc[0]
        loo_row = loo[(loo["gene"] == pair.gene) & (loo["lineage"] == pair.lineage)].iloc[0]
        group_map = {row.contrast_group: row for row in group_rows.itertuples(index=False)}
        perm_mask = (
            exact_summary["gene"].astype(str).str.upper().eq(str(pair.gene).upper())
            & exact_summary["lineage"].eq(pair.lineage)
        )
        perm_row = (
            exact_summary[perm_mask].iloc[0]
            if perm_mask.any()
            else None
        )
        low_count_groups = [
            name for name, row in group_map.items()
            if bool(row.group_mostly_zero)
            or (pd.notna(row.one_sample_supplies_majority_CPM) and bool(row.one_sample_supplies_majority_CPM))
        ]
        raw_driven = (
            "yes" if len(low_count_groups) == 2
            else "partly" if len(low_count_groups) == 1
            else "no"
        )
        survives_loo = (
            "yes"
            if loo_row["status"] == "tested"
            and bool(loo_row["leave_one_out_nominal_P_lt_0_05"])
            and loo_row["leave_one_out_log2FC"] > 0
            else "no"
            if loo_row["status"] == "tested"
            else "NOT VERIFIED"
        )
        count_pattern = bool(low_count_groups)
        prior_attenuates = (
            abs(fc_row["edgeR_logCPM_prior_count_1_log2FC"])
            < abs(fc_row["saved_log2FC"]) / 2
        )
        inflated = (
            "yes" if count_pattern and prior_attenuates
            else "no" if not count_pattern or not prior_attenuates
            else "NOT VERIFIED"
        )
        loo_sign_holds = (
            loo_row["status"] == "tested"
            and np.sign(loo_row["leave_one_out_log2FC"]) == np.sign(fc_row["saved_log2FC"])
        )
        loo_is_nominal = (
            loo_row["status"] == "tested"
            and bool(loo_row["leave_one_out_nominal_P_lt_0_05"])
        )
        low_abundance = float(fc_row["edgeR_AveLogCPM"]) < 1
        permutation_disagrees = (
            perm_row is not None
            and float(perm_row["exact_two_sided_permutation_P"]) >= 0.05
        )
        if not loo_sign_holds or not loo_is_nominal or all(
            bool(row.group_mostly_zero) for row in group_map.values()
        ):
            verdict = "ARTIFACT"
        elif (
            low_count_groups
            or low_abundance
            or inflated == "yes"
            or permutation_disagrees
        ):
            verdict = "FRAGILE"
        elif float(pair.BH_FDR_7hub_family) >= 0.05 or perm_row is None:
            verdict = "NOMINAL_ONLY_UNCONFIRMED"
        else:
            verdict = "ROBUST"
        if verdict == "ARTIFACT":
            manuscript_wording = "not supported after sensitivity checks"
        elif verdict == "FRAGILE":
            manuscript_wording = "higher in LV, but based on low counts / a few samples (interpret cautiously)"
        elif verdict == "NOMINAL_ONLY_UNCONFIRMED":
            manuscript_wording = (
                "nominal LV-minus-SV result only; BH correction failed and/or "
                "permutation test was not run"
            )
        else:
            manuscript_wording = "higher in LV with consistent counts across samples"
        eligible_lv = pair_sample[
            pair_sample["contrast_group"].eq("LV")
            & ~pair_sample["excluded_by_20_cell_rule"].astype(bool)
        ]
        max_cpm_sample = eligible_lv.loc[eligible_lv["CPM"].idxmax()]
        max_cpm_fraction = float(max_cpm_sample["CPM"] / eligible_lv["CPM"].sum())
        summaries.append(
            {
                "gene": pair.gene,
                "lineage": pair.lineage,
                "saved_log2FC": fc_row["saved_log2FC"],
                "saved_P": pair.P_value,
                "saved_BH_q": pair.BH_FDR_7hub_family,
                "edgeR_AveLogCPM": fc_row["edgeR_AveLogCPM"],
                "edgeR_QL_dispersion": fc_row["edgeR_QL_dispersion"],
                "zero_majority_groups": ";".join(
                    name for name, row in group_map.items() if row.group_mostly_zero
                ),
                "single_sample_majority_CPM_groups": ";".join(
                    name for name, row in group_map.items()
                    if row.one_sample_supplies_majority_CPM
                ),
                "groups_with_low_count_flags": ";".join(low_count_groups),
                "raw_count_driven": raw_driven,
                "leave_one_out_log2FC": loo_row["leave_one_out_log2FC"],
                "leave_one_out_P": loo_row["leave_one_out_P"],
                "leave_one_out_sample_id": loo_row["omitted_sample_id"],
                "leave_one_out_delta_log2FC": (
                    loo_row["leave_one_out_log2FC"] - fc_row["saved_log2FC"]
                ),
                "survives_leave_one_out": survives_loo,
                "highest_LV_CPM_sample_id": max_cpm_sample["sample_id"],
                "highest_LV_CPM_fraction": max_cpm_fraction,
                "permutation_empirical_p": (
                    perm_row["exact_two_sided_permutation_P"]
                    if perm_row is not None
                    else np.nan
                ),
                "permutation_count": (
                    int(perm_row["n_allocations"]) if perm_row is not None else 0
                ),
                "permutation_method": "exact" if perm_row is not None else "NOT RUN",
                "fold_change_inflated_by_low_counts": inflated,
                "verdict": verdict,
                "recommended_manuscript_wording": manuscript_wording,
                "n_samples_all_groups_included": len(pair_sample),
            }
        )
    audit = pd.DataFrame(summaries)
    display_order = {
        ("ECM1", "myeloid"): 0,
        ("ECM1", "mesothelial"): 1,
        ("ECM1", "epithelial"): 2,
        ("TNFSF15", "myeloid"): 3,
        ("EBI3", "myeloid"): 4,
        ("EBI3", "lymphoid"): 5,
    }
    audit["_display_order"] = audit.apply(
        lambda row: display_order.get((row["gene"], row["lineage"]), 99), axis=1
    )
    audit = audit.sort_values("_display_order").drop(columns="_display_order").reset_index(drop=True)
    controls = pd.read_csv(OUT / "sc7_negative_control_all_genes.csv")
    distribution_rows = []
    percentile_rows = []
    for pair in nominal.itertuples(index=False):
        fit_genes = controls[controls["lineage"].eq(pair.lineage)]
        other = fit_genes[fit_genes["gene"].astype(str).str.upper().ne(str(pair.gene).upper())]
        hrow = fit_genes[fit_genes["gene"].astype(str).str.upper().eq(str(pair.gene).upper())].iloc[0]
        for metric in ("logFC", "absolute_logFC", "PValue"):
            if metric == "absolute_logFC":
                vals = pd.to_numeric(other["logFC"], errors="coerce").abs().dropna()
                target = abs(float(hrow["logFC"]))
            else:
                vals = pd.to_numeric(other[metric], errors="coerce").dropna()
                target = float(hrow[metric])
            distribution_rows.append(
                {
                    "gene": pair.gene,
                    "lineage": pair.lineage,
                    "control_metric": metric,
                    "n_other_genes_tested": len(vals),
                    "min": vals.min(),
                    "q25": vals.quantile(0.25),
                    "median": vals.median(),
                    "q75": vals.quantile(0.75),
                    "max": vals.max(),
                }
            )
            percentile_rows.append(
                {
                    "gene": pair.gene,
                    "lineage": pair.lineage,
                    "metric": metric,
                    "hub_value": target,
                    "n_other_genes_tested": len(vals),
                    "ascending_percentile": 100 * float((vals <= target).mean()),
                    "interpretation": (
                        "fraction of other fitted genes with value <= hub value; "
                        "for PValue, smaller percentile means a smaller P value"
                    ),
                }
            )
    negative_distribution = pd.DataFrame(distribution_rows)
    negative_percentiles = pd.DataFrame(percentile_rows)
    write_csv(negative_distribution, "sc7_negative_control_distribution.csv")
    write_csv(negative_percentiles, "sc7_negative_control_percentiles.csv")
    abs_percentiles = negative_percentiles.loc[
        negative_percentiles["metric"].eq("absolute_logFC"),
        ["gene", "lineage", "ascending_percentile"],
    ].rename(columns={"ascending_percentile": "negative_control_abs_logFC_percentile"})
    audit = audit.merge(abs_percentiles, on=["gene", "lineage"], how="left")
    write_csv(audit, "sc7_pair_verdicts.csv")

    sample_lineage = pd.read_csv(OUT / "sc7_nominal_pair_sample_counts.csv")
    common_gene = str(
        nominal.loc[nominal["BH_FDR_7hub_family"].lt(0.05), "gene"].iloc[0]
    )
    common_lineages = list(
        nominal.loc[nominal["BH_FDR_7hub_family"].lt(0.05), "lineage"].unique()
    )
    bh_rows = sample_lineage[
        sample_lineage["gene"].eq(common_gene)
        & sample_lineage["lineage"].isin(common_lineages)
    ].copy()
    matrix = bh_rows.pivot(index="sample_id", columns="lineage", values="CPM")
    groups_by_sample = bh_rows.drop_duplicates("sample_id").set_index("sample_id")[
        "original_group"
    ]
    excluded = bh_rows.pivot(
        index="sample_id", columns="lineage", values="excluded_by_20_cell_rule"
    )
    matrix.columns = [f"{col}_CPM" for col in matrix.columns]
    excluded.columns = [f"{col}_excluded_by_20_cell_rule" for col in excluded.columns]
    matrix.insert(0, "original_group", groups_by_sample.reindex(matrix.index))
    matrix = matrix.join(excluded)
    matrix.insert(0, "sample_id", matrix.index)
    matrix.reset_index(drop=True, inplace=True)
    write_csv(matrix, "sc7_bh_hub_sample_lineage_cpm.csv")
    cpm_columns = [f"{lineage}_CPM" for lineage in common_lineages]
    raw_by_sample = bh_rows.pivot(index="sample_id", columns="lineage", values="raw_gene_count")
    expressed_each = raw_by_sample.reindex(columns=common_lineages).gt(0)
    overlap = pd.DataFrame(
        [
            {
                "gene": common_gene,
                "lineages": ";".join(common_lineages),
                "n_samples_in_any_lineage": len(raw_by_sample),
                "n_samples_expressing_all_three_lineages": int(expressed_each.all(axis=1).sum()),
                "n_samples_expressing_at_least_one_lineage": int(expressed_each.any(axis=1).sum()),
                "sample_ids_expressing_all_three": ";".join(
                    expressed_each.index[expressed_each.all(axis=1)].astype(str)
                ),
                "basis": "raw pseudobulk count > 0; all samples including <20-cell samples shown",
            }
        ]
    )
    write_csv(overlap, "sc7_bh_hub_expression_overlap.csv")
    return audit, negative_distribution, negative_percentiles, permutation_summary


def make_ecm1_subgroup_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    sample = pd.read_csv(OUT / "sc7_nominal_pair_sample_counts.csv")
    ecm1 = sample[sample["gene"].eq("ECM1")].copy()
    ecm1["excluded_by_20_cell_rule"] = ecm1["excluded_by_20_cell_rule"].astype(bool)
    ecm1["eligible_for_original_fit"] = ~ecm1["excluded_by_20_cell_rule"]
    sample_values = ecm1[
        [
            "lineage",
            "sample_id",
            "original_group",
            "n_qc_cells_in_lineage",
            "raw_gene_count",
            "total_library_size",
            "CPM",
            "excluded_by_20_cell_rule",
            "eligible_for_original_fit",
        ]
    ].copy()
    sample_values = sample_values.sort_values(
        ["lineage", "original_group", "sample_id"]
    ).reset_index(drop=True)
    write_csv(sample_values, "sc7_ecm1_subgroup_sample_values.csv")

    summaries = []
    for lineage in ("myeloid", "mesothelial", "epithelial"):
        for group in ("LV_UF", "LV_NOT_UF", "SV"):
            all_group = ecm1[
                ecm1["lineage"].eq(lineage) & ecm1["original_group"].eq(group)
            ]
            eligible = all_group[all_group["eligible_for_original_fit"]]
            summaries.append(
                {
                    "lineage": lineage,
                    "original_group": group,
                    "n_samples_total": int(len(all_group)),
                    "n_samples_eligible": int(len(eligible)),
                    "n_samples_excluded_by_20_cell_rule": int(
                        all_group["excluded_by_20_cell_rule"].sum()
                    ),
                    "n_zero_count_eligible": int((eligible["raw_gene_count"] == 0).sum()),
                    "mean_raw_gene_count_eligible": eligible["raw_gene_count"].mean(),
                    "mean_CPM_eligible": eligible["CPM"].mean(),
                }
            )
    summary = pd.DataFrame(summaries)
    write_csv(summary, "sc7_ecm1_subgroup_summary.csv")
    return sample_values, summary


def make_plots(permutation_summary: pd.DataFrame) -> list[str]:
    sample = pd.read_csv(OUT / "sc7_nominal_pair_sample_counts.csv")
    nominal = pd.read_csv(NETWORK / "sc7_hub_gene_lineage_results.csv")
    nominal = nominal[nominal["status"].eq("testable") & nominal["P_value"].lt(0.05)]
    group_colors = {"SV": "#2878A5", "LV_NOT_UF": "#D19A00", "LV_UF": "#C84B31"}
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    size_min, size_max = sample["total_library_size"].min(), sample["total_library_size"].max()
    for ax, row in zip(axes.flat, nominal.itertuples(index=False)):
        data = sample[(sample["gene"] == row.gene) & (sample["lineage"] == row.lineage)]
        sizes = 24 + 120 * (
            (data["total_library_size"] - size_min) / max(size_max - size_min, 1)
        )
        for group, color in group_colors.items():
            subset = data[data["original_group"].eq(group)]
            subset_sizes = sizes.loc[subset.index]
            ax.scatter(
                subset["original_group"].map({group: {"SV": 0, "LV_NOT_UF": 1, "LV_UF": 2}[group]}),
                subset["log2_CPM_plus_1"],
                s=subset_sizes,
                color=color,
                alpha=0.8,
                edgecolors="black",
                linewidths=0.35,
                label=group,
            )
        ax.set_xticks([0, 1, 2], ["SV", "LV_NOT_UF", "LV_UF"], rotation=25, ha="right")
        ax.set_ylabel("log2(CPM + 1)")
        ax.set_title(f"{row.gene} | {row.lineage}\nP={row.P_value:.3g}, BH q={row.BH_FDR_7hub_family:.3g}")
        ax.grid(axis="y", alpha=0.2)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Per-sample expression for the six nominal hub-gene pairs", y=1.04)
    save_plot(fig, "fig_sc7_nominal_sample_cpm")

    loo = pd.read_csv(OUT / "sc7_leave_one_out.csv")
    bh = nominal[nominal["BH_FDR_7hub_family"].lt(0.05)]
    fig, (ax_fc, ax_p) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    labels = []
    for row in bh.itertuples(index=False):
        labels.append(f"{row.gene} | {row.lineage}")
    ypos = np.arange(len(bh))
    for idx, row in enumerate(bh.itertuples(index=False)):
        loo_row = loo[(loo["gene"] == row.gene) & (loo["lineage"] == row.lineage)].iloc[0]
        ax_fc.scatter([row.log2FC_LV_vs_SV], [idx], color="#C0392B", label="Saved fit" if idx == 0 else None)
        ax_fc.scatter([loo_row["leave_one_out_log2FC"]], [idx], marker="x", s=80, color="#2878A5", label="Highest-count sample removed" if idx == 0 else None)
        ax_p.scatter([row.P_value], [idx], color="#C0392B", label="Saved fit" if idx == 0 else None)
        ax_p.scatter([loo_row["leave_one_out_P"]], [idx], marker="x", s=80, color="#2878A5", label="Highest-count sample removed" if idx == 0 else None)
    ax_fc.axvline(0, color="#555555", linewidth=0.8)
    ax_fc.set_yticks(ypos, labels)
    ax_fc.set_xlabel("log2FC (LV vs SV)")
    ax_fc.set_title("Effect estimate")
    ax_p.axvline(0.05, color="#555555", linestyle="--", linewidth=0.8)
    ax_p.set_yticks(ypos, labels)
    ax_p.set_xlabel("P value")
    ax_p.set_title("Nominal P-value sensitivity")
    ax_p.set_xscale("log")
    ax_p.invert_yaxis()
    ax_fc.invert_yaxis()
    ax_fc.legend(frameon=False, fontsize=8)
    ax_p.legend(frameon=False, fontsize=8)
    save_plot(fig, "fig_sc7_leave_one_out")

    matrix = pd.read_csv(OUT / "sc7_bh_hub_sample_lineage_cpm.csv")
    lineage_cols = [col for col in matrix.columns if col.endswith("_CPM")]
    heat = matrix.set_index("sample_id")[lineage_cols]
    fig, ax = plt.subplots(figsize=(8, max(5, len(heat) * 0.32)), constrained_layout=True)
    image = ax.imshow(np.log2(heat.to_numpy(dtype=float) + 1), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(lineage_cols)), [col[:-4] for col in lineage_cols], rotation=25, ha="right")
    ax.set_yticks(range(len(heat)), heat.index.astype(str))
    ax.set_title("Same hub gene across the three BH-significant lineages")
    ax.set_xlabel("Lineage; color = log2(CPM + 1)")
    fig.colorbar(image, ax=ax, label="log2(CPM + 1)")
    save_plot(fig, "fig_sc7_bh_hub_sample_lineage_heatmap")

    permutations = pd.read_csv(OUT / "sc7_permutation_results.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for ax, row in zip(axes.flat, bh.itertuples(index=False)):
        vals = permutations[
            (permutations["gene"] == row.gene)
            & (permutations["lineage"] == row.lineage)
        ]["permuted_log2FC"]
        ax.hist(vals, bins=20, color="#7AA6C2", edgecolor="white")
        ax.axvline(row.log2FC_LV_vs_SV, color="#C0392B", linewidth=2, label="Observed")
        ax.axvline(-row.log2FC_LV_vs_SV, color="#C0392B", linewidth=1, linestyle="--")
        ax.set_title(f"{row.gene} | {row.lineage}")
        ax.set_xlabel("Permuted log2FC")
        ax.set_ylabel("Permutation count")
        ax.legend(frameon=False)
    save_plot(fig, "fig_sc7_permutation_logfc")
    return [
        "fig_sc7_nominal_sample_cpm",
        "fig_sc7_leave_one_out",
        "fig_sc7_bh_hub_sample_lineage_heatmap",
        "fig_sc7_permutation_logfc",
    ]


def save_plot(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def create_report(
    audit: pd.DataFrame,
    negative_distribution: pd.DataFrame,
    negative_percentiles: pd.DataFrame,
    permutation_summary: pd.DataFrame,
    subgroup_summary: pd.DataFrame,
) -> None:
    fc = pd.read_csv(OUT / "sc7_fold_change_comparison.csv")
    groups = pd.read_csv(OUT / "sc7_pair_group_diagnostics.csv")
    loo = pd.read_csv(OUT / "sc7_leave_one_out.csv")
    confounding = pd.read_csv(OUT / "sc7_confounding_correlations.csv")
    three_group = pd.read_csv(OUT / "sc7_three_group_descriptive.csv")
    overlap = pd.read_csv(OUT / "sc7_bh_hub_expression_overlap.csv").iloc[0]
    versions = pd.read_csv(OUT / "edgeR_software_versions.csv")
    saved = pd.read_csv(NETWORK / "sc7_hub_gene_lineage_results.csv")
    nominal = saved[saved["status"].eq("testable") & saved["P_value"].lt(0.05)].copy()
    bh = nominal[nominal["BH_FDR_7hub_family"].lt(0.05)]
    no_bh = nominal[nominal["BH_FDR_7hub_family"].ge(0.05)]
    sample_counts = pd.read_csv(OUT / "sc7_nominal_pair_sample_counts.csv")
    exact_summary = pd.concat(
        [
            pd.read_csv(OUT / f"sc7_exact_permutation_summary_{lineage}.csv")
            for lineage in ("myeloid", "mesothelial", "epithelial")
        ],
        ignore_index=True,
    )
    decontx_contrasts = pd.read_csv(OUT / "sc7_decontx_contrast_results.csv")
    decontx_by_sample = pd.read_csv(OUT / "sc7_decontx_contamination_by_sample.csv")
    ambient_markers = pd.read_csv(OUT / "sc7_ecm1_ambient_marker_screen.csv")
    per_cell_ecm1 = pd.read_csv(OUT / "sc7_ecm1_cell_expression_by_sample_lineage.csv")
    pre_qc_ecm1 = pd.read_csv(OUT / "sc7_ecm1_pre_qc_excluded_summary.csv")

    lines = [
        "# GSE248762 seven-hub fold-change audit",
        "",
        "## Scope and reproducibility",
        "",
        "The analysis gene set was loaded from "
        "`results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`. "
        "Only the six nominally significant saved hub-gene/lineage pairs were audited. "
        "The 3 BH-significant pairs were prioritized. No gene symbols are used in narrative text.",
        "",
        f"The original edgeR quasi-likelihood and Monte Carlo results were generated in an earlier "
        f"full audit. This metadata-only refresh did not rerun edgeR, DecontX, or permutations. "
        f"The historical Monte Carlo outputs are retained below for provenance only and are "
        f"superseded by the separately completed exact enumeration. "
        f"R {versions.loc[versions.software.eq('R'), 'version'].iloc[0]} "
        f"and edgeR {versions.loc[versions.software.eq('edgeR'), 'version'].iloc[0]} were used.",
        "",
        "The original edgeR and Monte Carlo result files were produced in an earlier full audit "
        "run and reused unchanged. Exact permutation and DecontX sensitivity outputs were "
        "generated separately before this metadata-only refresh. SHA-256 checksums document "
        "current file contents and are not proof of which run generated reused outputs.",
        "",
        "Original eligibility was retained: samples needed at least 20 QC cells per lineage; "
        "a hub gene needed at least two samples per group with at least 10 raw pseudobulk UMIs. "
        "The edgeR design was `~0 + group`, with `LV - SV` contrast, `filterByExpr`, robust "
        "dispersion estimation, and robust `glmQLFit`, matching "
        "`04_pseudobulk_edger.R`.",
        "",
        "## Claim-by-claim findings",
        "",
        "| Claim | Finding | Source / calculation |",
        "|---|---|---|",
        f"| Saved comparisons audited | {len(nominal)} nominal pairs; {len(bh)} BH-significant and {len(no_bh)} nominal-only | `sc7_hub_gene_lineage_results.csv` |",
        f"| Saved edgeR output completeness | The saved raw file contains hub-gene summary rows, not the all-gene fit tables or per-gene quasi-likelihood dispersions; all-gene control fits and dispersions were generated during the earlier full audit, not this metadata-only refresh | `edger_gene_lineage_results_raw.csv`; `04_pseudobulk_edger.R` |",
        f"| BH-associated hub gene across three lineages | The same hub gene occurs in all three BH-significant rows; {int(overlap['n_samples_expressing_all_three_lineages'])} samples had nonzero raw pseudobulk counts in all three lineages | `sc7_bh_hub_expression_overlap.csv`; `sc7_bh_hub_sample_lineage_cpm.csv` |",
        "| Historical permutation outputs | Earlier Monte Carlo files are retained for provenance only and are superseded by exact enumeration for the ECM1 comparisons | `sc7_permutation_results.csv`; `sc7_exact_permutation_summary_<lineage>.csv` |",
        "",
        "The sample-count table includes all samples, including samples excluded from fitting by "
        "the 20-cell rule, and explicitly marks such samples. Group summaries use only samples "
        "meeting that rule. `sc7_pair_group_diagnostics.csv` contains zero/low/high count bins, "
        "raw-count and CPM summaries, library-size ranges, and maximum-sample contribution.",
        "",
        "### Pair-level verdict table",
        "",
        "| Hub | Lineage | Saved log2FC | EdgeR P | BH q | Exact permutation P (allocations) | Exact P - edgeR P; nominal agreement | LOO log2FC / P; omitted sample | Highest LV CPM sample (% LV CPM) | |log2FC| negative-control percentile | Verdict |",
        "|---|---|---:|---:|---:|---|---|---|---|---:|---|",
    ]
    for row in audit.itertuples(index=False):
        permutation_result = (
            f"{format_value(row.permutation_empirical_p)} ({row.permutation_count})"
            if row.permutation_count > 0
            else "NOT RUN"
        )
        if row.permutation_count > 0:
            permutation_delta = (
                f"{format_value(row.permutation_empirical_p - row.saved_P)}; "
                f"{'agree' if (row.permutation_empirical_p < 0.05) == (row.saved_P < 0.05) else 'disagree'} "
                "at P < 0.05"
            )
        else:
            permutation_delta = "NOT CALCULATED; no permutation test"
        lines.append(
            f"| {row.gene} | {row.lineage} | {format_value(row.saved_log2FC)} | "
            f"{format_value(row.saved_P)} | {format_value(row.saved_BH_q)} | "
            f"{permutation_result} | {permutation_delta} | "
            f"{format_value(row.leave_one_out_log2FC)} / {format_value(row.leave_one_out_P)}; "
            f"{row.leave_one_out_sample_id} | "
            f"{row.highest_LV_CPM_sample_id} ({format_value(100 * row.highest_LV_CPM_fraction)}%) | "
            f"{format_value(row.negative_control_abs_logFC_percentile)} | {row.verdict} |"
        )
    lines.extend(
        [
            "",
            "### Verdict rules applied by the current script",
            "",
            "The exact classification logic is:",
            "",
            "```python",
            "if not loo_sign_holds or not loo_is_nominal or all(",
            "    bool(row.group_mostly_zero) for row in group_map.values()",
            "):",
            '    verdict = "ARTIFACT"',
            "elif (",
            "    low_count_groups",
            "    or low_abundance",
            '    or inflated == "yes"',
            "    or permutation_disagrees",
            "):",
            '    verdict = "FRAGILE"',
            "elif float(pair.BH_FDR_7hub_family) >= 0.05 or perm_row is None:",
            '    verdict = "NOMINAL_ONLY_UNCONFIRMED"',
            "else:",
            '    verdict = "ROBUST"',
            "```",
            "",
            "`loo_is_nominal` means LOO P < 0.05; `low_abundance` means AveLogCPM < 1; "
            "`low_count_groups` includes groups with >50% zero counts or one sample contributing "
            ">50% of that group's CPM; `inflated` is yes when this count pattern exists and the "
            "absolute prior-count-1 logFC is <50% of the saved magnitude; and "
            "`permutation_disagrees` means an available permutation P is >=0.05. This ordering "
            "prevents a pair that failed BH or lacks a permutation test from being called ROBUST. "
            "The earlier saved table called the two nominal-only EBI3 pairs ROBUST; the current "
            "label for both is NOMINAL_ONLY_UNCONFIRMED.",
            "",
            "## LV subgroup check for the shared BH-associated hub",
            "",
            "Means below are arithmetic means of saved per-sample CPM and raw counts among "
            "samples meeting the original >=20-cell lineage rule. Zero counts are among eligible "
            "samples. These raw-CPM summaries do not refit the original subgroup model; "
            "separate DecontX-corrected contrasts are reported below. The per-sample source "
            "values, cell counts, and exclusion flags are in `sc7_ecm1_subgroup_sample_values.csv`.",
            "",
            "| Lineage | Original group | Eligible n | Excluded n | Mean raw count | Mean CPM | Zero counts |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in subgroup_summary.itertuples(index=False):
        lines.append(
            f"| {row.lineage} | {row.original_group} | {row.n_samples_eligible} | "
            f"{row.n_samples_excluded_by_20_cell_rule} | "
            f"{format_value(row.mean_raw_gene_count_eligible)} | "
            f"{format_value(row.mean_CPM_eligible)} | {row.n_zero_count_eligible} |"
        )
    lines.extend(
        [
            "",
            "LV_UF means exceed LV_NOT_UF and SV means in all three lineages, and the LV_UF CPM "
            "distributions include samples that dominate the pooled signal. Thus, the pooled LV "
            "versus SV contrast could be driven by the LV subgroup mix. The subgroup contrasts "
            "remain descriptive and underpowered; the CPM means are sensitive to individual "
            "samples.",
            "",
            "## Permutation-space feasibility (calculation only)",
            "",
            "The previous 200-draw tests sampled label assignments and used a plus-one "
            "correction. The approved exact test enumerated every fixed-size label allocation "
            "and used the inclusive two-sided tail count divided by the total number of "
            "allocations, without a plus-one correction.",
            "",
            "| Lineage | Eligible LV | Eligible SV | Distinct assignments | Exact extreme allocations | Exact P | EdgeR P | Exact minus edgeR P |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    ecm_samples = sample_counts[
        sample_counts["gene"].eq("ECM1")
        & ~sample_counts["excluded_by_20_cell_rule"].astype(bool)
    ]
    for lineage in ("myeloid", "mesothelial", "epithelial"):
        lineage_samples = ecm_samples[ecm_samples["lineage"].eq(lineage)]
        n_lv = int(lineage_samples["contrast_group"].eq("LV").sum())
        n_sv = int(lineage_samples["contrast_group"].eq("SV").sum())
        n_allocations = math.comb(n_lv + n_sv, n_lv)
        edge_p = float(
            saved.loc[
                saved["gene"].eq("ECM1") & saved["lineage"].eq(lineage),
                "P_value",
            ].iloc[0]
        )
        exact = exact_summary[exact_summary["lineage"].eq(lineage)].iloc[0]
        lines.append(
            f"| {lineage} | {n_lv} | {n_sv} | {n_allocations} | "
            f"{int(exact.n_abs_permuted_log2FC_ge_observed)} | "
            f"{format_value(exact.exact_two_sided_permutation_P)} | "
            f"{format_value(edge_p)} | "
            f"{format_value(exact.exact_two_sided_permutation_P - edge_p)} |"
        )
    lines.extend(
        [
            "",
            "These exact P values supersede the 200-draw Monte Carlo estimates for the "
            "sample-label permutation check. The edgeR P values remain model-based "
            "quasi-likelihood tests and use a different null calibration.",
            "",
            "## Raw-count and leave-one-out audit",
            "",
            "For each of the six pairs, `sc7_nominal_pair_sample_counts.csv` reports one row per "
            "sample with original and contrast group, lineage cell count, raw hub count, total "
            "library size, CPM, log2(CPM+1), and the 20-cell exclusion flag. Group-level "
            "zero/1–5/>5 bins and count/CPM summaries are in `sc7_pair_group_diagnostics.csv`.",
            "",
            "A group is flagged as mostly zero when more than half of its eligible samples have "
            "zero counts. A single sample is flagged as supplying a majority of signal when its "
            "CPM exceeds half the sum of eligible sample CPMs in that group. The verdict "
            "`raw_count_driven` is yes when both groups have at least one such flag, partly when "
            "one group does, and no otherwise. The low-count inflation flag is yes when a count "
            "pattern is present and the magnitude of the prior-count-1 logFC is less than half "
            "the saved logFC magnitude. These are audit decision rules, not biological "
            "significance thresholds.",
            "",
            "| Hub | Lineage | Omitted highest-count sample | Omitted group | Omitted count | LOO log2FC | LOO P |",
            "|---|---|---|---|---:|---:|---:|",
        ]
    )
    for row in loo.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {row.lineage} | {row.omitted_sample_id} | "
            f"{row.omitted_original_group} | {format_value(row.omitted_raw_count)} | "
            f"{format_value(row.leave_one_out_log2FC)} | {format_value(row.leave_one_out_P)} |"
        )
    lines.extend(
        [
            "",
            "## Fold-change scale and abundance",
            "",
            "The following table compares the saved edgeR log2FC with the unnormalized "
            "mean-CPM ratio using a 1-CPM pseudocount and with the mean group difference in "
            "edgeR logCPM from `cpm(prior.count=1)`. `edgeR_AveLogCPM` is the fitted average "
            "abundance from the edgeR quasi-likelihood table. Quasi-likelihood dispersion is "
            "reported from `glmQLFit`; it was not present in the saved raw hub-only results.",
            "",
            "| Hub | Lineage | Saved log2FC | Mean-CPM ratio log2FC (+1 CPM) | Prior-count-1 edgeR logCPM group difference | AveLogCPM | QL dispersion |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in fc.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {row.lineage} | {format_value(row.saved_log2FC)} | "
            f"{format_value(row.simple_CPM_ratio_log2FC_pseudocount_1)} | "
            f"{format_value(row.edgeR_logCPM_prior_count_1_log2FC)} | "
            f"{format_value(row.edgeR_AveLogCPM)} | {format_value(row.edgeR_QL_dispersion)} |"
        )
    lines.extend(
        [
            "",
            "## Library size, cell count and LV subtype checks",
            "",
            "`sc7_confounding_correlations.csv` reports Pearson correlation of the binary "
            "LV/SV indicator with eligible lineage cell count and library size, calculated "
            "within each pair's lineage. These descriptive correlations do not constitute "
            "a covariate-adjusted model.",
            "",
            "| Hub | Lineage | r(group, lineage cells) | r(group, library size) |",
            "|---|---|---:|---:|",
        ]
    )
    for row in confounding.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {row.lineage} | "
            f"{format_value(row.Pearson_r_group_vs_lineage_cell_count)} | "
            f"{format_value(row.Pearson_r_group_vs_library_size)} |"
        )
    lines.extend(
        [
            "",
            "The three-group descriptive edgeR contrasts compare SV separately with LV_NOT_UF "
            "and LV_UF. They use smaller groups and are underpowered; they are descriptive, "
            "not confirmatory. Full results, including non-significant P values and any "
            "non-testable contrasts, are in `sc7_three_group_descriptive.csv`.",
            "",
            "| Hub | Lineage | Contrast | n comparison group | n SV | log2FC | P | Status |",
            "|---|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in three_group.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {row.lineage} | {row.contrast} | "
            f"{format_value(row.n_comparison_group)} | {format_value(row.n_SV)} | "
            f"{format_value(row.log2FC)} | {format_value(row.P_value)} | {row.status} |"
        )
    lines.extend(
        [
            "",
            "The sample-by-lineage CPM matrix for the hub gene shared by the three BH results "
            "is in `sc7_bh_hub_sample_lineage_cpm.csv`; samples with fewer than 20 cells remain "
            "in that matrix but are flagged. Expression overlap is defined as raw count > 0. "
            "The overlap count above indicates whether the nominally repeated lineage results "
            "are driven by a shared set of samples; it does not establish statistical "
            "independence.",
            "",
            "## Negative controls: all fitted genes",
            "",
            "The saved edgeR raw file did not store all tested genes. The full three-lineage "
            "fit was therefore recomputed using the same filtering, robust dispersion "
            "estimation, and quasi-likelihood contrast as the source R script. The complete "
            "tested-gene distributions are in `sc7_negative_control_all_genes.csv`; quantiles "
            "and hub-gene percentile ranks are in the following tables.",
            "",
            "| Hub | Lineage | Metric | Other genes tested | Minimum | Q25 | Median | Q75 | Maximum |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in negative_distribution.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {row.lineage} | {row.control_metric} | "
            f"{row.n_other_genes_tested} | {format_value(row.min)} | {format_value(row.q25)} | "
            f"{format_value(row.median)} | {format_value(row.q75)} | {format_value(row.max)} |"
        )
    lines.extend(
        [
            "",
            "| Hub | Lineage | Metric | Hub value | Ascending percentile among other genes |",
            "|---|---|---|---:|---:|",
        ]
    )
    for row in negative_percentiles.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {row.lineage} | {row.metric} | "
            f"{format_value(row.hub_value)} | {format_value(row.ascending_percentile, 5)} |"
        )
    lines.extend(
        [
            "",
            "Percentiles are the fraction of other fitted genes with a value less than or equal "
            "to the hub gene's value. For raw P values, a lower percentile corresponds to a "
            "smaller P value. These genome-wide distributions are contextual negative controls, "
            "not independent validation.",
            "",
            "## Historical Monte Carlo permutation check (SUPERSEDED BY EXACT ENUMERATION)",
            "",
            "The following earlier-run results are retained for provenance only and must not be "
            "used in place of the exact-test results above. Labels were permuted within lineage "
            "while preserving the observed LV/SV sample "
            "counts, and each permuted dataset was refit using edgeR quasi-likelihood with the "
            "original filtering approach. The empirical two-sided P value is "
            "`(1 + number of |permuted log2FC| >= |observed log2FC|) / (1 + number of valid "
            "permutations)`. Complete permuted effects are in `sc7_permutation_results.csv`.",
            "",
            "| Hub | Lineage | Valid permutations | Observed log2FC | Permutation median | Empirical P |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in permutation_summary.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {row.lineage} | {row.n_permutations} | "
            f"{format_value(row.observed_log2FC)} | {format_value(row.permuted_log2FC_median)} | "
            f"{format_value(row.empirical_two_sided_p)} |"
        )
    lines.extend(
        [
            "",
            "## Verdict and recommended wording",
            "",
        ]
    )
    for row in audit.itertuples(index=False):
        permutation_result = (
            format_value(row.permutation_empirical_p)
            if row.permutation_count > 0
            else "NOT RUN"
        )
        lines.append(
            f"- **{row.gene} / {row.lineage}:** {row.recommended_manuscript_wording} "
            f"(raw-count driven: {row.raw_count_driven}; leave-one-out: "
            f"{row.survives_leave_one_out}; permutation P: "
            f"{permutation_result})."
        )
    lines.extend(
        [
            "",
            "## Ambient-RNA screen and DecontX sensitivity",
            "",
            "The extracted GEO directory contained the sample matrix, barcode, and feature files "
            "for the cell-called matrices, but no separate raw all-droplet or empty-droplet "
            "matrix. An empirical empty-droplet ambient profile is therefore NOT VERIFIED. "
            "The local `data/raw/GSE248762_series_matrix.txt.gz` file is 0 bytes and contains "
            "no readable series-matrix records or sample metadata. The official GEO download "
            "location is `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248762/matrix/GSE248762_series_matrix.txt.gz`; "
            "batch, donor, and collection-date fields are not inferred from this empty local copy. "
            "The raw matrices do contain pre-QC cell-called barcodes removed by the saved "
            "quality filter. Those excluded barcodes are not empty droplets, do not have saved "
            "lineage labels, and are summarized separately as a low-quality-cell check.",
            "DecontX was run on QC-retained cells separately by sample, using the saved cell-type "
            "labels as clusters and seed 20261019. No empty droplets were supplied. DecontX "
            "corrected expected counts were rounded as documented by the package, summed to "
            "sample-lineage pseudobulk, and tested with robust edgeR quasi-likelihood models. "
            "This model-based sensitivity analysis is not proof that removed signal was ambient.",
            "",
            "The per-cell expression export stores log1p(CPM10k). In the four LV_UF myeloid "
            "samples, hub-gene detection and expression were:",
            "",
            "| Sample | Myeloid cells | Positive cells | Percent positive | Mean log1p(CPM10k) |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    uf_myeloid = per_cell_ecm1[
        per_cell_ecm1["lineage"].eq("myeloid")
        & per_cell_ecm1["original_group"].eq("LV_UF")
    ].sort_values("sample_id")
    for row in uf_myeloid.itertuples(index=False):
        lines.append(
            f"| {row.sample_id} | {row.n_cells} | {row.n_expressing_cells} | "
            f"{format_sig(row.percent_expressing_cells)}% | "
            f"{format_sig(row.mean_log1p_CPM10k)} |"
        )
    lines.extend(
        [
            "",
            "In pre-QC cell-called barcodes excluded by the saved filter, no counts were "
            "observed in the first LV_UF sample; the two LV_UF samples with the strongest "
            "excluded-barcode signal had 12.0% and 19.2% positive barcodes, respectively. "
            "These observations do not estimate ambient RNA because the barcodes are not empty "
            "droplets and have no lineage assignments.",
            "",
            "| Sample | Group | Excluded barcodes | Positive (%) | Mean CPM | Maximum raw count |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in pre_qc_ecm1.sort_values("sample_id").itertuples(index=False):
        lines.append(
            f"| {row.sample_id} | {row.original_group} | "
            f"{row.qc_excluded_called_barcodes} | "
            f"{format_sig(row.ecm1_positive_percent_excluded)} | "
            f"{format_sig(row.ecm1_mean_cpm_excluded)} | "
            f"{row.ecm1_max_raw_count_excluded} |"
        )
    lines.extend(
        [
            "",
            "The strongest per-cell myeloid signal was concentrated in one LV_UF sample; that "
            "sample also had high expression in mesothelial and epithelial cells. The saved "
            "annotations assign cells by broad marker sets without using this hub gene, but "
            "cell labels cannot rule out ambient RNA or mixed-cell effects.",
            "",
            "An exact one-sided LV_UF-versus-other-samples marker screen used sample-level myeloid "
            "CPM and all 1,820 allocations of four LV_UF labels among 16 samples. This is a "
            "marker-pattern screen, not a decontamination test.",
            "",
            "| Gene | Mean CPM LV_UF | Mean CPM others | Exact P | BH q (8 genes) |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in ambient_markers.itertuples(index=False):
        lines.append(
            f"| {row.gene} | {format_sig(row.mean_CPM_LV_UF)} | "
            f"{format_sig(row.mean_CPM_other)} | "
            f"{format_sig(row.one_sided_exact_permutation_P)} | "
            f"{format_sig(row.BH_q_over_eight_screened_genes)} |"
        )
    lines.extend(
        [
            "",
            "The four broad fibroblast/stromal markers were not elevated in LV_UF in this screen. "
            "Mesothelial and keratin markers had higher LV_UF means, but none of the eight "
            "screened genes had BH q < 0.05. This pattern neither demonstrates nor excludes "
            "ambient contamination.",
            "",
            "DecontX returned per-cell contamination estimates and corrected-count contrasts:",
            "",
            "| Sample | Group | QC cells | Mean estimate | Median estimate |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in decontx_by_sample.sort_values("sample_id").itertuples(index=False):
        lines.append(
            f"| {row.sample_id} | {row.original_group} | {row.n_qc_cells} | "
            f"{format_sig(row.mean_contamination)} | {format_sig(row.median_contamination)} |"
        )
    lines.extend(
        [
            "",
            "| Lineage | Contrast | n contrast / SV | Before log2FC / P | After log2FC / P | Status |",
            "|---|---|---:|---:|---:|---|",
        ]
    )
    for row in decontx_contrasts.itertuples(index=False):
        lines.append(
            f"| {row.lineage} | {row.contrast} | "
            f"{row.n_contrast_samples} / {row.n_SV_samples} | "
            f"{format_sig(row.before_log2FC)} / {format_sig(row.before_P)} | "
            f"{format_sig(row.after_log2FC)} / {format_sig(row.after_P)} | "
            f"{row.after_status} |"
        )
    lines.extend(
        [
            "",
            "BH q-values for the DecontX contrasts are NOT VERIFIED and are not reported. "
            "No DecontX BH q-value table is included.",
        ]
    )
    lines.extend(
        [
            "",
            "Supporting files include the full sample-by-lineage cell-expression and pseudobulk "
            "CPM tables, per-sample myeloid marker CPMs, per-cell DecontX estimates, corrected "
            "pseudobulk counts, pre-QC excluded-barcode summaries, and per-lineage exact "
            "permutation results; see `FILE_INVENTORY.md`.",
        ]
    )
    lines.extend(
        [
            "",
            "All findings remain exploratory because group sizes are small (approximately "
            "6 SV and up to 10 LV samples), multiple lineages have limited eligible samples, "
            "and pseudobulk inference depends on sample-level rather than cell-level replication. "
            "Participant identifiers are not reported in the saved sample metadata, so samples "
            "may not be independent donors and donor-level independence cannot be verified.",
            "",
            "## Figures",
            "",
            "![Per-sample CPM by group; point size reflects library size](fig_sc7_nominal_sample_cpm.png)",
            "",
            "![Leave-one-highest-count-sample-out effect and P-value](fig_sc7_leave_one_out.png)",
            "",
            "![Sample-by-lineage CPM heatmap for the hub shared by the three BH results](fig_sc7_bh_hub_sample_lineage_heatmap.png)",
            "",
            "![Historical Monte Carlo log2FC distributions (superseded by exact enumeration; provenance only)](fig_sc7_permutation_logfc.png)",
            "",
            "## Output tables",
            "",
            "- `sc7_nominal_pair_sample_counts.csv`: per-sample counts, libraries, CPM, and exclusions.",
            "- `sc7_pair_group_diagnostics.csv`: group count bins and count/CPM/library summaries.",
            "- `sc7_fold_change_comparison.csv`: saved, simple CPM, prior-count logFC, AveLogCPM, dispersion.",
            "- `sc7_leave_one_out.csv`: highest-count-sample removal sensitivity for all nominal pairs.",
            "- `sc7_three_group_descriptive.csv`: LV_NOT_UF and LV_UF contrasts, kept descriptive.",
            "- `sc7_ecm1_subgroup_summary.csv`: descriptive CPM/count summaries by original LV subgroup.",
            "- `sc7_ecm1_subgroup_sample_values.csv`: per-sample expression and eligibility by lineage.",
            "- `sc7_confounding_correlations.csv`: group correlation with lineage cells and library size.",
            "- `sc7_bh_hub_sample_lineage_cpm.csv`: sample-by-lineage expression matrix for shared BH hub.",
            "- `sc7_negative_control_all_genes.csv`: all genes retained/tested in the three edgeR fits.",
            "- `sc7_negative_control_distribution.csv`, `sc7_negative_control_percentiles.csv`.",
            "- `sc7_permutation_results.csv`, `sc7_permutation_summary.csv`.",
            "- `sc7_exact_permutation_summary_<lineage>.csv`: exhaustive fixed-size label tests.",
            "- `sc7_ecm1_ambient_marker_screen.csv`: exact LV_UF sample-label marker screen.",
            "- `sc7_ecm1_pre_qc_excluded_summary.csv` and `sc7_ecm1_pre_qc_excluded_barcodes.csv`: low-quality cell-called barcodes removed by QC; not empty droplets.",
            "- `sc7_decontx_contrast_results.csv`: saved and DecontX-corrected contrasts.",
        ]
    )
    (OUT / "FC_AUDIT_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run or summarize the GSE248762 seven-hub fold-change audit."
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help=(
            "Regenerate reports, plots, inventory, and checksums from the saved original "
            "audit plus completed exact-permutation and DecontX outputs."
        ),
    )
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    required_inputs = INPUTS + [
        R_SCRIPT,
        EXACT_SCRIPT,
        DECONT_X_SCRIPT,
        AMBIENT_SCRIPT,
        PRE_QC_SCRIPT,
    ]
    if not args.metadata_only:
        required_inputs.append(RSCRIPT)
    missing = [str(path) for path in required_inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required audit inputs/tools:\n" + "\n".join(missing))
    if args.metadata_only:
        missing_saved_outputs = [
            str(OUT / name)
            for name in R_OUTPUTS
            if not (OUT / name).is_file()
        ]
        sensitivity_outputs = [
            *(OUT / f"sc7_exact_permutation_results_{lineage}.csv" for lineage in EXACT_LINEAGES),
            *(OUT / f"sc7_exact_permutation_summary_{lineage}.csv" for lineage in EXACT_LINEAGES),
            OUT / "sc7_decontx_contrast_results.csv",
            OUT / "sc7_decontx_contamination_by_sample.csv",
            OUT / "sc7_ecm1_ambient_marker_screen.csv",
            OUT / "sc7_ecm1_cell_expression_by_sample_lineage.csv",
            OUT / "sc7_ecm1_pseudobulk_cpm_by_sample.csv",
            OUT / "sc7_ecm1_ambient_marker_sample_cpm.csv",
            OUT / "sc7_ecm1_pre_qc_excluded_summary.csv",
            OUT / "sc7_ecm1_pre_qc_excluded_barcodes.csv",
        ]
        missing_saved_outputs.extend(
            str(path) for path in sensitivity_outputs if not path.is_file()
        )
        if missing_saved_outputs:
            raise FileNotFoundError(
                "Metadata-only mode requires these existing audit outputs:\n"
                + "\n".join(missing_saved_outputs)
            )
        for lineage in EXACT_LINEAGES:
            exact_results = pd.read_csv(
                OUT / f"sc7_exact_permutation_results_{lineage}.csv",
                usecols=["allocation_index"],
            )
            exact_summary = pd.read_csv(
                OUT / f"sc7_exact_permutation_summary_{lineage}.csv"
            ).iloc[0]
            if (
                len(exact_results) != int(exact_summary["n_allocations"])
                or int(exact_summary["n_fits_succeeded"]) != len(exact_results)
                or exact_results["allocation_index"].duplicated().any()
                or set(exact_results["allocation_index"].astype(int))
                != set(range(1, int(exact_summary["n_allocations"]) + 1))
            ):
                raise ValueError(f"{lineage}: exact allocation results are incomplete or duplicated.")
    else:
        run = subprocess.run(
            [
                str(RSCRIPT),
                str(R_SCRIPT),
                str(ROOT),
                str(OUT),
                str(SEED),
                str(N_PERMUTATIONS),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if run.stdout:
            print(run.stdout, end="")
        if run.stderr:
            print(run.stderr, end="")
        if run.returncode != 0:
            raise RuntimeError(f"edgeR sensitivity script failed with exit code {run.returncode}")

    audit, negative_distribution, negative_percentiles, permutation_summary = make_summaries()
    _, subgroup_summary = make_ecm1_subgroup_tables()
    plots = make_plots(permutation_summary)
    create_report(
        audit,
        negative_distribution,
        negative_percentiles,
        permutation_summary,
        subgroup_summary,
    )

    python_outputs = [
        OUT / "sc7_permutation_summary.csv",
        OUT / "sc7_pair_verdicts.csv",
        OUT / "sc7_negative_control_distribution.csv",
        OUT / "sc7_negative_control_percentiles.csv",
        OUT / "sc7_bh_hub_sample_lineage_cpm.csv",
        OUT / "sc7_bh_hub_expression_overlap.csv",
        OUT / "sc7_ecm1_subgroup_sample_values.csv",
        OUT / "sc7_ecm1_subgroup_summary.csv",
        OUT / "sc7_ecm1_cell_expression_by_sample_lineage.csv",
        OUT / "sc7_ecm1_pseudobulk_cpm_by_sample.csv",
        OUT / "sc7_ecm1_ambient_marker_screen.csv",
        OUT / "sc7_ecm1_ambient_marker_sample_cpm.csv",
        OUT / "sc7_decontx_contrast_results.csv",
        OUT / "sc7_decontx_contamination_per_cell.csv",
        OUT / "sc7_decontx_contamination_by_sample.csv",
        OUT / "FC_AUDIT_REPORT.md",
        MANUSCRIPT,
        *[OUT / f"{stem}.{ext}" for stem in plots for ext in ("png", "pdf")],
    ]
    analysis_outputs = [
        *OUT.glob("sc7_exact_permutation_results_*.csv"),
        *OUT.glob("sc7_exact_permutation_summary_*.csv"),
        *OUT.glob("sc7_decontx_*.csv"),
        *OUT.glob("sc7_ecm1_pre_qc_excluded_*.csv"),
        *OUT.glob("decontx_contamination_GSM*.csv"),
        *OUT.glob("decontx_GSM*.log"),
        *(OUT / "decontx_pseudobulk").glob("decontx_sample_counts_GSM*.csv"),
        *[OUT / f"sc7_decontx_counts_{lineage}.csv" for lineage in EXACT_LINEAGES],
    ]
    analysis_outputs.extend(
        [
            OUT / "run_exact_ecm1_permutations.R",
            EXACT_SCRIPT,
            DECONT_X_SCRIPT,
            AMBIENT_SCRIPT,
            PRE_QC_SCRIPT,
            Path(__file__).resolve(),
        ]
    )
    python_outputs.extend(analysis_outputs)
    if not args.metadata_only:
        python_outputs.extend(OUT / name for name in R_OUTPUTS)
    raw_inputs = [path for path in RAW_GSE248762.glob("*") if path.is_file()]
    read_paths = sorted(
        set(
            path.resolve()
            for path in INPUTS
            + [
                R_SCRIPT,
                Path(__file__).resolve(),
                EXACT_SCRIPT,
                OUT / "run_exact_ecm1_permutations.R",
                DECONT_X_SCRIPT,
                AMBIENT_SCRIPT,
                PRE_QC_SCRIPT,
                *raw_inputs,
                OUT / "sc7_ecm1_subgroup_sample_values.csv",
                OUT / "sc7_ecm1_cell_expression_by_sample_lineage.csv",
                OUT / "sc7_ecm1_pseudobulk_cpm_by_sample.csv",
                OUT / "sc7_ecm1_ambient_marker_screen.csv",
                OUT / "sc7_ecm1_ambient_marker_sample_cpm.csv",
                OUT / "sc7_ecm1_pre_qc_excluded_summary.csv",
                OUT / "sc7_ecm1_pre_qc_excluded_barcodes.csv",
                *[
                    OUT / f"sc7_exact_permutation_summary_{lineage}.csv"
                    for lineage in EXACT_LINEAGES
                ],
                OUT / "sc7_decontx_contrast_results.csv",
                OUT / "sc7_decontx_contamination_by_sample.csv",
                OUT / "sc7_decontx_contamination_per_cell.csv",
                MANUSCRIPT,
            ]
            + [OUT / name for name in PYTHON_READ_OUTPUTS]
        )
    )
    written_paths = sorted(
        set(path.resolve() for path in python_outputs if path.exists())
    )
    inventory_text = [
        "# File inventory",
        "",
        "Execution mode: metadata-only report refresh after the original audit, exact "
        "permutation, DecontX, ambient-marker, and pre-QC barcode analyses.",
        "",
        "All paths are repository-relative. Checksums and execution provenance are in "
        "`run_metadata.json`. Checksums record current content and do not prove which run "
        "generated reused outputs.",
        "",
        f"## Files read ({len(read_paths)})",
        "",
        *[f"- `{path.relative_to(ROOT)}`" for path in read_paths],
        "",
        f"## Files written ({len(written_paths) + 2})",
        "",
        *[f"- `{path.relative_to(ROOT)}`" for path in written_paths],
        "- `results\\revision\\sc7_fc_audit\\FILE_INVENTORY.md`",
        "- `results\\revision\\sc7_fc_audit\\run_metadata.json`",
        "",
        "## Earlier full-run outputs reused unchanged",
        "",
        *[f"- `results\\revision\\sc7_fc_audit\\{name}`" for name in R_OUTPUTS],
        "",
        "## Exact permutation outputs",
        "",
        *[
            f"- `{(OUT / f'sc7_exact_permutation_results_{lineage}.csv').relative_to(ROOT)}`"
            for lineage in EXACT_LINEAGES
        ],
        *[
            f"- `{(OUT / f'sc7_exact_permutation_summary_{lineage}.csv').relative_to(ROOT)}`"
            for lineage in EXACT_LINEAGES
        ],
        "- `results\\revision\\sc7_fc_audit\\run_exact_ecm1_permutations_parallel.R`",
        "",
        "## DecontX and ambient-screen outputs",
        "",
        *[
            f"- `{path.relative_to(ROOT)}`"
            for path in analysis_outputs
            if path.name.startswith("sc7_decontx_")
            or path.name.startswith("decontx_")
            or path.name.startswith("sc7_ecm1_")
        ],
        "- `results\\revision\\sc7_fc_audit\\run_decontx_ecm1.R`",
        "- `results\\revision\\sc7_fc_audit\\run_ambient_ecm1_screen.py`",
        "",
        "## Reports and metadata refreshed in this task",
        "",
        *[
            f"- `{path.relative_to(ROOT)}`"
            for path in written_paths
            if path.name in {
                "FC_AUDIT_REPORT.md",
                "MANUSCRIPT_WRITEUP_HUB_CONTEXT_AND_NETWORK.md",
                "FILE_INVENTORY.md",
                "run_metadata.json",
            }
        ],
        "",
    ]
    inventory_path = OUT / "FILE_INVENTORY.md"
    inventory_path.write_text("\n".join(inventory_text), encoding="utf-8")
    written_paths.extend([inventory_path.resolve()])
    package_versions = {}
    for package in ("matplotlib", "numpy", "pandas", "scipy"):
        try:
            package_versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            package_versions[package] = "NOT INSTALLED"
    version_table = pd.read_csv(OUT / "edgeR_software_versions.csv")
    package_versions["R"] = version_table.loc[version_table["software"].eq("R"), "version"].iloc[0]
    package_versions["edgeR"] = version_table.loc[version_table["software"].eq("edgeR"), "version"].iloc[0]
    r_version_query = subprocess.run(
        [
            str(RSCRIPT),
            "-e",
            "cat(paste(vapply(c('celda','Matrix','BiocManager'), function(pkg) as.character(packageVersion(pkg)), character(1)), collapse='|'))",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    for package, version in zip(
        ("celda", "Matrix", "BiocManager"),
        r_version_query.stdout.strip().split("|"),
        strict=True,
    ):
        package_versions[package] = version
    hash_paths = sorted(set(read_paths + written_paths))
    hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in hash_paths
        if path.is_file() and path.name != "run_metadata.json"
    }
    metadata = {
        "task": "GSE248762 seven-validated-hub pseudobulk fold-change audit",
        "execution_mode": "metadata-only",
        "analysis_modes": [
            "earlier full edgeR and Monte Carlo audit outputs reused unchanged",
            "exact fixed-size permutation enumeration run separately",
            "per-sample DecontX correction and corrected-count contrasts run separately",
            "per-cell marker screen and pre-QC called-barcode check run separately",
            "reports, inventory, and checksums regenerated in metadata-only mode",
        ],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "gene_set_source": "results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv",
        "random_seed": SEED,
        "permutations_per_BH_pair": N_PERMUTATIONS,
        "permutation_note": "Earlier Monte Carlo outputs (200 draws) were reused unchanged for provenance only and are superseded by exact enumeration for the three ECM1 comparisons.",
        "exact_permutation": {
            "runner": str(EXACT_SCRIPT.relative_to(ROOT)),
            "statistic": "absolute edgeR robust quasi-likelihood log2FC; full edgeR model refit for every fixed-size allocation",
            "tail": "inclusive two-sided; count(|permuted log2FC| >= |observed log2FC|) / number of allocations; no plus-one correction",
            "seed": None,
            "workers": 8,
            "summaries": [
                str((OUT / f"sc7_exact_permutation_summary_{lineage}.csv").relative_to(ROOT))
                for lineage in EXACT_LINEAGES
            ],
            "all_allocations_verified": True,
            "execution_note": "The first serial attempt was stopped after checkpointing; the parallel runner resumed the same allocation-index files and verified complete, unique allocation coverage.",
        },
        "decontx": {
            "runner": str(DECONT_X_SCRIPT.relative_to(ROOT)),
            "package": "celda::decontX",
            "seed": SEED,
            "unit": "per sample",
            "cell_clusters": "saved cell_type annotations",
            "background_empty_droplets_supplied": False,
            "corrected_count_handling": "Rounded decontXcounts per package documentation, then summed by sample and lineage.",
            "contrasts": "LV vs SV; LV_UF vs SV; LV_NOT_UF vs SV",
            "BH_q_values": "NOT VERIFIED; no DecontX BH q-values are reported.",
            "ambient_profile_status": "NOT VERIFIED; no separate raw all-droplet/empty-droplet matrix was present among the local GSE248762 extracted files.",
        },
        "pre_qc_barcode_check": {
            "runner": str(PRE_QC_SCRIPT.relative_to(ROOT)),
            "source": "local GSE248762 feature-barcode matrices and saved QC-retained cell IDs",
            "scope": "barcodes removed by the saved QC filter; these are cell-called barcodes, not empty droplets, and have no saved lineage labels",
            "outputs": [
                "results/revision/sc7_fc_audit/sc7_ecm1_pre_qc_excluded_summary.csv",
                "results/revision/sc7_fc_audit/sc7_ecm1_pre_qc_excluded_barcodes.csv",
            ],
        },
        "ambient_marker_screen": {
            "runner": str(AMBIENT_SCRIPT.relative_to(ROOT)),
            "contrast": "LV_UF vs all other samples, sample-level myeloid CPM",
            "test": "one-sided exact label permutation, all fixed-size allocations",
            "n_label_allocations": 1820,
            "multiplicity": "BH across eight screened markers",
        },
        "edgeR_model": {
            "design": "~0 + group",
            "contrast": "LV - SV",
            "min_cells_per_sample": 20,
            "min_umis_per_sample_gene": 10,
            "min_samples_per_group": 2,
            "filter": "edgeR::filterByExpr; original testable hub genes explicitly retained",
            "dispersion": "edgeR::estimateDisp(robust=TRUE)",
            "fit": "edgeR::glmQLFit(robust=TRUE)",
        },
        "package_versions": package_versions,
        "files_read": [str(path.relative_to(ROOT)) for path in read_paths],
        "files_written": [
            str(path.relative_to(ROOT))
            for path in written_paths + [(OUT / "run_metadata.json").resolve()]
        ],
        "provenance": {
            "original_edgeR_and_200_draw_permutation_outputs": {
                "origin": "Earlier full run of run_fc_audit.py; reused unchanged. The Monte Carlo results are superseded by exact enumeration for the ECM1 comparisons.",
                "files": [f"results/revision/sc7_fc_audit/{name}" for name in R_OUTPUTS],
            },
            "exact_permutation_outputs": {
                "origin": "Generated by the approved exact enumeration; the parallel runner resumed allocation-index checkpoints.",
                "files": [
                    str((OUT / f"sc7_exact_permutation_results_{lineage}.csv").relative_to(ROOT))
                    for lineage in EXACT_LINEAGES
                ]
                + [
                    str((OUT / f"sc7_exact_permutation_summary_{lineage}.csv").relative_to(ROOT))
                    for lineage in EXACT_LINEAGES
                ],
            },
            "decontx_outputs": {
                "origin": "Generated by the approved per-sample DecontX run using local QC-retained count matrices.",
                "files": [
                    str(path.relative_to(ROOT))
                    for path in written_paths
                    if path.name.startswith("sc7_decontx_")
                    or path.name.startswith("decontx_")
                ],
            },
            "ambient_and_pre_qc_outputs": {
                "origin": "Generated from saved per-cell exports, local GSE248762 cell-called matrices, and saved QC cell IDs; no empty droplets were available.",
                "files": [
                    str(path.relative_to(ROOT))
                    for path in written_paths
                    if path.name.startswith("sc7_ecm1_")
                    and (
                        "ambient" in path.name
                        or "cell_expression" in path.name
                        or "pre_qc_excluded" in path.name
                    )
                ],
            },
            "report_and_inventory_outputs": [
                str(path.relative_to(ROOT))
                for path in written_paths + [(OUT / "run_metadata.json").resolve()]
                if path.name in {
                    "FC_AUDIT_REPORT.md",
                    "MANUSCRIPT_WRITEUP_HUB_CONTEXT_AND_NETWORK.md",
                    "FILE_INVENTORY.md",
                    "run_metadata.json",
                }
            ],
            "checksum_interpretation": (
                "SHA-256 values record current file contents and are not proof of which run "
                "originally generated reused files."
            ),
        },
        "sha256_note": "run_metadata.json is omitted from its own checksum map to avoid self-reference.",
        "sha256": hashes,
    }
    metadata_path = OUT / "run_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "audit_pairs": len(audit),
                "BH_pairs": len(audit[audit["saved_BH_q"] < 0.05]),
                "permutation_rows": len(pd.read_csv(OUT / "sc7_permutation_results.csv")),
                "output_dir": str(OUT),
                "files_read": len(read_paths),
                "files_written": len(written_paths) + 1,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
