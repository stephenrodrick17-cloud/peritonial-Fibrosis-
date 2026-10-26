"""Apply global BH correction, export result lists, and make contextualization figures."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


OUT = Path(__file__).resolve().parent
HUB_ORDER = [
    "TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19",
    "NAT6", "ECM1", "CST7", "CCBE1",
]
LINEAGE_ORDER = [
    "lymphoid", "myeloid", "fibroblast/stromal", "mesothelial",
    "endothelial", "epithelial", "other/unresolved",
]
LINEAGE_COLORS = {
    "lymphoid": "#2878A5",
    "myeloid": "#C84B31",
    "fibroblast/stromal": "#218C74",
    "mesothelial": "#D19A00",
    "endothelial": "#7651A8",
    "epithelial": "#C05A9B",
    "other/unresolved": "#777777",
}
GROUP_COLORS = {"SV": "#2878A5", "LV_NOT_UF": "#D19A00", "LV_UF": "#C84B31"}


def save_figure(fig: plt.Figure, stem: str) -> None:
    for extension in ("png", "pdf"):
        fig.savefig(
            OUT / f"{stem}.{extension}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def make_umap_figures(annotations: pd.DataFrame) -> None:
    fig, axis = plt.subplots(figsize=(8.3, 6.8))
    for lineage in LINEAGE_ORDER:
        subset = annotations[annotations["lineage"] == lineage]
        if subset.empty:
            continue
        axis.scatter(
            subset["UMAP1"], subset["UMAP2"], s=0.55,
            c=LINEAGE_COLORS[lineage], label=lineage, alpha=0.75, rasterized=True,
        )
    axis.set(title="GSE248762 cells colored by marker-defined lineage", xlabel="UMAP1", ylabel="UMAP2")
    axis.legend(markerscale=8, frameon=False, fontsize=8, loc="best")
    save_figure(fig, "umap_by_lineage")

    fig, axis = plt.subplots(figsize=(8.3, 6.8))
    group_order = ["SV", "LV_NOT_UF", "LV_UF"]
    for group in group_order:
        subset = annotations[annotations["original_group"] == group]
        axis.scatter(
            subset["UMAP1"], subset["UMAP2"], s=0.55,
            c=GROUP_COLORS[group], label=group, alpha=0.75, rasterized=True,
        )
    axis.set(title="GSE248762 cells colored by original GEO group", xlabel="UMAP1", ylabel="UMAP2")
    axis.legend(markerscale=8, frameon=False, fontsize=8, loc="best")
    save_figure(fig, "umap_by_group")


def make_dotplot(profile: pd.DataFrame) -> None:
    lineages = [lineage for lineage in LINEAGE_ORDER if lineage in set(profile["lineage"])]
    profile = profile.copy()
    profile["gene"] = pd.Categorical(profile["gene"], categories=HUB_ORDER, ordered=True)
    profile["lineage"] = pd.Categorical(profile["lineage"], categories=lineages, ordered=True)
    profile = profile.sort_values(["gene", "lineage"])
    observed = profile[profile["feature_detected"] & profile["cells_in_lineage"].gt(0)]
    max_percent = max(float(observed["percent_expressing_cells"].max()), 1.0)
    sizes = 20 + 260 * observed["percent_expressing_cells"].astype(float) / max_percent

    fig, axis = plt.subplots(figsize=(max(8, len(lineages) * 1.35), 6.3), constrained_layout=True)
    points = axis.scatter(
        observed["lineage"].cat.codes,
        observed["gene"].cat.codes,
        s=sizes,
        c=observed["mean_log1p_CPM10k"].astype(float),
        cmap="viridis",
        edgecolors="#333333",
        linewidths=0.25,
    )
    axis.set_xticks(range(len(lineages)), lineages, rotation=30, ha="right")
    axis.set_yticks(range(len(HUB_ORDER)), HUB_ORDER)
    axis.set_xlabel("Marker-defined lineage")
    axis.set_ylabel("Project consensus hub gene")
    axis.set_title("Hub-gene expression across lineages")
    fig.colorbar(points, ax=axis, label="Mean log1p normalized expression (CP10k)")
    for percent in (10, 50, 90):
        axis.scatter([], [], s=20 + 260 * percent / max_percent, c="#777777", label=f"{percent}% expressing")
    axis.legend(title="Expressing cells", frameon=False, loc="upper left", bbox_to_anchor=(1.22, 1))
    save_figure(fig, "hub_gene_lineage_dotplot")


def make_feature_plot(
    annotations: pd.DataFrame, per_cell: pd.DataFrame, feature_presence: pd.DataFrame
) -> None:
    merged = annotations[["cell_id", "UMAP1", "UMAP2"]].merge(
        per_cell[["cell_id", *HUB_ORDER]], on="cell_id", validate="one_to_one"
    )
    present = feature_presence.set_index("gene")["feature_detected"].to_dict()
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 12), constrained_layout=True)
    for axis, gene in zip(axes.flat, HUB_ORDER):
        values = pd.to_numeric(merged[gene], errors="coerce").to_numpy(dtype=float)
        if not present.get(gene, False) or not np.isfinite(values).any():
            axis.text(0.5, 0.5, "Feature not detected", ha="center", va="center", transform=axis.transAxes)
            axis.set_title(gene)
            axis.set_axis_off()
            continue
        vmax = float(np.nanpercentile(values, 99.5))
        vmax = max(vmax, 1e-8)
        points = axis.scatter(
            merged["UMAP1"], merged["UMAP2"], c=values,
            s=0.45, cmap="magma", vmin=0, vmax=vmax, rasterized=True,
        )
        axis.set_title(gene)
        axis.set_xlabel("UMAP1")
        axis.set_ylabel("UMAP2")
        fig.colorbar(points, ax=axis, fraction=0.046, pad=0.02, label="log1p normalized expression")
    fig.suptitle("Consensus hub-gene feature plots")
    save_figure(fig, "umap_hub_gene_feature_plots")


def create_paragraph(
    audit: dict, membership: pd.DataFrame, comparison: pd.DataFrame
) -> tuple[str, dict]:
    testable = comparison[comparison["status"] == "testable"]
    nominal = testable[testable["P_value"] < 0.05].copy()
    fdr_hits = testable[testable["BH_FDR_global"] < 0.05].copy()
    nominal_only = nominal[nominal["BH_FDR_global"] >= 0.05].copy()
    tested_genes = sorted(testable["gene"].unique())
    not_testable_genes = sorted(set(membership["gene"]) - set(tested_genes))
    n_tests = len(testable)
    min_fdr = float(testable["BH_FDR_global"].min()) if n_tests else float("nan")

    original = audit["original_group_counts"]
    contrast = audit["primary_contrast_sample_counts"]
    if not fdr_hits.empty:
        hit_text = "; ".join(
            f"{row.gene} in {row.lineage} (log2FC = {row.log2FC_LV_vs_SV:.3f}, "
            f"P = {row.P_value:.3g}, FDR = {row.BH_FDR_global:.3g})"
            for row in fdr_hits.itertuples(index=False)
        )
        inference = f"Contrasts passing global FDR < 0.05 were {hit_text}."
    else:
        inference = (
            f"No gene-lineage contrast remained significant after correction across "
            f"the full set of {n_tests} tested gene-lineage-contrast combinations "
            f"(minimum FDR = {min_fdr:.3g}); nominal findings, if present, are "
            "directional trends only."
        )
    if nominal_only.empty:
        trend_text = "No nominal-only P < 0.05 trends were observed."
    else:
        trends = []
        for row in nominal_only.itertuples(index=False):
            direction = "higher" if row.log2FC_LV_vs_SV > 0 else "lower"
            trends.append(
                f"{row.gene} was directionally {direction} in LV than SV in "
                f"{row.lineage} (log2FC = {row.log2FC_LV_vs_SV:.3f}, "
                f"P = {row.P_value:.3g}, FDR = {row.BH_FDR_global:.3g})"
            )
        trend_text = (
            "Nominal-only P < 0.05 directional trends (not FDR-significant) were: "
            + "; ".join(trends)
            + "."
        )

    paragraph = (
        f"We examined GSE248762, a 10x single-cell RNA-sequencing dataset of "
        f"overnight peritoneal dialysis effluent, containing {audit['geo_sample_records']} "
        f"samples (SV n = {original.get('SV', 0)}, LV_NOT_UF n = "
        f"{original.get('LV_NOT_UF', 0)}, LV_UF n = {original.get('LV_UF', 0)}; "
        f"the pseudobulk contrast combined the long-vintage groups as LV "
        f"(n = {contrast.get('LV', 0)}) versus SV (n = {contrast.get('SV', 0)}). "
        f"After applying the prespecified cell filter, {audit['cells_after_recomputed_qc']:,} "
        f"of {audit['cells_before_qc']:,} cells were retained. The project shared-code "
        f"consensus hub panel comprised {len(membership)} genes "
        f"({', '.join(membership['gene'])}); {int(membership['in_five_gene_model_panel'].sum())} "
        f"were members of the five-gene model panel and "
        f"{int(membership['in_seven_gene_consensus_panel'].sum())} were members of the "
        f"seven-gene consensus panel. {len(tested_genes)} of the {len(membership)} hub "
        f"genes ({', '.join(tested_genes)}) were testable in at least one lineage under "
        f"the prespecified detection filter; {', '.join(not_testable_genes)} was not "
        f"testable because it was absent from the expression feature table. Across the "
        f"{n_tests} testable gene-lineage-contrast combinations, {len(nominal)} had "
        f"nominal P < 0.05 and {len(fdr_hits)} passed global FDR < 0.05 (minimum FDR = "
        f"{min_fdr:.3g}). {trend_text} {inference} "
        f"These pseudobulk comparisons use samples as "
        f"the replicate unit and provide an exploratory cellular contextualization, "
        f"not validation of the diagnostic models; the small number of samples per "
        f"group limits power."
    )
    summary = {
        "n_samples": int(audit["geo_sample_records"]),
        "original_group_counts": original,
        "primary_contrast_sample_counts": contrast,
        "cells_before_qc": int(audit["cells_before_qc"]),
        "cells_after_qc": int(audit["cells_after_recomputed_qc"]),
        "hub_panel_size": int(len(membership)),
        "hub_genes_testable_in_any_lineage": tested_genes,
        "hub_genes_not_testable_in_any_lineage": not_testable_genes,
        "n_tested_gene_lineage_contrasts": int(n_tests),
        "n_nominal_p_lt_0_05": int(len(nominal)),
        "n_global_fdr_lt_0_05": int(len(fdr_hits)),
        "minimum_global_fdr": min_fdr if np.isfinite(min_fdr) else None,
        "nominal_trends": nominal[
            ["gene", "lineage", "log2FC_LV_vs_SV", "P_value", "BH_FDR_global"]
        ].to_dict(orient="records"),
        "fdr_hits": fdr_hits[
            ["gene", "lineage", "log2FC_LV_vs_SV", "P_value", "BH_FDR_global"]
        ].to_dict(orient="records"),
    }
    return paragraph, summary


def main() -> None:
    audit = json.loads((OUT / "dataset_audit.json").read_text(encoding="utf-8"))
    annotations = pd.read_csv(OUT / "cell_type_annotation_table.csv")
    profile = pd.read_csv(OUT / "hub_gene_expression_by_lineage.csv")
    per_cell = pd.read_csv(OUT / "hub_gene_expression_per_cell.csv")
    membership = pd.read_csv(OUT / "hub_gene_panel_membership.csv")
    feature_presence = pd.read_csv(OUT / "hub_gene_feature_presence.csv")
    raw_results = pd.read_csv(OUT / "edger_gene_lineage_results_raw.csv")

    testable_mask = raw_results["status"] == "testable"
    raw_results["BH_FDR_global"] = np.nan
    if testable_mask.any():
        raw_results.loc[testable_mask, "BH_FDR_global"] = multipletests(
            raw_results.loc[testable_mask, "P_value"].to_numpy(dtype=float),
            method="fdr_bh",
        )[1]
    raw_results["nominal_P_lt_0_05"] = (
        testable_mask & raw_results["P_value"].lt(0.05)
    )
    raw_results["global_FDR_lt_0_05"] = (
        testable_mask & raw_results["BH_FDR_global"].lt(0.05)
    )
    raw_results.to_csv(OUT / "gene_lineage_pseudobulk_tests.csv", index=False)
    raw_results.loc[raw_results["nominal_P_lt_0_05"]].to_csv(
        OUT / "nominal_P_lt_0_05.csv", index=False
    )
    raw_results.loc[raw_results["global_FDR_lt_0_05"]].to_csv(
        OUT / "global_FDR_lt_0_05.csv", index=False
    )
    raw_results[[
        "gene", "lineage", "status", "not_testable_reason",
        "n_LV_samples_with_ge_20_cells", "n_SV_samples_with_ge_20_cells",
        "n_LV_samples_with_ge_10_umis", "n_SV_samples_with_ge_10_umis",
    ]].to_csv(OUT / "hub_gene_testability_by_lineage.csv", index=False)
    overall_testability = raw_results.groupby("gene", sort=False).agg(
        testable_in_any_lineage=("status", lambda values: bool((values == "testable").any())),
        testable_lineage_count=("status", lambda values: int((values == "testable").sum())),
        feature_absent_everywhere=(
            "not_testable_reason",
            lambda values: bool(
                len(values) > 0
                and all(value == "gene symbol absent from feature matrix" for value in values)
            ),
        ),
    ).reset_index()
    overall_testability.to_csv(OUT / "hub_gene_testability_overall.csv", index=False)

    make_umap_figures(annotations)
    make_dotplot(profile)
    make_feature_plot(annotations, per_cell, feature_presence)

    paragraph, summary = create_paragraph(audit, membership, raw_results)
    (OUT / "results_paragraph.txt").write_text(paragraph + "\n", encoding="utf-8")
    (OUT / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8"
    )
    print(paragraph)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
