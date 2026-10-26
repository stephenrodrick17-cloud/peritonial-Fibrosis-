"""Create publication figures from saved analysis outputs only.

No differential expression, WGCNA, enrichment, ROC, immune scoring, or
statistical testing is performed here.
"""

from __future__ import annotations

import hashlib
import io
import json
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "revision" / "figures_final"
SEED = 42
PANEL_V2 = ["FLT3LG", "EBI3", "LTB", "ADAM19", "CST7"]
ORIGINAL_PANEL = [
    "TNFSF15",
    "FLT3LG",
    "EBI3",
    "LTB",
    "ADAM19",
    "ECM1",
    "SERPINA10",
    "CST7",
]
COHORTS = ["GSE248762", "GSE130888", "GSE62928", "GSE92455"]
COLORS = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "purple": "#CC79A7",
    "yellow": "#E69F00",
    "grey": "#999999",
    "light_grey": "#E5E5E5",
    "dark": "#222222",
}
INPUT_HASHES: dict[str, str] = {}
LEGENDS: list[str] = []
MANIFEST: list[dict[str, str]] = []
SKIPPED: list[tuple[str, str]] = []


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def register_input(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"Required saved input is missing: {path}")
    relative = path.relative_to(ROOT).as_posix()
    if relative not in INPUT_HASHES:
        INPUT_HASHES[relative] = sha256(path)
    return relative


def read_csv(relative: str, **kwargs: Any) -> pd.DataFrame:
    path = ROOT / relative
    register_input(path)
    return pd.read_csv(path, **kwargs)


def read_json(relative: str) -> dict[str, Any]:
    path = ROOT / relative
    register_input(path)
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def source_csv(name: str, data: pd.DataFrame) -> str:
    path = OUT / name
    buffer = io.StringIO()
    data.to_csv(buffer, index=False)
    content = buffer.getvalue()
    path.write_text(content, encoding="utf-8", newline="")
    return path.name


def save_figure(fig: plt.Figure, stem: str, inputs: list[str]) -> list[str]:
    source_name = {
        "Fig1_study_workflow": "fig1_source_data.csv",
        "Fig2_discovery_volcano": "fig2_source_data.csv",
        "Fig3_overlap": "fig3_source_data.csv",
        "Fig4_three_model_union": "fig4_source_data.csv",
        "Fig5_hallmark_gsea": "fig5_source_data.csv",
        "Fig6_panel_cohorts": "fig6_source_data.csv",
        "Fig8_immune_context": "fig8_source_data.csv",
        "Fig9_historical_ppi_networks": "fig9_source_data.csv",
        "Fig10_single_cell_context": "fig10_source_data.csv",
        "FigS3_detectability_rules": "figS3_source_data.csv",
    }[stem]
    final_stem = f"{stem}_reviewed15"
    outputs = [f"{final_stem}.png", f"{final_stem}.pdf", source_name]
    for filename in outputs[:2]:
        if (OUT / filename).exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {OUT / filename}")
    fig.savefig(OUT / outputs[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / outputs[1], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    add_manifest(final_stem, inputs, outputs)
    return outputs


def add_manifest(figure: str, inputs: list[str], outputs: list[str]) -> None:
    for source in sorted(set(inputs)):
        MANIFEST.append(
            {
                "figure": figure,
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": "",
                "input_file": source,
                "input_sha256": INPUT_HASHES.get(source, ""),
                "output_files": ";".join(outputs),
                "status": "created",
            }
        )


def add_legend(text: str) -> None:
    LEGENDS.append(text)


def mm(value: float) -> float:
    return value / 25.4


def set_panel_title(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, loc="left", fontweight="bold", pad=5)


def flow_figure() -> None:
    files = [
        "results/tables/GSE125498_all_results.csv",
        "results/tables/GSE125498_sample_metadata.csv",
        "results/final_analysis/FINAL_DEG_Matrisomal_Genes.csv",
        "results/tables/GSE125498_ECM_intersection.csv",
        "results/revision/wgcna_modules_dynamic_taskB.csv",
        "results/revision/ml_624_intersection_results.csv",
        "results/revision/ml_38_ecm_wgcna_results.csv",
        "results/revision/panel_v2/panel_v2_rule_results.csv",
        "results/revision/panel_v2/panel_v2.json",
        "results/revision/effluent_validation/GSE248762/analysis_20261009/patient_celltype_composition.csv",
        "results/revision/effluent_validation/GSE130888/analysis_20261009/patient_celltype_composition.csv",
        "results/revision/gse62928_validation_8genes/validation_summary.json",
        "results/revision/panel_v2/GSE92455_20261011/panel_v2_score_summary.json",
    ]
    discovery = read_csv(files[0])
    metadata = read_csv(files[1])
    matrisome_degs = read_csv(files[2])
    measured_matrisome = read_csv(files[3])
    wgcna = read_csv(files[4])
    broad_pool = read_csv(files[5])
    ecm_pool = read_csv(files[6])
    rules = read_csv(files[7])
    panel_v2 = read_json(files[8])
    cohort_248 = read_csv(files[9])
    cohort_130 = read_csv(files[10])
    cohort_629 = read_json(files[11])
    cohort_924 = read_json(files[12])

    group_n = metadata["group"].value_counts().to_dict()
    status = discovery["DEG_Status"].value_counts().to_dict()
    non_grey_n = int(wgcna.loc[wgcna["Module"].astype(str).str.lower() != "grey", "Gene"].nunique())
    external_n = {
        "GSE62928": int(cohort_629["qc"]["sample_count"]),
        "GSE248762": int(len(cohort_248)),
        "GSE130888": int(len(cohort_130)),
    }
    labels = [
        (
            "Discovery",
            f"GSE125498\nLPD n={int(group_n['LPD'])}; SPD n={int(group_n['SPD'])}\n"
            f"{len(discovery):,} genes",
        ),
        (
            "DEG screen",
            f"|log2FC| >= 0.585; P < 0.05\n"
            f"{int(status.get('Upregulated', 0) + status.get('Downregulated', 0)):,} DEGs\n"
            f"{int(status.get('Upregulated', 0))} up / {int(status.get('Downregulated', 0))} down",
        ),
        (
            "Matrisome overlap",
            f"{len(measured_matrisome):,} measured\n{len(matrisome_degs)} ECM DEGs\n"
            f"{int((matrisome_degs['logFC'] > 0).sum())} up / {int((matrisome_degs['logFC'] < 0).sum())} down",
        ),
        (
            "WGCNA + 3-model ML",
            f"{non_grey_n:,} non-grey module genes\n"
            f"{len(broad_pool)} broad / {len(ecm_pool)} ECM pool\nLASSO, SVM-RFE, RF",
        ),
        ("Original fixed panel", f"{rules['gene'].nunique()} genes\nfixed 8-gene panel"),
        (
            "External cohorts",
            "GSE62928 n=8\nGSE248762 n=16\nGSE130888 n=13",
        ),
        (
            "Detectability rules",
            "8 genes across 3 cohorts\n(a) detected in >=2\n(b) direction concordant\n(c) BH q is a flag",
        ),
        (
            "Frozen panel v2",
            f"{len(panel_v2['genes_passing_selection_rules'])} genes "
            f"({panel_v2['date_defined']})\n"
            + ", ".join(panel_v2["genes_passing_selection_rules"]),
        ),
        (
            "GSE92455 test",
            f"Epithelioid n={cohort_924['n_epithelioid_patients']}\n"
            f"Non-epithelioid n={cohort_924['n_non_epithelioid_patients']}\nnot a vintage test",
        ),
        ("Exploratory context", "GSE125498; n=33\nsame data as discovery\nnot validation"),
    ]
    title_text = [
        "Discovery",
        "DEG screen",
        "Matrisome intersection",
        "WGCNA + 3-model ML",
        "Original fixed panel",
        "External cohorts",
        "Detectability rules",
        "Panel v2 freeze",
        "GSE92455 test",
        "Exploratory context",
    ]
    body_text = [
        f"LPD n={int(group_n['LPD'])} | SPD n={int(group_n['SPD'])}\n{len(discovery):,} genes",
        f"|log2FC| >= 0.585; P < 0.05\n"
        f"{int(status.get('Upregulated', 0) + status.get('Downregulated', 0))} DEGs\n"
        f"{int(status.get('Upregulated', 0))} up; {int(status.get('Downregulated', 0))} down",
        f"{len(measured_matrisome):,} measured\n{len(matrisome_degs)} ECM DEGs\n"
        f"{int((matrisome_degs['logFC'] > 0).sum())} up; "
        f"{int((matrisome_degs['logFC'] < 0).sum())} down",
        f"{non_grey_n:,} non-grey\n{len(broad_pool)} broad; {len(ecm_pool)} ECM\nLASSO / SVM-RFE / RF",
        f"{rules['gene'].nunique()}-gene fixed panel",
        f"GSE62928 n={external_n['GSE62928']}; GSE248762 n={external_n['GSE248762']}\n"
        f"GSE130888 n={external_n['GSE130888']}",
        f"{rules['gene'].nunique()} genes / 3 cohorts\n>=2 detected; direction concordant\nBH q = flag only",
        f"Frozen {panel_v2['date_defined']}\nFLT3LG, EBI3, LTB\nADAM19, CST7",
        f"Epithelioid n={cohort_924['n_epithelioid_patients']}; "
        f"non-epithelioid n={cohort_924['n_non_epithelioid_patients']}\nNot a vintage test",
        f"GSE125498 n={int(group_n['LPD']) + int(group_n['SPD'])}\n"
        "Same data as discovery\nNot validation",
    ]
    source = pd.DataFrame(
        [
            {
                "step": name,
                "title_text": title_text[i],
                "box_text": body_text[i],
                "order": i + 1,
            }
            for i, (name, _) in enumerate(labels)
        ]
    )
    source_name = source_csv("fig1_source_data.csv", source)

    fig, ax = plt.subplots(figsize=(mm(180), mm(190)))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    positions = [
        (0.27, 0.89),
        (0.73, 0.89),
        (0.73, 0.70),
        (0.27, 0.70),
        (0.27, 0.51),
        (0.73, 0.51),
        (0.73, 0.32),
        (0.27, 0.32),
        (0.27, 0.13),
        (0.73, 0.13),
    ]
    box_w, box_h = 0.40, 0.15
    for i, ((name, text), (x, y)) in enumerate(zip(labels, positions)):
        color = COLORS["green"] if i in (7, 8, 9) else COLORS["blue"]
        rect = plt.Rectangle(
            (x - box_w / 2, y - box_h / 2),
            box_w,
            box_h,
            linewidth=1.0,
            edgecolor=color,
            facecolor="white",
            zorder=2,
        )
        ax.add_patch(rect)
        ax.text(x, y + 0.040, title_text[i], ha="center", va="center", fontsize=8.0, fontweight="bold", linespacing=1.0)
        ax.text(x, y - 0.025, body_text[i], ha="center", va="center", fontsize=8.0, linespacing=1.08)
    for i in range(len(positions) - 1):
        x, y = positions[i]
        next_x, next_y = positions[i + 1]
        if np.isclose(y, next_y):
            direction = 1 if next_x > x else -1
            xy = (next_x - direction * box_w / 2, y)
            xytext = (x + direction * box_w / 2, y)
        else:
            xy = (next_x, next_y + box_h / 2)
            xytext = (x, y - box_h / 2)
        ax.annotate(
            "",
            xy=xy,
            xytext=xytext,
            arrowprops={"arrowstyle": "->", "lw": 1, "color": COLORS["dark"]},
        )
    ax.set_title(
        "Project workflow — saved-analysis overview; post hoc context is labeled",
        loc="left",
        fontsize=9,
        fontweight="bold",
    )
    outputs = save_figure(fig, "Fig1_study_workflow", files)
    add_legend(
        f"### Figure 1. Study workflow\n"
        f"Project flow from GSE125498 discovery through the fixed panel and its detectability-based "
        f"v2 definition, followed by the GSE92455 phenotype test and exploratory discovery-context work. "
        f"Discovery n={int(group_n['LPD'])} LPD and n={int(group_n['SPD'])} SPD; external sample sizes "
        f"are printed in the boxes. Counts are taken from the saved tables named in the manifest. "
        f"Exploratory context is the same data as discovery and is not validation. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def volcano_figure() -> None:
    files = [
        "results/tables/GSE125498_all_results.csv",
        "results/final_analysis/FINAL_DEG_Matrisomal_Genes.csv",
        "results/revision/panel_v2/panel_v2.json",
        "results/tables/GSE125498_sample_metadata.csv",
    ]
    data = read_csv(files[0])
    matrisome = read_csv(files[1])
    panel = read_json(files[2])["genes_passing_selection_rules"]
    metadata = read_csv(files[3])
    q_genes = set(matrisome.loc[matrisome["adj.P.Val"] < 0.05, "Gene.symbol"].dropna().astype(str))
    panel_symbols = set(panel)
    data["minus_log10_P"] = -np.log10(data["P.Value"].clip(lower=1e-300))
    data["highlight_type"] = ""
    data.loc[data["Gene.symbol"].isin(q_genes), "highlight_type"] = "matrisome_BH_q_lt_0.05"
    data.loc[data["Gene.symbol"].isin(panel_symbols), "highlight_type"] = data.loc[
        data["Gene.symbol"].isin(panel_symbols), "highlight_type"
    ].map(lambda item: item + ";panel_v2" if item else "panel_v2")
    data["plotted_color_group"] = data["DEG_Status"].map(
        {"Upregulated": "Up DEG", "Downregulated": "Down DEG"}
    ).fillna("Not significant")
    source = data[
        [
            "Gene.symbol",
            "logFC",
            "P.Value",
            "minus_log10_P",
            "DEG_Status",
            "plotted_color_group",
            "highlight_type",
        ]
    ].copy()
    source_name = source_csv("fig2_source_data.csv", source)
    fig, ax = plt.subplots(figsize=(mm(180), mm(125)))
    colors = {"Not significant": COLORS["light_grey"], "Up DEG": COLORS["blue"], "Down DEG": COLORS["orange"]}
    for group in ["Not significant", "Up DEG", "Down DEG"]:
        subset = source[source["plotted_color_group"] == group]
        ax.scatter(
            subset["logFC"],
            subset["minus_log10_P"],
            s=7 if group == "Not significant" else 10,
            c=colors[group],
            alpha=0.55 if group == "Not significant" else 0.85,
            linewidths=0,
            label=f"{group} (n={len(subset)})",
            rasterized=True,
        )
    ax.axvline(0.585, color=COLORS["dark"], linestyle="--", linewidth=0.7)
    ax.axvline(-0.585, color=COLORS["dark"], linestyle="--", linewidth=0.7)
    ax.axhline(-np.log10(0.05), color=COLORS["dark"], linestyle="--", linewidth=0.7)
    label_offsets = {
        "ADAM19": (10, 10),
        "CST7": (10, -12),
        "EBI3": (-8, -12),
        "FLT3LG": (10, 16),
        "LTB": (8, 0),
        "SEMA3E": (-48, -12),
        "TNFSF15": (-48, 12),
    }
    for symbol in sorted(q_genes | panel_symbols):
        points = source[source["Gene.symbol"] == symbol]
        if points.empty:
            continue
        row = points.iloc[0]
        ax.scatter(row["logFC"], row["minus_log10_P"], s=24, facecolors="none", edgecolors=COLORS["purple"], linewidths=0.9)
        ax.annotate(
            symbol,
            (row["logFC"], row["minus_log10_P"]),
            xytext=label_offsets.get(symbol, (3, 3)),
            textcoords="offset points",
            ha="left" if label_offsets.get(symbol, (3, 3))[0] >= 0 else "right",
            arrowprops={"arrowstyle": "-", "color": COLORS["dark"], "linewidth": 0.45},
            fontsize=8,
        )
    ax.set_xlabel("log2 fold change (LPD relative to SPD)")
    ax.set_ylabel("-log10 nominal P")
    ax.set_title(
        f"GSE125498 volcano — LPD n={sum(metadata['group'] == 'LPD')}, SPD n={sum(metadata['group'] == 'SPD')}; "
        "exploratory, same data as discovery; not validation",
        loc="left",
    )
    ax.legend(
        frameon=False,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.2),
    )
    fig.subplots_adjust(bottom=0.25)
    outputs = save_figure(fig, "Fig2_discovery_volcano", files)
    add_legend(
        f"### Figure 2. GSE125498 volcano plot\n"
        f"All {len(data):,} saved discovery genes are plotted with x=log2FC and y=-log10 nominal P. "
        f"Dashed thresholds are |log2FC| >= 0.585 and P < 0.05. Up/down colors follow the saved DEG_Status "
        f"field; no FDR cutoff was used for the primary screen. Matrisome genes passing BH q<0.05 and the "
        f"five panel-v2 genes are outlined and labeled. n={int(sum(metadata['group'] == 'LPD'))} LPD and "
        f"n={int(sum(metadata['group'] == 'SPD'))} SPD. Circular: same data as discovery; not validation. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def draw_upset_panel(
    fig: plt.Figure,
    slot: Any,
    patterns: list[str],
    counts: list[float],
    set_names: list[str],
    color: str,
    title: str,
    show_set_names: bool = True,
) -> None:
    inner = GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=slot,
        height_ratios=[3.2, 1.9],
        hspace=0.03,
    )
    bars = fig.add_subplot(inner[0])
    matrix = fig.add_subplot(inner[1], sharex=bars)
    x = np.arange(len(patterns))
    maximum = max(counts, default=1)
    bars.bar(x, counts, width=0.72, color=color, zorder=1)
    bars.set_xlim(-0.55, len(patterns) - 0.45)
    bars.set_ylim(0, maximum * 1.30)
    bars.set_ylabel("Genes", fontsize=8)
    bars.set_title(title, loc="left", fontsize=8, fontweight="bold")
    bars.tick_params(axis="x", bottom=False, labelbottom=False)
    bars.tick_params(axis="y", labelsize=8)
    for position, count in zip(x, counts):
        bars.annotate(
            f"{int(count)}",
            (position, count),
            xytext=(0, 2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    matrix_y = list(reversed(range(len(set_names))))
    for row, pattern in enumerate(patterns):
        selected = [index for index, member in enumerate(pattern) if member == "1"]
        if len(selected) > 1:
            selected_y = [matrix_y[index] for index in selected]
            matrix.vlines(
                row,
                min(selected_y),
                max(selected_y),
                color=COLORS["dark"],
                linewidth=0.65,
                zorder=1,
            )
        for index, row_y in enumerate(matrix_y):
            is_member = pattern[index] == "1"
            matrix.scatter(
                row,
                row_y,
                s=8 if is_member else 5,
                color=COLORS["dark"] if is_member else COLORS["grey"],
                zorder=2,
            )
    matrix.set_ylim(-0.5, len(set_names) - 0.5)
    matrix.set_yticks(matrix_y, set_names, fontsize=8 if show_set_names else 0)
    if show_set_names:
        matrix.tick_params(axis="y", length=0)
    else:
        matrix.set_yticklabels([])
        matrix.spines["left"].set_visible(False)
    matrix.tick_params(axis="x", bottom=False, labelbottom=False)
    matrix.spines[["top", "right", "bottom"]].set_visible(False)
    bars.spines[["top", "right"]].set_visible(False)


def overlap_figure() -> None:
    files = [
        "results/tables/GSE125498_all_results.csv",
        "results/tables/GSE125498_ECM_intersection.csv",
        "results/final_analysis/FINAL_DEG_Matrisomal_Genes.csv",
        "results/revision/wgcna_modules_dynamic_taskB.csv",
        "results/tables/phase3_counts.json",
    ]
    discovery = read_csv(files[0])
    measured = read_csv(files[1])
    matrisome_degs = read_csv(files[2])
    modules = read_csv(files[3])
    overlap_metrics = read_json(files[4])
    deg_set = set(
        discovery.loc[discovery["DEG_Status"].isin(["Upregulated", "Downregulated"]), "Gene.symbol"]
        .dropna()
        .astype(str)
    )
    measured_set = set(measured["Gene.symbol"].dropna().astype(str))
    wgcna_set = set(
        modules.loc[modules["Module"].astype(str).str.lower() != "grey", "Gene"].dropna().astype(str)
    )
    universe = deg_set | measured_set | wgcna_set
    intersections = []
    for mask in range(1, 8):
        members = [
            name
            for bit, name in [(1, "DEG"), (2, "Measured matrisome"), (4, "Non-grey WGCNA")]
            if mask & bit
        ]
        count = sum(
            (gene in deg_set) == bool(mask & 1)
            and (gene in measured_set) == bool(mask & 2)
            and (gene in wgcna_set) == bool(mask & 4)
            for gene in universe
        )
        pattern = "".join("1" if mask & bit else "0" for bit in (1, 2, 4))
        intersections.append(
            {
                "intersection": " & ".join(members),
                "membership_code": pattern,
                "gene_count": count,
                "plot": "a",
                "record_type": "intersection_count",
            }
        )
    up_n = int((discovery["DEG_Status"] == "Upregulated").sum())
    universe_n = int(discovery["Gene.symbol"].nunique())
    measured_n = int(measured["Gene.symbol"].nunique())
    observed_up_matrisome = int(
        ((matrisome_degs["logFC"] > 0) & matrisome_degs["Gene.symbol"].notna()).sum()
    )
    expected = up_n * measured_n / universe_n
    fold = float(overlap_metrics["fold_enrichment_up_ecm"])
    hyper_p = float(overlap_metrics["hypergeom_p_up_ecm"])
    permutation_p = float(overlap_metrics["permutation_ecm_p_10k"])
    overlap_source = pd.DataFrame(
        intersections
        + [
            {
                "intersection": "Observed up-matrisome overlap",
                "membership_code": np.nan,
                "gene_count": observed_up_matrisome,
                "plot": "b",
                "record_type": "enrichment_statistic",
            },
            {
                "intersection": "Expected overlap",
                "membership_code": np.nan,
                "gene_count": expected,
                "plot": "b",
                "record_type": "enrichment_statistic",
            },
            {
                "intersection": "Fold enrichment",
                "membership_code": np.nan,
                "gene_count": fold,
                "plot": "b",
                "record_type": "enrichment_statistic",
            },
            {
                "intersection": "Hypergeometric P",
                "membership_code": np.nan,
                "gene_count": hyper_p,
                "plot": "b",
                "record_type": "enrichment_statistic",
            },
            {
                "intersection": "Permutation empirical P",
                "membership_code": np.nan,
                "gene_count": permutation_p,
                "plot": "b",
                "record_type": "enrichment_statistic",
            },
        ]
    )
    set_size_rows = [
        {
            "intersection": f"{name} set size",
            "membership_code": "",
            "gene_count": size,
            "plot": "a",
            "record_type": "set_size",
        }
        for name, size in [
            ("DEG", len(deg_set)),
            ("Measured matrisome", len(measured_set)),
            ("Non-grey WGCNA", len(wgcna_set)),
        ]
    ]
    overlap_source = pd.concat(
        [pd.DataFrame(intersections + set_size_rows), overlap_source.iloc[len(intersections):]],
        ignore_index=True,
    )
    source_name = source_csv("fig3_source_data.csv", overlap_source)
    fig = plt.figure(figsize=(mm(180), mm(130)))
    grid = GridSpec(1, 2, width_ratios=[1.35, 1], wspace=0.38, figure=fig)
    upset = pd.DataFrame(intersections).sort_values("gene_count", ascending=False)
    draw_upset_panel(
        fig,
        grid[0, 0],
        upset["membership_code"].tolist(),
        upset["gene_count"].astype(float).tolist(),
        [
            f"DEG ({len(deg_set):,})",
            f"Matrisome ({len(measured_set):,})",
            f"WGCNA ({len(wgcna_set):,})",
        ],
        COLORS["blue"],
        "a  DEG / matrisome / WGCNA",
    )

    right_grid = GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=grid[0, 1],
        height_ratios=[3, 1],
        hspace=0.75,
    )
    ax2 = fig.add_subplot(right_grid[0])
    stats_ax = fig.add_subplot(right_grid[1])
    bars = [observed_up_matrisome, expected]
    ax2.bar(["Observed\nup matrisome", "Expected"], bars, color=[COLORS["orange"], COLORS["grey"]], width=0.62)
    ax2.set_ylabel("Genes")
    ax2.set_ylim(0, max(bars) * 1.55)
    ax2.text(0, observed_up_matrisome + 0.7, f"{observed_up_matrisome}", ha="center", fontsize=8)
    ax2.text(1, expected + 0.7, f"{expected:.2f}", ha="center", fontsize=8)
    stats_ax.axis("off")
    stats_ax.text(
        0.02,
        0.98,
        f"Fold enrichment = {fold:.4f}\n"
        f"Hypergeometric P = {hyper_p:.7g}\n"
        f"Permutation empirical P = {permutation_p:.7g}",
        transform=stats_ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        linespacing=1.3,
    )
    set_panel_title(ax2, "b  Up-matrisome enrichment")
    fig.suptitle(
        f"GSE125498 discovery, n=13 LPD / 20 SPD; exploratory, same data as discovery; not validation",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.80, bottom=0.16, left=0.33, right=0.98)
    outputs = save_figure(fig, "Fig3_overlap", files)
    add_legend(
        f"### Figure 3. Discovery overlap and up-matrisome enrichment\n"
        f"(a) UpSet intersections among {len(deg_set):,} discovery DEGs, {len(measured_set):,} measured "
        f"Matrisome genes, and {len(wgcna_set):,} non-grey WGCNA module genes. (b) The saved "
        f"{observed_up_matrisome} observed upregulated matrisome DEGs are compared "
        f"with the expectation from the saved measured-gene universe ({universe_n:,} genes), {up_n} up-DEGs, and "
        f"{measured_n} measured Matrisome genes. The exact expectation from these saved inputs is "
        f"{expected:.4f} (displayed as {expected:.2f}); this differs from the requested rounded 18.09. "
        f"Fold enrichment={fold:.4f}; one-sided hypergeometric P={hyper_p:.7g}; "
        f"10,000-draw permutation empirical P={permutation_p:.7g}. These P values use different null models. "
        f"n=13 LPD and n=20 SPD. Exploratory; same data as discovery; not validation. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def ml_figure() -> None:
    files = [
        "results/revision/ml_624_intersection_results.csv",
        "results/revision/ml_38_ecm_wgcna_results.csv",
        "results/revision/wgcna_modules_dynamic_taskB.csv",
        "results/tables/GSE125498_sample_metadata.csv",
    ]
    data_frames = []
    for pool, path in [("624-gene pool", files[0]), ("38-gene ECM/WGCNA pool", files[1])]:
        frame = read_csv(path).copy()
        frame["pool"] = pool
        frame["pattern"] = frame.apply(
            lambda row: "".join(
                str(int(bool(row[column])))
                for column in ["LASSO_selected", "SVM_RFE_selected", "RF_selected"]
            ),
            axis=1,
        )
        frame["vote_count"] = frame["pattern"].map(lambda item: item.count("1"))
        frame["majority_vote"] = frame["vote_count"] >= 2
        frame["three_of_three"] = frame["vote_count"] == 3
        data_frames.append(frame)
    all_votes = pd.concat(data_frames, ignore_index=True)
    metadata = read_csv(files[3])
    summary = (
        all_votes.groupby(["pool", "pattern"], as_index=False)
        .agg(gene_count=("Gene", "nunique"))
        .sort_values(["pool", "pattern"])
    )
    summary["record_type"] = "intersection_count"
    majority = all_votes.loc[all_votes["majority_vote"]].copy()
    majority["record_type"] = "majority_gene"
    vote_table = majority[
        ["pool", "pattern", "Gene", "vote_count", "three_of_three"]
    ].rename(columns={"Gene": "gene"})
    summary = summary[["pool", "pattern", "gene_count"]]
    source = pd.concat(
        [
            summary.assign(record_type="intersection_count"),
            vote_table.assign(record_type="majority_gene", gene_count=np.nan),
        ],
        ignore_index=True,
        sort=False,
    )
    source_name = source_csv("fig4_source_data.csv", source)
    fig = plt.figure(figsize=(mm(180), mm(170)))
    grid = GridSpec(1, 3, width_ratios=[1.2, 1.2, 1.5], wspace=0.38, figure=fig)
    for index, (pool, color, title) in enumerate(
        [
            ("624-gene pool", COLORS["blue"], "624-gene pool"),
            ("38-gene ECM/WGCNA pool", COLORS["orange"], "38-gene ECM/WGCNA pool"),
        ]
    ):
        patterns = [f"{i:03b}" for i in range(8)]
        counts = (
            summary.loc[summary["pool"] == pool]
            .set_index("pattern")["gene_count"]
            .reindex(patterns, fill_value=0)
        )
        draw_upset_panel(
            fig,
            grid[0, index],
            patterns,
            counts.astype(float).tolist(),
            ["LASSO", "SVM-RFE", "RF"],
            color,
            title,
            show_set_names=index == 0,
        )
    ax2 = fig.add_subplot(grid[0, 2])
    majority_display = majority.sort_values(["pool", "vote_count", "Gene"], ascending=[True, False, True])
    cell_text = [
        [
            {"624-gene pool": "624", "38-gene ECM/WGCNA pool": "38"}[row.pool],
            row.Gene,
            f"{int(row.vote_count)}/3",
            "yes" if row.three_of_three else "",
        ]
        for row in majority_display.itertuples()
    ]
    table = ax2.table(
        cellText=cell_text,
        colLabels=["Pool", "Gene", "Votes", "3/3"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.16, 0.48, 0.17, 0.13],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.18)
    ax2.axis("off")
    set_panel_title(ax2, "Majority-vote genes")
    gse_n = len(metadata)
    fig.suptitle(
        f"GSE125498 (n={gse_n}; 13 LPD, 20 SPD); three models only — no XGBoost\n"
        "Exploratory, same data as discovery; not validation",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.84, bottom=0.14, left=0.14, right=0.98)
    outputs = save_figure(fig, "Fig4_three_model_union", files)
    three_genes = majority.loc[
        (majority["pool"] == "38-gene ECM/WGCNA pool") & majority["three_of_three"], "Gene"
    ].tolist()
    add_legend(
        f"### Figure 4. Three-model machine-learning union\n"
        f"UpSet plots show all selection intersections and aligned membership matrices for the 624-gene "
        f"broad pool and 38-gene ECM/WGCNA pool using LASSO, SVM-RFE, and random forest only; XGBoost "
        f"is excluded. The adjacent "
        f"table lists genes with >=2/3 votes. In the 38-gene pool, the saved 3/3 list is "
        f"{', '.join(three_genes)}. GSE125498 n=13 LPD and n=20 SPD. WGCNA was not recomputed inside "
        f"the CV folds. Exploratory; same data as discovery; not validation. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def normalize_pathway(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).lower()).strip()


def gsea_figure() -> None:
    files = [
        "results/revision/panel_v2_context/hallmark_preranked_gsea.csv",
        "results/tables/GSE125498_sample_metadata.csv",
        "results/revision/clusterprofiler_GO_BP_all_taskB.csv",
        "scripts/revision_wgcna_clusterprofiler.R",
    ]
    gsea = read_csv(files[0])
    metadata = read_csv(files[1])
    read_csv(files[2])
    register_input(ROOT / files[3])
    pathways = gsea["pathway"].astype(str).map(normalize_pathway)
    targets = {
        "Epithelial Mesenchymal Transition": "epithelial mesenchymal transition",
        "TGF Beta Signaling": "tgf beta signaling",
        "Angiogenesis": "angiogenesis",
        "Coagulation": "coagulation",
        "Apical Junction": "apical junction",
    }
    target_rows = []
    for label, token in targets.items():
        found = gsea.loc[pathways.str.contains(token, regex=False)].copy()
        for _, row in found.iterrows():
            target_rows.append((label, row))
    target_pathways = {str(row["pathway"]) for _, row in target_rows}
    significant = gsea.loc[gsea["GSEA_FDR_q"] < 0.05].copy()
    plotted = pd.concat([significant, gsea[gsea["pathway"].isin(target_pathways)]], ignore_index=True)
    plotted = plotted.drop_duplicates("pathway").copy()
    plotted["plot_role"] = np.where(
        plotted["GSEA_FDR_q"] < 0.05, "BH FDR q < 0.05", "requested pathway; not significant"
    )
    plotted["panel_v2_leading_edge_genes"] = plotted["panel_v2_genes_in_leading_edge"].fillna("")
    plotted["minus_log10_FDR_q"] = -np.log10(plotted["GSEA_FDR_q"].clip(lower=1e-300))
    plotted = plotted.sort_values("NES")
    source = plotted[
        [
            "pathway",
            "NES",
            "nominal_P",
            "GSEA_FDR_q",
            "panel_v2_leading_edge_genes",
            "n_panel_v2_genes_in_leading_edge",
            "minus_log10_FDR_q",
            "plot_role",
        ]
    ].copy()
    source_name = source_csv("fig5_source_data.csv", source)
    fig = plt.figure(figsize=(mm(180), mm(242)))
    grid = GridSpec(1, 2, width_ratios=[1.6, 1], wspace=0.45, figure=fig)
    ax = fig.add_subplot(grid[0, 0])
    positions = np.arange(len(plotted))
    norm = matplotlib.colors.Normalize(
        vmin=float(plotted["minus_log10_FDR_q"].min()),
        vmax=float(plotted["minus_log10_FDR_q"].max()),
    )
    bar_colors = plt.cm.viridis(norm(plotted["minus_log10_FDR_q"].to_numpy()))
    significant_edges = plotted["GSEA_FDR_q"].lt(0.05).to_numpy()
    ax.barh(
        positions,
        plotted["NES"],
        color=bar_colors,
        edgecolor=[COLORS["dark"] if is_sig else "none" for is_sig in significant_edges],
        linewidth=0.45,
        height=0.68,
    )
    labels = []
    for _, row in plotted.iterrows():
        panel_genes = str(row["panel_v2_genes_in_leading_edge"] or "")
        suffix = f"  [{panel_genes}]" if panel_genes and panel_genes.lower() != "nan" else ""
        labels.append(f"{row['pathway']}{suffix}")
    ax.set_yticks(positions, labels, fontsize=8)
    ax.axvline(0, color=COLORS["dark"], linewidth=0.6)
    ax.set_xlabel("Normalized enrichment score (NES)")
    set_panel_title(ax, "Hallmark pathways (FDR-significant + requested terms)")
    ax.invert_yaxis()
    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap="viridis"),
        ax=ax,
        fraction=0.035,
        pad=0.025,
        label="-log10 BH FDR q",
    )

    side = fig.add_subplot(grid[0, 1])
    side.axis("off")
    side_text = ["Requested pathway", "NES | nominal P | FDR q"]
    for label, token in targets.items():
        found = gsea.loc[pathways.str.contains(token, regex=False)]
        if found.empty:
            side_text.extend([label, "missing from saved table", ""])
            continue
        row = found.iloc[0]
        side_text.extend(
            [
                label,
                f"{row['NES']:.3f} | {row['nominal_P']:.4g} | {row['GSEA_FDR_q']:.4g}",
                "panel v2 in leading edge: "
                + (str(row["panel_v2_genes_in_leading_edge"]) if pd.notna(row["panel_v2_genes_in_leading_edge"]) else "none"),
                "",
            ]
        )
    side.text(0.02, 0.98, "\n".join(side_text), va="top", fontsize=8, linespacing=1.35)
    fig.suptitle(
        f"GSE125498 Hallmark prerank (n=13 LPD, 20 SPD); exploratory, same data as discovery; not validation",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.91, bottom=0.08, left=0.44, right=0.98)
    outputs = save_figure(fig, "Fig5_hallmark_gsea", files)
    go_note = (
        "GO BP/CC/MF, KEGG, and DO ECM-DEG plots were not generated: the available saved GO BP file "
        "is per WGCNA module, not an enrichment table for the 44 ECM DEGs; no matching ECM-only "
        "CC/MF/KEGG/DO tables were found."
    )
    add_legend(
        f"### Figure 5. Hallmark preranked GSEA\n"
        f"Saved Hallmark pathways with BH FDR q<0.05 and the five requested pathways are shown; significant "
        f"bars are emphasized. Bracketed genes in pathway labels are panel-v2 genes recorded in that "
        f"leading edge. Positive NES corresponds to higher t statistics ranked descending; in the saved "
        f"GSE125498 LPD-versus-SPD contrast, positive direction is higher in LPD. n=13 LPD and n=20 SPD. "
        f"Exploratory, same data as discovery; not validation. {go_note} "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )
    SKIPPED.append(("Fig5b GO/KEGG/DO", go_note))


def panel_evaluation_figure() -> None:
    files = [
        "results/revision/panel_v2/panel_detectability_audit.csv",
        "results/revision/effluent_validation/GSE248762/analysis_20261009/primary_per_gene_results.csv",
        "results/revision/effluent_validation/GSE130888/analysis_20261009/primary_per_gene_results.csv",
        "results/revision/gse62928_validation_8genes/per_gene_primary_results.csv",
        "results/revision/panel_v2/GSE92455_20261011/panel_v2_per_gene_limma.csv",
        "results/revision/panel_v2/panel_v2_report.md",
        "results/revision/panel_v2/GSE92455_20261011/panel_v2_score_summary.json",
        "results/revision/effluent_validation/GSE248762/analysis_20261009/patient_celltype_composition.csv",
        "results/revision/effluent_validation/GSE130888/analysis_20261009/patient_celltype_composition.csv",
        "results/revision/gse62928_validation_8genes/validation_summary.json",
    ]
    detect = read_csv(files[0])
    data_248 = read_csv(files[1])
    data_130 = read_csv(files[2])
    data_629 = read_csv(files[3])
    data_924 = read_csv(files[4])
    composition_248 = read_csv(files[7])
    composition_130 = read_csv(files[8])
    cohort_629 = read_json(files[9])
    score_924 = read_json(files[6])
    (ROOT / files[5]).read_text(encoding="utf-8")
    register_input(ROOT / files[5])

    effect_lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for cohort, frame, gene_col, effect, lower, upper in [
        ("GSE248762", data_248, "gene", "log2FC", "ci_lower", "ci_upper"),
        ("GSE130888", data_130, "gene", "log2FC", "ci_lower", "ci_upper"),
        (
            "GSE62928",
            data_629,
            "Gene",
            "Validation_log2FC",
            "Validation_CI95_Lower",
            "Validation_CI95_Upper",
        ),
        (
            "GSE92455",
            data_924,
            "gene",
            "logFC_GEO_value_scale",
            "CI95_lower",
            "CI95_upper",
        ),
    ]:
        for _, row in frame.iterrows():
            effect_lookup[(cohort, str(row[gene_col]))] = {
                "effect": row[effect],
                "lower": row[lower],
                "upper": row[upper],
            }
    detection_lookup = {
        (str(row.gene), str(row.cohort)): bool(row.detected)
        for row in detect.itertuples()
    }
    cohort_n = {
        "GSE248762": (
            f"LV {int((composition_248['primary_group'] == 'LV').sum())} / "
            f"SV {int((composition_248['primary_group'] == 'SV').sum())}"
        ),
        "GSE130888": (
            f"LV {int((composition_130['primary_group'] == 'LV').sum())} / "
            f"SV {int((composition_130['primary_group'] == 'SV').sum())}; "
            f"{int((composition_130['primary_group'] == 'NORMAL').sum())} NORMAL excluded"
        ),
        "GSE62928": f"n={int(cohort_629['qc']['sample_count'])}",
        "GSE92455": (
            f"epithelioid {int(score_924['n_epithelioid_patients'])} / "
            f"non-epithelioid {int(score_924['n_non_epithelioid_patients'])}"
        ),
    }
    verdicts = {
        "GSE248762": "Not supported",
        "GSE130888": "Partially supportive",
        "GSE62928": "Not supported",
        "GSE92455": "Not a vintage test",
    }
    rows = []
    for cohort in COHORTS:
        for gene in ORIGINAL_PANEL:
            values = effect_lookup.get((cohort, gene), {})
            detection = detection_lookup.get((gene, cohort), None) if cohort != "GSE92455" else None
            raw_sum = np.nan
            if gene == "SERPINA10" and cohort in {"GSE248762", "GSE130888"}:
                match = detect.loc[(detect["gene"] == gene) & (detect["cohort"] == cohort)]
                if not match.empty:
                    raw_sum = match.iloc[0]["raw_counts_sum"]
            rows.append(
                {
                    "record_type": "gene_cohort",
                    "cohort": cohort,
                    "cohort_sample_n": cohort_n[cohort],
                    "gene": gene,
                    "log2FC_or_saved_scale_effect": values.get("effect", np.nan),
                    "CI_lower": values.get("lower", np.nan),
                    "CI_upper": values.get("upper", np.nan),
                    "detected_by_saved_rule": detection,
                    "raw_counts_sum_SERPINA10": raw_sum,
                    "frozen_verdict": verdicts[cohort],
                }
            )
    source = pd.DataFrame(rows)
    source_name = source_csv("fig6_source_data.csv", source)
    fig = plt.figure(figsize=(mm(180), mm(248)))
    outer = GridSpec(3, 1, height_ratios=[2.25, 1.35, 0.38], hspace=0.42, figure=fig)
    forest_grid = GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[0], hspace=0.40, wspace=0.32)
    for index, cohort in enumerate(COHORTS):
        ax = fig.add_subplot(forest_grid[index // 2, index % 2])
        frame = source[source["cohort"] == cohort].set_index("gene").reindex(ORIGINAL_PANEL)
        ypos = np.arange(len(ORIGINAL_PANEL))
        good = frame["log2FC_or_saved_scale_effect"].notna() & frame["CI_lower"].notna() & frame["CI_upper"].notna()
        effects = frame.loc[good, "log2FC_or_saved_scale_effect"].astype(float).to_numpy()
        lows = frame.loc[good, "CI_lower"].astype(float).to_numpy()
        highs = frame.loc[good, "CI_upper"].astype(float).to_numpy()
        ygood = ypos[good.to_numpy()]
        ax.errorbar(
            effects,
            ygood,
            xerr=np.vstack([effects - lows, highs - effects]),
            fmt="o",
            color=COLORS["blue"],
            ecolor=COLORS["grey"],
            capsize=2,
            markersize=3.5,
            linewidth=0.8,
        )
        ax.axvline(0, color=COLORS["dark"], linewidth=0.6, linestyle="--")
        ax.set_yticks(ypos, ORIGINAL_PANEL if index % 2 == 0 else [])
        ax.invert_yaxis()
        ax.set_xlabel("Saved log2FC / platform-scale estimate")
        set_panel_title(ax, f"{cohort} ({cohort_n[cohort]})")
    detect_ax = fig.add_subplot(outer[1])
    detect_matrix = np.full((len(ORIGINAL_PANEL), len(COHORTS)), np.nan)
    for gi, gene in enumerate(ORIGINAL_PANEL):
        for ci, cohort in enumerate(COHORTS):
            detected = detection_lookup.get((gene, cohort), None) if cohort != "GSE92455" else None
            if detected is not None:
                detect_matrix[gi, ci] = int(detected)
    cmap = matplotlib.colors.ListedColormap(["#F3B6B6", "#8CC7E8"])
    detect_ax.imshow(np.ma.masked_invalid(detect_matrix), aspect="auto", cmap=cmap, vmin=0, vmax=1)
    detect_ax.set_xticks(np.arange(len(COHORTS)), COHORTS)
    detect_ax.set_yticks(np.arange(len(ORIGINAL_PANEL)), ORIGINAL_PANEL)
    for gi, gene in enumerate(ORIGINAL_PANEL):
        for ci, cohort in enumerate(COHORTS):
            detected = detection_lookup.get((gene, cohort), None) if cohort != "GSE92455" else None
            label = "not assessed" if detected is None else ("yes" if detected else "no")
            if gene == "SERPINA10" and cohort in {"GSE248762", "GSE130888"}:
                raw_sum = source.loc[
                    (source["record_type"] == "gene_cohort")
                    & (source["gene"] == gene)
                    & (source["cohort"] == cohort),
                    "raw_counts_sum_SERPINA10",
                ].iloc[0]
                label += f"\nraw={int(raw_sum)}"
            detect_ax.text(ci, gi, label, ha="center", va="center", fontsize=8)
    detect_ax.set_title("Detection by saved operational rule; SERPINA10 raw counts shown for effluent cohorts", loc="left", fontsize=8)
    verdict_ax = fig.add_subplot(outer[2])
    verdict_ax.axis("off")
    verdict_text = "Frozen verdicts (unchanged): " + "  |  ".join(
        f"{cohort}: {verdicts[cohort]}" for cohort in COHORTS
    )
    verdict_ax.text(0.01, 0.5, verdict_text, va="center", fontsize=8)
    fig.suptitle(
        "Panel-gene effects and detectability across cohorts — exploratory display; frozen verdicts unchanged",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.94, bottom=0.04, left=0.18, right=0.98)
    outputs = save_figure(fig, "Fig6_panel_cohorts", files)
    add_legend(
        f"### Figure 6. Panel-gene effects, detection, and frozen verdicts\n"
        f"Upper forest panels show saved per-gene estimates and 95% CIs in the same original-panel gene order "
        f"for GSE248762 (LV 10/SV 6), GSE130888 (LV 4/SV 6; 3 NORMAL samples excluded), GSE62928 (n=8), "
        f"and GSE92455 (epithelioid n=9/non-epithelioid n=8). Platform-specific effect scales are not "
        f"directly comparable. Detection heatmap follows the saved rule for the first three cohorts; "
        f"GSE92455 had no matching detectability rule and is marked not assessed. SERPINA10 raw count sums "
        f"are shown: GSE248762=2 and GSE130888=0. Verdicts are taken from the saved panel report: "
        f"GSE248762 Not supported; GSE130888 Partially supportive; GSE62928 Not supported; "
        f"GSE92455 not a vintage test. No frozen verdict was changed. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def roc_figure_skip() -> None:
    reason = (
        "Skipped: missing saved ROC curve-point coordinates for the GSE125498 panel-v2 composite. "
        "The saved file contains AUC summaries and patient scores, but reconstructing ROC points would "
        "recompute a performance analysis, which this task forbids."
    )
    SKIPPED.append(("Fig7 honest performance", reason))
    MANIFEST.append(
        {
            "figure": "Fig7 honest performance",
            "script": "results/revision/figures_final/make_figures.py",
            "script_sha256": "",
            "input_file": "missing: saved ROC curve-point table",
            "input_sha256": "",
            "output_files": "",
            "status": reason,
        }
    )


def immune_figure() -> None:
    files = [
        "results/revision/panel_v2_context/immune_celltype_group_comparison.csv",
        "results/revision/panel_v2_context/panel_celltype_correlations.csv",
        "results/revision/panel_v2_context/composite_leukocyte_adjustment.csv",
    ]
    groups = read_csv(files[0])
    correlations = read_csv(files[1])
    adjustment = read_csv(files[2])
    corr = correlations[correlations["score"].astype(str).str.startswith("panel_v2_")].copy()
    corr["gene"] = corr["score"].str.replace("panel_v2_", "", regex=False)
    corr["gene"] = corr["gene"].replace({"composite": "Composite"})
    gene_order = PANEL_V2 + ["Composite"]
    corr["gene"] = pd.Categorical(corr["gene"], categories=gene_order, ordered=True)
    corr = corr.sort_values(["gene", "cell_type"])
    selected_model = adjustment[
        adjustment["composite"].astype(str).eq("panel_v2_composite")
    ].copy()
    selected_model["model_short"] = selected_model["model"].map(
        {
            "composite ~ group": "Before leukocyte adjustment",
            "composite ~ group + dominant leukocyte marker proxy": "After leukocyte adjustment",
        }
    )
    source_frames = [
        groups[
            ["cell_type", "LPD_n", "SPD_n", "LPD_minus_SPD_mean", "BH_q_across_cell_types"]
        ].assign(record_type="celltype_group_difference"),
        corr[
            ["gene", "cell_type", "Pearson_r", "Pearson_BH_q_within_score"]
        ].assign(record_type="panel_celltype_correlation"),
        selected_model[
            [
                "model_short",
                "group_beta_LPD_vs_SPD",
                "group_P",
                "group_CI95_low",
                "group_CI95_high",
                "group_beta_attenuation_percent",
            ]
        ].assign(record_type="composite_adjustment"),
    ]
    source = pd.concat(source_frames, ignore_index=True, sort=False)
    source_name = source_csv("fig8_source_data.csv", source)
    fig = plt.figure(figsize=(mm(180), mm(220)))
    grid = GridSpec(3, 1, height_ratios=[1.05, 1.25, 0.75], hspace=0.48, figure=fig)
    ax1 = fig.add_subplot(grid[0])
    ordered_group = groups.sort_values("LPD_minus_SPD_mean")
    colors = [
        COLORS["orange"] if q < 0.05 else COLORS["grey"]
        for q in ordered_group["BH_q_across_cell_types"]
    ]
    ax1.barh(ordered_group["cell_type"], ordered_group["LPD_minus_SPD_mean"], color=colors)
    ax1.axvline(0, color=COLORS["dark"], linewidth=0.6)
    ax1.set_xlabel("Mean ssGSEA score difference (LPD - SPD)")
    ax1.set_title("a  Immune signature estimates by group (LPD n=13, SPD n=20)", loc="left", fontsize=8)
    ax1.tick_params(axis="y", labelsize=8)
    for y, (_, row) in enumerate(ordered_group.iterrows()):
        ax1.annotate(
            f"q={row['BH_q_across_cell_types']:.3g}",
            (1.02, y),
            xycoords=ax1.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8,
            clip_on=False,
        )

    ax2 = fig.add_subplot(grid[1])
    pivot = corr.pivot(index="gene", columns="cell_type", values="Pearson_r").reindex(gene_order)
    image = ax2.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax2.set_xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=45, ha="right", fontsize=8)
    ax2.set_yticks(np.arange(len(pivot.index)), [str(x) for x in pivot.index])
    ax2.set_title("b  Pearson correlation of panel genes/composite with cell-type scores", loc="left", fontsize=8)
    fig.colorbar(image, ax=ax2, fraction=0.02, pad=0.02, label="Pearson r")

    ax3 = fig.add_subplot(grid[2])
    selected_model = selected_model.set_index("model_short").reindex(
        ["Before leukocyte adjustment", "After leukocyte adjustment"]
    )
    y = np.arange(len(selected_model))
    betas = selected_model["group_beta_LPD_vs_SPD"].astype(float).to_numpy()
    low = selected_model["group_CI95_low"].astype(float).to_numpy()
    high = selected_model["group_CI95_high"].astype(float).to_numpy()
    ax3.errorbar(
        betas,
        y,
        xerr=np.vstack([betas - low, high - betas]),
        fmt="o",
        color=COLORS["blue"],
        ecolor=COLORS["grey"],
        capsize=3,
    )
    ax3.set_yticks(y, selected_model.index)
    ax3.axvline(0, color=COLORS["dark"], linewidth=0.6, linestyle="--")
    ax3.set_xlabel("LPD group coefficient (95% CI)")
    fig.subplots_adjust(right=0.86)
    ax3.set_title("c  Composite coefficient before/after leukocyte-marker-proxy adjustment", loc="left", fontsize=8)
    for i, (_, row) in enumerate(selected_model.iterrows()):
        ax3.annotate(
            f"beta={row['group_beta_LPD_vs_SPD']:.3f}; P={row['group_P']:.3g}",
            (row["group_beta_LPD_vs_SPD"], i),
            xytext=(7, -12 if i == 1 else 6),
            textcoords="offset points",
            fontsize=8,
        )
    fig.suptitle(
        "GSE125498 immune context — exploratory, post hoc, same data as discovery; not validation",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.93, bottom=0.08, left=0.20, right=0.95)
    outputs = save_figure(fig, "Fig8_immune_context", files)
    adjusted = selected_model.loc["After leukocyte adjustment"]
    before = selected_model.loc["Before leukocyte adjustment"]
    add_legend(
        f"### Figure 8. Immune context\n"
        f"(a) Saved ssGSEA signature-score group differences with BH q values used only to color the bars "
        f"(orange q<0.05). (b) Pearson correlations for five panel-v2 genes and the panel composite against "
        f"the saved cell-type signature scores. (c) Saved composite group coefficients before and after "
        f"adjustment for a dominant-leukocyte marker-score proxy: beta {before['group_beta_LPD_vs_SPD']:.3f} "
        f"before and {adjusted['group_beta_LPD_vs_SPD']:.3f} after; adjusted P={adjusted['group_P']:.5g}; "
        f"attenuation={adjusted['group_beta_attenuation_percent']:.1f}%. Samples: LPD n=13, SPD n=20. "
        f"These are signature estimates/proxy scores, not cell counts or direct cell-fraction deconvolution. "
        f"Exploratory, post hoc, same data as discovery; not validation. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def ppi_figure() -> None:
    files = [
        "results/revision/legacy_11_hub_ppi/network_nodes.csv",
        "results/revision/legacy_11_hub_ppi/string_physical_interactions.tsv",
        "results/revision/legacy_11_hub_functional_500/network_nodes.csv",
        "results/revision/legacy_11_hub_functional_500/string_interactions.tsv",
    ]
    physical_nodes = read_csv(files[0])
    physical_edges = read_csv(files[1], sep="\t")
    functional_nodes = read_csv(files[2])
    functional_edges = read_csv(files[3], sep="\t")
    node_set = sorted(
        set(physical_nodes["gene"].astype(str))
        | set(functional_nodes["gene"].astype(str))
    )
    edge_rows = []
    for kind, frame in [("Physical", physical_edges), ("Functional", functional_edges)]:
        for _, row in frame.iterrows():
            edge_rows.append(
                {
                    "network": kind,
                    "gene_a": row["preferredName_A"],
                    "gene_b": row["preferredName_B"],
                    "score": row["score"],
                    "score_cutoff": 0.5,
                }
            )
    edges = pd.DataFrame(edge_rows)
    nodes = pd.DataFrame({"gene": node_set})
    source = pd.concat(
        [
            edges.assign(record_type="edge"),
            nodes.assign(record_type="node"),
        ],
        ignore_index=True,
        sort=False,
    )
    source_name = source_csv("fig9_source_data.csv", source)
    fig, axes = plt.subplots(1, 2, figsize=(mm(180), mm(112)))
    angles = np.linspace(0, 2 * np.pi, len(node_set), endpoint=False)
    coordinates = {
        gene: (np.cos(angle), np.sin(angle)) for gene, angle in zip(node_set, angles)
    }
    for ax, kind, edge_frame in [
        (axes[0], "Physical", physical_edges),
        (axes[1], "Functional", functional_edges),
    ]:
        ax.set_aspect("equal")
        for _, edge in edge_frame.iterrows():
            a, b = str(edge["preferredName_A"]), str(edge["preferredName_B"])
            x1, y1 = coordinates[a]
            x2, y2 = coordinates[b]
            ax.plot(
                [x1, x2],
                [y1, y2],
                color=COLORS["blue"] if kind == "Physical" else COLORS["orange"],
                linewidth=0.8 + 2.0 * float(edge["score"]),
                zorder=1,
            )
            ax.text((x1 + x2) / 2, (y1 + y2) / 2, f"{float(edge['score']):.3f}", fontsize=8, ha="center")
        for gene, (x, y) in coordinates.items():
            ax.scatter(x, y, s=23, color=COLORS["green"], zorder=2)
            ax.text(x * 1.12, y * 1.12, gene, ha="center", va="center", fontsize=8)
        ax.set_xlim(-1.45, 1.45)
        ax.set_ylim(-1.35, 1.35)
        ax.axis("off")
        set_panel_title(ax, f"{kind} STRING (cutoff 0.5)")
    fig.suptitle(
        "Historical 11-gene PPI networks — panel v2 was not queried",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.88, bottom=0.03, left=0.02, right=0.98, wspace=0.10)
    outputs = save_figure(fig, "Fig9_historical_ppi_networks", files)
    isolated = len(node_set) - len(set(edges["gene_a"].astype(str)) | set(edges["gene_b"].astype(str)))
    add_legend(
        f"### Figure 9. Saved STRING physical and functional networks\n"
        f"Only saved STRING edge tables are used; no live query was made. Both panels show the historical "
        f"11-gene list (not panel v2), with the saved score cutoff 0.5 and isolated nodes retained. "
        f"Physical: {len(physical_edges)} edge(s), including LTA-LTB score 0.998. Functional: "
        f"{len(functional_edges)} edges, including LTA-LTB 0.999 and LTA-FLT3LG 0.575. "
        f"{isolated} nodes are isolated in the union of displayed edges. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def single_cell_figure() -> None:
    files = [
        "results/revision/effluent_validation/GSE248762/analysis_20261009/patient_celltype_composition.csv",
        "results/revision/effluent_validation/GSE248762/analysis_20261009/gene_expression_by_cell_type.csv",
        "results/revision/effluent_validation/GSE130888/analysis_20261009/patient_celltype_composition.csv",
        "results/revision/effluent_validation/GSE130888/analysis_20261009/gene_expression_by_cell_type.csv",
        "results/revision/exploratory_posthoc_20261009/ecm1_posthoc_cell_concentration.csv",
    ]
    comp_248 = read_csv(files[0])
    expr_248 = read_csv(files[1])
    comp_130 = read_csv(files[2])
    expr_130 = read_csv(files[3])
    ecm1 = read_csv(files[4])
    composition_rows = []
    for cohort, frame in [("GSE248762", comp_248), ("GSE130888", comp_130)]:
        long = frame.melt(
            id_vars=["sample_id", "primary_group", "n_qc_cells"],
            value_vars=[column for column in frame if column.startswith("fraction_")],
            var_name="cell_type_column",
            value_name="fraction",
        )
        long["cohort"] = cohort
        long["cell_type"] = long["cell_type_column"].str.replace("fraction_", "", regex=False)
        composition_rows.append(long)
    composition = pd.concat(composition_rows, ignore_index=True)
    composition = composition[
        ["cohort", "sample_id", "primary_group", "n_qc_cells", "cell_type", "fraction"]
    ].copy()
    composition["record_type"] = "sample_celltype_composition"
    expression_rows = []
    for cohort, frame in [("GSE248762", expr_248), ("GSE130888", expr_130)]:
        selected = frame[frame["gene"].isin(PANEL_V2)].copy()
        selected["cohort"] = cohort
        selected["marker_size"] = 20 + 160 * selected["fraction_cells_nonzero"].astype(float)
        expression_rows.append(selected)
    expression = pd.concat(expression_rows, ignore_index=True)
    expression = expression[
        [
            "cohort",
            "gene",
            "cell_type",
            "fraction_cells_nonzero",
            "mean_log1p_normalized_expression",
            "marker_size",
        ]
    ].copy()
    expression["record_type"] = "panel_gene_celltype_expression"
    ecm1 = ecm1.copy()
    comp_frac = comp_248.set_index("sample_id")["fraction_fibroblast"]
    ecm1["fibroblast_cells_from_saved_fraction"] = (
        ecm1["sample_id"].map(comp_frac) * ecm1["qc_cells"]
    ).round().astype(int)
    ecm1["log1p_raw_counts_plotted"] = np.log1p(ecm1["ecm1_total_raw_counts"])
    ecm1 = ecm1[
        [
            "sample_id",
            "primary_group",
            "qc_cells",
            "ecm1_total_raw_counts",
            "ecm1_positive_cells",
            "fibroblast_cells_from_saved_fraction",
            "log1p_raw_counts_plotted",
        ]
    ].copy()
    ecm1["record_type"] = "ECM1_sample_raw_counts"
    source = pd.concat([composition, expression, ecm1], ignore_index=True, sort=False)
    source_name = source_csv("fig10_source_data.csv", source)
    fig = plt.figure(figsize=(mm(180), mm(242)))
    grid = GridSpec(3, 1, height_ratios=[1.15, 1.15, 0.9], hspace=0.56, figure=fig)
    comp_grid = GridSpecFromSubplotSpec(1, 2, subplot_spec=grid[0], wspace=0.22)
    for col, (cohort, frame) in enumerate([("GSE248762", comp_248), ("GSE130888", comp_130)]):
        ax = fig.add_subplot(comp_grid[0, col])
        stack_cols = [c for c in frame.columns if c.startswith("fraction_")]
        x = np.arange(len(frame))
        bottoms = np.zeros(len(frame))
        for i, column in enumerate(stack_cols):
            values = frame[column].fillna(0).astype(float).to_numpy()
            ax.bar(x, values, bottom=bottoms, width=0.8, label=column.replace("fraction_", ""))
            bottoms += values
        ax.set_xticks(x, frame["sample_id"], rotation=70, ha="right", fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Cell-type fraction")
        set_panel_title(ax, f"{cohort} (n={len(frame)} samples)")
        if col == 1:
            ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    expr_grid = GridSpecFromSubplotSpec(1, 2, subplot_spec=grid[1], wspace=0.32)
    cell_types = list(dict.fromkeys(expression["cell_type"].astype(str)))
    for col, cohort in enumerate(["GSE248762", "GSE130888"]):
        ax = fig.add_subplot(expr_grid[0, col])
        current = expression[expression["cohort"] == cohort]
        for gene in PANEL_V2:
            sub = current[current["gene"] == gene]
            ax.scatter(
                [gene] * len(sub),
                sub["cell_type"],
                s=sub["marker_size"],
                c=sub["mean_log1p_normalized_expression"],
                cmap="viridis",
                vmin=0,
                vmax=max(1.0, float(expression["mean_log1p_normalized_expression"].max())),
                edgecolors="black",
                linewidths=0.25,
            )
        ax.set_xticks(range(len(PANEL_V2)), PANEL_V2, rotation=35, ha="right")
        ax.set_yticks(range(len(cell_types)), cell_types)
        if col == 1:
            ax.set_yticklabels([])
        set_panel_title(ax, f"{cohort} panel-v2 expression")
    fig.text(0.92, 0.43, "Dot size = fraction cells nonzero\nColor = mean log1p normalized expression", fontsize=8, va="center")

    ax = fig.add_subplot(grid[2])
    ecm1 = ecm1.sort_values("ecm1_total_raw_counts")
    bars = ax.bar(
        ecm1["sample_id"],
        ecm1["log1p_raw_counts_plotted"],
        color=[COLORS["orange"] if sid == "GSM7919584" else COLORS["blue"] for sid in ecm1["sample_id"]],
    )
    ax.set_xticks(np.arange(len(ecm1)), ecm1["sample_id"], rotation=65, ha="right", fontsize=8)
    ax.set_ylabel("log(1 + raw ECM1 counts)")
    set_panel_title(ax, "ECM1 raw count concentration by GSE248762 sample")
    concern = ecm1.loc[ecm1["sample_id"] == "GSM7919584"].iloc[0]
    idx = list(ecm1["sample_id"]).index("GSM7919584")
    ax.annotate(
        f"{int(concern['ecm1_total_raw_counts']):,} counts\n"
        f"{int(concern['ecm1_positive_cells']):,} positive cells; "
        f"{int(concern['fibroblast_cells_from_saved_fraction'])} fibroblasts\n"
        f"{int(concern['qc_cells']):,} QC cells",
        (idx, bars[idx].get_height()),
        xytext=(5, 5),
        textcoords="offset points",
        fontsize=8,
    )
    fig.suptitle(
        "Single-cell context — supporting context, not validation; small n; no re-clustering",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.93, bottom=0.15, left=0.17, right=0.80)
    outputs = save_figure(fig, "Fig10_single_cell_context", files)
    add_legend(
        f"### Figure 10. Single-cell supporting context\n"
        f"(a) Saved per-sample cell-type composition for GSE248762 (n={len(comp_248)} samples) and "
        f"GSE130888 (n={len(comp_130)} samples). (b) Saved normalized expression of five panel-v2 genes "
        f"across annotated cell types; dot size is the saved fraction of cells nonzero and color is saved "
        f"mean log1p normalized expression. (c) Saved per-sample raw ECM1 counts in GSE248762. "
        f"GSM7919584 has {int(concern['ecm1_total_raw_counts']):,} counts in {int(concern['qc_cells']):,} QC cells; "
        f"the saved fibroblast fraction corresponds to {int(concern['fibroblast_cells_from_saved_fraction'])} "
        f"cells. This is supporting context, not validation; pseudobulk limits, small n, marker-based "
        f"annotations, and patient overlap unknown. No re-clustering or statistical testing was performed. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def detectability_supplement() -> None:
    files = [
        "results/revision/panel_v2/panel_v2_rule_results.csv",
        "results/revision/panel_v2/panel_detectability_audit.csv",
    ]
    rules = read_csv(files[0])
    detect = read_csv(files[1])
    rule_source = rules[
        [
            "gene",
            "rule_a_detected_in_at_least_2_of_3",
            "rule_b_concordant_in_every_detected_cohort",
            "discovery_BH_q_lt_0_05_flag",
            "passes_panel_v2_selection",
        ]
    ].copy()
    rule_source["record_type"] = "panel_rule_result"
    detection_rows = detect[
        ["gene", "cohort", "detected", "raw_counts_sum"]
    ].copy()
    detection_rows["record_type"] = "cohort_detectability"
    source = pd.concat([rule_source, detection_rows], ignore_index=True, sort=False)
    source_name = source_csv("figS3_source_data.csv", source)
    fig, axes = plt.subplots(1, 2, figsize=(mm(180), mm(118)), gridspec_kw={"width_ratios": [1.2, 1]})
    rule_cols = [
        "rule_a_detected_in_at_least_2_of_3",
        "rule_b_concordant_in_every_detected_cohort",
        "discovery_BH_q_lt_0_05_flag",
        "passes_panel_v2_selection",
    ]
    rule_labels = ["Detected >=2/3", "Concordant where detected", "BH q flag only", "Selected by a+b"]
    rule_matrix = np.zeros((len(rules), len(rule_cols)))
    for col_i, column in enumerate(rule_cols):
        rule_matrix[:, col_i] = rules[column].map({"PASS": 1, "FAIL": 0, True: 1, False: 0}).fillna(-1)
    ax = axes[0]
    image = ax.imshow(rule_matrix, aspect="auto", cmap=matplotlib.colors.ListedColormap(["#D9D9D9", "#F3B6B6", "#8CC7E8"]), vmin=-1, vmax=1)
    ax.set_yticks(np.arange(len(rules)), rules["gene"])
    ax.set_xticks(np.arange(len(rule_labels)), rule_labels, rotation=35, ha="right")
    for i in range(len(rules)):
        for j, col in enumerate(rule_cols):
            raw = rules.iloc[i][col]
            text = "pass" if str(raw).upper() == "PASS" or raw is True else "fail"
            ax.text(j, i, text, ha="center", va="center", fontsize=8)
    set_panel_title(ax, "Original 8-gene rule outcomes")

    detection_matrix = np.zeros((len(rules), 3), dtype=float)
    detection_columns = ["GSE62928", "GSE248762", "GSE130888"]
    for col_i, cohort in enumerate(detection_columns):
        cohort_calls = detect.loc[detect["cohort"] == cohort].set_index("gene")["detected"]
        detection_matrix[:, col_i] = rules["gene"].map(cohort_calls).map(
            {True: 1, False: 0}
        ).fillna(-1).astype(float)
    axes[1].imshow(
        detection_matrix,
        aspect="auto",
        cmap=matplotlib.colors.ListedColormap(["#D9D9D9", "#F3B6B6", "#8CC7E8"]),
        vmin=-1,
        vmax=1,
    )
    axes[1].set_yticks(np.arange(len(rules)), rules["gene"])
    axes[1].set_xticks(np.arange(3), ["GSE62928", "GSE248762", "GSE130888"], rotation=35, ha="right")
    for i in range(len(rules)):
        for j in range(3):
            value = detection_matrix[i, j]
            axes[1].text(j, i, "Y" if value == 1 else "N", ha="center", va="center", fontsize=8)
    set_panel_title(axes[1], "Saved detectability calls")
    fig.suptitle(
        "Supplementary Figure 3 — original panel detectability rules; n=8 genes across three cohorts",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0.01, 0.02, 0.99, 0.91])
    outputs = save_figure(fig, "FigS3_detectability_rules", files)
    add_legend(
        f"### Supplementary Figure 3. Original eight-gene detectability rules\n"
        f"Rule outcomes are taken directly from the saved panel-v2 rule-results table for the original "
        f"eight genes. The second panel shows saved detection calls in GSE62928 (n=8), GSE248762 (n=16), "
        f"and GSE130888 (n=10 primary-contrast samples; 13 total including 3 NORMAL samples). "
        f"BH q<0.05 remains a recorded flag, not a hard selection criterion. "
        f"Outputs: {', '.join(outputs)}; source data: {source_name}."
    )


def write_documentation() -> None:
    skipped_notes = [
        (
            "Supplementary Figure 1 — skipped: no saved panel-v2 miRNA/lncRNA regulatory result table "
            "was found. No database query was made."
        ),
        (
            "Supplementary Figure 2 — skipped: the saved nomogram tables are legacy and do not document "
            "a model for the current panel v2. No nomogram was built."
        ),
    ]
    for item in skipped_notes:
        SKIPPED.append((item.split(" — ")[0], item))
    skip_section = "\n".join(f"- **{name}:** {reason}" for name, reason in SKIPPED)
    legends_text = (
        "# Publication figure legends\n\n"
        "All figures were generated from saved result files only. No differential expression, WGCNA, "
        "machine learning, GSEA, immune scoring, ROC calculation, or statistical test was run. "
        "All exploratory results are labeled as such. PNGs are 300 dpi and PDFs are vector exports. "
        "Reviewed final renders use the `_reviewed15` filename suffix. Prior draft renders are retained "
        "and marked superseded in the manifest. The `_reviewed15` set is the only final figure set.\n\n"
        + "\n\n".join(LEGENDS)
        + "\n\n## Skipped panels and supplements\n\n"
        + skip_section
        + "\n"
    )
    legends_path = OUT / "figure_legends.md"
    manifest_path = OUT / "figure_manifest.csv"
    legends_path.write_text(legends_text, encoding="utf-8")

    script_hash = sha256(Path(__file__))
    draft_stems = [
        "Fig1_study_workflow",
        "Fig2_discovery_volcano",
        "Fig3_overlap",
        "Fig4_three_model_union",
        "Fig5_hallmark_gsea",
        "Fig6_panel_cohorts",
        "Fig8_immune_context",
        "Fig9_historical_ppi_networks",
        "Fig10_single_cell_context",
        "FigS3_detectability_rules",
    ]
    for stem in draft_stems:
        MANIFEST.append(
            {
                "figure": f"{stem} (superseded first-pass render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}.png;{stem}.pdf",
                "status": "Superseded after visual QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_final (superseded visual-QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_final.png;{stem}_final.pdf",
                "status": "Superseded after visual QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed (superseded visual-QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed.png;{stem}_reviewed.pdf",
                "status": "Superseded after visual QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed2 (superseded visual-QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed2.png;{stem}_reviewed2.pdf",
                "status": "Superseded after visual QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed3 (superseded partial render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed3.png;{stem}_reviewed3.pdf",
                "status": "Superseded after partial output during correction; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed4 (superseded visual-QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed4.png;{stem}_reviewed4.pdf",
                "status": "Superseded after UpSet layout QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed5 (superseded visual-QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed5.png;{stem}_reviewed5.pdf",
                "status": "Superseded after workflow layout QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed6 (superseded visual-QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed6.png;{stem}_reviewed6.pdf",
                "status": "Superseded after UpSet axis QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed7 (superseded visual-QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed7.png;{stem}_reviewed7.pdf",
                "status": "Superseded after exploratory-label QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed8 (superseded font-size QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed8.png;{stem}_reviewed8.pdf",
                "status": "Superseded after font-size QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed9 (superseded layout QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed9.png;{stem}_reviewed9.pdf",
                "status": "Superseded after final layout QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed10 (superseded label QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed10.png;{stem}_reviewed10.pdf",
                "status": "Superseded after label QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed11 (superseded data-label QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed11.png;{stem}_reviewed11.pdf",
                "status": "Superseded after workflow and statistics-panel QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed12 (superseded vote-table QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed12.png;{stem}_reviewed12.pdf",
                "status": "Superseded after vote-table layout QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed13 (superseded stats-layout QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed13.png;{stem}_reviewed13.pdf",
                "status": "Superseded after statistics and vote-table spacing QA; retained and not overwritten.",
            }
        )
        MANIFEST.append(
            {
                "figure": f"{stem}_reviewed14 (superseded vote-table width QA render)",
                "script": "results/revision/figures_final/make_figures.py",
                "script_sha256": script_hash,
                "input_file": "",
                "input_sha256": "",
                "output_files": f"{stem}_reviewed14.png;{stem}_reviewed14.pdf",
                "status": "Superseded after vote-table column width QA; retained and not overwritten.",
            }
        )
    for row in MANIFEST:
        row["script_sha256"] = script_hash
    manifest_rows = list(MANIFEST)
    pd.DataFrame(
        manifest_rows,
        columns=[
            "figure",
            "script",
            "script_sha256",
            "input_file",
            "input_sha256",
            "output_files",
            "status",
        ],
    ).to_csv(manifest_path, index=False)


def preflight() -> None:
    if OUT.exists():
        draft_stems = [
            "Fig1_study_workflow",
            "Fig2_discovery_volcano",
            "Fig3_overlap",
            "Fig4_three_model_union",
            "Fig5_hallmark_gsea",
            "Fig6_panel_cohorts",
            "Fig8_immune_context",
            "Fig9_historical_ppi_networks",
            "Fig10_single_cell_context",
            "FigS3_detectability_rules",
        ]
        allowed = {
            "make_figures.py",
            "figure_manifest.csv",
            "figure_legends.md",
            "fig1_source_data.csv",
            "fig2_source_data.csv",
            "fig3_source_data.csv",
            "fig4_source_data.csv",
            "fig5_source_data.csv",
            "fig6_source_data.csv",
            "fig8_source_data.csv",
            "fig9_source_data.csv",
            "fig10_source_data.csv",
            "figS3_source_data.csv",
        }
        allowed.update(f"{stem}.{extension}" for stem in draft_stems for extension in ("png", "pdf"))
        allowed.update(f"{stem}_final.{extension}" for stem in draft_stems for extension in ("png", "pdf"))
        existing = [
            p.name
            for p in OUT.iterdir()
            if p.name not in allowed
            and not any(
                re.fullmatch(
                    re.escape(stem) + r"(?:_final|_reviewed\d*)?\.(?:png|pdf)",
                    p.name,
                )
                for stem in draft_stems
            )
        ]
        if existing:
            raise FileExistsError(
                f"Output folder is not empty; refusing to overwrite anything: {existing}"
            )
    else:
        OUT.mkdir(parents=True)


def main() -> None:
    np.random.seed(SEED)
    preflight()
    flow_figure()
    volcano_figure()
    overlap_figure()
    ml_figure()
    gsea_figure()
    panel_evaluation_figure()
    roc_figure_skip()
    immune_figure()
    ppi_figure()
    single_cell_figure()
    detectability_supplement()
    write_documentation()


if __name__ == "__main__":
    main()
