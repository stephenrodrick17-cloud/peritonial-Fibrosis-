"""Plot ROC curves for every detectable gene in the consensus hub panel."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results" / "revision" / "ml_39_matrisome_wgcna" / "validation_auc"
PANEL_PATH = ROOT / "results" / "revision" / "ml_39_matrisome_wgcna" / "ml_39_gene_results.csv"
BOOTSTRAP_REPLICATES = 2_000
SEED = 42
OUTPUT_STEM = "all_detectable_consensus_hub_gene_ROC"

GENES = [
    "TNFSF15",
    "EBI3",
    "FLT3LG",
    "LTB",
    "ADAM19",
    "ECM1",
    "CST7",
    "NAT6",
    "CCBE1",
]
GENE_IDS = {
    "TNFSF15": "ENSG00000181634",
    "EBI3": "ENSG00000105246",
    "FLT3LG": "ENSG00000090554",
    "LTB": "ENSG00000227507",
    "ADAM19": "ENSG00000135074",
    "ECM1": "ENSG00000143369",
    "CST7": "ENSG00000077984",
    "NAT6": "ENSG00000243477",
    "CCBE1": "ENSG00000183287",
}
GENE_COLORS = {
    "TNFSF15": "#0072B2",
    "EBI3": "#E69F00",
    "FLT3LG": "#009E73",
    "LTB": "#D55E00",
    "ADAM19": "#CC79A7",
    "ECM1": "#56B4E9",
    "CST7": "#000000",
    "NAT6": "#F0E442",
    "CCBE1": "#882255",
}
GENE_MARKERS = {
    "TNFSF15": "o",
    "EBI3": "s",
    "FLT3LG": "^",
    "LTB": "D",
    "ADAM19": "v",
    "ECM1": "P",
    "CST7": "X",
    "NAT6": "*",
    "CCBE1": "h",
}
COHORTS = {
    "GSE248762": {
        "case": "LV",
        "control": "SV",
        "title": "LV vs SV",
        "contrast": "Long-vintage (LV) vs short-vintage (SV)",
        "n_pos_expected": 10,
        "n_neg_expected": 6,
    },
    "GSE130888": {
        "case": "LV",
        "control": "SV",
        "title": "Long-term PD vs short-term PD",
        "contrast": "Long-term PD vs short-term PD",
        "n_pos_expected": 4,
        "n_neg_expected": 6,
    },
}


def _is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _bootstrap_auc_ci(
    labels: np.ndarray,
    scores: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    positive_scores = scores[labels == 1]
    negative_scores = scores[labels == 0]
    positive_samples = rng.choice(
        positive_scores,
        size=(BOOTSTRAP_REPLICATES, len(positive_scores)),
        replace=True,
    )
    negative_samples = rng.choice(
        negative_scores,
        size=(BOOTSTRAP_REPLICATES, len(negative_scores)),
        replace=True,
    )
    pairwise_auc = (
        (positive_samples[:, :, None] > negative_samples[:, None, :]).astype(float)
        + 0.5 * (positive_samples[:, :, None] == negative_samples[:, None, :])
    ).mean(axis=(1, 2))
    low, high = np.percentile(pairwise_auc, [2.5, 97.5])
    return float(low), float(high)


def _delong_auc_ci(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float]:
    """Calculate the single-marker DeLong Wald interval from U-statistic placements."""
    positive_scores = scores[labels == 1]
    negative_scores = scores[labels == 0]
    comparison = (
        (positive_scores[:, None] > negative_scores[None, :]).astype(float)
        + 0.5 * (positive_scores[:, None] == negative_scores[None, :])
    )
    auc_value = float(comparison.mean())
    positive_placements = comparison.mean(axis=1)
    negative_placements = comparison.mean(axis=0)
    variance = (
        np.var(positive_placements, ddof=1) / len(positive_scores)
        + np.var(negative_placements, ddof=1) / len(negative_scores)
    )
    standard_error = float(np.sqrt(variance))
    return (
        max(0.0, auc_value - 1.96 * standard_error),
        min(1.0, auc_value + 1.96 * standard_error),
    )


def _load_cohort(
    cohort: str,
    settings: dict[str, object],
    panel_directions: dict[str, int],
    log_lines: list[str],
) -> dict[str, object]:
    cohort_dir = (
        ROOT
        / "results"
        / "revision"
        / "effluent_validation"
        / cohort
        / "analysis_20261009"
    )
    metadata = pd.read_csv(cohort_dir / "all_cell_pseudobulk_metadata.csv")
    metadata = metadata.loc[
        metadata["primary_group"].isin([settings["case"], settings["control"]])
    ].drop_duplicates("sample_id")
    metadata = metadata.set_index("sample_id")
    labels = (metadata["primary_group"] == settings["case"]).astype(int).to_numpy()
    n_pos = int(labels.sum())
    n_neg = int(len(labels) - n_pos)
    if (n_pos, n_neg) != (
        settings["n_pos_expected"],
        settings["n_neg_expected"],
    ):
        raise ValueError(f"{cohort} group sizes were {n_pos}/{n_neg}, contrary to expected values.")

    expression_path = cohort_dir / "edger_logCPM.csv"
    expression = pd.read_csv(expression_path, index_col=0)
    expression.index = expression.index.astype(str)
    missing_samples = set(metadata.index) - set(expression.columns)
    if missing_samples:
        raise ValueError(f"{cohort} logCPM matrix is missing samples: {sorted(missing_samples)}")

    detection = pd.read_csv(
        cohort_dir / "pseudobulk_gene_detection.csv",
        dtype={"gene_id": str},
    ).set_index("gene_id")
    edger = pd.read_csv(cohort_dir / "edger_all_gene_results.csv")
    edger = edger.loc[edger["contrast"] == "LV_vs_SV"].set_index("gene_id")

    log_lines.append(
        f"{cohort}: {len(metadata)} samples; n_pos={n_pos}; n_neg={n_neg}; "
        f"score matrix={expression_path.relative_to(ROOT)}"
    )
    log_lines.append(
        f"  TMM log2-CPM values used; observed panel range "
        f"{expression.loc[expression.index.intersection(GENE_IDS.values()), metadata.index].min().min():.4f} "
        f"to {expression.loc[expression.index.intersection(GENE_IDS.values()), metadata.index].max().max():.4f}"
    )

    gene_data: dict[str, dict[str, object]] = {}
    for gene in GENES:
        gene_id = GENE_IDS[gene]
        detected = gene_id in detection.index and _is_true(detection.loc[gene_id, "detected"])
        if not detected:
            gene_data[gene] = {
                "detected": False,
                "scores": None,
                "concordant": None,
                "n_pos": n_pos,
                "n_neg": n_neg,
            }
            log_lines.append(f"  {gene}: not detected; no ROC/AUC calculated")
            continue
        if gene_id not in expression.index:
            raise ValueError(f"{cohort} logCPM matrix is missing detected gene {gene} ({gene_id}).")
        if gene_id not in edger.index:
            raise ValueError(f"{cohort} edgeR table is missing primary-contrast gene {gene} ({gene_id}).")
        scores = expression.loc[gene_id, metadata.index].to_numpy(dtype=float)
        edge_log_fc = float(edger.loc[gene_id, "logFC"])
        concordant = bool(
            edge_log_fc != 0
            and np.sign(edge_log_fc) == panel_directions[gene]
        )
        gene_data[gene] = {
            "detected": True,
            "scores": scores,
            "concordant": concordant,
            "edge_log_fc": edge_log_fc,
            "n_pos": n_pos,
            "n_neg": n_neg,
        }
    return {
        "cohort": cohort,
        "settings": settings,
        "metadata": metadata,
        "labels": labels,
        "gene_data": gene_data,
    }


def _calculate_gene_rows(cohorts: dict[str, dict[str, object]]) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for cohort, data in cohorts.items():
        labels = data["labels"]
        for gene in GENES:
            gene_data = data["gene_data"][gene]
            row = {
                "gene": gene,
                "cohort": cohort,
                "auc": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "concordant": gene_data["concordant"],
                "n_pos": gene_data["n_pos"],
                "n_neg": gene_data["n_neg"],
                "detected": gene_data["detected"],
                "tied_positive_negative_pairs": np.nan,
                "n_tied_score_values": np.nan,
                "n_tied_score_observations": np.nan,
                "auc_delong_ci_low": np.nan,
                "auc_delong_ci_high": np.nan,
            }
            if gene_data["detected"]:
                scores = gene_data["scores"]
                positives = scores[labels == 1]
                negatives = scores[labels == 0]
                auc_value = float(roc_auc_score(labels, scores))
                ci_low, ci_high = _bootstrap_auc_ci(labels, scores, rng)
                delong_low, delong_high = _delong_auc_ci(labels, scores)
                score_counts = pd.Series(scores).value_counts()
                repeated_scores = score_counts.loc[score_counts > 1]
                row.update(
                    {
                        "auc": auc_value,
                        "ci_low": ci_low,
                        "ci_high": ci_high,
                        "auc_delong_ci_low": delong_low,
                        "auc_delong_ci_high": delong_high,
                        "n_tied_score_values": len(repeated_scores),
                        "n_tied_score_observations": int(repeated_scores.sum()),
                        "tied_positive_negative_pairs": int(
                            np.equal(positives[:, None], negatives[None, :]).sum()
                        ),
                    }
                )
            rows.append(row)
    return pd.DataFrame(rows)


def _format_auc(row: pd.Series) -> str:
    return f"AUC {row['auc']:.2f} [{row['ci_low']:.2f}-{row['ci_high']:.2f}]"


def _style_for(gene: str, concordant: bool) -> dict[str, object]:
    return {
        "color": GENE_COLORS[gene],
        "linestyle": "-" if concordant else "--",
        "linewidth": 1.8 if concordant else 1.25,
        "marker": GENE_MARKERS[gene],
        "markersize": 3.5,
        "markerfacecolor": "white",
        "markeredgewidth": 0.8,
        "clip_on": False,
    }


def _curve_handles(
    cohort: str,
    cohort_data: dict[str, object],
    result_by_gene: pd.DataFrame,
) -> list[Line2D]:
    handles = []
    ordered_genes = sorted(
        (
            gene
            for gene in GENES
            if cohort_data["gene_data"][gene]["detected"]
            and cohort_data["gene_data"][gene]["concordant"]
        )
    ) + sorted(
        (
            gene
            for gene in GENES
            if cohort_data["gene_data"][gene]["detected"]
            and not cohort_data["gene_data"][gene]["concordant"]
        )
    )
    for gene in ordered_genes:
        gene_data = cohort_data["gene_data"][gene]
        result = result_by_gene.loc[gene]
        concordant = bool(gene_data["concordant"])
        suffix = "" if concordant else " (discordant)"
        style = _style_for(gene, concordant)
        handles.append(
            Line2D(
                [0],
                [0],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=style["linewidth"],
                marker=style["marker"],
                markersize=style["markersize"],
                markerfacecolor="white",
                label=f"{gene}: {_format_auc(result)}{suffix}",
            )
        )
    return handles


def _draw_curve(
    ax: plt.Axes,
    labels: np.ndarray,
    gene: str,
    gene_data: dict[str, object],
) -> None:
    scores = gene_data["scores"]
    fpr, tpr, _ = roc_curve(labels, scores, pos_label=1)
    ax.plot(fpr, tpr, drawstyle="steps-post", **_style_for(gene, gene_data["concordant"]))


def _format_panel(
    ax: plt.Axes,
    cohort: str,
    settings: dict[str, object],
    cohort_data: dict[str, object],
    ylabel: bool,
) -> None:
    labels = cohort_data["labels"]
    gene_data = cohort_data["gene_data"]
    for gene in GENES:
        if gene_data[gene]["detected"]:
            _draw_curve(ax, labels, gene, gene_data[gene])
    ax.plot([0, 1], [0, 1], linestyle=":", color="#475569", linewidth=1.2, zorder=0, clip_on=False)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate" if ylabel else "")
    ax.set_title(
        f"{cohort}: {settings['title']} "
        f"(n={settings['n_pos_expected']} positive, n={settings['n_neg_expected']} negative)",
        fontsize=11,
        fontweight="bold",
    )
    ax.grid(color="#E2E8F0", linewidth=0.7)
    ax.set_axisbelow(True)
    missing = [gene for gene in GENES if not gene_data[gene]["detected"]]
    note = f"Not detected: {', '.join(missing)}" if missing else "All nine panel genes detected"
    ax.text(
        0.5,
        -0.12,
        note,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8.3,
        color="#475569",
    )


def _make_overview(
    cohorts: dict[str, dict[str, object]],
    results: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(17, 10))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[4.0, 1.55],
        left=0.07,
        right=0.98,
        top=0.87,
        bottom=0.17,
        wspace=0.08,
        hspace=0.33,
    )
    fig.suptitle(
        "Individual ROC curves for the nine consensus hub genes",
        fontsize=16,
        fontweight="bold",
        y=0.975,
    )
    fig.text(
        0.5,
        0.935,
        "Higher expression predicts LV / long-term PD. Solid = direction-concordant; dashed = discordant.",
        ha="center",
        fontsize=10,
        color="#475569",
    )
    for index, (cohort, cohort_data) in enumerate(cohorts.items()):
        ax = fig.add_subplot(grid[0, index])
        legend_ax = fig.add_subplot(grid[1, index])
        legend_ax.axis("off")
        _format_panel(ax, cohort, cohort_data["settings"], cohort_data, ylabel=index == 0)
        cohort_results = results.loc[results["cohort"] == cohort].set_index("gene")
        legend_ax.legend(
            handles=_curve_handles(cohort, cohort_data, cohort_results),
            loc="center",
            frameon=False,
            fontsize=9,
            ncol=2,
            title="Solid = direction-concordant; dashed = discordant. AUC [95% bootstrap CI]",
            title_fontsize=9,
            handlelength=2.1,
            labelspacing=0.8,
            columnspacing=1.6,
        )
    fig.text(
        0.5,
        0.055,
        "TMM-normalized log2-CPM; positive class = LV / long-term PD; direction not flipped; "
        "95% CIs use 2,000 stratified bootstrap resamples.\n"
        "Concordance is defined by the sign of the edgeR log2FC, so AUC and concordance can disagree "
        "(for example LTB in GSE130888).\n"
        "No correction for testing nine genes; estimates are exploratory.",
        ha="center",
        va="center",
        fontsize=8.5,
        color="#475569",
    )
    fig.savefig(RESULTS_DIR / f"{OUTPUT_STEM}.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS_DIR / f"{OUTPUT_STEM}.svg", bbox_inches="tight")
    plt.close(fig)


def _make_small_multiples(
    cohort: str,
    cohort_data: dict[str, object],
    results: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(
        3,
        3,
        figsize=(12, 12),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(top=0.91, bottom=0.20, left=0.08, right=0.98, hspace=0.34, wspace=0.25)
    fig.suptitle(
        f"{cohort}: individual consensus-gene ROC curves",
        fontsize=15,
        fontweight="bold",
        y=0.975,
    )
    labels = cohort_data["labels"]
    cohort_results = results.loc[results["cohort"] == cohort].set_index("gene")
    for ax, gene in zip(axes.flat, GENES):
        gene_data = cohort_data["gene_data"][gene]
        ax.plot(
            [0, 1],
            [0, 1],
            linestyle=":",
            color="#64748B",
            linewidth=1.0,
            zorder=0,
            clip_on=False,
        )
        if gene_data["detected"]:
            _draw_curve(ax, labels, gene, gene_data)
            result = cohort_results.loc[gene]
            suffix = "" if gene_data["concordant"] else " · discordant"
            title = f"{gene}: {_format_auc(result)}{suffix}"
        else:
            title = f"{gene}: not detected"
            ax.text(
                0.5,
                0.5,
                "Not detected",
                transform=ax.transAxes,
                ha="center",
                va="center",
                color="#64748B",
                fontsize=10,
            )
        ax.set_title(title, fontsize=8.7)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(color="#E2E8F0", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_xlabel("False-positive rate", fontsize=8)
        ax.set_ylabel("True-positive rate", fontsize=8)
    fig.text(
        0.5,
        0.04,
        "Positive class = LV / long-term PD; scores are TMM-normalized log2-CPM; "
        "solid = direction-concordant, dashed = discordant. "
        "Concordance follows edgeR log2FC sign, not AUC; CIs use 2,000 stratified bootstrap resamples.",
        ha="center",
        fontsize=8,
        color="#475569",
    )
    stem = f"{OUTPUT_STEM}_small_multiples_{cohort}"
    fig.savefig(RESULTS_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS_DIR / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(PANEL_PATH)
    consensus = panel.loc[panel["Consensus_At_Least_2_of_3"] == 1].set_index("Gene")
    if set(consensus.index) != set(GENES):
        raise ValueError(f"Current consensus panel differs from expected genes: {sorted(consensus.index)}")
    panel_directions = {
        gene: int(np.sign(float(consensus.loc[gene, "Discovery_log2FC"])))
        for gene in GENES
    }
    if any(direction == 0 for direction in panel_directions.values()):
        raise ValueError("A consensus gene has a zero discovery log2 fold change; direction is undefined.")

    log_lines = [
        "ROC score provenance:",
        "  Input values are the existing edger_logCPM.csv matrices.",
        "  The edgeR scripts create these as cpm(calcNormFactors(DGEList(counts)), log=TRUE, prior.count=2).",
        "  calcNormFactors uses edgeR's default TMM method; therefore the ROC scores are TMM-normalized log2-CPM.",
        "  Positive class: LV (GSE248762) or long-term PD (GSE130888); scores were not direction-flipped.",
        "  Discovery directions come from Discovery_log2FC in ml_39_gene_results.csv.",
        "  Concordance compares the sign of the validation edgeR LV_vs_SV logFC with discovery.",
        "  CIs: percentile intervals from 2,000 stratified bootstrap resamples; seed 42.",
        "  DeLong CIs: single-marker DeLong placement variance with normal 95% Wald limits, clipped to [0, 1].",
        "  Ties: number of positive-negative sample pairs with exactly equal log2-CPM; sklearn AUC gives each tie half credit.",
        "  CSV also reports unique score values occurring more than once and the number of observations with tied scores.",
    ]
    cohorts = {
        cohort: _load_cohort(cohort, settings, panel_directions, log_lines)
        for cohort, settings in COHORTS.items()
    }
    results = _calculate_gene_rows(cohorts)
    for column in (
        "tied_positive_negative_pairs",
        "n_tied_score_values",
        "n_tied_score_observations",
    ):
        results[column] = pd.array(results[column], dtype="Int64")
    results.to_csv(RESULTS_DIR / f"{OUTPUT_STEM}_auc.csv", index=False, float_format="%.6f")
    _make_overview(cohorts, results)
    for cohort, cohort_data in cohorts.items():
        _make_small_multiples(cohort, cohort_data, results)

    log_lines.append("Per-gene results:")
    for row in results.itertuples(index=False):
        if row.detected:
            log_lines.append(
                f"  {row.cohort} {row.gene}: AUC={row.auc:.6f}, "
                f"95% CI=[{row.ci_low:.6f}, {row.ci_high:.6f}], "
                f"DeLong 95% CI=[{row.auc_delong_ci_low:.6f}, {row.auc_delong_ci_high:.6f}], "
                f"concordant={row.concordant}, cross-class ties={row.tied_positive_negative_pairs}, "
                f"tied score values={row.n_tied_score_values}, "
                f"observations in tied score values={row.n_tied_score_observations}"
            )
        else:
            log_lines.append(f"  {row.cohort} {row.gene}: not detected")
    log_text = "\n".join(log_lines) + "\n"
    log_path = RESULTS_DIR / f"{OUTPUT_STEM}_analysis.log"
    log_path.write_text(log_text, encoding="utf-8")
    print(log_text, end="")
    print(f"Saved AUC table, overview and small-multiple figures under {RESULTS_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
