"""Create a publication-ready summary from the frozen seven-gene DSigDB results."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#333333",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)

OUT = Path(__file__).resolve().parent
INPUT = OUT / "validated7_dsigdb_enrichment_full.csv"
FDR_CUTOFF = 0.05
NOMINAL_CUTOFF = 0.05
MIN_OVERLAP = 2


def main() -> None:
    results = pd.read_csv(INPUT)
    required = {
        "dsigdb_term",
        "overlapping_genes",
        "overlap_size",
        "p_value",
        "adjusted_p_value",
        "combined_score",
    }
    missing = required.difference(results.columns)
    if missing:
        raise ValueError(f"DSigDB result table is missing columns: {sorted(missing)}")

    if results.empty:
        raise ValueError("DSigDB result table contains no enrichment terms.")

    candidates = results.loc[
        (results["p_value"] < NOMINAL_CUTOFF)
        & (results["overlap_size"] >= MIN_OVERLAP)
    ].copy()
    candidates = candidates.sort_values(
        ["p_value", "adjusted_p_value", "dsigdb_term"]
    ).reset_index(drop=True)
    candidates.to_csv(OUT / "validated7_dsigdb_multi_gene_nominal_candidates.csv", index=False)

    fdr_significant = int((results["adjusted_p_value"] < FDR_CUTOFF).sum())
    min_q = float(results["adjusted_p_value"].min())
    if candidates.empty:
        raise ValueError("No nominally enriched DSigDB terms overlap at least two query genes.")

    plot_data = candidates.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_data))
    neg_log10_q = -np.log10(plot_data["adjusted_p_value"].clip(lower=1e-300))
    colors = np.where(
        plot_data["adjusted_p_value"] < FDR_CUTOFF, "#C44E52", "#4C78A8"
    )
    sizes = 55 + 48 * plot_data["overlap_size"].to_numpy()

    fig_height = max(7.2, 0.46 * len(plot_data) + 2.2)
    fig, ax = plt.subplots(figsize=(13.2, fig_height), dpi=300)
    ax.scatter(
        neg_log10_q,
        y,
        s=sizes,
        c=colors,
        edgecolors="white",
        linewidths=0.7,
        zorder=3,
    )
    ax.axvline(
        -np.log10(FDR_CUTOFF),
        color="#C44E52",
        linestyle="--",
        linewidth=1.2,
        label="BH-FDR q = 0.05",
    )

    labels = [
        f"{row.dsigdb_term}\nGenes: {row.overlapping_genes.replace(';', ', ')}"
        for row in plot_data.itertuples(index=False)
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.1)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xlabel(r"$-\log_{10}$(BH-adjusted P)", fontsize=10.5, fontweight="bold")
    ax.set_title(
        "Exploratory DSigDB perturbation-signature enrichment for seven hub genes\n"
        f"No signatures pass BH-FDR < 0.05 (183 tested; minimum q = {min_q:.3f})",
        fontsize=12,
        fontweight="bold",
        pad=13,
    )
    ax.set_xlim(0, max(1.55, float(neg_log10_q.max()) + 0.63))
    ax.set_ylim(-0.7, len(plot_data) - 0.3)
    ax.grid(axis="x", linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)

    for idx, row in plot_data.iterrows():
        ax.annotate(
            f"P={row['p_value']:.3g}; q={row['adjusted_p_value']:.3f}",
            (neg_log10_q.iloc[idx], y[idx]),
            xytext=(8, 0),
            textcoords="offset points",
            va="center",
            fontsize=7.4,
            color="#333333",
        )

    fig.text(
        0.5,
        0.012,
        "Shown: nominal P < 0.05 and overlap with at least two of seven query genes. "
        "DSigDB chemical/cell-line signatures are hypothesis-generating associations, "
        "not evidence of therapeutic efficacy.",
        ha="center",
        va="bottom",
        fontsize=8.2,
        color="#444444",
        wrap=True,
    )
    fig.subplots_adjust(left=0.39, right=0.97, top=0.92, bottom=0.075)
    fig.savefig(OUT / "validated7_dsigdb_publication.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "validated7_dsigdb_publication.pdf", bbox_inches="tight")
    plt.close(fig)

    print(
        f"Saved publication plot and candidate table: {len(candidates)} multi-gene "
        f"nominal terms; {fdr_significant}/{len(results)} terms pass BH-FDR < "
        f"{FDR_CUTOFF}; minimum q={min_q:.6g}."
    )


if __name__ == "__main__":
    main()
