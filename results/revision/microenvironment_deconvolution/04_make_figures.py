"""Create boxplots and sample heatmaps for both signature variants."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns


OUT = Path(__file__).resolve().parent
VARIANTS = ("full_signature", "overlap_excluded_7gene")
GROUP_COLORS = {"SPD": "#2878A5", "LPD": "#C84B31"}


def main() -> None:
    scores = pd.read_csv(OUT / "sample_population_scores_long.csv")
    tests = pd.read_csv(OUT / "group_comparisons.csv")
    matrices: dict[str, pd.DataFrame] = {}
    for variant in VARIANTS:
        frame = scores[scores["signature_variant"] == variant].copy()
        matrix = pd.read_csv(
            OUT / f"sample_population_score_matrix_{variant}.csv", index_col="sample_id"
        )
        matrices[variant] = matrix

        populations = list(frame["population"].drop_duplicates())
        ncols = 3
        nrows = (len(populations) + ncols - 1) // ncols
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(14, 3.2 * nrows), squeeze=False, constrained_layout=True
        )
        for axis, population in zip(axes.flat, populations):
            sub = frame[frame["population"] == population]
            comparison = tests[
                (tests["signature_variant"] == variant)
                & (tests["population"] == population)
            ].iloc[0]
            sns.boxplot(
                data=sub, x="group", y="score_mean_z", order=["SPD", "LPD"],
                hue="group", hue_order=["SPD", "LPD"], palette=GROUP_COLORS,
                dodge=False, legend=False, width=0.58, fliersize=0, ax=axis,
            )
            sns.stripplot(
                data=sub, x="group", y="score_mean_z", order=["SPD", "LPD"],
                hue="group", palette=GROUP_COLORS, jitter=0.18, size=3.2,
                alpha=0.75, legend=False, ax=axis,
            )
            axis.set_title(
                f"{population}\nBH FDR = {comparison['BH_FDR_within_signature_variant']:.3g}"
            )
            axis.set_xlabel("")
            axis.set_ylabel("Mean marker z-score")
        for axis in axes.flat[len(populations):]:
            axis.remove()
        fig.suptitle(
            "Relative cell-population marker scores: " + variant.replace("_", " "),
            fontsize=15,
        )
        for ext in ("png", "pdf"):
            fig.savefig(OUT / f"cell_population_boxplots_{variant}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)

        ordered = matrix.assign(
            _group_order=matrix["group"].map({"SPD": 0, "LPD": 1})
        ).sort_values(["_group_order"]).drop(columns="_group_order")
        annotation = ordered.pop("group")
        fig_width = max(11, 0.37 * ordered.shape[0])
        fig_height = max(7, 0.45 * ordered.shape[1])
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
        sns.heatmap(
            ordered.T,
            cmap="vlag",
            center=0,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"label": "Mean marker z-score"},
            ax=ax,
        )
        for index, group in enumerate(annotation):
            ax.axvline(index, color=GROUP_COLORS[group], linewidth=0.5, alpha=0.85)
        ax.legend(
            handles=[
                Patch(facecolor=GROUP_COLORS[group], label=group)
                for group in ("SPD", "LPD")
            ],
            title="Sample group",
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
        )
        ax.set_title(
            "Sample-by-population relative scores: " + variant.replace("_", " ")
        )
        ax.set_xlabel("Samples (sorted: SPD, then LPD)")
        ax.set_ylabel("Population")
        ax.tick_params(axis="x", labelrotation=90, labelsize=7)
        for ext in ("png", "pdf"):
            fig.savefig(OUT / f"sample_population_heatmap_{variant}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
