"""Plot classic running-enrichment curves for the seven validated hub genes."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
GENES = ["TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "ECM1", "CST7"]
PATHWAY_SETS = OUT / "kegg_human_gene_sets.csv"
MAX_PATHWAYS_PER_HUB = 3
FDR_CUTOFF = 0.05
COLORS = ["#D62728", "#1F77B4", "#2CA02C"]

matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.75,
        "axes.edgecolor": "#333333",
        "figure.facecolor": "white",
    }
)


def running_profile(
    ranked_genes: list[str], scores: np.ndarray, members: set[str]
) -> tuple[float, list[int], np.ndarray]:
    hit_mask = np.isin(ranked_genes, list(members))
    hit_indices = np.flatnonzero(hit_mask).tolist()
    hit_count = len(hit_indices)
    if hit_count == 0 or hit_count == len(ranked_genes):
        raise ValueError("A pathway must overlap some, but not all, ranked genes.")
    hit_weights = np.abs(scores) * hit_mask
    if hit_weights.sum() == 0:
        raise ValueError("A pathway's hit genes have zero total ranking weight.")
    increments = hit_weights / hit_weights.sum()
    increments[~hit_mask] = -1.0 / (len(ranked_genes) - hit_count)
    running = np.cumsum(increments)
    maximum, minimum = float(running.max()), float(running.min())
    es = maximum if abs(maximum) > abs(minimum) else minimum
    return es, hit_indices, running


def load_pathway_sets() -> dict[str, set[str]]:
    table = pd.read_csv(PATHWAY_SETS)
    if table["pathway"].duplicated().any():
        raise ValueError("Cached KEGG pathway names must be unique.")
    return {
        str(row.pathway): set(str(row.genes).split(";")) - {""}
        for row in table.itertuples(index=False)
    }


def main() -> None:
    pathway_sets = load_pathway_sets()
    selected_rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []
    hub_data: list[dict[str, object]] = []
    for gene in GENES:
        ranking_path = OUT / f"validated7_kegg_spearman_ranking_{gene}.csv"
        results_path = OUT / f"validated7_kegg_gsea_{gene}.csv"
        ranking = pd.read_csv(ranking_path)
        results = pd.read_csv(results_path)
        if not {"gene", "spearman_rho"}.issubset(ranking.columns):
            raise ValueError(f"Unexpected ranking columns for {gene}.")
        if ranking["gene"].duplicated().any():
            raise ValueError(f"Ranking contains duplicate genes for {gene}.")
        ranking = ranking.dropna(subset=["gene", "spearman_rho"])
        ranking["gene"] = ranking["gene"].astype(str)
        ranking["spearman_rho"] = pd.to_numeric(
            ranking["spearman_rho"], errors="raise"
        )
        ranked_genes = ranking["gene"].tolist()
        scores = ranking["spearman_rho"].to_numpy(dtype=float)
        measured = set(ranked_genes)
        significant = results.loc[
            pd.to_numeric(results["FDR"], errors="coerce") < FDR_CUTOFF
        ].copy()
        significant["absolute_NES"] = pd.to_numeric(
            significant["NES"], errors="raise"
        ).abs()
        chosen = significant.sort_values(
            ["absolute_NES", "FDR", "pathway"],
            ascending=[False, True, True],
            kind="mergesort",
        ).head(MAX_PATHWAYS_PER_HUB)
        if chosen.empty:
            raise RuntimeError(f"No FDR-significant pathways available for {gene}.")

        curves: list[dict[str, object]] = []
        for pathway_index, result in enumerate(chosen.itertuples(index=False)):
            pathway = str(result.pathway)
            if pathway not in pathway_sets:
                raise KeyError(f"Cached pathway membership not found: {pathway}")
            members = pathway_sets[pathway] & measured
            expected_size = int(result.gene_set_size)
            if len(members) != expected_size:
                raise RuntimeError(
                    f"{gene}/{pathway}: cached set overlap {len(members)} "
                    f"does not match saved GSEA set size {expected_size}."
                )
            es, hits, running = running_profile(ranked_genes, scores, members)
            nes = float(result.NES)
            if np.sign(es) != np.sign(nes):
                raise RuntimeError(
                    f"Running ES direction does not match saved NES for {gene}/{pathway}."
                )
            color = COLORS[pathway_index]
            curve = {
                "pathway": pathway,
                "short_name": pathway.split(" - Homo sapiens")[0],
                "NES": nes,
                "FDR": float(result.FDR),
                "ES": es,
                "set_size": len(members),
                "ranked_genes": ranked_genes,
                "scores": scores,
                "hits": hits,
                "running": running,
                "color": color,
            }
            curves.append(curve)
            selected_rows.append(
                {
                    "hub_gene": gene,
                    "pathway": pathway,
                    "NES": nes,
                    "FDR": float(result.FDR),
                    "running_ES": es,
                    "gene_set_size": len(members),
                    "ranking_source": ranking_path.name,
                }
            )
            hit_set = set(hits)
            profile_rows.extend(
                {
                    "hub_gene": gene,
                    "pathway": pathway,
                    "rank": rank,
                    "ranked_gene": ranked_genes[rank - 1],
                    "spearman_rho": float(scores[rank - 1]),
                    "running_enrichment_score": float(running[rank - 1]),
                    "is_pathway_hit": int(rank - 1 in hit_set),
                }
                for rank in range(1, len(ranked_genes) + 1)
            )
        hub_data.append(
            {"gene": gene, "ranking": ranking, "curves": curves}
        )

    selected_table = pd.DataFrame(selected_rows)
    selected_table.to_csv(
        OUT / "validated7_kegg_running_gsea_selected_pathways.csv", index=False
    )
    pd.DataFrame(profile_rows).to_csv(
        OUT / "validated7_kegg_running_gsea_profiles.csv", index=False
    )

    fig = plt.figure(figsize=(13.2, 20.5), dpi=300)
    outer = gridspec.GridSpec(
        len(GENES),
        1,
        figure=fig,
        left=0.10,
        right=0.73,
        top=0.955,
        bottom=0.055,
        hspace=0.34,
    )
    for row_index, data in enumerate(hub_data):
        inner = outer[row_index].subgridspec(
            3, 1, height_ratios=[3.0, 0.48, 0.82], hspace=0.035
        )
        ax_es = fig.add_subplot(inner[0])
        ax_hits = fig.add_subplot(inner[1], sharex=ax_es)
        ax_rank = fig.add_subplot(inner[2], sharex=ax_es)

        ranking = data["ranking"]
        ranked_genes = ranking["gene"].tolist()
        scores = ranking["spearman_rho"].to_numpy(dtype=float)
        ranks = np.arange(1, len(ranked_genes) + 1)

        for curve_index, curve in enumerate(data["curves"]):
            ax_es.plot(
                ranks,
                curve["running"],
                color=curve["color"],
                linewidth=1.35,
                label=(
                    f"{curve['short_name']} "
                    f"(NES={curve['NES']:+.2f}; FDR={curve['FDR']:.3g})"
                ),
            )
            hit_row = len(data["curves"]) - curve_index
            ax_hits.vlines(
                np.asarray(curve["hits"]) + 1,
                hit_row - 0.33,
                hit_row + 0.33,
                color=curve["color"],
                linewidth=0.6,
            )

        ax_es.axhline(0, color="#666666", linewidth=0.6)
        ax_es.set_xlim(1, len(ranked_genes))
        ax_es.set_ylim(bottom=min(0.0, ax_es.get_ylim()[0]))
        ax_es.set_ylabel(
            f"{data['gene']}\nRunning ES", fontsize=8.4, fontweight="bold"
        )
        ax_es.grid(True, linestyle=":", color="#CBD5E1", alpha=0.55)
        ax_es.tick_params(axis="x", labelbottom=False, labelsize=7)

        ax_hits.set_ylim(0.4, len(data["curves"]) + 0.6)
        ax_hits.set_yticks([])
        ax_hits.set_ylabel("Hits", fontsize=7.3, fontweight="bold", labelpad=8)
        ax_hits.spines["left"].set_visible(False)
        ax_hits.spines["right"].set_visible(False)
        ax_hits.tick_params(axis="x", labelbottom=False)

        ax_rank.fill_between(
            ranks, scores, 0, color="#B8C7D9", alpha=0.72, linewidth=0
        )
        ax_rank.plot(ranks, scores, color="#536B87", linewidth=0.55)
        ax_rank.axhline(0, color="#536B87", linewidth=0.65)
        ax_rank.set_ylim(-1.05, 1.05)
        ax_rank.set_ylabel("Spearman\nrho", fontsize=7.2, fontweight="bold")
        ax_rank.set_xticks(np.linspace(1, len(ranked_genes), 5, dtype=int))
        ax_rank.tick_params(axis="both", labelsize=6.5)
        ax_rank.grid(True, axis="y", linestyle=":", color="#CBD5E1", alpha=0.5)
        if row_index == len(hub_data) - 1:
            ax_rank.set_xlabel(
                f"Rank in hub-specific ordered gene list ({len(ranked_genes):,} genes)",
                fontsize=8.4,
                fontweight="bold",
            )
        else:
            ax_rank.tick_params(axis="x", labelbottom=False)

        ax_es.legend(
            loc="center left",
            bbox_to_anchor=(1.005, 0.5),
            frameon=False,
            fontsize=7.2,
            handlelength=1.8,
            borderaxespad=0,
        )

    fig.suptitle(
        "Classic running-enrichment plots for the seven validated hub genes\n"
        "GSE125498 discovery cohort (n=33); top FDR-significant KEGG pathways "
        "per hub; 1,000 permutations",
        fontsize=13,
        fontweight="bold",
        y=0.99,
    )
    fig.savefig(
        OUT / "validated7_kegg_running_gsea.png", dpi=300, bbox_inches="tight"
    )
    fig.savefig(
        OUT / "validated7_kegg_running_gsea.pdf", bbox_inches="tight"
    )
    plt.close(fig)

    metadata = {
        "analysis": "Hub-specific preranked KEGG running-enrichment visualization",
        "cohort": "GSE125498",
        "sample_count": 33,
        "rank_metric": "Spearman correlation of each gene with the given hub across all samples",
        "pathway_set_source": str(PATHWAY_SETS.relative_to(ROOT)),
        "gsea_results_source": "validated7_kegg_gsea_<HUB>.csv",
        "ranked_gene_sources": "validated7_kegg_spearman_ranking_<HUB>.csv",
        "pathways_displayed_per_hub": MAX_PATHWAYS_PER_HUB,
        "selection": (
            "Three pathways with largest absolute NES among those passing "
            "within-hub BH-FDR < 0.05 (all significant pathways shown when fewer than three)"
        ),
        "permutations": 1000,
        "fdr_cutoff": FDR_CUTOFF,
        "hub_genes": GENES,
        "interpretation": (
            "Each plot shows pathway enrichment along that hub's gene-correlation "
            "ranking, not an LPD-versus-SPD pathway effect."
        ),
    }
    (OUT / "validated7_kegg_running_gsea_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(
        f"Saved seven hub-specific panels with "
        f"{len(selected_table)} validated pathway curves; all pathway hit counts "
        "and running ES directions match saved GSEA results."
    )


if __name__ == "__main__":
    main()
