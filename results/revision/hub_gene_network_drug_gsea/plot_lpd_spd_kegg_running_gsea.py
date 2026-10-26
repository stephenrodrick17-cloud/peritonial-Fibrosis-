"""Render classic running-enrichment curves for LPD-versus-SPD KEGG GSEA."""

from __future__ import annotations

import json
import warnings
from datetime import date
from pathlib import Path

import gseapy
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
RANKING_PATH = ROOT / "results" / "paper_style" / "Fig2_limma_gene_results_unique_symbols.csv"
GENE_SETS_PATH = OUT / "kegg_human_gene_sets.csv"
PERMUTATIONS = 1_000
SEED = 20261016
MIN_SIZE = 15
MAX_SIZE = 500
FDR_CUTOFF = 0.05

matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#333333",
        "figure.facecolor": "white",
    }
)


def load_inputs() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    ranking = pd.read_csv(RANKING_PATH)
    required = {"Gene.symbol", "t"}
    if not required.issubset(ranking.columns):
        raise ValueError(f"Limma ranking lacks required columns: {sorted(required - set(ranking.columns))}")
    ranking = ranking.loc[:, ["Gene.symbol", "t"]].dropna()
    ranking["Gene.symbol"] = ranking["Gene.symbol"].astype(str)
    ranking["t"] = pd.to_numeric(ranking["t"], errors="raise")
    if ranking["Gene.symbol"].duplicated().any():
        raise ValueError("Expected one limma moderated t-statistic per gene symbol.")
    if not np.isfinite(ranking["t"].to_numpy()).all():
        raise ValueError("Ranking contains non-finite moderated t-statistics.")
    ranking = ranking.sort_values(
        ["t", "Gene.symbol"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)

    set_table = pd.read_csv(GENE_SETS_PATH)
    if not {"pathway", "genes"}.issubset(set_table.columns):
        raise ValueError("Cached KEGG gene-set table lacks pathway or genes columns.")
    gene_sets = {
        str(row.pathway): sorted(set(str(row.genes).split(";")) - {""})
        for row in set_table.itertuples(index=False)
    }
    if len(gene_sets) != len(set_table):
        raise ValueError("Cached KEGG pathway names are not unique.")
    return ranking, gene_sets


def calculate_running_profile(
    genes: list[str], scores: np.ndarray, members: list[str]
) -> tuple[float, list[int], np.ndarray]:
    hit_mask = np.isin(genes, members)
    hit_indices = np.flatnonzero(hit_mask).tolist()
    hit_count = len(hit_indices)
    if hit_count == 0 or hit_count == len(genes):
        raise ValueError("A pathway must overlap some, but not all, ranked genes.")
    hit_weights = np.abs(scores) * hit_mask
    if hit_weights.sum() == 0:
        raise ValueError("A pathway's hit genes have zero total ranking weight.")
    increments = hit_weights / hit_weights.sum()
    increments[~hit_mask] = -1.0 / (len(genes) - hit_count)
    running = np.cumsum(increments)
    maximum, minimum = float(running.max()), float(running.min())
    es = maximum if abs(maximum) > abs(minimum) else minimum
    return es, hit_indices, running


def main() -> None:
    ranking, gene_sets = load_inputs()
    ranking.to_csv(OUT / "lpd_spd_kegg_gsea_ranked_genes.csv", index=False)
    ranked_genes = ranking["Gene.symbol"].tolist()
    ranked_scores = ranking["t"].to_numpy(dtype=float)
    measured = set(ranked_genes)
    eligible_sets = {
        name: sorted(measured.intersection(members))
        for name, members in gene_sets.items()
        if MIN_SIZE <= len(measured.intersection(members)) <= MAX_SIZE
    }
    if not eligible_sets:
        raise RuntimeError("No cached KEGG pathways meet the configured size limits.")

    warnings.filterwarnings("ignore", category=FutureWarning)
    result = gseapy.prerank(
        rnk=ranking[["Gene.symbol", "t"]],
        gene_sets=eligible_sets,
        min_size=MIN_SIZE,
        max_size=MAX_SIZE,
        permutation_num=PERMUTATIONS,
        seed=SEED,
        threads=1,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    table = result.res2d.rename(
        columns={
            "Term": "pathway",
            "ES": "ES",
            "NES": "NES",
            "NOM p-val": "p_value",
            "FDR q-val": "FDR",
            "Lead_genes": "leading_edge_genes",
        }
    )
    table["set_size_in_ranking"] = table["pathway"].map(
        lambda pathway: len(eligible_sets[str(pathway)])
    )
    table = table.sort_values(["FDR", "p_value"], kind="mergesort").reset_index(drop=True)
    table.to_csv(OUT / "lpd_spd_kegg_gsea_all_pathways.csv", index=False)
    positive = table.loc[
        (pd.to_numeric(table["NES"], errors="coerce") > 0)
        & (pd.to_numeric(table["FDR"], errors="coerce") < FDR_CUTOFF)
    ].copy()
    selected = positive.head(6).copy()
    if selected.empty:
        raise RuntimeError("No positively enriched KEGG pathway passed BH-FDR < 0.05.")
    selected.to_csv(OUT / "lpd_spd_kegg_gsea_selected_pathways.csv", index=False)

    colors = ["#E31A8A", "#08A5E8", "#E9B500", "#148A2C", "#9A4510", "#E32636"]
    profiles: list[dict[str, object]] = []
    curves: list[tuple[pd.Series, list[int], np.ndarray, str]] = []
    for color, row in zip(colors, selected.itertuples(index=False)):
        pathway = str(row.pathway)
        members = eligible_sets[pathway]
        es, hits, running = calculate_running_profile(ranked_genes, ranked_scores, members)
        if len(hits) != int(row.set_size_in_ranking) or len(running) != len(ranked_genes):
            raise RuntimeError(f"Running profile dimensions do not match for {pathway}.")
        curves.append((row, hits, running, color))
        hit_set = set(hits)
        for rank, score in enumerate(running, start=1):
            profiles.append(
                {
                    "pathway": pathway,
                    "rank": rank,
                    "running_enrichment_score": score,
                    "is_hit": int((rank - 1) in hit_set),
                    "ranked_gene": ranked_genes[rank - 1],
                    "moderated_t": ranked_scores[rank - 1],
                }
            )
        if not np.isclose(es, float(row.ES), atol=0.015):
            raise RuntimeError(
                f"Running ES ({es:.4f}) disagrees with GSEA result "
                f"({float(row.ES):.4f}) for {pathway}."
            )
    pd.DataFrame(profiles).to_csv(OUT / "lpd_spd_kegg_gsea_running_profiles.csv", index=False)

    fig = plt.figure(figsize=(11.4, 8.4), dpi=300)
    gs = gridspec.GridSpec(
        3,
        1,
        height_ratios=[4.8, 1.15, 1.4],
        hspace=0.06,
        left=0.105,
        right=0.98,
        top=0.90,
        bottom=0.12,
    )
    ax_es = fig.add_subplot(gs[0])
    ax_hits = fig.add_subplot(gs[1], sharex=ax_es)
    ax_rank = fig.add_subplot(gs[2], sharex=ax_es)

    x = np.arange(1, len(ranked_genes) + 1)
    for pathway_index, (row, hits, running, color) in enumerate(curves):
        q_value = float(row.FDR)
        q_label = f"FDR<0.001" if q_value < 0.001 else f"FDR={q_value:.3g}"
        label = (
            f"{row.pathway.split(' - Homo sapiens')[0]} "
            f"(NES={float(row.NES):+.2f}; {q_label})"
        )
        ax_es.plot(x, running, color=color, linewidth=1.75, label=label)
        row_number = len(curves) - pathway_index
        ax_hits.vlines(np.asarray(hits) + 1, row_number - 0.34, row_number + 0.34,
                       color=color, linewidth=0.68)

    ax_es.axhline(0, color="#666666", linewidth=0.65)
    ax_es.set_ylabel("Running Enrichment Score", fontsize=10, fontweight="bold")
    ax_es.set_xlim(1, len(ranked_genes))
    ax_es.set_ylim(bottom=min(0.0, ax_es.get_ylim()[0]))
    ax_es.grid(True, linestyle=":", color="#CBD5E1", alpha=0.6)
    ax_es.legend(loc="upper right", frameon=True, fontsize=7.6, ncol=1)
    ax_es.tick_params(axis="x", labelbottom=False)

    ax_hits.set_ylim(0.4, len(curves) + 0.6)
    ax_hits.set_yticks([])
    ax_hits.set_ylabel("Hits", fontsize=9, fontweight="bold", rotation=90, labelpad=12)
    ax_hits.spines["left"].set_visible(False)
    ax_hits.spines["right"].set_visible(False)
    ax_hits.grid(False)
    ax_hits.tick_params(axis="x", labelbottom=False)

    ax_rank.fill_between(x, ranked_scores, 0, color="#B8C7D9", alpha=0.68, linewidth=0)
    ax_rank.plot(x, ranked_scores, color="#536B87", linewidth=0.65)
    ax_rank.axhline(0, color="#536B87", linewidth=0.8)
    ax_rank.set_ylabel("Ranked Metric\n(moderated t)", fontsize=9, fontweight="bold")
    ax_rank.set_xlabel(f"Rank in Ordered Dataset ({len(ranked_genes):,} genes)", fontsize=10, fontweight="bold")
    ax_rank.grid(True, axis="y", linestyle=":", color="#CBD5E1", alpha=0.6)

    fig.suptitle(
        "KEGG pathways enriched in LPD versus SPD\n"
        f"GSE125498 (LPD n=13, SPD n=20); positive NES; BH-FDR < {FDR_CUTOFF}; "
        f"{PERMUTATIONS:,} permutations",
        fontsize=13,
        fontweight="bold",
        y=0.97,
    )
    fig.savefig(OUT / "lpd_spd_kegg_running_gsea.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "lpd_spd_kegg_running_gsea.pdf", bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "analysis": "Preranked KEGG GSEA, LPD versus SPD",
        "cohort": "GSE125498",
        "groups": {"LPD": 13, "SPD": 20},
        "ranked_metric": "limma moderated t-statistic, positive means higher in LPD",
        "ranking_source": str(RANKING_PATH.relative_to(ROOT)),
        "pathway_sets_source": str(GENE_SETS_PATH.relative_to(ROOT)),
        "pathway_set_source_release": "KEGG human, cached in project on 2026-10-16",
        "software": f"GSEApy {gseapy.__version__}",
        "permutations": PERMUTATIONS,
        "seed": SEED,
        "gene_set_size": [MIN_SIZE, MAX_SIZE],
        "pathways_tested": len(table),
        "positive_fdr_significant_pathways": int(len(positive)),
        "displayed_pathways": selected["pathway"].tolist(),
        "query_date": date.today().isoformat(),
        "interpretation": "Discovery-cohort pathway association; exploratory, not external validation.",
    }
    (OUT / "lpd_spd_kegg_running_gsea_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(
        f"Saved running GSEA plot for {len(selected)} pathways; ranked "
        f"{len(ranked_genes):,} genes; tested {len(table)} pathways; "
        f"{len(positive)} positive pathways passed BH-FDR < 0.05."
    )


if __name__ == "__main__":
    main()
