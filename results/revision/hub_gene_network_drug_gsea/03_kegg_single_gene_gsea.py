"""Run exploratory single-gene KEGG preranked GSEA on GSE125498."""

from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path

import gseapy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
EXPR_PATH = ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv"
META_PATH = ROOT / "results" / "tables" / "GSE125498_sample_metadata.csv"
GENES = ["FLT3LG", "TNFSF15", "LTB", "ADAM19", "SERPINA10", "EBI3", "CXCL14"]
SEED = 20261013
PERMUTATIONS = 1_000
MIN_SIZE = 15
MAX_SIZE = 500
FDR_CUTOFF = 0.05
KEGG = "https://rest.kegg.jp"


def fetch_kegg_sets() -> tuple[dict[str, set[str]], dict[str, str], str, str]:
    pathway_response = requests.get(f"{KEGG}/list/pathway/hsa", timeout=90)
    pathway_response.raise_for_status()
    pathway_names = {
        line.split("\t", 1)[0].removeprefix("path:"): line.split("\t", 1)[1]
        for line in pathway_response.text.splitlines()
        if "\t" in line
    }

    gene_response = requests.get(f"{KEGG}/list/hsa", timeout=120)
    gene_response.raise_for_status()
    gene_names: dict[str, str] = {}
    for line in gene_response.text.splitlines():
        if "\t" not in line:
            continue
        fields = line.split("\t")
        if len(fields) >= 4:
            gene_id, _, _, description = fields[:4]
        elif len(fields) == 3:
            gene_id, _, description = fields
        else:
            continue
        # KEGG gene records place the official symbol first in the description.
        symbol = description.split(",", 1)[0].split(";", 1)[0].strip()
        if symbol:
            gene_names[gene_id] = symbol

    links_response = requests.get(f"{KEGG}/link/pathway/hsa", timeout=120)
    links_response.raise_for_status()
    gene_sets: dict[str, set[str]] = {name: set() for name in pathway_names.values()}
    pathway_id_to_name = pathway_names
    for line in links_response.text.splitlines():
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        gene_id = parts[0].removeprefix("hsa:")
        pathway_id = parts[1].removeprefix("path:")
        symbol = gene_names.get(f"hsa:{gene_id}")
        pathway_name = pathway_id_to_name.get(pathway_id)
        if symbol and pathway_name:
            gene_sets[pathway_name].add(symbol)

    info_response = requests.get(f"{KEGG}/info/kegg", timeout=60)
    info_response.raise_for_status()
    pathway_release = "not reported"
    gene_release = "not reported"
    for line in info_response.text.splitlines():
        fields = line.split()
        if fields and len(fields) >= 3:
            if fields[0] == "pathway":
                pathway_release = fields[-1]
            elif fields[0] == "genes":
                gene_release = fields[-1]
    gene_sets = {name: members for name, members in gene_sets.items() if members}
    if not gene_sets:
        raise RuntimeError("KEGG REST returned no human pathway gene sets.")
    return gene_sets, pathway_names, pathway_release, gene_release


def normalized_result(
    result_table: pd.DataFrame,
    gene: str,
    gene_sets: dict[str, set[str]],
    ranked_genes: set[str],
) -> pd.DataFrame:
    if result_table.empty:
        return pd.DataFrame(
            columns=[
                "gene",
                "pathway",
                "NES",
                "p_value",
                "FDR",
                "gene_set_size",
                "leading_edge_genes",
                "empirical_floor_applied",
                "passes_fdr_0_05",
            ]
        )
    table = result_table.copy()
    table["gene"] = gene
    table = table.rename(
        columns={
            "Term": "pathway",
            "NES": "NES",
            "NOM p-val": "p_value",
            "FDR q-val": "FDR",
            "Lead_genes": "leading_edge_genes",
            "Tag %": "tag_fraction",
            "Gene %": "gene_fraction",
        }
    )
    table["gene_set_size"] = table["pathway"].map(
        lambda pathway: len(gene_sets.get(pathway, set()) & ranked_genes)
    )
    if "leading_edge_genes" not in table:
        table["leading_edge_genes"] = ""
    table["leading_edge_genes"] = table["leading_edge_genes"].fillna("").astype(str)
    table["empirical_floor_applied"] = (
        pd.to_numeric(table["p_value"], errors="coerce").eq(0)
        | pd.to_numeric(table["FDR"], errors="coerce").eq(0)
    )
    empirical_floor = 1.0 / (PERMUTATIONS + 1)
    table["p_value"] = pd.to_numeric(table["p_value"], errors="coerce").clip(
        lower=empirical_floor
    )
    table["FDR"] = pd.to_numeric(table["FDR"], errors="coerce").clip(
        lower=empirical_floor
    )
    table["passes_fdr_0_05"] = pd.to_numeric(table["FDR"], errors="coerce") < FDR_CUTOFF
    keep = [
        "gene",
        "pathway",
        "NES",
        "p_value",
        "FDR",
        "gene_set_size",
        "leading_edge_genes",
        "empirical_floor_applied",
        "passes_fdr_0_05",
    ]
    return table[[column for column in keep if column in table.columns]]


def main() -> None:
    expression_df = pd.read_csv(EXPR_PATH)
    metadata = pd.read_csv(META_PATH)
    if "Gene" not in expression_df or not {"sample_id", "stage_binary"}.issubset(
        metadata.columns
    ):
        raise ValueError("Discovery expression/metadata columns do not match project conventions.")
    if expression_df["Gene"].duplicated().any():
        raise ValueError("Expected unique gene symbols in the discovery expression matrix.")
    samples = metadata["sample_id"].astype(str).tolist()
    if len(samples) != 33 or len(set(samples)) != 33:
        raise ValueError(f"Expected 33 unique discovery samples, found {len(samples)}.")
    if not set(samples).issubset(expression_df.columns):
        raise ValueError("Discovery matrix does not contain all metadata samples.")
    matrix = expression_df.set_index("Gene")[samples].astype(float)
    if not np.isfinite(matrix.to_numpy()).all():
        raise ValueError("Discovery expression matrix contains non-finite values.")
    absent = sorted(set(GENES) - set(matrix.index))
    if absent:
        raise ValueError(f"Hub genes absent from discovery matrix: {absent}")

    gene_sets, pathway_names, pathway_release, gene_release = fetch_kegg_sets()
    pd.DataFrame(
        [
            {"pathway": name, "gene_set_size": len(members), "genes": ";".join(sorted(members))}
            for name, members in gene_sets.items()
        ]
    ).to_csv(OUT / "kegg_human_gene_sets.csv", index=False)

    all_results: list[pd.DataFrame] = []
    summary_rows = []
    for gene in GENES:
        scores = matrix.apply(
            lambda values: float(spearmanr(matrix.loc[gene].to_numpy(), values.to_numpy()).statistic),
            axis=1,
        )
        ranking = (
            pd.DataFrame({"gene": scores.index, "spearman_rho": scores.values})
            .dropna()
            .sort_values(["spearman_rho", "gene"], ascending=[False, True])
        )
        ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
        ranking.to_csv(OUT / f"kegg_spearman_ranking_{gene}.csv", index=False)

        prerank_result = gseapy.prerank(
            rnk=ranking[["gene", "spearman_rho"]],
            gene_sets=gene_sets,
            min_size=MIN_SIZE,
            max_size=MAX_SIZE,
            permutation_num=PERMUTATIONS,
            seed=SEED,
            threads=8,
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        table = normalized_result(
            prerank_result.res2d, gene, gene_sets, set(ranking["gene"])
        )
        table.to_csv(OUT / f"kegg_gsea_{gene}.csv", index=False)
        all_results.append(table)
        n_sig = int(table["passes_fdr_0_05"].sum()) if not table.empty else 0
        print(f"Processed {gene}: found {n_sig} significant KEGG pathways (FDR < 0.05)", flush=True)
        summary_rows.append(
            {
                "gene": gene,
                "samples": len(samples),
                "ranked_genes": len(ranking),
                "KEGG_pathways_tested": len(table),
                "significant_pathways_FDR_lt_0_05": n_sig,
                "min_gene_set_size": MIN_SIZE,
                "max_gene_set_size": MAX_SIZE,
                "permutations": PERMUTATIONS,
                "seed": SEED,
                "status": "significant pathways found" if n_sig else "no pathway passes FDR < 0.05",
            }
        )

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUT / "kegg_gsea_all_genes.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "kegg_gsea_summary.csv", index=False)

    plot_rows = combined[combined["passes_fdr_0_05"]].copy()
    selected_rows = []
    for gene in GENES:
        one_gene = plot_rows[plot_rows["gene"] == gene]
        up = one_gene[one_gene["NES"] > 0].nlargest(5, "NES")
        down = one_gene[one_gene["NES"] < 0].nsmallest(5, "NES")
        selected_rows.extend([up, down])
    figure_data = (
        pd.concat(selected_rows, ignore_index=True)
        if selected_rows
        else pd.DataFrame(columns=combined.columns)
    )
    figure_data.to_csv(OUT / "kegg_gsea_figure_data.csv", index=False)
    fig, axes = plt.subplots(len(GENES), 1, figsize=(13, max(12, 2.3 * len(GENES))))
    for ax, gene in zip(np.atleast_1d(axes), GENES):
        data = figure_data[figure_data["gene"] == gene].sort_values("NES")
        if data.empty:
            ax.text(0.5, 0.5, "No KEGG pathway passes FDR < 0.05", ha="center", va="center")
            ax.set_yticks([])
        else:
            colors = ["#3B82A0" if value > 0 else "#C05A47" for value in data["NES"]]
            ax.barh(data["pathway"], data["NES"], color=colors)
            ax.axvline(0, color="#333333", lw=0.8)
            ax.tick_params(axis="y", labelsize=7)
        ax.set_title(gene, loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel("Normalized enrichment score (NES)")
    fig.suptitle(
        "Exploratory single-gene KEGG preranked GSEA\n"
        "GSE125498; significant terms only (FDR < 0.05)",
        y=1.0,
    )
    fig.tight_layout()
    fig.savefig(OUT / "kegg_single_gene_gsea.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "kegg_single_gene_gsea.pdf", bbox_inches="tight")
    plt.close(fig)

    metadata_out = {
        "query_date": date.today().isoformat(),
        "expression_cohort": "GSE125498",
        "sample_count": len(samples),
        "LPD_sample_count": int((metadata["stage_binary"] == 1).sum()),
        "SPD_sample_count": int((metadata["stage_binary"] == 0).sum()),
        "expression_matrix": str(EXPR_PATH.relative_to(ROOT)),
        "rank_metric": "Spearman correlation of each hub gene against every gene across 33 samples; tied values are retained",
        "pathway_database": "KEGG human pathways via KEGG REST",
        "KEGG_pathway_database_updated": pathway_release,
        "KEGG_gene_database_updated": gene_release,
        "pathway_count_downloaded": len(pathway_names),
        "pathway_count_with_genes": len(gene_sets),
        "gsea_software": "GSEApy",
        "gseapy_version": gseapy.__version__,
        "method_note": "GSEApy prerank used because fgsea/clusterProfiler are unavailable in this Python-only environment. Many genes have tied Spearman values; GSEApy notes tie ordering may affect results.",
        "empirical_p_fdr_zero_handling": "GSEApy zero P/FDR estimates are floored to 1/(permutations + 1); empirical_floor_applied marks affected rows.",
        "permutations": PERMUTATIONS,
        "min_gene_set_size": MIN_SIZE,
        "max_gene_set_size": MAX_SIZE,
        "fdr_cutoff": FDR_CUTOFF,
        "seed": SEED,
        "query_genes": GENES,
    }
    (OUT / "kegg_gsea_metadata.json").write_text(
        json.dumps(metadata_out, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata_out, indent=2))


if __name__ == "__main__":
    main()
