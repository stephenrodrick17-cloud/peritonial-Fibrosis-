"""Build a STRING human physical-interaction network for Matrisome DEGs."""

from __future__ import annotations

import argparse
import csv
import io
import json
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = (
    ROOT
    / "results"
    / "revision"
    / "matrisome_volcano"
    / "GSE125498_matrisome_volcano_source_data.csv"
)
OUTPUT_DIR = ROOT / "results" / "revision" / "matrisome_string_ppi"
SPECIES = 9606
MIN_SCORE = 400
CALLER_IDENTITY = "peritoneal_fibrosis_matrisome_analysis"
API_ROOT = "https://string-db.org/api/tsv"
STRING_SYMBOL_OVERRIDES = {
    # The source table uses NAT6; its stable Ensembl gene ID ENSG00000243477
    # resolves to NAA80. STRING also maps the ambiguous alias NAT6 to HYAL1.
    "NAT6": "NAA80",
}


def post_string(endpoint: str, params: dict[str, str]) -> str:
    request = Request(
        f"{API_ROOT}/{endpoint}",
        data=urlencode(params).encode("utf-8"),
        headers={"User-Agent": f"{CALLER_IDENTITY}/1.0"},
    )
    with urlopen(request, timeout=90) as response:
        return response.read().decode("utf-8")


def read_tsv(payload: str, expected: set[str], endpoint: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(payload), delimiter="\t")
    if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
        raise ValueError(f"STRING {endpoint} response is missing expected columns.")
    return list(reader)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--network-type",
        choices=("physical", "functional"),
        default="functional",
        help="STRING network type (default: functional associations)",
    )
    args = parser.parse_args()
    network_type = args.network_type
    output_dir = OUTPUT_DIR / network_type

    if not SOURCE_PATH.is_file():
        raise FileNotFoundError(f"Missing MatrisomeAnalyzeR source table: {SOURCE_PATH}")

    source = pd.read_csv(SOURCE_PATH)
    required = {"Gene", "DEG_Status", "logFC", "P.Value"}
    if not required.issubset(source.columns):
        raise ValueError(f"Source table is missing columns: {sorted(required - set(source.columns))}")

    significant = source.loc[
        source["DEG_Status"].isin(["Upregulated", "Downregulated"])
    ].copy()
    if significant["Gene"].duplicated().any():
        raise ValueError("Significant MatrisomeAnalyzeR input contains duplicate gene symbols.")
    if len(significant) != 46:
        raise ValueError(f"Expected 46 MatrisomeAnalyzeR DEGs, found {len(significant)}.")
    direction_counts = significant["DEG_Status"].value_counts().to_dict()
    if direction_counts != {"Upregulated": 33, "Downregulated": 13}:
        raise ValueError(f"Expected 33 upregulated and 13 downregulated genes; found {direction_counts}.")

    genes = significant["Gene"].astype(str).tolist()
    query_symbols = {
        gene: STRING_SYMBOL_OVERRIDES.get(gene, gene) for gene in genes
    }
    identifiers = "\r\n".join(query_symbols.values())
    shared_params = {
        "identifiers": identifiers,
        "species": str(SPECIES),
        "caller_identity": CALLER_IDENTITY,
    }

    resolved_by_query: dict[str, dict[str, str]] = {}
    for gene in genes:
        query_symbol = query_symbols[gene]
        resolved = read_tsv(
            post_string("resolve", {**shared_params, "identifiers": query_symbol}),
            {"stringId", "preferredName", "ncbiTaxonId"},
            "resolve",
        )
        if not resolved:
            continue
        exact_matches = [
            row for row in resolved if row["preferredName"] == query_symbol
        ]
        resolved_by_query[gene] = exact_matches[0] if exact_matches else resolved[0]
    unresolved = sorted(set(genes) - set(resolved_by_query))
    if unresolved:
        raise ValueError(
            "STRING could not resolve the following input genes; refusing to silently "
            f"omit them: {', '.join(unresolved)}"
        )

    network_payload = post_string(
        "network",
        {
            **shared_params,
            "required_score": str(MIN_SCORE),
            "network_type": network_type,
        },
    )
    interactions = read_tsv(
        network_payload,
        {"preferredName_A", "preferredName_B", "score", "stringId_A", "stringId_B"},
        "network",
    )
    if not interactions:
        raise ValueError("STRING returned no physical interactions at the selected score threshold.")

    gene_by_string_id = {
        row["stringId"]: gene for gene, row in resolved_by_query.items()
    }
    mapped_interactions = []
    for edge in interactions:
        gene_a = gene_by_string_id.get(edge["stringId_A"])
        gene_b = gene_by_string_id.get(edge["stringId_B"])
        if gene_a is None or gene_b is None:
            raise ValueError("STRING returned an interaction outside the requested 46-gene set.")
        edge = dict(edge)
        edge["input_gene_A"] = gene_a
        edge["input_gene_B"] = gene_b
        edge["score"] = str(float(edge["score"]))
        mapped_interactions.append(edge)
    interactions = mapped_interactions

    output_dir.mkdir(parents=True, exist_ok=True)
    source.to_csv(output_dir / "matrisome_deg_source_table.csv", index=False)
    resolved_path = output_dir / "string_gene_mapping.csv"
    pd.DataFrame(
        [
            {
                "input_gene": gene,
                "query_symbol": query_symbols[gene],
                "string_id": resolved_by_query[gene]["stringId"],
                "string_preferred_name": resolved_by_query[gene]["preferredName"],
                "taxon_id": resolved_by_query[gene]["ncbiTaxonId"],
                "DEG_status": significant.set_index("Gene").loc[gene, "DEG_Status"],
            }
            for gene in genes
        ]
    ).to_csv(resolved_path, index=False)

    edge_path = output_dir / f"string_{network_type}_interactions.tsv"
    edge_fields = list(interactions[0].keys())
    with edge_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=edge_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(interactions)

    graph = nx.Graph()
    for _, row in significant.iterrows():
        graph.add_node(
            str(row["Gene"]),
            DEG_status=str(row["DEG_Status"]),
            logFC=float(row["logFC"]),
            p_value=float(row["P.Value"]),
        )
    for edge in interactions:
        gene_a = edge["input_gene_A"]
        gene_b = edge["input_gene_B"]
        score = float(edge["score"])
        if graph.has_edge(gene_a, gene_b):
            graph[gene_a][gene_b]["score"] = max(graph[gene_a][gene_b]["score"], score)
        else:
            graph.add_edge(gene_a, gene_b, score=score)

    nodes = []
    for gene, attrs in graph.nodes(data=True):
        nodes.append(
            {
                "gene": gene,
                "DEG_status": attrs["DEG_status"],
                "logFC": attrs["logFC"],
                "P_value": attrs["p_value"],
                "degree": graph.degree(gene),
                "STRING_mapped": True,
            }
        )
    pd.DataFrame(nodes).sort_values(
        ["degree", "gene"], ascending=[False, True]
    ).to_csv(output_dir / "string_network_nodes.csv", index=False)

    pos = nx.spring_layout(graph, seed=42, weight="score", k=2.0, iterations=400)
    fig, ax = plt.subplots(figsize=(18, 15), facecolor="white")
    ax.set_facecolor("white")

    edge_widths = [
        0.6 + 2.0 * float(attrs["score"]) / 1000
        for _, _, attrs in graph.edges(data=True)
    ]
    edge_collection = nx.draw_networkx_edges(
        graph,
        pos,
        ax=ax,
        width=edge_widths,
        edge_color="#64748B",
        alpha=0.50,
    )
    del edge_collection

    node_colors = [
        "#3182BD" if graph.nodes[gene]["DEG_status"] == "Upregulated" else "#C65353"
        for gene in graph.nodes
    ]
    node_sizes = [300 + 60 * graph.degree(gene) for gene in graph.nodes]
    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="white",
        linewidths=1.2,
    )
    label_pos = {gene: (x, y - 0.055) for gene, (x, y) in pos.items()}
    nx.draw_networkx_labels(
        graph,
        label_pos,
        ax=ax,
        font_size=7.8,
        font_weight="bold",
        font_color="#172B4D",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.5},
    )

    network_label = (
        "physical protein–protein interaction"
        if network_type == "physical"
        else "functional association"
    )
    fig.suptitle(f"STRING {network_label} network", fontsize=17, fontweight="bold", y=0.975)
    ax.set_title(
        "46 GSE125498 MatrisomeAnalyzeR DEGs | Human (taxon 9606) | "
        f"minimum STRING combined score = {MIN_SCORE / 1000:.1f}",
        fontsize=10,
        color="#475569",
        pad=15,
    )
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#3182BD",
               markeredgecolor="white", markersize=11, label="Upregulated (33)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#C65353",
               markeredgecolor="white", markersize=11, label="Downregulated (13)"),
    ]
    ax.legend(handles=legend, loc="upper left", frameon=False, fontsize=10)
    ax.text(
        0.5,
        -0.035,
        f"{graph.number_of_nodes()} queried genes | {graph.number_of_edges()} "
        f"unique {network_type} interactions | node size indicates degree. "
        "STRING confidence scores reflect evidence confidence, not binding strength.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=9,
        color="#475569",
    )
    ax.axis("off")
    fig.tight_layout(rect=(0.02, 0.04, 0.98, 0.94))
    figure_stem = f"GSE125498_46_matrisome_STRING_{network_type}"
    fig.savefig(output_dir / f"{figure_stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(output_dir / f"{figure_stem}.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "source_dataset": "GSE125498",
        "source_table": str(SOURCE_PATH.relative_to(ROOT)),
        "query_date": date.today().isoformat(),
        "string_api": "https://string-db.org/api",
        "species": "Homo sapiens",
        "taxon_id": SPECIES,
        "network_type": network_type,
        "minimum_string_combined_score": MIN_SCORE,
        "minimum_string_combined_score_scaled": MIN_SCORE / 1000,
        "input_gene_count": len(genes),
        "input_direction_counts": {
            "upregulated": int(direction_counts["Upregulated"]),
            "downregulated": int(direction_counts["Downregulated"]),
        },
        "resolved_gene_count": len(resolved_by_query),
        "unresolved_genes": unresolved,
        "STRING_symbol_overrides": {
            gene: {
                "query_symbol": query_symbols[gene],
                "basis": "GSE248762 Ensembl mapping: ENSG00000243477 maps to NAA80",
            }
            for gene in STRING_SYMBOL_OVERRIDES
            if gene in query_symbols
        },
        "network_node_count_including_isolates": graph.number_of_nodes(),
        "network_edge_count": graph.number_of_edges(),
        "connected_node_count": sum(graph.degree(gene) > 0 for gene in graph.nodes),
        "isolated_genes": sorted(gene for gene in graph.nodes if graph.degree(gene) == 0),
        "edge_definition_note": (
            "Direct/physical interactions only."
            if network_type == "physical"
            else "STRING functional associations combine physical and other evidence; "
            "they do not necessarily represent direct physical binding."
        ),
        "confidence_note": (
            "The STRING combined score represents confidence in an interaction, "
            "not interaction strength or direction."
        ),
        "outputs": {
            "png": str((output_dir / f"{figure_stem}.png").relative_to(ROOT)),
            "pdf": str((output_dir / f"{figure_stem}.pdf").relative_to(ROOT)),
            "edges": str(edge_path.relative_to(ROOT)),
            "nodes": str((OUTPUT_DIR / "string_network_nodes.csv").relative_to(ROOT)),
            "mapping": str(resolved_path.relative_to(ROOT)),
        },
    }
    (output_dir / "string_network_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
