"""Query GeneMANIA and export a seven-seed functional association network."""

from __future__ import annotations

import json
from collections import defaultdict
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import requests


OUT = Path(__file__).resolve().parent
GENES = ["FLT3LG", "TNFSF15", "LTB", "ADAM19", "SERPINA10", "EBI3", "CXCL14"]
GENEMANIA_URL = "https://genemania.org/json/search_results"
EXPECTED_TYPES = [
    "Co-expression",
    "Physical interaction",
    "Co-localization",
    "Pathway",
    "Genetic interaction",
    "Shared protein domains",
    "Predicted",
]


def main() -> None:
    payload = {
        "organism": 4,
        "genes": "\n".join(GENES),
        "weighting": "AUTOMATIC_SELECT",
        "resultsSize": 20,
        "attributeResultsSize": 10,
    }
    response = requests.post(GENEMANIA_URL, json=payload, timeout=120)
    response.raise_for_status()
    result = response.json()
    if result.get("error"):
        raise RuntimeError(f"GeneMANIA query failed: {result.get('error')}")

    (OUT / "genemania_raw_response.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    groups = result["resultNetworkGroups"]
    selected_types = {
        group["networkGroup"]["name"]: float(group["weight"]) for group in groups
    }
    pd.DataFrame(
        [
            {
                "interaction_type": name,
                "included_in_result": name in selected_types,
                "weight": selected_types.get(name, 0.0),
            }
            for name in EXPECTED_TYPES
        ]
    ).to_csv(OUT / "genemania_interaction_type_weights.csv", index=False)

    network_rows = []
    edge_rows = []
    for group in groups:
        group_name = group["networkGroup"]["name"]
        group_weight = float(group["weight"])
        for network in group["resultNetworks"]:
            network_name = network["network"]["name"]
            network_weight = float(network["weight"])
            network_rows.append(
                {
                    "interaction_type": group_name,
                    "source_network": network_name,
                    "interaction_type_weight": group_weight,
                    "source_network_weight": network_weight,
                    "combined_weight_contribution": group_weight * network_weight,
                }
            )
            for item in network["resultInteractions"]:
                source = item["fromGene"]
                target = item["toGene"]
                source_gene = source["gene"]["symbol"]
                target_gene = target["gene"]["symbol"]
                edge_rows.append(
                    {
                        "gene_a": source_gene,
                        "gene_b": target_gene,
                        "gene_a_is_query": bool(source["queryGene"]),
                        "gene_b_is_query": bool(target["queryGene"]),
                        "interaction_type": group_name,
                        "source_network": network_name,
                        "interaction_weight": float(item["interaction"]["weight"]),
                        "interaction_type_weight": group_weight,
                        "source_network_weight": network_weight,
                        "combined_weight_contribution": group_weight
                        * network_weight
                        * float(item["interaction"]["weight"]),
                    }
                )
    edges = pd.DataFrame(edge_rows)
    if edges.empty:
        raise RuntimeError("GeneMANIA returned no network edges.")
    edges.to_csv(OUT / "genemania_network_edges.csv", index=False)
    pd.DataFrame(network_rows).to_csv(
        OUT / "genemania_source_network_weights.csv", index=False
    )

    gene_rows = []
    for item in result["resultGenes"]:
        gene_rows.append(
            {
                "gene": item["gene"]["symbol"],
                "node_role": "query hub" if item["queryGene"] else "predicted interactor",
                "query_gene": bool(item["queryGene"]),
                "genemania_gene_score": float(item["score"]),
            }
        )
    nodes = pd.DataFrame(gene_rows)
    recognized = set(nodes.loc[nodes["query_gene"], "gene"])
    missing = sorted(set(GENES) - recognized)
    pd.DataFrame(
        [{"query_gene": gene, "recognized": gene in recognized} for gene in GENES]
    ).to_csv(OUT / "genemania_query_gene_recognition.csv", index=False)
    if missing:
        raise RuntimeError(f"GeneMANIA did not recognize query genes: {missing}")
    if len(nodes.loc[~nodes["query_gene"]]) != 20:
        raise RuntimeError(
            "GeneMANIA did not return exactly the requested 20 resultant genes."
        )
    nodes.to_csv(OUT / "genemania_network_nodes.csv", index=False)

    query_set = set(GENES)
    support_rows = []
    for row in edge_rows:
        if row["gene_a_is_query"] and not row["gene_b_is_query"]:
            hub, interactor = row["gene_a"], row["gene_b"]
        elif row["gene_b_is_query"] and not row["gene_a_is_query"]:
            hub, interactor = row["gene_b"], row["gene_a"]
        else:
            continue
        if hub not in query_set:
            continue
        support_rows.append(
            {
                "interactor_gene": interactor,
                "hub_gene": hub,
                "interaction_type": row["interaction_type"],
                "source_network": row["source_network"],
                "interaction_weight": row["interaction_weight"],
                "interaction_type_weight": row["interaction_type_weight"],
                "source_network_weight": row["source_network_weight"],
                "combined_weight_contribution": row[
                    "combined_weight_contribution"
                ],
            }
        )
    support = pd.DataFrame(support_rows)
    support.to_csv(OUT / "genemania_interactor_gene_weights.csv", index=False)
    if support.empty:
        raise RuntimeError("No direct query-gene to predicted-interactor edges returned.")

    # Plot only edges incident to the query genes to keep the seed/interactor
    # distinction readable; the full GeneMANIA edge list is exported above.
    graph = nx.Graph()
    graph.add_nodes_from(GENES, role="query hub")
    for _, row in support.iterrows():
        graph.add_edge(
            row["hub_gene"],
            row["interactor_gene"],
            weight=max(float(row["interaction_weight"]), 0.001),
            interaction_type=row["interaction_type"],
        )
    interactor_nodes = sorted(set(graph) - query_set)
    hubs_pos = nx.circular_layout(GENES, scale=1.0)
    interactor_pos = nx.spring_layout(
        graph.subgraph(interactor_nodes), seed=20261013, k=1.0
    ) if interactor_nodes else {}
    for index, gene in enumerate(interactor_nodes):
        x, y = interactor_pos[gene]
        interactor_pos[gene] = (x + 2.2, y)
    positions = {**hubs_pos, **interactor_pos}

    fig, ax = plt.subplots(figsize=(15, 10))
    widths = [
        0.6 + 2.0 * graph[u][v]["weight"]
        for u, v in graph.edges()
    ]
    nx.draw_networkx_edges(
        graph, positions, ax=ax, width=widths, edge_color="#8A94A6", alpha=0.55
    )
    nx.draw_networkx_nodes(
        graph, positions, nodelist=GENES, node_color="#F2A541",
        node_size=1500, edgecolors="#7A3E00", linewidths=2, ax=ax,
        label="Query hub gene",
    )
    if interactor_nodes:
        nx.draw_networkx_nodes(
            graph, positions, nodelist=interactor_nodes, node_color="#8EC9E8",
            node_size=900, edgecolors="#1E5A78", linewidths=1.3, ax=ax,
            label="GeneMANIA predicted interactor",
        )
    nx.draw_networkx_labels(graph, positions, font_size=8, ax=ax)
    ax.set_title(
        "GeneMANIA human functional-association network\n"
        "Seven query genes; up to 20 additional genes; automatic network weighting"
    )
    ax.legend(loc="upper left", frameon=False)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "genemania_network.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "genemania_network.pdf", bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "query_date": date.today().isoformat(),
        "database": "GeneMANIA",
        "organism": "Homo sapiens",
        "organism_id": 4,
        "query_genes": GENES,
        "recognized_query_genes": sorted(recognized),
        "unrecognized_query_genes": missing,
        "max_result_genes": payload["resultsSize"],
        "weighting": payload["weighting"],
        "result_gene_count_including_queries": len(nodes),
        "result_interactor_count": int((~nodes["query_gene"]).sum()),
        "interaction_types_and_weights": selected_types,
        "query_url": GENEMANIA_URL,
    }
    (OUT / "genemania_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
