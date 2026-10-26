"""Create a STRING physical-interaction network for the historical 11-gene hub list."""

from __future__ import annotations

import csv
import argparse
import io
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_BASE = ROOT / "results" / "revision"
GENES = [
    "FLT3LG",
    "TNFSF15",
    "ADAM19",
    "SERPINA10",
    "LTB",
    "CLEC4F",
    "SEMA3E",
    "SEMA4F",
    "EBI3",
    "LTA",
    "CXCL14",
]
SPECIES = 9606
API_URL = "https://string-db.org/api/tsv/network"


def fetch_interactions(network_type: str, min_score: int) -> list[dict[str, str]]:
    form = urlencode(
        {
            "identifiers": "\r\n".join(GENES),
            "species": str(SPECIES),
            "required_score": str(min_score),
            "network_type": network_type,
            "caller_identity": "peritoneal_fibrosis_analysis",
        }
    ).encode("ascii")
    request = Request(API_URL, data=form)
    with urlopen(request, timeout=60) as response:
        payload = response.read().decode("utf-8")

    reader = csv.DictReader(io.StringIO(payload), delimiter="\t")
    required_columns = {
        "preferredName_A",
        "preferredName_B",
        "score",
        "stringId_A",
        "stringId_B",
    }
    if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
        raise ValueError("STRING returned an unexpected network response.")
    rows = list(reader)
    for row in rows:
        if row["preferredName_A"] not in GENES or row["preferredName_B"] not in GENES:
            raise ValueError("STRING returned an interaction outside the requested gene list.")
    return rows


def save_tables(output_dir: Path, rows: list[dict[str, str]]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = output_dir / "string_interactions.tsv"
    raw_fields = [
        "stringId_A",
        "stringId_B",
        "preferredName_A",
        "preferredName_B",
        "ncbiTaxonId",
        "score",
        "nscore",
        "fscore",
        "pscore",
        "ascore",
        "escore",
        "dscore",
        "tscore",
    ]
    with raw_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=raw_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    degree = {gene: 0 for gene in GENES}
    for row in rows:
        degree[row["preferredName_A"]] += 1
        degree[row["preferredName_B"]] += 1

    node_path = output_dir / "network_nodes.csv"
    with node_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=["gene", "degree"])
        writer.writeheader()
        writer.writerows(
            {"gene": gene, "degree": degree[gene]} for gene in GENES
        )


def draw_network(
    output_dir: Path, rows: list[dict[str, str]], network_type: str, min_score: int
) -> None:
    linked_genes = {row["preferredName_A"] for row in rows}
    linked_genes.update(row["preferredName_B"] for row in rows)

    positions = {
        "LTA": (-2.05, 0.45),
        "LTB": (-1.15, 0.45),
    }
    remaining_genes = [gene for gene in GENES if gene not in positions]
    center_x, center_y, radius = 1.45, 0.0, 2.05
    for index, gene in enumerate(remaining_genes):
        angle = math.pi / 2 + (2 * math.pi * index / len(remaining_genes))
        positions[gene] = (
            center_x + radius * math.cos(angle),
            center_y + radius * math.sin(angle),
        )

    figure, axis = plt.subplots(figsize=(12, 8))
    network_label = (
        "physical-interaction"
        if network_type == "physical"
        else "functional-association"
    )
    axis.set_title(
        f"STRING {network_label} network: historical 11-gene hub list",
        fontsize=15,
        pad=18,
    )
    axis.text(
        0.5,
        1.01,
        f"Human (taxon {SPECIES}) | minimum STRING score: {min_score / 1000:.3f}",
        transform=axis.transAxes,
        ha="center",
        fontsize=10,
        color="#444444",
    )

    for row in rows:
        gene_a, gene_b = row["preferredName_A"], row["preferredName_B"]
        x_values = [positions[gene_a][0], positions[gene_b][0]]
        y_values = [positions[gene_a][1], positions[gene_b][1]]
        axis.plot(x_values, y_values, color="#536878", linewidth=2.5, zorder=1)
        midpoint = (
            (x_values[0] + x_values[1]) / 2,
            (y_values[0] + y_values[1]) / 2,
        )
        axis.text(
            midpoint[0],
            midpoint[1] + 0.18,
            f'{float(row["score"]):.3f}',
            ha="center",
            fontsize=9,
            color="#333333",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
            zorder=3,
        )

    for gene in GENES:
        x, y = positions[gene]
        color = "#3B82A0" if gene in linked_genes else "#B8C2CC"
        axis.scatter(
            [x],
            [y],
            s=2300,
            color=color,
            edgecolor="white",
            linewidth=1.8,
            zorder=2,
        )
        axis.text(
            x,
            y,
            gene,
            ha="center",
            va="center",
            color="white" if gene in linked_genes else "#344054",
            fontsize=8,
            weight="bold",
            zorder=4,
        )

    axis.text(
        0.5,
        -0.08,
        f"{len(GENES)} queried genes | {len(rows)} edge(s) at threshold",
        transform=axis.transAxes,
        ha="center",
        fontsize=10,
        color="#444444",
    )
    axis.text(
        0.5,
        -0.13,
        "STRING score reflects evidence confidence, not binding strength; isolated nodes had no qualifying edge.",
        transform=axis.transAxes,
        ha="center",
        fontsize=8.5,
        color="#555555",
    )
    axis.set_xlim(-2.8, 3.9)
    axis.set_ylim(-2.8, 2.8)
    axis.set_aspect("equal")
    axis.axis("off")
    figure.tight_layout()
    output_name = f"legacy_11_hub_{network_type}_network"
    figure.savefig(output_dir / f"{output_name}.png", dpi=300, bbox_inches="tight")
    figure.savefig(output_dir / f"{output_name}.pdf", bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--network-type",
        choices=("physical", "functional"),
        default="physical",
        help="STRING edge definition (default: physical)",
    )
    parser.add_argument(
        "--min-score",
        type=int,
        choices=(150, 400, 500, 700, 900),
        default=400,
        help="Minimum STRING confidence score on the 0-1000 scale (default: 400)",
    )
    args = parser.parse_args()
    output_name = (
        "legacy_11_hub_ppi"
        if args.network_type == "physical" and args.min_score == 400
        else f"legacy_11_hub_{args.network_type}_{args.min_score}"
    )
    output_dir = OUTPUT_BASE / output_name
    rows = fetch_interactions(args.network_type, args.min_score)
    save_tables(output_dir, rows)
    draw_network(output_dir, rows, args.network_type, args.min_score)
    print(
        f"Saved network for {len(GENES)} genes with {len(rows)} "
        f"{args.network_type} edge(s) to {output_dir}"
    )


if __name__ == "__main__":
    main()