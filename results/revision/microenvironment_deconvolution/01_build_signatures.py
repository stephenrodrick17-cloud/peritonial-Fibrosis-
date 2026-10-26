"""Build published human cell-marker signatures and audit panel overlap."""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
EXPR_PATH = ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv"
LIBRARY = "CellMarker_2024"
LIBRARY_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary"

# Each selected term is an exact human signature in Enrichr's CellMarker_2024.
# These are deliberately source terms, not hand-assembled marker lists.
SIGNATURE_TERMS = {
    "T cells": ["T Cell Peripheral Blood Human"],
    "B cells": ["B Cell Peripheral Blood Human"],
    "NK cells": ["Natural Killer Cell Peritoneal Fluid Human"],
    "Monocytes": ["Monocyte Peripheral Blood Human"],
    "M1 macrophages": ["M1 Macrophage Peripheral Blood Human"],
    "M2 macrophages": ["M2 Macrophage Peripheral Blood Human"],
    "Neutrophils": ["Neutrophil Peripheral Blood Human"],
    "Mast cells": ["Mast Cell Lung Human"],
    "Immature dendritic cells": ["Immature Dendritic Cell Stomach Human"],
    "Mature dendritic cells": ["Mature Dendritic Cell Liver Human"],
    "Fibroblasts": ["Fibroblast Epithelium Human"],
    "Myofibroblasts": ["Myofibroblast Epithelium Human"],
    "Mesothelial cells": ["Mesothelial Cell Epithelium Human"],
}
PANELS = {
    "five_gene_nomogram": {"TNFSF15", "ADAM19", "FLT3LG", "NAT6", "ECM1"},
    "seven_gene_consensus": {
        "FLT3LG", "TNFSF15", "LTB", "ADAM19", "SERPINA10", "EBI3", "CXCL14"
    },
}
MIN_DETECTABLE_MARKERS = 3


def download_library() -> dict[str, list[str]]:
    response = requests.get(
        LIBRARY_URL,
        params={"mode": "text", "libraryName": LIBRARY},
        timeout=90,
    )
    response.raise_for_status()
    if not response.text.strip():
        raise RuntimeError("Enrichr returned an empty CellMarker_2024 library.")
    library: dict[str, list[str]] = {}
    for line in response.text.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        term = fields[0].strip()
        genes = [gene.strip().upper() for gene in fields[2:] if gene.strip()]
        if term and genes:
            library[term] = genes
    (OUT / "cellmarker_2024_library_snapshot.txt").write_text(
        response.text, encoding="utf-8"
    )
    if not library:
        raise RuntimeError("Could not parse any signatures from CellMarker_2024.")
    return library


def main() -> None:
    if not EXPR_PATH.is_file():
        raise FileNotFoundError(f"Expression matrix not found: {EXPR_PATH}")
    expression = pd.read_csv(EXPR_PATH)
    if "Gene" not in expression:
        raise ValueError("Expression matrix has no Gene column.")
    expressed = set(expression["Gene"].astype(str).str.upper())
    if len(expression) != 11741 or expression["Gene"].duplicated().any():
        raise ValueError("Expected the existing 11,741-row unique-symbol expression matrix.")

    library = download_library()
    rows: list[dict[str, str | int]] = []
    summaries: list[dict[str, str | int]] = []
    all_overlaps: list[dict[str, str]] = []
    manifest: list[dict[str, str | int]] = []
    missing_terms: list[str] = []

    for population, terms in SIGNATURE_TERMS.items():
        source_genes: set[str] = set()
        for term in terms:
            if term not in library:
                missing_terms.append(term)
                continue
            source_genes.update(library[term])
        if not source_genes:
            continue
        present = sorted(source_genes & expressed)
        absent_count = len(source_genes - expressed)
        manifest.append(
            {
                "population": population,
                "database": "Enrichr CellMarker_2024",
                "source_signature_terms": "; ".join(terms),
                "marker_count_in_source": len(source_genes),
                "marker_count_in_expression_matrix": len(present),
                "source_markers_not_detected": absent_count,
                "run_date": date.today().isoformat(),
            }
        )
        summaries.append(
            {
                "population": population,
                "source_markers": len(source_genes),
                "detectable_markers": len(present),
                "below_3_detectable_markers": len(present) < MIN_DETECTABLE_MARKERS,
            }
        )
        for gene in present:
            rows.append(
                {
                    "population": population,
                    "gene": gene,
                    "source": f"{LIBRARY}: {'; '.join(terms)}",
                }
            )
        for panel_name, panel_genes in PANELS.items():
            overlap = sorted(set(present) & panel_genes)
            for gene in overlap:
                all_overlaps.append(
                    {"population": population, "gene": gene, "panel": panel_name}
                )

    if missing_terms:
        raise RuntimeError(f"Signature terms were missing from {LIBRARY}: {missing_terms}")
    signature_table = pd.DataFrame(rows)
    counts = pd.DataFrame(summaries)
    overlaps = pd.DataFrame(all_overlaps, columns=["population", "gene", "panel"])
    signature_table.to_csv(OUT / "cell_marker_signatures_detectable.csv", index=False)
    counts.to_csv(OUT / "signature_marker_counts.csv", index=False)
    overlaps.to_csv(OUT / "signature_panel_overlaps.csv", index=False)
    pd.DataFrame(manifest).to_csv(OUT / "signature_source_manifest.csv", index=False)
    metadata = {
        "run_date": date.today().isoformat(),
        "source_database": LIBRARY,
        "source_url": LIBRARY_URL,
        "species": "Homo sapiens",
        "expression_matrix": "results/tables/GSE125498_full_expression_matrix.csv",
        "expressed_gene_symbols": len(expressed),
        "minimum_detectable_markers": MIN_DETECTABLE_MARKERS,
        "score_method": "mean row-wise z-score across samples for each signature",
        "signature_terms": SIGNATURE_TERMS,
        "panels_checked": {key: sorted(value) for key, value in PANELS.items()},
        "overlap_genes_by_panel": {
            panel: sorted(
                set(overlaps.loc[overlaps["panel"] == panel, "gene"])
            ) if len(overlaps) else []
            for panel in PANELS
        },
        "overlap_population_gene_pairs_by_panel": {
            panel: overlaps.loc[overlaps["panel"] == panel].to_dict("records")
            for panel in PANELS
        },
        "overlapping_genes_excluded": False,
    }
    (OUT / "signature_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(counts.to_string(index=False))
    print("\nPanel overlaps:")
    print(overlaps.to_string(index=False) if len(overlaps) else "None")


if __name__ == "__main__":
    main()
