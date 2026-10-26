"""Query DSigDB through Enrichr for candidate drug-associated signatures."""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import requests


OUT = Path(__file__).resolve().parent
GENES = ["FLT3LG", "TNFSF15", "LTB", "ADAM19", "SERPINA10", "EBI3", "CXCL14"]
LIBRARY = "DSigDB"
FDR_CUTOFF = 0.05
BASE = "https://maayanlab.cloud/Enrichr"


def main() -> None:
    add_response = requests.post(
        f"{BASE}/addList",
        files={
            "list": (None, "\n".join(GENES)),
            "description": (None, "GSE125498 seven majority-vote diagnostic hub genes"),
        },
        timeout=60,
    )
    add_response.raise_for_status()
    user_list_id = add_response.json()["userListId"]
    enrich_response = requests.get(
        f"{BASE}/enrich",
        params={"userListId": user_list_id, "backgroundType": LIBRARY},
        timeout=120,
    )
    enrich_response.raise_for_status()
    enrich_json = enrich_response.json()
    if LIBRARY not in enrich_json:
        raise RuntimeError(f"Enrichr response did not include {LIBRARY}.")

    statistics_response = requests.get(f"{BASE}/datasetStatistics", timeout=60)
    statistics_response.raise_for_status()
    stats_payload = statistics_response.json()
    stats_rows = stats_payload.get("statistics", [])
    library_stats = next(
        (row for row in stats_rows if row.get("libraryName") == LIBRARY), None
    )
    if library_stats is None:
        raise RuntimeError(f"Enrichr did not report library metadata for {LIBRARY}.")

    rows = []
    for values in enrich_json[LIBRARY]:
        rank, term, p_value, z_score, combined_score, genes, adjusted_p, *_ = values
        overlap_genes = list(genes)
        rows.append(
            {
                "rank": int(rank),
                "dsigdb_term": term,
                "overlapping_genes": ";".join(overlap_genes),
                "overlap_size": len(overlap_genes),
                "p_value": float(p_value),
                "z_score": float(z_score),
                "combined_score": float(combined_score),
                "adjusted_p_value": float(adjusted_p),
                "significant_fdr_lt_0_05": float(adjusted_p) < FDR_CUTOFF,
            }
        )
    results = pd.DataFrame(rows).sort_values(
        ["adjusted_p_value", "p_value", "rank"]
    )
    results.to_csv(OUT / "dsigdb_enrichment_full.csv", index=False)
    significant = results[results["adjusted_p_value"] < FDR_CUTOFF].copy()
    significant.to_csv(OUT / "dsigdb_candidate_drugs_fdr_lt_0_05.csv", index=False)

    library_response = requests.get(
        f"{BASE}/geneSetLibrary",
        params={"mode": "text", "libraryName": LIBRARY},
        timeout=120,
    )
    library_response.raise_for_status()
    annotated_genes: set[str] = set()
    for line in library_response.text.splitlines():
        fields = line.rstrip("\r").split("\t")
        annotated_genes.update(value for value in fields[1:] if value)
    coverage = pd.DataFrame(
        [
            {
                "hub_gene": gene,
                "has_dsigdb_drug_signature_annotation": gene in annotated_genes,
            }
            for gene in GENES
        ]
    )
    coverage.to_csv(OUT / "dsigdb_gene_annotation_coverage.csv", index=False)
    coverage.loc[
        ~coverage["has_dsigdb_drug_signature_annotation"], ["hub_gene"]
    ].to_csv(OUT / "dsigdb_genes_without_drug_annotation.csv", index=False)

    plot_rows = (significant if not significant.empty else results).head(15).copy()
    plot_rows = plot_rows.sort_values("adjusted_p_value", ascending=False)
    labels = [
        f"{term} [{genes}]"
        for term, genes in zip(
            plot_rows["dsigdb_term"], plot_rows["overlapping_genes"]
        )
    ]
    fig_height = max(5.0, 0.38 * len(plot_rows) + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_height))
    colors = [
        "#C44E52" if value < FDR_CUTOFF else "#9AA0A6"
        for value in plot_rows["adjusted_p_value"]
    ]
    ax.barh(labels, -plot_rows["adjusted_p_value"].map(
        lambda value: __import__("math").log10(max(value, 1e-300))
    ), color=colors)
    ax.axvline(-__import__("math").log10(FDR_CUTOFF), color="#333333", ls="--", lw=1)
    ax.set_xlabel("-log10(Benjamini-Hochberg adjusted P)")
    ax.set_title(
        "Candidate DSigDB signatures associated with the hub-gene set"
        + ("" if not significant.empty else "\nNo term passes FDR < 0.05")
    )
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "dsigdb_top_candidate_drugs.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "dsigdb_top_candidate_drugs.pdf", bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "query_date": date.today().isoformat(),
        "database": "DSigDB",
        "library": LIBRARY,
        "library_version": "DSigDB v1.0 as listed by Enrichr",
        "library_statistics": library_stats,
        "query_genes": GENES,
        "enrichr_user_list_id": user_list_id,
        "fdr_cutoff": FDR_CUTOFF,
        "number_of_terms_returned": len(results),
        "number_significant_fdr_lt_0_05": len(significant),
        "hub_genes_without_drug_annotation": coverage.loc[
            ~coverage["has_dsigdb_drug_signature_annotation"], "hub_gene"
        ].tolist(),
    }
    (OUT / "dsigdb_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
