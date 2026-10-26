"""Score sample-level microenvironment signatures with marker z-score means."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
EXPR_PATH = ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv"
META_PATH = ROOT / "results" / "tables" / "GSE125498_sample_metadata.csv"
MIN_MARKERS = 3
FULL = "full_signature"
EXCLUDED = "overlap_excluded_7gene"
CONSENSUS_PANEL = {"FLT3LG", "TNFSF15", "LTB", "ADAM19", "SERPINA10", "EBI3", "CXCL14"}


def main() -> None:
    signatures = pd.read_csv(OUT / "cell_marker_signatures_detectable.csv")
    metadata = json.loads((OUT / "signature_metadata.json").read_text(encoding="utf-8"))
    expression = pd.read_csv(EXPR_PATH).set_index("Gene")
    samples = pd.read_csv(META_PATH)
    if expression.index.duplicated().any():
        raise ValueError("Expected unique gene symbols in the normalized matrix.")
    if samples["sample_id"].duplicated().any() or not set(samples["sample_id"]).issubset(expression.columns):
        raise ValueError("Sample metadata does not match the expression-matrix columns.")
    if len(samples) != 33 or (samples["group"] == "LPD").sum() != 13 or (samples["group"] == "SPD").sum() != 20:
        raise ValueError("Expected 33 samples: 13 LPD and 20 SPD.")
    expression = expression.loc[:, samples["sample_id"]].astype(float)
    z = expression.sub(expression.mean(axis=1), axis=0).div(
        expression.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )

    overlap_genes = set(metadata["overlap_genes_by_panel"]["seven_gene_consensus"])
    output_rows: list[dict[str, object]] = []
    detectable_counts = []
    for population, frame in signatures.groupby("population", sort=False):
        markers = sorted(set(frame["gene"]))
        markers_without_overlap = [gene for gene in markers if gene not in overlap_genes]
        for variant, genes in ((FULL, markers), (EXCLUDED, markers_without_overlap)):
            if len(genes) < MIN_MARKERS:
                detectable_counts.append(
                    {
                        "population": population,
                        "signature_variant": variant,
                        "detectable_markers": len(genes),
                        "included": False,
                        "drop_reason": f"fewer than {MIN_MARKERS} detectable markers",
                    }
                )
                continue
            for sample_id in samples["sample_id"]:
                output_rows.append(
                    {
                        "sample_id": sample_id,
                        "population": population,
                        "signature_variant": variant,
                        "score_mean_z": float(z.loc[genes, sample_id].mean()),
                    }
                )
            detectable_counts.append(
                {
                    "population": population,
                    "signature_variant": variant,
                    "detectable_markers": len(genes),
                    "included": True,
                    "drop_reason": "",
                }
            )

    long = pd.DataFrame(output_rows)
    long = long.merge(samples[["sample_id", "group"]], on="sample_id", validate="many_to_one")
    long.to_csv(OUT / "sample_population_scores_long.csv", index=False)
    count_table = pd.DataFrame(detectable_counts).merge(
        pd.read_csv(OUT / "signature_marker_counts.csv"),
        on="population",
        how="left",
        validate="many_to_one",
    )
    count_table.to_csv(OUT / "signature_marker_counts_by_variant.csv", index=False)

    output_signatures = []
    for _, row in signatures.iterrows():
        output_signatures.append(
            {
                "population": row["population"],
                "gene": row["gene"],
                "source": row["source"],
                "signature_variant": FULL,
                "included": True,
            }
        )
        output_signatures.append(
            {
                "population": row["population"],
                "gene": row["gene"],
                "source": row["source"],
                "signature_variant": EXCLUDED,
                "included": row["gene"] not in overlap_genes,
            }
        )
    pd.DataFrame(output_signatures).to_csv(
        OUT / "final_signature_table_both_variants.csv", index=False
    )

    for variant in (FULL, EXCLUDED):
        variant_scores = long[long["signature_variant"] == variant]
        matrix = variant_scores.pivot(
            index="sample_id", columns="population", values="score_mean_z"
        ).reindex(samples["sample_id"])
        matrix.insert(0, "group", samples.set_index("sample_id").loc[matrix.index, "group"])
        matrix.to_csv(OUT / f"sample_population_score_matrix_{variant}.csv", index_label="sample_id")
    print(
        f"Scored {long['population'].nunique()} populations for {len(samples)} samples "
        f"using both full and overlap-excluded signatures."
    )


if __name__ == "__main__":
    main()
