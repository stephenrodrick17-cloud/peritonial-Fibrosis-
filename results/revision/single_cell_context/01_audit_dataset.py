"""Audit the GEO sample metadata and reproduce the per-cell QC counts."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "processed" / "GSE248762_hubblind_allcells_qc.h5ad"
GEO_DIR = ROOT / "results" / "revision" / "effluent_validation" / "GSE248762"
METADATA = GEO_DIR / "sample_metadata.csv"
PREVIOUS_QC = GEO_DIR / "analysis_20261009" / "sample_qc.csv"
UMAP = GEO_DIR / "analysis_20261009" / "umap_coordinates.csv"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(METADATA, dtype=str, keep_default_na=False)
    previous = pd.read_csv(PREVIOUS_QC)
    coords = pd.read_csv(UMAP, usecols=["sample_id", "cell_id"])
    data = ad.read_h5ad(DATA, backed="r")
    obs = data.obs.copy()
    obs.index = obs.index.astype(str)
    data.file.close()

    if metadata["sample_id"].duplicated().any():
        raise ValueError("GEO sample metadata contains duplicate sample IDs.")
    if set(metadata["sample_id"]) != set(obs["gsm"].astype(str).unique()):
        raise ValueError("The GEO metadata sample set does not match the count object.")

    expected_original_groups = {"SV": 6, "LV_NOT_UF": 6, "LV_UF": 4}
    observed = metadata["group_from_metadata"].value_counts().to_dict()
    if observed != expected_original_groups or len(metadata) != 16:
        raise ValueError(f"Unexpected GSE248762 sample groups: {observed}")
    if sum(observed.values()) != len(metadata):
        raise ValueError("Original group counts do not sum to the sample total.")

    obs["pass_context_qc"] = (
        (obs["n_genes"].astype(float) >= 200)
        & (obs["pct_counts_mt"].astype(float) <= 20.0)
    )
    before = obs.groupby("gsm", observed=True).size().rename("cells_before_qc")
    after = (
        obs.loc[obs["pass_context_qc"]]
        .groupby("gsm", observed=True)
        .size()
        .rename("cells_after_qc_recomputed")
    )
    audit = (
        metadata[["sample_id", "title", "source_name", "group_from_metadata"]]
        .merge(previous[["sample_id", "cells_after_qc"]], on="sample_id", how="left", validate="one_to_one")
        .merge(before, left_on="sample_id", right_index=True, validate="one_to_one")
        .merge(after, left_on="sample_id", right_index=True, validate="one_to_one")
    )
    audit["primary_contrast_group"] = audit["group_from_metadata"].map(
        lambda group: "SV" if group == "SV" else "LV"
    )
    audit["difference_vs_saved_qc_cells"] = (
        audit["cells_after_qc_recomputed"] - audit["cells_after_qc"]
    )
    audit.to_csv(OUT / "sample_group_qc_audit.csv", index=False)

    coord_keys = {
        (str(row.sample_id), str(row.cell_id).split(":")[0])
        for row in coords.itertuples(index=False)
    }
    retained = obs.loc[obs["pass_context_qc"], ["gsm"]].copy()
    retained["barcode"] = retained.index.str.split("_", n=1).str[1]
    filter_keys = set(
        zip(retained["gsm"].astype(str), retained["barcode"].astype(str))
    )
    coord_only = coord_keys - filter_keys
    if len(filter_keys - coord_keys) != 0:
        raise ValueError("Filtered count-object cells are missing UMAP coordinates.")
    if len(coord_only) != 3:
        raise ValueError(
            f"Expected to audit three UMAP-only cells above the QC threshold; found {len(coord_only)}."
        )
    non_qc_umap = []
    for gsm, barcode in sorted(coord_only):
        cell = obs.loc[f"{gsm}_{barcode}"]
        non_qc_umap.append(
            {
                "sample_id": gsm,
                "barcode": barcode,
                "n_genes": int(cell["n_genes"]),
                "pct_counts_mt": float(cell["pct_counts_mt"]),
                "reason": "mitochondrial fraction > 20%",
            }
        )
    pd.DataFrame(non_qc_umap).to_csv(OUT / "umap_cells_excluded_by_recomputed_qc.csv", index=False)

    original_group_counts = metadata["group_from_metadata"].value_counts().to_dict()
    primary_counts = (
        metadata["group_from_metadata"]
        .map(lambda group: "SV" if group == "SV" else "LV")
        .value_counts()
        .to_dict()
    )
    local_download_timestamp = pd.Timestamp(
        pd.Timestamp.fromtimestamp(GEO_DIR.joinpath("GSE248762_RAW.tar").stat().st_mtime)
    ).isoformat()
    report = {
        "accession": "GSE248762",
        "data_type": "10x Genomics single-cell RNA-seq",
        "tissue": "overnight peritoneal dialysis effluent",
        "geo_sample_records": len(metadata),
        "original_group_counts": original_group_counts,
        "original_group_counts_sum_matches_total": sum(original_group_counts.values()) == len(metadata),
        "primary_contrast": "LV (LV_UF + LV_NOT_UF) versus SV",
        "primary_contrast_sample_counts": primary_counts,
        "primary_contrast_counts_sum_matches_total": sum(primary_counts.values()) == len(metadata),
        "cells_before_qc": int(len(obs)),
        "cells_after_recomputed_qc": int(obs["pass_context_qc"].sum()),
        "qc": {"minimum_detected_genes": 200, "maximum_mitochondrial_percent_inclusive": 20.0},
        "previous_saved_qc_cells_after": int(previous["cells_after_qc"].sum()),
        "recomputed_minus_previous_qc_cells": int(
            obs["pass_context_qc"].sum() - previous["cells_after_qc"].sum()
        ),
        "umap_only_cells_excluded": len(coord_only),
        "umap_only_cells_exclusion_reason": "pct_counts_mt > 20% in the current count object",
        "cached_geo_archive": str(GEO_DIR.joinpath("GSE248762_RAW.tar").relative_to(ROOT)),
        "cached_geo_archive_local_timestamp": local_download_timestamp,
        "original_geo_download_date_independently_logged": False,
    }
    (OUT / "dataset_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
