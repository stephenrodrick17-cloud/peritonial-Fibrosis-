"""Assign hub-blind broad lineages using reference marker-gene scores."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "processed" / "GSE248762_hubblind_allcells_qc.h5ad"
GEO_DIR = ROOT / "results" / "revision" / "effluent_validation" / "GSE248762"
UMAP_PATH = GEO_DIR / "analysis_20261009" / "umap_coordinates.csv"
MIN_GENES = 200
MAX_MT_PERCENT = 20.0

MARKERS = {
    "T/NK": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD74"],
    "Monocyte/macrophage": [
        "LST1", "TYROBP", "FCER1G", "CD14", "LILRB1", "C1QA", "C1QB"
    ],
    "Neutrophil": ["CSF3R", "FCGR3B", "S100A8", "S100A9", "CXCR2"],
    "Fibroblast/stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"],
    "Mesothelial": ["KRT8", "KRT18", "KRT19", "MSLN", "UPK3B", "KRT7"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "CDH5", "EMCN"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"],
}
SUBTYPE_TO_LINEAGE = {
    "T/NK": "lymphoid",
    "B": "lymphoid",
    "Monocyte/macrophage": "myeloid",
    "Neutrophil": "myeloid",
    "Fibroblast/stromal": "fibroblast/stromal",
    "Mesothelial": "mesothelial",
    "Endothelial": "endothelial",
    "Epithelial": "epithelial",
    "Other/unresolved": "other/unresolved",
}
MARKER_SOURCE = (
    "Project marker-score annotation panel; broad lineage marker reference: "
    "CellMarker 2.0 (Hu et al., Nucleic Acids Research 51(D1), D870-D876, "
    "2023; doi:10.1093/nar/gkac947; PMID:36300619)."
)


def main() -> None:
    audit = pd.read_json(OUT / "dataset_audit.json", typ="series")
    if int(audit["cells_after_recomputed_qc"]) <= 0:
        raise ValueError("Dataset audit reports no QC-retained cells.")

    coordinates = pd.read_csv(UMAP_PATH)
    coordinates["barcode"] = coordinates["cell_id"].astype(str).str.split(":").str[0]
    data = ad.read_h5ad(DATA, backed="r")
    obs = data.obs.copy()
    obs.index = obs.index.astype(str)
    obs["object_cell_id"] = obs.index
    obs["barcode"] = obs.index.str.split("_", n=1).str[1]
    obs["pass_context_qc"] = (
        (obs["n_genes"].astype(float) >= MIN_GENES)
        & (obs["pct_counts_mt"].astype(float) <= MAX_MT_PERCENT)
    )
    obs = obs.loc[obs["pass_context_qc"]].copy()
    coords = coordinates.merge(
        obs,
        left_on=["sample_id", "barcode"],
        right_on=["gsm", "barcode"],
        how="inner",
        validate="one_to_one",
        suffixes=("", "_object"),
    )
    if len(coords) != int(audit["cells_after_recomputed_qc"]):
        data.file.close()
        raise ValueError("UMAP coordinates do not cover all QC-retained cells.")
    if coords["object_cell_id"].duplicated().any():
        data.file.close()
        raise ValueError("Cell-to-coordinate mapping is not one-to-one.")

    symbol_names = pd.Index(data.var_names.astype(str)).str.upper()
    marker_columns: dict[str, list[int]] = {
        marker: np.flatnonzero(symbol_names == marker).tolist()
        for marker in sorted({gene for genes in MARKERS.values() for gene in genes})
    }
    available = {
        cell_type: [gene for gene in genes if marker_columns[gene]]
        for cell_type, genes in MARKERS.items()
    }
    missing = {
        cell_type: [gene for gene in genes if not marker_columns[gene]]
        for cell_type, genes in MARKERS.items()
    }

    subtype_assignments = pd.Series(index=coords.index, dtype="object")
    for sample_id, subset in coords.groupby("sample_id", sort=False):
        object_rows = data.obs_names.get_indexer(subset["object_cell_id"])
        if np.any(object_rows < 0):
            data.file.close()
            raise ValueError(f"{sample_id}: annotation references missing count rows.")
        counts = data.X[object_rows, :]
        if not sparse.issparse(counts):
            counts = sparse.csr_matrix(counts)
        else:
            counts = counts.tocsr()
        library_size = np.asarray(counts.sum(axis=1)).ravel().astype(float)
        if np.any(library_size <= 0):
            data.file.close()
            raise ValueError(f"{sample_id}: zero-library cells passed QC.")

        type_scores = []
        for cell_type in MARKERS:
            gene_scores = []
            for marker in available[cell_type]:
                marker_counts = counts[:, marker_columns[marker]]
                if marker_counts.shape[1] > 1:
                    marker_counts = marker_counts.sum(axis=1)
                values = np.asarray(marker_counts.toarray()).ravel()
                gene_scores.append(np.log1p(values * (10_000.0 / library_size)))
            if gene_scores:
                type_scores.append(np.mean(np.vstack(gene_scores), axis=0))
            else:
                type_scores.append(np.zeros(len(subset), dtype=float))

        score_matrix = np.column_stack(type_scores)
        max_score = score_matrix.max(axis=1)
        top = np.isclose(score_matrix, max_score[:, None], rtol=1e-12, atol=1e-12)
        unique = top.sum(axis=1) == 1
        winner = np.argmax(score_matrix, axis=1)
        assigned = np.full(len(subset), "Other/unresolved", dtype=object)
        for position, cell_type in enumerate(MARKERS):
            assigned[unique & (max_score > 0) & (winner == position)] = cell_type
        subtype_assignments.loc[subset.index] = assigned
    data.file.close()

    coords["cell_type"] = subtype_assignments.to_numpy()
    coords["lineage"] = coords["cell_type"].map(SUBTYPE_TO_LINEAGE)
    coords["primary_group"] = np.where(
        coords["original_group"].astype(str) == "SV", "SV", "LV"
    )
    if coords["lineage"].isna().any():
        raise ValueError("A cell-type assignment did not map to a broad lineage.")
    keep_columns = [
        "object_cell_id", "sample_id", "original_group", "primary_group",
        "cell_type", "lineage", "n_counts", "n_genes", "pct_counts_mt",
        "UMAP1", "UMAP2",
    ]
    annotations = coords[keep_columns].rename(columns={"object_cell_id": "cell_id"})
    annotations.to_csv(OUT / "cell_type_annotation_table.csv", index=False)
    annotations[["cell_id", "sample_id", "original_group", "primary_group", "lineage", "UMAP1", "UMAP2"]].to_csv(
        OUT / "umap_coordinates_annotated.csv", index=False
    )
    marker_rows = []
    for cell_type, genes in MARKERS.items():
        for marker in genes:
            marker_rows.append(
                {
                    "cell_type": cell_type,
                    "lineage": SUBTYPE_TO_LINEAGE[cell_type],
                    "marker_gene": marker,
                    "detected_in_feature_table": bool(marker_columns[marker]),
                    "source": MARKER_SOURCE,
                }
            )
    pd.DataFrame(marker_rows).to_csv(OUT / "lineage_marker_genes_and_sources.csv", index=False)

    lineage_order = list(dict.fromkeys(SUBTYPE_TO_LINEAGE.values()))
    sample_groups = annotations[["sample_id", "original_group", "primary_group"]].drop_duplicates()
    counts = (
        annotations.groupby(["sample_id", "lineage"], observed=True)
        .size()
        .rename("n_qc_cells")
        .reset_index()
    )
    sample_lineage = (
        sample_groups.assign(_key=1)
        .merge(pd.DataFrame({"lineage": lineage_order, "_key": 1}), on="_key")
        .drop(columns="_key")
        .merge(counts, on=["sample_id", "lineage"], how="left", validate="one_to_one")
        .fillna({"n_qc_cells": 0})
    )
    sample_lineage["n_qc_cells"] = sample_lineage["n_qc_cells"].astype(int)
    sample_lineage.to_csv(OUT / "sample_lineage_cell_counts.csv", index=False)

    print(
        f"Annotated {len(annotations)} QC-retained cells across "
        f"{annotations['sample_id'].nunique()} samples; lineages: "
        f"{annotations['lineage'].value_counts().to_dict()}"
    )
    print("Markers unavailable by subtype:", missing)


if __name__ == "__main__":
    main()
