"""Profile consensus hubs by lineage and aggregate raw counts by sample."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "processed" / "GSE248762_hubblind_allcells_qc.h5ad"
HUBS = [
    "TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19",
    "NAT6", "ECM1", "CST7", "CCBE1",
]
MODEL_PANEL = {"TNFSF15", "ADAM19", "FLT3LG", "NAT6", "ECM1"}
SEVEN_GENE_PANEL = {
    "FLT3LG", "TNFSF15", "LTB", "ADAM19", "SERPINA10", "EBI3", "CXCL14"
}


def as_csr(matrix: object) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(matrix)


def main() -> None:
    annotation = pd.read_csv(OUT / "cell_type_annotation_table.csv")
    sample_audit = pd.read_csv(OUT / "sample_group_qc_audit.csv")
    lineages = pd.read_csv(OUT / "sample_lineage_cell_counts.csv")["lineage"].drop_duplicates().tolist()
    membership = pd.DataFrame(
        {
            "gene": HUBS,
            "in_project_9_gene_hubs": True,
            "in_five_gene_model_panel": [gene in MODEL_PANEL for gene in HUBS],
            "in_seven_gene_consensus_panel": [gene in SEVEN_GENE_PANEL for gene in HUBS],
        }
    )
    membership.to_csv(OUT / "hub_gene_panel_membership.csv", index=False)

    data = ad.read_h5ad(DATA, backed="r")
    symbols = pd.Index(data.var_names.astype(str))
    normalized_symbols = symbols.str.upper()
    gene_ids = data.var["gene_ids"].astype(str).to_numpy()
    if pd.Index(gene_ids).duplicated().any():
        data.file.close()
        raise ValueError("Feature IDs are duplicated in the processed GEO object.")
    feature_map = pd.DataFrame(
        {"feature_id": gene_ids, "gene_symbol": symbols.to_numpy()}
    )
    feature_map.to_csv(OUT / "gene_id_symbol_map.csv", index=False)
    symbol_indices = {
        gene: np.flatnonzero(normalized_symbols == gene).tolist() for gene in HUBS
    }
    feature_presence = pd.DataFrame(
        {
            "gene": HUBS,
            "feature_detected": [bool(symbol_indices[gene]) for gene in HUBS],
            "feature_ids": [
                ";".join(gene_ids[index] for index in symbol_indices[gene])
                for gene in HUBS
            ],
            "feature_symbols": [
                ";".join(symbols[index] for index in symbol_indices[gene])
                for gene in HUBS
            ],
        }
    )
    feature_presence.to_csv(OUT / "hub_gene_feature_presence.csv", index=False)

    cell_to_row = pd.Series(
        np.arange(data.n_obs, dtype=np.int64), index=data.obs_names.astype(str)
    )
    if annotation["cell_id"].duplicated().any():
        data.file.close()
        raise ValueError("Cell annotation table contains duplicate cell identifiers.")
    annotation["object_row"] = annotation["cell_id"].map(cell_to_row)
    if annotation["object_row"].isna().any():
        data.file.close()
        raise ValueError("Some annotated cells are absent from the expression object.")
    annotation["object_row"] = annotation["object_row"].astype(np.int64)

    profile_accumulator = {
        lineage: {
            gene: {"cells": 0, "detected": 0, "log_expression_sum": 0.0}
            for gene in HUBS
        }
        for lineage in lineages
    }
    pseudobulk: dict[str, dict[str, np.ndarray]] = {
        lineage: {
            sample_id: np.zeros(data.n_vars, dtype=np.int64)
            for sample_id in sample_audit["sample_id"].astype(str)
        }
        for lineage in lineages
    }
    expression_rows = []

    for sample_id, sample_cells in annotation.groupby("sample_id", sort=False):
        rows = np.sort(sample_cells["object_row"].to_numpy(dtype=np.int64))
        sample_counts = as_csr(data.X[rows, :])
        if sample_counts.data.size and (
            np.any(sample_counts.data < 0)
            or not np.all(np.equal(sample_counts.data, np.floor(sample_counts.data)))
        ):
            data.file.close()
            raise ValueError(f"{sample_id}: input feature matrix is not raw integer counts.")
        library_size = np.asarray(sample_counts.sum(axis=1)).ravel().astype(float)
        if np.any(library_size <= 0):
            data.file.close()
            raise ValueError(f"{sample_id}: a QC-retained cell has zero library size.")

        sample_cells = sample_cells.sort_values("object_row").reset_index(drop=True)
        if not np.array_equal(sample_cells["object_row"].to_numpy(), rows):
            data.file.close()
            raise ValueError(f"{sample_id}: annotation/count row order mismatch.")
        group = str(sample_cells["primary_group"].iloc[0])
        hub_expression: dict[str, np.ndarray] = {}
        for gene in HUBS:
            indices = symbol_indices[gene]
            if not indices:
                values = np.full(len(sample_cells), np.nan, dtype=float)
            else:
                raw_gene = sample_counts[:, indices].sum(axis=1)
                raw_gene = np.asarray(raw_gene).ravel().astype(float)
                values = np.log1p(raw_gene * (10_000.0 / library_size))
            hub_expression[gene] = values

        per_cell = sample_cells[["cell_id", "sample_id", "original_group", "primary_group", "lineage"]].copy()
        for gene in HUBS:
            per_cell[gene] = hub_expression[gene]
        expression_rows.append(per_cell)

        for lineage in lineages:
            local_positions = np.flatnonzero(sample_cells["lineage"].to_numpy() == lineage)
            if local_positions.size == 0:
                continue
            lineage_counts = sample_counts[local_positions, :]
            pseudobulk[lineage][str(sample_id)] = (
                np.asarray(lineage_counts.sum(axis=0)).ravel().astype(np.int64)
            )
            for gene in HUBS:
                values = hub_expression[gene][local_positions]
                stats = profile_accumulator[lineage][gene]
                stats["cells"] += int(local_positions.size)
                if np.isfinite(values).all():
                    raw_counts = sample_counts[local_positions, :][:, symbol_indices[gene]]
                    detected = np.asarray(raw_counts.sum(axis=1)).ravel() > 0
                    stats["detected"] += int(detected.sum())
                    stats["log_expression_sum"] += float(values.sum())
    data.file.close()

    profile_rows = []
    for lineage in lineages:
        for gene in HUBS:
            stats = profile_accumulator[lineage][gene]
            n_cells = int(stats["cells"])
            present = bool(symbol_indices[gene])
            profile_rows.append(
                {
                    "gene": gene,
                    "lineage": lineage,
                    "feature_detected": present,
                    "cells_in_lineage": n_cells,
                    "expressing_cells_raw_count_gt_0": (
                        int(stats["detected"]) if present else pd.NA
                    ),
                    "percent_expressing_cells": (
                        100.0 * stats["detected"] / n_cells
                        if present and n_cells else pd.NA
                    ),
                    "mean_log1p_CPM10k": (
                        stats["log_expression_sum"] / n_cells
                        if present and n_cells else pd.NA
                    ),
                }
            )
    pd.DataFrame(profile_rows).to_csv(OUT / "hub_gene_expression_by_lineage.csv", index=False)
    pd.concat(expression_rows, ignore_index=True).to_csv(
        OUT / "hub_gene_expression_per_cell.csv", index=False
    )

    count_dir = OUT / "pseudobulk_counts"
    count_dir.mkdir(exist_ok=True)
    for lineage in lineages:
        safe_lineage = (
            lineage.lower().replace("/", "_").replace(" ", "_").replace("-", "_")
        )
        matrix = pd.DataFrame(
            {sample_id: values for sample_id, values in pseudobulk[lineage].items()}
        )
        matrix.insert(0, "gene_symbol", symbols.to_numpy())
        matrix.insert(0, "feature_id", gene_ids)
        matrix.to_csv(count_dir / f"counts_{safe_lineage}.csv", index=False)
    print(
        f"Generated expression profiles and pseudobulk matrices for {len(HUBS)} "
        f"hub genes, {len(lineages)} lineages, and {len(sample_audit)} samples."
    )


if __name__ == "__main__":
    main()
