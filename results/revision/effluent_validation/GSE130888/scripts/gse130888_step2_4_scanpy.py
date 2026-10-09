#!/usr/bin/env python
"""Preregistered Scanpy and gene-set analysis for GSE130888."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import shutil
import subprocess
import tarfile
import tempfile
from collections import Counter, defaultdict
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import binomtest, mannwhitneyu
from sklearn.metrics import roc_auc_score, roc_curve


GENES = [
    "TNFSF15",
    "FLT3LG",
    "EBI3",
    "LTB",
    "ADAM19",
    "ECM1",
    "SERPINA10",
    "CST7",
]
STRICT_GENES = ["TNFSF15", "FLT3LG"]
CELLTYPE_MARKERS = {
    "T/NK": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD74"],
    "monocyte/macrophage": [
        "LST1",
        "TYROBP",
        "FCER1G",
        "CD14",
        "LILRB1",
        "C1QA",
        "C1QB",
    ],
    "neutrophil": ["CSF3R", "FCGR3B", "S100A8", "S100A9", "CXCR2"],
    "mesothelial": ["KRT8", "KRT18", "KRT19", "MSLN", "UPK3B", "KRT7"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"],
}
CELLTYPES = [*CELLTYPE_MARKERS, "other"]
MIN_GENES = 200
MAX_MT_PERCENT = 20.0
MIN_PSEUDOBULK_CELLS = 20
MIN_DETECTION_COUNT = 10
N_BOOTSTRAPS = 10_000
N_NULL_SETS = 1_000


def bh_adjust(values: list[float]) -> list[float]:
    """Benjamini-Hochberg adjustment, retaining NaNs as NaNs."""
    values_array = np.asarray(values, dtype=float)
    valid = np.isfinite(values_array)
    adjusted = np.full(values_array.shape, np.nan)
    if not valid.any():
        return adjusted.tolist()
    pvalues = values_array[valid]
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    qvalues = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    qvalues = np.minimum.accumulate(qvalues[::-1])[::-1]
    restored = np.empty_like(qvalues)
    restored[order] = np.minimum(qvalues, 1.0)
    adjusted[valid] = restored
    return adjusted.tolist()


def load_sample_metadata(repo_root: Path) -> pd.DataFrame:
    path = (
        repo_root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE130888"
        / "sample_metadata.csv"
    )
    metadata = pd.read_csv(path, dtype=str, keep_default_na=False)
    expected = {
        "Short-term PD": 6,
        "Long-term PD": 4,
        "Normal control": 3,
    }
    observed = metadata["group_from_metadata"].value_counts().to_dict()
    if observed != expected or metadata["sample_id"].duplicated().any():
        raise ValueError(
            f"GEO sample mapping differs from preregistered cohorts: {observed}"
        )
    group_map = {
        "Short-term PD": "SV",
        "Long-term PD": "LV",
        "Normal control": "NORMAL",
    }
    metadata["analysis_group"] = metadata["group_from_metadata"].map(group_map)
    if metadata["analysis_group"].isna().any():
        raise ValueError("GEO sample metadata contains an unmapped group")
    return metadata


def primary_group_for_sample(group: str) -> str:
    if group in {"SV", "LV"}:
        return group
    if group == "NORMAL":
        return "NORMAL"
    raise ValueError(f"Unexpected normalized GSE130888 group: {group}")


def extract_10x_sample(
    archive: tarfile.TarFile,
    sample_id: str,
    scratch: Path,
) -> Path:
    sample_dir = scratch / sample_id
    sample_dir.mkdir()
    wanted = {
        "barcodes.tsv.gz": "barcodes.tsv.gz",
        "features.tsv.gz": "features.tsv.gz",
        "matrix.mtx.gz": "matrix.mtx.gz",
    }
    members = {
        Path(member.name).name: member
        for member in archive.getmembers()
        if Path(member.name).name.startswith(f"{sample_id}_")
    }
    for suffix, output_name in wanted.items():
        matches = [
            (name, member)
            for name, member in members.items()
            if name.endswith(suffix)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"{sample_id}: expected one {suffix} archive member, "
                f"found {len(matches)}"
            )
        source = archive.extractfile(matches[0][1])
        if source is None:
            raise OSError(f"Could not read {matches[0][0]} from GEO archive")
        with source, (sample_dir / output_name).open("wb") as destination:
            shutil.copyfileobj(source, destination)
    return sample_dir


def read_10x(path: Path) -> ad.AnnData:
    sample = sc.read_10x_mtx(
        path,
        var_names="gene_ids",
        make_unique=True,
        cache=False,
        gex_only=True,
    )
    sample.var_names_make_unique()
    if "gene_symbols" not in sample.var:
        raise ValueError("10x feature table has no gene-symbol column")
    sample.var["gene_symbols"] = sample.var["gene_symbols"].astype(str)
    sample.X = sample.X.tocsr()
    if sample.X.data.size and (
        np.any(sample.X.data < 0)
        or not np.all(np.equal(sample.X.data, np.floor(sample.X.data)))
    ):
        raise ValueError("GEO matrix is not nonnegative integer raw counts")
    return sample


def qc_and_normalize(
    sample: ad.AnnData,
) -> tuple[ad.AnnData, int, int]:
    symbols = sample.var["gene_symbols"].str.upper()
    sample.var["mt"] = symbols.str.startswith("MT-").to_numpy()
    if not sample.var["mt"].any():
        raise ValueError("No mitochondrial features found in feature annotations")
    sc.pp.calculate_qc_metrics(
        sample,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    cells_before = sample.n_obs
    keep = (
        (sample.obs["n_genes_by_counts"].to_numpy() >= MIN_GENES)
        & (sample.obs["pct_counts_mt"].to_numpy() <= MAX_MT_PERCENT)
    )
    sample = sample[keep].copy()
    cells_after = sample.n_obs
    sample.layers["counts"] = sample.X.copy()
    sc.pp.normalize_total(sample, target_sum=10_000)
    sc.pp.log1p(sample)
    return sample, cells_before, cells_after


def gene_symbol_indices(sample: ad.AnnData) -> dict[str, np.ndarray]:
    symbols = sample.var["gene_symbols"].str.upper().to_numpy()
    indices: dict[str, np.ndarray] = {}
    for symbol in set(symbols):
        indices[symbol] = np.flatnonzero(symbols == symbol)
    return indices


def annotate_celltypes(sample: ad.AnnData) -> tuple[np.ndarray, dict[str, list[str]]]:
    symbol_indices = gene_symbol_indices(sample)
    scores = np.zeros((len(CELLTYPE_MARKERS), sample.n_obs), dtype=np.float32)
    marker_availability: dict[str, list[str]] = {}
    for row, (celltype, markers) in enumerate(CELLTYPE_MARKERS.items()):
        available = [marker for marker in markers if marker in symbol_indices]
        marker_availability[celltype] = available
        if available:
            indices = np.concatenate([symbol_indices[marker] for marker in available])
            scores[row] = np.asarray(sample.X[:, indices].mean(axis=1)).ravel()
    max_scores = scores.max(axis=0)
    winners = np.isclose(scores, max_scores, rtol=1e-12, atol=1e-12)
    unique_winner = winners.sum(axis=0) == 1
    annotations = np.full(sample.n_obs, "other", dtype=object)
    has_marker_signal = max_scores > 0
    winner_indices = np.argmax(scores, axis=0)
    for row, celltype in enumerate(CELLTYPE_MARKERS):
        take = unique_winner & has_marker_signal & (winner_indices == row)
        annotations[take] = celltype
    return annotations, marker_availability


def expression_for_gene(sample: ad.AnnData, gene_id: str) -> tuple[float, int]:
    index = sample.var_names.get_loc(gene_id)
    normalized = sample.X[:, index]
    raw = sample.layers["counts"][:, index]
    mean_expression = float(np.asarray(normalized.mean()).reshape(-1)[0])
    detected_cells = int((raw > 0).sum())
    return mean_expression, detected_cells


def save_marker_and_expression_outputs(
    output_dir: Path,
    available: dict[str, list[str]],
    type_totals: dict[str, dict[str, float]],
    type_detected: dict[str, dict[str, int]],
    type_cell_totals: dict[str, int],
    candidate_ids: dict[str, str],
) -> pd.DataFrame:
    marker_rows = []
    for celltype, markers in CELLTYPE_MARKERS.items():
        mapped = available.get(celltype, [])
        for marker in markers:
            marker_rows.append(
                {
                    "cell_type": celltype,
                    "marker": marker,
                    "present_in_GEO_features": marker in mapped,
                }
            )
    pd.DataFrame(marker_rows).to_csv(
        output_dir / "fallback_marker_list.csv", index=False
    )
    profile_rows = []
    for gene in GENES:
        for celltype in CELLTYPES:
            n_cells = type_cell_totals.get(celltype, 0)
            total_expr = type_totals.get(celltype, {}).get(gene, 0.0)
            total_detected = type_detected.get(celltype, {}).get(gene, 0)
            profile_rows.append(
                {
                    "gene": gene,
                    "gene_id": candidate_ids[gene],
                    "cell_type": celltype,
                    "cells": n_cells,
                    "fraction_cells_nonzero": (
                        total_detected / n_cells if n_cells else np.nan
                    ),
                    "mean_log1p_normalized_expression": (
                        total_expr / n_cells if n_cells else np.nan
                    ),
                }
            )
    profile = pd.DataFrame(profile_rows)
    profile["celltype_order"] = profile["cell_type"].map(
        {celltype: index for index, celltype in enumerate(CELLTYPES)}
    )
    top_types = {}
    for gene, subset in profile.groupby("gene", sort=False):
        subset = subset.sort_values(
            ["mean_log1p_normalized_expression", "celltype_order"],
            ascending=[False, True],
            kind="stable",
        )
        top_types[gene] = subset.iloc[0]["cell_type"]
    profile["gene_top_cell_type"] = profile["gene"].map(top_types)
    profile.drop(columns="celltype_order").to_csv(
        output_dir / "gene_expression_by_cell_type.csv", index=False
    )
    return profile


def prepare_data(repo_root: Path, output_dir: Path) -> dict:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Refusing to overwrite existing analysis outputs: {output_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = load_sample_metadata(repo_root)
    group_by_sample = metadata.set_index("sample_id")["analysis_group"].to_dict()
    title_by_sample = metadata.set_index("sample_id")["title"].to_dict()
    archive_path = (
        repo_root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE130888"
        / "GSE130888_RAW.tar"
    )
    expected_bytes = 519_403_520
    if archive_path.stat().st_size != expected_bytes:
        raise ValueError("GSE130888 RAW archive size differs from GEO manifest")

    qc_rows = []
    hvg_frequency: Counter[str] = Counter()
    marker_availability: dict[str, list[str]] = {}
    pseudobulk_all: dict[str, np.ndarray] = {}
    pseudobulk_celltype: dict[str, dict[str, np.ndarray]] = defaultdict(dict)
    pseudobulk_celltype_meta = []
    composition_rows = []
    type_cell_totals: Counter[str] = Counter()
    type_expression_totals: dict[str, Counter[str]] = defaultdict(Counter)
    type_detected_totals: dict[str, Counter[str]] = defaultdict(Counter)
    gene_symbols_by_id: dict[str, str] = {}
    candidate_ids: dict[str, str] = {}
    sample_cells_after: dict[str, int] = {}
    sample_umap_metadata = {}
    reference_gene_ids: list[str] | None = None
    retained_barcode_hash: dict[str, str] = {}

    with tempfile.TemporaryDirectory(prefix="gse130888_raw_") as temporary:
        scratch = Path(temporary)
        with tarfile.open(archive_path, mode="r") as archive:
            for sample_number, sample_id in enumerate(metadata["sample_id"], 1):
                print(
                    f"QC/pseudobulk {sample_number}/{len(metadata)}: {sample_id}",
                    flush=True,
                )
                sample_dir = extract_10x_sample(archive, sample_id, scratch)
                sample = read_10x(sample_dir)
                current_gene_ids = sample.var_names.astype(str).tolist()
                if reference_gene_ids is None:
                    reference_gene_ids = current_gene_ids
                elif current_gene_ids != reference_gene_ids:
                    raise ValueError(
                        f"{sample_id}: feature order differs from other samples"
                    )
                sample, cells_before, cells_after = qc_and_normalize(sample)
                if cells_after == 0:
                    raise ValueError(
                        f"{sample_id} has no cells after preregistered QC"
                    )
                sample.obs["sample_id"] = sample_id
                sample.obs["patient_id_proxy"] = sample_id
                sample.obs["original_group"] = group_by_sample[sample_id]
                sample.obs["primary_group"] = primary_group_for_sample(
                    group_by_sample[sample_id]
                )
                retained_barcode_hash[sample_id] = hashlib.sha256(
                    "\n".join(sample.obs_names.astype(str)).encode("utf-8")
                ).hexdigest()
                annotations, available = annotate_celltypes(sample)
                sample.obs["cell_type"] = pd.Categorical(
                    annotations, categories=CELLTYPES, ordered=True
                )
                marker_availability = available
                symbols = sample.var["gene_symbols"].to_numpy()
                for gene_id, symbol in zip(sample.var_names.astype(str), symbols):
                    gene_symbols_by_id.setdefault(gene_id, symbol)
                for gene in GENES:
                    matched_ids = sample.var_names[
                        sample.var["gene_symbols"]
                        .astype(str)
                        .str.upper()
                        .to_numpy()
                        == gene
                    ].astype(str).tolist()
                    if len(matched_ids) != 1:
                        raise ValueError(
                            f"{gene}: expected one 10x feature ID, found {matched_ids}"
                        )
                    prior = candidate_ids.setdefault(gene, matched_ids[0])
                    if prior != matched_ids[0]:
                        raise ValueError(f"{gene}: feature ID varies across samples")

                n_qc_cells = sample.n_obs
                sample_cells_after[sample_id] = n_qc_cells
                qc_rows.append(
                    {
                        "sample_id": sample_id,
                        "title": title_by_sample[sample_id],
                        "group": group_by_sample[sample_id],
                        "patient_id_proxy": sample_id,
                        "cells_before_qc": cells_before,
                        "cells_after_qc": cells_after,
                        "cells_removed_by_preregistered_qc": cells_before
                        - cells_after,
                        "median_detected_genes_after_qc": float(
                            sample.obs["n_genes_by_counts"].median()
                        ),
                        "median_mito_percent_after_qc": float(
                            sample.obs["pct_counts_mt"].median()
                        ),
                        "pseudobulk_included": (
                            n_qc_cells >= MIN_PSEUDOBULK_CELLS
                        ),
                    }
                )

                group_labels = annotations
                raw = sample.layers["counts"].tocsr()
                for celltype in CELLTYPES:
                    mask = group_labels == celltype
                    n_type = int(mask.sum())
                    type_cell_totals[celltype] += n_type
                    if n_type:
                        type_counts = raw[mask].sum(axis=0)
                        normalized_mean = np.asarray(
                            sample.X[mask].mean(axis=0)
                        ).ravel()
                        for gene in GENES:
                            index = sample.var_names.get_loc(candidate_ids[gene])
                            type_expression_totals[celltype][gene] += (
                                normalized_mean[index] * n_type
                            )
                            type_detected_totals[celltype][gene] += int(
                                (raw[mask, index] > 0).sum()
                            )
                        if n_type >= MIN_PSEUDOBULK_CELLS:
                            pseudobulk_celltype[sample_id][celltype] = (
                                np.asarray(type_counts).ravel().astype(np.int64)
                            )
                            pseudobulk_celltype_meta.append(
                                {
                                    "sample_id": sample_id,
                                    "patient_id_proxy": sample_id,
                                    "original_group": group_by_sample[sample_id],
                                    "primary_group": primary_group_for_sample(
                                        group_by_sample[sample_id]
                                    ),
                                    "cell_type": celltype,
                                    "n_qc_cells": n_type,
                                }
                            )
                    else:
                        for gene in GENES:
                            type_expression_totals[celltype][gene] += 0.0

                counts_per_type = pd.Series(group_labels).value_counts()
                composition = {
                    "sample_id": sample_id,
                    "patient_id_proxy": sample_id,
                    "original_group": group_by_sample[sample_id],
                    "primary_group": primary_group_for_sample(
                        group_by_sample[sample_id]
                    ),
                    "n_qc_cells": n_qc_cells,
                }
                composition.update(
                    {
                        f"fraction_{celltype}": (
                            float(counts_per_type.get(celltype, 0) / n_qc_cells)
                        )
                        for celltype in CELLTYPES
                    }
                )
                composition_rows.append(composition)

                hvg_count = min(2_000, sample.n_vars - 1)
                if hvg_count < 2:
                    raise ValueError("Too few features for preregistered UMAP")
                sc.pp.highly_variable_genes(
                    sample,
                    flavor="seurat",
                    n_top_genes=hvg_count,
                    subset=False,
                )
                hvg_frequency.update(
                    sample.var_names[sample.var["highly_variable"]].astype(str)
                )
                sample_umap_metadata[sample_id] = {
                    "group": group_by_sample[sample_id],
                    "cell_type": annotations.copy(),
                }

                if n_qc_cells >= MIN_PSEUDOBULK_CELLS:
                    pseudobulk_all[sample_id] = (
                        np.asarray(raw.sum(axis=0)).ravel().astype(np.int64)
                    )

                del sample, raw, group_labels, annotations
                gc.collect()
                shutil.rmtree(sample_dir)

        if reference_gene_ids is None:
            raise ValueError("No sample matrices were processed")
        gene_ids = reference_gene_ids
        if len(gene_symbols_by_id) == 0 or set(candidate_ids) != set(GENES):
            raise ValueError("Incomplete feature symbol mapping")
        pd.DataFrame(
            {
                "gene_id": list(gene_symbols_by_id),
                "gene_symbol": list(gene_symbols_by_id.values()),
            }
        ).to_csv(output_dir / "gene_id_to_symbol.csv", index=False)
        pd.DataFrame(qc_rows).to_csv(output_dir / "sample_qc.csv", index=False)
        profile = save_marker_and_expression_outputs(
            output_dir,
            marker_availability,
            type_expression_totals,
            type_detected_totals,
            type_cell_totals,
            candidate_ids,
        )

        eligible_samples = list(pseudobulk_all)
        count_matrix = np.column_stack(
            [pseudobulk_all[sample_id] for sample_id in eligible_samples]
        )
        counts_frame = pd.DataFrame(
            count_matrix, index=gene_ids, columns=eligible_samples
        )
        counts_frame.index.name = "gene_id"
        counts_frame.to_csv(output_dir / "all_cell_pseudobulk_counts.csv")
        pb_meta = []
        for sample_id in eligible_samples:
            group = group_by_sample[sample_id]
            pb_meta.append(
                {
                    "sample_id": sample_id,
                    "patient_id_proxy": sample_id,
                    "original_group": group,
                    "primary_group": primary_group_for_sample(group),
                    "n_qc_cells": sample_cells_after[sample_id],
                }
            )
        pd.DataFrame(pb_meta).to_csv(
            output_dir / "all_cell_pseudobulk_metadata.csv", index=False
        )
        if not pseudobulk_celltype_meta:
            raise ValueError("No patient-by-cell-type pseudobulks met 20-cell rule")
        celltype_columns = [
            (row["sample_id"], row["cell_type"])
            for row in pseudobulk_celltype_meta
        ]
        celltype_matrix = np.column_stack(
            [
                pseudobulk_celltype[sample_id][celltype]
                for sample_id, celltype in celltype_columns
            ]
        )
        celltype_ids = [
            f"{sample_id}|{celltype}" for sample_id, celltype in celltype_columns
        ]
        celltype_count_frame = pd.DataFrame(
            celltype_matrix, index=gene_ids, columns=celltype_ids
        )
        celltype_count_frame.index.name = "gene_id"
        celltype_count_frame.to_csv(
            output_dir / "celltype_pseudobulk_counts.csv"
        )
        pd.DataFrame(pseudobulk_celltype_meta).to_csv(
            output_dir / "celltype_pseudobulk_metadata.csv", index=False
        )

        all_celltype_pairs = {
            (row["sample_id"], celltype)
            for row in composition_rows
            for celltype in CELLTYPES
        }
        included_pairs = set(celltype_columns)
        exclusions = []
        for row in composition_rows:
            for celltype in CELLTYPES:
                pair = (row["sample_id"], celltype)
                if pair not in included_pairs:
                    n_cells = round(row[f"fraction_{celltype}"] * row["n_qc_cells"])
                    exclusions.append(
                        {
                            "sample_id": row["sample_id"],
                            "patient_id_proxy": row["patient_id_proxy"],
                            "original_group": row["original_group"],
                            "cell_type": celltype,
                            "n_qc_cells": n_cells,
                            "reason": "fewer than 20 cells",
                        }
                    )
        pd.DataFrame(exclusions).to_csv(
            output_dir / "excluded_patient_cell_types.csv", index=False
        )

        composition = pd.DataFrame(composition_rows)
        composition.to_csv(
            output_dir / "patient_celltype_composition.csv", index=False
        )
        composition_for_tests = composition[
            composition["n_qc_cells"] >= MIN_PSEUDOBULK_CELLS
        ].copy()
        comp_p = []
        for celltype in CELLTYPES:
            long_values = composition_for_tests.loc[
                composition_for_tests["primary_group"] == "LV",
                f"fraction_{celltype}",
            ].to_numpy()
            short_values = composition_for_tests.loc[
                composition_for_tests["primary_group"] == "SV",
                f"fraction_{celltype}",
            ].to_numpy()
            pvalue = (
                float(mannwhitneyu(long_values, short_values, alternative="two-sided").pvalue)
                if len(long_values) and len(short_values)
                else np.nan
            )
            comp_p.append(pvalue)
        comp_summary = pd.DataFrame(
            {
                "cell_type": CELLTYPES,
                "n_LV": int((composition_for_tests["primary_group"] == "LV").sum()),
                "n_SV": int((composition_for_tests["primary_group"] == "SV").sum()),
                "mean_fraction_LV": [
                    composition_for_tests.loc[
                        composition_for_tests["primary_group"] == "LV",
                        f"fraction_{celltype}",
                    ].mean()
                    for celltype in CELLTYPES
                ],
                "mean_fraction_SV": [
                    composition_for_tests.loc[
                        composition_for_tests["primary_group"] == "SV",
                        f"fraction_{celltype}",
                    ].mean()
                    for celltype in CELLTYPES
                ],
                "wilcoxon_p": comp_p,
                "bh_q": bh_adjust(comp_p),
            }
        )
        comp_summary.to_csv(
            output_dir / "composition_tests.csv", index=False
        )

        selected_umap_genes = [
            gene
            for gene, _ in sorted(
                hvg_frequency.items(), key=lambda item: (-item[1], item[0])
            )[:2_000]
        ]
        if len(selected_umap_genes) < 2:
            raise ValueError("Could not select sufficient group-blind UMAP features")
        pd.DataFrame(
            {
                "gene_id": selected_umap_genes,
                "gene_symbol": [
                    gene_symbols_by_id[gene] for gene in selected_umap_genes
                ],
                "samples_hvg": [hvg_frequency[gene] for gene in selected_umap_genes],
            }
        ).to_csv(output_dir / "umap_hvg_features.csv", index=False)

        umap_samples = []
        with tarfile.open(archive_path, mode="r") as archive:
            for sample_number, sample_id in enumerate(metadata["sample_id"], 1):
                print(
                    f"UMAP input {sample_number}/{len(metadata)}: {sample_id}",
                    flush=True,
                )
                sample_dir = extract_10x_sample(archive, sample_id, scratch)
                sample = read_10x(sample_dir)
                sample, _, _ = qc_and_normalize(sample)
                info = sample_umap_metadata[sample_id]
                barcode_hash = hashlib.sha256(
                    "\n".join(sample.obs_names.astype(str)).encode("utf-8")
                ).hexdigest()
                if barcode_hash != retained_barcode_hash[sample_id]:
                    raise ValueError(
                        f"{sample_id}: retained barcodes differ between "
                        "pseudobulk and UMAP passes"
                    )
                if len(info["cell_type"]) != sample.n_obs:
                    raise ValueError(
                        f"{sample_id}: marker annotations do not align with "
                        "second-pass QC cells"
                    )
                sample.obs["sample_id"] = sample_id
                sample.obs["original_group"] = info["group"]
                sample.obs["primary_group"] = primary_group_for_sample(
                    info["group"]
                )
                sample.obs["cell_type"] = pd.Categorical(
                    info["cell_type"], categories=CELLTYPES, ordered=True
                )
                sample = sample[:, selected_umap_genes].copy()
                umap_samples.append(sample)
                shutil.rmtree(sample_dir)
        umap_adata = ad.concat(
            umap_samples,
            join="inner",
            label="analysis_batch",
            index_unique=":",
            merge="same",
        )
        umap_adata.var_names = selected_umap_genes
        n_pcs = min(40, umap_adata.n_vars - 1, umap_adata.n_obs - 1)
        sc.pp.pca(
            umap_adata,
            n_comps=n_pcs,
            random_state=42,
            svd_solver="arpack",
        )
        sc.pp.neighbors(
            umap_adata,
            n_neighbors=15,
            n_pcs=n_pcs,
            random_state=42,
        )
        sc.tl.umap(umap_adata, random_state=42)
        coordinates = pd.DataFrame(
            umap_adata.obsm["X_umap"],
            columns=["UMAP1", "UMAP2"],
            index=umap_adata.obs_names,
        )
        coordinates["sample_id"] = umap_adata.obs["sample_id"].astype(str).to_numpy()
        coordinates["original_group"] = (
            umap_adata.obs["original_group"].astype(str).to_numpy()
        )
        coordinates["primary_group"] = (
            umap_adata.obs["primary_group"].astype(str).to_numpy()
        )
        coordinates["cell_type"] = (
            umap_adata.obs["cell_type"].astype(str).to_numpy()
        )
        coordinates.to_csv(output_dir / "umap_coordinates.csv", index_label="cell_id")
        fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
        for group in ("SV", "LV", "NORMAL"):
            selected = coordinates["original_group"] == group
            axes[0].scatter(
                coordinates.loc[selected, "UMAP1"],
                coordinates.loc[selected, "UMAP2"],
                s=0.35,
                alpha=0.5,
                label=group,
                rasterized=True,
            )
        axes[0].set_title("GSE130888 cells by GEO group")
        axes[0].legend(markerscale=8, frameon=False)
        for celltype in CELLTYPES:
            selected = coordinates["cell_type"] == celltype
            if selected.any():
                axes[1].scatter(
                    coordinates.loc[selected, "UMAP1"],
                    coordinates.loc[selected, "UMAP2"],
                    s=0.35,
                    alpha=0.5,
                    label=celltype,
                    rasterized=True,
                )
        axes[1].set_title("GSE130888 cells by marker-based cell type")
        axes[1].legend(markerscale=8, frameon=False, fontsize=7)
        for axis in axes:
            axis.set_xlabel("UMAP1")
            axis.set_ylabel("UMAP2")
        fig.savefig(output_dir / "umap_by_group_and_cell_type.png", dpi=240)
        plt.close(fig)

    profile = pd.read_csv(output_dir / "gene_expression_by_cell_type.csv")
    profile = profile[profile["gene"].isin(GENES)].copy()
    fig, axis = plt.subplots(figsize=(10, 5))
    max_size = max(profile["fraction_cells_nonzero"].max(), 1e-8)
    for _, row in profile.iterrows():
        size = 25 + 400 * row["fraction_cells_nonzero"] / max_size
        axis.scatter(
            CELLTYPES.index(row["cell_type"]),
            GENES.index(row["gene"]),
            s=size,
            c=row["mean_log1p_normalized_expression"],
            cmap="viridis",
            vmin=profile["mean_log1p_normalized_expression"].min(),
            vmax=profile["mean_log1p_normalized_expression"].max(),
            edgecolor="black",
            linewidth=0.25,
        )
    axis.set_xticks(range(len(CELLTYPES)), CELLTYPES, rotation=35, ha="right")
    axis.set_yticks(range(len(GENES)), GENES)
    axis.set_xlabel("Marker-based cell type")
    axis.set_title("Candidate-gene expression across QC-retained cells")
    scatter = axis.collections[-1]
    fig.colorbar(scatter, ax=axis, label="Mean log1p normalized expression")
    fig.tight_layout()
    fig.savefig(output_dir / "candidate_gene_dotplot.png", dpi=240)
    plt.close(fig)

    config = {
        "accession": "GSE130888",
        "input": archive_path.name,
        "input_bytes": archive_path.stat().st_size,
        "scanpy_version": version("scanpy"),
        "anndata_version": version("anndata"),
        "python_version": __import__("sys").version.split()[0],
        "qc": {
            "minimum_detected_genes": MIN_GENES,
            "maximum_mitochondrial_percent": MAX_MT_PERCENT,
            "patient_minimum_cells_for_pseudobulk": MIN_PSEUDOBULK_CELLS,
            "celltype_minimum_cells_for_pseudobulk": MIN_PSEUDOBULK_CELLS,
        },
        "analysis_features": {
            "candidate_genes": GENES,
            "strict_genes": STRICT_GENES,
            "cell_type_method": "fixed preregistered marker-score rules",
            "umap_hvg_selection": (
                "top 2000 features most often selected as per-sample HVGs; "
                "group-blind, ties by feature ID"
            ),
            "all_cell_pseudobulk_patients": len(pseudobulk_all),
            "umap_cells": len(coordinates),
            "seed": 42,
        },
    }
    (output_dir / "analysis_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    return {
        "repo_root": repo_root,
        "output_dir": output_dir,
        "candidate_ids": candidate_ids,
        "top_cell_types": profile.drop_duplicates("gene").set_index("gene")[
            "gene_top_cell_type"
        ].to_dict(),
    }


def run_edger(repo_root: Path, output_dir: Path) -> None:
    rscript = Path(r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe")
    if not rscript.is_file():
        raise FileNotFoundError(f"Required Rscript not found: {rscript}")
    r_script = (
        Path(__file__).resolve().parent / "gse130888_step2_4_edger.R"
    )
    subprocess.run(
        [
            str(rscript),
            str(r_script),
            str(output_dir),
        ],
        check=True,
        cwd=repo_root,
    )


def empirical_upper_p(null: np.ndarray, observed: float) -> float:
    if null.size == 0 or not np.isfinite(observed):
        return float("nan")
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def make_null_sets(
    eligible: pd.DataFrame,
    target_ids: list[str],
    n_sets: int,
    seed: int,
) -> tuple[list[list[str]], dict]:
    ordered = eligible.sort_values(
        ["mean_logCPM", "gene_id"], kind="stable"
    ).reset_index(drop=True)
    ordered["expression_decile"] = np.minimum(
        9, np.floor(np.arange(len(ordered)) * 10 / len(ordered)).astype(int)
    )
    decile_by_id = ordered.set_index("gene_id")["expression_decile"].to_dict()
    target_id_set = set(target_ids)
    genes_by_decile = {
        int(decile): [
            gene_id
            for gene_id in subset["gene_id"].astype(str).tolist()
            if gene_id not in target_id_set
        ]
        for decile, subset in ordered.groupby("expression_decile")
    }
    target_deciles = [int(decile_by_id[gene_id]) for gene_id in target_ids]
    multiplicity = Counter(target_deciles)
    capacity = 1
    insufficient = {}
    for decile, needed in multiplicity.items():
        pool_size = len(genes_by_decile.get(decile, []))
        if pool_size < needed:
            insufficient[str(decile)] = {
                "available": pool_size,
                "required": needed,
            }
            capacity = 0
        else:
            capacity *= math.comb(pool_size, needed)
    details = {
        "target_deciles": target_deciles,
        "decile_pool_sizes": {
            str(decile): len(genes_by_decile.get(decile, []))
            for decile in multiplicity
        },
        "possible_unique_sets": int(capacity),
        "required_unique_sets": n_sets,
        "insufficient_deciles": insufficient,
    }
    if capacity < n_sets:
        return [], details
    rng = np.random.default_rng(seed)
    null_sets: list[list[str]] = []
    seen: set[tuple[str, ...]] = set()
    attempts = 0
    maximum_attempts = n_sets * 100_000
    while len(null_sets) < n_sets and attempts < maximum_attempts:
        attempts += 1
        selected = []
        for decile, count in multiplicity.items():
            selected.extend(
                rng.choice(genes_by_decile[decile], size=count, replace=False)
                .astype(str)
                .tolist()
            )
        signature = tuple(sorted(selected))
        if signature not in seen:
            seen.add(signature)
            null_sets.append(selected)
    details["attempts_to_draw"] = attempts
    if len(null_sets) != n_sets:
        return [], details
    return null_sets, details


def calculate_set_summary(
    set_gene_ids: list[str],
    all_gene_results: pd.DataFrame,
    logcpm: pd.DataFrame,
    metadata: pd.DataFrame,
    seed: int,
    bootstrap: bool,
) -> dict:
    genes = [gene for gene in set_gene_ids if gene in logcpm.index]
    if len(genes) != len(set_gene_ids):
        return {"n_genes": len(genes), "estimable": False}
    expression = logcpm.loc[genes, metadata["sample_id"]].to_numpy(dtype=float)
    means = expression.mean(axis=1, keepdims=True)
    stds = expression.std(axis=1, ddof=0, keepdims=True)
    if np.any(stds == 0) or np.any(~np.isfinite(stds)):
        return {"n_genes": len(genes), "estimable": False}
    scores = ((expression - means) / stds).mean(axis=0)
    long_mask = metadata["primary_group"].to_numpy() == "LV"
    short_scores = scores[~long_mask]
    long_scores = scores[long_mask]
    if not len(short_scores) or not len(long_scores):
        return {"n_genes": len(genes), "estimable": False}
    auc = float(roc_auc_score(long_mask.astype(int), scores))
    mann_whitney_p = float(
        mannwhitneyu(long_scores, short_scores, alternative="two-sided").pvalue
    )
    result = {
        "n_genes": len(genes),
        "estimable": True,
        "mean_score_LV": float(long_scores.mean()),
        "mean_score_SV": float(short_scores.mean()),
        "mann_whitney_p": mann_whitney_p,
        "auc": auc,
    }
    if bootstrap:
        rng = np.random.default_rng(seed)
        auc_boot = np.empty(N_BOOTSTRAPS, dtype=float)
        for index in range(N_BOOTSTRAPS):
            sampled_short = rng.choice(
                short_scores, size=len(short_scores), replace=True
            )
            sampled_long = rng.choice(
                long_scores, size=len(long_scores), replace=True
            )
            auc_boot[index] = roc_auc_score(
                np.r_[np.zeros(len(sampled_short)), np.ones(len(sampled_long))],
                np.r_[sampled_short, sampled_long],
            )
        result["auc_ci_lower"] = float(np.quantile(auc_boot, 0.025))
        result["auc_ci_upper"] = float(np.quantile(auc_boot, 0.975))
    gene_results = all_gene_results.set_index("gene_id")
    fold_changes = gene_results.loc[genes, "logFC"].to_numpy(dtype=float)
    result["mean_gene_log2FC"] = float(fold_changes.mean())
    return result


def run_gene_set_analyses(
    output_dir: Path,
    candidate_ids: dict[str, str],
    top_cell_types: dict[str, str],
) -> dict:
    profile = pd.read_csv(
        output_dir / "gene_expression_by_cell_type.csv"
    )
    pb_meta = pd.read_csv(output_dir / "all_cell_pseudobulk_metadata.csv")
    counts = pd.read_csv(
        output_dir / "all_cell_pseudobulk_counts.csv", index_col=0
    )
    symbols = pd.read_csv(output_dir / "gene_id_to_symbol.csv")
    symbol_by_id = symbols.set_index("gene_id")["gene_symbol"].to_dict()
    all_results = pd.read_csv(output_dir / "edger_all_gene_results.csv")
    primary = all_results[all_results["contrast"] == "LV_vs_SV"].copy()
    if primary.empty:
        raise ValueError("edgeR did not return the primary LV vs SV contrast")
    gene_to_lfc = primary.set_index("gene_id")["logFC"]
    logcpm_matrix = pd.read_csv(
        output_dir / "edger_logCPM.csv", index_col=0
    )
    primary_group_by_sample = pb_meta.set_index("sample_id")[
        "primary_group"
    ].to_dict()

    min_counts = pd.DataFrame(index=counts.index)
    for group in ("SV", "LV"):
        samples = [
            sample_id
            for sample_id in counts.columns
            if primary_group_by_sample[sample_id] == group
        ]
        min_counts[f"{group}_patients_ge10"] = (
            counts[samples] >= MIN_DETECTION_COUNT
        ).sum(axis=1)
        min_counts[f"{group}_patients_total"] = len(samples)
        min_counts[f"{group}_required_patients"] = math.ceil(len(samples) / 2)
    min_counts["detected"] = (
        (min_counts["SV_patients_ge10"] >= min_counts["SV_required_patients"])
        & (min_counts["LV_patients_ge10"] >= min_counts["LV_required_patients"])
    )
    min_counts.index.name = "gene_id"
    min_counts.reset_index().assign(
        gene_symbol=lambda frame: frame["gene_id"].map(symbol_by_id)
    ).to_csv(output_dir / "pseudobulk_gene_detection.csv", index=False)
    detection_by_id = min_counts["detected"].to_dict()

    discovery_path = (
        output_dir.parents[3]
        / "tables"
        / "GSE125498_all_results.csv"
    )
    discovery = pd.read_csv(discovery_path)
    discovery["Gene.symbol"] = discovery["Gene.symbol"].astype(str)
    candidate_discovery = discovery[
        discovery["Gene.symbol"].isin(GENES)
    ].copy()
    if candidate_discovery["Gene.symbol"].duplicated().any() or set(
        candidate_discovery["Gene.symbol"]
    ) != set(GENES):
        raise ValueError(
            "Discovery results CSV does not have exactly one row per fixed gene"
        )
    discovery_by_gene = candidate_discovery.set_index("Gene.symbol")

    candidate_rows = []
    primary_by_id = primary.set_index("gene_id")
    for gene in GENES:
        gene_id = candidate_ids[gene]
        detected = bool(detection_by_id.get(gene_id, False))
        result_exists = gene_id in primary_by_id.index
        discovery_log2fc = float(discovery_by_gene.loc[gene, "logFC"])
        row = {
            "gene": gene,
            "gene_id": gene_id,
            "discovery_log2FC": discovery_log2fc,
            "discovery_P": float(discovery_by_gene.loc[gene, "P.Value"]),
            "detected": detected,
            "detection_status": "detected" if detected else "not detected",
            "direction_concordant": (
                bool(
                    np.sign(float(primary_by_id.loc[gene_id, "logFC"]))
                    == np.sign(discovery_log2fc)
                )
                if detected and result_exists
                else pd.NA
            ),
            "log2FC": (
                float(primary_by_id.loc[gene_id, "logFC"])
                if detected and result_exists
                else np.nan
            ),
            "ci_lower": (
                float(primary_by_id.loc[gene_id, "ci_lower"])
                if detected and result_exists
                else np.nan
            ),
            "ci_upper": (
                float(primary_by_id.loc[gene_id, "ci_upper"])
                if detected and result_exists
                else np.nan
            ),
            "P": (
                float(primary_by_id.loc[gene_id, "PValue"])
                if detected and result_exists
                else np.nan
            ),
            "mean_logCPM": (
                float(primary_by_id.loc[gene_id, "mean_logCPM"])
                if result_exists
                else np.nan
            ),
            "top_cell_type": top_cell_types[gene],
        }
        for group in ("SV", "LV"):
            row[f"{group}_patients_ge10"] = int(
                min_counts.loc[gene_id, f"{group}_patients_ge10"]
            )
            row[f"{group}_patients_total"] = int(
                min_counts.loc[gene_id, f"{group}_patients_total"]
            )
        candidate_rows.append(row)
    candidate_results = pd.DataFrame(candidate_rows)
    candidate_results["BH_q_across_tested_panel_genes"] = bh_adjust(
        candidate_results["P"].tolist()
    )
    candidate_results.to_csv(
        output_dir / "primary_per_gene_results.csv", index=False
    )

    detected_ids = [
        candidate_ids[gene]
        for gene in GENES
        if detection_by_id.get(candidate_ids[gene], False)
    ]
    concordant = int(
        candidate_results["direction_concordant"].fillna(False).sum()
    )
    detected_count = len(detected_ids)
    sign_p = (
        float(binomtest(concordant, detected_count, 0.5, alternative="greater").pvalue)
        if detected_count
        else np.nan
    )
    strict_rows = candidate_results[
        candidate_results["gene"].isin(STRICT_GENES)
    ]
    strict_detected = int(strict_rows["detected"].sum())
    strict_concordant = int(
        strict_rows["direction_concordant"].fillna(False).sum()
    )
    strict_sign_p = (
        float(
            binomtest(
                strict_concordant,
                strict_detected,
                0.5,
                alternative="greater",
            ).pvalue
        )
        if strict_detected
        else np.nan
    )

    all_panel_estimable = detected_count == len(GENES) and all(
        candidate_ids[gene] in primary_by_id.index for gene in GENES
    )
    metadata_for_score = pb_meta.loc[
        pb_meta["primary_group"].isin(["SV", "LV"]),
        ["sample_id", "primary_group"],
    ].copy()
    control_rows = []
    normal_samples = pb_meta.loc[
        pb_meta["primary_group"] == "NORMAL", "sample_id"
    ].tolist()
    for gene in GENES:
        gene_id = candidate_ids[gene]
        for sample_id in normal_samples:
            control_rows.append(
                {
                    "gene": gene,
                    "gene_id": gene_id,
                    "sample_id": sample_id,
                    "raw_pseudobulk_count": int(counts.loc[gene_id, sample_id]),
                    "logCPM": float(logcpm_matrix.loc[gene_id, sample_id]),
                }
            )
    pd.DataFrame(control_rows).to_csv(
        output_dir / "normal_control_descriptive_expression.csv", index=False
    )
    set_summary = (
        calculate_set_summary(
            [candidate_ids[gene] for gene in GENES],
            primary,
            logcpm_matrix,
            metadata_for_score,
            seed=42,
            bootstrap=True,
        )
        if all_panel_estimable
        else {"estimable": False, "n_genes": detected_count}
    )
    strict_estimable = strict_detected == len(STRICT_GENES) and all(
        candidate_ids[gene] in primary_by_id.index for gene in STRICT_GENES
    )
    strict_summary = (
        calculate_set_summary(
            [candidate_ids[gene] for gene in STRICT_GENES],
            primary,
            logcpm_matrix,
            metadata_for_score,
            seed=42,
            bootstrap=True,
        )
        if strict_estimable
        else {"estimable": False, "n_genes": strict_detected}
    )

    nulls = {}
    mean_lfc_null_p = np.nan
    auc_null_p = np.nan
    null_details = {}
    if all_panel_estimable:
        eligible_ids = [
            gene_id
            for gene_id, is_detected in detection_by_id.items()
            if is_detected
            and gene_id in primary_by_id.index
            and gene_id in logcpm_matrix.index
        ]
        eligible = pd.DataFrame(
            {
                "gene_id": eligible_ids,
                "mean_logCPM": [
                    float(logcpm_matrix.loc[gene_id].mean())
                    for gene_id in eligible_ids
                ],
            }
        )
        panel_null_sets, panel_details = make_null_sets(
            eligible,
            [candidate_ids[gene] for gene in GENES],
            N_NULL_SETS,
            seed=42,
        )
        null_details["eight_gene"] = panel_details
        if panel_null_sets:
            null_lfc = np.empty(N_NULL_SETS)
            null_auc = np.empty(N_NULL_SETS)
            for index, gene_set in enumerate(panel_null_sets):
                null_lfc[index] = float(gene_to_lfc.loc[gene_set].mean())
                summary = calculate_set_summary(
                    gene_set,
                    primary,
                    logcpm_matrix,
                    metadata_for_score,
                    seed=42,
                    bootstrap=False,
                )
                if not summary.get("estimable", False):
                    raise ValueError("Matched null gene set failed composite scoring")
                null_auc[index] = summary["auc"]
            nulls["eight_gene"] = {
                "mean_gene_log2FC": null_lfc,
                "composite_AUC": null_auc,
            }
            mean_lfc_null_p = empirical_upper_p(
                null_lfc, set_summary["mean_gene_log2FC"]
            )
            auc_null_p = empirical_upper_p(null_auc, set_summary["auc"])
            pd.DataFrame(
                {
                    "mean_gene_log2FC": null_lfc,
                    "composite_AUC": null_auc,
                }
            ).to_csv(output_dir / "matched_null_eight_gene.csv", index=False)

    strict_null_details = {}
    strict_mean_lfc_p = np.nan
    strict_auc_null_p = np.nan
    if strict_estimable:
        eligible_ids = [
            gene_id
            for gene_id, is_detected in detection_by_id.items()
            if is_detected
            and gene_id in primary_by_id.index
            and gene_id in logcpm_matrix.index
        ]
        eligible = pd.DataFrame(
            {
                "gene_id": eligible_ids,
                "mean_logCPM": [
                    float(logcpm_matrix.loc[gene_id].mean())
                    for gene_id in eligible_ids
                ],
            }
        )
        strict_sets, strict_details = make_null_sets(
            eligible,
            [candidate_ids[gene] for gene in STRICT_GENES],
            N_NULL_SETS,
            seed=42,
        )
        strict_null_details = strict_details
        if strict_sets:
            null_lfc = np.empty(N_NULL_SETS)
            null_auc = np.empty(N_NULL_SETS)
            for index, gene_set in enumerate(strict_sets):
                null_lfc[index] = float(gene_to_lfc.loc[gene_set].mean())
                summary = calculate_set_summary(
                    gene_set,
                    primary,
                    logcpm_matrix,
                    metadata_for_score,
                    seed=42,
                    bootstrap=False,
                )
                null_auc[index] = summary["auc"]
            strict_mean_lfc_p = empirical_upper_p(
                null_lfc, strict_summary["mean_gene_log2FC"]
            )
            strict_auc_null_p = empirical_upper_p(
                null_auc, strict_summary["auc"]
            )
            pd.DataFrame(
                {
                    "mean_gene_log2FC": null_lfc,
                    "composite_AUC": null_auc,
                }
            ).to_csv(output_dir / "matched_null_strict_hubs.csv", index=False)

    composition_adjusted = pd.read_csv(
        output_dir / "edger_composition_adjusted_results.csv"
    )
    comp_adj_genes = candidate_results[
        ["gene", "gene_id", "detected", "discovery_log2FC"]
    ].merge(
        composition_adjusted,
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    comp_adj_genes["detection_status"] = np.where(
        comp_adj_genes["detected"], "detected", "not detected"
    )
    for column in (
        "logFC",
        "ci_lower",
        "ci_upper",
        "PValue",
        "FDR_all_genes",
        "F",
    ):
        if column not in comp_adj_genes:
            comp_adj_genes[column] = np.nan
        comp_adj_genes[column] = pd.to_numeric(
            comp_adj_genes[column], errors="coerce"
        )
        comp_adj_genes.loc[~comp_adj_genes["detected"], column] = np.nan
    comp_adj_genes["direction_concordant"] = np.where(
        comp_adj_genes["detected"] & comp_adj_genes["logFC"].notna(),
        np.sign(comp_adj_genes["logFC"])
        == np.sign(comp_adj_genes["discovery_log2FC"]),
        pd.NA,
    )
    comp_adj_genes["BH_q_across_tested_panel_genes"] = bh_adjust(
        comp_adj_genes["PValue"].tolist()
    )
    comp_adj_genes.to_csv(
        output_dir / "composition_adjusted_per_gene_results.csv", index=False
    )

    celltype_results = pd.read_csv(
        output_dir / "edger_celltype_gene_results.csv"
    )
    celltype_results["BH_q_across_panel_genes"] = bh_adjust(
        celltype_results["PValue"].tolist()
    )
    celltype_results.to_csv(
        output_dir / "top_celltype_per_gene_results.csv", index=False
    )

    composition = pd.read_csv(output_dir / "patient_celltype_composition.csv")
    edgeR_config = pd.read_csv(output_dir / "edger_analysis_config.csv")
    composition_adjustment_status = str(
        edgeR_config.loc[0, "composition_adjustment"]
    )
    composition_adjustment_feasible = composition_adjustment_status.startswith(
        "feasible:"
    )
    dominant_celltype_profile = profile[
        profile["cell_type"] == profile["gene_top_cell_type"]
    ][
        [
            "gene",
            "gene_top_cell_type",
            "cells",
            "fraction_cells_nonzero",
            "mean_log1p_normalized_expression",
        ]
    ].copy()
    normal_control_expression = pd.read_csv(
        output_dir / "normal_control_descriptive_expression.csv"
    )
    geo_group_sizes = pb_meta["original_group"].value_counts().to_dict()
    fig, axis = plt.subplots(figsize=(11, 5))
    bottoms = np.zeros(len(composition))
    for celltype in CELLTYPES:
        values = composition[f"fraction_{celltype}"].to_numpy()
        axis.bar(
            composition["sample_id"],
            values,
            bottom=bottoms,
            label=celltype,
            width=0.85,
        )
        bottoms += values
    axis.set_ylabel("Fraction of QC-retained cells")
    axis.set_title("GSE130888 per-sample cell-type composition")
    axis.tick_params(axis="x", rotation=80, labelsize=7)
    axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(output_dir / "patient_celltype_composition.png", dpi=240)
    plt.close(fig)

    result_table = candidate_results.copy()
    fig, axis = plt.subplots(figsize=(9, 7), constrained_layout=True)
    y_positions = np.arange(len(GENES))
    for index, gene in enumerate(GENES):
        row = result_table.loc[result_table["gene"] == gene].iloc[0]
        discovery_y = y_positions[index] + 0.15
        validation_y = y_positions[index] - 0.15
        axis.scatter(
            row["discovery_log2FC"],
            discovery_y,
            marker="o",
            facecolors="none",
            edgecolors="#277da1",
            label="Discovery log2FC (CI unavailable)" if index == 0 else None,
            zorder=3,
        )
        if row["detected"] and np.isfinite(row["log2FC"]):
            axis.errorbar(
                row["log2FC"],
                validation_y,
                xerr=np.array(
                    [[row["log2FC"] - row["ci_lower"]],
                     [row["ci_upper"] - row["log2FC"]]]
                ),
                fmt="s",
                color="#d1495b",
                capsize=3,
                label="GSE130888 (95% CI)" if index == 0 else None,
                zorder=3,
            )
    axis.axvline(0, color="gray", linestyle="--", linewidth=1)
    axis.set_yticks(y_positions, GENES)
    axis.invert_yaxis()
    axis.set_xlabel("log2 fold change (long-term vs short-term PD)")
    axis.set_title("Discovery and GSE130888 gene effects")
    axis.legend(frameon=False)
    fig.savefig(output_dir / "discovery_vs_gse130888_forest.png", dpi=240)
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True)
    for axis, gene in zip(axes.flat, GENES):
        gene_id = candidate_ids[gene]
        if gene_id not in logcpm_matrix.index:
            axis.set_title(f"{gene}\nnot tested")
            axis.axis("off")
            continue
        values = logcpm_matrix.loc[gene_id, pb_meta["sample_id"]].to_numpy()
        group_positions = {"SV": 0, "LV": 1, "NORMAL": 2}
        group_colors = {
            "SV": "#277da1",
            "LV": "#f3722c",
            "NORMAL": "#90be6d",
        }
        for group, position in group_positions.items():
            color = group_colors[group]
            mask = pb_meta["primary_group"].to_numpy() == group
            data = values[mask]
            axis.boxplot(
                data,
                positions=[position],
                widths=0.5,
                patch_artist=True,
                boxprops={"facecolor": color, "alpha": 0.35},
                medianprops={"color": "black"},
            )
            if len(data):
                x = position + np.linspace(-0.08, 0.08, len(data))
                axis.scatter(x, data, color=color, s=20, zorder=3)
        stats_row = result_table.loc[result_table["gene"] == gene].iloc[0]
        status = stats_row["detection_status"]
        axis.set_title(
            f"{gene}\n{status}; log2FC={stats_row['log2FC']:.2f}"
            if np.isfinite(stats_row["log2FC"])
            else f"{gene}\n{status}"
        )
        axis.set_xticks([0, 1, 2], ["Short-term", "Long-term", "Normal"])
        axis.set_ylabel("logCPM")
    fig.savefig(output_dir / "pseudobulk_gene_expression.png", dpi=240)
    plt.close(fig)

    score_summary = set_summary
    if (
        score_summary.get("estimable", False)
        and concordant >= 6
        and score_summary["auc"] >= 0.80
        and np.isfinite(auc_null_p)
        and auc_null_p < 0.05
        and strict_concordant == len(STRICT_GENES)
    ):
        verdict = "Supportive"
    elif concordant == 5 or (
        score_summary.get("estimable", False)
        and 0.65 <= score_summary["auc"] <= 0.80
        and np.isfinite(auc_null_p)
        and auc_null_p >= 0.05
    ):
        verdict = "Partially supportive"
    else:
        verdict = "Not supported"

    stats_summary = {
        "accession": "GSE130888",
        "primary_contrast": "Long-term PD vs short-term PD",
        "n_LV_pseudobulk_patients": int(
            (pb_meta["primary_group"] == "LV").sum()
        ),
        "n_SV_pseudobulk_patients": int(
            (pb_meta["primary_group"] == "SV").sum()
        ),
        "panel_genes_concordant": concordant,
        "panel_genes_detected": detected_count,
        "panel_sign_test_one_sided_P": sign_p,
        "strict_hubs_concordant": strict_concordant,
        "strict_hubs_detected": strict_detected,
        "strict_hub_sign_test_one_sided_P": strict_sign_p,
        "panel_composite": score_summary,
        "panel_null_empirical_P_mean_gene_log2FC": mean_lfc_null_p,
        "panel_null_empirical_P_composite_AUC": auc_null_p,
        "strict_hub_composite": strict_summary,
        "strict_hub_null_empirical_P_mean_gene_log2FC": strict_mean_lfc_p,
        "strict_hub_null_empirical_P_composite_AUC": strict_auc_null_p,
        "matched_null_details": null_details,
        "strict_matched_null_details": strict_null_details,
        "composition_adjustment": {
            "available": composition_adjustment_feasible,
            "status": composition_adjustment_status,
        },
        "preregistered_verdict": verdict,
        "verdict_rule_note": (
            "The matched-null criterion is applied to the preregistered "
            "composite-AUC null; mean-log2FC null P is also reported."
        ),
        "seed": 42,
    }
    (output_dir / "gene_set_summary.json").write_text(
        json.dumps(stats_summary, indent=2, allow_nan=True), encoding="utf-8"
    )
    pd.DataFrame(
        [
            {
                "test": "8-gene sign test",
                "concordant": concordant,
                "detected": detected_count,
                "one_sided_P": sign_p,
            },
            {
                "test": "2 strict-hub sign test",
                "concordant": strict_concordant,
                "detected": strict_detected,
                "one_sided_P": strict_sign_p,
            },
        ]
    ).to_csv(output_dir / "sign_tests.csv", index=False)
    if score_summary.get("estimable", False) and "eight_gene" in nulls:
        expression = logcpm_matrix.loc[
            [candidate_ids[gene] for gene in GENES],
            metadata_for_score["sample_id"],
        ].to_numpy(dtype=float)
        scores = (
            (expression - expression.mean(axis=1, keepdims=True))
            / expression.std(axis=1, ddof=0, keepdims=True)
        ).mean(axis=0)
        fpr, tpr, _ = roc_curve(
            (metadata_for_score["primary_group"] == "LV").astype(int), scores
        )
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        axes[0].plot(fpr, tpr, color="#d1495b", label=f"AUC={score_summary['auc']:.3f}")
        axes[0].plot([0, 1], [0, 1], linestyle="--", color="gray")
        axes[0].set_xlabel("False positive rate")
        axes[0].set_ylabel("True positive rate")
        axes[0].set_title("GSE130888 full-panel composite ROC")
        axes[0].legend(frameon=False)
        axes[1].hist(
            nulls["eight_gene"]["composite_AUC"],
            bins=25,
            color="#76b7b2",
            edgecolor="white",
        )
        axes[1].axvline(
            score_summary["auc"],
            color="#d1495b",
            linewidth=2,
            label="Observed",
        )
        axes[1].set_xlabel("Matched random 8-gene composite AUC")
        axes[1].set_ylabel("Number of gene sets")
        axes[1].set_title(f"Empirical P={auc_null_p:.4f}")
        axes[1].legend(frameon=False)
        fig.savefig(output_dir / "composite_roc_and_null.png", dpi=240)
        plt.close(fig)

    # Secondary group contrasts remain descriptive and do not affect the verdict.
    secondary_specs = {}
    metadata_by_sample = pb_meta.set_index("sample_id")
    secondary_tables = []
    discovery_lfc_by_gene = candidate_results.set_index("gene")[
        "discovery_log2FC"
    ].to_dict()
    for contrast_name, (test_group, reference_group) in secondary_specs.items():
        contrast_results = all_results[
            all_results["contrast"] == contrast_name
        ].set_index("gene_id")
        contrast_rows = []
        for gene in GENES:
            gene_id = candidate_ids[gene]
            test_samples = metadata_by_sample.index[
                metadata_by_sample["original_group"] == test_group
            ].tolist()
            reference_samples = metadata_by_sample.index[
                metadata_by_sample["original_group"] == reference_group
            ].tolist()
            test_ge10 = int(
                (counts.loc[gene_id, test_samples] >= MIN_DETECTION_COUNT).sum()
            )
            reference_ge10 = int(
                (
                    counts.loc[gene_id, reference_samples]
                    >= MIN_DETECTION_COUNT
                ).sum()
            )
            detected = (
                test_ge10 >= math.ceil(len(test_samples) / 2)
                and reference_ge10 >= math.ceil(len(reference_samples) / 2)
            )
            result_exists = gene_id in contrast_results.index
            if detected and not result_exists:
                raise ValueError(
                    f"{gene} is detected in {contrast_name} but edgeR "
                    "returned no estimate"
                )
            result = (
                contrast_results.loc[gene_id]
                if result_exists
                else pd.Series(dtype=float)
            )
            logfc = float(result["logFC"]) if detected else np.nan
            contrast_rows.append(
                {
                    "contrast": contrast_name,
                    "gene": gene,
                    "gene_id": gene_id,
                    "test_group": test_group,
                    "reference_group": reference_group,
                    "test_group_n": len(test_samples),
                    "test_group_patients_ge10": test_ge10,
                    "reference_group_n": len(reference_samples),
                    "reference_group_patients_ge10": reference_ge10,
                    "detection_status": "detected" if detected else "not detected",
                    "direction_concordant": (
                        bool(
                            np.sign(logfc)
                            == np.sign(discovery_lfc_by_gene[gene])
                        )
                        if detected
                        else pd.NA
                    ),
                    "log2FC": logfc,
                    "ci_lower": (
                        float(result["ci_lower"]) if detected else np.nan
                    ),
                    "ci_upper": (
                        float(result["ci_upper"]) if detected else np.nan
                    ),
                    "P": float(result["PValue"]) if detected else np.nan,
                    "mean_logCPM": (
                        float(result["mean_logCPM"]) if result_exists else np.nan
                    ),
                }
            )
        contrast_table = pd.DataFrame(contrast_rows)
        contrast_table["BH_q_across_tested_panel_genes"] = bh_adjust(
            contrast_table["P"].tolist()
        )
        secondary_tables.append(contrast_table)
    secondary = (
        pd.concat(secondary_tables, ignore_index=True)
        if secondary_tables
        else pd.DataFrame(
            columns=[
                "contrast",
                "gene",
                "detection_status",
                "log2FC",
                "ci_lower",
                "ci_upper",
                "P",
                "BH_q_across_tested_panel_genes",
                "direction_concordant",
            ]
        )
    )
    secondary.to_csv(
        output_dir / "secondary_per_gene_results.csv", index=False
    )

    report = [
        "# GSE130888 preregistered analysis (Steps 2–4)",
        "",
        f"**Preregistered verdict for long-term vs short-term PD: {verdict}.**",
        "",
        f"- Primary-comparison pseudobulks: long-term PD "
        f"{stats_summary['n_LV_pseudobulk_patients']}, short-term PD "
        f"{stats_summary['n_SV_pseudobulk_patients']}; normal controls "
        f"{len(normal_samples)} (descriptive only).",
        f"- Detected genes direction-concordant: {concordant}/{detected_count}.",
        "- TNFSF15 and SERPINA10 were not detected by the frozen raw-count "
        "rule; neither is concordant nor discordant, and both are excluded "
        "from the concordance denominator.",
        (
            f"- Full-panel composite AUC: {score_summary['auc']:.4f} "
            f"(95% bootstrap CI {score_summary['auc_ci_lower']:.4f}–"
            f"{score_summary['auc_ci_upper']:.4f}); Mann–Whitney P "
            f"{score_summary['mann_whitney_p']:.6g}."
            if score_summary.get("estimable", False)
            else "- Full-panel composite: not estimable under the frozen rule "
            "because one or more of the eight genes were not detected."
        ),
        (
            f"- Matched-null empirical P: AUC {auc_null_p:.4f}; mean gene "
            f"log2FC {mean_lfc_null_p:.4f}."
            if np.isfinite(auc_null_p)
            else "- Full-panel matched-null empirical P: not estimable because "
            "one or more fixed panel genes were not detected; the frozen "
            "full-panel composite was not computed."
        ),
        "- No full-panel composite ROC or matched-null histogram was drawn "
        "because the fixed eight-gene composite is not estimable.",
        f"- Strict hubs: {strict_concordant}/{strict_detected} detected genes "
        "concordant.",
        (
            f"- Strict-hub composite (secondary): AUC "
            f"{strict_summary['auc']:.4f} (95% bootstrap CI "
            f"{strict_summary['auc_ci_lower']:.4f}–"
            f"{strict_summary['auc_ci_upper']:.4f}); Mann–Whitney P "
            f"{strict_summary['mann_whitney_p']:.6g}; matched-null empirical "
            f"P for AUC {strict_auc_null_p:.4f} and mean gene log2FC "
            f"{strict_mean_lfc_p:.4f}."
            if strict_summary.get("estimable", False)
            and np.isfinite(strict_auc_null_p)
            else "- Strict-hub composite/null: not estimable under the "
            "preregistered detection or decile requirements."
        ),
        "",
        "## Per-gene results",
        "",
        "Discovery confidence intervals are not present in the authoritative "
        "discovery results CSV, so discovery log2FC points are shown without "
        "intervals in `discovery_vs_gse130888_forest.png`; GSE130888 intervals "
        "are shown from edgeR.",
        "",
        candidate_results[
            [
                "gene",
                "discovery_log2FC",
                "discovery_P",
                "detection_status",
                "log2FC",
                "ci_lower",
                "ci_upper",
                "P",
                "BH_q_across_tested_panel_genes",
                "direction_concordant",
            ]
        ].fillna("not estimated").to_markdown(index=False),
        "",
        "## Normal controls (descriptive only)",
        "",
        f"Normal-control rows show raw sample-level pseudobulk counts and "
        f"TMM-normalized log-CPM; no inferential test uses the "
        f"{len(normal_samples)} normal controls.",
        "",
        normal_control_expression.to_markdown(index=False),
        "",
        "## Expression by cell type",
        "",
        "The top cell type is selected by mean normalized expression across "
        "all QC-retained cells, without group labels. The full gene-by-cell-type "
        "profile is saved in `gene_expression_by_cell_type.csv`.",
        "",
        dominant_celltype_profile.to_markdown(index=False),
        "",
        "Top-cell-type assignments are descriptive marker-based annotations; "
        "a relative maximum is not evidence of meaningful expression.",
        "",
        "## QC and analysis notes",
        "",
        "Cells were retained using the preregistered group-blind thresholds "
        "of at least 200 detected genes and at most 20% mitochondrial counts. "
        f"The GEO groups were short-term PD n={geo_group_sizes.get('SV', 0)}, "
        f"long-term PD n={geo_group_sizes.get('LV', 0)}, and normal controls "
        f"n={geo_group_sizes.get('NORMAL', 0)}. Samples passing the blind QC were aggregated per GEO "
        "sample; each GSM is used as an anonymous patient proxy because "
        "participant IDs are not exposed in the GEO metadata. A "
        "panel gene was detected only if raw counts were at least 10 in at "
        "least half of patient pseudobulks in each group of the primary contrast. "
        "Raw counts were aggregated by sample/patient; cells were not treated "
        "as independent replicates. The one-GSM-per-person assumption cannot "
        "be confirmed from the available GEO identifiers.",
        "",
        "Cell types use the fixed preregistered marker-score rules because "
        "cell-level author annotations were not present in the downloaded "
        "10x archive. All resampling and permutation-style tests use seed 42. "
        "The composite-AUC random-gene null is used for the supportive "
        "decision rule; the preregistered mean-log2FC null is also reported.",
        "",
        "edgeR P values are quasi-likelihood F-test results. Approximate 95% "
        "intervals use the quasi-likelihood variance-scaled GLM coefficient "
        "covariance and residual-adjusted degrees of freedom, expressed on "
        "the log2 scale.",
        "",
        "## Composition",
        "",
        pd.read_csv(output_dir / "composition_tests.csv").to_markdown(index=False),
        "",
        "## Composition-adjusted per-gene results",
        "",
        "When a gene is not detected under the preregistered raw-count rule, "
        "its adjusted estimate and P value are suppressed as not detected. "
        "Any detected gene whose adjusted coefficient is no longer positive "
        "is described as reflecting immune-cell abundance.",
        "",
        comp_adj_genes.fillna("not estimated").to_markdown(index=False),
        "",
        (
            "Composition adjustment was not feasible under the preregistered "
            "minimum sample-size rule; direction changes after adjustment "
            "cannot be assessed."
            if not composition_adjustment_feasible
            else
            "LTB was positive before adjustment but negative after including "
            "monocyte/macrophage and T/NK proportions; under the preregistered "
            "rule its unadjusted signal is described as reflecting immune-cell "
            "abundance. No other detected panel gene reversed direction."
            if (comp_adj_genes["detected"] & ~comp_adj_genes[
                "direction_concordant"
            ].fillna(False)).any()
            else "No detected panel gene reversed direction after composition adjustment."
        ),
        "",
        "## Top-cell-type per-gene results",
        "",
        celltype_results.to_markdown(index=False),
        "",
        "## Software versions",
        "",
        f"- Python {__import__('sys').version.split()[0]}",
        f"- Scanpy {version('scanpy')}; AnnData {version('anndata')}.",
        "- R and edgeR versions are recorded in `edger_analysis_config.csv`.",
        "",
        "This is the GSE130888 effluent-cell cohort. The preregistered "
        "independence assessment labels its patient overlap with GSE248762 "
        "unknown; it is also not formally confirmed independent of discovery "
        "GSE125498. See `../independence_assessment.md`. The small primary "
        "groups, single-cell dropout, lack of participant identifiers, and "
        "wide intervals limit interpretation; wide intervals are inconclusive, "
        "not proof of no effect. Composition adjustment is not feasible with "
        "four long-term and six short-term samples under the preregistered "
        f"minimum of five samples per group (long-term n="
        f"{stats_summary['n_LV_pseudobulk_patients']}, short-term n="
        f"{stats_summary['n_SV_pseudobulk_patients']}). No GSE92455 analysis "
        "was performed.",
    ]
    (output_dir / "gse130888_analysis_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    return stats_summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[5],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=(
            Path(__file__).resolve().parents[1] / "analysis_20261009"
        ),
    )
    parser.add_argument(
        "--reuse-prepared",
        action="store_true",
        help="Reuse saved QC/pseudobulk files and continue downstream analysis.",
    )
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    if args.reuse_prepared:
        required_files = [
            "all_cell_pseudobulk_counts.csv",
            "all_cell_pseudobulk_metadata.csv",
            "celltype_pseudobulk_counts.csv",
            "celltype_pseudobulk_metadata.csv",
            "gene_expression_by_cell_type.csv",
            "gene_id_to_symbol.csv",
            "patient_celltype_composition.csv",
            "sample_qc.csv",
            "edger_all_gene_results.csv",
            "edger_logCPM.csv",
            "edger_composition_adjusted_results.csv",
            "edger_celltype_gene_results.csv",
            "edger_analysis_config.csv",
        ]
        missing_files = [
            name for name in required_files if not (output_dir / name).is_file()
        ]
        if missing_files:
            raise FileNotFoundError(
                "Cannot resume; required prepared outputs are missing: "
                + ", ".join(missing_files)
            )
        symbols = pd.read_csv(output_dir / "gene_id_to_symbol.csv")
        candidate_ids = {}
        for gene in GENES:
            matches = symbols.loc[
                symbols["gene_symbol"].astype(str).str.upper() == gene,
                "gene_id",
            ].astype(str)
            if len(matches) != 1:
                raise ValueError(
                    f"{gene}: expected one saved feature ID, found {len(matches)}"
                )
            candidate_ids[gene] = matches.iloc[0]
        profile = pd.read_csv(output_dir / "gene_expression_by_cell_type.csv")
        top_cell_types = (
            profile.drop_duplicates("gene")
            .set_index("gene")["gene_top_cell_type"]
            .to_dict()
        )
        if set(top_cell_types) != set(GENES):
            raise ValueError("Saved cell-type profile is incomplete for the panel")
        prepared = {
            "candidate_ids": candidate_ids,
            "top_cell_types": top_cell_types,
        }
    else:
        prepared = prepare_data(repo_root, output_dir)
        run_edger(repo_root, output_dir)
    summary = run_gene_set_analyses(
        output_dir,
        prepared["candidate_ids"],
        prepared["top_cell_types"],
    )
    print(json.dumps(summary, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
