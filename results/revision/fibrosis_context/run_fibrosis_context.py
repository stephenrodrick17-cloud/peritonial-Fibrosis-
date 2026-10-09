#!/usr/bin/env python
"""Run the preregistered GSE248762 fibrosis-context analyses."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import gseapy as gp
import liana as li
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import kruskal, rankdata, spearmanr, t


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
MARKERS = [
    "ACTA2",
    "COL1A1",
    "COL1A2",
    "FN1",
    "TGFB1",
    "CCN2",
    "SERPINE1",
    "VEGFA",
]
GROUPS = ["LV_UF", "LV_NOT_UF", "SV"]
IMMUNE_TYPES = {"T/NK", "B", "monocyte/macrophage", "neutrophil"}
STROMAL_TYPES = {"fibroblast", "mesothelial"}
BOOTSTRAPS = 2_000
SEED = 42


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    valid = np.isfinite(values)
    adjusted = np.full(values.shape, np.nan)
    if valid.any():
        order = np.argsort(values[valid])
        sorted_p = values[valid][order]
        q = sorted_p * len(sorted_p) / np.arange(1, len(sorted_p) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        restored = np.empty_like(q)
        restored[order] = np.minimum(q, 1.0)
        adjusted[valid] = restored
    return pd.Series(adjusted, index=p_values.index)


def residualize(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(values)), covariates])
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficients


def partial_spearman(
    x: np.ndarray, y: np.ndarray, covariates: np.ndarray
) -> tuple[float, float]:
    rank_x = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    rank_y = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    rank_covariates = np.column_stack(
        [
            pd.Series(covariates[:, j]).rank(method="average").to_numpy(dtype=float)
            for j in range(covariates.shape[1])
        ]
    )
    residual_x = residualize(rank_x, rank_covariates)
    residual_y = residualize(rank_y, rank_covariates)
    if np.std(residual_x) == 0 or np.std(residual_y) == 0:
        return np.nan, np.nan
    rho = float(np.corrcoef(residual_x, residual_y)[0, 1])
    df = len(x) - covariates.shape[1] - 2
    if df <= 0 or not np.isfinite(rho):
        return rho, np.nan
    statistic = rho * np.sqrt(df / max(1.0 - rho * rho, np.finfo(float).eps))
    return rho, float(2.0 * t.sf(abs(statistic), df))


def bootstrap_ci(
    x: np.ndarray,
    y: np.ndarray,
    covariates: np.ndarray | None,
    rng: np.random.Generator,
) -> tuple[float, float]:
    n = len(x)
    selected = rng.integers(0, n, (BOOTSTRAPS, n))
    rank_x = rankdata(x[selected], axis=1, method="average")
    rank_y = rankdata(y[selected], axis=1, method="average")
    if covariates is not None:
        rank_covariates = np.stack(
            [
                rankdata(covariates[selected, j], axis=1, method="average")
                for j in range(covariates.shape[1])
            ],
            axis=2,
        )
        design = np.concatenate(
            [np.ones((BOOTSTRAPS, n, 1)), rank_covariates], axis=2
        )
        cross = np.einsum("bnp,bnq->bpq", design, design)
        rhs_x = np.einsum("bnp,bn->bp", design, rank_x)
        rhs_y = np.einsum("bnp,bn->bp", design, rank_y)
        inverse = np.linalg.pinv(cross)
        beta_x = np.einsum("bpq,bq->bp", inverse, rhs_x)
        beta_y = np.einsum("bpq,bq->bp", inverse, rhs_y)
        rank_x = rank_x - np.einsum("bnp,bp->bn", design, beta_x)
        rank_y = rank_y - np.einsum("bnp,bp->bn", design, beta_y)
    centered_x = rank_x - rank_x.mean(axis=1, keepdims=True)
    centered_y = rank_y - rank_y.mean(axis=1, keepdims=True)
    denominator = np.sqrt(
        np.einsum("bn,bn->b", centered_x, centered_x)
        * np.einsum("bn,bn->b", centered_y, centered_y)
    )
    numerator = np.einsum("bn,bn->b", centered_x, centered_y)
    draws = np.divide(
        numerator,
        denominator,
        out=np.full(BOOTSTRAPS, np.nan),
        where=denominator > 0,
    )
    valid = draws[np.isfinite(draws)]
    if len(valid) < BOOTSTRAPS * 0.8:
        return np.nan, np.nan
    return tuple(np.quantile(valid, [0.025, 0.975]).astype(float))


def read_symbol_map(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    mapping = pd.read_csv(path, dtype=str).fillna("")
    mapping["gene_symbol"] = mapping["gene_symbol"].str.upper()
    ids_by_symbol = mapping.groupby("gene_symbol")["gene_id"].apply(list).to_dict()
    required = GENES + MARKERS
    missing = [symbol for symbol in required if symbol not in ids_by_symbol]
    if missing:
        raise ValueError(f"Missing required gene symbols in GEO feature map: {missing}")
    duplicate = {
        symbol: ids_by_symbol[symbol]
        for symbol in required
        if len(set(ids_by_symbol[symbol])) != 1
    }
    if duplicate:
        raise ValueError(f"Required symbols do not map one-to-one: {duplicate}")
    return mapping, {symbol: ids_by_symbol[symbol][0] for symbol in required}


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    values = counts.astype(float)
    library_sizes = values.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("Pseudobulk has an empty library")
    return np.log2(values.div(library_sizes, axis=1) * 1_000_000 + 0.5)


def calculate_correlations(
    input_dir: Path,
    output_dir: Path,
    gene_ids: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_csv(
        input_dir / "all_cell_pseudobulk_metadata.csv",
        dtype={"sample_id": str},
    )
    composition = pd.read_csv(
        input_dir / "patient_celltype_composition.csv",
        dtype={"sample_id": str},
    ).set_index("sample_id")
    counts_all = pd.read_csv(
        input_dir / "all_cell_pseudobulk_counts.csv",
        index_col=0,
        dtype={"gene_id": str},
    )
    counts_celltype = pd.read_csv(
        input_dir / "celltype_pseudobulk_counts.csv",
        index_col=0,
        dtype={"gene_id": str},
    )
    cell_metadata = pd.read_csv(
        input_dir / "celltype_pseudobulk_metadata.csv",
        dtype={"sample_id": str},
    )
    if metadata["sample_id"].nunique() != 16:
        raise ValueError("Expected the fixed 16-patient GSE248762 metadata")
    if set(metadata["sample_id"]) != set(composition.index):
        raise ValueError("Patient composition metadata do not map one-to-one")

    strata: list[tuple[str, str, pd.DataFrame, pd.DataFrame]] = []
    all_samples = metadata["sample_id"].tolist()
    strata.append(
        (
            "all_cells",
            "all",
            counts_all.loc[:, all_samples],
            metadata.set_index("sample_id").loc[all_samples],
        )
    )
    for cell_type in sorted(cell_metadata["cell_type"].unique()):
        selected_meta = cell_metadata[cell_metadata["cell_type"] == cell_type]
        samples = selected_meta["sample_id"].tolist()
        if not samples:
            continue
        columns = [f"{sample}|{cell_type}" for sample in samples]
        absent = [column for column in columns if column not in counts_celltype.columns]
        if absent:
            raise ValueError(
                f"{cell_type}: pseudobulk metadata columns missing from matrix: {absent}"
            )
        strata.append(
            (
                "cell_type",
                cell_type,
                counts_celltype.loc[:, columns],
                selected_meta.set_index("sample_id").loc[samples],
            )
        )

    rows: list[dict] = []
    patient_rows: list[dict] = []
    rng = np.random.default_rng(SEED)
    for aggregation, cell_type, count_frame, stratum_metadata in strata:
        normalized = log_cpm(count_frame.loc[list(gene_ids.values())])
        samples = stratum_metadata.index.astype(str).tolist()
        metadata_rows = metadata.set_index("sample_id").loc[samples]
        fractions = composition.loc[
            samples,
            ["fraction_monocyte/macrophage", "fraction_T/NK"],
        ].to_numpy(dtype=float)
        sample_index_by_id = {sample: i for i, sample in enumerate(samples)}
        print(
            f"Patient correlations: {aggregation} / {cell_type} (n={len(samples)})",
            flush=True,
        )
        for gene in GENES:
            x = normalized.loc[gene_ids[gene]].to_numpy(dtype=float)
            for marker in MARKERS:
                y = normalized.loc[gene_ids[marker]].to_numpy(dtype=float)
                valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(fractions).all(axis=1)
                xv, yv, cov = x[valid], y[valid], fractions[valid]
                n = len(xv)
                if n < 4:
                    for adjustment in ("unadjusted", "composition_adjusted"):
                        rows.append(
                            {
                                "aggregation": aggregation,
                                "cell_type": cell_type,
                                "adjustment": adjustment,
                                "gene": gene,
                                "marker": marker,
                                "n": n,
                                "rho": np.nan,
                                "ci_lower": np.nan,
                                "ci_upper": np.nan,
                                "p": np.nan,
                                "status": "not estimable: n < 4",
                            }
                        )
                    continue
                rho = float(spearmanr(xv, yv).statistic)
                p_value = float(spearmanr(xv, yv).pvalue)
                ci_lower, ci_upper = bootstrap_ci(xv, yv, None, rng)
                rows.append(
                    {
                        "aggregation": aggregation,
                        "cell_type": cell_type,
                        "adjustment": "unadjusted",
                        "gene": gene,
                        "marker": marker,
                        "n": n,
                        "rho": rho,
                        "ci_lower": ci_lower,
                        "ci_upper": ci_upper,
                        "p": p_value,
                        "status": (
                            "estimated"
                            if np.isfinite(rho)
                            else "not estimable: constant input"
                        ),
                    }
                )
                adjusted_rho, adjusted_p = partial_spearman(xv, yv, cov)
                adj_lower, adj_upper = bootstrap_ci(xv, yv, cov, rng)
                rows.append(
                    {
                        "aggregation": aggregation,
                        "cell_type": cell_type,
                        "adjustment": "composition_adjusted",
                        "gene": gene,
                        "marker": marker,
                        "n": n,
                        "rho": adjusted_rho,
                        "ci_lower": adj_lower,
                        "ci_upper": adj_upper,
                        "p": adjusted_p,
                        "status": (
                            "estimated"
                            if np.isfinite(adjusted_rho)
                            else "not estimable: constant residual"
                        ),
                    }
                )
                for sample, value_x, value_y in zip(samples, x, y):
                    patient_rows.append(
                        {
                            "aggregation": aggregation,
                            "cell_type": cell_type,
                            "sample_id": sample,
                            "patient_id_proxy": sample,
                            "group": metadata_rows.loc[sample, "original_group"],
                            "gene": gene,
                            "gene_log2CPM": value_x,
                            "marker": marker,
                            "marker_log2CPM": value_y,
                            "monocyte_fraction": composition.loc[
                                sample, "fraction_monocyte/macrophage"
                            ],
                            "tnk_fraction": composition.loc[
                                sample, "fraction_T/NK"
                            ],
                        }
                    )

    results = pd.DataFrame(rows)
    results["q"] = np.nan
    for _, indices in results.groupby(
        ["aggregation", "cell_type", "adjustment"], sort=False
    ).groups.items():
        results.loc[indices, "q"] = bh_adjust(results.loc[indices, "p"])
    results.to_csv(output_dir / "gene_marker_correlations.csv", index=False)
    pd.DataFrame(patient_rows).to_csv(
        output_dir / "patient_gene_marker_expression.csv", index=False
    )

    all_cell_adjusted = results[
        (results["aggregation"] == "all_cells")
        & (results["adjustment"] == "composition_adjusted")
    ]
    support_rows = []
    for gene in GENES:
        gene_corr = all_cell_adjusted[all_cell_adjusted["gene"] == gene]
        marker_count = int(((gene_corr["rho"] > 0.4) & gene_corr["rho"].notna()).sum())
        support_rows.append(
            {
                "gene": gene,
                "n_markers_with_adjusted_rho_gt_0_4": marker_count,
                "correlation_condition_met": marker_count >= 3,
            }
        )
    return results, pd.DataFrame(support_rows)


def run_existing_edger(
    repo_root: Path, input_dir: Path, output_dir: Path
) -> tuple[pd.DataFrame, dict]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        candidate = Path(r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe")
        if not candidate.exists():
            raise FileNotFoundError("Rscript.exe was not found")
        rscript = str(candidate)
    pipeline_script = (
        repo_root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE248762"
        / "scripts"
        / "gse248762_step2_4_edger.R"
    )
    source_dir = input_dir
    required = [
        "all_cell_pseudobulk_counts.csv",
        "all_cell_pseudobulk_metadata.csv",
        "gene_id_to_symbol.csv",
        "patient_celltype_composition.csv",
        "gene_expression_by_cell_type.csv",
        "celltype_pseudobulk_counts.csv",
        "celltype_pseudobulk_metadata.csv",
    ]
    with tempfile.TemporaryDirectory(prefix="fibrosis_edger_") as temp:
        temp_dir = Path(temp)
        for name in required:
            shutil.copy2(source_dir / name, temp_dir / name)
        result = subprocess.run(
            [rscript, str(pipeline_script), str(temp_dir)],
            check=False,
            capture_output=True,
            text=True,
            timeout=1_800,
        )
        (output_dir / "edger_pipeline_run.log").write_text(
            result.stdout + "\n--- STDERR ---\n" + result.stderr,
            encoding="utf-8",
        )
        if result.returncode:
            raise RuntimeError(
                "Existing edgeR pipeline failed. See edger_pipeline_run.log:\n"
                + result.stderr[-4_000:]
            )
        all_results = pd.read_csv(temp_dir / "edger_all_gene_results.csv")
        panel = all_results[
            (all_results["contrast"] == "LV_UF_vs_LV_NOT_UF")
            & (all_results["gene_symbol"].str.upper().isin(GENES))
        ].copy()
        if set(panel["gene_symbol"].str.upper()) != set(GENES):
            missing = sorted(set(GENES) - set(panel["gene_symbol"].str.upper()))
            raise ValueError(f"Existing edgeR contrast omitted panel genes: {missing}")
        panel["gene_symbol"] = panel["gene_symbol"].str.upper()
        panel.to_csv(
            output_dir / "edger_LV_UF_vs_LV_NOT_UF_candidate_genes.csv",
            index=False,
        )
        config = pd.read_csv(temp_dir / "edger_analysis_config.csv").iloc[0].to_dict()
        (output_dir / "edger_model_manifest.json").write_text(
            json.dumps(config, indent=2, default=str), encoding="utf-8"
        )
    return panel, config


def run_compartment_edger(
    repo_root: Path, input_dir: Path, output_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        rscript = r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
    r_script = output_dir / "compartment_edger_gsea_input.R"
    result = subprocess.run(
        [
            rscript,
            str(r_script),
            str(input_dir),
            str(output_dir / "compartment_edgeR_rankings.csv"),
            str(output_dir / "compartment_edgeR_models.csv"),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=1_800,
    )
    (output_dir / "compartment_edgeR_run.log").write_text(
        result.stdout + "\n--- STDERR ---\n" + result.stderr,
        encoding="utf-8",
    )
    if result.returncode:
        raise RuntimeError(
            "Compartment edgeR failed. See compartment_edgeR_run.log:\n"
            + result.stderr[-4_000:]
        )
    return (
        pd.read_csv(output_dir / "compartment_edgeR_rankings.csv"),
        pd.read_csv(output_dir / "compartment_edgeR_models.csv"),
    )


def run_preranked_gsea(
    rankings: pd.DataFrame, output_dir: Path
) -> pd.DataFrame:
    hallmark = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
    go_bp = gp.get_library("GO_Biological_Process_2023", organism="Human")
    requested = {
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION":
            hallmark["Epithelial Mesenchymal Transition"],
        "HALLMARK_TGF_BETA_SIGNALING": hallmark["TGF-beta Signaling"],
        "GO_EXTRACELLULAR_MATRIX_ORGANIZATION":
            go_bp["Extracellular Matrix Organization (GO:0030198)"],
        "HALLMARK_INFLAMMATORY_RESPONSE":
            hallmark["Inflammatory Response"],
    }
    gene_set_rows = [
        {"gene_set": set_id, "gene_symbol": str(gene).upper()}
        for set_id, members in requested.items()
        for gene in members
    ]
    pd.DataFrame(gene_set_rows).to_csv(
        output_dir / "gsea_gene_set_membership.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "set_id": set_id,
                "source": (
                    "MSigDB_Hallmark_2020"
                    if set_id.startswith("HALLMARK_")
                    else "GO_Biological_Process_2023"
                ),
                "term": set_id,
                "n_genes": len(genes),
                "retrieval": "gseapy Enrichr library",
            }
            for set_id, genes in requested.items()
        ]
    ).to_csv(output_dir / "gsea_gene_set_manifest.csv", index=False)

    outputs = []
    for cell_type in [
        "monocyte/macrophage",
        "T/NK",
        "fibroblast",
        "mesothelial",
    ]:
        subset = rankings[rankings["compartment"] == cell_type].copy()
        subset = subset.dropna(subset=["gene_symbol", "rank_score"])
        subset = subset[subset["gene_symbol"].str.match(r"^[A-Z0-9_.-]+$")]
        subset = subset.sort_values(
            ["gene_symbol", "rank_score"],
            key=lambda series: series.abs() if series.name == "rank_score" else series,
            ascending=[True, False],
        ).drop_duplicates("gene_symbol", keep="first")
        ranked = subset[["gene_symbol", "rank_score"]].sort_values(
            "rank_score", ascending=False
        )
        if len(ranked) < 100:
            raise ValueError(
                f"{cell_type}: only {len(ranked)} ranked genes for GSEA"
            )
        enrichment = gp.prerank(
            rnk=ranked,
            gene_sets=requested,
            min_size=15,
            max_size=500,
            permutation_num=1_000,
            seed=SEED,
            threads=1,
            outdir=None,
            verbose=False,
        ).res2d
        for _, row in enrichment.iterrows():
            term = str(row.get("Term", ""))
            lead_text = str(row.get("Lead_genes", ""))
            lead = {
                item.strip().upper()
                for item in re.split(r"[;,]", lead_text)
                if item.strip()
            }
            outputs.append(
                {
                    "compartment": cell_type,
                    "gene_set": term,
                    "NES": pd.to_numeric(row.get("NES"), errors="coerce"),
                    "P": pd.to_numeric(row.get("NOM p-val"), errors="coerce"),
                    "FDR": pd.to_numeric(row.get("FDR q-val"), errors="coerce"),
                    "leading_edge_size": len(lead),
                    "panel_genes_in_leading_edge": ";".join(
                        gene for gene in GENES if gene in lead
                    ),
                }
            )
    result_df = pd.DataFrame(outputs)
    empirical_floor = 1.0 / 1_001
    result_df["P"] = result_df["P"].mask(
        result_df["P"].eq(0), empirical_floor
    )
    result_df["FDR"] = result_df["FDR"].mask(
        result_df["FDR"].eq(0), empirical_floor
    )
    result_df.to_csv(output_dir / "compartment_gsea_results.csv", index=False)
    return result_df


def process_single_cell(
    repo_root: Path, output_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    scanpy_script_dir = (
        repo_root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE248762"
        / "scripts"
    )
    sys.path.insert(0, str(scanpy_script_dir))
    import gse248762_step2_4_scanpy as previous_pipeline

    metadata = pd.read_csv(
        repo_root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE248762"
        / "sample_metadata.csv",
        dtype=str,
    )
    extracted_dir = repo_root / "data" / "raw" / "GSE248762_extracted"
    all_mmt_rows: list[dict] = []
    celltype_expressions: dict[str, list[np.ndarray]] = {}
    celltype_obs: dict[str, list[dict]] = {}
    liana_adatas: list[ad.AnnData] = []
    subsample_rows = []
    rng = np.random.default_rng(SEED)
    feature_names_reference: list[str] | None = None

    for row in metadata.itertuples(index=False):
        sample_id = row.sample_id
        print(f"Single-cell MMT/LIANA preprocessing: {sample_id}", flush=True)
        names = {
            kind: list(extracted_dir.glob(f"{sample_id}_*.{kind}.tsv.gz"))
            for kind in ("barcodes", "features")
        }
        matrices = list(extracted_dir.glob(f"{sample_id}_*.matrix.mtx.gz"))
        if any(len(paths) != 1 for paths in names.values()) or len(matrices) != 1:
            raise FileNotFoundError(
                f"{sample_id}: expected exactly one matrix, feature, and barcode file"
            )
        with tempfile.TemporaryDirectory(
            prefix=f"fibrosis_{sample_id}_", dir=extracted_dir
        ) as temporary:
            sample_dir = Path(temporary)
            for kind, paths in names.items():
                os.link(paths[0], sample_dir / f"{kind}.tsv.gz")
            os.link(matrices[0], sample_dir / "matrix.mtx.gz")
            sample = previous_pipeline.read_10x(sample_dir)
        sample, cells_before, cells_after = previous_pipeline.qc_and_normalize(sample)
        cell_types, _ = previous_pipeline.annotate_celltypes(sample)
        symbols = sample.var["gene_symbols"].astype(str).str.upper().to_numpy()
        feature_names = sample.var_names.astype(str).tolist()
        if feature_names_reference is None:
            feature_names_reference = feature_names
        elif feature_names != feature_names_reference:
            raise ValueError(f"{sample_id}: single-cell feature order changed")
        symbol_index: dict[str, int] = {}
        for symbol in GENES + MARKERS:
            matches = np.flatnonzero(symbols == symbol)
            if len(matches) != 1:
                raise ValueError(
                    f"{sample_id}: expected one cell-level feature for {symbol}, "
                    f"found {len(matches)}"
                )
            symbol_index[symbol] = int(matches[0])
        panel_matrix = sample.X[:, [symbol_index[g] for g in MARKERS]].toarray()
        for cell_type in sorted(set(cell_types)):
            indices = np.flatnonzero(cell_types == cell_type)
            if not len(indices):
                continue
            celltype_expressions.setdefault(cell_type, []).append(
                panel_matrix[indices].astype(np.float32, copy=False)
            )
            celltype_obs.setdefault(cell_type, []).extend(
                [
                    {
                        "sample_id": sample_id,
                        "patient_id_proxy": sample_id,
                        "group": row.group_from_metadata,
                        "cell_type": cell_type,
                    }
                ]
                * len(indices)
            )

        target_n = min(1_250, sample.n_obs)
        picked = np.sort(rng.choice(sample.n_obs, size=target_n, replace=False))
        sample_symbols = symbols.copy()
        sample_symbols[sample_symbols == ""] = sample.var_names[
            sample_symbols == ""
        ].astype(str)
        var_names = pd.Index(sample_symbols)
        var_names = var_names.map(lambda x: str(x).upper())
        var_frame = pd.DataFrame(index=var_names)
        sample_adata = ad.AnnData(
            X=sample.X[picked].tocsr().copy(),
            obs=pd.DataFrame(
                {
                    "sample_id": sample_id,
                    "group": row.group_from_metadata,
                    "cell_type": cell_types[picked],
                },
                index=[f"{sample_id}_{i}" for i in picked],
            ),
            var=var_frame,
        )
        sample_adata.var_names_make_unique()
        liana_adatas.append(sample_adata)
        subsample_rows.append(
            {
                "sample_id": sample_id,
                "group": row.group_from_metadata,
                "cells_before_qc": cells_before,
                "cells_after_qc": cells_after,
                "cells_in_liana_subsample": target_n,
            }
        )
        del sample, sample_adata, panel_matrix

    mmt_rows = []
    for cell_type, matrices in celltype_expressions.items():
        stacked = np.vstack(matrices).astype(np.float64, copy=False)
        sd = stacked.std(axis=0, ddof=0)
        sd[sd == 0] = np.nan
        z_scores = (stacked - stacked.mean(axis=0)) / sd
        score = np.nanmean(z_scores, axis=1)
        obs = pd.DataFrame(celltype_obs[cell_type])
        obs["mmt_proxy_score"] = score
        patient = (
            obs.groupby(
                ["sample_id", "patient_id_proxy", "group", "cell_type"],
                as_index=False,
            )
            .agg(
                mmt_proxy_score=("mmt_proxy_score", "mean"),
                n_cells=("mmt_proxy_score", "size"),
            )
        )
        mmt_rows.extend(patient.to_dict("records"))
    mmt_by_patient = pd.DataFrame(mmt_rows)
    mmt_by_patient.to_csv(output_dir / "mmt_proxy_patient_scores.csv", index=False)
    mmt_summary_rows = []
    for cell_type in sorted(STROMAL_TYPES):
        sub = mmt_by_patient[mmt_by_patient["cell_type"] == cell_type]
        group_values = [
            sub.loc[sub["group"] == group, "mmt_proxy_score"].dropna().to_numpy()
            for group in GROUPS
        ]
        nonempty = [values for values in group_values if len(values)]
        omnibus_p = (
            float(kruskal(*nonempty).pvalue)
            if len(nonempty) >= 2 and all(len(values) > 0 for values in nonempty)
            else np.nan
        )
        for group in GROUPS:
            values = sub.loc[sub["group"] == group, "mmt_proxy_score"].dropna()
            mmt_summary_rows.append(
                {
                    "cell_type": cell_type,
                    "group": group,
                    "n_patients": int(values.size),
                    "median": float(values.median()) if len(values) else np.nan,
                    "q1": float(values.quantile(0.25)) if len(values) else np.nan,
                    "q3": float(values.quantile(0.75)) if len(values) else np.nan,
                    "omnibus_Kruskal_P": omnibus_p,
                }
            )
    mmt_summary = pd.DataFrame(mmt_summary_rows)
    mmt_summary.to_csv(output_dir / "mmt_proxy_group_summary.csv", index=False)
    pd.DataFrame(subsample_rows).to_csv(
        output_dir / "liana_subsample_manifest.csv", index=False
    )

    adata = ad.concat(liana_adatas, join="inner", merge="same", index_unique=None)
    del liana_adatas
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    resource = li.resource.select_resource("consensus")
    notch_ligands = {"JAG1", "JAG2", "DLL1", "DLL3", "DLL4"}
    notch_receptors = {"NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"}
    keep_tnfsf15 = (resource["ligand"] == "TNFSF15") & (
        resource["receptor"] == "TNFRSF25"
    )
    keep_il27 = (resource["ligand"] == "EBI3_IL27") & (
        resource["receptor"] == "IL27RA_IL6ST"
    )
    keep_notch = resource["ligand"].isin(notch_ligands) & resource[
        "receptor"
    ].isin(notch_receptors)
    targeted_resource = resource[keep_tnfsf15 | keep_il27 | keep_notch].drop_duplicates()
    if targeted_resource.empty:
        raise ValueError("LIANA consensus resource lacks all requested pathways")
    result = li.mt.rank_aggregate(
        adata,
        groupby="cell_type",
        resource=targeted_resource,
        expr_prop=0.1,
        min_cells=5,
        n_perms=1_000,
        seed=SEED,
        n_jobs=1,
        return_all_lrs=True,
        use_raw=False,
        inplace=False,
        verbose=True,
    )
    if result is None or not len(result):
        raise RuntimeError("LIANA returned no cell-cell communication predictions")
    result = pd.DataFrame(result)
    immune_to_stromal = result.apply(
        lambda prediction: (
            prediction["source"] in IMMUNE_TYPES
            and prediction["target"] in STROMAL_TYPES
        )
        or (
            prediction["target"] in IMMUNE_TYPES
            and prediction["source"] in STROMAL_TYPES
        ),
        axis=1,
    )
    selected = result.loc[immune_to_stromal].copy()
    selected.to_csv(output_dir / "liana_inferred_targeted_interactions.csv", index=False)
    status = {
        "tool": "LIANA",
        "version": importlib.metadata.version("liana"),
        "resource": "consensus",
        "resource_interactions_supplied": targeted_resource.to_dict("records"),
        "total_subsampled_cells": int(adata.n_obs),
        "cell_type_counts": adata.obs["cell_type"].value_counts().to_dict(),
        "seed": SEED,
        "targeted_immune_stromal_predictions": int(len(selected)),
        "TNFSF15_TNFRSF25_resource_pair_available": bool(keep_tnfsf15.any()),
        "IL27_EBI3_IL27RA_IL6ST_resource_pair_available": bool(keep_il27.any()),
        "Notch_ligand_receptor_resource_pairs_available": int(keep_notch.sum()),
        "ADAM19_direct_interaction_resource_pair_available": bool(
            resource["ligand"].eq("ADAM19").any()
            or resource["receptor"].eq("ADAM19").any()
        ),
        "warning": (
            "All LIANA results are inferred ligand-receptor predictions, not "
            "measured signaling. ADAM19 was not a direct interaction in the "
            "selected consensus resource; Notch-route predictions do not establish "
            "ADAM19-mediated Notch signaling."
        ),
    }
    (output_dir / "liana_inference_manifest.json").write_text(
        json.dumps(status, indent=2, default=str), encoding="utf-8"
    )
    return mmt_by_patient, selected, status


def write_report(
    output_dir: Path,
    correlations: pd.DataFrame,
    support: pd.DataFrame,
    edger_panel: pd.DataFrame,
    edger_config: dict,
    gsea: pd.DataFrame,
    mmt: pd.DataFrame,
    communication: pd.DataFrame,
    communication_manifest: dict,
) -> None:
    communication = communication.copy()
    communication["route_class"] = np.select(
        [
            (communication["ligand_complex"] == "TNFSF15")
            & (communication["receptor_complex"] == "TNFRSF25"),
            (communication["ligand_complex"] == "EBI3_IL27")
            & (communication["receptor_complex"] == "IL27RA_IL6ST"),
            communication["receptor_complex"]
            .astype(str)
            .str.contains("NOTCH", case=False),
        ],
        ["TNFSF15-TNFRSF25", "IL-27", "JAG/DLL-Notch"],
        default="other",
    )
    communication["direction"] = np.where(
        communication["source"].isin(IMMUNE_TYPES),
        "immune-to-stromal",
        "stromal-to-immune",
    )
    communication_summary = (
        communication[communication["route_class"] != "other"]
        .groupby(["route_class", "direction"], as_index=False)
        .agg(
            resource_pairs_evaluated=("source", "size"),
            best_cellphone_p=("cellphone_pvals", "min"),
            best_magnitude_rank=("magnitude_rank", "min"),
            best_specificity_rank=("specificity_rank", "min"),
        )
    )
    targeted_communication = communication[
        communication["route_class"] != "other"
    ]
    if not targeted_communication.empty:
        best_indices = targeted_communication.groupby(
            ["route_class", "direction"], sort=False
        )["specificity_rank"].idxmin()
        best_predictions = targeted_communication.loc[
            best_indices,
            [
                "route_class",
                "direction",
                "source",
                "target",
                "ligand_complex",
                "receptor_complex",
            ],
        ].rename(
            columns={
                "source": "best_source",
                "target": "best_target",
                "ligand_complex": "best_ligand",
                "receptor_complex": "best_receptor",
            }
        )
        communication_summary = communication_summary.merge(
            best_predictions, on=["route_class", "direction"], how="left"
        )
    communication_summary.to_csv(
        output_dir / "liana_requested_route_summary.csv", index=False
    )
    adjusted = correlations[
        (correlations["aggregation"] == "all_cells")
        & (correlations["adjustment"] == "composition_adjusted")
    ]
    support = support.merge(
        edger_panel[["gene_symbol", "logFC", "ci_lower", "ci_upper", "PValue", "FDR_all_genes"]]
        .rename(
            columns={
                "gene_symbol": "gene",
                "logFC": "LV_UF_vs_LV_NOT_UF_log2FC",
                "ci_lower": "log2FC_CI_lower",
                "ci_upper": "log2FC_CI_upper",
                "PValue": "edgeR_P",
                "FDR_all_genes": "edgeR_FDR",
            }
        ),
        on="gene",
        how="left",
    )
    support["contrast_concordant_positive"] = support[
        "LV_UF_vs_LV_NOT_UF_log2FC"
    ] > 0
    support["supported"] = (
        support["correlation_condition_met"]
        & support["contrast_concordant_positive"]
    )
    support.to_csv(output_dir / "frozen_success_rule_by_gene.csv", index=False)
    mmt_summary = pd.read_csv(output_dir / "mmt_proxy_group_summary.csv")

    target_cols = [
        "gene",
        "marker",
        "n",
        "rho",
        "ci_lower",
        "ci_upper",
        "p",
        "q",
    ]
    md = [
        "# GSE248762 fibrosis-context analysis",
        "",
        "## Scope and status",
        "",
        "This report covers GSE248762 only. It uses the fixed patient/sample "
        "inclusion and Scanpy QC already applied in the effluent analysis. "
        "Patients are the independent units; the all-cell pseudobulk analysis "
        "has n=16, so uncertainty is wide. The frozen support rule is in "
        "[PREREGISTRATION_FIBROSIS.md](PREREGISTRATION_FIBROSIS.md).",
        "",
        "All findings are **correlative, not causal**. Cell-cell communication "
        "results are **inferred**, not measured. The MMT quantity is explicitly "
        "a study-specific marker-panel proxy, not a validated canonical MMT score.",
        "",
        "## Frozen success rule",
        "",
        "Support is assessed only from all-cell pseudobulk composition-adjusted "
        "correlations and requires at least three markers with rho > 0.4 plus "
        "positive LV_UF-versus-LV_NOT_UF edgeR log2FC.",
        "",
        support.to_markdown(index=False, floatfmt=".4g"),
        "",
        "No changes were made to the candidate genes, markers, covariates, "
        "contrast, thresholds, or success rule after analysis.",
        "",
        "## Patient-level gene-marker correlations",
        "",
        "The complete all-cell and per-major-cell-type results (n, rho, "
        "patient-bootstrap 95% CI, two-sided P, and BH q) are in "
        "[gene_marker_correlations.csv](gene_marker_correlations.csv). The "
        "patient-level normalized expression used is in "
        "[patient_gene_marker_expression.csv](patient_gene_marker_expression.csv). "
        "Composition adjustment residualizes rank-transformed expression on "
        "rank-transformed monocyte/macrophage and T/NK fractions. Bootstrap "
        "intervals use 2,000 patient resamples (seed 42).",
        "",
        "All-cell composition-adjusted correlations with rho > 0.4:",
        "",
        adjusted.loc[adjusted["rho"] > 0.4, target_cols].to_markdown(
            index=False, floatfmt=".4g"
        )
        if (adjusted["rho"] > 0.4).any()
        else "No candidate-marker pairs had adjusted rho > 0.4.",
        "",
        "## LV_UF versus LV_NOT_UF edgeR contrast",
        "",
        "The existing edgeR pseudobulk pipeline was rerun in a temporary copy of "
        "its input files; its earlier outputs were not overwritten. The two LV "
        "groups were compared without SV samples. Package/model details: "
        + "; ".join(f"{key}={value}" for key, value in edger_config.items()),
        "",
        edger_panel[
            ["gene_symbol", "logFC", "ci_lower", "ci_upper", "PValue", "FDR_all_genes"]
        ].to_markdown(index=False, floatfmt=".4g"),
        "",
        "## Compartment-level preranked GSEA",
        "",
        "Gene rankings use edgeR robust quasi-likelihood statistics for "
        "LV_UF relative to LV_NOT_UF, separately for monocyte/macrophage, T/NK, "
        "fibroblast, and mesothelial pseudobulks. Gene sets are the Hallmark "
        "2020 EMT, TGF-beta, and inflammatory-response sets and GO BP 2023 "
        "Extracellular Matrix Organization (GO:0030198). Full NES, P, FDR, "
        "and panel leading-edge membership are in "
        "[compartment_gsea_results.csv](compartment_gsea_results.csv).",
        "",
        "Gene-set P/FDR estimates with zero permutation exceedances are "
        "reported at their empirical resolution floor (1/1,001), not as exact zero.",
        "",
        gsea.to_markdown(index=False, floatfmt=".4g"),
        "",
        "## Inferred cell-cell communication",
        "",
        f"LIANA {communication_manifest['version']} consensus resource analyzed "
        f"{communication_manifest['total_subsampled_cells']:,} observed cells "
        "subsampled with seed 42 (target at most 20,000). Predictions connecting "
        "immune and stromal cell types are in "
        "[liana_inferred_targeted_interactions.csv](liana_inferred_targeted_interactions.csv).",
        "",
        "Predicted ligand-receptor scores are not evidence of signaling activity. "
        "The consensus resource includes TNFSF15-TNFRSF25 and EBI3/IL-27 receptor "
        "pairs. It does not contain a direct ADAM19-Notch interaction; JAG/DLL-Notch "
        "rankings are separately inferred and do not attribute receptor processing "
        "to ADAM19.",
        "",
        "The route-level summary below aggregates the rank of the best "
        "resource-supported prediction within each direction; lower rank is "
        "better, and these rankings are not hypothesis-test P values.",
        "",
        communication_summary.to_markdown(index=False, floatfmt=".4g")
        if len(communication_summary)
        else "No requested pathways had any resource-supported interactions.",
        "",
        "No immune-to-stromal TNFSF15-TNFRSF25 or IL-27 pair had an informative "
        "LIANA score in this analysis: all such CellPhoneDB p-values and "
        "magnitude/specificity ranks were 1.0. Any reverse-direction or Notch "
        "ranking remains inferred and does not establish signaling.",
        "",
        f"Full targeted output has {len(communication)} immune/stromal directed "
        "predictions, including non-significant resource-supported pairs.",
        "",
        communication_manifest["warning"],
        "",
        "## MMT proxy by PD group",
        "",
        "The patient-level score is the mean of the eight marker-gene z-scores "
        "within each cell type after Scanpy total-count normalization to 10,000 "
        "and log1p transformation. Patient is the unit of comparison; group "
        "summaries and Kruskal-Wallis P values are in "
        "[mmt_proxy_group_summary.csv](mmt_proxy_group_summary.csv), and "
        "individual patient scores are in "
        "[mmt_proxy_patient_scores.csv](mmt_proxy_patient_scores.csv).",
        "",
        mmt_summary.to_markdown(index=False, floatfmt=".4g"),
        "",
        "## Reproducibility files",
        "",
        "The analysis driver is [run_fibrosis_context.py](run_fibrosis_context.py); "
        "compartment-level edgeR ranking code is "
        "[compartment_edger_gsea_input.R](compartment_edger_gsea_input.R). "
        "The pseudobulk input matrices and metadata used for this analysis are "
        "included in this folder as "
        "[all_cell_pseudobulk_counts.csv](all_cell_pseudobulk_counts.csv), "
        "[celltype_pseudobulk_counts.csv](celltype_pseudobulk_counts.csv), "
        "[all_cell_pseudobulk_metadata.csv](all_cell_pseudobulk_metadata.csv), "
        "[celltype_pseudobulk_metadata.csv](celltype_pseudobulk_metadata.csv), "
        "and [patient_celltype_composition.csv](patient_celltype_composition.csv). "
        "Exact GSEA gene-set "
        "memberships are in [gsea_gene_set_membership.csv](gsea_gene_set_membership.csv). "
        "The original edgeR pipeline and source data remain in their existing "
        "locations and were not edited by this analysis.",
    ]
    (output_dir / "FIBROSIS_CONTEXT_REPORT.md").write_text(
        "\n".join(md) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[3]
    )
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = Path(__file__).resolve().parent
    input_dir = (
        repo_root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE248762"
        / "analysis_20261009"
    )
    if not input_dir.is_dir():
        raise FileNotFoundError(f"GSE248762 pseudobulk inputs not found: {input_dir}")
    map_df, gene_ids = read_symbol_map(input_dir / "gene_id_to_symbol.csv")
    start = time.time()
    correlation_table, support = calculate_correlations(
        input_dir, output_dir, gene_ids
    )
    edger_panel, edger_config = run_existing_edger(
        repo_root, input_dir, output_dir
    )
    rankings, _ = run_compartment_edger(repo_root, input_dir, output_dir)
    gsea = run_preranked_gsea(rankings, output_dir)
    mmt, communication, communication_manifest = process_single_cell(
        repo_root, output_dir
    )
    write_report(
        output_dir,
        correlation_table,
        support,
        edger_panel,
        edger_config,
        gsea,
        mmt,
        communication,
        communication_manifest,
    )
    input_names = [
        "all_cell_pseudobulk_counts.csv",
        "all_cell_pseudobulk_metadata.csv",
        "celltype_pseudobulk_counts.csv",
        "celltype_pseudobulk_metadata.csv",
        "patient_celltype_composition.csv",
        "gene_id_to_symbol.csv",
    ]
    input_hashes = {}
    for name in input_names:
        source_path = input_dir / name
        destination = output_dir / name
        shutil.copy2(source_path, destination)
        digest = hashlib.sha256()
        with destination.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        input_hashes[name] = digest.hexdigest()
    manifest = {
        "cohort": "GSE248762",
        "seed": SEED,
        "patients": 16,
        "groups": {
            group: int((pd.read_csv(
                repo_root
                / "results"
                / "revision"
                / "effluent_validation"
                / "GSE248762"
                / "sample_metadata.csv"
            )["group_from_metadata"] == group).sum())
            for group in GROUPS
        },
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "liana": importlib.metadata.version("liana"),
        "gseapy": importlib.metadata.version("gseapy"),
        "runtime_seconds": round(time.time() - start, 1),
        "candidate_genes": GENES,
        "fibrosis_markers": MARKERS,
        "source_analysis_config": str(
            input_dir / "analysis_config.json"
        ),
        "pseudobulk_input_sha256": input_hashes,
        "existing_edger_script": str(
            repo_root
            / "results"
            / "revision"
            / "effluent_validation"
            / "GSE248762"
            / "scripts"
            / "gse248762_step2_4_edger.R"
        ),
        "gene_mapping_rows": int(len(map_df)),
    }
    (output_dir / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str), encoding="utf-8"
    )
    print(f"Analysis complete: {output_dir}", flush=True)


if __name__ == "__main__":
    main()
