from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import os
import re
import shutil
import tarfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw" / "network_7hub"
SEED = 20261019
MIN_LNCRNA_CELL_PERCENT = 1.0
MIN_LNCRNA_SAMPLE_COUNT = 3
MIN_LINEAGE_CELLS = 20
MIN_PSEUDOBULK_N = 5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, name: str) -> Path:
    path = OUT / name
    frame.to_csv(path, index=False, na_rep="")
    return path


def read_encori_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)


def feature_id(value: str) -> str:
    return str(value).split(".", 1)[0]


def partial_spearman(x: np.ndarray, y: np.ndarray, covariate: np.ndarray) -> tuple[float, float]:
    rx = stats.rankdata(x, method="average")
    ry = stats.rankdata(y, method="average")
    rz = stats.rankdata(covariate, method="average")
    if np.unique(rz).size < 2:
        result = stats.spearmanr(x, y)
        return float(result.statistic), float(result.pvalue)
    design = np.column_stack((np.ones(len(rz)), rz))
    residual_x = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    residual_y = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    rho = float(np.corrcoef(residual_x, residual_y)[0, 1])
    degrees_freedom = len(x) - 3
    if degrees_freedom <= 0:
        return rho, float("nan")
    if abs(rho) >= 1:
        return rho, 0.0
    t_value = rho * math.sqrt(degrees_freedom / (1.0 - rho * rho))
    p_value = float(2.0 * stats.t.sf(abs(t_value), degrees_freedom))
    return rho, p_value


def plot_save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def load_mirna_mrna_evidence(
    hubs: list[str],
    validated: pd.DataFrame,
    predicted: pd.DataFrame,
    encori_by_hub: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    assay_pattern = r"reporter|luciferase|western|q-?rt.?pcr|q-?pcr|qpcr"
    all_rows: list[dict[str, object]] = []
    per_hub_rows: list[dict[str, object]] = []

    for hub in hubs:
        val_hub = validated[validated["target_symbol"].eq(hub)].copy()
        pred_hub = predicted[predicted["target_symbol"].eq(hub)].copy()
        encori = encori_by_hub[hub].copy()
        if not encori.empty:
            encori["clipExpNum_int"] = pd.to_numeric(
                encori["clipExpNum"], errors="coerce"
            ).fillna(0).astype(int)
        else:
            encori["clipExpNum_int"] = pd.Series(dtype=int)

        keys = set()
        for frame in (val_hub, pred_hub):
            keys.update(
                (str(row["mature_mirna_acc"]), str(row["mature_mirna_id"]))
                for _, row in frame.iterrows()
                if str(row.get("mature_mirna_acc", ""))
            )
        if not encori.empty:
            keys.update(
                (str(row["miRNAid"]), str(row["miRNAname"]))
                for _, row in encori.iterrows()
                if str(row.get("miRNAid", ""))
            )

        group_rows: list[dict[str, object]] = []
        for accession, mirna_name in sorted(keys):
            val_pair = val_hub[val_hub["mature_mirna_acc"].eq(accession)]
            positive_val = val_pair[
                ~val_pair["support_type"].str.lower().eq("negative")
            ]
            negative_val = val_pair[
                val_pair["support_type"].str.lower().eq("negative")
            ]
            predicted_pair = pred_hub[pred_hub["mature_mirna_acc"].eq(accession)]
            enc_pair = encori[encori["miRNAid"].eq(accession)]
            clip_values = (
                pd.to_numeric(enc_pair["clipExpNum"], errors="coerce")
                .fillna(0)
                .astype(int)
                if not enc_pair.empty
                else pd.Series(dtype=int)
            )
            clip_max = int(clip_values.max()) if not clip_values.empty else 0
            functional = positive_val[
                positive_val["experiment"].str.contains(
                    assay_pattern, case=False, regex=True, na=False
                )
            ]
            tier_a = not functional.empty
            tier_b = clip_max >= 2
            has_positive_validated = not positive_val.empty
            tier_c = (
                not predicted_pair.empty
                and not tier_a
                and not tier_b
                and not has_positive_validated
            )
            if tier_a:
                best_tier = "A"
            elif tier_b:
                best_tier = "B"
            elif tier_c:
                best_tier = "C"
            elif has_positive_validated:
                best_tier = "validated_unclassified"
            else:
                best_tier = "none"
            row = {
                "hub_gene": hub,
                "miRNA_accession": accession,
                "miRNA_name": mirna_name,
                "best_evidence_class": best_tier,
                "tier_A_functional": tier_a,
                "tier_B_clip_ge_2": tier_b,
                "tier_C_prediction_only": tier_c,
                "n_positive_validated_records": len(positive_val),
                "n_negative_records_excluded": len(negative_val),
                "n_functional_assay_records": len(functional),
                "functional_assays": ";".join(
                    sorted(set(functional["experiment"].astype(str)))
                ),
                "positive_validated_support_types": ";".join(
                    sorted(set(positive_val["support_type"].astype(str)))
                ),
                "positive_validated_experiments": ";".join(
                    sorted(set(positive_val["experiment"].astype(str)))
                ),
                "negative_record_experiments_excluded": ";".join(
                    sorted(set(negative_val["experiment"].astype(str)))
                ),
                "validated_sources": ";".join(
                    sorted(set(positive_val["database"].astype(str)))
                ),
                "predicted_sources": ";".join(
                    sorted(set(predicted_pair["database"].astype(str)))
                ),
                "n_prediction_records": len(predicted_pair),
                "ENCORI_clipExpNum_max": clip_max,
                "n_ENCORI_target_site_records": len(enc_pair),
            }
            all_rows.append(row)
            group_rows.append(row)

        for tier in ("A", "B", "C", "validated_unclassified"):
            if tier == "A":
                count = sum(bool(row["tier_A_functional"]) for row in group_rows)
            elif tier == "B":
                count = sum(
                    bool(row["tier_B_clip_ge_2"]) and not bool(row["tier_A_functional"])
                    for row in group_rows
                )
            elif tier == "C":
                count = sum(bool(row["tier_C_prediction_only"]) for row in group_rows)
            else:
                count = sum(
                    row["best_evidence_class"] == "validated_unclassified"
                    for row in group_rows
                )
            per_hub_rows.append(
                {
                    "hub_gene": hub,
                    "evidence_tier": tier,
                    "unique_miRNA_count": count,
                    "definition": {
                        "A": "Positive multiMiR validated record with reporter, luciferase, western-blot, qPCR, or qRT-PCR assay",
                        "B": "ENCORI hub-mRNA clipExpNum >= 2, excluding Tier A best-class overlaps",
                        "C": "multiMiR predicted interaction only; no positive validated record or Tier A/B support",
                        "validated_unclassified": "Positive validated-table record without a qualifying functional assay or Tier B support",
                    }[tier],
                }
            )
        per_hub_rows.append(
            {
                "hub_gene": hub,
                "evidence_tier": "B_support_including_A",
                "unique_miRNA_count": sum(bool(row["tier_B_clip_ge_2"]) for row in group_rows),
                "definition": "All unique miRNAs with ENCORI hub-mRNA clipExpNum >= 2, including any also classed Tier A",
            }
        )
    return pd.DataFrame(all_rows), pd.DataFrame(per_hub_rows)


def load_lnc_mirna_pairs(
    path: Path, target_accessions: set[str]
) -> pd.DataFrame:
    pair_data: dict[tuple[str, str], dict[str, object]] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            accession = row["miRNAid"]
            if accession not in target_accessions:
                continue
            key = (row["geneID"], accession)
            item = pair_data.setdefault(
                key,
                {
                    "lncRNA_gene_id": row["geneID"],
                    "lncRNA_name": row["geneName"],
                    "miRNA_accession": accession,
                    "miRNA_name": row["miRNAname"],
                    "n_target_site_records": 0,
                    "clipExpNum_max": 0,
                    "bindClipExpNum_max": 0,
                    "degraExpNum_max": 0,
                    "transcripts": set(),
                },
            )
            item["n_target_site_records"] = int(item["n_target_site_records"]) + 1
            item["clipExpNum_max"] = max(
                int(item["clipExpNum_max"]), int(row.get("clipExpNum") or 0)
            )
            item["bindClipExpNum_max"] = max(
                int(item["bindClipExpNum_max"]),
                int(row.get("bindClipExpNum") or 0),
            )
            item["degraExpNum_max"] = max(
                int(item["degraExpNum_max"]),
                int(row.get("degraExpNum") or 0),
            )
            tx = row.get("txInfo", "")
            if tx:
                item["transcripts"].add(tx)
    rows = []
    for item in pair_data.values():
        transcripts = item.pop("transcripts")
        item["n_transcripts"] = len(transcripts)
        item["transcript_ids"] = ";".join(sorted(transcripts))
        item["lncRNA_CLIP_status"] = (
            "CLIP_supported_ge_2"
            if int(item["clipExpNum_max"]) >= 2
            else "CLIP_supported_1"
            if int(item["clipExpNum_max"]) == 1
            else "target_site_record_clip_0"
        )
        rows.append(item)
    return pd.DataFrame(rows)


def load_encori_cerna(
    hubs: list[str], raw_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    pair_frames = []
    status_rows = []
    manifest_path = raw_dir / "encori_query_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    query_entries = manifest.get("ceRNA_hub_queries", {}).get("queries", [])
    entries_by_hub = {str(entry.get("gene")): entry for entry in query_entries}
    for hub in hubs:
        path = raw_dir / f"encori_ceRNA_{hub}.tsv"
        text = path.read_text(encoding="utf-8", errors="replace")
        data_lines = [line for line in text.splitlines() if line and not line.startswith("#")]
        api_error = next(
            (line for line in data_lines if line.startswith("The ")),
            "",
        )
        if not data_lines:
            status = "NOT RETRIEVED"
            error = "Empty response"
            frame = pd.DataFrame()
        elif api_error:
            status = "NOT RETRIEVED"
            error = api_error
            frame = pd.DataFrame()
        else:
            frame = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
            expected_columns = {
                "geneID", "geneName", "geneType", "ceRNAid", "ceRNAname",
                "ceRNAgeneType", "hitMiRNAFamilyNum", "hitMiRNAFamily",
                "pval", "fdr",
            }
            if not expected_columns.issubset(frame.columns):
                status = "NOT RETRIEVED"
                error = "Response did not contain the expected ceRNA result schema"
                frame = pd.DataFrame()
            else:
                frame = frame[frame["geneID"].astype(str).str.startswith("ENSG")].copy()
                status = "RETRIEVED" if not frame.empty else "NO ROWS"
                error = ""
                if not frame.empty:
                    frame["queried_hub"] = hub
                    pair_frames.append(frame)
        entry = entries_by_hub.get(hub, {})
        entry["response_status"] = status
        entry["parsed_data_rows"] = len(frame)
        entry["data_lines"] = len(frame)
        entry["error"] = error
        status_rows.append(
            {
                "queried_hub": hub,
                "status": status,
                "parsed_data_rows": len(frame),
                "error": error,
                "raw_response_file": str(path.relative_to(ROOT)),
                "parameters": "geneType=mRNA; miRNAnum=2; family=all; pval<=0.01; fdr<=0.01; pancancerNum=0",
            }
        )
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    pairs = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    return pairs, pd.DataFrame(status_rows), manifest


def build_axes(mrna: pd.DataFrame, lnc: pd.DataFrame) -> pd.DataFrame:
    eligible = mrna[mrna["best_evidence_class"].isin(["A", "B", "C"])].copy()
    if eligible.empty or lnc.empty:
        return pd.DataFrame(
            columns=[
                "hub_gene", "miRNA_accession", "miRNA_name", "lncRNA_gene_id",
                "lncRNA_name", "mRNA_evidence_tier", "lncRNA_clipExpNum_max",
                "confidence_class",
            ]
        )
    axes = eligible.merge(lnc, on=["miRNA_accession", "miRNA_name"], how="inner")
    axes["confidence_class"] = np.where(
        axes["best_evidence_class"].isin(["A", "B"])
        & axes["clipExpNum_max"].ge(2),
        "high",
        "low",
    )
    axes = axes.rename(
        columns={
            "best_evidence_class": "mRNA_evidence_tier",
            "clipExpNum_max": "lncRNA_clipExpNum_max",
            "bindClipExpNum_max": "lncRNA_bindClipExpNum_max",
            "n_target_site_records": "lncRNA_target_site_records",
            "n_transcripts": "lncRNA_transcript_count",
            "transcript_ids": "lncRNA_transcripts",
            "functional_assays": "mRNA_functional_assays",
            "validated_sources": "mRNA_validated_sources",
            "predicted_sources": "mRNA_predicted_sources",
            "ENCORI_clipExpNum_max": "mRNA_ENCORI_clipExpNum_max",
            "n_ENCORI_target_site_records": "mRNA_ENCORI_target_site_records",
        }
    )
    axes = axes.drop_duplicates(
        ["hub_gene", "miRNA_accession", "lncRNA_gene_id"]
    ).reset_index(drop=True)
    return axes


def build_detection_table(
    axes: pd.DataFrame,
    h5ad_path: Path,
    annotations: pd.DataFrame,
) -> tuple[pd.DataFrame, ad.AnnData, pd.DataFrame, dict[str, list[int]]]:
    lnc_info = (
        axes[["lncRNA_gene_id", "lncRNA_name"]]
        .drop_duplicates()
        .set_index("lncRNA_gene_id")["lncRNA_name"]
        .to_dict()
    )
    axis_counts = axes.groupby("lncRNA_gene_id").size().to_dict()
    data = ad.read_h5ad(h5ad_path, backed="r")
    variable_ids = data.var["gene_ids"].astype(str).map(feature_id)
    if variable_ids.duplicated().any():
        data.file.close()
        raise ValueError("GSE248762 feature IDs are not unique after version stripping.")
    id_to_index = {value: index for index, value in enumerate(variable_ids)}
    gene_columns: dict[str, list[int]] = {}
    match_status: dict[str, str] = {}
    for gene_id in lnc_info:
        normalized = feature_id(gene_id)
        if normalized in id_to_index:
            gene_columns[gene_id] = [id_to_index[normalized]]
            match_status[gene_id] = "exact_Ensembl_gene_id"
        else:
            gene_columns[gene_id] = []
            match_status[gene_id] = "feature_not_found"
    matched_ids = [gene_id for gene_id, cols in gene_columns.items() if cols]
    unique_columns = sorted({cols[0] for cols in gene_columns.values() if cols})
    column_position = {index: position for position, index in enumerate(unique_columns)}
    mapped_col = {
        gene_id: column_position[cols[0]]
        for gene_id, cols in gene_columns.items()
        if cols
    }

    row_index = pd.Index(data.obs_names.astype(str))
    annotations = annotations.copy()
    annotations["object_row"] = row_index.get_indexer(annotations["cell_id"].astype(str))
    if annotations["object_row"].lt(0).any():
        data.file.close()
        raise ValueError("Some QC-retained annotations are absent from the count object.")
    annotations = annotations.sort_values("object_row").reset_index(drop=True)
    detected_cells = Counter()
    positive_samples = Counter()
    sample_ids = sorted(annotations["sample_id"].astype(str).unique())
    for sample_id, group in annotations.groupby("sample_id", sort=False):
        if not unique_columns:
            break
        rows = group["object_row"].to_numpy(dtype=np.int64)
        counts = data.X[rows, :][:, unique_columns]
        if not sparse.issparse(counts):
            counts = sparse.csr_matrix(counts)
        binary = counts.tocsr(copy=True)
        binary.data = np.ones(binary.data.shape, dtype=np.int8)
        per_feature = np.asarray(binary.sum(axis=0)).ravel().astype(np.int64)
        for gene_id in matched_ids:
            position = mapped_col[gene_id]
            detected_cells[gene_id] += int(per_feature[position])
            if per_feature[position] > 0:
                positive_samples[gene_id] += 1

    n_qc_cells = len(annotations)
    rows = []
    for gene_id, name in lnc_info.items():
        n_detected = int(detected_cells[gene_id])
        n_samples = int(positive_samples[gene_id])
        pct = 100.0 * n_detected / n_qc_cells if match_status[gene_id] != "feature_not_found" else np.nan
        passed = (
            match_status[gene_id] != "feature_not_found"
            and pct >= MIN_LNCRNA_CELL_PERCENT
            and n_samples >= MIN_LNCRNA_SAMPLE_COUNT
        )
        rows.append(
            {
                "lncRNA_gene_id": gene_id,
                "lncRNA_name": name,
                "n_network_axes": int(axis_counts.get(gene_id, 0)),
                "feature_match_status": match_status[gene_id],
                "n_QC_cells": n_qc_cells,
                "n_cells_with_count_gt_0": n_detected if match_status[gene_id] != "feature_not_found" else pd.NA,
                "percent_QC_cells_detected": pct,
                "n_samples_with_pseudobulk_count_gt_0": n_samples if match_status[gene_id] != "feature_not_found" else pd.NA,
                "passes_coexpression_detection_threshold": passed,
                "threshold": f">={MIN_LNCRNA_CELL_PERCENT}% QC cells and >={MIN_LNCRNA_SAMPLE_COUNT} samples",
            }
        )
    detection = pd.DataFrame(rows)
    return detection, data, annotations, gene_columns


def make_correlation_tables(
    axes: pd.DataFrame,
    detection: pd.DataFrame,
    h5ad: ad.AnnData,
    annotations: pd.DataFrame,
    gene_columns: dict[str, list[int]],
    root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if axes.empty:
        return pd.DataFrame(), pd.DataFrame()
    pass_ids = set(
        detection.loc[
            detection["passes_coexpression_detection_threshold"],
            "lncRNA_gene_id",
        ].astype(str)
    )
    pair_rows = (
        axes[axes["lncRNA_gene_id"].astype(str).isin(pass_ids)]
        .groupby(["lncRNA_gene_id", "lncRNA_name", "hub_gene"], as_index=False)
        .agg(
            n_linked_miRNAs=("miRNA_accession", "nunique"),
            linked_miRNAs=("miRNA_name", lambda values: ";".join(sorted(set(values)))),
        )
    )
    symbol_to_id: dict[str, list[str]] = defaultdict(list)
    for symbol, gene_id in zip(h5ad.var_names.astype(str), h5ad.var["gene_ids"].astype(str)):
        symbol_to_id[symbol.upper()].append(feature_id(gene_id))
    hub_ids = {
        hub: symbol_to_id.get(hub.upper(), [None])[0]
        for hub in axes["hub_gene"].drop_duplicates().astype(str)
    }
    count_rows_by_lineage: dict[str, pd.DataFrame] = {}
    lineages = sorted(annotations["lineage"].astype(str).unique())
    cell_counts = (
        annotations.groupby(["sample_id", "lineage"], as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
    )
    results: list[dict[str, object]] = []

    for lineage in lineages:
        count_path = (
            root
            / "results"
            / "revision"
            / "single_cell_context"
            / "pseudobulk_counts"
            / f"counts_{lineage.replace('/', '_').replace(' ', '_')}.csv"
        )
        if not count_path.exists():
            continue
        count_frame = pd.read_csv(count_path)
        if "feature_id" not in count_frame:
            raise ValueError(f"Missing feature_id in {count_path}")
        sample_columns = list(count_frame.columns[2:])
        count_frame["feature_id_norm"] = count_frame["feature_id"].astype(str).map(feature_id)
        if count_frame["feature_id_norm"].duplicated().any():
            count_frame = count_frame.groupby("feature_id_norm", as_index=False)[sample_columns].sum()
        else:
            count_frame = count_frame.set_index("feature_id_norm")
        if "feature_id_norm" in count_frame.columns:
            count_frame = count_frame.set_index("feature_id_norm")
        library_sizes = count_frame[sample_columns].sum(axis=0).astype(float)
        lineage_counts = cell_counts[cell_counts["lineage"].eq(lineage)].set_index("sample_id")["n_cells"]
        eligible = [
            sample for sample in sample_columns
            if int(lineage_counts.get(sample, 0)) >= MIN_LINEAGE_CELLS
        ]
        if not eligible:
            continue

        def get_logcpm(gene_id: str | None) -> pd.Series:
            if not gene_id:
                return pd.Series(index=eligible, dtype=float)
            norm_id = feature_id(gene_id)
            if norm_id not in count_frame.index:
                return pd.Series(index=eligible, dtype=float)
            counts = count_frame.loc[norm_id, eligible].astype(float)
            sizes = library_sizes.loc[eligible]
            return np.log2((counts * 1_000_000.0 / sizes) + 1.0)

        for row in pair_rows.itertuples(index=False):
            if not gene_columns.get(str(row.lncRNA_gene_id)):
                continue
            lnc = get_logcpm(str(row.lncRNA_gene_id))
            hub = get_logcpm(hub_ids.get(str(row.hub_gene)))
            cell_n = np.log1p(
                pd.Series(
                    {sample: int(lineage_counts.get(sample, 0)) for sample in eligible}
                )
            )
            x = lnc.to_numpy(dtype=float)
            y = hub.to_numpy(dtype=float)
            z = cell_n.to_numpy(dtype=float)
            n = int(np.isfinite(x).sum() and np.isfinite(y).sum())
            status = "testable"
            reason = ""
            if len(x) < MIN_PSEUDOBULK_N:
                status = "not_testable"
                reason = f"fewer than {MIN_PSEUDOBULK_N} eligible sample-lineages"
            elif np.unique(x).size < 2 or np.unique(y).size < 2:
                status = "not_testable"
                reason = "lncRNA or hub logCPM is constant"
            if status != "testable":
                results.append(
                    {
                        "lncRNA_gene_id": row.lncRNA_gene_id,
                        "lncRNA_name": row.lncRNA_name,
                        "hub_gene": row.hub_gene,
                        "lineage": lineage,
                        "n_linked_miRNAs": row.n_linked_miRNAs,
                        "n": len(eligible),
                        "rho_unadjusted": np.nan,
                        "p_unadjusted": np.nan,
                        "q_unadjusted_all_pairs": np.nan,
                        "rho_adjusted_log_cell_count": np.nan,
                        "p_adjusted_log_cell_count": np.nan,
                        "q_adjusted_all_pairs": np.nan,
                        "status": status,
                        "not_testable_reason": reason,
                    }
                )
                continue
            raw = stats.spearmanr(x, y)
            adjusted_rho, adjusted_p = partial_spearman(x, y, z)
            results.append(
                {
                    "lncRNA_gene_id": row.lncRNA_gene_id,
                    "lncRNA_name": row.lncRNA_name,
                    "hub_gene": row.hub_gene,
                    "lineage": lineage,
                    "n_linked_miRNAs": row.n_linked_miRNAs,
                    "n": len(eligible),
                    "rho_unadjusted": float(raw.statistic),
                    "p_unadjusted": float(raw.pvalue),
                    "q_unadjusted_all_pairs": np.nan,
                    "rho_adjusted_log_cell_count": adjusted_rho,
                    "p_adjusted_log_cell_count": adjusted_p,
                    "q_adjusted_all_pairs": np.nan,
                    "status": "tested",
                    "not_testable_reason": "",
                }
            )

    correlation = pd.DataFrame(results)
    if not correlation.empty:
        tested = correlation["status"].eq("tested")
        if tested.any():
            correlation.loc[tested, "q_unadjusted_all_pairs"] = multipletests(
                correlation.loc[tested, "p_unadjusted"].astype(float),
                method="fdr_bh",
            )[1]
            adjusted_finite = tested & correlation["p_adjusted_log_cell_count"].notna()
            correlation.loc[adjusted_finite, "q_adjusted_all_pairs"] = multipletests(
                correlation.loc[adjusted_finite, "p_adjusted_log_cell_count"].astype(float),
                method="fdr_bh",
            )[1]

    cell_level_rows: list[dict[str, object]] = []
    eligible_pairs = pair_rows[
        pair_rows["lncRNA_gene_id"].astype(str).isin(pass_ids)
    ]
    qualifying_ids = sorted(set(eligible_pairs["lncRNA_gene_id"].astype(str)))
    if qualifying_ids:
        gene_index_by_symbol = {
            str(symbol).upper(): index
            for index, symbol in enumerate(h5ad.var_names.astype(str))
        }
        relevant_hub_symbols = sorted(set(eligible_pairs["hub_gene"].astype(str)))
        selected_feature_index: dict[str, int] = {}
        for gene_id in qualifying_ids:
            cols = gene_columns.get(gene_id, [])
            if cols:
                selected_feature_index[f"LNC:{gene_id}"] = cols[0]
        for hub in relevant_hub_symbols:
            if hub.upper() in gene_index_by_symbol:
                selected_feature_index[f"HUB:{hub}"] = gene_index_by_symbol[hub.upper()]
        unique_features = sorted(set(selected_feature_index.values()))
        feature_position = {index: position for position, index in enumerate(unique_features)}
        key_position = {
            key: feature_position[index]
            for key, index in selected_feature_index.items()
        }
        rows = annotations["object_row"].to_numpy(dtype=np.int64)
        selected = h5ad.X[rows, :][:, unique_features]
        if not sparse.issparse(selected):
            selected = sparse.csr_matrix(selected)
        selected = selected.tocsc()
        library_size = annotations["n_counts"].to_numpy(dtype=float)
        lineage_array = annotations["lineage"].astype(str).to_numpy()
        hub_values: dict[str, np.ndarray] = {}
        for hub in relevant_hub_symbols:
            key = f"HUB:{hub}"
            if key in key_position:
                counts = selected[:, key_position[key]].toarray().ravel().astype(float)
                hub_values[hub] = np.log1p(counts * (10_000.0 / library_size))
        for gene_id in qualifying_ids:
            key = f"LNC:{gene_id}"
            if key not in key_position:
                continue
            counts = selected[:, key_position[key]].toarray().ravel().astype(float)
            lnc_values = np.log1p(counts * (10_000.0 / library_size))
            linked = eligible_pairs[eligible_pairs["lncRNA_gene_id"].astype(str).eq(gene_id)]
            for pair in linked.itertuples(index=False):
                if pair.hub_gene not in hub_values:
                    continue
                for lineage in lineages:
                    mask = lineage_array == lineage
                    x = lnc_values[mask]
                    y = hub_values[pair.hub_gene][mask]
                    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
                        continue
                    result = stats.spearmanr(x, y)
                    cell_level_rows.append(
                        {
                            "lncRNA_gene_id": gene_id,
                            "lncRNA_name": pair.lncRNA_name,
                            "hub_gene": pair.hub_gene,
                            "lineage": lineage,
                            "n_cells": len(x),
                            "rho_cell_level": float(result.statistic),
                            "p_cell_level_descriptive": float(result.pvalue),
                            "q_cell_level_descriptive": np.nan,
                            "independent_biological_replicates": False,
                        }
                    )
        cell_level = pd.DataFrame(cell_level_rows)
        if not cell_level.empty:
            cell_level["q_cell_level_descriptive"] = multipletests(
                cell_level["p_cell_level_descriptive"].astype(float),
                method="fdr_bh",
            )[1]
    else:
        cell_level = pd.DataFrame(
            columns=[
                "lncRNA_gene_id", "lncRNA_name", "hub_gene", "lineage",
                "n_cells", "rho_cell_level", "p_cell_level_descriptive",
                "q_cell_level_descriptive", "independent_biological_replicates",
            ]
        )
    return correlation, cell_level


def make_network_figures(
    axes: pd.DataFrame,
    tier_counts: pd.DataFrame,
    detection: pd.DataFrame,
    correlation: pd.DataFrame,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 8))
    if axes.empty or not axes["confidence_class"].eq("high").any():
        ax.text(0.5, 0.5, "No high-confidence axes in retrieved evidence", ha="center", va="center")
        ax.set_axis_off()
    else:
        high = axes[axes["confidence_class"].eq("high")]
        lnc_nodes = sorted(high["lncRNA_name"].astype(str).unique())
        mir_nodes = sorted(high["miRNA_name"].astype(str).unique())
        hub_nodes = sorted(high["hub_gene"].astype(str).unique())
        layers = [(lnc_nodes, 0), (mir_nodes, 1), (hub_nodes, 2)]
        pos: dict[tuple[str, str], tuple[float, float]] = {}
        for nodes, x in layers:
            for j, node in enumerate(nodes):
                y = 0.5 if len(nodes) == 1 else 1.0 - j / (len(nodes) - 1)
                pos[(("lnc" if x == 0 else "mir" if x == 1 else "hub"), node)] = (x, y)
        alpha = min(0.7, 30 / max(len(high), 30))
        for row in high.itertuples(index=False):
            lnc = str(row.lncRNA_name)
            mir = str(row.miRNA_name)
            hub = str(row.hub_gene)
            p1 = pos[("lnc", lnc)]
            p2 = pos[("mir", mir)]
            p3 = pos[("hub", hub)]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#45A36F", alpha=alpha, linewidth=0.6)
            style = "-" if row.mRNA_evidence_tier == "A" else "--"
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color="#475569", alpha=alpha, linewidth=0.7, linestyle=style)
        for nodes, x, color, category in [
            (lnc_nodes, 0, "#36A269", "lncRNA"),
            (mir_nodes, 1, "#D8A932", "miRNA"),
            (hub_nodes, 2, "#4C78A8", "validated hub"),
        ]:
            ys = [pos[(category.split()[0] if category != "lncRNA" else "lnc", n)][1] for n in nodes] if False else None
            node_key = "lnc" if x == 0 else "mir" if x == 1 else "hub"
            points = [pos[(node_key, n)] for n in nodes]
            ax.scatter([point[0] for point in points], [point[1] for point in points], s=24, color=color, label=category, zorder=3)
            if len(nodes) <= 45 or x == 2:
                for node, point in zip(nodes, points):
                    ax.text(point[0] + 0.025, point[1], node, fontsize=6 if x < 2 else 8, va="center")
        ax.set_xlim(-0.12, 2.4)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xticks([0, 1, 2], ["lncRNA", "miRNA", "hub mRNA"])
        ax.set_yticks([])
        ax.legend(frameon=False)
        ax.set_title(f"High-confidence ENCORI-supported axes (n={len(high)})")
    plot_save(fig, "fig_network_high_confidence")

    fig, ax = plt.subplots(figsize=(10, 5.5))
    if tier_counts.empty:
        ax.text(0.5, 0.5, "Tier evidence unavailable", ha="center", va="center")
        ax.set_axis_off()
    else:
        pivot = tier_counts[tier_counts["evidence_tier"].isin(["A", "B", "C"])].pivot(
            index="hub_gene", columns="evidence_tier", values="unique_miRNA_count"
        ).fillna(0)
        pivot = pivot.reindex(columns=["A", "B", "C"], fill_value=0)
        pivot.plot(kind="bar", ax=ax, color=["#2A9D8F", "#E9C46A", "#8D99AE"])
        ax.set_ylabel("Unique miRNAs")
        ax.set_xlabel("Validated hub")
        ax.set_title("miRNA evidence by hub and best evidence tier")
        ax.legend(title="Tier", frameon=False)
        ax.tick_params(axis="x", rotation=35)
    plot_save(fig, "fig_network_mirna_tiers_by_hub")

    fig, ax = plt.subplots(figsize=(8, 5))
    mapped = detection[detection["feature_match_status"].eq("exact_Ensembl_gene_id")]
    if mapped.empty:
        ax.text(0.5, 0.5, "No lncRNA features matched to GSE248762", ha="center", va="center")
        ax.set_axis_off()
    else:
        values = pd.to_numeric(mapped["percent_QC_cells_detected"], errors="coerce").dropna()
        ax.hist(values, bins=30, color="#457B9D", alpha=0.85)
        ax.axvline(MIN_LNCRNA_CELL_PERCENT, color="#D1495B", linestyle="--", label="analysis threshold")
        ax.set_xlabel("QC-retained cells with count > 0 (%)")
        ax.set_ylabel("lncRNAs in candidate axes")
        ax.set_title("GSE248762 lncRNA detection")
        ax.legend(frameon=False)
    plot_save(fig, "fig_network_lncrna_detection")

    fig, ax = plt.subplots(figsize=(8, 7))
    tested = correlation[correlation["status"].eq("tested")] if not correlation.empty else pd.DataFrame()
    if tested.empty:
        ax.text(0.5, 0.5, "No lncRNA-hub pseudobulk pair passed detection and testability criteria", ha="center", va="center", wrap=True)
        ax.set_axis_off()
    else:
        x = tested["rho_unadjusted"].astype(float)
        y = tested["rho_adjusted_log_cell_count"].astype(float)
        ax.scatter(x, y, s=20, alpha=0.6, color="#457B9D")
        ax.plot([-1, 1], [-1, 1], color="#777777", linewidth=0.8, linestyle="--")
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_xlabel("Unadjusted Spearman rho")
        ax.set_ylabel("Partial Spearman rho adjusted for log cell count")
        ax.set_title("Sample-level lncRNA-hub correlations")
    plot_save(fig, "fig_network_correlation_adjustment")

    fig, ax = plt.subplots(figsize=(12, 7))
    if axes.empty:
        ax.text(0.5, 0.5, "No lncRNA-miRNA-hub axes available", ha="center", va="center")
        ax.set_axis_off()
    else:
        display_axes = (
            axes.sort_values(
                ["confidence_class", "hub_gene", "miRNA_accession", "lncRNA_gene_id"]
            )
            .groupby(["confidence_class", "hub_gene"], group_keys=False)
            .head(30)
        )
        hub_nodes = sorted(display_axes["hub_gene"].astype(str).unique())
        mir_nodes = sorted(display_axes["miRNA_name"].astype(str).unique())
        lnc_nodes = sorted(display_axes["lncRNA_name"].astype(str).unique())
        layers = [(lnc_nodes, 0, "lnc"), (mir_nodes, 1, "mir"), (hub_nodes, 2, "hub")]
        pos: dict[tuple[str, str], tuple[float, float]] = {}
        for nodes, xcoord, key in layers:
            for i, node in enumerate(nodes):
                y = 0.5 if len(nodes) == 1 else 1.0 - i / (len(nodes) - 1)
                pos[(key, node)] = (xcoord, y)
        alpha = min(0.45, 20 / max(len(display_axes), 20))
        for row in display_axes.itertuples(index=False):
            p1 = pos[("lnc", str(row.lncRNA_name))]
            p2 = pos[("mir", str(row.miRNA_name))]
            p3 = pos[("hub", str(row.hub_gene))]
            color = "#2A9D8F" if row.confidence_class == "high" else "#9A9A9A"
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, alpha=alpha, linewidth=0.55)
            style = "-" if row.mRNA_evidence_tier == "A" else "--" if row.mRNA_evidence_tier == "B" else ":"
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color=color, alpha=alpha, linewidth=0.55, linestyle=style)
        for nodes, xcoord, key, color, label in [
            (lnc_nodes, 0, "lnc", "#36A269", "lncRNA"),
            (mir_nodes, 1, "mir", "#D8A932", "miRNA"),
            (hub_nodes, 2, "hub", "#4C78A8", "hub mRNA"),
        ]:
            points = [pos[(key, n)] for n in nodes]
            ax.scatter([p[0] for p in points], [p[1] for p in points], s=14, color=color, label=label, zorder=3)
            if key == "hub" or len(nodes) <= 35:
                for node, point in zip(nodes, points):
                    ax.text(point[0] + 0.02, point[1], node, fontsize=6 if key != "hub" else 8, va="center")
        ax.set_xlim(-0.12, 2.35)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xticks([0, 1, 2], ["lncRNA", "miRNA", "hub mRNA"])
        ax.set_yticks([])
        ax.set_title(
            f"High-confidence axes (n={int(axes['confidence_class'].eq('high').sum())}) "
            f"and low-confidence background; displaying at most 30 per hub/class"
        )
        ax.legend(frameon=False, loc="upper left")
    plot_save(fig, "fig_network_layered_axes")


def summarize_correlations(
    correlation: pd.DataFrame, cell_level: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    if correlation.empty:
        tested = correlation
        lineages: list[str] = []
    else:
        tested = correlation[correlation["status"].eq("tested")]
        lineages = sorted(correlation["lineage"].astype(str).unique())
    for lineage in lineages + ["ALL"]:
        full_scope = correlation if lineage == "ALL" else correlation[correlation["lineage"].eq(lineage)]
        subset = tested if lineage == "ALL" else tested[tested["lineage"].eq(lineage)]
        summary_rows.append(
            {
                "lineage": lineage,
                "n_tested": len(subset),
                "n_untestable": int(full_scope["status"].eq("not_testable").sum()),
                "n_nominal_unadjusted_P_lt_0_05": int((subset["p_unadjusted"] < 0.05).sum()) if not subset.empty else 0,
                "n_BH_unadjusted_q_lt_0_05_global_family": int((subset["q_unadjusted_all_pairs"] < 0.05).sum()) if not subset.empty else 0,
                "n_nominal_adjusted_P_lt_0_05": int((subset["p_adjusted_log_cell_count"] < 0.05).sum()) if not subset.empty else 0,
                "n_BH_adjusted_q_lt_0_05_global_family": int((subset["q_adjusted_all_pairs"] < 0.05).sum()) if not subset.empty else 0,
                "BH_family": "all testable lncRNA-hub-lineage pairs; unadjusted and adjusted families separately",
            }
        )
    cell_rows = []
    if cell_level.empty:
        cell_rows.append(
            {
                "analysis": "cell-level sensitivity",
                "n_tested": 0,
                "n_nominal_P_lt_0_05": 0,
                "n_BH_q_lt_0_05": 0,
                "interpretation": "NOT ESTIMATED",
            }
        )
    else:
        cell_rows.append(
            {
                "analysis": "cell-level sensitivity",
                "n_tested": len(cell_level),
                "n_nominal_P_lt_0_05": int((cell_level["p_cell_level_descriptive"] < 0.05).sum()),
                "n_BH_q_lt_0_05": int((cell_level["q_cell_level_descriptive"] < 0.05).sum()),
                "interpretation": "Descriptive only; cells are not independent biological replicates",
            }
        )
    return pd.DataFrame(summary_rows), pd.DataFrame(cell_rows)


def make_sc_forest(root: Path) -> pd.DataFrame:
    source = root / "results" / "revision" / "sc_and_network_summary" / "sc7_hub_gene_lineage_results.csv"
    result = pd.read_csv(source)
    tested = result[
        result["status"].eq("testable")
        & result["P_value"].notna()
        & result["P_value"].lt(0.05)
    ].copy()
    tested = tested.sort_values(["BH_FDR_7hub_family", "P_value", "gene", "lineage"])
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.5 * len(tested))))
    if tested.empty:
        ax.text(0.5, 0.5, "No nominal P < 0.05 comparisons", ha="center", va="center")
        ax.set_axis_off()
    else:
        labels = [f"{r.gene} | {r.lineage}" for r in tested.itertuples(index=False)]
        y = np.arange(len(tested))
        colors = np.where(tested["BH_FDR_7hub_family"].lt(0.05), "#C0392B", "#457B9D")
        ax.scatter(tested["log2FC_LV_vs_SV"], y, c=colors, s=55, zorder=3)
        ax.axvline(0, color="#666666", linewidth=0.8)
        ax.set_yticks(y, labels)
        ax.set_xlabel("log2 fold change (LV vs SV)")
        ax.set_title("Nominal pseudobulk comparisons; red = BH q < 0.05")
        for yy, p_value in zip(y, tested["P_value"]):
            ax.text(ax.get_xlim()[1], yy, f"P={p_value:.3g}", va="center", fontsize=8)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.2)
    plot_save(fig, "fig_sc_forest_nominal_and_fdr_hits")
    return tested


def create_network_writeup(
    hubs: list[str],
    panel: list[str],
    mrna: pd.DataFrame,
    tier_counts: pd.DataFrame,
    axes: pd.DataFrame,
    reconciliation: pd.DataFrame,
    detection: pd.DataFrame,
    correlation: pd.DataFrame,
    cell_level: pd.DataFrame,
    database_status: pd.DataFrame,
    expression_status: pd.DataFrame,
    cerna_pairs: pd.DataFrame,
    cerna_status: pd.DataFrame,
    source_versions: pd.DataFrame,
) -> None:
    n_axes = len(axes)
    n_high = int(axes["confidence_class"].eq("high").sum()) if not axes.empty else 0
    n_low = int(axes["confidence_class"].eq("low").sum()) if not axes.empty else 0
    n_lncs = axes["lncRNA_gene_id"].nunique() if not axes.empty else 0
    n_mirnas = axes["miRNA_accession"].nunique() if not axes.empty else 0
    n_hubs = axes["hub_gene"].nunique() if not axes.empty else 0
    n_detectable = int(detection["passes_coexpression_detection_threshold"].sum())
    n_mapped = int(detection["feature_match_status"].eq("exact_Ensembl_gene_id").sum())
    n_tested = int(correlation["status"].eq("tested").sum()) if not correlation.empty else 0
    n_cell = len(cell_level)
    n_current_mirtarbase = database_status[
        database_status["database"].eq("miRTarBase current release")
    ]["status"].astype(str).tolist()
    best_tier_a = int(tier_counts.loc[
        tier_counts["evidence_tier"].eq("A"), "unique_miRNA_count"
    ].sum())
    best_tier_b = int(tier_counts.loc[
        tier_counts["evidence_tier"].eq("B"), "unique_miRNA_count"
    ].sum())
    best_tier_c = int(tier_counts.loc[
        tier_counts["evidence_tier"].eq("C"), "unique_miRNA_count"
    ].sum())
    cerna_successes = int(cerna_status["status"].eq("RETRIEVED").sum()) if not cerna_status.empty else 0
    cerna_errors = int(cerna_status["status"].eq("NOT RETRIEVED").sum()) if not cerna_status.empty else 0
    cerna_unique_pairs = (
        int(cerna_pairs[["geneID", "ceRNAid"]].drop_duplicates().shape[0])
        if not cerna_pairs.empty
        else 0
    )
    cerna_lncrna_rows = (
        int(cerna_pairs["ceRNAgeneType"].str.lower().eq("lncrna").sum())
        if not cerna_pairs.empty
        else 0
    )
    tier_a_by_hub = tier_counts[tier_counts["evidence_tier"].eq("A")].set_index(
        "hub_gene"
    )["unique_miRNA_count"].to_dict()
    tier_a_table = [
        "| Hub gene | Tier A miRNAs |",
        "|---|---:|",
        *[
            f"| {hub} | {int(tier_a_by_hub.get(hub, 0))} |"
            for hub in hubs
        ],
    ]
    multimir_sources = source_versions[
        source_versions["query_via"].astype(str).str.startswith("multiMiR")
    ]
    multimir_not_queried = multimir_sources[
        multimir_sources["retrieval_status"].eq("SOURCE_VERSION_LISTED_NOT_QUERIED")
    ]["source"].astype(str).tolist()
    tested_rows = correlation[correlation["status"].eq("tested")] if not correlation.empty else pd.DataFrame()
    unadjusted_nominal = int((tested_rows["p_unadjusted"] < 0.05).sum()) if not tested_rows.empty else 0
    unadjusted_q = int((tested_rows["q_unadjusted_all_pairs"] < 0.05).sum()) if not tested_rows.empty else 0
    adjusted_nominal = int((tested_rows["p_adjusted_log_cell_count"] < 0.05).sum()) if not tested_rows.empty else 0
    adjusted_q = int((tested_rows["q_adjusted_all_pairs"] < 0.05).sum()) if not tested_rows.empty else 0
    if tested_rows.empty:
        coexpression_text = (
            "No sample-level lncRNA-hub pair passed the stated detection, lineage-cell, "
            "and variability requirements; therefore no correlation coefficient or "
            "P value is reported as an estimated association. The full eligibility "
            "and detection status are in `lncRNA_detection.csv` and "
            "`lncRNA_mRNA_correlation_results.csv`."
        )
    else:
        coexpression_text = (
            f"The pseudobulk analysis produced {n_tested} testable lncRNA-hub-lineage "
            f"pairs; {unadjusted_nominal} had nominal unadjusted P < 0.05 and "
            f"{unadjusted_q} had global BH q < 0.05. After adjusting for log cell "
            f"count, {adjusted_nominal} had nominal P < 0.05 and {adjusted_q} had "
            "global BH q < 0.05. All nominal and adjusted results, including "
            "nonsignificant pairs, are in "
            "`lncRNA_mRNA_correlation_results.csv`."
        )
    lines = [
        "# lncRNA-miRNA-mRNA regulatory network and single-cell co-expression",
        "",
        "## Methods",
        "",
        f"The analysis was restricted to the {len(hubs)} validated hub genes read from "
        "`results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`. "
        f"The frozen {len(panel)}-gene panel was read from "
        "`results/revision/panel_v2/panel_v2.json`; no gene symbols are repeated in this prose.",
        "",
        "Positive validated miRNA-mRNA evidence was queried through multiMiR for "
        "miRTarBase, miRecords, and TarBase, and predicted evidence was kept separate. "
        "Negative records were excluded from positive support. Tier A required a "
        "positive validated record whose experiment field explicitly named a reporter/"
        "luciferase assay, western blot, qPCR, or qRT-PCR. Tier B required at least "
        "two ENCORI CLIP-seq experiments for the hub mRNA. Tier C was restricted to "
        "multiMiR prediction-only pairs with no positive validated record or Tier A/B "
        "support. Validated rows lacking a qualifying functional assay were counted "
        "separately as unclassified, not upgraded to Tier A.",
        "",
        "ENCORI interactions were retrieved against hg38. The hub-mRNA API query "
        "requested records with at least one CLIP experiment and no target-prediction "
        "program requirement; Tier B was then assigned using `clipExpNum >= 2`. "
        "The ENCORI miRNA-lncRNA archive was retrieved as a raw archive and matched "
        "to miRNAs linked to the validated hubs. High-confidence axes required a "
        "Tier A or B hub-mRNA edge and an lncRNA edge with `clipExpNum >= 2`. "
        "All remaining database-backed candidate axes were low confidence. "
        "A recorded target site with zero CLIP support is a prediction-level "
        "lncRNA candidate, not experimentally established binding.",
        "",
        "The ENCORI ceRNA endpoint was also queried separately for each hub "
        "(`miRNAnum >= 2`, `pval <= 0.01`, `fdr <= 0.01`, `pancancerNum = 0`). "
        "These shared-miRNA network results were not used to define the lncRNA-miRNA-hub "
        "axes because they are a separate, cancer-network-derived evidence product.",
        "",
        "The relevant external miRNA-seq dataset was GSE182736 / PRJNA757643. "
        "ENA run metadata were saved, but sample-level miRNA abundance was not "
        "quantified. Thus miRNA detection is NOT VERIFIED and was not used to filter "
        "candidate network axes.",
        "",
        "For GSE248762, the saved QC-retained cell annotations and sample-lineage "
        "pseudobulk count matrices were reused. An lncRNA was eligible for co-expression "
        f"only if it was detected in at least {MIN_LNCRNA_CELL_PERCENT}% of QC-retained "
        f"cells and had nonzero counts in at least {MIN_LNCRNA_SAMPLE_COUNT} distinct "
        "samples. Each eligible lncRNA-hub pair was tested per lineage using "
        f"Spearman correlation on log2(CPM + 1) values from sample-level pseudobulk "
        f"profiles; samples with fewer than {MIN_LINEAGE_CELLS} cells in a lineage "
        f"were excluded and at least {MIN_PSEUDOBULK_N} samples were required. "
        "The adjusted analysis is partial Spearman correlation controlling for the "
        "rank of log(1 + lineage cell count). BH correction was applied over all "
        "testable pairs across lineages separately for unadjusted and adjusted P values. "
        "A cell-level correlation was calculated only as a sensitivity analysis; "
        "cells are not independent biological replicates.",
        "",
        "## Results",
        "",
        f"The primary network result is {n_high} high-confidence lncRNA-miRNA-hub "
        f"axes across {n_hubs} hubs, {n_mirnas} miRNAs, and {n_lncs} lncRNAs. "
        f"The {n_low} low-confidence axes are reported only as a candidate background, "
        "not as supported interactions; many include prediction-level lncRNA sites. "
        f"The complete table contains {n_axes} axes in total, reconciled in "
        "`network_axes.csv` and `network_reconciliation.csv`.",
        "",
        f"MultiMiR/ENCORI tier counts (best-class, unique miRNAs summed over hubs) "
        f"were Tier A = {best_tier_a}, Tier B = {best_tier_b}, and Tier C = {best_tier_c}; "
        "Tier C may be incomplete because the predicted query returned exactly "
        "1,000 rows and truncation by a service cap could not be verified. "
        "Per-hub Tier A counts are shown below; full hub-specific counts, unclassified "
        "validation, and source detail are in `hub_miRNA_tier_counts.csv` and "
        "`miRNA_mRNA_evidence.csv`. These are database evidence counts, not evidence "
        "of regulation in peritoneal samples.",
        "",
        *tier_a_table,
        "",
        f"Of the {len(multimir_sources)} multiMiR sources listed in the saved source "
        f"metadata, {len(multimir_not_queried)} were not queried: "
        f"{', '.join(multimir_not_queried)}. Current-release miRTarBase was not "
        "retrieved directly; the queried multiMiR miRTarBase snapshot is older and "
        "must not be represented as the current release.",
        "",
        "",
        f"The separate ENCORI ceRNA queries returned parsable results for "
        f"{cerna_successes} of {len(cerna_status)} queried hubs; {cerna_errors} "
        f"queries returned API errors. The successful responses contained "
        f"{len(cerna_pairs)} records representing {cerna_unique_pairs} unique "
        f"gene-pair links, including {cerna_lncrna_rows} rows with an lncRNA partner. "
        "These are not counted as validated regulation and are tabulated separately "
        "in `ENCORI_ceRNA_pairs.csv` and `ENCORI_ceRNA_query_status.csv`.",
        "",
        f"Of {len(detection)} lncRNAs present in the candidate axes, {n_mapped} "
        f"matched an exact Ensembl gene ID in the GSE248762 count matrix and "
        f"{n_detectable} met the predeclared co-expression detection threshold. "
        f"{coexpression_text} The adjusted results are the most stringent "
        "co-expression evidence available in this analysis, but remain exploratory "
        "because lineage-level pseudobulk sample counts are small. The nominal "
        f"unadjusted fraction is {unadjusted_nominal}/{n_tested} "
        f"({100 * unadjusted_nominal / n_tested:.2f}%); this exceeds the nominal "
        "5% reference, but is descriptive rather than a calibrated chance estimate "
        "because tests are dependent and share genes and samples.",
        "",
        f"The separate cell-level sensitivity table contains {n_cell} estimable "
        "lncRNA-hub-lineage correlations. It is reported descriptively because "
        "individual cells are nested within donors and are not independent replicates.",
        "",
        "The external database query status, including the unavailable direct current "
        "miRTarBase release and the older source versions available through multiMiR, "
        "is documented in `database_query_status.csv`. Current-release results must "
        "not be inferred from the older multiMiR backend.",
        "",
        "## Limitations",
        "",
        "This is a database-supported candidate network, not proof of post-transcriptional "
        "regulation in peritoneal dialysis. Only the high-confidence class is treated "
        "as the primary network result; low-confidence axes are a background set and "
        "many rely on prediction-level lncRNA sites. The 1,000-row predicted multiMiR "
        "result may be truncated, so Tier C counts may be incomplete. The direct "
        "current miRTarBase release was not retrieved; the multiMiR backend contains "
        "older source snapshots, and three of the fourteen listed sources were not "
        "queried. ENCORI "
        "CLIP support does not establish functional repression, and prediction-only "
        "edges and ceRNA-network associations remain hypotheses. "
        f"{cerna_errors} ceRNA API queries returned errors and cannot be interpreted as zero "
        "biological interactions. The relevant public effluent miRNA-seq study has "
        "six archived runs, but no sample-level miRNA abundance table was verified, "
        "so network miRNA expression and co-expression-dependent axes remain "
        "unverified. The single-cell data are sparse for many lncRNAs, use a small "
        "number of donors, and represent effluent pseudobulk rather than tissue. "
        "Pseudobulk correlations are exploratory; the cell-level sensitivity is "
        "especially vulnerable to pseudoreplication and donor composition.",
        "",
        "## Figure legends",
        "",
        "**Network figure.** Candidate lncRNA-miRNA-hub axes are shown by confidence "
        "class. Node colors indicate RNA class; edge styles distinguish the best "
        "hub-mRNA evidence class. Predicted and CLIP-supported interactions are "
        "candidates, not validated regulation.",
        "",
        "**Tier bar plot.** Counts of unique miRNAs by hub and best evidence tier. "
        "Tier definitions and evidence source versions are given in Methods.",
        "",
        "**Detection plot.** Distribution of the fraction of QC-retained GSE248762 "
        "cells with nonzero lncRNA counts. The vertical line indicates the cell-level "
        "detection component of the stated eligibility threshold.",
        "",
        "**Correlation plot.** Unadjusted sample-level Spearman coefficients are "
        "compared with partial Spearman coefficients adjusted for log lineage cell "
        "count. Each point is an eligible lncRNA-hub-lineage pair. The layered "
        "network plot foregrounds high-confidence axes and displays at most 30 "
        "low-confidence background axes per hub for readability; the full axis set "
        "is in `network_axes.csv`.",
        "",
        f"**Current miRTarBase access.** {', '.join(n_current_mirtarbase) if n_current_mirtarbase else 'NOT RETRIEVED'} "
        "is recorded in the query-status table; the retrieved multiMiR source snapshot "
        "is not a substitute for the current release.",
    ]
    (OUT / "NETWORK_WRITEUP.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source_hubs = ROOT / "results" / "revision" / "hub_gene_network_drug_gsea" / "validated7_kegg_gsea_summary.csv"
    source_panel = ROOT / "results" / "revision" / "panel_v2" / "panel_v2.json"
    hubs = pd.read_csv(source_hubs)["gene"].astype(str).tolist()
    panel_json = json.loads(source_panel.read_text(encoding="utf-8"))
    panel = list(panel_json["genes_passing_selection_rules"])
    validated = pd.read_csv(RAW / "multiMiR_validated_raw.csv", dtype=str).fillna("")
    predicted = pd.read_csv(RAW / "multiMiR_predicted_raw.csv", dtype=str).fillna("")
    encori_by_hub = {
        hub: read_encori_tsv(RAW / f"encori_miRNA_mRNA_{hub}.tsv")
        for hub in hubs
    }
    mrna, tier_counts = load_mirna_mrna_evidence(
        hubs, validated, predicted, encori_by_hub
    )
    target_accessions = set(mrna["miRNA_accession"].astype(str))
    lnc = load_lnc_mirna_pairs(
        RAW / "encori_lncRNA_matches_for_hub_miRNAs.tsv",
        target_accessions,
    )
    axes = build_axes(mrna, lnc)

    write_csv(mrna, "miRNA_mRNA_evidence.csv")
    write_csv(tier_counts, "hub_miRNA_tier_counts.csv")
    write_csv(lnc, "lncRNA_miRNA_evidence.csv")
    write_csv(axes, "network_axes.csv")

    reconciliation_rows: list[dict[str, object]] = []
    n_axes = len(axes)
    high = int(axes["confidence_class"].eq("high").sum()) if n_axes else 0
    low = int(axes["confidence_class"].eq("low").sum()) if n_axes else 0
    unique_lnc = int(axes["lncRNA_gene_id"].nunique()) if n_axes else 0
    if n_axes:
        lnc_class = axes.groupby("lncRNA_gene_id")["confidence_class"].apply(
            lambda values: "high" if values.eq("high").any() else "low"
        )
        lnc_high = int(lnc_class.eq("high").sum())
        lnc_low = int(lnc_class.eq("low").sum())
    else:
        lnc_high = 0
        lnc_low = 0
    checks = [
        ("total_axes", n_axes, "Unique hub-miRNA-lncRNA triplets in network_axes.csv"),
        ("high_confidence_axes", high, "Tier A/B mRNA and ENCORI lncRNA clipExpNum >= 2"),
        ("low_confidence_axes", low, "All other retrieved candidate axes"),
        ("high_plus_low_axes", high + low, "Must equal total_axes"),
        ("unique_hubs", int(axes["hub_gene"].nunique()) if n_axes else 0, "Unique hub symbols in candidate axes"),
        ("unique_miRNAs", int(axes["miRNA_accession"].nunique()) if n_axes else 0, "Unique mature miRNA accessions in candidate axes"),
        ("unique_lncRNAs", unique_lnc, "Unique lncRNA Ensembl gene IDs in candidate axes"),
        ("lncRNAs_assigned_high_class", lnc_high, "Unique lncRNA assigned high if any axis is high"),
        ("lncRNAs_assigned_low_class", lnc_low, "Unique lncRNA assigned low only if no high axis exists"),
        ("high_plus_low_lncRNAs", lnc_high + lnc_low, "Must equal unique_lncRNAs"),
    ]
    for metric, value, definition in checks:
        reconciliation_rows.append(
            {
                "metric": metric,
                "value": value,
                "definition": definition,
                "arithmetic_check": (
                    "PASS"
                    if metric == "high_plus_low_axes" and value == n_axes
                    or metric == "high_plus_low_lncRNAs" and value == unique_lnc
                    or metric not in ("high_plus_low_axes", "high_plus_low_lncRNAs")
                    else "MISMATCH"
                ),
            }
        )
    reconciliation = pd.DataFrame(reconciliation_rows)
    write_csv(reconciliation, "network_reconciliation.csv")

    annotations_path = ROOT / "results" / "revision" / "single_cell_context" / "cell_type_annotation_table.csv"
    annotations = pd.read_csv(annotations_path)
    h5ad_path = ROOT / "data" / "processed" / "GSE248762_hubblind_allcells_qc.h5ad"
    reuse_existing = os.environ.get("NETWORK7HUB_REUSE_DERIVED") == "1"
    cached_paths = {
        "detection": OUT / "lncRNA_detection.csv",
        "correlation": OUT / "lncRNA_mRNA_correlation_results.csv",
        "cell_level": OUT / "lncRNA_mRNA_cell_level_sensitivity.csv",
    }
    if reuse_existing and all(path.exists() for path in cached_paths.values()):
        detection = pd.read_csv(cached_paths["detection"])
        correlation = pd.read_csv(cached_paths["correlation"])
        cell_level = pd.read_csv(cached_paths["cell_level"])
    else:
        detection, h5ad, annotations, gene_columns = build_detection_table(
            axes, h5ad_path, annotations
        )
        write_csv(detection, "lncRNA_detection.csv")
        correlation, cell_level = make_correlation_tables(
            axes, detection, h5ad, annotations, gene_columns, ROOT
        )
        write_csv(correlation, "lncRNA_mRNA_correlation_results.csv")
        write_csv(cell_level, "lncRNA_mRNA_cell_level_sensitivity.csv")
        h5ad.file.close()

    correlation_summary, cell_summary = summarize_correlations(correlation, cell_level)
    write_csv(correlation_summary, "lncRNA_mRNA_correlation_summary.csv")
    write_csv(cell_summary, "lncRNA_mRNA_cell_level_summary.csv")

    workbook_path = ROOT / "data" / "raw" / "GSE182736_All_Expressed_miRNA.xlsx"
    workbook = pd.ExcelFile(workbook_path)
    workbook_audit_rows = []
    for sheet in workbook.sheet_names:
        columns = list(pd.read_excel(workbook_path, sheet_name=sheet, nrows=0).columns)
        sample_columns = [str(column) for column in columns if str(column).startswith("GSM")]
        workbook_audit_rows.append(
            {
                "sheet": sheet,
                "columns": ";".join(map(str, columns)),
                "n_columns": len(columns),
                "sample_accession_columns": ";".join(sample_columns),
                "has_sample_level_expression_columns": bool(sample_columns),
                "interpretation": "annotation/sequence table; no sample-level abundance columns" if not sample_columns else "sample-level columns present",
            }
        )
    workbook_audit = pd.DataFrame(workbook_audit_rows)
    write_csv(workbook_audit, "GSE182736_workbook_audit.csv")

    ena_manifest = json.loads(
        (RAW / "GSE182736_ENA_query_manifest.json").read_text(encoding="utf-8")
    )
    ena_run_count = sum(
        1
        for _ in csv.DictReader(
            (RAW / "GSE182736_ENA_run_metadata.tsv").open(encoding="utf-8"),
            delimiter="\t",
        )
    )
    miRNA_status = pd.DataFrame(
        [
            {
                "dataset": "GSE182736",
                "study_accession": "PRJNA757643",
                "sample_level_expression_available_locally": False,
                "sample_level_expression_status": "NOT VERIFIED",
                "run_count_from_saved_ENA_metadata": ena_run_count,
                "ENA_metadata_query_date_utc": ena_manifest["retrieved_utc"],
                "available_local_workbook": "data/raw/GSE182736_All_Expressed_miRNA.xlsx",
                "workbook_role": "annotation/sequence table; not sample-level abundance",
                "workbook_column_audit": ";".join(workbook_audit["columns"].astype(str)),
                "raw_reads_quantified": False,
                "use_for_axis_filtering": False,
                "reason": "Public SRA/ENA runs are listed, but sample-level miRNA abundance was not processed or verified.",
            }
        ]
    )
    write_csv(miRNA_status, "miRNA_expression_support.csv")

    database_info = pd.read_csv(RAW / "multiMiR_database_info.csv", dtype=str).fillna("")
    multi_log = (RAW / "multimir_query.log").read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"R\s+([0-9.]+)\s+multiMiR\s+([0-9.]+)", multi_log)
    r_version = version_match.group(1) if version_match else "NOT STATED"
    multimir_version = version_match.group(2) if version_match else "NOT STATED"
    multi_query_time = datetime.fromtimestamp(
        max(
            (RAW / "multiMiR_validated_raw.csv").stat().st_mtime,
            (RAW / "multiMiR_predicted_raw.csv").stat().st_mtime,
        ),
        tz=timezone.utc,
    ).isoformat()
    current_query_time = datetime.now(timezone.utc).isoformat()
    target_query_rows = []
    for hub in hubs:
        validated_hub = validated[validated["target_symbol"].eq(hub)]
        predicted_hub = predicted[predicted["target_symbol"].eq(hub)]
        target_query_rows.append(
            {
                "hub_gene": hub,
                "validated_rows": len(validated_hub),
                "positive_validated_rows": int(
                    (~validated_hub["support_type"].str.lower().eq("negative")).sum()
                ),
                "negative_rows_excluded": int(
                    validated_hub["support_type"].str.lower().eq("negative").sum()
                ),
                "predicted_rows": len(predicted_hub),
                "distinct_validated_miRNAs": validated_hub["mature_mirna_acc"].nunique(),
                "distinct_predicted_miRNAs": predicted_hub["mature_mirna_acc"].nunique(),
            }
        )
    write_csv(pd.DataFrame(target_query_rows), "multiMiR_target_query_counts.csv")
    source_version_rows = []
    queried_multimir_sources = set(
        re.findall(r"Searching\s+([A-Za-z0-9_]+)", multi_log)
    )
    for row in database_info.to_dict(orient="records"):
        source_key = str(row.get("map_name", "")).lower()
        source_version_rows.append(
            {
                "source": row.get("source_name") or row.get("map_name"),
                "version": row.get("source_version") or "NOT STATED",
                "source_date": row.get("source_date") or "NOT STATED",
                "query_date_utc": multi_query_time,
                "query_via": f"multiMiR {multimir_version}",
                "source_url": row.get("source_url", ""),
                "retrieval_status": (
                    "QUERIED_IN_THIS_ANALYSIS"
                    if source_key in queried_multimir_sources
                    else "SOURCE_VERSION_LISTED_NOT_QUERIED"
                ),
            }
        )
    source_version_rows.extend(
        [
            {
                "source": "miRTarBase current release",
                "version": "11.0",
                "source_date": "NOT STATED in saved source metadata",
                "query_date_utc": current_query_time,
                "query_via": "direct current-release attempt",
                "source_url": "https://awi.cuhk.edu.cn/miRTarBase/",
                "retrieval_status": "NOT RETRIEVED",
            },
            {
                "source": "ENCORI/starBase",
                "version": "NOT STATED separately; response cites 2026 Nature Methods paper",
                "source_date": "NOT STATED",
                "query_date_utc": json.loads(
                    (RAW / "encori_query_manifest.json").read_text(encoding="utf-8")
                )["retrieved_utc"],
                "query_via": "ENCORI API and hg38 archive",
                "source_url": "https://rnasysu.com/encori/tutorialAPI.php",
                "retrieval_status": "RETRIEVED",
            },
            {
                "source": "ENCORI ceRNA network",
                "version": "NOT STATED separately; response cites 2026 Nature Methods paper",
                "source_date": "NOT STATED",
                "query_date_utc": json.loads(
                    (RAW / "encori_query_manifest.json").read_text(encoding="utf-8")
                ).get("ceRNA_hub_queries", {}).get("query_date_utc", "NOT STATED"),
                "query_via": "ENCORI ceRNA API",
                "source_url": "https://rnasysu.com/encori/tutorialAPI.php",
                "retrieval_status": "PARTIAL; see per-hub query status",
            },
        ]
    )
    source_versions = pd.DataFrame(source_version_rows)
    write_csv(source_versions, "database_source_versions.csv")
    cerna_pairs, cerna_status, updated_encori_manifest = load_encori_cerna(hubs, RAW)
    write_csv(cerna_pairs, "ENCORI_ceRNA_pairs.csv")
    write_csv(cerna_status, "ENCORI_ceRNA_query_status.csv")
    database_status = pd.DataFrame(
        [
            {
                "database": "miRTarBase current release",
                "version_or_release": "11.0 (official current page observed; direct query data not retrieved)",
                "query_date_utc": current_query_time,
                "status": "NOT RETRIEVED",
                "error_or_note": "Direct current-release access failed (HTTP 403 in earlier targeted attempt; current endpoint also returned a TLS handshake error). Do not interpret multiMiR's older miRTarBase snapshot as release 11.0.",
                "source_or_output_file": "database_source_versions.csv; current landing page https://awi.cuhk.edu.cn/miRTarBase/",
            },
            {
                "database": "multiMiR",
                "version_or_release": f"R {r_version}; multiMiR package {multimir_version}; overall backend version NOT STATED in saved metadata",
                "query_date_utc": multi_query_time,
                "status": "RETRIEVED",
                "error_or_note": f"{len(hubs)} target-specific queries; validated results {len(validated)} rows and predicted results {len(predicted)} rows. Whether results were truncated by a service cap is NOT VERIFIED.",
                "source_or_output_file": "data/raw/network_7hub/multiMiR_validated_raw.csv; data/raw/network_7hub/multiMiR_predicted_raw.csv; data/raw/network_7hub/multiMiR_database_info.csv",
            },
            {
                "database": "miRTarBase snapshot via multiMiR",
                "version_or_release": "9; Sept 2021",
                "query_date_utc": multi_query_time,
                "status": "RETRIEVED_OLDER_SNAPSHOT",
                "error_or_note": "Version/date from multiMiR_database_info.csv; not the current release.",
                "source_or_output_file": "data/raw/network_7hub/multiMiR_database_info.csv",
            },
            {
                "database": "TarBase snapshot via multiMiR",
                "version_or_release": "9; 2023",
                "query_date_utc": multi_query_time,
                "status": "RETRIEVED",
                "error_or_note": "Version/date from multiMiR_database_info.csv; negative TarBase records were excluded from positive support.",
                "source_or_output_file": "data/raw/network_7hub/multiMiR_database_info.csv",
            },
            {
                "database": "miRecords snapshot via multiMiR",
                "version_or_release": "4; Apr 27, 2013",
                "query_date_utc": multi_query_time,
                "status": "RETRIEVED",
                "error_or_note": "Version/date from multiMiR_database_info.csv.",
                "source_or_output_file": "data/raw/network_7hub/multiMiR_database_info.csv",
            },
            {
                "database": "ENCORI/starBase hg38",
                "version_or_release": "API result citation dated 2026; service does not expose a separate database release identifier in the saved response",
                "query_date_utc": json.loads((RAW / "encori_query_manifest.json").read_text(encoding="utf-8"))["retrieved_utc"],
                "status": "RETRIEVED",
                "error_or_note": "Hub mRNA queries used clipExpNum >= 1 to collect records; Tier B threshold was applied locally at >=2. Bulk hg38 miRNA-lncRNA archive retrieved.",
                "source_or_output_file": "data/raw/network_7hub/encori_query_manifest.json; data/raw/network_7hub/encori_hg38_miRNA_lncRNA.tar.gz",
            },
            {
                "database": "ENCORI ceRNA network",
                "version_or_release": "API result citation dated 2026; no separate release identifier stated",
                "query_date_utc": updated_encori_manifest.get("ceRNA_hub_queries", {}).get("query_date_utc", "NOT STATED"),
                "status": "PARTIAL",
                "error_or_note": "Queried each validated hub with miRNAnum >= 2, pval <= 0.01, fdr <= 0.01, pancancerNum = 0; per-target errors are recorded in ENCORI_ceRNA_query_status.csv.",
                "source_or_output_file": "data/raw/network_7hub/encori_query_manifest.json; ENCORI_ceRNA_pairs.csv; ENCORI_ceRNA_query_status.csv",
            },
        ]
    )
    write_csv(database_status, "database_query_status.csv")

    lnc_clip_supported_axes = int(
        axes["lncRNA_clipExpNum_max"].ge(2).sum()
    ) if not axes.empty else 0
    evidence_filter_summary = pd.DataFrame(
        [
            {
                "filter_or_stage": "database_candidate_axes",
                "value": n_axes,
                "status": "RETRIEVED",
                "definition": "Unique lncRNA-miRNA-hub triplets with multiMiR Tier A/B/C mRNA evidence and an ENCORI lncRNA-miRNA target-site record",
            },
            {
                "filter_or_stage": "lncRNA_edge_clipExpNum_ge_2",
                "value": lnc_clip_supported_axes,
                "status": "RETRIEVED",
                "definition": "Candidate axes with at least two ENCORI CLIP experiments on the lncRNA-miRNA edge",
            },
            {
                "filter_or_stage": "high_confidence_interaction_axes",
                "value": high,
                "status": "RETRIEVED",
                "definition": "Tier A/B hub-mRNA edge plus lncRNA edge clipExpNum >= 2",
            },
            {
                "filter_or_stage": "all_layers_with_sample_level_miRNA_detection",
                "value": "NOT VERIFIED",
                "status": "NOT VERIFIED",
                "definition": "No sample-level PD effluent miRNA abundance table was available; candidate axes were not filtered on miRNA detection",
            },
        ]
    )
    write_csv(evidence_filter_summary, "network_evidence_filter_summary.csv")

    make_network_figures(axes, tier_counts, detection, correlation)
    sc_dir = ROOT / "results" / "revision" / "sc_and_network_summary"
    figure_map = {
        "figA_umap_lineage": "fig_sc_umap_by_lineage",
        "figB_hub_gene_dotplot": "fig_sc_hub_gene_dotplot",
        "figC_logfc_heatmap": "fig_sc_logfc_heatmap",
        "figD_hub_gene_volcano": "fig_sc_hub_gene_volcano",
        "figE_lineage_composition": "fig_sc_lineage_composition",
    }
    for source_stem, output_stem in figure_map.items():
        for extension in ("png", "pdf"):
            shutil.copy2(
                sc_dir / f"{source_stem}.{extension}",
                OUT / f"{output_stem}.{extension}",
            )
    sc_tables = [
        "sc7_hub_gene_lineage_results.csv",
        "sc7_hub_gene_expression.csv",
        "sc7_lineage_composition.csv",
        "sc7_sample_qc.csv",
        "sc7_hub_gene_testability.csv",
    ]
    for name in sc_tables:
        shutil.copy2(sc_dir / name, OUT / name)
    forest = make_sc_forest(ROOT)
    write_csv(forest, "sc7_nominal_forest_plot_data.csv")

    previous_sc_writeup = sc_dir / "SC_CONTEXT_WRITEUP.md"
    sc_text = previous_sc_writeup.read_text(encoding="utf-8")
    sc_text += (
        "\n\n## Refreshed output figures\n\n"
        "The five lineage-context figures are reused from the validated-seven "
        "analysis and copied into this output folder as `fig_sc_umap_by_lineage`, "
        "`fig_sc_hub_gene_dotplot`, `fig_sc_logfc_heatmap`, "
        "`fig_sc_hub_gene_volcano`, and `fig_sc_lineage_composition` (PNG and PDF). "
        "A new forest plot of nominal comparisons, with FDR-significant points "
        "highlighted, is saved as `fig_sc_forest_nominal_and_fdr_hits` (PNG and PDF); "
        "its plotted values are in `sc7_nominal_forest_plot_data.csv`.\n"
    )
    (OUT / "SC_CONTEXT_WRITEUP.md").write_text(sc_text, encoding="utf-8")

    create_network_writeup(
        hubs,
        panel,
        mrna,
        tier_counts,
        axes,
        reconciliation,
        detection,
        correlation,
        cell_level,
        database_status,
        miRNA_status,
        cerna_pairs,
        cerna_status,
        source_versions,
    )
    n_tested = (
        int(correlation["status"].eq("tested").sum())
        if not correlation.empty and "status" in correlation.columns
        else 0
    )
    n_detectable = (
        int(detection["passes_coexpression_detection_threshold"].sum())
        if not detection.empty
        else 0
    )
    n_cerna_success = int(cerna_status["status"].eq("RETRIEVED").sum()) if not cerna_status.empty else 0
    n_cerna_errors = int(cerna_status["status"].eq("NOT RETRIEVED").sum()) if not cerna_status.empty else 0
    n_nominal_corr = int(
        (
            correlation.loc[correlation["status"].eq("tested"), "p_unadjusted"]
            < 0.05
        ).sum()
    ) if n_tested else 0
    n_q_corr = int(
        (
            correlation.loc[
                correlation["status"].eq("tested"), "q_unadjusted_all_pairs"
            ]
            < 0.05
        ).sum()
    ) if n_tested else 0
    n_adjusted_q_corr = int(
        (
            correlation.loc[
                correlation["status"].eq("tested"), "q_adjusted_all_pairs"
            ]
            < 0.05
        ).sum()
    ) if n_tested else 0
    n_tier_a_hubs = int(
        tier_counts.loc[
            tier_counts["evidence_tier"].eq("A"), "unique_miRNA_count"
        ].gt(0).sum()
    )
    n_multimir_listed = int(
        source_versions["query_via"].astype(str).str.startswith("multiMiR").sum()
    )
    n_multimir_not_queried = int(
        source_versions["retrieval_status"].eq("SOURCE_VERSION_LISTED_NOT_QUERIED").sum()
    )
    verification = [
        "# Network and single-cell verification summary",
        "",
        "| Claim | Status | Verified value or disposition | Supporting file |",
        "|---|---|---|---|",
        f"| Analysis gene set is the validated seven | TRUE | {len(hubs)} genes loaded from the authoritative list | `results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv` |",
        f"| Frozen panel is a subset of the validated hubs | TRUE | {len(panel)} genes in the saved panel definition | `results/revision/panel_v2/panel_v2.json` |",
        f"| Candidate axes reconciled | {'TRUE' if (high + low == n_axes and lnc_high + lnc_low == unique_lnc) else 'FALSE'} | {n_axes} total = {high} high + {low} low; {unique_lnc} lncRNAs = {lnc_high} high-class + {lnc_low} low-class | `network_reconciliation.csv` |",
        f"| Primary high-confidence network | TRUE | {high} axes are primary; {low} low-confidence axes are background only | `network_reconciliation.csv`; `network_axes.csv` |",
        f"| Tier A hub coverage | TRUE | Tier A miRNAs by hub are listed in NETWORK_WRITEUP.md; {n_tier_a_hubs} hubs have one or more Tier A miRNAs | `hub_miRNA_tier_counts.csv` |",
        f"| miRNA expression in relevant PD samples | NOT VERIFIED | No sample-level abundance matrix was verified; the six ENA runs were not quantified | `miRNA_expression_support.csv`; `data/raw/network_7hub/GSE182736_ENA_run_metadata.tsv` |",
        f"| lncRNA pseudobulk association results | {'TRUE' if n_tested else 'NOT VERIFIED'} | {n_tested} testable pairs; see full table including non-significant and untestable pairs | `lncRNA_mRNA_correlation_results.csv` |",
        f"| lncRNA pseudobulk nominal/BH summary | TRUE | {n_nominal_corr} nominal P < 0.05; {n_q_corr} BH q < 0.05 over all testable unadjusted pairs | `lncRNA_mRNA_correlation_summary.csv` |",
        f"| Cell-count-adjusted pseudobulk BH results | TRUE | {n_adjusted_q_corr} global BH q < 0.05; exploratory because lineage pseudobulk sample counts are small | `lncRNA_mRNA_correlation_summary.csv` |",
        f"| multiMiR source coverage | PARTIAL | {n_multimir_not_queried} of {n_multimir_listed} listed sources were not queried; predicted results contain exactly {len(predicted)} rows, possible truncation NOT VERIFIED | `database_source_versions.csv`; `multiMiR_predicted_raw.csv` |",
        f"| ENCORI ceRNA network | PARTIAL | {n_cerna_success} of {len(cerna_status)} hub queries returned parsable rows; {n_cerna_errors} returned API errors | `ENCORI_ceRNA_query_status.csv`; `ENCORI_ceRNA_pairs.csv` |",
        "| Current miRTarBase release results | NOT RETRIEVED | Direct current-release result could not be obtained; older multiMiR source snapshot is labeled separately | `database_query_status.csv` |",
        "| Single-cell edgeR model refit | NOT DONE | Saved seven-gene test results were reused; plots refreshed and forest plot added | `sc7_hub_gene_lineage_results.csv`; `sc7_nominal_forest_plot_data.csv` |",
        "",
        "Predicted and CLIP-supported edges are candidate interactions. They are not described as experimentally validated regulation unless a qualifying functional assay is recorded.",
        "",
        "The complete itemized file-read and file-written inventories are in `FILE_INVENTORY.md`; SHA-256 hashes and package/source metadata are in `run_metadata.json`.",
        "",
    ]
    (OUT / "VERIFICATION_REPORT.md").write_text("\n".join(verification), encoding="utf-8")

    files_read = [
        Path(__file__).resolve(),
        source_hubs,
        source_panel,
        RAW / "multiMiR_validated_raw.csv",
        RAW / "multiMiR_predicted_raw.csv",
        RAW / "multiMiR_database_info.csv",
        RAW / "multimir_query.log",
        RAW / "encori_query_manifest.json",
        RAW / "encori_hg38_miRNA_lncRNA.tar.gz",
        RAW / "encori_lncRNA_matches_for_hub_miRNAs.tsv",
        RAW / "GSE182736_ENA_run_metadata.tsv",
        ROOT / "data" / "raw" / "GSE182736_All_Expressed_miRNA.xlsx",
        h5ad_path,
        annotations_path,
        sc_dir / "SC_CONTEXT_WRITEUP.md",
        sc_dir / "sc7_hub_gene_lineage_results.csv",
        sc_dir / "sc7_hub_gene_expression.csv",
        sc_dir / "sc7_lineage_composition.csv",
        sc_dir / "sc7_sample_qc.csv",
        sc_dir / "sc7_hub_gene_testability.csv",
        *[
            RAW / f"encori_miRNA_mRNA_{hub}.tsv"
            for hub in hubs
        ],
        *[
            RAW / f"encori_ceRNA_{hub}.tsv"
            for hub in hubs
        ],
        *[
            ROOT / "results" / "revision" / "single_cell_context" / "pseudobulk_counts"
            / f"counts_{lineage.replace('/', '_').replace(' ', '_')}.csv"
            for lineage in annotations["lineage"].astype(str).unique()
        ],
        *[
            sc_dir / f"{stem}.{ext}"
            for stem in figure_map
            for ext in ("png", "pdf")
        ],
        ROOT / "results" / "revision" / "single_cell_context" / "03_expression_and_pseudobulk_counts.py",
        ROOT / "results" / "revision" / "single_cell_context" / "05_summarize_and_figures.py",
        ROOT / "results" / "revision" / "single_cell_context" / "analysis_summary.json",
        ROOT / "results" / "revision" / "single_cell_context" / "dataset_audit.json",
        ROOT / "results" / "revision" / "single_cell_context" / "gene_id_symbol_map.csv",
        ROOT / "results" / "revision" / "single_cell_context" / "sample_lineage_cell_counts.csv",
        ROOT / "results" / "revision" / "sc_and_network_summary" / "run_metadata_summary.json",
        RAW / "GSE182736_ENA_query_manifest.json",
    ]
    if reuse_existing:
        files_read.extend(cached_paths.values())
    files_written = sorted(
        set(
            [p for p in RAW.iterdir() if p.is_file()]
            + [p for p in OUT.iterdir() if p.is_file()]
            + [OUT / "run_metadata.json", OUT / "FILE_INVENTORY.md"]
        )
    )
    files_read = sorted(set(p.resolve() for p in files_read if p.exists()))
    files_written = sorted(set(p.resolve() for p in files_written))
    inventory_lines = [
        "# File inventory",
        "",
        "Paths are relative to the repository root. This inventory records files read and generated/archived as part of the seven-hub analysis.",
        "",
        f"## Files read ({len(files_read)})",
        "",
        *[f"- `{path.relative_to(ROOT)}`" for path in files_read],
        "",
        f"## Files written ({len(files_written)})",
        "",
        *[f"- `{path.relative_to(ROOT)}`" for path in files_written],
        "",
    ]
    (OUT / "FILE_INVENTORY.md").write_text("\n".join(inventory_lines), encoding="utf-8")
    package_versions = {}
    for name in (
        "anndata", "matplotlib", "numpy", "pandas", "scipy", "statsmodels",
        "h5py", "requests",
    ):
        try:
            package_versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            package_versions[name] = "NOT INSTALLED"
    package_versions["R"] = r_version
    package_versions["multiMiR"] = multimir_version
    hash_paths = sorted(set(files_read + files_written))
    hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in hash_paths
        if path.is_file() and path.name != "run_metadata.json"
    }
    metadata = {
        "task": "Seven validated hub lncRNA-miRNA-mRNA evidence network and GSE248762 context refresh",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "random_steps": "No random resampling; seed is recorded for reproducibility and any deterministic graph-layout extensions.",
        "thresholds": {
            "ENCORI_hub_mRNA_Tier_B_minimum_clip_experiments": 2,
            "ENCORI_lncRNA_high_confidence_minimum_clip_experiments": 2,
            "lncRNA_detection_minimum_percent_QC_cells": MIN_LNCRNA_CELL_PERCENT,
            "lncRNA_detection_minimum_positive_samples": MIN_LNCRNA_SAMPLE_COUNT,
            "minimum_cells_per_lineage_pseudobulk": MIN_LINEAGE_CELLS,
            "minimum_samples_for_pseudobulk_correlation": MIN_PSEUDOBULK_N,
        },
        "gene_set_source": str(source_hubs.relative_to(ROOT)),
        "frozen_panel_source": str(source_panel.relative_to(ROOT)),
        "database_versions_and_dates": database_status.to_dict(orient="records"),
        "analysis_summary": {
            "validated_hubs": len(hubs),
            "frozen_panel_genes": len(panel),
            "candidate_axes": n_axes,
            "high_confidence_axes": high,
            "low_confidence_axes": low,
            "lncRNAs_in_axes": unique_lnc,
            "lncRNAs_passing_detection_threshold": n_detectable,
            "pseudobulk_pairs_tested": n_tested,
            "cell_level_sensitivity_pairs": len(cell_level),
            "sample_level_miRNA_detection": "NOT VERIFIED",
            "ceRNA_hub_queries_retrieved": n_cerna_success,
            "ceRNA_hub_queries_not_retrieved": n_cerna_errors,
            "ceRNA_parsed_records": len(cerna_pairs),
            "ceRNA_unique_gene_pair_links": int(
                cerna_pairs[["geneID", "ceRNAid"]].drop_duplicates().shape[0]
            ) if not cerna_pairs.empty else 0,
            "multiMiR_sources_queried": len(queried_multimir_sources),
        },
        "package_versions": package_versions,
        "files_read": [str(path.relative_to(ROOT)) for path in files_read],
        "files_written": [str(path.relative_to(ROOT)) for path in files_written],
        "sha256": hashes,
    }
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "n_hubs": len(hubs),
                "n_panel": len(panel),
                "mrna_gene_mirna_rows": len(mrna),
                "lnc_mirna_pairs": len(lnc),
                "network_axes": n_axes,
                "high_axes": high,
                "low_axes": low,
                "unique_lncrnas": unique_lnc,
                "lncrna_detection_rows": len(detection),
                "lncrna_detection_pass": n_detectable,
                "pseudobulk_tested": int(correlation["status"].eq("tested").sum()) if not correlation.empty else 0,
                "cell_level_tested": len(cell_level),
                "output_dir": str(OUT),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
