"""Audit fixed-panel detectability and derive an exploratory panel v2."""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
REVISION = ROOT / "results" / "revision"
OUTPUT = REVISION / "panel_v2"
GENES = ["TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "ECM1", "SERPINA10", "CST7"]
DISCOVERY_FILE = ROOT / "results" / "tables" / "GSE125498_all_results.csv"
BOOTSTRAPS = 10_000
NULL_SETS = 1_000
SEED = 42


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required input file(s) missing:\n" + "\n".join(missing))


def read_series_matrix(path: Path) -> pd.DataFrame:
    lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as source:
        for line in source:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                lines.append(line)
    if not lines:
        raise ValueError(f"No series-matrix expression table found in {path}")
    table = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    table.columns = [str(column).strip('"') for column in table.columns]
    table.iloc[:, 0] = table.iloc[:, 0].astype(str).str.strip('"')
    table = table.set_index(table.columns[0])
    table = table.apply(lambda column: pd.to_numeric(column.str.strip('"'), errors="raise"))
    if table.index.duplicated().any() or table.columns.duplicated().any():
        raise ValueError(f"Duplicate probe or sample IDs in {path}")
    return table


def is_true(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def as_float(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if np.isfinite(result) else float("nan")


def load_discovery() -> tuple[pd.DataFrame, pd.DataFrame]:
    discovery = pd.read_csv(DISCOVERY_FILE)
    if "Gene" not in discovery or "adj.P.Val" not in discovery:
        raise ValueError(f"Unexpected discovery table schema: {DISCOVERY_FILE}")
    if len(discovery) != 11_741 or discovery["Gene"].nunique() != 11_741:
        raise ValueError(
            "Expected 11,741 unique genes in the discovery universe; "
            f"found {len(discovery)} rows and {discovery['Gene'].nunique()} unique."
        )
    matches = discovery[discovery["Gene"].isin(GENES)].copy()
    if matches["Gene"].duplicated().any() or set(matches["Gene"]) != set(GENES):
        raise ValueError("The discovery table does not contain exactly one row per panel gene.")
    matches["BH_q_lt_0_05"] = matches["adj.P.Val"] < 0.05
    evidence = matches.rename(
        columns={
            "ID": "probe_id",
            "logFC": "logFC",
            "P.Value": "nominal_P",
            "adj.P.Val": "BH_adjusted_P_11741_gene_universe",
        }
    )[
        [
            "Gene",
            "probe_id",
            "logFC",
            "nominal_P",
            "BH_adjusted_P_11741_gene_universe",
            "BH_q_lt_0_05",
        ]
    ].rename(columns={"Gene": "gene"})
    return discovery, evidence


def load_gse62928() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = REVISION / "gse62928_validation_8genes"
    probe_info = pd.read_csv(base / "probe_selection_and_detection.csv")
    sample_groups = pd.read_csv(base / "sample_group_mapping.csv")
    matrix = read_series_matrix(base / "GSE62928_series_matrix.txt.gz")
    if set(sample_groups["sample_id"]) != set(matrix.columns):
        raise ValueError("GSE62928 sample mapping does not match series-matrix columns.")
    primary_probe = probe_info.set_index("Gene")["Primary_Highest_Mean_Probe"].to_dict()
    if set(primary_probe) != set(GENES):
        raise ValueError("GSE62928 primary-probe table does not cover the fixed panel.")
    rows: list[dict[str, Any]] = []
    expression_rows: list[pd.Series] = []
    for gene in GENES:
        probe = str(primary_probe[gene])
        if probe not in matrix.index:
            raise ValueError(f"GSE62928 primary probe {probe} for {gene} is missing.")
        values = matrix.loc[probe].astype(float)
        mean_expression = float(values.mean())
        samples_above = int((values > 5.0).sum())
        group_counts = {
            group: int((values[sample_groups.set_index("sample_id")["primary_group"] == group] > 5.0).sum())
            for group in ("EPS", "Non-EPS")
        }
        detected = mean_expression > 5.0
        rows.append(
            {
                "gene": gene,
                "cohort": "GSE62928",
                "contrast": "EPS vs combined non-EPS",
                "detection_rule": "highest-mean probe mean log2 expression > 5.0 across all 8 samples",
                "detected": detected,
                "annotation_present": True,
                "removal_filter": "none" if detected else "below expression cutoff",
                "expression_measure": "primary-probe RMA log2 expression",
                "mean_expression_or_count": mean_expression,
                "cutoff_value": 5.0,
                "samples_above_cutoff_n": samples_above,
                "samples_total_n": int(len(values)),
                "group1": "EPS",
                "group1_samples_above_cutoff_n": group_counts["EPS"],
                "group1_samples_n": 4,
                "group2": "Non-EPS",
                "group2_samples_above_cutoff_n": group_counts["Non-EPS"],
                "group2_samples_n": 4,
                "min_cell_filter_removed_samples_n": 0,
                "min_cell_rule": "not applicable to bulk biopsy microarray",
                "probe_or_gene_id": probe,
            }
        )
        expression_rows.append(values.rename(gene))
    expression = pd.DataFrame(expression_rows)
    expression.columns.name = "sample_id"
    return pd.DataFrame(rows), expression, sample_groups


def load_effluent_cohort(accession: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    analysis = REVISION / "effluent_validation" / accession / "analysis_20261009"
    counts = pd.read_csv(analysis / "all_cell_pseudobulk_counts.csv", index_col=0)
    counts.index = counts.index.astype(str)
    metadata = pd.read_csv(analysis / "all_cell_pseudobulk_metadata.csv")
    symbols = pd.read_csv(analysis / "gene_id_to_symbol.csv", dtype=str)
    primary = pd.read_csv(analysis / "primary_per_gene_results.csv")
    detection = pd.read_csv(analysis / "pseudobulk_gene_detection.csv")
    if metadata["sample_id"].duplicated().any() or set(metadata["sample_id"]) != set(counts.columns):
        raise ValueError(f"{accession}: sample metadata does not match raw pseudobulk columns.")
    if symbols["gene_id"].duplicated().any():
        raise ValueError(f"{accession}: duplicate gene IDs in symbol annotation.")
    if set(primary["gene"]) != set(GENES) or primary["gene"].duplicated().any():
        raise ValueError(f"{accession}: primary panel result does not cover the fixed panel.")
    if not set(GENES).issubset(set(detection["gene_symbol"])):
        raise ValueError(f"{accession}: detection table does not cover the fixed panel.")
    symbol_map = symbols.set_index("gene_id")["gene_symbol"].to_dict()
    gene_ids = primary.set_index("gene")["gene_id"].astype(str).to_dict()
    annotation_ids_by_symbol: dict[str, list[str]] = {}
    for gene_id, symbol in symbol_map.items():
        annotation_ids_by_symbol.setdefault(str(symbol), []).append(str(gene_id))
    groups = ("LV", "SV")
    contrast_meta = metadata[metadata["primary_group"].isin(groups)].copy()
    group_sample_ids = {
        group: contrast_meta.loc[contrast_meta["primary_group"] == group, "sample_id"].astype(str).tolist()
        for group in groups
    }
    if not group_sample_ids["LV"] or not group_sample_ids["SV"]:
        raise ValueError(f"{accession}: primary contrast sample groups are empty.")
    min_cells = 20
    bad_cell_samples = contrast_meta.loc[contrast_meta["n_qc_cells"] < min_cells, "sample_id"].tolist()
    if bad_cell_samples:
        raise ValueError(
            f"{accession}: contrast pseudobulk sample(s) below frozen {min_cells}-cell rule: "
            f"{bad_cell_samples}"
        )

    result_rows: list[dict[str, Any]] = []
    expression_rows: list[pd.Series] = []
    primary_index = primary.set_index("gene")
    detection_index = detection.set_index("gene_symbol")
    logcpm = pd.read_csv(analysis / "edger_logCPM.csv", index_col=0)
    for gene in GENES:
        gene_id = gene_ids[gene]
        annotation_present = gene_id in symbol_map and symbol_map[gene_id] == gene
        if gene_id not in counts.index:
            raise ValueError(f"{accession}: raw count row for annotated {gene} ({gene_id}) is missing.")
        vector = counts.loc[gene_id, contrast_meta["sample_id"]].astype(float)
        sample_sum = float(vector.sum())
        mean_count = float(vector.mean())
        counts_ge10 = {
            group: int((counts.loc[gene_id, group_sample_ids[group]] >= 10).sum())
            for group in groups
        }
        required = {group: math.ceil(len(group_sample_ids[group]) / 2) for group in groups}
        detected = all(counts_ge10[group] >= required[group] for group in groups)
        if not annotation_present:
            removal_filter = "not in annotation"
        elif sample_sum == 0:
            removal_filter = "zero counts"
        elif detected:
            removal_filter = "none"
        elif max((counts.loc[gene_id, contrast_meta["sample_id"]] >= 10).sum(), 0) == 0:
            removal_filter = "below expression cutoff"
        else:
            removal_filter = "filtered by min-sample rule"

        saved_detection = is_true(detection_index.loc[gene, "detected"])
        if saved_detection != detected:
            raise ValueError(f"{accession}: recomputed frozen detection differs for {gene}.")

        result_rows.append(
            {
                "gene": gene,
                "cohort": accession,
                "contrast": "LV vs SV",
                "detection_rule": (
                    "raw pseudobulk count >= 10 in at least half (ceiling) "
                    "of samples in both LV and SV groups"
                ),
                "detected": detected,
                "annotation_present": annotation_present,
                "removal_filter": removal_filter,
                "expression_measure": "raw pseudobulk counts per sample",
                "mean_expression_or_count": mean_count,
                "cutoff_value": 10,
                "samples_above_cutoff_n": int((vector >= 10).sum()),
                "samples_total_n": len(vector),
                "group1": "LV",
                "group1_samples_above_cutoff_n": counts_ge10["LV"],
                "group1_samples_n": len(group_sample_ids["LV"]),
                "group2": "SV",
                "group2_samples_above_cutoff_n": counts_ge10["SV"],
                "group2_samples_n": len(group_sample_ids["SV"]),
                "group1_required_samples_n": required["LV"],
                "group2_required_samples_n": required["SV"],
                "min_cell_filter_removed_samples_n": 0,
                "min_cell_rule": f"all contrast samples had >= {min_cells} QC cells",
                "probe_or_gene_id": gene_id,
                "raw_counts_sum": sample_sum,
                "saved_mean_logCPM": as_float(primary_index.loc[gene, "mean_logCPM"]),
            }
        )
        expression_rows.append(
            pd.Series(
                logcpm.loc[gene_id, metadata["sample_id"]]
                .astype(float)
                .to_numpy(),
                index=metadata["sample_id"].astype(str),
                name=gene,
            )
        )
    expression = pd.DataFrame(expression_rows)
    expression.columns.name = "sample_id"
    return pd.DataFrame(result_rows), expression, metadata, primary


def build_discovery_evidence(
    discovery: pd.DataFrame, output: pd.DataFrame
) -> pd.DataFrame:
    # Output is already selected from the full 11,741-gene discovery table.
    return output.copy()


def make_rule_table(
    audit: pd.DataFrame,
    discovery_evidence: pd.DataFrame,
    cohort_results: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    audit_index = audit.set_index(["gene", "cohort"])
    discovery_index = discovery_evidence.set_index("gene")
    row_by_cohort = {
        accession: result.set_index("gene")
        for accession, result in cohort_results.items()
    }
    rows: list[dict[str, Any]] = []
    for gene in GENES:
        per_cohort_detected = {
            accession: bool(audit_index.loc[(gene, accession), "detected"])
            for accession in ("GSE62928", "GSE248762", "GSE130888")
        }
        n_detected = sum(per_cohort_detected.values())
        rule_a = n_detected >= 2
        concordance: dict[str, Any] = {}
        detected_cohort_concordances: list[bool] = []
        for accession in ("GSE62928", "GSE248762", "GSE130888"):
            result = row_by_cohort[accession].loc[gene]
            if per_cohort_detected[accession]:
                concordant = is_true(
                    result.get(
                        "direction_concordant",
                        result.get("Direction_Concordant", False),
                    )
                )
                logfc = as_float(
                    result.get("logFC", result.get("Validation_log2FC", np.nan))
                )
                detected_cohort_concordances.append(concordant)
            else:
                concordant = None
                logfc = as_float(
                    result.get("logFC", result.get("Validation_log2FC", np.nan))
                )
            concordance[accession] = (concordant, logfc)
        rule_b = all(detected_cohort_concordances)
        rule_b_status = (
            "NOT APPLICABLE (no detected cohort)"
            if not detected_cohort_concordances
            else "PASS"
            if rule_b
            else "FAIL"
        )
        bh_flag = bool(discovery_index.loc[gene, "BH_q_lt_0_05"])
        rows.append(
            {
                "gene": gene,
                "GSE62928_detected": per_cohort_detected["GSE62928"],
                "GSE248762_detected": per_cohort_detected["GSE248762"],
                "GSE130888_detected": per_cohort_detected["GSE130888"],
                "detected_validation_cohorts_n": n_detected,
                "rule_a_detected_in_at_least_2_of_3": "PASS" if rule_a else "FAIL",
                "rule_b_detected_cohorts_n": len(detected_cohort_concordances),
                "GSE62928_direction_concordant_if_detected": concordance["GSE62928"][0],
                "GSE62928_validation_logFC": concordance["GSE62928"][1],
                "GSE248762_direction_concordant_if_detected": concordance["GSE248762"][0],
                "GSE248762_validation_logFC": concordance["GSE248762"][1],
                "GSE130888_direction_concordant_if_detected": concordance["GSE130888"][0],
                "GSE130888_validation_logFC": concordance["GSE130888"][1],
                "rule_b_concordant_in_every_detected_cohort": rule_b_status,
                "discovery_BH_q_lt_0_05_flag": bh_flag,
                "rule_c_status": "FLAG ONLY; NOT A HARD FILTER",
                "passes_panel_v2_selection": bool(rule_a and rule_b),
            }
        )
    return pd.DataFrame(rows)


def stratified_score(
    expression: pd.DataFrame,
    genes: list[str],
    metadata: pd.DataFrame,
    positive_group: str,
    negative_group: str,
    seed: int,
) -> tuple[dict[str, Any], pd.Series]:
    sample_ids = metadata["sample_id"].astype(str).tolist()
    groups = metadata.set_index("sample_id").loc[sample_ids, "primary_group"]
    keep = groups.isin([positive_group, negative_group]).to_numpy()
    selected = expression.loc[genes, sample_ids].to_numpy(dtype=float)[:, keep]
    selected_groups = groups.iloc[np.flatnonzero(keep)].to_numpy()
    labels = selected_groups == positive_group
    if selected.shape[1] != len(selected_groups) or not labels.any() or labels.all():
        raise ValueError("Cannot compute v2 score: contrast is missing a group.")
    means = selected.mean(axis=1, keepdims=True)
    stds = selected.std(axis=1, ddof=0, keepdims=True)
    if np.any(stds == 0) or np.any(~np.isfinite(stds)):
        raise ValueError("Cannot compute v2 score: at least one gene has no variation.")
    scores = ((selected - means) / stds).mean(axis=0)
    auc = float(roc_auc_score(labels.astype(int), scores))
    positive = scores[labels]
    negative = scores[~labels]
    mw_p = float(mannwhitneyu(positive, negative, alternative="two-sided").pvalue)
    rng = np.random.default_rng(seed)
    bootstrap_auc = np.empty(BOOTSTRAPS)
    for index in range(BOOTSTRAPS):
        positive_sample = rng.choice(positive, size=len(positive), replace=True)
        negative_sample = rng.choice(negative, size=len(negative), replace=True)
        bootstrap_auc[index] = roc_auc_score(
            np.r_[np.zeros(len(negative_sample)), np.ones(len(positive_sample))],
            np.r_[negative_sample, positive_sample],
        )
    complete_scores = pd.Series(np.nan, index=sample_ids, dtype=float)
    complete_scores.loc[np.asarray(sample_ids)[keep]] = scores
    return (
        {
            "n_genes": len(genes),
            "auc": auc,
            "auc_ci95_lower": float(np.quantile(bootstrap_auc, 0.025)),
            "auc_ci95_upper": float(np.quantile(bootstrap_auc, 0.975)),
            "mann_whitney_p_two_sided": mw_p,
            "n_positive": int(labels.sum()),
            "n_negative": int((~labels).sum()),
            "mean_score_positive": float(positive.mean()),
            "mean_score_negative": float(negative.mean()),
        },
        complete_scores,
    )


def matched_null_p(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    eligible_gene_ids: list[str],
    id_to_symbol: dict[str, str],
    primary: pd.DataFrame,
    positive_group: str,
    negative_group: str,
    observed_auc: float,
    seed: int,
) -> tuple[float, int, str]:
    logcpm_path = expression
    contrast_samples = metadata.loc[
        metadata["primary_group"].isin([positive_group, negative_group]), "sample_id"
    ].astype(str).tolist()
    rows = []
    for gene_id in eligible_gene_ids:
        if gene_id not in logcpm_path.index:
            continue
        symbol = id_to_symbol.get(gene_id, "")
        if not symbol:
            continue
        values = logcpm_path.loc[gene_id, contrast_samples].to_numpy(dtype=float)
        rows.append({"gene_id": gene_id, "symbol": symbol, "mean_expression": values.mean()})
    eligible = pd.DataFrame(rows).drop_duplicates("gene_id")
    if len(eligible) < len(genes) + 1:
        return float("nan"), 0, "not enough detected, annotated null genes"
    eligible = eligible.sort_values(["mean_expression", "gene_id"], kind="stable").reset_index(drop=True)
    eligible["decile"] = np.minimum(
        9, np.floor(np.arange(len(eligible)) * 10 / len(eligible)).astype(int)
    )
    target_deciles: list[int] = []
    symbols_to_ids = {
        str(row["symbol"]): str(row["gene_id"]) for row in rows
    }
    if any(gene not in symbols_to_ids for gene in genes):
        return float("nan"), 0, "target gene absent from matched-null expression universe"
    eligible_index = eligible.set_index("gene_id")
    target_ids = [symbols_to_ids[gene] for gene in genes]
    target_deciles = [int(eligible_index.loc[gene_id, "decile"]) for gene_id in target_ids]
    excluded_ids = set(target_ids)
    pools = {
        decile: subset.loc[~subset["gene_id"].isin(excluded_ids), "gene_id"].astype(str).tolist()
        for decile, subset in eligible.groupby("decile")
    }
    multiplicity = pd.Series(target_deciles).value_counts().to_dict()
    capacity = math.prod(
        math.comb(len(pools.get(int(decile), [])), int(required))
        if len(pools.get(int(decile), [])) >= int(required)
        else 0
        for decile, required in multiplicity.items()
    )
    if capacity < NULL_SETS:
        return float("nan"), 0, f"only {capacity} unique decile-matched sets available"
    rng = np.random.default_rng(seed)
    null_auc: list[float] = []
    seen: set[tuple[str, ...]] = set()
    max_attempts = NULL_SETS * 100_000
    attempts = 0
    groups = metadata.set_index("sample_id").loc[contrast_samples, "primary_group"].to_numpy()
    labels = groups == positive_group
    while len(null_auc) < NULL_SETS and attempts < max_attempts:
        attempts += 1
        selected_ids: list[str] = []
        for decile, needed in multiplicity.items():
            selected_ids.extend(
                rng.choice(pools[int(decile)], size=int(needed), replace=False).astype(str).tolist()
            )
        signature = tuple(sorted(selected_ids))
        if signature in seen:
            continue
        seen.add(signature)
        matrix = logcpm_path.loc[selected_ids, contrast_samples].to_numpy(dtype=float)
        gene_std = matrix.std(axis=1, ddof=0, keepdims=True)
        if np.any(gene_std == 0):
            continue
        zscore = (matrix - matrix.mean(axis=1, keepdims=True)) / gene_std
        null_auc.append(float(roc_auc_score(labels.astype(int), zscore.mean(axis=0))))
    if len(null_auc) != NULL_SETS:
        return float("nan"), len(null_auc), "could not draw 1,000 unique matched sets"
    p_value = float((1 + np.count_nonzero(np.asarray(null_auc) >= observed_auc)) / (NULL_SETS + 1))
    return p_value, len(null_auc), "1,000 detected-gene expression-decile matched sets"


def load_and_score_effluent(
    accession: str, genes: list[str]
) -> tuple[dict[str, Any], pd.DataFrame]:
    analysis = REVISION / "effluent_validation" / accession / "analysis_20261009"
    metadata = pd.read_csv(analysis / "all_cell_pseudobulk_metadata.csv")
    primary_all = pd.read_csv(analysis / "edger_all_gene_results.csv")
    primary = primary_all[primary_all["contrast"] == "LV_vs_SV"].copy()
    logcpm = pd.read_csv(analysis / "edger_logCPM.csv", index_col=0)
    detection = pd.read_csv(analysis / "pseudobulk_gene_detection.csv")
    symbols = pd.read_csv(analysis / "gene_id_to_symbol.csv", dtype=str)
    gene_id_for_symbol = (
        pd.read_csv(analysis / "primary_per_gene_results.csv")
        .set_index("gene")["gene_id"]
        .astype(str)
        .to_dict()
    )
    gene_ids = [gene_id_for_symbol[gene] for gene in genes]
    expression = logcpm.copy()
    score_metadata = metadata[metadata["primary_group"].isin(["LV", "SV"])].copy()
    score, sample_scores = stratified_score(
        expression,
        gene_ids,
        score_metadata,
        positive_group="LV",
        negative_group="SV",
        seed=SEED,
    )
    detected_ids = detection.loc[detection["detected"].map(is_true), "gene_id"].astype(str).tolist()
    id_to_symbol = symbols.set_index("gene_id")["gene_symbol"].to_dict()
    null_p, null_sets, null_note = matched_null_p(
        expression,
        metadata,
        genes,
        detected_ids,
        id_to_symbol,
        primary,
        "LV",
        "SV",
        score["auc"],
        SEED,
    )
    score.update(
        {
            "matched_null_auc_empirical_p": null_p,
            "matched_null_sets": null_sets,
            "matched_null_note": null_note,
            "analysis_label": "post hoc, exploratory, not preregistered",
        }
    )
    return score, sample_scores


def existing_score_summaries() -> dict[str, dict[str, Any]]:
    summaries = json.loads(
        (REVISION / "exploratory_posthoc_20261009" / "detected_gene_composites_summary.json").read_text(
            encoding="utf-8"
        )
    )
    return {entry["accession"]: entry for entry in summaries}


def original_verdict_rows() -> dict[str, dict[str, Any]]:
    biopsy = json.loads(
        (REVISION / "gse62928_validation_8genes" / "validation_summary.json").read_text(
            encoding="utf-8"
        )
    )
    effluent: dict[str, dict[str, Any]] = {}
    for accession in ("GSE248762", "GSE130888"):
        effluent[accession] = json.loads(
            (
                REVISION
                / "effluent_validation"
                / accession
                / "analysis_20261009"
                / "gene_set_summary.json"
            ).read_text(encoding="utf-8")
        )
    eight = next(item for item in biopsy["gene_set_results"] if item["Set"] == "Eight-gene panel")
    return {
        "GSE62928": {
            "verdict": "Not supported",
            "estimable": True,
            "n_genes": 8,
            "auc": eight["Composite_AUC"],
            "ci_lower": eight["Composite_AUC_Bootstrap_CI95_Lower"],
            "ci_upper": eight["Composite_AUC_Bootstrap_CI95_Upper"],
            "mann_whitney_p": eight["Composite_Mann_Whitney_Exact_Two_Sided_P"],
            "matched_null_p": eight["Null_AUC_Empirical_P"],
            "matched_null_note": "saved preregistered eight-gene result",
        },
        "GSE248762": {
            "verdict": effluent["GSE248762"]["preregistered_verdict"],
            "estimable": False,
            "n_genes": 8,
            "auc": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "mann_whitney_p": np.nan,
            "matched_null_p": np.nan,
            "matched_null_note": "full-panel composite not estimable under frozen rule",
        },
        "GSE130888": {
            "verdict": effluent["GSE130888"]["preregistered_verdict"],
            "estimable": False,
            "n_genes": 8,
            "auc": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "mann_whitney_p": np.nan,
            "matched_null_p": np.nan,
            "matched_null_note": "full-panel composite not estimable under frozen rule",
        },
    }


def write_outputs(
    audit: pd.DataFrame,
    discovery_evidence: pd.DataFrame,
    rules: pd.DataFrame,
    panel_genes: list[str],
    scores: dict[str, dict[str, Any]],
    score_rows: list[pd.DataFrame],
    original: dict[str, dict[str, Any]],
    existing_posthoc: dict[str, dict[str, Any]],
    current_date: str,
) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    audit.to_csv(OUTPUT / "panel_detectability_audit.csv", index=False)
    discovery_evidence.to_csv(OUTPUT / "panel_discovery_evidence.csv", index=False)
    rules.to_csv(OUTPUT / "panel_v2_rule_results.csv", index=False)

    panel_json = {
        "panel_name": "panel_v2",
        "date_defined": current_date,
        "status": "defined before any new test data",
        "defined_before_any_new_test_data": True,
        "genes_passing_selection_rules": panel_genes,
        "rule_text": {
            "a": "Detected in at least 2 of the 3 validation cohorts (GSE62928, GSE248762, GSE130888).",
            "b": "Concordant direction with discovery in every cohort where detected.",
            "c": "Record discovery BH q < 0.05 as a flag only; it is not a hard selection filter.",
        },
        "selection_logic": "A gene is included only when rules (a) AND (b) pass; rule (c) is recorded but does not filter.",
        "status_note": (
            f"This rule-derived panel is defined on {current_date} using existing saved cohort results. "
            "It is not preregistered; all v2 scores on these same existing cohorts are post hoc, "
            "exploratory, and not validation claims. The status means the panel is fixed before "
            "any future new test data."
        ),
    }
    (OUTPUT / "panel_v2.json").write_text(
        json.dumps(panel_json, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    sensitivity_rows: list[dict[str, Any]] = []
    for accession in ("GSE62928", "GSE248762", "GSE130888"):
        old = original[accession]
        prior_posthoc = existing_posthoc.get(accession, {})
        v2 = scores[accession]
        sensitivity_rows.append(
            {
                "cohort": accession,
                "original_8_gene_preregistered_verdict_unchanged": old["verdict"],
                "original_8_gene_composite_estimable": old["estimable"],
                "original_8_gene_n": old["n_genes"],
                "original_8_gene_AUC": old["auc"],
                "original_8_gene_AUC_CI95_lower": old["ci_lower"],
                "original_8_gene_AUC_CI95_upper": old["ci_upper"],
                "original_8_gene_Mann_Whitney_P": old["mann_whitney_p"],
                "original_8_gene_matched_null_P": old["matched_null_p"],
                "original_8_gene_matched_null_note": old["matched_null_note"],
                "existing_detected_gene_posthoc_AUC": prior_posthoc.get("composite_auc", np.nan),
                "existing_detected_gene_posthoc_CI95_lower": (
                    prior_posthoc.get("auc_bootstrap_ci_95", [np.nan, np.nan])[0]
                ),
                "existing_detected_gene_posthoc_CI95_upper": (
                    prior_posthoc.get("auc_bootstrap_ci_95", [np.nan, np.nan])[1]
                ),
                "existing_detected_gene_posthoc_Mann_Whitney_P": prior_posthoc.get(
                    "mann_whitney_p", np.nan
                ),
                "existing_detected_gene_posthoc_matched_null_P": prior_posthoc.get(
                    "matched_null_auc_empirical_p", np.nan
                ),
                "panel_v2_genes": ";".join(panel_genes),
                "panel_v2_score_label": "post hoc, exploratory, not preregistered",
                "panel_v2_n_genes": v2["n_genes"],
                "panel_v2_AUC": v2["auc"],
                "panel_v2_AUC_CI95_lower": v2["auc_ci95_lower"],
                "panel_v2_AUC_CI95_upper": v2["auc_ci95_upper"],
                "panel_v2_Mann_Whitney_P": v2["mann_whitney_p_two_sided"],
                "panel_v2_matched_null_P": v2["matched_null_auc_empirical_p"],
                "panel_v2_matched_null_n_sets": v2["matched_null_sets"],
                "panel_v2_matched_null_note": v2["matched_null_note"],
                "preregistered_verdict_changed": False,
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(OUTPUT / "panel_v2_score_sensitivity.csv", index=False)
    if score_rows:
        pd.concat(score_rows, ignore_index=True).to_csv(
            OUTPUT / "panel_v2_sample_scores.csv", index=False
        )

    counts_by_cohort = audit.groupby("cohort")["detected"].sum().astype(int).to_dict()
    report = [
        "# Panel detectability audit and exploratory panel v2",
        "",
        f"**Generated:** {current_date}",
        "",
        "## Detectability audit",
        "",
        "The audit reads saved cohort tables and the deposited GSE62928 series matrix; it does not rerun the preregistered analyses.",
        "",
        "| Cohort | Panel genes detected | Detection basis |",
        "|---|---:|---|",
        "| GSE62928 | "
        + f"{counts_by_cohort.get('GSE62928', 0)}/8 | Highest-mean primary probe mean log2 expression > 5 across 8 samples |",
        "| GSE248762 | "
        + f"{counts_by_cohort.get('GSE248762', 0)}/8 | Raw pseudobulk count ≥10 in at least half the samples in both LV and SV |",
        "| GSE130888 | "
        + f"{counts_by_cohort.get('GSE130888', 0)}/8 | Raw pseudobulk count ≥10 in at least half the samples in both LV and SV |",
        "",
        "### SERPINA10",
        "",
    ]
    serpina = audit[audit["gene"] == "SERPINA10"]
    for row in serpina.to_dict(orient="records"):
        report.append(
            f"- **{row['cohort']}:** detected={row['detected']}; filter={row['removal_filter']}; "
            f"annotation present={row['annotation_present']}; mean={row['mean_expression_or_count']:.6g} "
            f"({row['expression_measure']}); samples above cutoff={row['samples_above_cutoff_n']}/"
            f"{row['samples_total_n']}."
        )
    report.extend(
        [
            "",
            "The saved effluent annotations contain SERPINA10 in both cohorts. Its GSE248762 raw counts are nonzero in the saved pseudobulk matrix but remain below the count threshold; GSE130888's saved raw pseudobulk counts are all zero. Thus it is not missing from annotation, but its undetected status has cohort-specific causes.",
            "",
            "See `panel_detectability_audit.csv` for all genes, cohort-level means, sample counts above cutoff, and per-group counts. The biopsy cutoff is an operational reporting threshold, not a biological detection limit.",
            "",
        ]
    )
    (OUTPUT / "panel_detectability_audit.md").write_text("\n".join(report), encoding="utf-8")

    score_report = [
        "# Panel v2 rule and score report",
        "",
        f"**Panel date:** {current_date}. **Status:** defined before any new test data.",
        "",
        "Panel inclusion is determined mechanically by rule (a) AND rule (b). BH q<0.05 is only a flag (rule c), not a filter. All v2 scores below reuse existing cohorts and are labeled post hoc, exploratory, and not preregistered; they are not new validation claims.",
        "",
        "## Rule-derived genes",
        "",
        "| Gene | Detected in ≥2/3 | Direction concordant wherever detected | Discovery BH q<0.05 flag | Included |",
        "|---|---|---|---|---|",
    ]
    for row in rules.to_dict(orient="records"):
        score_report.append(
            f"| {row['gene']} | {row['rule_a_detected_in_at_least_2_of_3']} | "
            f"{row['rule_b_concordant_in_every_detected_cohort']} | "
            f"{'Yes' if row['discovery_BH_q_lt_0_05_flag'] else 'No'} | "
            f"{'Yes' if row['passes_panel_v2_selection'] else 'No'} |"
        )
    score_report.extend(
        [
            "",
            f"**Panel v2 ({len(panel_genes)} genes):** {', '.join(panel_genes) if panel_genes else 'none'}.",
            "",
            "## Original verdicts and v2 exploratory scores",
            "",
            "| Cohort | Original 8-gene verdict (unchanged) | Original 8-gene AUC (95% CI) | Original MW P | Original matched-null P | V2 exploratory AUC (95% CI) | V2 MW P | V2 matched-null P |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in sensitivity_rows:
        def fmt(value: Any) -> str:
            return "not estimable/not computed" if pd.isna(value) else f"{float(value):.4g}"

        old_auc = (
            "not estimable"
            if pd.isna(row["original_8_gene_AUC"])
            else (
                f"{fmt(row['original_8_gene_AUC'])} "
                f"({fmt(row['original_8_gene_AUC_CI95_lower'])}–"
                f"{fmt(row['original_8_gene_AUC_CI95_upper'])})"
            )
        )
        new_auc = (
            f"{fmt(row['panel_v2_AUC'])} "
            f"({fmt(row['panel_v2_AUC_CI95_lower'])}–{fmt(row['panel_v2_AUC_CI95_upper'])})"
        )
        score_report.append(
            f"| {row['cohort']} | {row['original_8_gene_preregistered_verdict_unchanged']} | "
            f"{old_auc} | {fmt(row['original_8_gene_Mann_Whitney_P'])} | "
            f"{fmt(row['original_8_gene_matched_null_P'])} | {new_auc} | "
            f"{fmt(row['panel_v2_Mann_Whitney_P'])} | {fmt(row['panel_v2_matched_null_P'])} |"
        )
    score_report.extend(
        [
            "",
            "Original effluent 8-gene composites were not estimable under their frozen detection rule. Existing detected-gene exploratory scores are separately retained in `panel_v2_score_sensitivity.csv`; they are not the original fixed-panel test. The historical GSE130888 detected-gene AUC of 0.8056 included normal controls; this audit uses the corrected saved result of 0.7083 for LV vs SV only. The biopsy v2 genes remain below the operational expression threshold, so its v2 matched-null P is not calculated. For effluent v2 score nulls, 1,000 detected-gene expression-decile matched sets and the plus-one empirical-P formula were used.",
            "",
            "The v2 score uses the mean of per-gene z-scored expression across the primary contrast samples; AUC confidence intervals use 10,000 stratified bootstrap resamples (seed 42). GSE130888 normal controls are excluded from the primary LV-vs-SV score and test.",
            "",
            "## Integrity",
            "",
            "No preregistered result, verdict, threshold, sample inclusion, or file under the previous cohort analysis folders was modified. The panel is rule-derived; no gene was selected by manual judgment.",
            "",
        ]
    )
    (OUTPUT / "panel_v2_report.md").write_text("\n".join(score_report), encoding="utf-8")


def main() -> None:
    required = [
        DISCOVERY_FILE,
        REVISION / "gse62928_validation_8genes" / "probe_selection_and_detection.csv",
        REVISION / "gse62928_validation_8genes" / "sample_group_mapping.csv",
        REVISION / "gse62928_validation_8genes" / "GSE62928_series_matrix.txt.gz",
        REVISION / "gse62928_validation_8genes" / "validation_summary.json",
        REVISION / "exploratory_posthoc_20261009" / "detected_gene_composites_summary.json",
    ]
    for accession in ("GSE248762", "GSE130888"):
        analysis = REVISION / "effluent_validation" / accession / "analysis_20261009"
        required.extend(
            analysis / filename
            for filename in (
                "all_cell_pseudobulk_counts.csv",
                "all_cell_pseudobulk_metadata.csv",
                "gene_id_to_symbol.csv",
                "primary_per_gene_results.csv",
                "pseudobulk_gene_detection.csv",
                "edger_logCPM.csv",
                "edger_all_gene_results.csv",
                "gene_set_summary.json",
            )
        )
    require_files(required)
    if OUTPUT.exists() and any(OUTPUT.iterdir()):
        raise FileExistsError(
            f"{OUTPUT} already contains files; refusing to overwrite prior panel_v2 outputs."
        )

    discovery, evidence = load_discovery()
    audit_frames: list[pd.DataFrame] = []
    gse62928_audit, gse62928_expression, gse62928_metadata = load_gse62928()
    audit_frames.append(gse62928_audit)
    cohort_results: dict[str, pd.DataFrame] = {
        "GSE62928": pd.read_csv(
            REVISION / "gse62928_validation_8genes" / "per_gene_primary_results.csv"
        ).rename(
            columns={
                "Gene": "gene",
                "Direction_Concordant": "direction_concordant",
                "Validation_log2FC": "logFC",
            }
        )
    }
    for accession in ("GSE248762", "GSE130888"):
        cohort_audit, _, _, cohort_primary = load_effluent_cohort(accession)
        audit_frames.append(cohort_audit)
        cohort_results[accession] = cohort_primary
    audit = pd.concat(audit_frames, ignore_index=True)
    if len(audit) != len(GENES) * 3:
        raise ValueError(f"Expected {len(GENES) * 3} detectability rows, got {len(audit)}.")

    rules = make_rule_table(audit, evidence, cohort_results)
    panel_genes = rules.loc[rules["passes_panel_v2_selection"], "gene"].tolist()

    scores: dict[str, dict[str, Any]] = {}
    sample_score_tables: list[pd.DataFrame] = []

    # Score using the saved RMA probe values; controls are the combined non-EPS group.
    gse62928_score, gse62928_score_values = stratified_score(
        gse62928_expression,
        panel_genes,
        gse62928_metadata,
        positive_group="EPS",
        negative_group="Non-EPS",
        seed=SEED,
    )
    gse62928_score.update(
        {
            "matched_null_auc_empirical_p": float("nan"),
            "matched_null_sets": 0,
            "matched_null_note": (
                "not computed: all v2 genes are below the operational mean-expression cutoff"
            ),
            "analysis_label": "post hoc, exploratory, not preregistered",
        }
    )
    scores["GSE62928"] = gse62928_score
    biopsy_score_table = gse62928_metadata[["sample_id", "primary_group"]].copy()
    biopsy_score_table["panel_v2_composite"] = gse62928_score_values.reindex(
        biopsy_score_table["sample_id"]
    ).to_numpy()
    biopsy_score_table["analysis_label"] = "post hoc, exploratory, not preregistered"
    biopsy_score_table["accession"] = "GSE62928"
    sample_score_tables.append(biopsy_score_table)

    for accession in ("GSE248762", "GSE130888"):
        score, sample_scores = load_and_score_effluent(accession, panel_genes)
        scores[accession] = score
        analysis = REVISION / "effluent_validation" / accession / "analysis_20261009"
        metadata = pd.read_csv(analysis / "all_cell_pseudobulk_metadata.csv")
        table = metadata[["sample_id", "primary_group"]].copy()
        table["panel_v2_composite"] = sample_scores.reindex(table["sample_id"]).to_numpy()
        table["analysis_label"] = "post hoc, exploratory, not preregistered"
        table["accession"] = accession
        sample_score_tables.append(table)

    existing_posthoc = existing_score_summaries()
    original = original_verdict_rows()
    write_outputs(
        audit,
        build_discovery_evidence(discovery, evidence),
        rules,
        panel_genes,
        scores,
        sample_score_tables,
        original,
        existing_posthoc,
        date.today().isoformat(),
    )
    print(f"Panel v2 genes: {', '.join(panel_genes)}")
    for accession, values in scores.items():
        print(
            f"{accession}: AUC={values['auc']:.4f}, "
            f"95% CI=({values['auc_ci95_lower']:.4f}, {values['auc_ci95_upper']:.4f}), "
            f"MW P={values['mann_whitney_p_two_sided']:.4g}, "
            f"matched-null P={values['matched_null_auc_empirical_p']}"
        )
    print(f"New audit and panel-v2 outputs saved under {OUTPUT}")


if __name__ == "__main__":
    main()
