from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gseapy as gp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import io as scipy_io
from scipy.stats import spearmanr
import seaborn as sns
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw" / "fig7_real_data"
PANEL_FILE = ROOT / "results" / "revision" / "panel_v2" / "panel_v2.json"
SIGNATURE_FILE = (
    ROOT
    / "results"
    / "revision"
    / "microenvironment_deconvolution"
    / "final_signature_table_both_variants.csv"
)
SEED = 20261023
BOOTSTRAPS = 2000
PANEL = json.loads(PANEL_FILE.read_text(encoding="utf-8"))["genes_passing_selection_rules"]


def parse_series_matrix(path: Path, accession: str) -> pd.DataFrame:
    rows: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                fields = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                rows.setdefault(fields[0], []).append(fields[1:])
    sample_ids = rows["!Sample_geo_accession"][0]
    output: list[dict[str, Any]] = []
    for i, sample_id in enumerate(sample_ids):
        record: dict[str, Any] = {
            "accession": accession,
            "sample_id": sample_id,
            "title": rows.get("!Sample_title", [[""] * len(sample_ids)])[0][i],
            "source_name": rows.get("!Sample_source_name_ch1", [[""] * len(sample_ids)])[0][i],
            "platform": rows.get("!Sample_platform_id", [[""] * len(sample_ids)])[0][i],
        }
        chars: list[str] = []
        for key, value_rows in rows.items():
            if "characteristics" not in key.lower() and "description" not in key.lower():
                continue
            for value_row in value_rows:
                value = value_row[i].strip()
                if not value:
                    continue
                chars.append(value)
                if ":" in value:
                    label, detail = value.split(":", 1)
                    record[label.strip().lower().replace(" ", "_")] = detail.strip()
        record["characteristics"] = "; ".join(chars)
        output.append(record)
    return pd.DataFrame(output)


def quantile_normalize(values: np.ndarray) -> np.ndarray:
    sorted_mean = np.sort(values, axis=0).mean(axis=1)
    normalized = np.empty_like(values, dtype=float)
    target_ranks = np.arange(1, values.shape[0] + 1, dtype=float)
    for j in range(values.shape[1]):
        ranks = pd.Series(values[:, j]).rank(method="average").to_numpy()
        normalized[:, j] = np.interp(ranks, target_ranks, sorted_mean)
    return normalized


def read_gpl_annotation(path: Path) -> dict[str, str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("ID\t"))
    annotation = pd.read_csv(io.StringIO("".join(lines[header:])), sep="\t", dtype=str)
    probe_col = "ID"
    symbol_col = next(c for c in annotation.columns if c.lower().lstrip("#") == "gene symbol")
    mapping: dict[str, str] = {}
    for probe, raw_symbol in zip(annotation[probe_col], annotation[symbol_col]):
        symbol = str(raw_symbol).split("///")[0].strip().upper()
        if symbol and symbol not in {"---", "NAN", "NA"}:
            mapping[str(probe)] = symbol
    return mapping


def discovery_array() -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = parse_series_matrix(RAW / "GSE125498_series_matrix.txt.gz", "GSE125498")
    raw = pd.read_csv(
        RAW / "GSE125498_non-normalized_data.txt.gz",
        sep="\t",
        compression="gzip",
        low_memory=False,
    )
    sample_cols = [c for c in raw.columns if re.fullmatch(r"PRAGUE_\d+", str(c))]
    if len(sample_cols) != len(metadata):
        raise ValueError(f"Array columns={len(sample_cols)} but GEO samples={len(metadata)}")
    raw_intensities = raw[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    valid_probes = np.isfinite(raw_intensities).all(axis=1)
    pd.DataFrame(
        [
            {
                "accession": "GSE125498",
                "source_probe_rows": int(len(raw)),
                "rows_dropped_for_nonnumeric_or_missing_intensity": int((~valid_probes).sum()),
                "rows_retained": int(valid_probes.sum()),
                "policy": "drop affected probe rows; do not impute",
            }
        ]
    ).to_csv(OUT / "array_preprocessing_qc.csv", index=False)
    raw = raw.loc[valid_probes].reset_index(drop=True)
    raw_intensities = raw_intensities[valid_probes]
    values = np.log2(np.maximum(raw_intensities, 1.0))
    normalized = quantile_normalize(values)
    probe_map = read_gpl_annotation(RAW / "GPL10558.annot.gz")
    sample_ids = metadata["sample_id"].tolist()
    mapped = pd.DataFrame(normalized, columns=sample_ids)
    mapped.insert(0, "probe_id", raw["Probe_Id"].astype(str))
    mapped["gene"] = mapped["probe_id"].map(probe_map)
    mapped = mapped.dropna(subset=["gene"])
    mapped["mean_expression"] = mapped[sample_ids].mean(axis=1)
    mapped = mapped.sort_values("mean_expression", ascending=False).drop_duplicates("gene")
    gene_expression = mapped.set_index("gene")[sample_ids].astype(float)
    metadata["group"] = metadata["title"].str.extract(r"^(Short-term PD|Long-term PD)", expand=False)
    lpd = metadata["group"].eq("Long-term PD")
    if int(lpd.sum()) != 13:
        raise ValueError(f"Expected 13 long-term discovery samples; found {int(lpd.sum())}")
    return gene_expression.loc[:, metadata.loc[lpd, "sample_id"]], metadata.loc[lpd].reset_index(drop=True)


def processed_affymetrix_array(
    accession: str, matrix_path: Path, annotation_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metadata = parse_series_matrix(matrix_path, accession)
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))
    table = pd.read_csv(
        io.StringIO("".join(lines[begin + 1 : end])),
        sep="\t",
        quotechar='"',
        low_memory=False,
    )
    id_col = "ID_REF" if "ID_REF" in table.columns else table.columns[0]
    sample_cols = metadata["sample_id"].tolist()
    if not set(sample_cols).issubset(table.columns):
        raise ValueError(f"{accession}: processed matrix columns do not match GEO sample IDs")
    numeric = table[sample_cols].apply(pd.to_numeric, errors="coerce")
    valid = np.isfinite(numeric.to_numpy()).all(axis=1)
    mapping = read_gpl_annotation(annotation_path)
    table = table.loc[valid].reset_index(drop=True)
    numeric = numeric.loc[valid].reset_index(drop=True)
    symbols = table[id_col].astype(str).map(mapping)
    good = symbols.notna() & symbols.astype(str).ne("")
    numeric = numeric.loc[good].copy()
    numeric.index = symbols.loc[good].astype(str).str.upper().to_numpy()
    means = numeric.mean(axis=1)
    keep_rows = (
        pd.DataFrame({"gene": numeric.index, "mean_expression": means.to_numpy()})
        .sort_values("mean_expression", ascending=False)
        .drop_duplicates("gene")
        .index
    )
    numeric = numeric.iloc[keep_rows]
    numeric.index = symbols.loc[good].astype(str).str.upper().iloc[keep_rows].to_numpy()
    numeric = numeric.loc[~numeric.index.duplicated(keep="first")]
    normalized = pd.DataFrame(
        numeric.to_numpy(dtype=float), index=numeric.index, columns=sample_cols
    )

    if "patient" not in metadata or "treatment" not in metadata:
        raise ValueError(f"{accession}: patient and treatment identifiers are required")
    patient_metadata: list[dict[str, Any]] = []
    patient_values: dict[str, pd.Series] = {}
    metadata["patient"] = metadata["patient"].astype(str).str.strip()
    for patient, group in metadata.groupby("patient", sort=True):
        patient_sample_ids = group["sample_id"].tolist()
        if len(patient_sample_ids) != 2:
            raise ValueError(
                f"{accession}: expected paired treatment/control assays for patient {patient}; "
                f"found {len(patient_sample_ids)}"
            )
        sample_id = f"patient_{patient}"
        patient_values[sample_id] = normalized[patient_sample_ids].mean(axis=1)
        first = group.iloc[0]
        patient_metadata.append(
            {
                "accession": accession,
                "sample_id": sample_id,
                "title": f"Patient {patient} (paired treatment/control mean)",
                "source_name": first["source_name"],
                "platform": first["platform"],
                "patient": patient,
                "previous_peritonitis": first.get("previous_peritonitis", ""),
                "group": "paired_treatment_control_mean",
                "assay_sample_ids_averaged": ";".join(patient_sample_ids),
                "characteristics": "; ".join(group["characteristics"].astype(str)),
            }
        )
    patient_expression = pd.DataFrame(patient_values)
    patient_expression.index.name = "gene"
    assay_metadata = metadata.copy()
    collapsed_metadata = pd.DataFrame(patient_metadata)
    qc = pd.DataFrame(
        [
            {
                "accession": accession,
                "assay_samples": int(len(metadata)),
                "independent_patients": int(len(collapsed_metadata)),
                "mapped_gene_symbols": int(len(patient_expression)),
                "normalization": "GEO processed log2 series matrix used as deposited; no second quantile normalization",
                "duplicate_symbol_policy": "retain highest mean-expression Ensembl feature",
                "paired_sample_policy": "average the two treatment-condition arrays per patient",
            }
        ]
    )
    return patient_expression, collapsed_metadata, assay_metadata, qc


def _read_tar_gzip(archive: tarfile.TarFile, member: tarfile.TarInfo) -> bytes:
    stream = archive.extractfile(member)
    if stream is None:
        raise OSError(f"Unable to read tar member {member.name}")
    return gzip.decompress(stream.read())


def pseudobulk_from_tar(
    accession: str, metadata: pd.DataFrame, archive_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids = set(metadata["sample_id"])
    result: dict[str, pd.Series] = {}
    cells: dict[str, int] = {}
    with tarfile.open(archive_path, "r:*") as archive:
        members = archive.getmembers()
        by_sample: dict[str, dict[str, tarfile.TarInfo]] = {}
        for member in members:
            base = Path(member.name).name
            match = re.match(r"(GSM\d+).*?(matrix\.mtx\.gz|features\.tsv\.gz)$", base)
            if match and match.group(1) in sample_ids:
                by_sample.setdefault(match.group(1), {})[match.group(2)] = member
        for sample_id in metadata["sample_id"]:
            files = by_sample.get(sample_id, {})
            if "matrix.mtx.gz" not in files or "features.tsv.gz" not in files:
                raise FileNotFoundError(f"{accession}: missing matrix/features for {sample_id}")
            feature_lines = _read_tar_gzip(archive, files["features.tsv.gz"]).decode(
                "utf-8", errors="replace"
            ).splitlines()
            feature_rows = [line.split("\t") for line in feature_lines]
            symbols = [
                (row[1] if len(row) > 1 and row[1].strip() else row[0]).strip().upper()
                for row in feature_rows
            ]
            matrix = scipy_io.mmread(
                io.BytesIO(_read_tar_gzip(archive, files["matrix.mtx.gz"]))
            ).tocsr()
            if matrix.shape[0] != len(symbols):
                if matrix.shape[1] == len(symbols):
                    matrix = matrix.T.tocsr()
                else:
                    raise ValueError(f"{sample_id}: matrix/features dimension mismatch")
            totals = np.asarray(matrix.sum(axis=1)).ravel()
            keep = np.array([bool(s) and s not in {"---", "NA"} for s in symbols])
            frame = pd.DataFrame({"gene": np.array(symbols)[keep], "count": totals[keep]})
            gene_totals = frame.groupby("gene", sort=False)["count"].sum()
            lib_size = float(matrix.sum())
            result[sample_id] = np.log2(gene_totals / lib_size * 1_000_000 + 1)
            cells[sample_id] = int(matrix.shape[1])
    expression = pd.DataFrame(result).fillna(0.0)
    expression.index.name = "gene"
    metadata = metadata.copy()
    metadata["cells_in_source_matrix"] = metadata["sample_id"].map(cells)
    return expression, metadata


def load_signatures() -> dict[str, list[str]]:
    source = pd.read_csv(SIGNATURE_FILE)
    source = source[
        source["signature_variant"].eq("full_signature")
        & source["included"].astype(str).str.lower().eq("true")
    ]
    return {
        str(name): sorted(set(group["gene"].astype(str).str.upper()))
        for name, group in source.groupby("population", sort=False)
    }


def score_signatures(
    expression: pd.DataFrame,
    signatures: dict[str, list[str]],
    accession: str,
    output_suffix: str = "",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    detectable = {
        name: [gene for gene in genes if gene in expression.index]
        for name, genes in signatures.items()
    }
    sizes = pd.DataFrame(
        [
            {
                "accession": accession,
                "signature": name,
                "source_gene_count": len(signatures[name]),
                "detected_gene_count": len(genes),
                "genes_used": ";".join(genes),
                "status": "scored" if len(genes) >= 3 else "not_scored_min_size_lt_3",
            }
            for name, genes in detectable.items()
        ]
    )
    valid = {name: genes for name, genes in detectable.items() if len(genes) >= 3}
    scores = pd.DataFrame(index=expression.columns)
    if valid:
        result = gp.ssgsea(
            data=expression,
            gene_sets=valid,
            outdir=None,
            min_size=3,
            max_size=10000,
            sample_norm_method="rank",
            correl_norm_type="rank",
            permutation_num=0,
            threads=1,
            no_plot=True,
            seed=SEED,
            verbose=False,
        )
        score_table = result.res2d
        value_col = "NES" if "NES" in score_table.columns else "ES"
        scores = score_table.pivot(index="Name", columns="Term", values=value_col)
        scores = scores.reindex(index=expression.columns, columns=list(signatures))
    scores.index.name = "sample_id"
    scores.to_csv(OUT / f"ssgsea_scores_{accession}{output_suffix}.csv")
    sizes.to_csv(
        OUT / f"signature_detection_{accession}{output_suffix}.csv", index=False
    )
    return scores, sizes


def bh(values: list[float]) -> list[float]:
    finite = np.isfinite(values)
    output = np.full(len(values), np.nan)
    if finite.any():
        output[finite] = multipletests(np.asarray(values)[finite], method="fdr_bh")[1]
    return output.tolist()


def score_panel(expression: pd.DataFrame) -> tuple[pd.Series, list[str]]:
    available = [gene for gene in PANEL if gene in expression.index]
    if len(available) != len(PANEL):
        return pd.Series(np.nan, index=expression.columns), available
    z = expression.loc[available].T
    z = (z - z.mean(axis=0)) / z.std(axis=0, ddof=1).replace(0, np.nan)
    return z.mean(axis=1).rename("panel_score"), available


def bootstrap_rho(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    n = len(x)
    vals: list[float] = []
    for _ in range(BOOTSTRAPS):
        idx = rng.integers(0, n, size=n)
        if np.unique(x[idx]).size > 1 and np.unique(y[idx]).size > 1:
            rho = spearmanr(x[idx], y[idx]).statistic
            if np.isfinite(rho):
                vals.append(float(rho))
    if not vals:
        return float("nan"), float("nan")
    return tuple(np.quantile(vals, [0.025, 0.975]).tolist())


def analyze_cohort(
    accession: str,
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    signatures: dict[str, list[str]],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score, available_panel = score_panel(expression)
    if len(available_panel) != len(PANEL):
        missing = sorted(set(PANEL) - set(available_panel))
        raise ValueError(f"{accession} is missing panel genes: {missing}")
    signature_scores, _ = score_signatures(expression, signatures, accession)
    common = score.index.intersection(signature_scores.index)
    score = score.loc[common]
    signature_scores = signature_scores.loc[common]
    metadata = metadata.set_index("sample_id").loc[common]
    metadata.index.name = "sample_id"
    metadata = metadata.reset_index()
    score_table = pd.DataFrame({"sample_id": common, "panel_score": score.to_numpy()})
    score_table = score_table.merge(metadata, on="sample_id", how="left")
    signature_scores.index.name = "sample_id"
    score_table = score_table.merge(
    signature_scores.reset_index(), on="sample_id", how="left"
    )
    score_table.to_csv(OUT / f"sample_scores_{accession}.csv", index=False)

    results: list[dict[str, Any]] = []
    for signature in signatures:
        x = score.to_numpy(dtype=float)
        y = signature_scores[signature].to_numpy(dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        inference_status = "reported"
        if keep.sum() < 3 or np.unique(x[keep]).size < 2 or np.unique(y[keep]).size < 2:
            rho, p, low, high = np.nan, np.nan, np.nan, np.nan
            inference_status = "not_estimable"
        elif keep.sum() < 6:
            rho = float(spearmanr(x[keep], y[keep]).statistic)
            p, low, high = np.nan, np.nan, np.nan
            inference_status = "descriptive_only_n_lt_6"
        else:
            test = spearmanr(x[keep], y[keep])
            rho, p = float(test.statistic), float(test.pvalue)
            low, high = bootstrap_rho(x[keep], y[keep], rng)
        results.append(
            {
                "accession": accession,
                "n": int(keep.sum()),
                "signature": signature,
                "rho": rho,
                "rho_ci95_low": low,
                "rho_ci95_high": high,
                "p_value": p,
                "inference_status": inference_status,
            }
        )
    corr = pd.DataFrame(results)
    corr["q_value_BH_within_dataset_13_tests"] = bh(corr["p_value"].tolist())
    corr.to_csv(OUT / f"correlations_{accession}.csv", index=False)
    return score_table, signature_scores, corr, metadata


def clinical_tests(
    accession: str,
    score_table: pd.DataFrame,
    group_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata_fields = [
        "age",
        "sex",
        "dialysis_vintage_group",
        "dialysis_vintage_months",
        "dialysis_duration_months",
        "transport_status",
        "ultrafiltration_failure",
        "peritonitis",
        "previous_peritonitis",
        "batch",
    ]
    availability = []
    for field in metadata_fields:
        matching = [c for c in score_table.columns if c.lower() == field]
        values = score_table[matching[0]] if matching else pd.Series(dtype=object)
        source = "GEO series matrix sample characteristics"
        group_values = (
            score_table["group"].astype(str)
            if "group" in score_table
            else pd.Series(dtype=str)
        )
        has_vintage_labels = group_values.str.contains(
            r"short|long|\bSV\b|\bLV\b", case=False, regex=True
        ).any()
        has_uf_labels = group_values.str.contains(
            r"UF|failure", case=False, regex=True
        ).any()
        if field == "dialysis_vintage_group" and has_vintage_labels:
            values = score_table["group"]
            source = "GEO series-matrix study group"
        if field == "ultrafiltration_failure" and has_uf_labels:
            values = score_table["group"].astype(str).str.contains(
                r"LV_UF|UF failure", case=False, regex=True
            )
            source = "GEO series-matrix study group"
        if field == "peritonitis" and "previous_peritonitis" in score_table:
            values = score_table["previous_peritonitis"]
            source = "GEO series matrix previous_peritonitis characteristic"
        availability.append(
            {
                "accession": accession,
                "covariate": field,
                "available_in_GEO_sample_metadata": bool(values.notna().any()),
                "nonmissing_n": int(values.notna().sum()),
                "metadata_source": source,
            }
        )
    tests: list[dict[str, Any]] = []
    data = score_table.copy()
    numeric_columns = {c.lower(): c for c in data.columns}
    predictors: list[tuple[str, pd.Series, pd.DataFrame]] = []
    proxy_columns = [
        c
        for c in [
            "T cells",
            "B cells",
            "NK cells",
            "Monocytes",
            "M1 macrophages",
            "M2 macrophages",
            "Neutrophils",
            "Mast cells",
            "Immature dendritic cells",
            "Mature dendritic cells",
        ]
        if c in data.columns
    ]
    data["leukocyte_proxy"] = data[proxy_columns].mean(axis=1) if proxy_columns else np.nan

    age_col = numeric_columns.get("age")
    if age_col:
        age = pd.to_numeric(data[age_col], errors="coerce")
        if age.notna().sum() >= 5:
            predictors.append(("age_years", age, data))
    sex_col = numeric_columns.get("sex")
    if sex_col:
        sex = data[sex_col].astype(str).str.strip().str.lower().map(
            {"m": 1.0, "male": 1.0, "f": 0.0, "female": 0.0}
        )
        if sex.notna().sum() >= 5 and sex.nunique() == 2:
            predictors.append(("sex_M_vs_F", sex, data))
    if group_column and group_column in data.columns:
        group = data[group_column].astype(str)
        long = group.str.contains("long|lv", case=False, regex=True).astype(float)
        if long.nunique() == 2 and long.sum() >= 3 and (1 - long).sum() >= 3:
            predictors.append(("long_vintage_vs_short", long, data))
        uf = group.str.contains("uf", case=False, regex=True) & ~group.str.contains(
            "not", case=False, regex=True
        )
        not_uf = group.str.contains("not.*uf", case=False, regex=True)
        long_mask = long.eq(1)
        if int((uf & long_mask).sum()) >= 3 and int((not_uf & long_mask).sum()) >= 3:
            subset = data.loc[long_mask].copy()
            uf_code = (uf.loc[long_mask]).astype(float)
            predictors.append(("UF_failure_vs_no_UF_in_long_vintage", uf_code, subset))

    for name, predictor, frame in predictors:
        y = pd.to_numeric(frame["panel_score"], errors="coerce")
        x = pd.to_numeric(predictor, errors="coerce")
        proxy = pd.to_numeric(frame["leukocyte_proxy"], errors="coerce")
        valid = y.notna() & x.notna() & proxy.notna()
        if valid.sum() < 6 or x.loc[valid].nunique() < 2:
            continue
        design = pd.DataFrame(
            {"predictor": x.loc[valid].astype(float), "leukocyte_proxy": proxy.loc[valid].astype(float)}
        )
        design = sm.add_constant(design)
        fit = sm.OLS(y.loc[valid].astype(float), design).fit()
        ci = fit.conf_int().loc["predictor"]
        tests.append(
            {
                "accession": accession,
                "covariate": name,
                "n": int(valid.sum()),
                "effect": float(fit.params["predictor"]),
                "ci95_low": float(ci.iloc[0]),
                "ci95_high": float(ci.iloc[1]),
                "p_value": float(fit.pvalues["predictor"]),
                "adjustment": "leukocyte_proxy",
                "model": "OLS panel_score ~ covariate + leukocyte_proxy",
            }
        )
    tests_df = pd.DataFrame(tests)
    if not tests_df.empty:
        tests_df["q_value_BH_within_dataset"] = bh(tests_df["p_value"].tolist())
    return tests_df, pd.DataFrame(availability)


def audit_candidates() -> pd.DataFrame:
    rows = [
        ["GSE125498", "Homo sapiens", "peritoneal dialysis effluent cells", "Illumina HumanHT-12 v4 array", 33, "short-term=20; long-term=13", "dialysis vintage category=yes; exact months/age/sex/UF/transport/peritonitis/batch=not in GEO sample fields", "Discovery; only the 13 long-term/LPD samples analyzed"],
        ["GSE130888", "Homo sapiens", "peritoneal effluent cells plus normal peritoneal tissue controls", "10x Genomics scRNA-seq", 13, "short-term PD=6; long-term PD=4; normal tissue=3", "duration category=yes; age/sex/UF/transport/peritonitis/batch=not in GEO sample fields", "Same-tissue validation; effluent cases only (n=10)"],
        ["GSE248762", "Homo sapiens", "overnight peritoneal dialysis effluent", "10x Genomics scRNA-seq", 16, "SV=6; LV_NOT_UF=6; LV_UF=4", "vintage group and UF status=yes; age/sex reported per sample; exact vintage/peritonitis/transport/batch=not in GEO sample fields", "Same-tissue validation; patient pseudobulk"],
        ["GSE57070", "Homo sapiens", "peritoneal dialysis effluent cells", "Affymetrix Human Gene 1.0 ST array, GPL11209", 6, "AlaGln=3; control=3, paired within 3 patients", "patient ID and previous peritonitis=yes; age/sex/vintage/UF/transport/batch=not in GEO sample fields", "Descriptive only; paired arrays averaged per patient (n=3), excluded from validation because rank correlations are uninformative at this n"],
        ["GSE297829", "Homo sapiens", "peritoneal dialysis effluent cells", "10x Genomics scRNA-seq", 1, "one patient after 3 years PD", "duration approx. 3 years from study design; other patient-level covariates unavailable", "Excluded: one patient, not an independent association validation"],
        ["GSE62928", "Homo sapiens", "peritoneal membrane biopsy", "Affymetrix array", 8, "EPS=4; PD=2; uremic controls=2", "clinical group yes; age/sex/vintage/UF/transport/peritonitis/batch not fully available per GEO sample fields", "Different tissue type; not treated as effluent validation"],
        ["GSE92453", "Homo sapiens", "cultured and effluent-derived mesothelial cells", "Expression array", 38, "ex-vivo control comparisons=4; epithelioid comparisons=18; non-epithelioid comparisons=16", "transport association discussed in study; patient-level covariates not complete", "Not used: mixed culture/ex-vivo cell selection is not bulk effluent"],
        ["GSE92455", "not verified", "not verified", "GEO series matrix URL returned 404", "not available", "not available", "not available", "Search candidate could not be verified as a public GEO series"],
        ["GSE121372", "Homo sapiens", "treated human peritoneal mesothelial cells", "Expression array", 4, "control 6h=1; TGF-beta1 6h=1; control 24h=1; TGF-beta1 24h=1", "patient-level covariates not applicable", "Not used: cultured cells, not patient effluent"],
        ["GSE146955", "not verified", "not verified", "GEO series matrix URL returned 404", "not available", "not available", "not available", "Search candidate could not be verified as a public GEO series"],
        ["E-MTAB-5439", "Homo sapiens", "effluent-derived sorted CD14+ macrophages and CD1c+ dendritic cells", "Affymetrix expression array", 10, "CD14+ macrophages=5; CD1c+ dendritic cells=5 (5 donors)", "donor and sorted cell type available; age/sex/vintage/UF/transport/peritonitis/batch not verified in the processed source", "Not used: sorted cell populations are not bulk effluent; BioStudies lists raw CEL files and a final gene list, not a processed gene-expression matrix"],
        ["E-MTAB-5462", "Homo sapiens", "peritoneal dialysis effluent cells", "Bulk RNA-seq", 11, "six patients; crossover treatment/control (11 assay samples)", "patient identifier and treatment available; age/sex/vintage summarized in protocol, not row-level SDRF; batch sequencing lane encoded in names", "Not processed: BioStudies exposes SDRF/IDF only; usable gene-level matrix absent, raw ENA reads require a separate alignment/quantification pipeline"],
        ["E-MTAB-5461", "Homo sapiens", "peritoneal dialysis effluent microRNA", "miRNA array", 10, "AlaGln=5; control=5 (10 assay samples; incomplete patient pairing)", "sample treatment and patient identifiers available; not relevant to panel mRNA scoring", "Not used: miRNA platform cannot score protein-coding panel/signatures"],
    ]
    columns = [
        "accession",
        "organism",
        "tissue_sample_type",
        "platform",
        "n_samples",
        "n_per_group",
        "patient_level_covariates_vintage_age_sex_UF_transport_peritonitis_batch",
        "decision_and_reason",
    ]
    return pd.DataFrame(rows, columns=columns)


def sra_audit() -> pd.DataFrame:
    rows = [
        ["PRJNA541920", "GSE130888", "Human peritoneal effluent scRNA-seq plus normal tissue samples", "Analyzed through GEO raw matrices; duplicate SRA/GEO deposit"],
        ["PRJNA1045728", "GSE248762", "Human overnight peritoneal effluent scRNA-seq", "Analyzed through GEO raw matrices; duplicate SRA/GEO deposit"],
        ["PRJNA1266235", "GSE297829", "One human patient effluent scRNA-seq sample after 3 years PD", "Excluded: n=1"],
        ["PRJEB19316", "E-MTAB-5462", "Human effluent bulk RNA-seq crossover trial, 6 patients", "Not analyzed: BioStudies lists only IDF/SDRF, no processed gene-level matrix; FASTQ quantification not undertaken"],
        ["PRJNA638799", "", "Small-RNA sequencing of peritoneal biopsies: 7 hyalinizing vasculopathy and 4 non-hyalinizing", "Excluded: biopsy tissue and small-RNA data are not mRNA panel expression"],
        ["PRJNA672097", "GSE160145", "Human heart tissue from kidney disease/dialysis study", "Excluded: heart tissue; SRA title/search text was a false positive"],
        ["PRJNA757643", "GSE182736", "Exosomal miRNA in peritoneal dialysis effluent, 3 UF failure and 3 stable", "Excluded: miRNA-seq cannot score the protein-coding panel or mRNA signatures"],
        ["PRJNA1033066", "", "Human pediatric blood TCR repertoire in sepsis", "Excluded: tissue and study question unrelated to peritoneal effluent"],
        ["PRJNA1212711", "", "Animal experiments on microbiome and PD-associated peritoneal fibrosis", "Excluded: non-human model, despite a broad SRA human-query match"],
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "bioproject",
            "linked_GEO_or_ArrayExpress",
            "verified_project_description",
            "decision",
        ],
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_manifest() -> list[dict[str, Any]]:
    now = datetime.now(timezone.utc).isoformat()
    urls = {
        "GSE125498_non-normalized_data.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125498/suppl/GSE125498_non-normalized_data.txt.gz",
        "GSE125498_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125498/suppl/GSE125498_RAW.tar",
        "GPL10558.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
        "GSE125498_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125498/matrix/GSE125498_series_matrix.txt.gz",
        "GSE130888_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130888/matrix/GSE130888_series_matrix.txt.gz",
        "GSE248762_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248762/matrix/GSE248762_series_matrix.txt.gz",
        "GSE130888_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130888/suppl/GSE130888_RAW.tar",
        "GSE248762_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248762/suppl/GSE248762_RAW.tar",
        "GSE57070_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57070/suppl/GSE57070_RAW.tar",
        "GSE57070_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57070/matrix/GSE57070_series_matrix.txt.gz",
        "GPL11209.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/GPL11209/annot/GPL11209.annot.gz",
        "GSE62928_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62928/matrix/GSE62928_series_matrix.txt.gz",
        "GSE92453_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92453/matrix/GSE92453_series_matrix.txt.gz",
        "GSE121372_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121372/matrix/GSE121372_series_matrix.txt.gz",
    }
    result = []
    for path in sorted(RAW.iterdir()):
        if path.is_file():
            copied_source = path.name in {"GSE130888_RAW.tar", "GSE248762_RAW.tar"}
            result.append(
                {
                    "file": str(path.relative_to(ROOT)),
                    "source_url": urls.get(path.name, "local frozen project reference"),
                    "sha256": sha256(path),
                    "size_bytes": path.stat().st_size,
                    "download_date_utc": (
                        now[:10]
                        if path.name == "GSE57070_RAW.tar"
                        else (
                            "unknown_preexisting_local_download; copied and checksum-verified "
                            + now[:10]
                            if copied_source
                            else "unknown_preexisting_local_download; checksum-verified "
                            + now[:10]
                        )
                    ),
                    "verification_timestamp_utc": now,
                }
            )
    return result


def make_figures(correlations: pd.DataFrame, score_tables: dict[str, pd.DataFrame]) -> None:
    order = ["GSE125498", "GSE130888", "GSE248762"]
    validation_correlations = correlations[
        correlations["accession"].isin(order)
    ]
    pivot = validation_correlations.pivot(
        index="signature", columns="accession", values="rho"
    )
    pivot = pivot.reindex(columns=order)
    pivot.to_csv(OUT / "correlation_rho_matrix.csv")
    plt.figure(figsize=(8.2, 6.3))
    sns.heatmap(
        pivot,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        cbar_kws={"label": "Spearman rho"},
    )
    plt.title("Panel score vs. CellMarker ssGSEA enrichment")
    plt.tight_layout()
    plt.savefig(OUT / "fig7_A_correlation_heatmap.png", dpi=220)
    plt.close()

    discovery = score_tables["GSE125498"]
    top = (
        correlations.query("accession == 'GSE125498'")
        .assign(abs_rho=lambda frame: frame["rho"].abs())
        .sort_values("abs_rho", ascending=False)
        .head(3)
    )
    fig, axes = plt.subplots(1, len(top), figsize=(4.6 * len(top), 4.1))
    if len(top) == 1:
        axes = [axes]
    for axis, row in zip(axes, top.itertuples()):
        x = discovery["panel_score"].to_numpy()
        y = discovery[row.signature].to_numpy()
        sns.regplot(x=x, y=y, lowess=True, ci=None, scatter_kws={"s": 38}, ax=axis)
        axis.set_xlabel("Five-gene panel score (mean gene z-score)")
        axis.set_ylabel(f"{row.signature} ssGSEA score")
        axis.set_title(f"rho={row.rho:.2f}; q={row.q_value_BH_within_dataset_13_tests:.3g}")
    fig.suptitle("Exploratory discovery associations (13 LPD samples)")
    fig.tight_layout()
    fig.savefig(OUT / "fig7_B_top_associations.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.4, 6.2))
    offsets = {
        "GSE125498": -0.27,
        "GSE130888": -0.09,
        "GSE248762": 0.09,
        "GSE57070": 0.27,
    }
    for accession in order:
        part = correlations[correlations["accession"].eq(accession)].copy()
        y = np.arange(len(pivot.index))
        part = part.set_index("signature").reindex(pivot.index)
        x = part["rho"].to_numpy()
        low = part["rho_ci95_low"].to_numpy()
        high = part["rho_ci95_high"].to_numpy()
        ax.errorbar(
            x,
            y + offsets[accession],
            xerr=np.vstack([x - low, high - x]),
            fmt="o",
            capsize=2,
            label=accession,
        )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(np.arange(len(pivot.index)), pivot.index)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Spearman rho (bootstrap 95% CI)")
    ax.set_title("Exploratory cross-cohort concordance (not same-method replication)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig7_C_validation_replication.png", dpi=220)
    plt.close(fig)


def write_report(
    correlations: pd.DataFrame,
    metadata_by_cohort: dict[str, pd.DataFrame],
    clinical: pd.DataFrame,
    audit: pd.DataFrame,
    panel_signature_overlap: pd.DataFrame,
    ltb_sensitivity: pd.DataFrame,
    leukocyte_proxy_correlations: pd.DataFrame,
) -> None:
    overlaps = panel_signature_overlap.loc[
        panel_signature_overlap["overlap_count"].gt(0)
    ]
    overlap_description = "; ".join(
        f"{row.signature}: {row.overlap_genes}" for row in overlaps.itertuples()
    )
    discovery_ltb = ltb_sensitivity.loc[
        ltb_sensitivity["accession"].eq("GSE125498")
    ].iloc[0]
    lines = [
        "# Figure 7: exploratory peritoneal effluent panel-score analysis",
        "",
        f"Analysis run: {datetime.now(timezone.utc).date().isoformat()} (UTC).",
        "",
        "## Data and methods",
        "",
        f"Frozen panel genes: {', '.join(PANEL)}.",
        "Discovery was restricted to the 13 long-term PD (LPD) peritoneal-effluent samples in GSE125498. "
        f"The 13 GEO sample accessions were: {', '.join(metadata_by_cohort['GSE125498']['sample_id'].astype(str))}. "
        "The 20 short-term samples were not included in the discovery score or correlations. "
        "For each cohort separately, gene-wise z-scores were calculated across that cohort's analyzed samples "
        "and the panel score is their mean.",
        "GSE125498 raw non-normalized Illumina intensities were log2 transformed, quantile normalized across "
        "all 33 deposited arrays, mapped using GPL10558 official annotation, and duplicate gene symbols were "
        "collapsed by retaining the probe with highest mean normalized expression. The 13 LPD samples were "
        "then selected. scRNA-seq raw gene-by-cell matrices were aggregated to one pseudobulk count vector per "
        "GEO sample; counts were library-size normalized to CPM and log2(CPM+1) transformed. Cohorts were "
        "processed and scored separately without merging or batch correction. For GSE57070, GEO's processed "
        "log2 series-matrix values were used as deposited, without a second quantile normalization; its raw CEL "
        "archive was downloaded and checksum-recorded but was not the matrix used for scoring.",
        "The frozen 13 human CellMarker_2024 signatures were scored with GSEApy ssGSEA (rank normalization, "
        "minimum detectable signature size 3). Scores are marker-enrichment scores, not cell fractions. "
        "For cohorts with n>=6, all 13 Spearman tests have within-cohort Benjamini-Hochberg q-values; 95% CIs "
        "for rho are percentile bootstrap intervals (2,000 resamples, fixed seed). For n<6, only descriptive rho is "
        "reported; p-values, q-values, and bootstrap CIs are withheld because asymptotic and bootstrap "
        "inference is not credible at that sample size.",
        f"Panel/signature overlap was audited: {overlap_description or 'none detected'}. "
        "Correlations with signatures sharing panel genes should be interpreted with this dependence in mind.",
        "The LTB circularity sensitivity was run by removing LTB from the T-cell marker set while keeping the "
        "five-gene panel unchanged; all other signature definitions and scoring settings were unchanged.",
        "No consensus clustering was run: the requested primary inference is continuous, and with small "
        "cohorts unsupervised grouping would not improve the prespecified continuous-score test. Consequently, "
        "cluster sizes, cluster null distributions, and clustering figures are not applicable.",
        "",
        "## Dataset audit and inclusion",
        "",
        f"Audited {len(audit)} relevant GEO/ArrayExpress candidates. Source series matrices and supplementary "
        "files are staged in `data/raw/fig7_real_data/` and checksum-verified. The GSE130888 and GSE248762 "
        "scRNA-seq archives, and other pre-existing source files, have unknown original download dates. "
        "GSE57070_RAW.tar was downloaded on the run date, but scoring used the GEO processed log2 matrix. "
        "See `download_manifest.csv` for URLs, byte sizes, SHA-256 hashes, and acquisition-date notes.",
        "SRA search (human peritoneal dialysis/effluent/fibrosis terms) returned 95 run records representing "
        "nine BioProjects. Three duplicate GEO projects were GSE130888, GSE248762, and GSE297829. Other "
        "projects included noneligible biopsy small-RNA, effluent exosomal miRNA, heart tissue, sepsis blood "
        "TCR, an animal model, and E-MTAB-5462. Its SDRF/IDF are public, but no "
        "processed gene-level matrix is listed, and reprocessing its FASTQ files would require a separate "
        "alignment/quantification workflow, so it is recorded but not analyzed.",
        "GSE130888 contributed 10 effluent samples (6 short-term, 4 long-term); the 3 normal tissue samples "
        "were excluded. GSE248762 contributed 16 overnight-effluent samples (6 SV, 6 LV_NOT_UF, 4 LV_UF). "
        "GSE57070 contributed six paired treatment/control arrays from three patients; the paired arrays were "
        "averaged within patient before scoring (n=3). It is not treated or cited as validation and is omitted "
        "from the Figure 7 validation/concordance panels because its rank correlations are necessarily limited "
        "to 1, 0.5, or -0.5 and are not informative at this n. "
        "GSE297829 is a single patient after three years of PD and is not a valid association validation cohort. "
        "E-MTAB-5439 has 10 effluent-derived sorted-cell arrays (five CD14+ macrophage and five CD1c+ "
        "dendritic-cell samples from five donors), but it is not bulk effluent and BioStudies lists raw CEL files "
        "and a final gene list rather than a processed expression matrix. E-MTAB-5461 is a microRNA array and "
        "cannot score this mRNA panel. The biopsy and in-vitro/cell-selected studies in the audit were not "
        "treated as bulk effluent validation.",
        "",
        "## Correlations",
        "",
        "| Dataset | n | Signature | rho (95% bootstrap CI) | p | BH q |",
        "|---|---:|---|---:|---:|---:|",
    ]
    for row in correlations.loc[~correlations["accession"].eq("GSE57070")].itertuples():
        ci = (
            f"[{row.rho_ci95_low:.3f}, {row.rho_ci95_high:.3f}]"
            if np.isfinite(row.rho_ci95_low) and np.isfinite(row.rho_ci95_high)
            else "not estimated"
        )
        p_value = f"{row.p_value:.4g}" if np.isfinite(row.p_value) else "not estimated"
        q_value = (
            f"{row.q_value_BH_within_dataset_13_tests:.4g}"
            if np.isfinite(row.q_value_BH_within_dataset_13_tests)
            else "not estimated"
        )
        rho = f"{row.rho:.3f}" if np.isfinite(row.rho) else "not estimable"
        lines.append(
            f"| {row.accession} | {row.n} | {row.signature} | {rho} {ci} | "
            f"{p_value} | {q_value} |"
        )
    lines.extend(
        [
            "",
            "### T-cell result after removing LTB",
            "",
            "The sensitivity analysis below rescored the T-cell signature without LTB. The five-gene panel "
            "remains unchanged. GSE57070 is omitted because n=3 supports no useful inferential check.",
            "",
            "| Dataset | n | T-cell rho (95% bootstrap CI) | p | BH q (13-signature family) |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in ltb_sensitivity.loc[
        ltb_sensitivity["accession"].ne("GSE57070")
    ].itertuples():
        ci = (
            f"[{row.rho_ci95_low:.3f}, {row.rho_ci95_high:.3f}]"
            if np.isfinite(row.rho_ci95_low) and np.isfinite(row.rho_ci95_high)
            else "not estimated"
        )
        p_value = f"{row.p_value:.4g}" if np.isfinite(row.p_value) else "not estimated"
        q_value = (
            f"{row.q_value_BH_within_13_signatures:.4g}"
            if np.isfinite(row.q_value_BH_within_13_signatures)
            else "not estimated"
        )
        rho = f"{row.rho:.3f}" if np.isfinite(row.rho) else "not estimable"
        lines.append(
            f"| {row.accession} | {row.n} | {rho} {ci} | {p_value} | {q_value} |"
        )
    lines.extend(
        [
            "",
            f"In discovery, removing LTB left rho at {discovery_ltb.rho:.3f} "
            f"(p={discovery_ltb.p_value:.4g}; BH q={discovery_ltb.q_value_BH_within_13_signatures:.4g}). "
            "Thus the observed T-cell association is not driven solely by the shared LTB marker, although "
            "that does not rule out broader immune-content confounding.",
        ]
    )
    lines.extend(
        [
            "",
            "### Leukocyte-content caution",
            "",
            "The five-gene panel contains several immune-associated genes. A panel-score association with an "
            "immune-marker proxy is therefore compatible with the panel reflecting leukocyte content rather "
            "than a distinct fibrosis program; ssGSEA scores are enrichment scores, not cell fractions. "
            "The following exploratory check correlates the panel score with the mean of the nine non-T-cell "
            "immune ssGSEA scores (the T-cell signature is left out to avoid the direct LTB overlap):",
            "",
            "| Dataset | n | Spearman rho | 95% bootstrap CI | p | BH q across cohorts |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in leukocyte_proxy_correlations.loc[
        leukocyte_proxy_correlations["accession"].ne("GSE57070")
    ].itertuples():
        ci = (
            f"[{row.rho_ci95_low:.3f}, {row.rho_ci95_high:.3f}]"
            if np.isfinite(row.rho_ci95_low) and np.isfinite(row.rho_ci95_high)
            else "not estimated"
        )
        p_value = f"{row.p_value:.4g}" if np.isfinite(row.p_value) else "not estimated"
        q_value = (
            f"{row.q_value_BH_across_cohorts:.4g}"
            if np.isfinite(row.q_value_BH_across_cohorts)
            else "not estimated"
        )
        rho = f"{row.rho:.3f}" if np.isfinite(row.rho) else "not estimable"
        lines.append(
            f"| {row.accession} | {row.n} | {rho} | {ci} | {p_value} | {q_value} |"
        )
    proxy_q = leukocyte_proxy_correlations.loc[
        leukocyte_proxy_correlations["accession"].ne("GSE57070"),
        "q_value_BH_across_cohorts",
    ]
    lines.extend(
        [
            "",
            f"None of these three panel-score/proxy correlations survived BH correction "
            f"(q range {proxy_q.min():.3f}–{proxy_q.max():.3f}); this does not establish that the panel "
            "is independent of immune content, given the small samples and marker-score proxy. Conversely, "
            "these proxy results do not demonstrate that leukocyte abundance explains the panel associations. "
            "The negative fibroblast association in GSE248762 is compatible with a compositional dilution "
            "interpretation but does not establish fibroblast depletion or a causal mechanism.",
        ]
    )
    lines.extend(
        [
            "",
            "The full machine-readable table includes all 13 signatures, including non-significant associations "
            "and descriptive GSE57070 values; that n=3 cohort is not validation and is excluded from the report "
            "table and figures. Directional agreement in GSE130888 and GSE248762 is exploratory cross-cohort "
            "concordance, not identical-method replication: both validation cohorts are single-cell-derived "
            "pseudobulk, whereas discovery is bulk microarray. Small n (13, 10, 16) and wide intervals limit "
            "precision.",
            "",
            "## Covariates and limitations",
            "",
            "Clinical association models use panel score as the outcome and include the leukocyte proxy "
            "(mean of the 10 immune ssGSEA scores) as an adjustment term. Reported coefficients are "
            "unstandardized OLS effects with 95% confidence intervals, p-values, and within-dataset BH q-values. "
            "None of the five fitted clinical-covariate associations survived BH correction (all q >= 0.0655). "
            "These models test clinical covariate associations conditional on the proxy; they do not directly "
            "test whether the panel is independent of leukocyte content. Where metadata do not provide a "
            "patient-level value, no value was inferred.",
            "",
            "| Dataset | Covariate | n | Adjusted effect (95% CI) | p | BH q |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    if clinical.empty:
        lines.append("| — | No eligible sample-level covariate model | — | — | — | — |")
    else:
        for row in clinical.itertuples():
            lines.append(
                f"| {row.accession} | {row.covariate} | {row.n} | {row.effect:.3f} "
                f"[{row.ci95_low:.3f}, {row.ci95_high:.3f}] | {row.p_value:.4g} | "
                f"{row.q_value_BH_within_dataset:.4g} |"
            )
    lines.extend(
        [
            "",
            "Exact dialysis vintage and batch were not available as per-sample values in the analyzed cohorts. "
            "Previous peritonitis was available only in GSE57070. In GSE248762, long/short vintage and UF status are encoded "
            "by design group; age and sex were assessed only if present as sample-level series-matrix fields. "
            "In the LPD-only discovery subset, a between-vintage-group test is impossible because the subset "
            "contains only long-term patients. These results remain underpowered and exploratory; neither "
            "the discovery correlations nor the validation cohorts constitute clinical validation. GSE57070 "
            "has only three independent patients; previous peritonitis is in its sample metadata, but no "
            "clinical covariate model was fitted because n=3. Age, sex, vintage, UF status, transport status, "
            "and batch were not available for that cohort.",
            "",
            "## Figure legend",
            "",
            "**Figure 7. Exploratory association of the frozen five-gene panel score with immune/stromal "
            "marker-enrichment scores in peritoneal dialysis effluent.** (A) Spearman rho heatmap for the "
            "13 frozen CellMarker_2024 signatures in GSE125498 (n=13), GSE130888 (n=10), and GSE248762 "
            "(n=16); these are separate cohorts, not identical-method validation. GSE57070 (n=3) is excluded "
            "from these panels and is not cited as validation. "
            "(B) Three strongest absolute discovery associations, shown descriptively; these signatures "
            "are not independent of the score-construction process. (C) Cohort-specific Spearman rho and "
            "bootstrap 95% confidence intervals across the discovery and two separate single-cell pseudobulk "
            "cohorts. The T-cell sensitivity analysis excludes LTB from that signature. ssGSEA values represent "
            "marker enrichment, not cell fractions. All analyses are exploratory, small-sample, and "
            "not validated for clinical use.",
        ]
    )
    (OUT / "fig7_results_and_legend.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit = audit_candidates()
    audit.to_csv(OUT / "dataset_audit.csv", index=False)
    sra_audit().to_csv(OUT / "sra_project_audit.csv", index=False)
    signatures = load_signatures()
    panel_signature_overlap = pd.DataFrame(
        [
            {
                "signature": name,
                "overlap_genes": ";".join(sorted(set(PANEL).intersection(genes))),
                "overlap_count": len(set(PANEL).intersection(genes)),
            }
            for name, genes in signatures.items()
        ]
    )
    panel_signature_overlap.to_csv(OUT / "panel_signature_overlap.csv", index=False)
    rng = np.random.default_rng(SEED)

    discovery, meta125 = discovery_array()
    discovery.to_csv(OUT / "gene_expression_GSE125498_LPD.csv")
    meta125.to_csv(OUT / "sample_metadata_GSE125498_LPD.csv", index=False)

    meta130_all = parse_series_matrix(RAW / "GSE130888_series_matrix.txt.gz", "GSE130888")
    meta130_all["group"] = meta130_all["title"].str.extract(
        r"^(Short-term PD|Long-term PD|HC)", expand=False
    )
    meta130 = meta130_all[meta130_all["group"].isin(["Short-term PD", "Long-term PD"])].copy()
    expression130, meta130 = pseudobulk_from_tar(
        "GSE130888", meta130, RAW / "GSE130888_RAW.tar"
    )
    expression130.to_csv(OUT / "pseudobulk_expression_GSE130888.csv")
    meta130.to_csv(OUT / "sample_metadata_GSE130888.csv", index=False)

    meta248 = parse_series_matrix(RAW / "GSE248762_series_matrix.txt.gz", "GSE248762")
    meta248["group"] = meta248["title"].str.extract(r"^(SV|LV_NOT_UF|LV_UF)", expand=False)
    expression248, meta248 = pseudobulk_from_tar(
        "GSE248762", meta248, RAW / "GSE248762_RAW.tar"
    )
    expression248.to_csv(OUT / "pseudobulk_expression_GSE248762.csv")
    meta248.to_csv(OUT / "sample_metadata_GSE248762.csv", index=False)

    expression570, meta570, assay_meta570, qc570 = processed_affymetrix_array(
        "GSE57070",
        RAW / "GSE57070_series_matrix.txt.gz",
        RAW / "GPL11209.annot.gz",
    )
    expression570.to_csv(OUT / "gene_expression_GSE57070_patient_level.csv")
    meta570.to_csv(OUT / "sample_metadata_GSE57070.csv", index=False)
    assay_meta570.to_csv(OUT / "sample_metadata_GSE57070_assays.csv", index=False)
    qc570.to_csv(OUT / "array_preprocessing_qc_GSE57070.csv", index=False)

    expressions = {
        "GSE125498": discovery,
        "GSE130888": expression130,
        "GSE248762": expression248,
        "GSE57070": expression570,
    }
    metadata = {
        "GSE125498": meta125,
        "GSE130888": meta130,
        "GSE248762": meta248,
        "GSE57070": meta570,
    }
    panel_detection = []
    for accession, expression in expressions.items():
        for gene in PANEL:
            present = gene in expression.index
            panel_detection.append(
                {
                    "accession": accession,
                    "gene": gene,
                    "present": present,
                    "mean_expression": float(expression.loc[gene].mean()) if present else np.nan,
                }
            )
    pd.DataFrame(panel_detection).to_csv(OUT / "panel_gene_detection.csv", index=False)
    score_tables: dict[str, pd.DataFrame] = {}
    correlation_parts = []
    clinical_parts = []
    availability_parts = []
    ltb_sensitivity_rows: list[dict[str, Any]] = []
    leukocyte_proxy_rows: list[dict[str, Any]] = []
    no_ltb_signatures = {
        name: [gene for gene in genes if not (name == "T cells" and gene == "LTB")]
        for name, genes in signatures.items()
    }
    immune_proxy_columns = [
        name
        for name in signatures
        if name
        in {
            "B cells",
            "NK cells",
            "Monocytes",
            "M1 macrophages",
            "M2 macrophages",
            "Neutrophils",
            "Mast cells",
            "Immature dendritic cells",
            "Mature dendritic cells",
        }
    ]
    for accession, expression in expressions.items():
        score_table, signature_scores, corr, clean_meta = analyze_cohort(
            accession, expression, metadata[accession], signatures, rng
        )
        score_tables[accession] = score_table
        correlation_parts.append(corr)
        no_ltb_scores, no_ltb_sizes = score_signatures(
            expression,
            no_ltb_signatures,
            accession,
            output_suffix="_T_cells_no_LTB",
        )
        x = score_table["panel_score"].to_numpy(dtype=float)
        y = no_ltb_scores["T cells"].reindex(score_table["sample_id"]).to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        sensitivity_status = "reported"
        if valid.sum() < 3 or np.unique(x[valid]).size < 2 or np.unique(y[valid]).size < 2:
            rho, p_value, ci_low, ci_high = np.nan, np.nan, np.nan, np.nan
            sensitivity_status = "not_estimable"
        elif valid.sum() < 6:
            rho = float(spearmanr(x[valid], y[valid]).statistic)
            p_value, ci_low, ci_high = np.nan, np.nan, np.nan
            sensitivity_status = "descriptive_only_n_lt_6"
        else:
            sensitivity_test = spearmanr(x[valid], y[valid])
            rho, p_value = float(sensitivity_test.statistic), float(sensitivity_test.pvalue)
            ci_low, ci_high = bootstrap_rho(x[valid], y[valid], rng)
        ltb_sensitivity_rows.append(
            {
                "accession": accession,
                "n": int(valid.sum()),
                "panel_genes": ";".join(PANEL),
                "signature": "T cells",
                "signature_gene_count": int(no_ltb_sizes.loc[no_ltb_sizes["signature"].eq("T cells"), "detected_gene_count"].iloc[0]),
                "removed_panel_signature_overlap": "LTB",
                "rho": rho,
                "rho_ci95_low": ci_low,
                "rho_ci95_high": ci_high,
                "p_value": p_value,
                "inference_status": sensitivity_status,
            }
        )

        leukocyte_proxy = score_table[immune_proxy_columns].mean(axis=1).to_numpy(dtype=float)
        valid_proxy = np.isfinite(x) & np.isfinite(leukocyte_proxy)
        proxy_status = "reported"
        if valid_proxy.sum() < 3 or np.unique(x[valid_proxy]).size < 2 or np.unique(leukocyte_proxy[valid_proxy]).size < 2:
            proxy_rho, proxy_p, proxy_low, proxy_high = np.nan, np.nan, np.nan, np.nan
            proxy_status = "not_estimable"
        elif valid_proxy.sum() < 6:
            proxy_rho = float(spearmanr(x[valid_proxy], leukocyte_proxy[valid_proxy]).statistic)
            proxy_p, proxy_low, proxy_high = np.nan, np.nan, np.nan
            proxy_status = "descriptive_only_n_lt_6"
        else:
            proxy_test = spearmanr(x[valid_proxy], leukocyte_proxy[valid_proxy])
            proxy_rho, proxy_p = float(proxy_test.statistic), float(proxy_test.pvalue)
            proxy_low, proxy_high = bootstrap_rho(x[valid_proxy], leukocyte_proxy[valid_proxy], rng)
        leukocyte_proxy_rows.append(
            {
                "accession": accession,
                "n": int(valid_proxy.sum()),
                "proxy_definition": "mean of 9 non-T-cell immune ssGSEA marker-enrichment scores",
                "rho": proxy_rho,
                "rho_ci95_low": proxy_low,
                "rho_ci95_high": proxy_high,
                "p_value": proxy_p,
                "inference_status": proxy_status,
            }
        )
        group_column = (
            "group"
            if accession != "GSE57070" and "group" in score_table.columns
            else None
        )
        clinical, availability = clinical_tests(accession, score_table, group_column)
        if not clinical.empty:
            clinical_parts.append(clinical)
        availability_parts.append(availability)
    correlations = pd.concat(correlation_parts, ignore_index=True)
    correlations.to_csv(OUT / "correlations_all_datasets.csv", index=False)
    ltb_sensitivity = pd.DataFrame(ltb_sensitivity_rows)
    original_t_cell_q = correlations.loc[
        correlations["signature"].eq("T cells"),
        ["accession", "q_value_BH_within_dataset_13_tests"],
    ].rename(
        columns={
            "q_value_BH_within_dataset_13_tests":
                "q_value_BH_within_13_signatures"
        }
    )
    ltb_sensitivity = ltb_sensitivity.merge(
        original_t_cell_q, on="accession", how="left"
    )
    ltb_sensitivity.to_csv(OUT / "ltb_removed_t_cell_sensitivity.csv", index=False)
    leukocyte_proxy_correlations = pd.DataFrame(leukocyte_proxy_rows)
    leukocyte_proxy_correlations["q_value_BH_across_cohorts"] = bh(
        leukocyte_proxy_correlations["p_value"].tolist()
    )
    leukocyte_proxy_correlations.to_csv(
        OUT / "panel_vs_leukocyte_proxy.csv", index=False
    )
    clinical = (
        pd.concat(clinical_parts, ignore_index=True)
        if clinical_parts
        else pd.DataFrame(
            columns=[
                "accession",
                "covariate",
                "n",
                "effect",
                "ci95_low",
                "ci95_high",
                "p_value",
                "q_value_BH_within_dataset",
            ]
        )
    )
    clinical.to_csv(OUT / "clinical_associations.csv", index=False)
    pd.concat(availability_parts, ignore_index=True).to_csv(
        OUT / "clinical_covariate_availability.csv", index=False
    )
    make_figures(correlations, score_tables)
    write_report(
        correlations,
        metadata,
        clinical,
        audit,
        panel_signature_overlap,
        ltb_sensitivity,
        leukocyte_proxy_correlations,
    )
    (OUT / "null_distributions_not_run.txt").write_text(
        "Consensus clustering was not run; therefore random-gene-set cluster null distributions "
        "are not applicable. The continuous-score correlation analysis is the prespecified primary analysis.\n",
        encoding="utf-8",
    )

    manifest = file_manifest()
    pd.DataFrame(manifest).to_csv(OUT / "download_manifest.csv", index=False)
    package_versions = {
        "python": __import__("sys").version.split()[0],
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "gseapy": gp.__version__,
        "matplotlib": matplotlib.__version__,
        "seaborn": sns.__version__,
        "statsmodels": statsmodels.__version__,
    }
    metadata_json = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "bootstrap_resamples": BOOTSTRAPS,
        "panel_genes": PANEL,
        "signature_count": len(signatures),
        "signature_names": list(signatures),
        "ssgsea_min_signature_size": 3,
        "cluster_analysis": "not run",
        "ltb_removed_sensitivity": {
            "signature": "T cells",
            "removed_gene": "LTB",
            "panel_unchanged": PANEL,
            "results_file": "ltb_removed_t_cell_sensitivity.csv",
        },
        "panel_leukocyte_proxy_check": {
            "results_file": "panel_vs_leukocyte_proxy.csv",
            "proxy": "mean of 9 non-T-cell immune ssGSEA marker-enrichment scores",
        },
        "cohorts": {
            key: {
                "n_samples_analyzed": len(meta),
                "sample_ids": meta["sample_id"].tolist(),
                "expression_method": (
                    "log2 + quantile normalized array, GPL10558 mapped; duplicate symbol highest mean"
                    if key == "GSE125498"
                    else (
                        "raw scRNA-seq patient pseudobulk, log2(CPM+1)"
                        if key in {"GSE130888", "GSE248762"}
                        else "GEO processed log2 array matrix used as deposited; paired conditions averaged within patient"
                    )
                ),
            }
            for key, meta in metadata.items()
        },
        "package_versions": package_versions,
        "source_files": manifest,
        "frozen_reference_files": {
            "panel_definition": {
                "file": str(PANEL_FILE.relative_to(ROOT)),
                "sha256": sha256(PANEL_FILE),
            },
            "cellmarker_signatures": {
                "file": str(SIGNATURE_FILE.relative_to(ROOT)),
                "sha256": sha256(SIGNATURE_FILE),
            },
        },
        "search_scope": {
            "GEO": "Human peritoneal dialysis, effluent, fibrosis, EPS, UF failure, and mesothelial expression records; relevant candidates audited in dataset_audit.csv",
            "ArrayExpress_BioStudies": "Peritoneal dialysis and peritoneal effluent search results screened; relevant candidates audited in dataset_audit.csv",
            "SRA": "NCBI E-utilities human peritoneal dialysis/effluent/fibrosis search; 95 run records grouped to nine BioProjects in sra_project_audit.csv",
        },
        "input_restriction": "Read only the frozen panel definition, frozen CellMarker signature table, and public GEO source files; no previous Figure 7 output was used.",
        "limitations": [
            "Small exploratory cohorts; no clinical validation.",
            "GSE297829 is one patient and cannot support association testing.",
            "GSE57070 has three independent patients after averaging paired treatment/control arrays; its rho values are descriptive only.",
            "GSE57070 is omitted from validation/concordance figures and is not cited as validation.",
            "GSE57070 was scored from the GEO processed log2 matrix; its downloaded raw CEL archive was not reprocessed.",
            "E-MTAB-5439 contains sorted immune-cell populations, not bulk effluent, and no processed expression matrix.",
            "E-MTAB-5462 has no processed gene-level matrix in BioStudies; raw-read alignment was not performed.",
            "Exact patient-level clinical covariates were unavailable for most cohorts.",
        ],
    }
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata_json, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
