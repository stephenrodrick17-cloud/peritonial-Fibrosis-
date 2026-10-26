from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[3]
RAW_ROOT = ROOT / "data" / "raw" / "GSE248762_extracted"
SINGLE_CELL_ROOT = ROOT / "results" / "revision" / "single_cell_context"
OUTPUT_ROOT = ROOT / "results" / "revision" / "sc7_fc_audit"


def main() -> None:
    qc_cells = pd.read_csv(
        SINGLE_CELL_ROOT / "cell_type_annotation_table.csv",
        usecols=["cell_id", "sample_id"],
        dtype={"cell_id": str, "sample_id": str},
    )
    qc_by_sample = {
        sample_id: set(group["cell_id"])
        for sample_id, group in qc_cells.groupby("sample_id", sort=False)
    }
    sample_audit = pd.read_csv(
        SINGLE_CELL_ROOT / "sample_group_qc_audit.csv",
        dtype={"sample_id": str, "group_from_metadata": str},
    ).set_index("sample_id")

    summary_rows: list[dict[str, object]] = []
    barcode_rows: list[dict[str, object]] = []
    matrix_files = sorted(RAW_ROOT.glob("*.matrix.mtx.gz"))
    if len(matrix_files) != 16:
        raise ValueError(f"Expected 16 sample matrices; found {len(matrix_files)}")

    for matrix_path in matrix_files:
        sample_id = matrix_path.name.split("_", maxsplit=1)[0]
        prefix = matrix_path.name.removesuffix(".matrix.mtx.gz")
        feature_path = RAW_ROOT / f"{prefix}.features.tsv.gz"
        barcode_path = RAW_ROOT / f"{prefix}.barcodes.tsv.gz"

        features = pd.read_csv(feature_path, sep="\t", header=None, dtype=str)
        symbols = features.iloc[:, 1].fillna("").str.upper()
        gene_rows = symbols.index[symbols.eq("ECM1")].to_list()
        if not gene_rows:
            raise ValueError(f"ECM1 feature is absent from {feature_path.name}")

        with gzip.open(barcode_path, "rt") as handle:
            barcodes = [line.rstrip("\r\n") for line in handle]
        with gzip.open(matrix_path, "rb") as handle:
            matrix = mmread(handle).tocsr()
        if matrix.shape != (len(features), len(barcodes)):
            raise ValueError(f"Feature/barcode dimensions do not match {prefix}")

        qc_ids = qc_by_sample[sample_id]
        retained = [
            index
            for index, barcode in enumerate(barcodes)
            if f"{sample_id}_{barcode}" in qc_ids
        ]
        excluded = [
            index
            for index, barcode in enumerate(barcodes)
            if f"{sample_id}_{barcode}" not in qc_ids
        ]
        audit = sample_audit.loc[sample_id]
        expected_before = int(audit["cells_before_qc"])
        expected_after = int(audit["cells_after_qc_recomputed"])
        if len(barcodes) != expected_before or len(retained) != expected_after:
            raise ValueError(
                f"QC barcode reconciliation failed for {sample_id}: "
                f"matrix={len(barcodes)}, retained={len(retained)}, "
                f"saved={expected_before}/{expected_after}"
            )

        gene_counts = matrix[gene_rows, :].sum(axis=0).A1
        library_sizes = matrix.sum(axis=0).A1
        excluded_gene_counts = gene_counts[excluded]
        excluded_libraries = library_sizes[excluded]
        cpm = (
            excluded_gene_counts / excluded_libraries * 1_000_000
            if len(excluded)
            else pd.Series(dtype=float)
        )
        original_group = str(audit["group_from_metadata"])

        summary_rows.append(
            {
                "sample_id": sample_id,
                "original_group": original_group,
                "raw_matrix_barcodes": len(barcodes),
                "qc_retained_barcodes": len(retained),
                "qc_excluded_called_barcodes": len(excluded),
                "ecm1_positive_excluded_barcodes": int(
                    (excluded_gene_counts > 0).sum()
                ),
                "ecm1_positive_percent_excluded": (
                    float((excluded_gene_counts > 0).mean() * 100)
                    if len(excluded)
                    else 0.0
                ),
                "ecm1_mean_raw_count_excluded": (
                    float(excluded_gene_counts.mean()) if len(excluded) else 0.0
                ),
                "ecm1_median_raw_count_excluded": (
                    float(pd.Series(excluded_gene_counts).median())
                    if len(excluded)
                    else 0.0
                ),
                "ecm1_mean_cpm_excluded": (
                    float(cpm.mean()) if len(excluded) else 0.0
                ),
                "ecm1_median_cpm_excluded": (
                    float(pd.Series(cpm).median()) if len(excluded) else 0.0
                ),
                "ecm1_max_raw_count_excluded": (
                    int(excluded_gene_counts.max()) if len(excluded) else 0
                ),
            }
        )
        barcode_rows.extend(
            {
                "sample_id": sample_id,
                "original_group": original_group,
                "barcode": barcodes[index],
                "ecm1_raw_count": int(gene_counts[index]),
                "total_library_size": int(library_sizes[index]),
                "ecm1_cpm": float(
                    gene_counts[index] / library_sizes[index] * 1_000_000
                )
                if library_sizes[index]
                else 0.0,
            }
            for index in excluded
        )
        print(
            f"{sample_id}: {len(excluded)} QC-excluded called barcodes; "
            f"{int((excluded_gene_counts > 0).sum())} with nonzero counts",
            flush=True,
        )

    pd.DataFrame(summary_rows).to_csv(
        OUTPUT_ROOT / "sc7_ecm1_pre_qc_excluded_summary.csv", index=False
    )
    pd.DataFrame(barcode_rows).to_csv(
        OUTPUT_ROOT / "sc7_ecm1_pre_qc_excluded_barcodes.csv", index=False
    )


if __name__ == "__main__":
    main()
