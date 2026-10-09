#!/usr/bin/env python
"""Post hoc / exploratory ECM1 raw-count localization and concentration audit."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import tarfile
import tempfile
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd


MARKER_GENES = (
    "MSLN",
    "UPK3B",
    "S100A8",
    "S100A9",
    "MPO",
    "FCGR3B",
    "LCN2",
    "CSF3R",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    script_dir = (
        root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE248762"
        / "scripts"
    )
    sys.path.insert(0, str(script_dir))
    import gse248762_step2_4_scanpy as pipeline

    cohort_dir = root / "results" / "revision" / "effluent_validation" / "GSE248762"
    archive_path = cohort_dir / "GSE248762_RAW.tar"
    metadata = pipeline.load_sample_metadata(root)
    saved_qc = pd.read_csv(
        cohort_dir / "analysis_20261009" / "sample_qc.csv"
    ).set_index("sample_id")["cells_after_qc"]
    saved_celltypes = pd.read_csv(
        cohort_dir / "analysis_20261009" / "celltype_pseudobulk_metadata.csv"
    )
    saved_celltypes = saved_celltypes.set_index(["sample_id", "cell_type"])
    saved_pseudobulk = pd.read_csv(
        cohort_dir / "analysis_20261009" / "all_cell_pseudobulk_counts.csv",
        index_col=0,
    )

    celltype_rows: list[dict] = []
    marker_rows: list[dict] = []
    distribution_rows: list[dict] = []
    celltypes = [*pipeline.CELLTYPE_MARKERS, "other"]

    with tempfile.TemporaryDirectory(prefix="gse248762_posthoc_ecm1_") as temp:
        scratch = Path(temp)
        with tarfile.open(archive_path, mode="r") as archive:
            for _, sample_meta in metadata.iterrows():
                sample_id = str(sample_meta["sample_id"])
                sample_dir = pipeline.extract_10x_sample(
                    archive, sample_id, scratch
                )
                sample = pipeline.read_10x(sample_dir)
                sample, _, cells_after = pipeline.qc_and_normalize(sample)
                if cells_after != int(saved_qc.loc[sample_id]):
                    raise ValueError(
                        f"{sample_id}: QC-retained cell count differs from "
                        "the saved primary analysis"
                    )
                annotations, _ = pipeline.annotate_celltypes(sample)
                raw = sample.layers["counts"].tocsr()
                symbols = sample.var["gene_symbols"].str.upper().to_numpy()

                symbol_indices: dict[str, np.ndarray] = {}
                for symbol in ("ECM1", *MARKER_GENES):
                    symbol_indices[symbol] = np.flatnonzero(symbols == symbol)
                if len(symbol_indices["ECM1"]) != 1:
                    raise ValueError(
                        f"{sample_id}: expected one ECM1 feature, found "
                        f"{len(symbol_indices['ECM1'])}"
                    )
                ecm1_index = int(symbol_indices["ECM1"][0])
                ecm1_counts = np.asarray(raw[:, ecm1_index].toarray()).ravel()
                ecm1_positive = ecm1_counts > 0
                positive_count = int(ecm1_positive.sum())
                sorted_counts = np.sort(ecm1_counts)[::-1]
                top_1_n = max(1, int(np.ceil(len(sorted_counts) * 0.01)))
                top_10_n = max(1, int(np.ceil(len(sorted_counts) * 0.10)))
                total_ecm1 = int(ecm1_counts.sum())

                distribution_rows.append(
                    {
                        "sample_id": sample_id,
                        "original_group": sample_meta["group_from_metadata"],
                        "primary_group": (
                            "SV"
                            if sample_meta["group_from_metadata"] == "SV"
                            else "LV"
                        ),
                        "qc_cells": int(sample.n_obs),
                        "ecm1_total_raw_counts": total_ecm1,
                        "ecm1_positive_cells": positive_count,
                        "ecm1_positive_fraction": positive_count / sample.n_obs,
                        "ecm1_max_count_per_cell": int(sorted_counts[0]),
                        "ecm1_median_count_per_cell": float(
                            np.median(ecm1_counts)
                        ),
                        "ecm1_p99_count_per_cell": float(
                            np.quantile(ecm1_counts, 0.99)
                        ),
                        "top_1pct_cells_ecm1_count_share": (
                            float(sorted_counts[:top_1_n].sum() / total_ecm1)
                            if total_ecm1
                            else np.nan
                        ),
                        "top_10pct_cells_ecm1_count_share": (
                            float(sorted_counts[:top_10_n].sum() / total_ecm1)
                            if total_ecm1
                            else np.nan
                        ),
                    }
                )

                for celltype in celltypes:
                    cell_mask = annotations == celltype
                    n_cells = int(cell_mask.sum())
                    cell_counts = ecm1_counts[cell_mask]
                    ecm1_count_type = int(cell_counts.sum())
                    ecm1_positive_type = int((cell_counts > 0).sum())
                    celltype_rows.append(
                        {
                            "sample_id": sample_id,
                            "original_group": sample_meta[
                                "group_from_metadata"
                            ],
                            "primary_group": (
                                "SV"
                                if sample_meta["group_from_metadata"] == "SV"
                                else "LV"
                            ),
                            "cell_type": celltype,
                            "qc_cells": n_cells,
                            "ecm1_raw_count_sum": ecm1_count_type,
                            "ecm1_positive_cells": ecm1_positive_type,
                            "ecm1_positive_fraction": (
                                ecm1_positive_type / n_cells if n_cells else 0.0
                            ),
                            "ecm1_mean_raw_count_per_cell": (
                                ecm1_count_type / n_cells if n_cells else 0.0
                            ),
                        }
                    )
                    if n_cells >= pipeline.MIN_PSEUDOBULK_CELLS:
                        saved_n = saved_celltypes.loc[
                            (sample_id, celltype), "n_qc_cells"
                        ]
                        if int(saved_n) != n_cells:
                            raise ValueError(
                                f"{sample_id}/{celltype}: cell count differs "
                                "from saved pseudobulk metadata"
                            )

                    for marker in MARKER_GENES:
                        indices = symbol_indices[marker]
                        if len(indices):
                            marker_counts = np.asarray(
                                raw[cell_mask, :][:, indices].sum(axis=1)
                            ).ravel()
                            marker_positive = marker_counts > 0
                            both_positive = marker_positive & (cell_counts > 0)
                            count_value: int | None = int(marker_counts.sum())
                            n_positive: int | None = int(marker_positive.sum())
                            n_both: int | None = int(both_positive.sum())
                        else:
                            count_value = None
                            n_positive = None
                            n_both = None
                        marker_rows.append(
                            {
                                "sample_id": sample_id,
                                "original_group": sample_meta[
                                    "group_from_metadata"
                                ],
                                "cell_type": celltype,
                                "marker": marker,
                                "marker_detected_in_feature_table": bool(
                                    len(indices)
                                ),
                                "qc_cells": n_cells,
                                "marker_raw_count_sum": count_value,
                                "marker_positive_cells": n_positive,
                                "ecm1_marker_double_positive_cells": n_both,
                            }
                        )

                if int(ecm1_counts.sum()) != int(
                    saved_pseudobulk.loc[
                        sample.var_names[ecm1_index], sample_id
                    ]
                ):
                    raise ValueError(
                        f"{sample_id}: ECM1 raw-count sum differs from "
                        "the saved all-cell pseudobulk matrix"
                    )
                del sample
                shutil.rmtree(sample_dir)

    celltype_table = pd.DataFrame(celltype_rows)
    marker_table = pd.DataFrame(marker_rows)
    distribution_table = pd.DataFrame(distribution_rows)
    celltype_table.to_csv(output / "ecm1_posthoc_by_sample_celltype.csv", index=False)
    marker_table.to_csv(output / "ecm1_posthoc_marker_comparison.csv", index=False)
    distribution_table.to_csv(output / "ecm1_posthoc_cell_concentration.csv", index=False)

    config = {
        "analysis_label": "post hoc / exploratory",
        "accession": "GSE248762",
        "input": str(archive_path.relative_to(root)),
        "input_bytes": archive_path.stat().st_size,
        "python_version": sys.version.split()[0],
        "scanpy_version": version("scanpy"),
        "anndata_version": version("anndata"),
        "qc_and_annotation": (
            "Unchanged functions and thresholds from the saved primary "
            "GSE248762 Scanpy pipeline; all 16 GEO samples retained."
        ),
        "qc_minimum_detected_genes": pipeline.MIN_GENES,
        "qc_maximum_mitochondrial_percent": pipeline.MAX_MT_PERCENT,
        "celltype_pseudobulk_minimum_cells": pipeline.MIN_PSEUDOBULK_CELLS,
        "marker_comparison": list(MARKER_GENES),
        "seed": 42,
        "ambient_rna_limitation": (
            "Filtered-cell matrices do not include empty-droplet/background "
            "profiles; descriptive localization is not an ambient-RNA test."
        ),
    }
    (output / "ecm1_posthoc_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    print(
        distribution_table.loc[
            distribution_table["sample_id"] == "GSM7919584"
        ].to_string(index=False)
    )
    print(
        celltype_table.loc[
            celltype_table["sample_id"] == "GSM7919584"
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
