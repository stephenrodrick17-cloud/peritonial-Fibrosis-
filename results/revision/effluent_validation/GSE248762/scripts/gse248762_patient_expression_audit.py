#!/usr/bin/env python
"""Save per-patient raw-count and detected-cell audit for selected genes."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import tarfile
import tempfile
from pathlib import Path

import pandas as pd

from gse248762_step2_4_scanpy import (
    annotate_celltypes,
    extract_10x_sample,
    load_sample_metadata,
    qc_and_normalize,
    read_10x,
)


AUDIT_GENES = ("TNFSF15", "SERPINA10", "ECM1")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()

    metadata = load_sample_metadata(repo_root)
    archive_path = (
        repo_root
        / "results"
        / "revision"
        / "effluent_validation"
        / "GSE248762"
        / "GSE248762_RAW.tar"
    )
    all_counts = pd.read_csv(
        output_dir / "all_cell_pseudobulk_counts.csv", index_col=0
    )
    expected_cells = pd.read_csv(output_dir / "sample_qc.csv").set_index(
        "sample_id"
    )["cells_after_qc"]
    rows = []
    with tempfile.TemporaryDirectory(prefix="gse248762_audit_") as temp:
        scratch = Path(temp)
        with tarfile.open(archive_path, mode="r") as archive:
            for _, sample_meta in metadata.iterrows():
                sample_id = sample_meta["sample_id"]
                sample_dir = extract_10x_sample(archive, sample_id, scratch)
                sample = read_10x(sample_dir)
                sample, _, _ = qc_and_normalize(sample)
                if sample.n_obs != int(expected_cells.loc[sample_id]):
                    raise ValueError(
                        f"{sample_id}: QC-retained cells differ from main analysis"
                    )
                barcode_hash = hashlib.sha256(
                    "\n".join(sample.obs_names.astype(str)).encode("utf-8")
                ).hexdigest()
                annotations, _ = annotate_celltypes(sample)
                symbols = sample.var["gene_symbols"].astype(str).str.upper()
                raw = sample.layers["counts"].tocsr()
                gene_ids = {}
                for gene in AUDIT_GENES:
                    matched = sample.var_names[symbols.to_numpy() == gene].astype(str)
                    if len(matched) != 1:
                        raise ValueError(
                            f"{sample_id}: {gene} has {len(matched)} feature IDs"
                        )
                    gene_ids[gene] = matched[0]
                    if int(all_counts.loc[matched[0], sample_id]) < 0:
                        raise ValueError("Invalid negative pseudobulk count")

                fibroblast_mask = annotations == "fibroblast"
                ecm1_index = sample.var_names.get_loc(gene_ids["ECM1"])
                rows.append(
                    {
                        "sample_id": sample_id,
                        "title": sample_meta["title"],
                        "group": sample_meta["group_from_metadata"],
                        "primary_group": (
                            "SV"
                            if sample_meta["group_from_metadata"] == "SV"
                            else "LV"
                        ),
                        "patient_id_proxy": sample_id,
                        "qc_cells": sample.n_obs,
                        "raw_pseudobulk_TNFSF15": int(
                            all_counts.loc[gene_ids["TNFSF15"], sample_id]
                        ),
                        "cells_expressing_TNFSF15": int(
                            (raw[:, sample.var_names.get_loc(gene_ids["TNFSF15"])] > 0)
                            .sum()
                        ),
                        "raw_pseudobulk_SERPINA10": int(
                            all_counts.loc[gene_ids["SERPINA10"], sample_id]
                        ),
                        "cells_expressing_SERPINA10": int(
                            (
                                raw[
                                    :,
                                    sample.var_names.get_loc(gene_ids["SERPINA10"]),
                                ]
                                > 0
                            ).sum()
                        ),
                        "raw_pseudobulk_ECM1": int(
                            all_counts.loc[gene_ids["ECM1"], sample_id]
                        ),
                        "fibroblast_cells": int(fibroblast_mask.sum()),
                        "fibroblasts_expressing_ECM1": int(
                            (
                                raw[
                                    fibroblast_mask,
                                    ecm1_index,
                                ]
                                > 0
                            ).sum()
                        ),
                        "all_cells_expressing_ECM1": int(
                            (raw[:, ecm1_index] > 0).sum()
                        ),
                        "retained_barcode_sha256": barcode_hash,
                    }
                )
                del sample
                shutil.rmtree(sample_dir)

    result = pd.DataFrame(rows)
    if set(result["sample_id"]) != set(metadata["sample_id"]):
        raise ValueError("Patient-level audit is missing samples")
    result.to_csv(output_dir / "patient_expression_audit.csv", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
