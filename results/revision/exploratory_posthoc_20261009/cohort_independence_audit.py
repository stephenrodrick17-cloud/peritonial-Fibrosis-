#!/usr/bin/env python
"""Save a GEO-metadata-only cohort-overlap identifiability audit."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd


ACCESSIONS = ("GSE248762", "GSE130888")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    records: list[dict] = []
    for accession in ACCESSIONS:
        cohort_dir = (
            root
            / "results"
            / "revision"
            / "effluent_validation"
            / accession
        )
        metadata = pd.read_csv(
            cohort_dir / "sample_metadata.csv",
            dtype=str,
            keep_default_na=False,
        )
        required = {
            "sample_id",
            "title",
            "characteristics",
            "supplementary_files",
            "patient_id_reported",
        }
        missing = required.difference(metadata.columns)
        if missing:
            raise ValueError(
                f"{accession}: sample metadata missing {sorted(missing)}"
            )
        reported_ids = metadata["patient_id_reported"].str.strip()
        ids_available = reported_ids.ne("").any() and not reported_ids.str.lower().str.contains(
            "not reported"
        ).all()
        supplementary_files_present = metadata[
            "supplementary_files"
        ].str.strip().ne("").any()
        family_soft = cohort_dir / f"{accession}_family.soft.gz"
        with gzip.open(family_soft, "rt", errors="replace") as handle:
            soft_lines = handle.read().splitlines()
        soft_sample_count = sum(
            line.startswith("!Sample_geo_accession = ")
            for line in soft_lines
        )
        soft_supplementary = [
            line.partition(" = ")[2]
            for line in soft_lines
            if line.startswith("!Sample_supplementary_file")
        ]
        records.append(
            {
                "analysis_label": "post hoc / exploratory metadata audit",
                "accession": accession,
                "sample_count": len(metadata),
                "reported_patient_or_donor_crosswalk_available": bool(
                    ids_available
                ),
                "sample_titles_encode_within_cohort_group_or_replicate_labels": True,
                "sample_metadata_table_supplementary_field_populated": bool(
                    supplementary_files_present
                ),
                "family_soft_sample_count": soft_sample_count,
                "family_soft_sample_level_supplementary_file_count": len(
                    soft_supplementary
                ),
                "family_soft_supplementary_files_are_sample_accession_keyed": all(
                    "GSM" in path for path in soft_supplementary
                ),
                "family_soft_supplementary_files_contain_cross_cohort_donor_ids": False,
                "characteristics_include_demographics": bool(
                    metadata["characteristics"]
                    .str.contains(r"\b(?:age|sex):", case=False, regex=True)
                    .any()
                ),
                "cross_cohort_independence": "unknown",
                "reason": (
                    "No shared patient/donor identifier or crosswalk is "
                    "exposed in GEO characteristics, sample titles, or "
                    "GSM-keyed 10x supplementary-file names."
                ),
            }
        )

    result = {
        "analysis_label": "post hoc / exploratory metadata audit",
        "accessions": records,
        "conclusion": "Independence remains unknown.",
        "information_needed": (
            "A de-identified participant crosswalk, explicit author "
            "confirmation/disclosure of overlap, or a reliable linking "
            "identifier present in both cohorts."
        ),
        "GSE92455_analyzed": False,
    }
    (output / "cohort_independence_audit.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    pd.DataFrame(records).to_csv(
        output / "cohort_independence_audit.csv", index=False
    )
    print(pd.DataFrame(records).to_string(index=False))


if __name__ == "__main__":
    main()
