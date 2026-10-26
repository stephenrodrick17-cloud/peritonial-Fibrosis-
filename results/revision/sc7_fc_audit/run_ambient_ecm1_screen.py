from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SINGLE_CELL = ROOT / "results" / "revision" / "single_cell_context"


def bh_adjust(values: list[float]) -> list[float]:
    p_values = np.asarray(values, dtype=float)
    order = np.argsort(p_values)
    adjusted_sorted = np.minimum.accumulate(
        (p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1))[::-1]
    )[::-1]
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return adjusted.tolist()


def main() -> None:
    expression = pd.read_csv(
        SINGLE_CELL / "hub_gene_expression_per_cell.csv",
        usecols=["sample_id", "original_group", "lineage", "ECM1"],
    )
    expression["ECM1"] = pd.to_numeric(expression["ECM1"], errors="raise")
    cell_summary = (
        expression.groupby(["sample_id", "original_group", "lineage"], as_index=False)
        .agg(
            n_cells=("ECM1", "size"),
            n_expressing_cells=("ECM1", lambda values: int((values > 0).sum())),
            percent_expressing_cells=("ECM1", lambda values: 100.0 * (values > 0).mean()),
            mean_log1p_CPM10k=("ECM1", "mean"),
        )
    )
    cell_summary.to_csv(OUT / "sc7_ecm1_cell_expression_by_sample_lineage.csv", index=False)

    saved_sample_values = pd.read_csv(OUT / "sc7_ecm1_subgroup_sample_values.csv")
    ecm1_cpm = saved_sample_values.copy()
    cpm_wide = ecm1_cpm.pivot(index=["sample_id", "original_group"], columns="lineage", values="CPM")
    cells_wide = ecm1_cpm.pivot(
        index=["sample_id", "original_group"],
        columns="lineage",
        values="n_qc_cells_in_lineage",
    )
    excluded_wide = ecm1_cpm.pivot(
        index=["sample_id", "original_group"],
        columns="lineage",
        values="excluded_by_20_cell_rule",
    )
    cpm_wide.columns = [f"{name}_CPM" for name in cpm_wide.columns]
    cells_wide.columns = [f"{name}_n_qc_cells" for name in cells_wide.columns]
    excluded_wide.columns = [f"{name}_excluded_by_20_cell_rule" for name in excluded_wide.columns]
    sample_matrix = cpm_wide.join(cells_wide).join(excluded_wide).reset_index()
    sample_matrix.to_csv(OUT / "sc7_ecm1_pseudobulk_cpm_by_sample.csv", index=False)

    count_table = pd.read_csv(
        SINGLE_CELL / "pseudobulk_counts" / "counts_myeloid.csv",
        usecols=lambda name: name in {"feature_id", "gene_symbol"} or name.startswith("GSM"),
    )
    sample_ids = [column for column in count_table.columns if column.startswith("GSM")]
    symbols = count_table["gene_symbol"].fillna("").astype(str).str.upper()
    symbols = symbols.where(symbols.ne(""), count_table["feature_id"].astype(str))
    counts = count_table[sample_ids].apply(pd.to_numeric, errors="raise")
    collapsed = counts.groupby(symbols, sort=False).sum()
    library_sizes = counts.sum(axis=0)
    cpm = collapsed.div(library_sizes, axis=1) * 1_000_000

    group_table = pd.read_csv(SINGLE_CELL / "sample_group_qc_audit.csv")
    group_by_sample = group_table.set_index("sample_id")["group_from_metadata"].astype(str)
    group_by_sample = group_by_sample.reindex(sample_ids)
    if group_by_sample.isna().any():
        raise ValueError("Some myeloid pseudobulk samples lack saved group metadata.")
    uf_samples = group_by_sample.index[group_by_sample.eq("LV_UF")].tolist()
    other_samples = group_by_sample.index[group_by_sample.ne("LV_UF")].tolist()
    if len(uf_samples) != 4 or len(other_samples) != 12:
        raise ValueError("Expected four LV_UF and twelve comparison samples.")
    marker_genes = ["ECM1", "COL1A1", "COL1A2", "DCN", "LUM", "MSLN", "KRT8", "KRT18"]
    missing = [gene for gene in marker_genes if gene not in cpm.index]
    if missing:
        raise ValueError(f"Marker symbols missing from myeloid matrix: {missing}")

    sample_rows = []
    summary_rows = []
    allocations = list(itertools.combinations(sample_ids, len(uf_samples)))
    for gene in marker_genes:
        gene_cpm = cpm.loc[gene]
        uf_values = gene_cpm.loc[uf_samples]
        other_values = gene_cpm.loc[other_samples]
        observed_difference = float(uf_values.mean() - other_values.mean())
        extreme_count = 0
        for selected in allocations:
            selected_set = set(selected)
            selected_values = gene_cpm.loc[list(selected)]
            complement = gene_cpm.loc[[sample for sample in sample_ids if sample not in selected_set]]
            if float(selected_values.mean() - complement.mean()) >= observed_difference - 1e-12:
                extreme_count += 1
        p_value = extreme_count / len(allocations)
        summary_rows.append(
            {
                "gene": gene,
                "n_LV_UF": len(uf_values),
                "n_other": len(other_values),
                "mean_CPM_LV_UF": float(uf_values.mean()),
                "mean_CPM_other": float(other_values.mean()),
                "median_CPM_LV_UF": float(uf_values.median()),
                "median_CPM_other": float(other_values.median()),
                "mean_difference_LV_UF_minus_other": observed_difference,
                "one_sided_exact_permutation_P": p_value,
                "n_label_allocations": len(allocations),
            }
        )
        for sample_id in sample_ids:
            sample_rows.append(
                {
                    "gene": gene,
                    "sample_id": sample_id,
                    "original_group": group_by_sample.loc[sample_id],
                    "myeloid_library_size": int(library_sizes.loc[sample_id]),
                    "myeloid_CPM": float(gene_cpm.loc[sample_id]),
                }
            )
    q_values = bh_adjust([row["one_sided_exact_permutation_P"] for row in summary_rows])
    for row, q_value in zip(summary_rows, q_values, strict=True):
        row["BH_q_over_eight_screened_genes"] = q_value
    pd.DataFrame(summary_rows).to_csv(OUT / "sc7_ecm1_ambient_marker_screen.csv", index=False)
    pd.DataFrame(sample_rows).to_csv(OUT / "sc7_ecm1_ambient_marker_sample_cpm.csv", index=False)

    print(
        f"Wrote sample-by-lineage ECM1 expression for {len(cell_summary)} sample-lineage "
        f"groups; tested {len(marker_genes)} myeloid markers across "
        f"{len(allocations)} exact LV_UF label allocations."
    )


if __name__ == "__main__":
    main()
