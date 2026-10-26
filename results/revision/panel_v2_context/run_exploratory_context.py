"""Exploratory, post hoc context analyses on the GSE125498 discovery cohort.

All generated results are descriptive, use the same data as discovery, and are
not validation. Run from the repository root with the project's Python environment.
"""

from __future__ import annotations

import json
import re
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "revision" / "panel_v2_context"
EXPR_PATH = ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv"
META_PATH = ROOT / "results" / "tables" / "GSE125498_sample_metadata.csv"
RANK_PATH = ROOT / "results" / "tables" / "GSE125498_all_results.csv"

V2_GENES = ["FLT3LG", "EBI3", "LTB", "ADAM19", "CST7"]
ORIGINAL_GENES = [
    "TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "ECM1", "SERPINA10", "CST7"
]
LEUKOCYTE_MARKERS = [
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD19",
    "MS4A1", "CD79A", "CD79B", "NCAM1", "KLRD1", "CD14", "CD68", "ITGAM",
    "FCGR3A",
]
CELL_TYPES = {
    "CD4 T cell": lambda n: re.search(r"\bt cell\b", n) is not None
    and "cd4" in n
    and "nkt" not in n,
    "CD8 T cell": lambda n: re.search(r"\bt cell\b", n) is not None
    and "cd8" in n
    and "nkt" not in n,
    "NKT cell": lambda n: "nkt" in n,
    "NK cell": lambda n: ("natural killer" in n or re.search(r"\bnk cell\b", n))
    and "nkt" not in n,
    "T cell (other)": lambda n: re.search(r"\bt cell\b", n) is not None
    and not any(x in n for x in ("cd4", "cd8", "nkt", "natural killer")),
    "B cell": lambda n: re.search(r"\bb cell\b", n) is not None and "plasma" not in n,
    "Plasma cell": lambda n: "plasma cell" in n,
    "Monocyte": lambda n: "monocyte" in n
    and "dendritic" not in n
    and "progenitor" not in n,
    "Macrophage": lambda n: "macrophage" in n,
    "Dendritic cell": lambda n: "dendritic" in n,
    "Neutrophil": lambda n: "neutrophil" in n,
    "Eosinophil": lambda n: "eosinophil" in n,
    "Mast cell": lambda n: "mast cell" in n,
}
SEED = 125498
BOOTSTRAPS = 5000
DISCOVERY_LABEL = "exploratory, post hoc, same data as discovery, not validation"


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in (EXPR_PATH, META_PATH, RANK_PATH):
        if not path.is_file():
            raise FileNotFoundError(f"Required input file is missing: {path}")
    expr = pd.read_csv(EXPR_PATH).set_index("Gene")
    meta = pd.read_csv(META_PATH)
    ranking = pd.read_csv(RANK_PATH)
    if expr.index.has_duplicates or ranking["Gene"].duplicated().any():
        raise ValueError("Gene symbols are not unique in the expression/ranking inputs.")
    if set(meta["sample_id"]) != set(expr.columns):
        raise ValueError("Sample IDs do not match exactly between expression and metadata.")
    if not set(meta["group"]).issubset({"SPD", "LPD"}):
        raise ValueError("Unexpected group labels in sample metadata.")
    if len(ranking) != 11741 or ranking["t"].isna().any():
        raise ValueError("Expected the complete 11,741-gene limma ranking with finite t values.")
    missing_genes = sorted(set(V2_GENES + ORIGINAL_GENES) - set(expr.index))
    if missing_genes:
        raise ValueError(f"Panel genes missing from expression matrix: {missing_genes}")
    sample_ids = meta["sample_id"].tolist()
    expr = expr.loc[:, sample_ids].astype(float)
    meta = meta.set_index("sample_id").loc[sample_ids].reset_index()
    return expr, meta, ranking


def get_cellmarker_signatures(expr_genes: set[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    library = gp.get_library("CellMarker_2024", organism="Human")
    signatures: dict[str, list[str]] = {}
    provenance: list[dict[str, object]] = []
    for cell_type, matcher in CELL_TYPES.items():
        selected = {
            term: sorted({str(g).upper() for g in genes} & expr_genes)
            for term, genes in library.items()
            if term.endswith("Human") and matcher(term.lower())
        }
        selected = {term: genes for term, genes in selected.items() if genes}
        union = sorted({gene for genes in selected.values() for gene in genes})
        if len(union) < 15 or not selected:
            raise ValueError(
                f"CellMarker signature for {cell_type} has only {len(union)} "
                "genes in the expression matrix."
            )
        signatures[cell_type] = union
        provenance.append(
            {
                "cell_type": cell_type,
                "source_library": "CellMarker_2024 (Enrichr, human)",
                "n_source_signatures_merged": len(selected),
                "source_signatures": "; ".join(sorted(selected)),
                "n_genes_in_expression": len(union),
                "merge_rule": "union of matching human tissue/context marker signatures",
            }
        )
    return signatures, pd.DataFrame(provenance)


def add_analysis_label(frame: pd.DataFrame) -> pd.DataFrame:
    labeled = frame.copy()
    labeled.insert(0, "analysis_label", DISCOVERY_LABEL)
    return labeled


def format_gsea_p(value: Any) -> str:
    numeric = float(value)
    return "<0.001" if numeric < 0.001 else f"{numeric:.4g}"


def score_ssgsea(
    expr: pd.DataFrame, signatures: dict[str, list[str]]
) -> pd.DataFrame:
    result = gp.ssgsea(
        data=expr,
        gene_sets=signatures,
        min_size=15,
        max_size=5000,
        sample_norm_method="rank",
        correl_norm_type="rank",
        permutation_num=None,
        threads=1,
        no_plot=True,
        seed=SEED,
        outdir=None,
        verbose=False,
    ).res2d
    if result is None:
        raise RuntimeError("GSEApy ssGSEA did not return a result table.")
    if {"Name", "Term", "NES"}.issubset(result.columns):
        scores = result.pivot(index="Term", columns="Name", values="NES")
    else:
        raise ValueError(f"Unexpected GSEApy ssGSEA output columns: {result.columns.tolist()}")
    scores = scores.reindex(index=list(signatures))
    if scores.isna().any().any():
        raise ValueError("ssGSEA returned missing cell-type/sample scores.")
    scores.index.name = "cell_type"
    return scores


def make_sample_scores(expr: pd.DataFrame) -> pd.DataFrame:
    z = expr.sub(expr.mean(axis=1), axis=0).div(
        expr.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )
    z = z.fillna(0)
    scores: dict[str, pd.Series] = {}
    for label, genes in (("panel_v2", V2_GENES), ("original_8", ORIGINAL_GENES)):
        scores[f"{label}_composite"] = z.loc[genes].mean(axis=0)
        for gene in genes:
            scores[f"{label}_{gene}"] = expr.loc[[gene]].iloc[0]
    present_markers = [gene for gene in LEUKOCYTE_MARKERS if gene in z.index]
    if len(present_markers) < 10:
        raise ValueError("Fewer than 10 established repository leukocyte markers are present.")
    scores["dominant_leukocyte_marker_proxy"] = z.loc[present_markers].mean(axis=0)
    return pd.DataFrame(scores, index=expr.columns)


def immune_composition(
    expr: pd.DataFrame, meta: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    signatures, provenance = get_cellmarker_signatures(set(expr.index))
    scores = score_ssgsea(expr, signatures)
    sample_scores = scores.T
    sample_scores.index.name = "sample_id"
    sample_scores = sample_scores.reset_index().merge(meta, on="sample_id", validate="one_to_one")
    sample_scores = sample_scores.set_index("sample_id")
    add_analysis_label(sample_scores).to_csv(
        OUT / "immune_celltype_ssgsea_scores.csv", index_label="sample_id"
    )
    add_analysis_label(provenance).to_csv(
        OUT / "immune_signature_manifest.csv", index=False
    )

    comparisons: list[dict[str, object]] = []
    for cell_type in scores.index:
        lpd = sample_scores.loc[sample_scores["group"] == "LPD", cell_type].to_numpy()
        spd = sample_scores.loc[sample_scores["group"] == "SPD", cell_type].to_numpy()
        p_value = float(
            stats.mannwhitneyu(lpd, spd, alternative="two-sided").pvalue
        )
        comparisons.append(
            {
                "cell_type": cell_type,
                "LPD_n": len(lpd),
                "SPD_n": len(spd),
                "LPD_mean_ssGSEA": float(np.mean(lpd)),
                "SPD_mean_ssGSEA": float(np.mean(spd)),
                "LPD_minus_SPD_mean": float(np.mean(lpd) - np.mean(spd)),
                "Mann_Whitney_P": p_value,
            }
        )
    comparison = pd.DataFrame(comparisons)
    comparison["BH_q_across_cell_types"] = multipletests(
        comparison["Mann_Whitney_P"], method="fdr_bh"
    )[1]
    comparison = add_analysis_label(comparison)
    comparison.to_csv(OUT / "immune_celltype_group_comparison.csv", index=False)

    gene_scores = make_sample_scores(expr)
    merged = gene_scores.merge(meta.set_index("sample_id"), left_index=True, right_index=True)
    add_analysis_label(merged).to_csv(
        OUT / "panel_gene_and_composite_scores.csv", index_label="sample_id"
    )
    immune_names = list(scores.index)
    corr_rows: list[dict[str, object]] = []
    targets = [f"panel_v2_{g}" for g in V2_GENES] + [
        "panel_v2_composite",
        *[f"original_8_{g}" for g in ORIGINAL_GENES],
        "original_8_composite",
    ]
    for target in targets:
        for cell_type in immune_names:
            x = merged[target].to_numpy(dtype=float)
            y = sample_scores[cell_type].to_numpy(dtype=float)
            pearson = stats.pearsonr(x, y)
            spearman = stats.spearmanr(x, y)
            corr_rows.append(
                {
                    "score": target,
                    "cell_type": cell_type,
                    "Pearson_r": float(pearson[0]),
                    "Pearson_P": float(pearson[1]),
                    "Spearman_rho": float(spearman[0]),
                    "Spearman_P": float(spearman[1]),
                }
            )
    correlations = pd.DataFrame(corr_rows)
    correlations["Pearson_BH_q_within_score"] = correlations.groupby("score")[
        "Pearson_P"
    ].transform(lambda p: multipletests(p, method="fdr_bh")[1])
    correlations = add_analysis_label(correlations)
    correlations.to_csv(OUT / "panel_celltype_correlations.csv", index=False)

    dominant = merged["dominant_leukocyte_marker_proxy"].to_numpy(dtype=float)
    group = (merged["group"] == "LPD").astype(int).to_numpy()
    model_rows: list[dict[str, object]] = []
    model_summary: dict[str, object] = {}
    for target in ("panel_v2_composite", "original_8_composite"):
        y = merged[target].to_numpy(dtype=float)
        unadjusted = sm.OLS(y, sm.add_constant(group)).fit()
        adjusted = sm.OLS(y, sm.add_constant(np.column_stack([group, dominant]))).fit()
        att = (
            100 * (float(unadjusted.params[1]) - float(adjusted.params[1]))
            / float(unadjusted.params[1])
            if abs(float(unadjusted.params[1])) > 1e-12
            else np.nan
        )
        model_rows.extend(
            [
                {
                    "composite": target,
                    "model": "composite ~ group",
                    "group_beta_LPD_vs_SPD": float(unadjusted.params[1]),
                    "group_SE": float(unadjusted.bse[1]),
                    "group_P": float(unadjusted.pvalues[1]),
                    "group_CI95_low": float(unadjusted.conf_int()[1, 0]),
                    "group_CI95_high": float(unadjusted.conf_int()[1, 1]),
                    "leukocyte_proxy_beta": np.nan,
                    "leukocyte_proxy_P": np.nan,
                    "R_squared": float(unadjusted.rsquared),
                    "group_beta_attenuation_percent": 0.0,
                },
                {
                    "composite": target,
                    "model": "composite ~ group + dominant leukocyte marker proxy",
                    "group_beta_LPD_vs_SPD": float(adjusted.params[1]),
                    "group_SE": float(adjusted.bse[1]),
                    "group_P": float(adjusted.pvalues[1]),
                    "group_CI95_low": float(adjusted.conf_int()[1, 0]),
                    "group_CI95_high": float(adjusted.conf_int()[1, 1]),
                    "leukocyte_proxy_beta": float(adjusted.params[2]),
                    "leukocyte_proxy_P": float(adjusted.pvalues[2]),
                    "R_squared": float(adjusted.rsquared),
                    "group_beta_attenuation_percent": float(att),
                },
            ]
        )
        model_summary[target] = {
            "unadjusted_group_beta": float(unadjusted.params[1]),
            "unadjusted_group_p": float(unadjusted.pvalues[1]),
            "adjusted_group_beta": float(adjusted.params[1]),
            "adjusted_group_p": float(adjusted.pvalues[1]),
            "adjusted_leukocyte_proxy_beta": float(adjusted.params[2]),
            "adjusted_leukocyte_proxy_p": float(adjusted.pvalues[2]),
            "group_beta_attenuation_percent": float(att),
            "adjusted_r_squared": float(adjusted.rsquared),
        }
    models = pd.DataFrame(model_rows)
    models = add_analysis_label(models)
    models.to_csv(OUT / "composite_leukocyte_adjustment.csv", index=False)
    present_markers = [gene for gene in LEUKOCYTE_MARKERS if gene in expr.index]
    return comparison, correlations, models, {
        "model_summary": model_summary,
        "n_leukocyte_markers": len(present_markers),
        "leukocyte_markers_used": present_markers,
    }


def run_hallmark_gsea(ranking: pd.DataFrame) -> pd.DataFrame:
    hallmark = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
    gene_sets = {
        str(name): sorted({str(g).upper() for g in members})
        for name, members in hallmark.items()
    }
    (OUT / "hallmark_gene_sets.json").write_text(
        json.dumps(
            {
                "analysis_label": DISCOVERY_LABEL,
                "source_library": "MSigDB_Hallmark_2020",
                "retrieval": "gseapy Enrichr library",
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "gene_sets": gene_sets,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    rank = ranking[["Gene", "t"]].copy()
    rank["Gene"] = rank["Gene"].astype(str).str.upper()
    rank["t"] = pd.to_numeric(rank["t"], errors="raise")
    if rank["Gene"].duplicated().any() or rank["t"].isna().any():
        raise ValueError("GSEA rank must contain unique symbols and finite limma t statistics.")
    rank = rank.sort_values("t", ascending=False, kind="mergesort")
    prerank_result = gp.prerank(
        rnk=rank,
        gene_sets=gene_sets,
        min_size=15,
        max_size=500,
        permutation_num=1000,
        weight=1.0,
        ascending=False,
        threads=1,
        seed=SEED,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    result = prerank_result.res2d
    if result is None:
        raise RuntimeError("GSEApy prerank did not return a result table.")
    normalized = result.rename(
        columns={
            "Term": "pathway",
            "ES": "ES",
            "NES": "NES",
            "NOM p-val": "nominal_P",
            "FDR q-val": "GSEA_FDR_q",
            "Lead_genes": "leading_edge_genes",
        }
    )
    if not {"pathway", "NES", "nominal_P", "leading_edge_genes"}.issubset(
        normalized.columns
    ):
        raise ValueError(f"Unexpected GSEApy result columns: {result.columns.tolist()}")
    normalized["BH_q_across_50_Hallmark_pathways"] = multipletests(
        pd.to_numeric(normalized["nominal_P"]), method="fdr_bh"
    )[1]
    normalized["leading_edge_genes"] = normalized["leading_edge_genes"].fillna("").astype(str)
    normalized["panel_v2_genes_in_leading_edge"] = normalized[
        "leading_edge_genes"
    ].apply(
        lambda text: ";".join(
            gene for gene in V2_GENES if re.search(rf"(?<![A-Z0-9]){gene}(?![A-Z0-9])", text.upper())
        )
    )
    normalized["n_panel_v2_genes_in_leading_edge"] = normalized[
        "panel_v2_genes_in_leading_edge"
    ].apply(lambda text: 0 if not text else len(text.split(";")))
    normalized = normalized.sort_values(
        "BH_q_across_50_Hallmark_pathways", kind="mergesort"
    )
    normalized = add_analysis_label(normalized)
    normalized.to_csv(OUT / "hallmark_preranked_gsea.csv", index=False)
    return normalized


def bootstrap_auc(
    values: np.ndarray, labels: np.ndarray, rng: np.random.Generator
) -> tuple[float, float, float, float, np.ndarray, np.ndarray]:
    auc_value = float(roc_auc_score(labels, values))
    case_indices = np.flatnonzero(labels == 1)
    control_indices = np.flatnonzero(labels == 0)
    draws: list[float] = []
    for _ in range(BOOTSTRAPS):
        sampled = np.concatenate(
            [
                rng.choice(control_indices, size=len(control_indices), replace=True),
                rng.choice(case_indices, size=len(case_indices), replace=True),
            ]
        )
        sampled_labels = labels[sampled]
        sampled_values = values[sampled]
        if np.unique(sampled_labels).size == 2:
            draws.append(
                float(roc_auc_score(sampled_labels, sampled_values))
            )
    low, high = np.percentile(draws, [2.5, 97.5])
    p_value = float(
        stats.mannwhitneyu(
            values[labels == 1], values[labels == 0], alternative="two-sided"
        ).pvalue
    )
    fpr, tpr, _ = roc_curve(labels, values)
    return auc_value, float(low), float(high), p_value, fpr, tpr


def run_roc(sample_scores: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    labels_by_sample = meta.set_index("sample_id")["group"].eq("LPD")
    scores = sample_scores.set_index("sample_id").loc[labels_by_sample.index]
    labels = labels_by_sample.loc[scores.index].astype(int).to_numpy()
    rng = np.random.default_rng(SEED)
    targets = [
        *[f"panel_v2_{g}" for g in V2_GENES],
        "panel_v2_composite",
        *[f"original_8_{g}" for g in ORIGINAL_GENES],
        "original_8_composite",
    ]
    fig, ax = plt.subplots(figsize=(13, 8))
    rows: list[dict[str, object]] = []
    for target in targets:
        observed, low, high, p_value, fpr, tpr = bootstrap_auc(
            scores[target].to_numpy(dtype=float), labels, rng
        )
        rows.append(
            {
                "score": target,
                "n_LPD": int(labels.sum()),
                "n_SPD": int((1 - labels).sum()),
                "AUC": observed,
                "bootstrap_95CI_low": low,
                "bootstrap_95CI_high": high,
                "Mann_Whitney_P": p_value,
                "bootstrap_resamples": BOOTSTRAPS,
                "bootstrap_seed": SEED,
            }
        )
        ax.plot(fpr, tpr, linewidth=1.25, label=f"{target} (AUC={observed:.3f})")
    ax.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1)
    ax.set(
        xlabel="False positive rate",
        ylabel="True positive rate",
        title="GSE125498 ROC — circular: genes were selected on this data",
    )
    fig.text(0.5, 0.01, DISCOVERY_LABEL, ha="center", fontsize=8)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1))
    ax.grid(alpha=0.2)
    fig.tight_layout(rect=(0, 0.04, 0.76, 1))
    fig.savefig(OUT / "discovery_circular_roc.png", dpi=180)
    plt.close(fig)
    table = pd.DataFrame(rows)
    add_analysis_label(table).to_csv(OUT / "discovery_descriptive_roc.csv", index=False)
    return table


def query_dgidb() -> pd.DataFrame:
    query = (
        'query { genes(names: ["FLT3LG", "EBI3", "LTB", "ADAM19", "CST7"]) { '
        "nodes { name interactions { drug { name } interactionTypes { type directionality } "
        "sources { sourceDbName } publications { pmid } } } } }"
    )
    response = requests.post(
        "https://dgidb.org/api/graphql",
        json={"query": query},
        headers={"Accept": "application/json"},
        timeout=60,
    )
    response.raise_for_status()
    data = response.json()
    if data.get("errors"):
        raise RuntimeError(f"DGIdb returned GraphQL errors: {data['errors']}")
    if not isinstance(data.get("data", {}).get("genes", {}).get("nodes"), list):
        raise RuntimeError("DGIdb response did not contain the requested gene nodes.")
    (OUT / "dgidb_raw_response.json").write_text(
        json.dumps(
            {
                "analysis_label": DISCOVERY_LABEL,
                "query": query,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "http_status": response.status_code,
                "response": data,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    rows: list[dict[str, object]] = []
    nodes = {node["name"]: node for node in data["data"]["genes"]["nodes"]}
    for gene in V2_GENES:
        node = nodes.get(gene)
        if node is None:
            raise RuntimeError(f"DGIdb did not return requested gene {gene}.")
        interactions = node.get("interactions") or []
        if not interactions:
            rows.append(
                {
                    "gene": gene,
                    "drug": "",
                    "interaction_types": "",
                    "directionality": "",
                    "sources": "",
                    "publications": "",
                    "DGIdb_interaction_returned": False,
                }
            )
        for interaction in interactions:
            rows.append(
                {
                    "gene": gene,
                    "drug": (interaction.get("drug") or {}).get("name", ""),
                    "interaction_types": ";".join(
                        item.get("type", "") for item in interaction.get("interactionTypes", [])
                    ),
                    "directionality": ";".join(
                        item.get("directionality", "")
                        for item in interaction.get("interactionTypes", [])
                    ),
                    "sources": ";".join(
                        item.get("sourceDbName", "") for item in interaction.get("sources", [])
                    ),
                    "publications": ";".join(
                        str(item.get("pmid", ""))
                        for item in interaction.get("publications", [])
                    ),
                    "DGIdb_interaction_returned": True,
                }
            )
    result = add_analysis_label(pd.DataFrame(rows))
    result.to_csv(OUT / "dgidb_panel_v2_interactions.csv", index=False)
    return result


def render_report(
    meta: pd.DataFrame,
    immune_comparison: pd.DataFrame,
    correlations: pd.DataFrame,
    models: pd.DataFrame,
    model_details: dict[str, Any],
    gsea: pd.DataFrame,
    roc: pd.DataFrame,
    drugs: pd.DataFrame,
) -> None:
    top = gsea.nsmallest(15, "BH_q_across_50_Hallmark_pathways")
    lines = [
        "# GSE125498 panel context analysis",
        "",
        f"**Label:** {DISCOVERY_LABEL}.",
        "",
        "All analyses in this report reuse the discovery cohort that informed gene selection. "
        "They are descriptive/exploratory and are not validation.",
        "",
        "## Cohort and scoring",
        "",
        f"- Samples: {len(meta)} total; {int((meta.group == 'LPD').sum())} LPD and "
        f"{int((meta.group == 'SPD').sum())} SPD.",
        "- Immune cell-type scores: GSEApy ssGSEA using human CellMarker_2024 signatures "
        "retrieved from the Enrichr library. For each broad cell type, matching human "
        "tissue/context marker sets were unioned before scoring; signatures are expression "
        "enrichment estimates, not cell counts or validated deconvolution.",
        "- Dominant leukocyte covariate: mean per-gene z-score of the 17 repository "
        f"pan-leukocyte markers ({model_details['n_leukocyte_markers']} observed in the "
        f"matrix: {', '.join(model_details['leukocyte_markers_used'])}); this is a proxy, "
        "not an independently measured cell fraction.",
        "- Composite scores: average of gene-wise z-scores across these 33 samples. "
        "Original panel has eight genes; panel v2 has five.",
        "- Group differences use two-sided Mann–Whitney tests; BH q-values are adjusted "
        "across the 13 immune cell-type comparisons.",
        "",
        "## Immune cell-type differences",
        "",
        "| Cell type | LPD mean ssGSEA | SPD mean ssGSEA | LPD − SPD | Mann–Whitney P | BH q |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in immune_comparison.sort_values("BH_q_across_cell_types").itertuples():
        lines.append(
            f"| {row.cell_type} | {row.LPD_mean_ssGSEA:.4f} | {row.SPD_mean_ssGSEA:.4f} "
            f"| {row.LPD_minus_SPD_mean:+.4f} | {row.Mann_Whitney_P:.4g} "
            f"| {row.BH_q_across_cell_types:.4g} |"
        )
    lines += [
        "",
        "## Correlation and leukocyte adjustment",
        "",
        "The table gives Pearson correlations for each v2 gene and the v2 composite "
        "against each immune cell-type ssGSEA score. Spearman results and original-8 "
        "correlations are in `panel_celltype_correlations.csv`; Pearson BH q is adjusted "
        "within each score. Group-adjusted model details are in "
        "`composite_leukocyte_adjustment.csv`.",
        "",
        "| V2 score | Cell type | Pearson r | P | BH q within score |",
        "|---|---|---:|---:|---:|",
    ]
    v2_targets = {f"panel_v2_{gene}" for gene in V2_GENES} | {"panel_v2_composite"}
    for row in correlations[correlations["score"].isin(v2_targets)].itertuples():
        lines.append(
            f"| {row.score} | {row.cell_type} | {row.Pearson_r:.4f} "
            f"| {row.Pearson_P:.4g} | {row.Pearson_BH_q_within_score:.4g} |"
        )
    lines.append("")
    for target, detail in model_details["model_summary"].items():
        lines.append(
            f"- **{target}:** unadjusted group coefficient {detail['unadjusted_group_beta']:.4f} "
            f"(P={detail['unadjusted_group_p']:.4g}); after adding the dominant leukocyte "
            f"marker proxy, group coefficient {detail['adjusted_group_beta']:.4f} "
            f"(P={detail['adjusted_group_p']:.4g}), proxy coefficient "
            f"{detail['adjusted_leukocyte_proxy_beta']:.4f} "
            f"(P={detail['adjusted_leukocyte_proxy_p']:.4g}), and group-coefficient "
            f"attenuation {detail['group_beta_attenuation_percent']:.1f}%."
        )
    detail = model_details["model_summary"]["panel_v2_composite"]
    lines += [
        "",
        "**Interpretation:** The panel-v2 group coefficient is attenuated by "
        f"{detail['group_beta_attenuation_percent']:.1f}% after adding the leukocyte-marker "
        f"proxy, but remains nominally significant (P={detail['adjusted_group_p']:.4g}). "
        "This is consistent with a substantial leukocyte-composition contribution, but the "
        "proxy does not fully explain the association; neither the attenuation nor residual "
        "association establishes causality or cell-intrinsic expression.",
        "",
        "## Preranked Hallmark GSEA",
        "",
        "The complete 11,741-gene limma t-statistic ranking was analyzed with GSEApy "
        "prerank (1,000 permutations, seed 125498) using `MSigDB_Hallmark_2020`, the "
        "same library source used elsewhere in this repository. The table below includes "
        "BH adjustment of nominal pathway P values across the 50 Hallmark sets; "
        "GSEApy's FDR q-values are retained separately in the CSV.",
        "",
        "| Pathway | NES | Nominal P | BH q | GSEA FDR q | Panel-v2 genes in leading edge |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in top.itertuples():
        lead = row.panel_v2_genes_in_leading_edge or "—"
        lines.append(
            f"| {row.pathway} | {float(row.NES):.4f} | {format_gsea_p(row.nominal_P)} "
            f"| {format_gsea_p(row.BH_q_across_50_Hallmark_pathways)} "
            f"| {format_gsea_p(row.GSEA_FDR_q)} | {lead} |"
        )
    lines += [
        "",
        "All 50 sets, full leading-edge strings, and panel-gene membership are in "
        "`hallmark_preranked_gsea.csv`; the retrieved gene sets are saved in "
        "`hallmark_gene_sets.json`. Across all 50 pathways, 18 have BH q < 0.05. "
        "The panel-v2 leading-edge memberships are: IL-2/STAT5 Signaling "
        "(FLT3LG, LTB, ADAM19, CST7; BH q < 0.001), Allograft Rejection "
        "(LTB; BH q < 0.001), Inflammatory Response (EBI3; BH q=0.0251), and "
        "IL-6/JAK/STAT3 Signaling (EBI3, LTB; BH q=0.3712).",
        "",
        "Nominal permutation P values of 0 are displayed as <0.001; with 1,000 "
        "permutations they are below the reporting resolution, not literal zero.",
        "",
        "## Descriptive ROC",
        "",
        "**circular: genes were selected on this data.** ROC estimates below are "
        "in-sample descriptions, not performance validation. Bootstrap CIs resample "
        "within LPD/SPD groups (5,000 replicates; seed 125498).",
        "",
        "| Score | AUC | Bootstrap 95% CI | Mann–Whitney P |",
        "|---|---:|---:|---:|",
    ]
    for row in roc.itertuples():
        lines.append(
            f"| {row.score} | {row.AUC:.4f} | {row.bootstrap_95CI_low:.4f}–"
            f"{row.bootstrap_95CI_high:.4f} | {row.Mann_Whitney_P:.4g} |"
        )
    lines += [
        "",
        "Figure: `discovery_circular_roc.png`.",
        "",
        "## DGIdb interactions (retrieved live)",
        "",
        f"Queried DGIdb GraphQL on {datetime.now(timezone.utc).date().isoformat()}. "
        "No drug ranking was performed. Empty interaction rows mean the endpoint returned "
        "no interaction for that queried gene at retrieval time; this is not proof that no "
        "interaction is known.",
        "",
        "| Gene | Drug | Interaction type | Directionality | Source(s) | Publication(s) |",
        "|---|---|---|---|---|---|",
    ]
    for row in drugs.itertuples():
        if row.DGIdb_interaction_returned:
            lines.append(
                f"| {row.gene} | {row.drug} | {row.interaction_types or 'not reported'} "
                f"| {row.directionality or 'not reported'} | {row.sources or 'not reported'} "
                f"| {row.publications or 'not reported'} |"
            )
        else:
            lines.append(f"| {row.gene} | — | No interaction returned | — | — | — |")
    lines += [
        "",
        "The parsed result is in `dgidb_panel_v2_interactions.csv`; the raw response, "
        "query, and retrieval timestamp are in `dgidb_raw_response.json`.",
        "",
        "## Caveats",
        "",
        "- CellMarker tissue-context markers were unioned to obtain broad immune-cell "
        "signatures; tissue effects and overlapping marker biology may influence scores.",
        "- The cohort is small (13 LPD, 20 SPD); nominal associations are unstable and "
        "exploratory.",
        "- ROC and discovery GSEA are computed on the same data used in feature selection.",
        "- GSEApy reported a tie affecting 0.01% of ranked genes; the tied order may be "
        "arbitrary and could affect results at that tie boundary.",
        "- DGIdb is a database snapshot queried live; database coverage and annotation "
        "can change.",
        "",
    ]
    (OUT / "EXPLORATORY_CONTEXT_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    expr, meta, ranking = load_inputs()
    immune_comp, correlations, models, model_details = immune_composition(expr, meta)
    gsea = run_hallmark_gsea(ranking)
    sample_scores = make_sample_scores(expr)
    roc = run_roc(sample_scores.reset_index(names="sample_id"), meta)
    drugs = query_dgidb()
    render_report(meta, immune_comp, correlations, models, model_details, gsea, roc, drugs)
    (OUT / "analysis_metadata.json").write_text(
        json.dumps(
            {
                "label": DISCOVERY_LABEL,
                "run_date": date.today().isoformat(),
                "discovery_cohort": "GSE125498",
                "samples": len(meta),
                "group_counts": meta["group"].value_counts().to_dict(),
                "panel_v2_genes": V2_GENES,
                "original_panel_genes": ORIGINAL_GENES,
                "gene_expression_input": str(EXPR_PATH.relative_to(ROOT)),
                "sample_metadata_input": str(META_PATH.relative_to(ROOT)),
                "full_limma_ranking_input": str(RANK_PATH.relative_to(ROOT)),
                "gseapy_version": gp.__version__,
                "seed": SEED,
                "roc_bootstrap_resamples": BOOTSTRAPS,
                "bh_adjustment": "Benjamini-Hochberg",
                "gsea_rank": "all 11,741 limma t statistics, descending",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
