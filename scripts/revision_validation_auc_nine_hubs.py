"""Create exploratory external-cohort ROC/AUC figures for the nine-gene hub panel."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_curve, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
PANEL_PATH = ROOT / "results" / "revision" / "ml_39_matrisome_wgcna" / "ml_39_gene_results.csv"
OUTPUT_DIR = ROOT / "results" / "revision" / "ml_39_matrisome_wgcna" / "validation_auc"
BOOTSTRAP_REPLICATES = 5_000
SEED = 42

PANEL_GENES = [
    "TNFSF15",
    "FLT3LG",
    "EBI3",
    "LTB",
    "ADAM19",
    "NAT6",
    "ECM1",
    "CST7",
    "CCBE1",
]
COMMON_GENES = [
    "FLT3LG",
    "EBI3",
    "LTB",
    "ADAM19",
    "NAT6",
    "ECM1",
    "CST7",
    "CCBE1",
]
GENE_IDS = {
    "TNFSF15": "ENSG00000181634",
    "FLT3LG": "ENSG00000090554",
    "EBI3": "ENSG00000105246",
    "LTB": "ENSG00000227507",
    "ADAM19": "ENSG00000135074",
    "NAT6": "ENSG00000243477",
    "ECM1": "ENSG00000143369",
    "CST7": "ENSG00000077984",
    "CCBE1": "ENSG00000183287",
}
COHORTS = {
    "GSE248762": {
        "contrast": "LV (n=10) vs SV (n=6)",
        "case": "LV",
        "control": "SV",
        "n_cases_expected": 10,
        "n_controls_expected": 6,
    },
    "GSE130888": {
        "contrast": "Long-term PD (n=4) vs short-term PD (n=6)",
        "case": "LV",
        "control": "SV",
        "n_cases_expected": 4,
        "n_controls_expected": 6,
    },
}


def bootstrap_auc(
    y: np.ndarray, scores: np.ndarray, rng: np.random.Generator
) -> tuple[float, float]:
    case_idx = np.flatnonzero(y == 1)
    control_idx = np.flatnonzero(y == 0)
    sampled_cases = rng.choice(
        case_idx, size=(BOOTSTRAP_REPLICATES, len(case_idx)), replace=True
    )
    sampled_controls = rng.choice(
        control_idx, size=(BOOTSTRAP_REPLICATES, len(control_idx)), replace=True
    )
    case_scores = scores[sampled_cases]
    control_scores = scores[sampled_controls]
    pairwise = (
        (case_scores[:, :, None] > control_scores[:, None, :]).astype(float)
        + 0.5 * (case_scores[:, :, None] == control_scores[:, None, :])
    )
    aucs = pairwise.mean(axis=(1, 2))
    return tuple(float(v) for v in np.percentile(aucs, [2.5, 97.5]))


def load_cohort(cohort: str, settings: dict) -> dict:
    cohort_dir = ROOT / "results" / "revision" / "effluent_validation" / cohort / "analysis_20261009"
    metadata = pd.read_csv(cohort_dir / "all_cell_pseudobulk_metadata.csv")
    metadata = metadata.loc[metadata["primary_group"].isin([settings["case"], settings["control"]])].copy()
    metadata = metadata.drop_duplicates("sample_id").set_index("sample_id")
    labels = (metadata["primary_group"] == settings["case"]).astype(int)
    y = labels.to_numpy(dtype=int)
    n_cases = int(y.sum())
    n_controls = int(len(y) - n_cases)
    if (n_cases, n_controls) != (
        settings["n_cases_expected"],
        settings["n_controls_expected"],
    ):
        raise ValueError(f"{cohort} group sizes were {n_cases}/{n_controls}, contrary to expected values.")

    expression = pd.read_csv(cohort_dir / "edger_logCPM.csv", index_col=0)
    expression.index = expression.index.astype(str)
    missing_samples = set(metadata.index) - set(expression.columns)
    if missing_samples:
        raise ValueError(f"{cohort} expression matrix is missing samples: {sorted(missing_samples)}")

    detection = pd.read_csv(cohort_dir / "pseudobulk_gene_detection.csv", dtype={"gene_id": str})
    detection = detection.set_index("gene_id")
    detected_genes = set()
    for gene in PANEL_GENES:
        gene_id = GENE_IDS[gene]
        if gene_id not in detection.index:
            continue
        if bool(detection.loc[gene_id, "detected"]):
            detected_genes.add(gene)
    if cohort == "GSE248762" and len(detected_genes) != 9:
        raise ValueError(f"GSE248762 should have all 9 genes detected; found {sorted(detected_genes)}.")
    if cohort == "GSE130888" and detected_genes != set(COMMON_GENES):
        raise ValueError(f"GSE130888 detected panel set differs from expected: {sorted(detected_genes)}.")

    values: dict[str, np.ndarray] = {}
    for gene in detected_genes:
        gene_id = GENE_IDS[gene]
        if gene_id not in expression.index:
            raise ValueError(f"{cohort} logCPM matrix is missing {gene} ({gene_id}).")
        values[gene] = expression.loc[gene_id, metadata.index].to_numpy(dtype=float)
    return {
        "cohort": cohort,
        "dir": cohort_dir,
        "metadata": metadata,
        "y": y,
        "values": values,
        "detected_genes": detected_genes,
        "n_cases": n_cases,
        "n_controls": n_controls,
        "settings": settings,
    }


def calculate_results(cohorts: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rng = np.random.default_rng(SEED)
    gene_rows = []
    composite_rows = []
    score_data = {}
    for cohort_name, data in cohorts.items():
        y = data["y"]
        cohort_gene_scores = {}
        for gene in PANEL_GENES:
            if gene not in data["detected_genes"]:
                gene_rows.append(
                    {
                        "Cohort": cohort_name,
                        "Gene": gene,
                        "Detected": False,
                        "N_case": data["n_cases"],
                        "N_control": data["n_controls"],
                        "AUC_higher_expression_predicts_LV": np.nan,
                        "Bootstrap_CI95_low": np.nan,
                        "Bootstrap_CI95_high": np.nan,
                        "Direction_concordant_with_discovery": "not assessed",
                    }
                )
                continue
            scores = data["values"][gene]
            value_auc = float(roc_auc_score(y, scores))
            ci_low, ci_high = bootstrap_auc(y, scores, rng)
            edge = pd.read_csv(
                ROOT / "results" / "revision" / "effluent_validation" / cohort_name
                / "analysis_20261009" / "edger_all_gene_results.csv"
            )
            edge = edge.loc[
                (edge["contrast"] == "LV_vs_SV")
                & (edge["gene_id"] == GENE_IDS[gene])
            ]
            if len(edge) != 1:
                raise ValueError(f"Expected one LV_vs_SV edgeR row for {cohort_name}:{gene}, found {len(edge)}.")
            concordant = bool(float(edge.iloc[0]["logFC"]) > 0)
            cohort_gene_scores[gene] = scores
            gene_rows.append(
                {
                    "Cohort": cohort_name,
                    "Gene": gene,
                    "Detected": True,
                    "N_case": data["n_cases"],
                    "N_control": data["n_controls"],
                    "AUC_higher_expression_predicts_LV": value_auc,
                    "Bootstrap_CI95_low": ci_low,
                    "Bootstrap_CI95_high": ci_high,
                    "Direction_concordant_with_discovery": concordant,
                    "edgeR_log2FC": float(edge.iloc[0]["logFC"]),
                    "edgeR_P": float(edge.iloc[0]["PValue"]),
                    "edgeR_FDR": float(edge.iloc[0]["FDR_all_genes"]),
                }
            )

        common_values = np.vstack([data["values"][gene] for gene in COMMON_GENES]).T
        common_z = (common_values - common_values.mean(axis=0)) / common_values.std(axis=0, ddof=0)
        composite_score = common_z.mean(axis=1)
        composite_auc = float(roc_auc_score(y, composite_score))
        ci_low, ci_high = bootstrap_auc(y, composite_score, rng)
        composite_rows.append(
            {
                "Cohort": cohort_name,
                "Score": "mean_of_8_common_gene_z_scores",
                "Genes": ";".join(COMMON_GENES),
                "N_genes": len(COMMON_GENES),
                "N_case": data["n_cases"],
                "N_control": data["n_controls"],
                "AUC_higher_score_predicts_LV": composite_auc,
                "Bootstrap_CI95_low": ci_low,
                "Bootstrap_CI95_high": ci_high,
            }
        )
        score_data[cohort_name] = {
            "y": y,
            "score": composite_score,
            "auc": composite_auc,
            "genes": list(COMMON_GENES),
        }

    gene_df = pd.DataFrame(gene_rows)
    composite_df = pd.DataFrame(composite_rows)
    summary = {
        "analysis_type": "Exploratory external-cohort, expression-score ROC analysis",
        "panel_genes": PANEL_GENES,
        "cohort_contrasts": {
            cohort: {
                "contrast": data["settings"]["contrast"],
                "n_case": data["n_cases"],
                "n_control": data["n_controls"],
                "detected_panel_genes": sorted(data["detected_genes"]),
                "not_detected_panel_genes": sorted(set(PANEL_GENES) - data["detected_genes"]),
            }
            for cohort, data in cohorts.items()
        },
        "per_gene_auc_definition": (
            "ROC AUC of the gene's validation logCPM, with higher expression scored as predicting LV; "
            "directions are not flipped to maximize AUC."
        ),
        "composite_definition": (
            "Equal-weight mean of within-cohort z-scores for the same eight genes detected in both cohorts; "
            "no classifier was fitted. Standardization did not use phenotype labels."
        ),
        "ci_method": f"Stratified percentile bootstrap, {BOOTSTRAP_REPLICATES} resamples per score; seed {SEED}.",
        "limitations": [
            "Small validation sample sizes make AUC estimates imprecise.",
            "These are exploratory score evaluations, not fitted logistic-regression performance estimates.",
            "Panel discovery used GSE125498; these plots use additional effluent cohorts.",
            "Independence of GSE248762 and GSE130888 patient samples is not formally established.",
            "The shared eight-gene composite omits TNFSF15, which was not detected in GSE130888.",
        ],
    }
    return gene_df, composite_df, {"summary": summary, "score_data": score_data}


def make_figure(gene_df: pd.DataFrame, composite_df: pd.DataFrame, details: dict) -> None:
    colors = {"GSE248762": "#2676A5", "GSE130888": "#D47727"}
    fig, (ax_auc, ax_roc) = plt.subplots(
        1, 2, figsize=(15, 8), gridspec_kw={"width_ratios": [1.5, 1.0]}
    )
    fig.suptitle(
        "Exploratory ROC analysis of the nine consensus genes in effluent cohorts",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.935,
        "GSE248762 (LV 10 vs SV 6) and GSE130888 (LV 4 vs SV 6); "
        "gene AUCs use higher expression as the LV score",
        ha="center",
        fontsize=9.5,
        color="#475569",
    )

    gene_order = PANEL_GENES
    positions = np.arange(len(gene_order))
    offsets = {"GSE248762": -0.12, "GSE130888": 0.12}
    for cohort_name, color in colors.items():
        subset = gene_df.loc[gene_df["Cohort"] == cohort_name].set_index("Gene")
        for i, gene in enumerate(gene_order):
            if gene not in subset.index or not bool(subset.loc[gene, "Detected"]):
                continue
            row = subset.loc[gene]
            y_pos = i + offsets[cohort_name]
            ax_auc.errorbar(
                row["AUC_higher_expression_predicts_LV"],
                y_pos,
                xerr=[
                    [row["AUC_higher_expression_predicts_LV"] - row["Bootstrap_CI95_low"]],
                    [row["Bootstrap_CI95_high"] - row["AUC_higher_expression_predicts_LV"]],
                ],
                fmt="o",
                color=color,
                ecolor=color,
                capsize=2.5,
                markersize=6,
                label=cohort_name if i == 0 else None,
            )
    ax_auc.axvline(0.5, color="#64748B", linestyle="--", linewidth=1)
    ax_auc.set_yticks(positions)
    ax_auc.set_yticklabels(gene_order, fontstyle="italic")
    ax_auc.invert_yaxis()
    ax_auc.set_xlim(0, 1)
    ax_auc.set_xlabel("AUC (higher expression predicts long-vintage / long-term PD)")
    ax_auc.set_title("A. Per-gene AUC (95% stratified bootstrap CI)", loc="left", fontweight="bold")
    ax_auc.grid(axis="x", color="#E2E8F0", linewidth=0.8)
    ax_auc.legend(frameon=False, loc="lower right")

    ax_roc.plot([0, 1], [0, 1], linestyle="--", color="#94A3B8", linewidth=1, label="Chance (AUC = 0.50)")
    for cohort_name, color in colors.items():
        score = details["score_data"][cohort_name]
        fpr, tpr, _ = roc_curve(score["y"], score["score"])
        composite_row = composite_df.loc[composite_df["Cohort"] == cohort_name].iloc[0]
        ax_roc.plot(
            fpr,
            tpr,
            color=color,
            linewidth=2.2,
            label=(
                f"{cohort_name}, common 8-gene score "
                f"(AUC {composite_row['AUC_higher_score_predicts_LV']:.2f})"
            ),
        )
    ax_roc.set_xlim(0, 1)
    ax_roc.set_ylim(0, 1)
    ax_roc.set_aspect("equal", adjustable="box")
    ax_roc.set_xlabel("False-positive rate")
    ax_roc.set_ylabel("True-positive rate")
    ax_roc.set_title("B. Fixed common 8-gene score ROC", loc="left", fontweight="bold")
    ax_roc.grid(color="#E2E8F0", linewidth=0.8)
    ax_roc.legend(frameon=False, loc="lower right", fontsize=8.5)

    fig.text(
        0.5,
        0.025,
        "Common score: equal-weight mean of within-cohort z-scores for FLT3LG, EBI3, LTB, ADAM19, "
        "NAT6/NAA80, ECM1, CST7 and CCBE1. TNFSF15 was unavailable in GSE130888. "
        "Exploratory; small cohorts and independence not formally confirmed.",
        ha="center",
        fontsize=8,
        color="#475569",
    )
    fig.tight_layout(rect=(0.03, 0.07, 0.98, 0.91))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / "nine_hub_external_validation_auc.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "nine_hub_external_validation_auc.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    panel = pd.read_csv(PANEL_PATH)
    selected = set(
        panel.loc[panel["Consensus_At_Least_2_of_3"] == 1, "Gene"].astype(str)
    )
    if selected != set(PANEL_GENES):
        raise ValueError(f"Saved nine-gene panel differs from the expected panel: {sorted(selected)}")

    cohorts = {
        cohort: load_cohort(cohort, settings)
        for cohort, settings in COHORTS.items()
    }
    gene_df, composite_df, details = calculate_results(cohorts)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gene_df.to_csv(OUTPUT_DIR / "per_gene_validation_auc.csv", index=False)
    composite_df.to_csv(OUTPUT_DIR / "common_8_gene_composite_auc.csv", index=False)
    make_figure(gene_df, composite_df, details)
    (OUTPUT_DIR / "validation_auc_summary.json").write_text(
        json.dumps(details["summary"], indent=2), encoding="utf-8"
    )
    print(gene_df.to_string(index=False))
    print("\nCommon 8-gene composite:")
    print(composite_df.to_string(index=False))
    print(f"\nSaved figures and tables under {OUTPUT_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
