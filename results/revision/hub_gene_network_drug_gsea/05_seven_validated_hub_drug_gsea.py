"""
Potential drug prediction (DSigDB/Enrichr) and single-gene KEGG GSEA for the
seven validated consensus hub genes (Version A):

    TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, CST7

This script re-runs the drug-signature enrichment and the preranked single-gene
KEGG GSEA on the GSE125498 discovery cohort (33 samples; 13 LPD / 20 SPD) for
the seven validated hubs, then renders two publication-quality figures:

  * validated7_dsigdb_drug_candidates.png / .pdf
  * validated7_kegg_gsea_panel.png / .pdf

All outputs are written into this directory; the earlier seven-gene exploratory
panel (which included SERPINA10 and CXCL14) is left untouched.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import gseapy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import requests
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr

matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#333333",
        "figure.facecolor": "white",
    }
)

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
EXPR_PATH = ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv"
META_PATH = ROOT / "results" / "tables" / "GSE125498_sample_metadata.csv"

# The seven validated consensus hub genes (Version A).
GENES = ["TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "ECM1", "CST7"]

# DSigDB settings
ENRICHR = "https://maayanlab.cloud/Enrichr"
DSIGDB_LIBRARY = "DSigDB"
FDR_CUTOFF = 0.05

# KEGG GSEA settings
KEGG = "https://rest.kegg.jp"
SEED = 20261013
PERMUTATIONS = 1_000
MIN_SIZE = 15
MAX_SIZE = 500

# Figure palette
COLOR_ACTIVATED = "#C0392B"   # red - activated (NES > 0)
COLOR_SUPPRESSED = "#2E5FA3"  # blue - suppressed (NES < 0)
COLOR_SIG = "#C44E52"
COLOR_NS = "#9AA0A6"


def run_dsigdb() -> pd.DataFrame:
    """Query Enrichr DSigDB for the seven validated hub genes."""
    add_response = requests.post(
        f"{ENRICHR}/addList",
        files={
            "list": (None, "\n".join(GENES)),
            "description": (None, "GSE125498 seven validated consensus hub genes (Version A)"),
        },
        timeout=60,
    )
    add_response.raise_for_status()
    user_list_id = add_response.json()["userListId"]

    enrich_response = requests.get(
        f"{ENRICHR}/enrich",
        params={"userListId": user_list_id, "backgroundType": DSIGDB_LIBRARY},
        timeout=180,
    )
    enrich_response.raise_for_status()
    enrich_json = enrich_response.json()
    if DSIGDB_LIBRARY not in enrich_json:
        raise RuntimeError(f"Enrichr response did not include {DSIGDB_LIBRARY}.")

    rows = []
    for values in enrich_json[DSIGDB_LIBRARY]:
        rank, term, p_value, z_score, combined_score, genes, adjusted_p, *_ = values
        overlap_genes = list(genes)
        rows.append(
            {
                "rank": int(rank),
                "dsigdb_term": term,
                "overlapping_genes": ";".join(overlap_genes),
                "overlap_size": len(overlap_genes),
                "p_value": float(p_value),
                "z_score": float(z_score),
                "combined_score": float(combined_score),
                "adjusted_p_value": float(adjusted_p),
                "significant_fdr_lt_0_05": float(adjusted_p) < FDR_CUTOFF,
            }
        )
    results = pd.DataFrame(rows).sort_values(
        ["adjusted_p_value", "p_value", "rank"]
    )
    results.to_csv(OUT / "validated7_dsigdb_enrichment_full.csv", index=False)
    significant = results[results["adjusted_p_value"] < FDR_CUTOFF].copy()
    significant.to_csv(
        OUT / "validated7_dsigdb_candidate_drugs_fdr_lt_0_05.csv", index=False
    )

    metadata = {
        "query_date": date.today().isoformat(),
        "database": "DSigDB (via Enrichr)",
        "library_version": "DSigDB v1.0 as listed by Enrichr",
        "query_genes": GENES,
        "enrichr_user_list_id": user_list_id,
        "fdr_cutoff": FDR_CUTOFF,
        "number_of_terms_returned": len(results),
        "number_significant_fdr_lt_0_05": int(significant.shape[0]),
    }
    (OUT / "validated7_dsigdb_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(
        f"[DSigDB] {len(results)} signatures returned; "
        f"{significant.shape[0]} pass FDR < {FDR_CUTOFF}",
        flush=True,
    )
    return results


def fetch_kegg_sets() -> tuple[dict[str, set[str]], int, str, str]:
    """Download human KEGG pathway definitions and gene-to-pathway links."""
    pathway_response = requests.get(f"{KEGG}/list/pathway/hsa", timeout=90)
    pathway_response.raise_for_status()
    pathway_names = {
        line.split("\t", 1)[0].removeprefix("path:"): line.split("\t", 1)[1]
        for line in pathway_response.text.splitlines()
        if "\t" in line
    }

    gene_response = requests.get(f"{KEGG}/list/hsa", timeout=120)
    gene_response.raise_for_status()
    gene_names: dict[str, str] = {}
    for line in gene_response.text.splitlines():
        if "\t" not in line:
            continue
        fields = line.split("\t")
        if len(fields) >= 4:
            gene_id, _, _, description = fields[:4]
        elif len(fields) == 3:
            gene_id, _, description = fields
        else:
            continue
        symbol = description.split(",", 1)[0].split(";", 1)[0].strip()
        if symbol:
            gene_names[gene_id] = symbol

    links_response = requests.get(f"{KEGG}/link/pathway/hsa", timeout=120)
    links_response.raise_for_status()
    gene_sets: dict[str, set[str]] = {name: set() for name in pathway_names.values()}
    for line in links_response.text.splitlines():
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        gene_id = parts[0].removeprefix("hsa:")
        pathway_id = parts[1].removeprefix("path:")
        symbol = gene_names.get(f"hsa:{gene_id}")
        pathway_name = pathway_names.get(pathway_id)
        if symbol and pathway_name:
            gene_sets[pathway_name].add(symbol)

    info_response = requests.get(f"{KEGG}/info/kegg", timeout=60)
    info_response.raise_for_status()
    pathway_release = "not reported"
    gene_release = "not reported"
    for line in info_response.text.splitlines():
        fields = line.split()
        if fields and len(fields) >= 3:
            if fields[0] == "pathway":
                pathway_release = fields[-1]
            elif fields[0] == "genes":
                gene_release = fields[-1]

    gene_sets = {name: members for name, members in gene_sets.items() if members}
    if not gene_sets:
        raise RuntimeError("KEGG REST returned no human pathway gene sets.")
    return gene_sets, len(pathway_names), pathway_release, gene_release


def run_kegg_gsea() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run preranked single-gene KEGG GSEA on the discovery cohort."""
    expression_df = pd.read_csv(EXPR_PATH)
    metadata = pd.read_csv(META_PATH)
    if "Gene" not in expression_df or not {"sample_id", "stage_binary"}.issubset(
        metadata.columns
    ):
        raise ValueError("Discovery expression/metadata columns do not match conventions.")
    if expression_df["Gene"].duplicated().any():
        raise ValueError("Expected unique gene symbols in the expression matrix.")

    samples = metadata["sample_id"].astype(str).tolist()
    if len(samples) != 33 or len(set(samples)) != 33:
        raise ValueError(f"Expected 33 unique discovery samples, found {len(samples)}.")
    matrix = expression_df.set_index("Gene")[samples].astype(float)
    if not np.isfinite(matrix.to_numpy()).all():
        raise ValueError("Discovery expression matrix contains non-finite values.")
    absent = sorted(set(GENES) - set(matrix.index))
    if absent:
        raise ValueError(f"Validated hub genes absent from discovery matrix: {absent}")

    gene_sets, n_pathways, pathway_release, gene_release = fetch_kegg_sets()
    print(
        f"[KEGG] {n_pathways} pathways listed, {len(gene_sets)} usable gene sets",
        flush=True,
    )

    all_results: list[pd.DataFrame] = []
    summary_rows = []
    for gene in GENES:
        scores = matrix.apply(
            lambda values: float(
                spearmanr(matrix.loc[gene].to_numpy(), values.to_numpy()).statistic
            ),
            axis=1,
        )
        ranking = (
            pd.DataFrame({"gene": scores.index, "spearman_rho": scores.values})
            .dropna()
            .sort_values(["spearman_rho", "gene"], ascending=[False, True])
        )
        ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
        ranking.to_csv(OUT / f"validated7_kegg_spearman_ranking_{gene}.csv", index=False)

        prerank_result = gseapy.prerank(
            rnk=ranking[["gene", "spearman_rho"]],
            gene_sets=gene_sets,
            min_size=MIN_SIZE,
            max_size=MAX_SIZE,
            permutation_num=PERMUTATIONS,
            seed=SEED,
            threads=8,
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        table = prerank_result.res2d.copy()
        table["gene"] = gene
        table = table.rename(
            columns={
                "Term": "pathway",
                "NES": "NES",
                "NOM p-val": "p_value",
                "FDR q-val": "FDR",
                "Lead_genes": "leading_edge_genes",
            }
        )
        table["gene_set_size"] = table["pathway"].map(
            lambda pathway: len(gene_sets.get(pathway, set()) & set(ranking["gene"]))
        )
        table["leading_edge_genes"] = (
            table.get("leading_edge_genes", pd.Series(dtype=str)).fillna("").astype(str)
        )
        empirical_floor = 1.0 / (PERMUTATIONS + 1)
        table["NES"] = pd.to_numeric(table["NES"], errors="coerce")
        table["p_value"] = pd.to_numeric(table["p_value"], errors="coerce").clip(
            lower=empirical_floor
        )
        table["FDR"] = pd.to_numeric(table["FDR"], errors="coerce").clip(
            lower=empirical_floor
        )
        table["empirical_floor_applied"] = (
            pd.to_numeric(prerank_result.res2d["NOM p-val"], errors="coerce").eq(0)
            | pd.to_numeric(prerank_result.res2d["FDR q-val"], errors="coerce").eq(0)
        )
        table["passes_fdr_0_05"] = table["FDR"] < FDR_CUTOFF
        table = table[
            [
                "gene",
                "pathway",
                "NES",
                "p_value",
                "FDR",
                "gene_set_size",
                "leading_edge_genes",
                "empirical_floor_applied",
                "passes_fdr_0_05",
            ]
        ]
        table.to_csv(OUT / f"validated7_kegg_gsea_{gene}.csv", index=False)
        all_results.append(table)

        n_sig = int(table["passes_fdr_0_05"].sum())
        print(f"[KEGG] {gene}: {n_sig} pathways pass FDR < 0.05", flush=True)
        summary_rows.append(
            {
                "gene": gene,
                "samples": len(samples),
                "ranked_genes": len(ranking),
                "KEGG_pathways_tested": len(table),
                "significant_pathways_FDR_lt_0_05": n_sig,
                "permutations": PERMUTATIONS,
                "seed": SEED,
                "status": "significant pathways found"
                if n_sig
                else "no pathway passes FDR < 0.05",
            }
        )

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUT / "validated7_kegg_gsea_all_genes.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "validated7_kegg_gsea_summary.csv", index=False)

    metadata_out = {
        "query_date": date.today().isoformat(),
        "expression_cohort": "GSE125498",
        "sample_count": len(samples),
        "LPD_sample_count": int((metadata["stage_binary"] == 1).sum()),
        "SPD_sample_count": int((metadata["stage_binary"] == 0).sum()),
        "rank_metric": "Spearman correlation of each hub gene against every gene; ties retained",
        "pathway_database": "KEGG human pathways via KEGG REST",
        "KEGG_pathway_database_updated": pathway_release,
        "KEGG_gene_database_updated": gene_release,
        "gsea_software": f"GSEApy {gseapy.__version__}",
        "permutations": PERMUTATIONS,
        "min_gene_set_size": MIN_SIZE,
        "max_gene_set_size": MAX_SIZE,
        "fdr_cutoff": FDR_CUTOFF,
        "seed": SEED,
        "query_genes": GENES,
    }
    (OUT / "validated7_kegg_gsea_metadata.json").write_text(
        json.dumps(metadata_out, indent=2), encoding="utf-8"
    )
    return combined, summary, metadata_out


def plot_dsigdb(results: pd.DataFrame) -> None:
    """Publication bar chart of top DSigDB candidate drug signatures."""
    plot_rows = results.head(15).copy()
    plot_rows = plot_rows.sort_values("combined_score", ascending=True)
    log10_combined = np.log10(plot_rows["combined_score"].clip(lower=1.0))

    def clean(term: str) -> str:
        """Turn 'raloxifene CTD 00007367' into 'raloxifene (CTD 00007367)'."""
        match = re.match(
            r"^(.*?)\s+(CTD\s+\d+|PC3\s+(?:UP|DOWN)|HL60\s+(?:UP|DOWN)"
            r"|MCF7\s+(?:UP|DOWN)|BOSS)$",
            term.strip(),
        )
        if match:
            return f"{match.group(1)} ({match.group(2)})"
        return term.strip()

    labels = [f"{clean(t)}\n[{g}]" for t, g in zip(plot_rows["dsigdb_term"], plot_rows["overlapping_genes"])]
    colors = [
        COLOR_SIG if sig else "#7AA6D6"
        for sig in plot_rows["significant_fdr_lt_0_05"]
    ]

    fig, ax = plt.subplots(figsize=(11.5, 7.6), dpi=300)
    y_pos = np.arange(len(plot_rows))
    ax.barh(y_pos, log10_combined, color=colors, height=0.7, edgecolor="none", alpha=0.92)

    for i, (_, row) in enumerate(plot_rows.iterrows()):
        ax.text(
            log10_combined.iloc[i] + 0.03,
            i,
            f"P = {row['p_value']:.2e} | FDR = {row['adjusted_p_value']:.2f}",
            va="center",
            fontsize=8,
            color="#333333",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.6)
    ax.set_xlabel(r"$\log_{10}$(DSigDB combined score)", fontsize=11, fontweight="bold")
    ax.set_title(
        "Candidate drug signatures (DSigDB) for the seven validated hub genes\n"
        "TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, CST7",
        fontsize=11.5,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlim(0, max(log10_combined) * 1.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle="--", alpha=0.35)

    handles = [
        Rectangle((0, 0), 1, 1, color=COLOR_SIG),
        Rectangle((0, 0), 1, 1, color="#7AA6D6"),
    ]
    ax.legend(
        handles,
        [f"FDR < {FDR_CUTOFF}", "Nominal P only (exploratory)"],
        loc="lower right",
        frameon=True,
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(OUT / "validated7_dsigdb_drug_candidates.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "validated7_dsigdb_drug_candidates.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[plot] validated7_dsigdb_drug_candidates.png / .pdf written", flush=True)


def _shorten(pathway: str) -> str:
    name = pathway.split(" - Homo sapiens")[0]
    return name if len(name) <= 42 else name[:39] + "..."


def plot_kegg_gsea(combined: pd.DataFrame) -> None:
    """Per-gene KEGG GSEA panel: top activated and suppressed pathways."""
    fig, axes = plt.subplots(
        nrows=len(GENES), ncols=2,
        figsize=(15.5, 2.1 * len(GENES)),
        gridspec_kw={"width_ratios": [1, 1]},
    )

    for row_idx, gene in enumerate(GENES):
        sig = combined[(combined["gene"] == gene) & (combined["passes_fdr_0_05"])].copy()
        activated = sig[sig["NES"] > 0].nlargest(6, "NES").sort_values("NES")
        suppressed = sig[sig["NES"] < 0].nsmallest(6, "NES").sort_values("NES")

        for col_idx, (data, color, title) in enumerate(
            [
                (activated, COLOR_ACTIVATED, "Activated (NES > 0)"),
                (suppressed, COLOR_SUPPRESSED, "Suppressed (NES < 0)"),
            ]
        ):
            ax = axes[row_idx, col_idx]
            if data.empty:
                ax.text(0.5, 0.5, "No pathway passes FDR < 0.05",
                        ha="center", va="center", fontsize=8.5, color="#666666")
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                y = np.arange(len(data))
                ax.barh(y, data["NES"], color=color, alpha=0.9, height=0.68, edgecolor="none")
                ax.set_yticks(y)
                ax.set_yticklabels([_shorten(p) for p in data["pathway"]], fontsize=7.4)
                for i, (_, r) in enumerate(data.iterrows()):
                    offset = 0.05 if r["NES"] > 0 else -0.05
                    ha = "left" if r["NES"] > 0 else "right"
                    ax.text(r["NES"] + offset, i, f"{r['NES']:.2f}", va="center",
                            ha=ha, fontsize=6.9, color="#333333")
                ax.axvline(0, color="#333333", lw=0.8)
                xmax = max(abs(data["NES"]).max() * 1.28, 1.0)
                ax.set_xlim(-xmax, xmax)
                ax.grid(axis="x", linestyle=":", alpha=0.35)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.tick_params(axis="x", labelsize=7)
            ax.set_title(f"{gene} — {title}", loc="left", fontsize=8.8, fontweight="bold")

    fig.suptitle(
        "Single-gene KEGG preranked GSEA for the seven validated hub genes\n"
        f"GSE125498 discovery (33 samples); significant terms only (FDR < {FDR_CUTOFF}); "
        f"{PERMUTATIONS} permutations",
        fontsize=12,
        fontweight="bold",
        y=0.995,
    )
    fig.subplots_adjust(left=0.135, right=0.975, top=0.93, bottom=0.035,
                        hspace=0.85, wspace=0.62)
    fig.savefig(OUT / "validated7_kegg_gsea_panel.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "validated7_kegg_gsea_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[plot] validated7_kegg_gsea_panel.png / .pdf written", flush=True)


def main() -> None:
    dsigdb_results = run_dsigdb()
    plot_dsigdb(dsigdb_results)

    combined, summary, metadata_out = run_kegg_gsea()
    plot_kegg_gsea(combined)

    print("\n=== KEGG GSEA summary ===")
    print(summary.to_string(index=False))
    print("\nAll validated-seven outputs written to:", OUT)


if __name__ == "__main__":
    main()
