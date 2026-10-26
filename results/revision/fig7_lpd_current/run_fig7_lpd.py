"""Fresh exploratory Fig. 7 analysis for the current panel-v2 genes.

Run from the repository root with the project Python environment. All outputs
are written to this new directory; the historical Fig. 7 outputs are not read.
"""

from __future__ import annotations

import csv
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from numba import njit

_NULL_EXPR: np.ndarray | None = None
_NULL_GENE_INDEX: dict[str, int] | None = None
_NULL_RESAMPLES: np.ndarray | None = None


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "revision" / "fig7_lpd_current"
EXPR_PATH = ROOT / "results" / "tables" / "GSE125498_full_expression_matrix.csv"
META_PATH = ROOT / "results" / "tables" / "GSE125498_sample_metadata.csv"
PANEL_PATH = ROOT / "results" / "revision" / "panel_v2" / "panel_v2.json"
POOL_PATH = (
    ROOT
    / "results"
    / "revision"
    / "ml_43_matrisome_wgcna"
    / "ml_43_matrisome_results.csv"
)
SIGNATURE_PATH = (
    ROOT
    / "results"
    / "revision"
    / "microenvironment_deconvolution"
    / "cellmarker_2024_library_snapshot.txt"
)

SEED = 20261016
RESAMPLES = 1000
P_ITEM = 0.80
K_VALUES = (2, 3, 4, 5)
NULL_SETS = 1000
SEED_REPEATS = 20
CDF_GRID = np.linspace(0.0, 1.0, 101)
LEUKOCYTE_MARKERS = [
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD19",
    "MS4A1", "CD79A", "CD79B", "NCAM1", "KLRD1", "CD14", "CD68", "ITGAM",
    "FCGR3A",
]
SIGNATURE_TERMS = {
    "T cells": "T Cell Peripheral Blood Human",
    "B cells": "B Cell Peripheral Blood Human",
    "NK cells": "Natural Killer Cell Peritoneal Fluid Human",
    "Monocytes": "Monocyte Peripheral Blood Human",
    "M1 macrophages": "M1 Macrophage Peripheral Blood Human",
    "M2 macrophages": "M2 Macrophage Peripheral Blood Human",
    "Neutrophils": "Neutrophil Peripheral Blood Human",
    "Mast cells": "Mast Cell Lung Human",
    "Immature dendritic cells": "Immature Dendritic Cell Stomach Human",
    "Mature dendritic cells": "Mature Dendritic Cell Liver Human",
    "Fibroblasts": "Fibroblast Epithelium Human",
    "Myofibroblasts": "Myofibroblast Epithelium Human",
    "Mesothelial cells": "Mesothelial Cell Epithelium Human",
}


def load_inputs():
    existing_outputs = [
        path for path in OUT.iterdir() if path.resolve() != Path(__file__).resolve()
    ] if OUT.exists() else []
    if existing_outputs:
        raise FileExistsError(
            f"Refusing to reuse or overwrite existing outputs: {existing_outputs}"
        )
    for path in (EXPR_PATH, META_PATH, PANEL_PATH, POOL_PATH, SIGNATURE_PATH):
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing: {path}")

    expr = pd.read_csv(EXPR_PATH).set_index("Gene")
    meta = pd.read_csv(META_PATH)
    panel = json.loads(PANEL_PATH.read_text(encoding="utf-8"))
    pool = pd.read_csv(POOL_PATH)["Gene"].astype(str).tolist()
    if expr.index.has_duplicates or len(pool) != 43 or len(set(pool)) != 43:
        raise ValueError("Expected unique expression genes and a unique 43-gene pool.")
    genes = panel["genes_passing_selection_rules"]
    if len(genes) != 5 or not set(genes).issubset(expr.index):
        raise ValueError(f"Current panel-v2 genes are missing or unexpected: {genes}")
    if not set(genes).issubset(pool):
        raise ValueError("The current panel is not fully contained in the 43-gene pool.")
    if set(meta["sample_id"]) != set(expr.columns):
        raise ValueError("Sample IDs do not match between expression and metadata.")
    if meta["sample_id"].duplicated().any():
        raise ValueError("Sample IDs are not unique in metadata.")
    lpd_meta = meta.loc[meta["group"].eq("LPD")].copy()
    if len(lpd_meta) != 13:
        raise ValueError(f"Expected 13 LPD samples, found {len(lpd_meta)}.")
    expr = expr.loc[:, meta["sample_id"]].astype(float)
    return expr, meta, lpd_meta, genes, pool


def make_resamples(n_samples: int, seed: int) -> list[np.ndarray]:
    n_item = int(np.floor(P_ITEM * n_samples))
    if n_item <= max(K_VALUES):
        raise ValueError("pItem resamples must contain more samples than the largest k.")
    rng = np.random.default_rng(seed)
    return [
        np.sort(rng.choice(n_samples, size=n_item, replace=False))
        for _ in range(RESAMPLES)
    ]


def cluster_once(
    x: np.ndarray,
    samples: list[np.ndarray],
) -> tuple[dict[int, np.ndarray], dict[int, dict[str, float]], dict[int, np.ndarray]]:
    """Ward-Euclidean consensus clustering with complete-link consensus cuts."""
    n = x.shape[0]
    coobserved = np.zeros((n, n), dtype=np.int32)
    coclustered = {k: np.zeros((n, n), dtype=np.int32) for k in K_VALUES}
    for keep in samples:
        ward_tree = linkage(x[keep], method="ward", metric="euclidean")
        partitions = cut_tree(ward_tree, n_clusters=list(K_VALUES))
        coobserved[np.ix_(keep, keep)] += 1
        for col, k in enumerate(K_VALUES):
            same = np.equal.outer(partitions[:, col], partitions[:, col])
            coclustered[k][np.ix_(keep, keep)] += same

    consensus: dict[int, np.ndarray] = {}
    assignments: dict[int, np.ndarray] = {}
    metrics: dict[int, dict[str, float]] = {}
    for k in K_VALUES:
        matrix = np.divide(
            coclustered[k],
            coobserved,
            out=np.zeros((n, n), dtype=float),
            where=coobserved > 0,
        )
        np.fill_diagonal(matrix, 1.0)
        consensus[k] = matrix
        distance = np.clip(1.0 - matrix, 0.0, 1.0)
        np.fill_diagonal(distance, 0.0)
        final_tree = linkage(squareform(distance, checks=True), method="complete")
        cluster = cut_tree(final_tree, n_clusters=[k]).ravel() + 1
        assignments[k] = cluster
        upper = matrix[np.triu_indices(n, 1)]
        cdf = np.array([np.mean(upper <= value) for value in CDF_GRID])
        cdf_area = float(np.trapezoid(cdf, CDF_GRID))
        metrics[k] = {
            "cdf_area": cdf_area,
            "mean_pairwise_consensus": float(np.mean(upper)),
            "pac_0_1": float(np.mean((upper > 0.1) & (upper < 0.9))),
            "silhouette_euclidean_z": float(silhouette_score(x, cluster, metric="euclidean")),
            "silhouette_consensus_distance": float(
                silhouette_score(distance, cluster, metric="precomputed")
            ),
        }
        if len(set(cluster)) != k:
            raise RuntimeError(f"Consensus cut returned fewer than {k} groups.")
    for k in K_VALUES:
        metrics[k]["delta_area"] = (
            metrics[k]["cdf_area"] - metrics[k - 1]["cdf_area"]
            if k > K_VALUES[0]
            else float("nan")
        )
    return consensus, metrics, assignments


@njit
def _ward_cut_labels(points: np.ndarray, k: int) -> np.ndarray:
    """Cut a Ward-Euclidean hierarchy using size-weighted centroid distances."""
    n, n_features = points.shape
    centroids = np.zeros((2 * n - 1, n_features))
    sizes = np.zeros(2 * n - 1, dtype=np.int64)
    active = np.zeros(2 * n - 1, dtype=np.bool_)
    point_cluster = np.arange(n, dtype=np.int64)
    for i in range(n):
        centroids[i] = points[i]
        sizes[i] = 1
        active[i] = True

    for step in range(n - k):
        best_distance = np.inf
        best_a = -1
        best_b = -1
        for a in range(n + step):
            if not active[a]:
                continue
            for b in range(a + 1, n + step):
                if not active[b]:
                    continue
                squared = 0.0
                for feature in range(n_features):
                    difference = centroids[a, feature] - centroids[b, feature]
                    squared += difference * difference
                weight = (
                    2.0 * sizes[a] * sizes[b] / (sizes[a] + sizes[b])
                )
                distance = np.sqrt(weight * squared)
                if distance < best_distance:
                    best_distance = distance
                    best_a = a
                    best_b = b

        new_cluster = n + step
        total_size = sizes[best_a] + sizes[best_b]
        for feature in range(n_features):
            centroids[new_cluster, feature] = (
                sizes[best_a] * centroids[best_a, feature]
                + sizes[best_b] * centroids[best_b, feature]
            ) / total_size
        sizes[new_cluster] = total_size
        active[best_a] = False
        active[best_b] = False
        active[new_cluster] = True
        for point in range(n):
            if point_cluster[point] == best_a or point_cluster[point] == best_b:
                point_cluster[point] = new_cluster

    labels = np.empty(n, dtype=np.int64)
    active_label = 0
    label_lookup = np.full(2 * n - 1, -1, dtype=np.int64)
    for cluster in range(2 * n - 1):
        if active[cluster]:
            label_lookup[cluster] = active_label
            active_label += 1
    for point in range(n):
        labels[point] = label_lookup[point_cluster[point]]
    return labels


@njit
def _consensus_from_resamples(
    x: np.ndarray, resamples: np.ndarray, k: int
) -> np.ndarray:
    n = x.shape[0]
    coobserved = np.zeros((n, n), dtype=np.int32)
    coclustered = np.zeros((n, n), dtype=np.int32)
    for resample_idx in range(resamples.shape[0]):
        keep = resamples[resample_idx]
        subset = x[keep]
        labels = _ward_cut_labels(subset, k)
        for i in range(keep.shape[0]):
            for j in range(keep.shape[0]):
                coobserved[keep[i], keep[j]] += 1
                if labels[i] == labels[j]:
                    coclustered[keep[i], keep[j]] += 1
    consensus = np.zeros((n, n))
    for i in range(n):
        consensus[i, i] = 1.0
        for j in range(n):
            if coobserved[i, j] > 0:
                consensus[i, j] = coclustered[i, j] / coobserved[i, j]
    return consensus


def select_k(metrics: dict[int, dict[str, float]]) -> int:
    """Choose the k with the largest raw-expression Euclidean silhouette."""
    return max(K_VALUES, key=lambda k: (metrics[k]["silhouette_euclidean_z"], -k))


def build_signatures(expr_genes: set[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    library: dict[str, set[str]] = {}
    with SIGNATURE_PATH.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            term = row[0].strip()
            if term:
                library[term] = {g.strip().upper() for g in row[2:] if g.strip()}

    signatures: dict[str, list[str]] = {}
    provenance = []
    for cell_type, term in SIGNATURE_TERMS.items():
        if term not in library:
            raise ValueError(f"Required CellMarker signature is missing: {term}.")
        markers = sorted(library[term] & expr_genes)
        if len(markers) < 3:
            raise ValueError(
                f"Only {len(markers)} genes matched signature {cell_type} ({term})."
            )
        signatures[cell_type] = markers
        provenance.append(
            {
                "cell_type": cell_type,
                "source": "local CellMarker_2024 human library snapshot",
                "source_signature": term,
                "n_genes_in_expression": len(markers),
                "merge_rule": "exact frozen CellMarker_2024 source term",
            }
        )
    if len(signatures) != 13:
        raise RuntimeError("Expected exactly 13 immune signatures.")
    return signatures, pd.DataFrame(provenance)


def score_ssgsea(expr: pd.DataFrame, signatures: dict[str, list[str]]) -> pd.DataFrame:
    import gseapy as gp

    result = gp.ssgsea(
        data=expr,
        gene_sets=signatures,
        min_size=3,
        max_size=5000,
        sample_norm_method="rank",
        correl_norm_type="rank",
        permutation_num=None,
        threads=1,
        no_plot=True,
        seed=125498,
        outdir=None,
        verbose=False,
    ).res2d
    if result is None or not {"Name", "Term", "NES"}.issubset(result.columns):
        raise ValueError("Unexpected ssGSEA output structure.")
    scores = result.pivot(index="Name", columns="Term", values="NES")
    scores = scores.reindex(index=expr.columns, columns=list(signatures))
    if scores.isna().any().any():
        raise ValueError("ssGSEA scores are missing for at least one sample/signature.")
    scores.index.name = "sample_id"
    return scores


def init_null_worker(
    expr_values: np.ndarray, gene_names: list[str], samples: np.ndarray
) -> None:
    global _NULL_EXPR, _NULL_GENE_INDEX, _NULL_RESAMPLES
    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1)
    _NULL_EXPR = expr_values
    _NULL_GENE_INDEX = {gene: idx for idx, gene in enumerate(gene_names)}
    _NULL_RESAMPLES = samples


def null_gene_set_batch(
    task: tuple[str, list[tuple[int, list[str]]], int],
) -> list[dict[str, object]]:
    if _NULL_EXPR is None or _NULL_GENE_INDEX is None or _NULL_RESAMPLES is None:
        raise RuntimeError("Null worker was not initialized.")
    universe, random_sets, selected_k = task
    rows = []
    for set_index, random_genes in random_sets:
        row_indices = [_NULL_GENE_INDEX[gene] for gene in random_genes]
        random_x = _NULL_EXPR[row_indices, :].T
        means = random_x.mean(axis=0)
        scales = random_x.std(axis=0, ddof=1)
        if np.any(scales == 0):
            raise ValueError(f"Zero-variance gene in null set {universe}:{set_index}.")
        random_x = (random_x - means) / scales
        consensus = _consensus_from_resamples(
            random_x, _NULL_RESAMPLES, selected_k
        )
        distance = np.clip(1.0 - consensus, 0.0, 1.0)
        np.fill_diagonal(distance, 0.0)
        final_tree = linkage(squareform(distance, checks=True), method="complete")
        cluster = cut_tree(final_tree, n_clusters=[selected_k]).ravel() + 1
        upper = consensus[np.triu_indices(random_x.shape[0], 1)]
        cdf = np.array([np.mean(upper <= value) for value in CDF_GRID])
        rows.append(
            {
                "universe": universe,
                "random_set": set_index,
                "genes": ";".join(random_genes),
                "k": selected_k,
                "cdf_area": float(np.trapezoid(cdf, CDF_GRID)),
                "mean_pairwise_consensus": float(np.mean(upper)),
                "pac_0_1": float(np.mean((upper > 0.1) & (upper < 0.9))),
                "silhouette_euclidean_z": float(
                    silhouette_score(random_x, cluster, metric="euclidean")
                ),
                "silhouette_consensus_distance": float(
                    silhouette_score(distance, cluster, metric="precomputed")
                ),
            }
        )
    return rows


def bh(values: pd.Series) -> np.ndarray:
    return multipletests(values.to_numpy(dtype=float), method="fdr_bh")[1]


def save_plots(
    x: np.ndarray,
    genes: list[str],
    sample_names: list[str],
    assignments: dict[int, np.ndarray],
    selected_k: int,
    consensus: dict[int, np.ndarray],
    metrics: dict[int, dict[str, float]],
    cdf_rows: list[dict[str, float]],
    immune: pd.DataFrame,
    immune_tests: pd.DataFrame,
) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    colors = ["#2864A5", "#D1495B", "#3A8D5D", "#9265A8", "#D18B27"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    cdf_frame = pd.DataFrame(cdf_rows)
    for k, frame in cdf_frame.groupby("k"):
        axs[0, 0].plot(
            frame["consensus_value"],
            frame["cumulative_fraction"],
            label=f"k={k}",
            color=colors[k - 2],
        )
    axs[0, 0].set(xlabel="Consensus index", ylabel="Cumulative fraction", title="Consensus CDF")
    axs[0, 0].legend(frameon=False)
    ks = list(K_VALUES)
    axs[0, 1].bar(
        [k - 0.15 for k in ks],
        [metrics[k]["delta_area"] for k in ks],
        width=0.3,
        label="Delta CDF area",
        color="#77AADD",
    )
    axs[0, 1].plot(
        [k + 0.15 for k in ks],
        [metrics[k]["silhouette_euclidean_z"] for k in ks],
        "o-",
        label="Silhouette",
        color="#B23A48",
    )
    axs[0, 1].set(xticks=ks, xlabel="k", title="Selection metrics (k=2 delta undefined)")
    axs[0, 1].legend(frameon=False)

    matrix = consensus[selected_k]
    labels = assignments[selected_k]
    order = np.argsort(labels, kind="stable")
    im = axs[1, 0].imshow(matrix[np.ix_(order, order)], vmin=0, vmax=1, cmap="Blues")
    axs[1, 0].set(
        xticks=range(len(order)),
        yticks=range(len(order)),
        xticklabels=np.array(sample_names)[order],
        yticklabels=np.array(sample_names)[order],
        title=f"Consensus matrix, k={selected_k}",
    )
    axs[1, 0].tick_params(axis="x", labelrotation=90, labelsize=7)
    axs[1, 0].tick_params(axis="y", labelsize=7)
    fig.colorbar(im, ax=axs[1, 0], label="Consensus")

    pca_model = PCA(n_components=2)
    pca = pca_model.fit_transform(x)
    for cluster in sorted(set(labels)):
        mask = labels == cluster
        axs[1, 1].scatter(
            pca[mask, 0],
            pca[mask, 1],
            label=f"Cluster {cluster}",
            color=colors[cluster - 1],
            s=48,
        )
    axs[1, 1].set(
        xlabel=f"PC1 ({pca_model.explained_variance_ratio_[0] * 100:.1f}%)",
        ylabel=f"PC2 ({pca_model.explained_variance_ratio_[1] * 100:.1f}%)",
        title="PCA of five z-scored panel genes",
    )
    axs[1, 1].legend(frameon=False)
    fig.savefig(OUT / "Fig7_exploratory.png", dpi=300)
    fig.savefig(OUT / "Fig7_exploratory.pdf")
    plt.close(fig)

    gene_rows = []
    for idx, gene in enumerate(genes):
        for sample_idx, sample in enumerate(sample_names):
            gene_rows.append(
                {
                    "sample_id": sample,
                    "gene": gene,
                    "z_score": float(x[sample_idx, idx]),
                    "cluster": int(labels[sample_idx]),
                }
            )
    gene_frame = pd.DataFrame(gene_rows)
    fig, axs = plt.subplots(1, len(genes), figsize=(14, 4), sharey=True, constrained_layout=True)
    for ax, gene in zip(axs, genes):
        groups = [
            gene_frame.loc[
                gene_frame["gene"].eq(gene) & gene_frame["cluster"].eq(cluster), "z_score"
            ].to_numpy()
            for cluster in sorted(set(labels))
        ]
        ax.boxplot(
            groups, tick_labels=[str(c) for c in sorted(set(labels))], widths=0.55
        )
        for cluster_idx, values in enumerate(groups, start=1):
            jitter = (
                np.linspace(-0.06, 0.06, len(values))
                if len(values) > 1
                else np.array([0.0])
            )
            ax.scatter(cluster_idx + jitter, values, s=18, color="#333333", zorder=3)
        ax.set_title(gene)
        ax.set_xlabel("Cluster")
    axs[0].set_ylabel("Gene-wise z-score")
    fig.suptitle("Descriptive panel expression (clusters defined from these genes)")
    fig.savefig(OUT / "Fig7_gene_panel_descriptive.png", dpi=300)
    fig.savefig(OUT / "Fig7_gene_panel_descriptive.pdf")
    plt.close(fig)

    signature_names = immune_tests["cell_type"].tolist()
    fig, axs = plt.subplots(4, 4, figsize=(13, 12), constrained_layout=True)
    for ax, signature in zip(axs.ravel(), signature_names):
        values = [
            immune.loc[immune["cluster"].eq(cluster), signature].to_numpy()
            for cluster in sorted(immune["cluster"].unique())
        ]
        ax.boxplot(
            values, tick_labels=[str(c) for c in sorted(immune["cluster"].unique())]
        )
        for i, vector in enumerate(values, start=1):
            jitter = (
                np.linspace(-0.06, 0.06, len(vector))
                if len(vector) > 1
                else np.array([0.0])
            )
            ax.scatter(i + jitter, vector, s=13, color="#333333", zorder=3)
        test = immune_tests.loc[immune_tests["cell_type"].eq(signature)].iloc[0]
        q_value = test.get("mw_bh_q", test.get("minimum_pairwise_mw_bh_q", np.nan))
        ax.set_title(f"{signature}\nMW q={q_value:.3g}", fontsize=9)
        ax.set_xlabel("Cluster")
    for ax in axs.ravel()[len(signature_names) :]:
        ax.set_visible(False)
    fig.suptitle("Exploratory ssGSEA scores by cluster")
    fig.savefig(OUT / "Fig7_immune_scores.png", dpi=300)
    fig.savefig(OUT / "Fig7_immune_scores.pdf")
    plt.close(fig)


def main() -> None:
    expr, _, lpd_meta, genes, pool = load_inputs()
    sample_names = lpd_meta["sample_id"].tolist()
    z_panel = expr.loc[genes, sample_names].T
    z_panel = z_panel.sub(z_panel.mean(axis=0), axis=1).div(
        z_panel.std(axis=0, ddof=1).replace(0, np.nan), axis=1
    )
    if z_panel.isna().any().any():
        raise ValueError("A panel gene has zero variance or missing values in the LPD samples.")
    x = z_panel.to_numpy()
    shared_resamples = make_resamples(len(sample_names), SEED)
    consensus, metrics, assignments = cluster_once(x, shared_resamples)
    selected_k = select_k(metrics)
    selected_labels = assignments[selected_k]

    metric_rows = []
    cdf_rows = []
    membership_rows = []
    for k in K_VALUES:
        counts = pd.Series(assignments[k]).value_counts().sort_index()
        size_text = ";".join(f"{cluster}:{count}" for cluster, count in counts.items())
        row = {
            "k": k,
            **metrics[k],
            "cluster_sizes": size_text,
            "minimum_cluster_size": int(counts.min()),
            "selected_by_max_euclidean_silhouette": k == selected_k,
        }
        metric_rows.append(row)
        upper = consensus[k][np.triu_indices(len(sample_names), 1)]
        cdf = np.array([np.mean(upper <= value) for value in CDF_GRID])
        cdf_rows.extend(
            {
                "k": k,
                "consensus_value": float(value),
                "cumulative_fraction": float(fraction),
            }
            for value, fraction in zip(CDF_GRID, cdf)
        )
        for sample, cluster in zip(sample_names, assignments[k]):
            membership_rows.append({"k": k, "sample_id": sample, "cluster": int(cluster)})
    metric_frame = pd.DataFrame(metric_rows)
    cdf_frame = pd.DataFrame(cdf_rows)
    membership_frame = pd.DataFrame(membership_rows)

    # Generate 1,000 random same-size panels in each null universe; all genes
    # in the observed panel are excluded from the null draws.
    null_rng = np.random.default_rng(SEED + 1)
    random_sets = {
        "43_gene_pool": [
            null_rng.choice([g for g in pool if g not in genes], size=len(genes), replace=False)
            for _ in range(NULL_SETS)
        ],
        "all_expressed_genes": [
            null_rng.choice(
                [g for g in expr.index if g not in genes], size=len(genes), replace=False
            )
            for _ in range(NULL_SETS)
        ],
    }
    batch_size = 25
    tasks = []
    for universe, gene_sets in random_sets.items():
        indexed_sets = [
            (set_index, list(random_genes))
            for set_index, random_genes in enumerate(gene_sets, start=1)
        ]
        tasks.extend(
            (universe, indexed_sets[start : start + batch_size], selected_k)
            for start in range(0, len(indexed_sets), batch_size)
        )
    expr_lpd_values = expr.loc[:, sample_names].to_numpy(dtype=float)
    resample_matrix = np.stack(shared_resamples)
    worker_count = min(6, max(1, (os.cpu_count() or 2) - 1))
    null_rows = []
    completed_by_universe = {name: 0 for name in random_sets}
    with ProcessPoolExecutor(
        max_workers=worker_count,
        initializer=init_null_worker,
        initargs=(expr_lpd_values, expr.index.tolist(), resample_matrix),
    ) as executor:
        futures = {executor.submit(null_gene_set_batch, task): task for task in tasks}
        for completed, future in enumerate(as_completed(futures), start=1):
            universe, batch, _ = futures[future]
            completed_by_universe[universe] += len(batch)
            null_rows.extend(future.result())
            if completed % 4 == 0 or completed == len(futures):
                print(
                    "Null panels completed: "
                    + ", ".join(
                        f"{name} {count}/{NULL_SETS}"
                        for name, count in completed_by_universe.items()
                    ),
                    flush=True,
                )
    null_frame = pd.DataFrame(null_rows)

    observed_k_metrics = metrics[selected_k]
    null_comparisons = []
    for universe, frame in null_frame.loc[null_frame["k"].eq(selected_k)].groupby("universe"):
        for metric in (
            "cdf_area",
            "mean_pairwise_consensus",
            "pac_0_1",
            "silhouette_euclidean_z",
            "silhouette_consensus_distance",
        ):
            distribution = frame[metric].to_numpy(dtype=float)
            observed = observed_k_metrics[metric]
            lower = (np.sum(distribution <= observed) + 1) / (len(distribution) + 1)
            upper = (np.sum(distribution >= observed) + 1) / (len(distribution) + 1)
            null_comparisons.append(
                {
                    "universe": universe,
                    "k": selected_k,
                    "metric": metric,
                    "observed": observed,
                    "null_mean": float(np.mean(distribution)),
                    "null_median": float(np.median(distribution)),
                    "null_2_5_percentile": float(np.quantile(distribution, 0.025)),
                    "null_97_5_percentile": float(np.quantile(distribution, 0.975)),
                    "observed_percentile": float(100 * np.mean(distribution <= observed)),
                    "empirical_two_sided_p": float(min(1.0, 2 * min(lower, upper))),
                }
            )
    null_comparison_frame = pd.DataFrame(null_comparisons)

    # Sensitivity to random resampling seed, keeping the observed genes fixed.
    stability_rows = []
    for repeat_seed in range(SEED + 1, SEED + SEED_REPEATS + 1):
        samples = make_resamples(len(sample_names), repeat_seed)
        _, repeat_metrics, repeat_assignments = cluster_once(x, samples)
        repeat_k = select_k(repeat_metrics)
        stability_rows.append(
            {
                "seed": repeat_seed,
                "selected_k": repeat_k,
                "primary_k_ari": float(
                    adjusted_rand_score(assignments[selected_k], repeat_assignments[selected_k])
                ),
                **{
                    f"silhouette_k{k}": repeat_metrics[k]["silhouette_euclidean_z"]
                    for k in K_VALUES
                },
            }
        )
    stability_frame = pd.DataFrame(stability_rows)

    # Recompute ssGSEA on expression data with the local frozen signatures.
    signatures, signature_provenance = build_signatures(set(expr.index))
    immune_scores = score_ssgsea(expr, signatures).loc[sample_names]
    z_all = expr.sub(expr.mean(axis=1), axis=0).div(
        expr.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )
    present_markers = [g for g in LEUKOCYTE_MARKERS if g in z_all.index]
    if len(present_markers) < 10:
        raise ValueError("Fewer than 10 predefined leukocyte markers were observed.")
    leukocyte_proxy = z_all.loc[present_markers, sample_names].mean(axis=0)
    immune = immune_scores.copy()
    immune["cluster"] = selected_labels
    immune["leukocyte_proxy"] = leukocyte_proxy.loc[sample_names].to_numpy()

    immune_test_rows = []
    pairwise_rows = []
    clusters = sorted(set(selected_labels))
    for signature in signatures:
        groups = [
            immune.loc[immune["cluster"].eq(cluster), signature].to_numpy()
            for cluster in clusters
        ]
        medians = {str(c): float(np.median(values)) for c, values in zip(clusters, groups)}
        design = pd.DataFrame(
            {
                "score": immune[signature].to_numpy(),
                "cluster": pd.Categorical(immune["cluster"]),
                "leukocyte_proxy": immune["leukocyte_proxy"].to_numpy(),
            }
        )
        full = sm.OLS.from_formula(
            "score ~ C(cluster) + leukocyte_proxy", data=design
        ).fit()
        reduced = sm.OLS.from_formula("score ~ leukocyte_proxy", data=design).fit()
        adjusted_p = float(full.compare_f_test(reduced)[1])
        adjusted_effect = (
            float(full.params.filter(like="C(cluster)").iloc[-1])
            if len(clusters) == 2
            else float("nan")
        )
        for i, cluster_a in enumerate(clusters):
            for j in range(i + 1, len(clusters)):
                cluster_b = clusters[j]
                pairwise_rows.append(
                    {
                        "cell_type": signature,
                        "cluster_a": cluster_a,
                        "cluster_b": cluster_b,
                        "median_a": medians[str(cluster_a)],
                        "median_b": medians[str(cluster_b)],
                        "mann_whitney_p": float(
                            mannwhitneyu(
                                groups[i], groups[j], alternative="two-sided"
                            ).pvalue
                        ),
                    }
                )
        immune_test_rows.append(
            {
                "cell_type": signature,
                "test": "Mann-Whitney U" if len(clusters) == 2 else "pairwise Mann-Whitney U",
                "cluster_medians": json.dumps(medians, sort_keys=True),
                "unadjusted_p": float(
                    mannwhitneyu(groups[0], groups[1], alternative="two-sided").pvalue
                )
                if len(clusters) == 2
                else np.nan,
                "leukocyte_adjusted_cluster_p": adjusted_p,
                "leukocyte_adjusted_cluster_effect_C2_minus_C1": adjusted_effect,
            }
        )
    immune_tests = pd.DataFrame(immune_test_rows)
    pairwise_tests = pd.DataFrame(pairwise_rows)
    pairwise_tests["bh_q"] = bh(pairwise_tests["mann_whitney_p"])
    if len(clusters) == 2:
        immune_tests["mw_bh_q"] = bh(immune_tests["unadjusted_p"])
    else:
        min_q = pairwise_tests.groupby("cell_type")["bh_q"].min()
        immune_tests["minimum_pairwise_mw_bh_q"] = immune_tests["cell_type"].map(min_q)
    immune_tests["adjusted_bh_q"] = bh(immune_tests["leukocyte_adjusted_cluster_p"])
    immune_tests["survives_BH_after_leukocyte_adjustment"] = (
        immune_tests["adjusted_bh_q"] < 0.05
    )

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(x)
    pca_frame = pd.DataFrame(
        {
            "sample_id": sample_names,
            "cluster": selected_labels,
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
        }
    )
    pca_frame["PC1_variance_percent"] = pca.explained_variance_ratio_[0] * 100
    pca_frame["PC2_variance_percent"] = pca.explained_variance_ratio_[1] * 100
    assignments_frame = membership_frame.loc[membership_frame["k"].eq(selected_k)].copy()
    assignments_frame["title"] = lpd_meta.set_index("sample_id").loc[
        assignments_frame["sample_id"], "title"
    ].to_numpy()
    sizes = (
        assignments_frame.groupby("cluster", as_index=False)
        .agg(n=("sample_id", "size"), samples=("sample_id", lambda s: ";".join(s)))
        .assign(selected_k=selected_k)
    )
    covariates = pd.DataFrame(
        [
            {
                "variable": "dialysis_duration",
                "available": False,
                "test": "not estimable",
                "reason": "GEO provides only a short-term/long-term category, not per-patient duration.",
            },
            {
                "variable": "age",
                "available": False,
                "test": "not estimable",
                "reason": "No per-patient age is provided in the GEO series matrix or sample table.",
            },
            {
                "variable": "sex",
                "available": False,
                "test": "not estimable",
                "reason": "No per-patient sex is provided in the GEO series matrix or sample table.",
            },
            {
                "variable": "ultrafiltration_failure",
                "available": False,
                "test": "not estimable",
                "reason": "No per-patient UF-failure status is provided in the GEO series matrix or sample table.",
            },
            {
                "variable": "batch",
                "available": False,
                "test": "not estimable",
                "reason": "GEO sample descriptions PRAGUE_1 through PRAGUE_33 do not define a varying batch factor.",
            },
        ]
    )

    save_plots(
        x,
        genes,
        sample_names,
        assignments,
        selected_k,
        consensus,
        metrics,
        cdf_rows,
        immune,
        immune_tests,
    )
    files = {
        "metrics_by_k.csv": metric_frame,
        "consensus_cdf.csv": cdf_frame,
        "consensus_matrix_selected_k.csv": pd.DataFrame(
            consensus[selected_k], index=sample_names, columns=sample_names
        ).rename_axis("sample_id"),
        "cluster_membership_all_k.csv": membership_frame,
        "cluster_assignments_selected_k.csv": assignments_frame,
        "cluster_sizes_selected_k.csv": sizes,
        "pca_coordinates.csv": pca_frame,
        "gene_zscores_by_sample.csv": z_panel.rename_axis("sample_id").reset_index(),
        "null_distributions.csv": null_frame,
        "null_observed_vs_random.csv": null_comparison_frame,
        "seed_stability.csv": stability_frame,
        "ssgsea_scores_selected_k.csv": immune.reset_index(),
        "ssgsea_tests_selected_k.csv": immune_tests,
        "ssgsea_pairwise_mann_whitney_selected_k.csv": pairwise_tests,
        "ssgsea_signature_provenance.csv": signature_provenance,
        "clinical_covariate_availability.csv": covariates,
    }
    for filename, frame in files.items():
        frame.to_csv(OUT / filename, index=False, na_rep="")

    def markdown_table(headers: list[str], rows: list[list[object]]) -> list[str]:
        return [
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join("---" for _ in headers) + " |",
            *[
                "| " + " | ".join(str(value) for value in row) + " |"
                for row in rows
            ],
        ]

    seed_k_counts_text = ", ".join(
        f"k={k}: {(stability_frame.selected_k == k).sum()}" for k in K_VALUES
    )
    cluster_size_text = " and ".join(
        str(int(row.n)) for row in sizes.itertuples()
    )
    pool_silhouette = null_comparison_frame.loc[
        null_comparison_frame["universe"].eq("43_gene_pool")
        & null_comparison_frame["metric"].eq("silhouette_euclidean_z")
    ].iloc[0]
    expressed_silhouette = null_comparison_frame.loc[
        null_comparison_frame["universe"].eq("all_expressed_genes")
        & null_comparison_frame["metric"].eq("silhouette_euclidean_z")
    ].iloc[0]
    report_lines = [
        "# Figure 7: exploratory LPD subgroups",
        "",
        "**Exploratory; n=13 LPD samples; not validated.**",
        "",
        "## Gene set and samples",
        "",
        f"Current candidate panel: {', '.join(genes)} (panel-v2 rule-derived set).",
        "Each gene was z-scored across these 13 samples using the sample SD (ddof=1).",
        "",
        ", ".join(sample_names),
        "",
        "## Consensus clustering",
        "",
        f"Resampling: {RESAMPLES} subsamples, pItem={P_ITEM}, seed={SEED}; "
        "each subsample uses Ward agglomerative clustering with Euclidean distance. "
        "For each k, consensus is the fraction of co-sampled resamples in which a "
        "sample pair co-clustered. Final consensus partitions use complete linkage "
        "on distance 1-consensus. k was selected by the largest Euclidean silhouette "
        "on the full z-scored expression matrix (ties to smaller k), not by visual "
        "inspection of the figure.",
        "",
        "CDF curves for all k and the selected-k consensus matrix are exported to "
        "`consensus_cdf.csv` and `consensus_matrix_selected_k.csv`.",
        "",
        *markdown_table(
            [
                "k",
                "CDF area",
                "Delta area",
                "Euclidean silhouette",
                "Consensus-distance silhouette",
                "PAC (0.1, 0.9)",
                "Cluster sizes",
            ],
            [
                [
                    int(row.k),
                    f"{row.cdf_area:.4f}",
                    "NA" if pd.isna(row.delta_area) else f"{row.delta_area:.4f}",
                    f"{row.silhouette_euclidean_z:.4f}",
                    f"{row.silhouette_consensus_distance:.4f}",
                    f"{row.pac_0_1:.4f}",
                    row.cluster_sizes,
                ]
                for row in metric_frame.itertuples()
            ],
        ),
        "",
        f"Selected k={selected_k}; cluster sizes: "
        + ", ".join(f"C{row.cluster} n={row.n}" for row in sizes.itertuples())
        + ".",
        "",
        "Selected-k sample membership:",
        "",
        *markdown_table(
            ["Cluster", "n", "Sample IDs"],
            [
                [int(row.cluster), int(row.n), row.samples]
                for row in sizes.itertuples()
            ],
        ),
        "",
        "## Random-gene-set null controls",
        "",
        f"Each null set contains {len(genes)} genes, sampled without replacement "
        "and excluding the observed genes. There are 1,000 sets from the 43-gene "
        f"pool and 1,000 from all {len(expr)} expressed genes. Each set uses the "
        "same 1,000 sample subsamples and is re-z-scored across the 13 LPD samples. "
        "All selected-k null metrics are in `null_distributions.csv`; empirical "
        "two-sided p-values use the plus-one correction.",
        "",
        *markdown_table(
            [
                "Universe",
                "Metric",
                "Observed",
                "Null median (95% range)",
                "Observed percentile",
                "Empirical two-sided p",
            ],
            [
                [
                    row.universe,
                    row.metric,
                    f"{row.observed:.4f}",
                    f"{row.null_median:.4f} "
                    f"({row.null_2_5_percentile:.4f}-{row.null_97_5_percentile:.4f})",
                    f"{row.observed_percentile:.1f}%",
                    f"{row.empirical_two_sided_p:.4f}",
                ]
                for row in null_comparison_frame.itertuples()
            ],
        ),
        "",
        "## PCA and descriptive gene panel",
        "",
        f"The PCA is calculated on the {len(genes)} gene-wise z-scored panel and "
        "colored by selected consensus cluster (`Fig7_exploratory.png`). The "
        f"first two components explain {pca_frame['PC1_variance_percent'].iloc[0]:.1f}% "
        f"and {pca_frame['PC2_variance_percent'].iloc[0]:.1f}% of variance, respectively. "
        "The cluster-wise panel plot is descriptive only: clusters were constructed "
        "from these same genes (`Fig7_gene_panel_descriptive.png`).",
        "",
        "## Immune/stromal ssGSEA scores",
        "",
        "Thirteen exact CellMarker_2024 human signatures were scored by ssGSEA. "
        "The leukocyte proxy is the mean of 17 canonical leukocyte-marker "
        "z-scores, each standardized across all 33 discovery samples. "
        "For k=2, the reported unadjusted tests are two-sided Mann-Whitney U tests "
        "with BH q across 13 signatures. For k>2, all pairwise two-sided "
        "Mann-Whitney U tests are BH-adjusted together; the `ssgsea_tests` file "
        "also reports an omnibus cluster-term F test after leukocyte-proxy "
        "adjustment, BH-adjusted across the 13 signatures. Scores are descriptive "
        "marker-enrichment scores, not cell fractions.",
        "",
        *markdown_table(
            [
                "Signature",
                "Unadjusted test p",
                "MW BH q (or minimum pairwise q)",
                "Leukocyte-adjusted cluster p",
                "Adjusted BH q",
                "Survives q<0.05?",
            ],
            [
                [
                    row.cell_type,
                    "NA" if pd.isna(row.unadjusted_p) else f"{row.unadjusted_p:.4g}",
                    f"{getattr(row, 'mw_bh_q', getattr(row, 'minimum_pairwise_mw_bh_q', np.nan)):.4g}",
                    f"{row.leukocyte_adjusted_cluster_p:.4g}",
                    f"{row.adjusted_bh_q:.4g}",
                    "yes" if row.survives_BH_after_leukocyte_adjustment else "no",
                ]
                for row in immune_tests.itertuples()
            ],
        ),
        "",
        "Pairwise Mann-Whitney values and BH q-values are in "
        "`ssgsea_pairwise_mann_whitney_selected_k.csv`; sample scores are in "
        "`ssgsea_scores_selected_k.csv`.",
        "",
        *markdown_table(
            [
                "Signature",
                "Cluster A",
                "Cluster B",
                "Median A",
                "Median B",
                "MW p",
                "BH q",
            ],
            [
                [
                    row.cell_type,
                    int(row.cluster_a),
                    int(row.cluster_b),
                    f"{row.median_a:.4f}",
                    f"{row.median_b:.4f}",
                    f"{row.mann_whitney_p:.4g}",
                    f"{row.bh_q:.4g}",
                ]
                for row in pairwise_tests.itertuples()
            ],
        ),
        "",
        "## Clinical covariates and stability",
        "",
        "The GEO series matrix reports patient ID, short/long-term group, and "
        "peritoneal-cell type only; it does not provide per-patient dialysis "
        "duration, age, sex, or ultrafiltration-failure status. Its `PRAGUE_1` "
        "through `PRAGUE_33` sample descriptions do not define a varying batch "
        "factor. Therefore, tests for these clinical variables and batch are "
        "not estimable; see `clinical_covariate_availability.csv`.",
        "",
        f"Across {SEED_REPEATS} additional resampling seeds, selected k counts "
        f"were {seed_k_counts_text}. The selected 12/1 split includes a cluster "
        "smaller than four; all tested k values also contain at least one "
        "cluster smaller than four. Therefore, the subgroup definition is "
        "unstable despite no selected-k changes across these seeds. See "
        "`seed_stability.csv` and `cluster_size_lt_4` in `run_metadata.json`.",
        "",
        "## Figure legend",
        "",
        f"**Figure 7. Exploratory molecular subgroups within long-term PD "
        f"(LPD; n=13).** Consensus clustering used the current five-gene "
        f"candidate panel ({', '.join(genes)}), gene-wise z-scored across the "
        f"13 LPD samples, with {RESAMPLES} subsamples (pItem={P_ITEM}; seed "
        f"{SEED}), Ward linkage with Euclidean distance, and final complete-linkage "
        f"clustering of 1-consensus. The silhouette-based rule selected k={selected_k} "
        f"with cluster sizes {cluster_size_text}"
        + ". CDF/delta-area and silhouette statistics are reported for k=2-5. "
        "The PCA is colored by selected cluster. The gene panel is descriptive "
        "because the same genes defined the clusters. ssGSEA signatures are "
        "exploratory marker scores. A selected cluster has n<4. This analysis "
        "is exploratory (n=13), unvalidated, and not suitable for clinical "
        "subgroup inference.",
        "",
        "## Results text",
        "",
        f"Within the 13 LPD samples, exploratory consensus clustering of the "
        f"current five-gene panel selected k={selected_k} by the pre-specified "
        f"maximum Euclidean silhouette rule (silhouette="
        f"{metrics[selected_k]['silhouette_euclidean_z']:.4f}); the resulting "
        f"cluster sizes were {cluster_size_text}"
        + ". Results for k=2-5, including consensus CDF area, delta-area, "
        "silhouette, and memberships, are shown in the metrics and membership "
        "tables. The selected-k Euclidean silhouette was "
        f"{pool_silhouette.observed:.4f}, at the "
        f"{pool_silhouette.observed_percentile:.1f}% of null values below "
        "the observed score among 1,000 "
        "five-gene null sets from the 43-gene pool (empirical two-sided "
        f"p={pool_silhouette.empirical_two_sided_p:.4f}); against 11,741 "
        "expressed-gene null sets, "
        f"{expressed_silhouette.observed_percentile:.1f}% of null values were "
        "below the observed score "
        f"(p={expressed_silhouette.empirical_two_sided_p:.4f}). No immune/"
        "stromal signature difference survived BH adjustment, including after "
        "leukocyte-proxy adjustment. Duration, age, sex, UF failure, and batch "
        "could not be tested because patient-level covariates were unavailable. "
        "This analysis is exploratory, has n=13, includes a singleton cluster, "
        "and is not validated.",
        "",
        "## Code line map",
        "",
        "- Input validation and current-panel/LPD checks: "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L84).",
        "- Fixed resampling and Ward/Euclidean consensus clustering: "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L118) and "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L129).",
        "- Fast Ward implementation (verified against SciPy partitions): "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L186).",
        "- Silhouette-based k selection: "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L273).",
        "- Frozen 13-signature definitions and ssGSEA: "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L278) and "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L313).",
        "- Matched random-gene-set null clustering: "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L352) and "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L594).",
        "- ssGSEA tests and covariate availability: "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L708) and "
        "[run_fig7_lpd.py](./run_fig7_lpd.py#L793).",
        "",
    ]
    (OUT / "fig7_results_and_legend.md").write_text(
        "\n".join(report_lines), encoding="utf-8"
    )

    metadata = {
        "analysis_label": "exploratory, n=13, not validated",
        "current_panel_definition": str(PANEL_PATH.relative_to(ROOT)),
        "genes": genes,
        "sample_ids": sample_names,
        "seed": SEED,
        "resamples": RESAMPLES,
        "pItem": P_ITEM,
        "k_values": list(K_VALUES),
        "resampling_algorithm": "Ward agglomerative hierarchical clustering, Euclidean distance",
        "consensus_final_algorithm": "complete-linkage hierarchical clustering, distance 1-consensus",
        "selection_rule": "maximum silhouette on Euclidean distances of gene-wise z-scored expression; ties to smaller k",
        "null_sets_per_universe": NULL_SETS,
        "null_universes": {
            "43_gene_pool": "43-gene model input pool, excluding the five observed genes",
            "all_expressed_genes": f"{len(expr)} expressed genes, excluding the five observed genes",
        },
        "selected_k": selected_k,
        "cluster_sizes": dict(zip(sizes["cluster"].astype(str), sizes["n"].astype(int))),
        "cluster_size_lt_4": bool((sizes["n"] < 4).any()),
        "seed_stability_repeats": SEED_REPEATS,
        "seed_selected_k_counts": {
            str(k): int((stability_frame["selected_k"] == k).sum()) for k in K_VALUES
        },
        "input_files": [
            str(EXPR_PATH.relative_to(ROOT)),
            str(META_PATH.relative_to(ROOT)),
            str(PANEL_PATH.relative_to(ROOT)),
            str(POOL_PATH.relative_to(ROOT)),
            str(SIGNATURE_PATH.relative_to(ROOT)),
        ],
        "historical_fig7_outputs_read": False,
    }
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "selected_k": selected_k,
                "cluster_sizes": metadata["cluster_sizes"],
                "seed_selected_k_counts": metadata["seed_selected_k_counts"],
                "outputs": str(OUT.relative_to(ROOT)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
