import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "results" / "revision" / "ml_39_matrisome_wgcna"
RESULTS_PATH = OUTPUT_DIR / "ml_39_gene_results.csv"
SUMMARY_PATH = OUTPUT_DIR / "ml_39_gene_summary.json"
PNG_PATH = OUTPUT_DIR / "ml_39_gene_selection_plot.png"
PDF_PATH = OUTPUT_DIR / "ml_39_gene_selection_plot.pdf"

results = pd.read_csv(RESULTS_PATH)
summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
required = {
    "Gene",
    "Matrisome_DEG_Status",
    "LASSO_selected",
    "SVM_RFE_selected",
    "RF_selected",
    "Total_Votes_out_of_3",
    "Consensus_At_Least_2_of_3",
    "Discovery_log2FC",
}
if not required.issubset(results.columns):
    raise ValueError(f"Results table is missing required columns: {sorted(required - set(results.columns))}")
if results["Gene"].duplicated().any():
    raise ValueError("Results table contains duplicate genes.")
if len(results) != 39 or int(results["Consensus_At_Least_2_of_3"].sum()) != 9:
    raise ValueError("Results do not match the verified 39-input, 9-consensus analysis.")
if int((results["Total_Votes_out_of_3"] == 3).sum()) != 2:
    raise ValueError("Expected two strict 3-of-3 consensus genes.")

results = results.sort_values(
    ["Total_Votes_out_of_3", "Discovery_log2FC", "Gene"],
    ascending=[False, False, True],
).reset_index(drop=True)
genes = results["Gene"].tolist()
y = np.arange(len(results))
consensus = results["Consensus_At_Least_2_of_3"].astype(bool).to_numpy()
strict = results["Total_Votes_out_of_3"].eq(3).to_numpy()

method_columns = [
    ("LASSO", "LASSO_selected", "#2F6B9A"),
    ("SVM-RFE", "SVM_RFE_selected", "#D17A22"),
    ("RF", "RF_selected", "#4B8B62"),
]

fig = plt.figure(figsize=(10.5, 12.2), facecolor="white")
grid = fig.add_gridspec(
    nrows=1,
    ncols=3,
    width_ratios=[1.55, 1.15, 2.15],
    left=0.12,
    right=0.97,
    top=0.84,
    bottom=0.10,
    wspace=0.04,
)
gene_ax = fig.add_subplot(grid[0, 0])
selection_ax = fig.add_subplot(grid[0, 1], sharey=gene_ax)
effect_ax = fig.add_subplot(grid[0, 2], sharey=gene_ax)

for ax in (gene_ax, selection_ax, effect_ax):
    ax.set_ylim(len(results) - 0.5, -0.5)
    ax.set_yticks(y)
    ax.tick_params(axis="y", length=0, labelleft=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E8EDF2", linewidth=0.55)

gene_ax.set_xlim(0, 1)
gene_ax.set_xticks([])
gene_ax.set_title("Gene", loc="left", fontsize=11, fontweight="bold", pad=12)
for idx, row in results.iterrows():
    if consensus[idx]:
        gene_ax.axhspan(idx - 0.48, idx + 0.48, color="#FFF4D6", zorder=0)
    label = row["Gene"] + ("  ★" if strict[idx] else "")
    gene_ax.text(
        0.02,
        idx,
        label,
        ha="left",
        va="center",
        fontsize=9.2,
        fontstyle="italic",
        fontweight="bold" if consensus[idx] else "normal",
        color="#172B4D",
    )

selection_ax.set_xlim(-0.45, 3.45)
selection_ax.set_xticks([0, 1, 2, 3])
selection_ax.set_xticklabels(["LASSO", "SVM-\nRFE", "RF", "Votes"], fontsize=8.5)
selection_ax.xaxis.tick_top()
selection_ax.tick_params(axis="x", length=0, pad=5)
selection_ax.set_title("Method selection", loc="left", fontsize=11, fontweight="bold", pad=30)
selection_ax.set_yticklabels([])
for idx, row in results.iterrows():
    if consensus[idx]:
        selection_ax.axhspan(idx - 0.48, idx + 0.48, color="#FFF4D6", zorder=0)
    for x_pos, (_, column, color) in enumerate(method_columns):
        selected = bool(row[column])
        selection_ax.scatter(
            x_pos,
            idx,
            s=47 if selected else 25,
            color=color if selected else "#D7DEE6",
            edgecolors="white" if selected else "none",
            linewidths=0.45,
            zorder=3,
        )
    votes = int(row["Total_Votes_out_of_3"])
    selection_ax.text(
        3,
        idx,
        f"{votes}/3",
        ha="center",
        va="center",
        fontsize=8.8,
        fontweight="bold" if consensus[idx] else "normal",
        color="#172B4D" if consensus[idx] else "#64748B",
    )
selection_ax.set_xticks([0, 1, 2, 3])
selection_ax.set_xticklabels(["LASSO", "SVM-\nRFE", "RF", "Votes"], fontsize=8.5)
selection_ax.axvline(2.55, color="#CBD5E1", linewidth=0.8)

effect_ax.set_xlim(-1.55, 1.9)
effect_ax.set_xticks([-1, 0, 1])
effect_ax.set_xlabel("Discovery log₂ fold change\n(long-term vs short-term PD)", fontsize=8.5, labelpad=8)
effect_ax.xaxis.tick_top()
effect_ax.xaxis.set_label_position("top")
effect_ax.tick_params(axis="x", length=0, pad=5, labelsize=8.5)
effect_ax.set_title("GSE125498 expression effect", loc="left", fontsize=11, fontweight="bold", pad=12)
effect_ax.set_yticklabels([])
effect_ax.axvline(0, color="#64748B", linewidth=0.9, zorder=1)
for idx, row in results.iterrows():
    if consensus[idx]:
        effect_ax.axhspan(idx - 0.48, idx + 0.48, color="#FFF4D6", zorder=0)
    log_fc = float(row["Discovery_log2FC"])
    color = "#238B57" if log_fc >= 0 else "#B84A4A"
    effect_ax.barh(
        idx,
        log_fc,
        height=0.56,
        left=0,
        color=color,
        alpha=0.92,
        zorder=2,
    )

fig.suptitle(
    "Machine-learning feature selection from the WGCNA–Matrisome intersection",
    fontsize=14,
    fontweight="bold",
    y=0.965,
)
fig.text(
    0.5,
    0.925,
    "39 input genes (28 upregulated, 11 downregulated)  |  "
    "Consensus: 9 genes selected by ≥2 of 3 methods  |  Strict 3/3: 2",
    ha="center",
    va="center",
    fontsize=9.7,
    color="#334155",
)

legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="none",
           markersize=7, label=label)
    for label, _, color in method_columns
]
legend_handles += [
    Patch(facecolor="#238B57", label="Positive log₂FC"),
    Patch(facecolor="#B84A4A", label="Negative log₂FC"),
    Patch(facecolor="#FFF4D6", edgecolor="#E7C56A", label="Consensus hub (≥2/3)"),
    Line2D([0], [0], marker="$★$", color="none", label="Strict 3/3 hub", markersize=10),
]
fig.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.045),
    ncol=4,
    frameon=False,
    fontsize=8.5,
    columnspacing=1.6,
    handletextpad=0.5,
)
fig.text(
    0.5,
    0.015,
    "Feature selection and expression effects are from the same discovery dataset; "
    "this figure does not represent independent validation.",
    ha="center",
    va="center",
    fontsize=8,
    color="#64748B",
)

fig.savefig(PNG_PATH, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(PDF_PATH, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"Saved {PNG_PATH.relative_to(ROOT)}")
print(f"Saved {PDF_PATH.relative_to(ROOT)}")
print(f"Plotted {len(genes)} genes; consensus={int(consensus.sum())}; strict={int(strict.sum())}")
