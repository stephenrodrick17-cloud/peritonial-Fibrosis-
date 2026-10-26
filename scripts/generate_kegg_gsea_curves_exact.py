"""
Generate mathematically exact, fully authentic GSEA Running Enrichment Score Plots:
Left: Top Upregulated KEGG Pathways in LPD
Right: Top Downregulated KEGG Pathways in LPD (Upregulated in SPD)

Refined layout: legend on right panel moved to upper right (clean and free from curves).
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

OUT_DIR = "results/revision/wgcna_paper_figures"
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load data
rank_df = pd.read_csv("results/paper_style/Fig2_DEG_moderated_t_ranking_Entrez.csv")
kegg_up = pd.read_csv("results/paper_style/Fig2D_KEGG_up_LPD_FDR05.csv")
kegg_down = pd.read_csv("results/paper_style/Fig2E_KEGG_up_SPD_FDR05.csv")

with open("results/revision/wgcna_paper_figures/top12_kegg_pathway_entrez.json") as f:
    pw_genes = json.load(f)

ranked_entrez = list(rank_df["ENTREZID"].astype(str))
ranked_metric = rank_df["moderated_t"].values
N_genes = len(ranked_entrez)

# Exact GSEA running sum function
def compute_exact_gsea(gene_set, ranked_genes, weights, p=1.0):
    indicator = np.array([1 if g in gene_set else 0 for g in ranked_genes])
    nh = indicator.sum()
    if nh == 0:
        return np.zeros(len(ranked_genes)), []
    hit_weights = (np.abs(weights) ** p) * indicator
    nr = hit_weights.sum()
    p_hit = np.cumsum(hit_weights) / nr
    p_miss = np.cumsum(1 - indicator) / (len(ranked_genes) - nh)
    es = p_hit - p_miss
    hit_indices = np.where(indicator == 1)[0]
    return es, hit_indices

# Top 6 pathways
up_pathways = kegg_up.head(6).copy()
down_pathways = kegg_down.head(6).copy()

up_colors = ["#E63385", "#00A0E9", "#E6B800", "#008000", "#8B4513", "#D92121"]
down_colors = ["#E63385", "#00A0E9", "#E6B800", "#008000", "#8B4513", "#D92121"]

fig = plt.figure(figsize=(15.5, 8.2), dpi=300)
gs = gridspec.GridSpec(3, 2, height_ratios=[3.6, 1.1, 1.2], hspace=0.08, wspace=0.18)

x_ranks = np.arange(N_genes)

# ==================== LEFT: UPREGULATED IN LPD ====================
ax_up_curve = fig.add_subplot(gs[0, 0])
ax_up_hits = fig.add_subplot(gs[1, 0], sharex=ax_up_curve)
ax_up_rank = fig.add_subplot(gs[2, 0], sharex=ax_up_curve)

for idx, (_, row) in enumerate(up_pathways.iterrows()):
    col = up_colors[idx]
    pid = row["ID"]
    gset = set(pw_genes[pid])
    es_curve, hit_indices = compute_exact_gsea(gset, ranked_entrez, ranked_metric)
    
    lbl = f"{row['Description']} (NES={row['NES']:+.2f})"
    ax_up_curve.plot(x_ranks, es_curve, color=col, lw=2.0, label=lbl)
    ax_up_hits.vlines(hit_indices, idx, idx + 0.85, color=col, lw=0.9, alpha=0.85)

ax_up_curve.set_ylabel("Running Enrichment Score", fontsize=11, fontweight="bold")
ax_up_curve.axhline(0, color="#64748B", lw=0.8, linestyle=":")
ax_up_curve.legend(loc="upper right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#E2E8F0", fontsize=8.2)
ax_up_curve.grid(True, linestyle=":", alpha=0.35)
ax_up_curve.tick_params(labelbottom=False)
ax_up_curve.set_title("Upregulated KEGG Pathways in LPD (n=13) vs SPD (n=20)", fontsize=12.5, fontweight="bold", pad=10)
ax_up_curve.set_ylim(-0.05, 0.80)

ax_up_hits.set_yticks([])
ax_up_hits.set_ylim(-0.2, len(up_pathways))
ax_up_hits.tick_params(labelbottom=False)
ax_up_hits.grid(False)

ax_up_rank.fill_between(x_ranks, 0, ranked_metric, color="#94A3B8", alpha=0.45)
ax_up_rank.plot(x_ranks, ranked_metric, color="#475569", lw=1.0)
ax_up_rank.axhline(0, color="#1E293B", lw=0.8)
ax_up_rank.set_ylabel("Ranked Metric\n(moderated t)", fontsize=9.5, fontweight="bold")
ax_up_rank.set_xlabel("Rank in Ordered Dataset (11,139 Entrez Genes)", fontsize=10.5, fontweight="bold")
ax_up_rank.set_xlim(0, N_genes)
ax_up_rank.set_ylim(-6.0, 7.5)
ax_up_rank.grid(True, linestyle=":", alpha=0.35)

# ==================== RIGHT: DOWNREGULATED IN LPD (UP IN SPD) ====================
ax_down_curve = fig.add_subplot(gs[0, 1])
ax_down_hits = fig.add_subplot(gs[1, 1], sharex=ax_down_curve)
ax_down_rank = fig.add_subplot(gs[2, 1], sharex=ax_down_curve)

for idx, (_, row) in enumerate(down_pathways.iterrows()):
    col = down_colors[idx]
    pid = row["ID"]
    gset = set(pw_genes[pid])
    es_curve, hit_indices = compute_exact_gsea(gset, ranked_entrez, ranked_metric)
    
    lbl = f"{row['Description']} (NES={row['NES']:+.2f})"
    ax_down_curve.plot(x_ranks, es_curve, color=col, lw=2.0, label=lbl)
    ax_down_hits.vlines(hit_indices, idx, idx + 0.85, color=col, lw=0.9, alpha=0.85)

ax_down_curve.set_ylabel("Running Enrichment Score", fontsize=11, fontweight="bold")
ax_down_curve.axhline(0, color="#64748B", lw=0.8, linestyle=":")
# In the right panel, curves dip negative into [-0.45, 0.0] and trough around rank 7,500-8,500.
# Upper right (or upper middle) has maximum blank space. Let's use upper right with y-limit up to +0.08
ax_down_curve.legend(loc="upper right", frameon=True, framealpha=0.92, facecolor="white", edgecolor="#E2E8F0", fontsize=8.2)
ax_down_curve.grid(True, linestyle=":", alpha=0.35)
ax_down_curve.tick_params(labelbottom=False)
ax_down_curve.set_title("Downregulated KEGG Pathways in LPD (n=13) vs SPD (n=20)", fontsize=12.5, fontweight="bold", pad=10)
ax_down_curve.set_ylim(-0.52, 0.08)

ax_down_hits.set_yticks([])
ax_down_hits.set_ylim(-0.2, len(down_pathways))
ax_down_hits.tick_params(labelbottom=False)
ax_down_hits.grid(False)

ax_down_rank.fill_between(x_ranks, 0, ranked_metric, color="#94A3B8", alpha=0.45)
ax_down_rank.plot(x_ranks, ranked_metric, color="#475569", lw=1.0)
ax_down_rank.axhline(0, color="#1E293B", lw=0.8)
ax_down_rank.set_ylabel("Ranked Metric\n(moderated t)", fontsize=9.5, fontweight="bold")
ax_down_rank.set_xlabel("Rank in Ordered Dataset (11,139 Entrez Genes)", fontsize=10.5, fontweight="bold")
ax_down_rank.set_xlim(0, N_genes)
ax_down_rank.set_ylim(-6.0, 7.5)
ax_down_rank.grid(True, linestyle=":", alpha=0.35)

for ext in ["png", "pdf", "svg"]:
    out_file = os.path.join(OUT_DIR, f"KEGG_GSEA_running_score_LPD_vs_SPD_exact.{ext}")
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_file}")

plt.close()
print("Exact GSEA curve generation complete!")
