"""
Generate clean, beautifully spaced STRING PPI network visualization for the 43 Matrisome-WGCNA genes.
Uses kamada_kawai / spring layout with collision avoidance, distinct community spacing,
and clear labels.
"""

import os
import json
import urllib.request
import urllib.parse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

ROOT = r"d:/Peritoneal Project"
OUT_DIR = os.path.join(ROOT, "results", "revision", "matrisome_string_ppi")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load data & edges
df_43 = pd.read_csv(os.path.join(ROOT, "results", "revision", "ml_43_matrisome_wgcna", "ml_43_matrisome_results.csv"))
genes43 = df_43["Gene"].tolist()
gene_dict = df_43.set_index("Gene").to_dict(orient="index")

df_edges_func = pd.read_csv(os.path.join(OUT_DIR, "string_43_functional_edges.csv"))
df_edges_phys = pd.read_csv(os.path.join(OUT_DIR, "string_43_physical_edges.csv"))

def generate_clean_network(df_edges, title, filename_prefix, net_desc):
    G = nx.Graph()
    for g in genes43:
        info = gene_dict[g]
        G.add_node(
            g,
            status=info["DEG_Status"],
            logFC=info["logFC"],
            hub=(info["Model_Vote_Count"] >= 2),
            votes=info["Model_Vote_Count"]
        )
    for _, r in df_edges.iterrows():
        G.add_edge(r["gene_A"], r["gene_B"], weight=r["score"])

    degrees = dict(G.degree())
    
    # Separate connected and isolated
    connected = [n for n in G.nodes() if degrees[n] > 0]
    isolated = [n for n in G.nodes() if degrees[n] == 0]

    # Layout for connected component: spring with large k and expansion
    subG = G.subgraph(connected)
    pos_conn = nx.spring_layout(subG, k=0.9, iterations=500, seed=42, weight="weight")
    
    # Scale connected nodes to fit inside a circle of radius 0.65
    coords = np.array(list(pos_conn.values()))
    c_min = coords.min(axis=0)
    c_max = coords.max(axis=0)
    c_center = (c_min + c_max) / 2.0
    c_span = max((c_max - c_min).max() / 2.0, 0.01)

    pos = {}
    for n in connected:
        norm_coord = (pos_conn[n] - c_center) / c_span * 0.58
        pos[n] = norm_coord

    # Place isolated nodes in a regular outer ring (radius 0.88)
    n_iso = len(isolated)
    for i, n in enumerate(isolated):
        angle = 2.0 * np.pi * i / n_iso
        pos[n] = np.array([0.88 * np.cos(angle), 0.88 * np.sin(angle)])

    fig, ax = plt.subplots(figsize=(15.5, 15.5), facecolor="#FFFFFF")
    ax.set_facecolor("#FFFFFF")
    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-1.08, 1.08)
    ax.axis("off")

    # Draw edges
    edges = list(G.edges(data=True))
    if edges:
        weights = [d["weight"] for _, _, d in edges]
        widths = [1.2 + 2.5 * (w - 0.4) / 0.6 for w in weights]
        edge_colors = ["#94A3B8" if w < 0.70 else "#1E293B" for w in weights]
        edge_alphas = [0.55 if w < 0.70 else 0.85 for w in weights]
        
        for (u, v, d), w_val, col, a_val in zip(edges, widths, edge_colors, edge_alphas):
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color=col, lw=w_val, alpha=a_val, zorder=2)

    # Draw nodes
    for n in G.nodes():
        x, y = pos[n]
        is_up = (G.nodes[n]["status"] == "Upregulated")
        is_hub = G.nodes[n]["hub"]
        deg = degrees[n]
        
        # Color & size
        fill_color = "#E11D48" if is_up else "#2563EB"
        r_size = 520 + 75 * deg + (220 if is_hub else 0)
        
        # Halo for ML hubs
        if is_hub:
            ax.scatter([x], [y], s=r_size + 380, color="#F59E0B", alpha=0.9, zorder=3)
            ax.scatter([x], [y], s=r_size, color=fill_color, edgecolor="#FFFFFF", linewidth=2.0, zorder=4)
        else:
            ax.scatter([x], [y], s=r_size, color=fill_color, edgecolor="#FFFFFF", linewidth=1.5, zorder=4)

        # Label placement: above node or inside node depending on size
        ax.text(x, y, n, fontsize=8.8 if is_hub else 7.8, fontweight="bold",
                ha="center", va="center", color="#FFFFFF", zorder=5)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor="#E11D48", edgecolor="#FFFFFF", label="Upregulated (n=33)"),
        mpatches.Patch(facecolor="#2563EB", edgecolor="#FFFFFF", label="Downregulated (n=10)"),
        Line2D([0], [0], marker="o", color="w", markeredgecolor="#F59E0B", markeredgewidth=3.0,
               markerfacecolor="#E11D48", markersize=14, label="ML Hub Gene (Votes >= 2, n=10)"),
        Line2D([0], [0], color="#1E293B", lw=3.0, label="High Confidence Edge (Score >= 0.70)"),
        Line2D([0], [0], color="#94A3B8", lw=1.5, label="Medium Confidence Edge (0.40 - 0.70)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", frameon=True,
              facecolor="#F8FAFC", edgecolor="#CBD5E1", fontsize=9.8, title="Node & Edge Attributes", title_fontsize=11.0)

    # Title & summary stats
    plt.title(title, fontsize=15.0, fontweight="bold", pad=20, color="#0F172A")
    top_deg = sorted([(n, degrees[n]) for n in G.nodes()], key=lambda x: x[1], reverse=True)[:5]
    top_str = ", ".join([f"{n} (k={d})" for n, d in top_deg if d > 0])
    
    summary_text = (
        f"STRING PPI Summary ({net_desc})\n"
        f"• Network Nodes: 43 (33 Upregulated, 10 Downregulated)\n"
        f"• Connected Subnetwork: {len(connected)} genes | Isolated: {len(isolated)} genes\n"
        f"• Total Network Edges: {len(df_edges)}\n"
        f"• Interaction Score Cutoff: >= 0.400 (medium-to-high confidence)\n"
        f"• Top Connected Hubs: {top_str}"
    )
    ax.text(0.02, 0.02, summary_text, transform=ax.transAxes, fontsize=9.5,
            verticalalignment="bottom", linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.6", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.96))

    out_png = os.path.join(OUT_DIR, f"{filename_prefix}.png")
    out_pdf = os.path.join(OUT_DIR, f"{filename_prefix}.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Generated clean network: {out_png}")

generate_clean_network(
    df_edges_func,
    "STRING Protein-Protein Interaction Network (43 Matrisome-WGCNA Genes)\nFunctional Associations (Interaction Score >= 0.400)",
    "STRING_PPI_43_genes_functional",
    "Functional Associations"
)

generate_clean_network(
    df_edges_phys,
    "STRING Protein-Protein Interaction Network (43 Matrisome-WGCNA Genes)\nPhysical Protein Interactions (Interaction Score >= 0.400)",
    "STRING_PPI_43_genes_physical",
    "Physical Interactions"
)
