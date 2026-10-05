"""
scripts/31_replot_string_with_legend.py
Re-draws results/figures/D_string_ppi_network.{png,pdf} from the saved STRING edge table
(results/tables/D_string_ppi_network_source_data.csv) - no new API query - and adds a legend.
"""
import sys
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append("scripts")
from string_legend import add_string_legend

HUBS = ["ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"]
edges = pd.read_csv("results/tables/D_string_ppi_network_source_data.csv")
cent = pd.read_csv("results/tables/D_hub_centrality.csv")

G = nx.Graph()
G.add_nodes_from(HUBS)
for _, r in edges.iterrows():
    G.add_edge(r["Node_A"], r["Node_B"], weight=float(r["STRING_Functional_Score"]))

wdeg = dict(G.degree(weight="weight"))
top = max(wdeg, key=wdeg.get)
size_fn = lambda d: 500 + 350 * d
width_fn = lambda s: s * 3.5

fig, ax = plt.subplots(figsize=(12.5, 8), facecolor="white")
pos = nx.spring_layout(G, seed=42, k=0.95)
nx.draw_networkx_edges(G, pos, width=[width_fn(G[u][v]["weight"]) for u, v in G.edges()],
                       alpha=0.6, edge_color="#64748b", ax=ax)
nx.draw_networkx_nodes(G, pos, node_size=[size_fn(G.degree(n)) for n in G.nodes()],
                       node_color=["#ef4444" if n == top else ("#3b82f6" if G.degree(n) > 0 else "#94a3b8") for n in G.nodes()],
                       edgecolors="#1e293b", linewidths=1.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=11, font_weight="bold", ax=ax)
legs = add_string_legend(ax, top, wdeg[top], size_fn, width_fn)
ax.set_title("Protein-Protein Interaction Network (STRING v12.5)", fontsize=13, fontweight="bold", pad=15)
iso = [n for n in G.nodes() if G.degree(n) == 0]
fig.text(0.40, 0.02, f"STRING v12.5 functional network, 11 hub genes, {G.number_of_edges()} edges (combined score >= 0.400). "
         f"Top edge: {edges.iloc[0]['Node_A']}-{edges.iloc[0]['Node_B']} ({edges.iloc[0]['STRING_Functional_Score']:.3f}). "
         f"Isolated: {', '.join(iso)}.", ha="center", fontsize=9, color="#334155")
ax.axis("off")
plt.tight_layout(rect=(0, 0.04, 1, 1))
plt.savefig("results/figures/D_string_ppi_network.png", dpi=300, bbox_inches="tight", bbox_extra_artists=legs)
plt.savefig("results/figures/D_string_ppi_network.pdf", bbox_inches="tight", bbox_extra_artists=legs)
print(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()} top_weighted_degree={top} ({wdeg[top]:.3f}) isolated={iso}")
print("Check vs D_hub_centrality.csv weighted degree:", dict(zip(cent.Gene, cent.Weighted_Degree))[top])
print("Saved results/figures/D_string_ppi_network.png/.pdf with legend")
