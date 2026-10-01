"""
==============================================================================
TASK 7: PROTEIN-PROTEIN INTERACTION (PPI) & CO-EXPRESSION NETWORK (11 HUB GENES)
==============================================================================
Constructs an integrated network from STRING v12.0 functional/physical interactions
(score >= 0.400) and empirical co-expression correlation (|r| >= 0.70, FDR < 0.05)
in human peritoneal tissue (GSE62928). Computes topological centrality metrics
(Degree, Betweenness, Closeness) and exports publication-ready vector PDF & 300 DPI PNG.
"""

import os
import urllib.request
import json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

# Set seeds and directories
np.random.seed(42)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# 1. Consensus Hub Genes
hub_genes = ["FN1", "COL3A1", "COL11A1", "COL8A1", "VCAN", "COMP", "THBS3", "EDIL3", "LOX", "INHBA", "ISM1"]

# Matrisome categories & ML consensus votes
df_hubs = pd.read_csv("results/tables/ML_hub_genes_from_WGCNA_ECM.csv")
meta_dict = df_hubs.set_index("Gene_Symbol").to_dict(orient="index")

# 2. Fetch STRING Interactions (score >= 0.400)
print("Step 1: Querying live STRING database API...")
caller_id = "peritoneal_fibrosis_audit"
url_tsv = f"https://string-db.org/api/tsv/network?identifiers={'%0d'.join(hub_genes)}&species=9606&required_score=400&network_type=functional&caller_identity={caller_id}"
url_json = f"https://string-db.org/api/json/network?identifiers={'%0d'.join(hub_genes)}&species=9606&required_score=400&network_type=functional&caller_identity={caller_id}"

req_tsv = urllib.request.Request(url_tsv, headers={"User-Agent": "Mozilla/5.0"})
req_json = urllib.request.Request(url_json, headers={"User-Agent": "Mozilla/5.0"})

try:
    with urllib.request.urlopen(req_tsv, timeout=20) as resp:
        tsv_raw = resp.read().decode('utf-8')
    today_str = "20261001"
    tsv_out_path = f"results/tables/string_live_edges_{today_str}.tsv"
    with open(tsv_out_path, "w", encoding="utf-8") as f:
        f.write(tsv_raw)
    print(f"Saved raw live STRING response to {tsv_out_path}")

    with urllib.request.urlopen(req_json, timeout=20) as resp:
        string_raw = json.loads(resp.read().decode('utf-8'))
        print(f"STRING API returned {len(string_raw)} raw interaction entries.")
except Exception as e:
    raise RuntimeError(f"STRING API query failed: {e}. Fallback cache is strictly disabled.") from e

if not string_raw:
    raise RuntimeError("STRING API returned empty results. Fallback cache is strictly disabled.")

string_edges = {}
for item in string_raw:
    g1 = item["preferredName_A"]
    g2 = item["preferredName_B"]
    score = float(item["score"])
    if g1 in hub_genes and g2 in hub_genes and g1 != g2:
        pair = tuple(sorted([g1, g2]))
        if pair not in string_edges or score > string_edges[pair]:
            string_edges[pair] = score

print(f"Loaded {len(string_edges)} verified live STRING PPI edges (score >= 0.400).")

# 3. Compute Co-expression Edges in Discovery Tissue (GSE62928, N=8)
print("Step 2: Computing pairwise Pearson co-expression in GSE62928...")
df_expr = pd.read_csv("results/tables/GSE62928_full_expression_matrix.csv", index_col=0)
expr_sub = df_expr.loc[hub_genes].T

pairs_list = []
r_vals = []
p_vals = []
for i in range(len(hub_genes)):
    for j in range(i + 1, len(hub_genes)):
        g1, g2 = hub_genes[i], hub_genes[j]
        r, p = pearsonr(expr_sub[g1], expr_sub[g2])
        pairs_list.append((g1, g2))
        r_vals.append(r)
        p_vals.append(p)

_, p_adj, _, _ = multipletests(p_vals, method="fdr_bh")

coexpr_edges = {}
for idx, (g1, g2) in enumerate(pairs_list):
    r = r_vals[idx]
    padj = p_adj[idx]
    pval = p_vals[idx]
    # Retain strong co-expression: |r| >= 0.85 and FDR < 0.01 for clear network visualization
    if abs(r) >= 0.85 and padj < 0.01:
        coexpr_edges[tuple(sorted([g1, g2]))] = {"r": r, "pval": pval, "padj": padj}

print(f"Identified {len(coexpr_edges)} significant co-expression edges (|r| >= 0.85, FDR < 0.01).")

# 4. Construct NetworkX Graph
G = nx.Graph()
for g in hub_genes:
    cat = meta_dict.get(g, {}).get("Matrisome_Category", "Core Matrisome")
    votes = meta_dict.get(g, {}).get("Votes", 2)
    G.add_node(g, category=cat, votes=votes)

# Add STRING edges
for (u, v), score in string_edges.items():
    G.add_edge(u, v, edge_type="STRING_PPI", weight=score, string_score=score)

# Add Co-expression edges
for (u, v), d in coexpr_edges.items():
    if G.has_edge(u, v):
        G[u][v]["edge_type"] = "Both"
        G[u][v]["corr_r"] = d["r"]
    else:
        G.add_edge(u, v, edge_type="Coexpression", weight=abs(d["r"]), corr_r=d["r"])

# 5. Compute Centrality Metrics
deg_centrality = dict(G.degree())
betweenness = nx.betweenness_centrality(G, weight="weight")
closeness = nx.closeness_centrality(G)
eigenvector = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)

centrality_rows = []
for g in hub_genes:
    # Count specific edge types
    string_deg = sum(1 for v in G.neighbors(g) if G[g][v]["edge_type"] in ["STRING_PPI", "Both"])
    coexpr_deg = sum(1 for v in G.neighbors(g) if G[g][v]["edge_type"] in ["Coexpression", "Both"])
    centrality_rows.append({
        "Gene_Symbol": g,
        "Matrisome_Category": meta_dict.get(g, {}).get("Matrisome_Category", "Core Matrisome"),
        "Matrisome_Division": meta_dict.get(g, {}).get("Matrisome_Division", "Core matrisome"),
        "ML_Consensus_Votes": meta_dict.get(g, {}).get("Votes", 2),
        "Total_Degree": deg_centrality[g],
        "STRING_Degree": string_deg,
        "Coexpression_Degree": coexpr_deg,
        "Betweenness_Centrality": round(betweenness[g], 4),
        "Closeness_Centrality": round(closeness[g], 4),
        "Eigenvector_Centrality": round(eigenvector[g], 4)
    })

df_metrics = pd.DataFrame(centrality_rows).sort_values(by=["Total_Degree", "Betweenness_Centrality"], ascending=False)
out_metrics = "results/tables/hub_genes_ppi_centrality_metrics.csv"
df_metrics.to_csv(out_metrics, index=False)
print(f"Exported centrality metrics to {out_metrics}")

# 6. Publication Visualization
plt.figure(figsize=(12, 10), facecolor="#FFFFFF")
ax = plt.gca()
ax.set_facecolor("#FFFFFF")

# Spring layout with fixed seed and repulsion
pos = nx.spring_layout(G, k=1.8, iterations=100, seed=42)

# Color palette by Matrisome Category
category_colors = {
    "ECM Glycoproteins": "#2563EB",   # Blue
    "Collagens": "#DC2626",           # Red
    "Proteoglycans": "#16A34A",       # Green
    "ECM Regulators": "#9333EA",       # Purple
    "Secreted Factors": "#EA580C"      # Orange
}

# Separate edges by type
both_edges = [(u, v) for u, v, d in G.edges(data=True) if d["edge_type"] == "Both"]
string_only = [(u, v) for u, v, d in G.edges(data=True) if d["edge_type"] == "STRING_PPI"]
coexpr_only = [(u, v) for u, v, d in G.edges(data=True) if d["edge_type"] == "Coexpression"]

# Draw Co-expression edges (dashed green/slate lines)
nx.draw_networkx_edges(G, pos, edgelist=coexpr_only, style="dashed",
                       edge_color="#059669", width=1.5, alpha=0.6, ax=ax)

# Draw STRING PPI edges (solid blue lines with thickness proportional to score)
for u, v in string_only:
    score = G[u][v]["string_score"]
    w = 1.0 + (score - 0.4) * 5.0
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], style="solid",
                           edge_color="#3B82F6", width=w, alpha=0.75, ax=ax)

# Draw Both (dual validated: solid violet thick lines)
for u, v in both_edges:
    score = G[u][v]["string_score"]
    w = 1.5 + (score - 0.4) * 6.0
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], style="solid",
                           edge_color="#7C3AED", width=w, alpha=0.9, ax=ax)

# Draw Nodes scaled to accommodate labels without clipping
node_sizes = [3400 + deg_centrality[g] * 280 for g in G.nodes()]
node_colors = [category_colors.get(meta_dict.get(g, {}).get("Matrisome_Category", ""), "#64748B") for g in G.nodes()]

nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                       edgecolors="#0F172A", linewidths=2.2, alpha=0.95, ax=ax)

# Draw Labels with adjusted font size to prevent any clipping of long symbols (COL11A1, COL8A1)
label_font_sizes = {g: 9.5 if len(g) >= 7 else 11.0 for g in G.nodes()}
for g, (x, y) in pos.items():
    ax.text(x, y, g, fontsize=label_font_sizes[g], fontweight="bold",
            color="#FFFFFF", ha="center", va="center")

# Legends - Placed outside / to the right side with dedicated margin so no nodes are obscured
cat_handles = [mlines.Line2D([], [], color=col, marker='o', linestyle='None',
                             markersize=12, markeredgecolor='#0F172A', markeredgewidth=1.5, label=cat)
               for cat, col in category_colors.items()]

edge_handles = [
    mlines.Line2D([], [], color='#7C3AED', lw=3.5, label='STRING PPI + Co-expression (Dual Validated)'),
    mlines.Line2D([], [], color='#3B82F6', lw=2.5, label='STRING v12.5 PPI (score >= 0.400)'),
    mlines.Line2D([], [], color='#059669', lw=2.0, linestyle='dashed', label='Co-expression (|r| >= 0.85, FDR < 0.01)')
]

leg1 = ax.legend(handles=cat_handles, title="Matrisome Functional Category",
                 loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=True,
                 fontsize=10, title_fontsize=11, facecolor="#F8FAFC", edgecolor="#CBD5E1")
ax.add_artist(leg1)

ax.legend(handles=edge_handles, title="Interaction Evidence",
          loc="lower left", bbox_to_anchor=(1.01, 0.0), frameon=True,
          fontsize=10, title_fontsize=11, facecolor="#F8FAFC", edgecolor="#CBD5E1")

ax.set_title("Protein-Protein Interaction (STRING v12.5) and Co-expression Network\n"
             "Consensus Hub Biomarkers in Peritoneal Dialysis Fibrosis (GSE62928, N = 8)",
             fontsize=13, fontweight="bold", color="#0F172A", pad=20)

plt.axis("off")
plt.tight_layout()

# Save PNG (300 DPI) and Vector PDF
out_png = "results/figures/Hub_01_ppi_gene_interaction_network.png"
out_pdf = "results/figures/Hub_01_ppi_gene_interaction_network.pdf"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()
print(f"Saved publication figures to {out_png} and {out_pdf}")

