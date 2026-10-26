"""
Plot GeneMANIA network matching the exact Cytoscape/GeneMANIA publication layout.
Features:
- Bold panel label "A" in upper left.
- Concentric circular arrangement:
  - Outer ring: 20 predicted interactors
  - Inner cluster: 7 diagnostic query hub genes
- Query hub genes: larger nodes with cross-hatched texture and functional slices
- Interactor genes: circular nodes with functional category pie slices
- Edges: multi-network colored lines with GeneMANIA palette and transparency
- Dual bottom legend: "Networks" on left, "Functions" on right with colored markers
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Output directories
out_dir = r"d:\Peritoneal Project\results\revision\hub_gene_network_drug_gsea"
art_dir = r"C:\Users\asus\.gemini\antigravity-ide\brain\a3bd40f9-59d0-450e-883c-c58b2c00e17e"
os.makedirs(out_dir, exist_ok=True)
os.makedirs(art_dir, exist_ok=True)

# Load data
nodes_df = pd.read_csv(os.path.join(out_dir, "genemania_network_nodes.csv"))
edges_df = pd.read_csv(os.path.join(out_dir, "genemania_network_edges.csv"))

with open(os.path.join(out_dir, "genemania_raw_response.json"), "r", encoding="utf-8") as f:
    raw_data = json.load(f)

# 7 Selected Enriched Functions from GeneMANIA ontology
function_palette = [
    ("cellular response to tumor necrosis factor", "#3B7EBD"),  # Blue
    ("positive regulation of lymphocyte activation", "#D39E36"),  # Ochre/Amber
    ("positive regulation of cell-cell adhesion", "#8C59A4"),  # Purple
    ("mononuclear cell proliferation", "#4E9B58"),  # Green
    ("regulation of endopeptidase activity", "#4FB5B5"),  # Teal
    ("NIK/NF-kappaB signaling", "#C24C44"),  # Coral Red
    ("cytokine receptor binding", "#B8558A"),  # Magenta/Rose
]

# Map category descriptions to ID
desc_to_color = {f[0]: f[1] for f in function_palette}
id_to_desc = {}
for c in raw_data.get("resultOntologyCategories", []):
    desc = c["ontologyCategory"]["description"]
    if desc in desc_to_color:
        id_to_desc[c["id"]] = desc

# Map genes to their active function colors
gene_func_colors = {}
for g in raw_data.get("resultGenes", []):
    sym = g["gene"]["symbol"]
    cat_ids = [cat["id"] for cat in g.get("resultOntologyCategories", [])]
    funcs = [id_to_desc[cid] for cid in cat_ids if cid in id_to_desc]
    unique_funcs = list(dict.fromkeys(funcs))
    gene_func_colors[sym] = [desc_to_color[fn] for fn in unique_funcs]

# Separate query hub genes and interactors
query_genes = nodes_df[nodes_df["query_gene"] == True]["gene"].tolist()
interactor_genes = nodes_df[nodes_df["query_gene"] == False]["gene"].tolist()

# Define positions
pos = {}
center_y = 0.66  # Shifted up to give plenty of room for legends
center_x = 0.50

r_outer = 0.35
r_inner = 0.15

# Place outer 20 interactors in a circle
n_outer = len(interactor_genes)
for i, g in enumerate(interactor_genes):
    # rotate clockwise
    angle = np.pi / 2 - 2 * np.pi * i / n_outer
    pos[g] = (center_x + r_outer * np.cos(angle), center_y + r_outer * np.sin(angle))

# Place inner query genes
n_inner = len(query_genes)
for j, g in enumerate(query_genes):
    angle = np.pi / 2 - 2 * np.pi * j / n_inner
    pos[g] = (center_x + r_inner * np.cos(angle), center_y + r_inner * np.sin(angle))

# Set up figure
fig = plt.figure(figsize=(9, 11), dpi=300, facecolor="white")
ax = fig.add_axes([0, 0, 1, 1], facecolor="white")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.15)
ax.axis("off")

# Add panel "A"
ax.text(0.06, 1.08, "A", fontsize=28, fontweight="bold", fontfamily="sans-serif", color="#222222")

# Edge colors
network_edge_styles = {
    "Co-expression": {"color": "#C288D9", "alpha": 0.45, "lw": 1.1},
    "Shared protein domains": {"color": "#D5C67A", "alpha": 0.38, "lw": 1.2},
    "Physical Interactions": {"color": "#E78282", "alpha": 0.45, "lw": 1.1},
    "Predicted": {"color": "#F3B05C", "alpha": 0.45, "lw": 1.1},
    "Co-localization": {"color": "#7BAAF7", "alpha": 0.45, "lw": 1.1},
    "Genetic Interactions": {"color": "#88CF88", "alpha": 0.45, "lw": 1.1},
    "Pathway": {"color": "#64D2D2", "alpha": 0.45, "lw": 1.1},
}

# Draw edges
for _, row in edges_df.iterrows():
    u = row["gene_a"]
    v = row["gene_b"]
    itype = row["interaction_type"]
    if u in pos and v in pos:
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        style = network_edge_styles.get(itype, {"color": "#C288D9", "alpha": 0.4, "lw": 1.0})
        ax.plot([x1, x2], [y1, y2], color=style["color"], alpha=style["alpha"], lw=style["lw"], zorder=2)

# Node radii
outer_node_radius = 0.035
inner_node_radius = 0.048

# Draw nodes exactly matching GeneMANIA Cytoscape
def draw_node(gene_name, x, y, radius, is_query):
    colors = gene_func_colors.get(gene_name, [])
    base_gray = "#484848"
    
    # 1. Base dark circle
    base_circle = patches.Circle((x, y), radius, facecolor=base_gray, edgecolor="#2A2A2A", lw=0.8, zorder=4)
    ax.add_patch(base_circle)
    
    # 2. If gene has functional categories, draw wedges
    if colors:
        n_slices = len(colors)
        slice_angle = 360.0 / n_slices
        for s_idx, col in enumerate(colors):
            theta1 = 90.0 + s_idx * slice_angle
            theta2 = theta1 + slice_angle
            wedge = patches.Wedge((x, y), radius, theta1, theta2, facecolor=col, edgecolor="#2A2A2A", lw=0.5, zorder=5)
            ax.add_patch(wedge)
            
    # 3. For query genes: GeneMANIA adds distinct diagonal cross-hatching
    if is_query:
        hatch_circle = patches.Circle((x, y), radius, facecolor="none", edgecolor="#DDDDDD", hatch="////", lw=0, zorder=6, alpha=0.55)
        ax.add_patch(hatch_circle)
        border = patches.Circle((x, y), radius, facecolor="none", edgecolor="#222222", lw=1.6, zorder=7)
        ax.add_patch(border)
    else:
        border = patches.Circle((x, y), radius, facecolor="none", edgecolor="#222222", lw=1.0, zorder=7)
        ax.add_patch(border)
        
    # 4. Gene label in the center
    # In Cytoscape/GeneMANIA, text has crisp readability inside the node
    font_size = 7.5 if is_query else 6.2
    # Slight contrast halo
    ax.text(x, y, gene_name, color="#FFFFFF", fontsize=font_size, fontweight="bold",
            ha="center", va="center", zorder=8, fontfamily="sans-serif",
            bbox=dict(boxstyle="round,pad=0.10", facecolor="#181818", edgecolor="none", alpha=0.72))

# Draw outer nodes
for g in interactor_genes:
    x, y = pos[g]
    draw_node(g, x, y, outer_node_radius, is_query=False)

# Draw inner query hub nodes
for g in query_genes:
    x, y = pos[g]
    draw_node(g, x, y, inner_node_radius, is_query=True)

# ----------------- LEGENDS AT BOTTOM -----------------
legend_y_start = 0.19
legend_line_h = 0.022

# Left column: Networks
net_x = 0.08
ax.text(net_x, legend_y_start, "Networks", fontsize=13, fontweight="bold", fontfamily="sans-serif", color="#111111")

all_network_items = [
    ("Physical Interactions", "#E78282"),
    ("Co-expression", "#C288D9"),
    ("Predicted", "#F3B05C"),
    ("Co-localization", "#7BAAF7"),
    ("Genetic Interactions", "#88CF88"),
    ("Pathway", "#64D2D2"),
    ("Shared protein domains", "#D5C67A"),
]

for idx, (net_name, net_col) in enumerate(all_network_items):
    y_pos = legend_y_start - (idx + 1) * legend_line_h
    sq = patches.Rectangle((net_x, y_pos - 0.0035), 0.015, 0.009, facecolor=net_col, edgecolor="none")
    ax.add_patch(sq)
    ax.text(net_x + 0.026, y_pos, net_name, fontsize=9.5, fontfamily="sans-serif", color="#333333", va="center")

# Right column: Functions
func_x = 0.50
ax.text(func_x, legend_y_start, "Functions", fontsize=13, fontweight="bold", fontfamily="sans-serif", color="#111111")

for idx, (func_name, func_col) in enumerate(function_palette):
    y_pos = legend_y_start - (idx + 1) * legend_line_h
    sq = patches.Rectangle((func_x, y_pos - 0.0035), 0.015, 0.009, facecolor=func_col, edgecolor="none")
    ax.add_patch(sq)
    ax.text(func_x + 0.026, y_pos, func_name, fontsize=9.5, fontfamily="sans-serif", color="#333333", va="center")

# Save outputs to both results dir and artifacts dir
png_path = os.path.join(out_dir, "genemania_network.png")
pdf_path = os.path.join(out_dir, "genemania_network.pdf")
art_png = os.path.join(art_dir, "Figure_A_GeneMANIA_Network.png")
art_pdf = os.path.join(art_dir, "Figure_A_GeneMANIA_Network.pdf")

for p in [png_path, art_png]:
    plt.savefig(p, dpi=300, bbox_inches="tight")
for p in [pdf_path, art_pdf]:
    plt.savefig(p, dpi=300, bbox_inches="tight")

plt.close()
print("GeneMANIA publication figure updated and synced!")
