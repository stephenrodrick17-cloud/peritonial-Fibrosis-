"""
scripts/string_legend.py
Legend for the STRING PPI network figure (node colour, node size = degree, edge width = score).
Used by scripts/07_module_d_string_ppi.py and scripts/31_replot_string_with_legend.py.
"""
from matplotlib.lines import Line2D


def add_string_legend(ax, top_node, top_wdeg, node_size_fn, edge_width_fn,
                      degrees=(1, 4, 7), scores=(0.4, 0.7, 0.95)):
    """node_size_fn(degree) -> marker area (pt^2) used by nx.draw_networkx_nodes;
    edge_width_fn(score) -> linewidth used by nx.draw_networkx_edges."""
    colour = [
        Line2D([], [], ls="", marker="o", ms=12, mfc="#ef4444", mec="#1e293b",
               label=f"{top_node} (max weighted degree, {top_wdeg:.2f})"),
        Line2D([], [], ls="", marker="o", ms=12, mfc="#3b82f6", mec="#1e293b", label="Connected hub gene"),
        Line2D([], [], ls="", marker="o", ms=12, mfc="#94a3b8", mec="#1e293b", label="Isolated hub gene (degree 0)"),
    ]
    # scatter area s (pt^2) -> marker diameter ms (pt) = sqrt(s)
    size = [Line2D([], [], ls="", marker="o", ms=node_size_fn(d) ** 0.5 * 0.5, mfc="white", mec="#1e293b",
                   label=f"degree {d}") for d in degrees]
    width = [Line2D([], [], color="#64748b", alpha=0.6, lw=edge_width_fn(s),
                    label=f"score {s:.2f}") for s in scores]
    leg1 = ax.legend(handles=colour, title="Node colour", loc="upper left", bbox_to_anchor=(1.0, 1.0),
                     frameon=False, fontsize=9, title_fontsize=10)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=size, title="Node size = degree\n(markers drawn at 0.5x)", loc="upper left",
                     bbox_to_anchor=(1.0, 0.78), frameon=False, fontsize=9, title_fontsize=10,
                     labelspacing=2.6, borderpad=1.2, handletextpad=1.5)
    ax.add_artist(leg2)
    leg3 = ax.legend(handles=width, title="Edge width =\nSTRING combined score", loc="upper left",
                     bbox_to_anchor=(1.0, 0.30), frameon=False, fontsize=9, title_fontsize=10)
    return [leg1, leg2, leg3]
