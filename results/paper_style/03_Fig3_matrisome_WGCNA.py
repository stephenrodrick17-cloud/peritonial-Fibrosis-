"""
Figure 3: Screening of Matrisome-Related Genes in Peritoneal Dialysis Effluent
Replicating the exact figure layout and visual style of Sun et al. (Electron J Biotechnol 2026, Fig. 3).
Panels:
  A: WGCNA scale independence (left) and mean connectivity (right) vs soft-thresholding power
  B: Gene dendrogram and dynamic tree cut module colors
  C: Module-trait relationships heatmap (PD duration)
  D: Three-way Venn diagram (DEGs, WGCNA non-grey modules, Matrisome)
  E: Functional enrichment analysis (GO BP, GO CC, GO MF, KEGG)
  F: Protein-protein interaction (PPI) network of matrisome-related genes
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from PIL import Image

# Output directory
OUT_DIR = os.path.join("results", "paper_style")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load Data
soft_df = pd.read_csv(os.path.join(OUT_DIR, "Fig3A_soft_threshold_source.csv"))
trait_df = pd.read_csv(os.path.join(OUT_DIR, "Fig3C_module_trait_source.csv"))
venn_counts = pd.read_csv(os.path.join(OUT_DIR, "Fig3D_three_way_intersection_counts.csv"))
enrich_bp = pd.read_csv(os.path.join("results", "revision", "wgcna_paper_figures", "GO_BP_enrichment_39_non_grey_intersection_FDR05.csv"))
enrich_all = pd.read_csv(os.path.join("results", "tables", "T03_enrichment_all_terms.csv"))

# Dendrogram and STRING PPI images
dendro_img_path = os.path.join("results", "revision", "wgcna_paper_figures", "Fig3B_gene_dendrogram_modules.png")
string_img_path = os.path.join("results", "figures", "STRING_PPI_44_ECM_DEGs.png")

# Figure Setup (14 x 16.5 inches)
fig = plt.figure(figsize=(14.0, 16.8), facecolor="#FFFFFF")
gs = gridspec.GridSpec(3, 2, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.15, 1.25], hspace=0.32, wspace=0.26)

# ---------------------------------------------------------
# PANEL A: WGCNA Soft Threshold Power (Scale Free R^2 & Mean Connectivity)
# ---------------------------------------------------------
gs_a = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0, 0], wspace=0.28)
ax_a1 = fig.add_subplot(gs_a[0, 0])
ax_a2 = fig.add_subplot(gs_a[0, 1])

# A1: Scale independence
powers = soft_df["power"].values
r2_vals = soft_df["scale_free_fit_R2"].values
conn_vals = soft_df["mean_connectivity"].values

ax_a1.scatter(powers, r2_vals, color="#DC2626", s=18, zorder=3)
for p, r in zip(powers, r2_vals):
    ax_a1.text(p, r + 0.035, f"{int(p)}", color="#DC2626", fontsize=7.5, ha="center", va="bottom", fontweight="bold")

ax_a1.axhline(0.85, color="#DC2626", linestyle="-", lw=1.2)
ax_a1.set_xlabel("Soft Threshold (power)", fontsize=9.5)
ax_a1.set_ylabel("Scale Free Topology Model Fit, signed R^2", fontsize=9)
ax_a1.set_title("Scale independence", fontsize=10.5, fontweight="bold", pad=8)
ax_a1.set_ylim(0.0, 1.08)
ax_a1.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
ax_a1.text(-0.25, 1.12, "A", transform=ax_a1.transAxes, fontsize=17, fontweight="bold", va="top")

# A2: Mean connectivity
ax_a2.scatter(powers, conn_vals, color="#DC2626", s=18, zorder=3)
for p, c in zip(powers, conn_vals):
    ax_a2.text(p, c + 35, f"{int(p)}", color="#DC2626", fontsize=7.5, ha="center", va="bottom", fontweight="bold")

ax_a2.set_xlabel("Soft Threshold (power)", fontsize=9.5)
ax_a2.set_ylabel("Mean Connectivity", fontsize=9.5)
ax_a2.set_title("Mean connectivity", fontsize=10.5, fontweight="bold", pad=8)
ax_a2.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")

# ---------------------------------------------------------
# PANEL B: Gene Dendrogram and Module Colors
# ---------------------------------------------------------
ax_b = fig.add_subplot(gs[1, 0])
if os.path.exists(dendro_img_path):
    img_b = Image.open(dendro_img_path)
    ax_b.imshow(img_b)
    ax_b.axis("off")
else:
    ax_b.text(0.5, 0.5, "Gene Dendrogram & Module Colors", ha="center", va="center")
ax_b.set_title("Gene dendrogram and module colors", fontsize=11, fontweight="bold", pad=8)
ax_b.text(-0.10, 1.05, "B", transform=ax_b.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL C: Module-Trait Relationships
# ---------------------------------------------------------
ax_c = fig.add_subplot(gs[2, 0])

# Prepare module-trait matrix
mod_names = trait_df["Module"].tolist()
corrs = trait_df["Correlation"].values
pvals = trait_df["P_value"].values

# Sort by correlation
sort_idx = np.argsort(corrs)[::-1]
mod_sorted = [mod_names[i] for i in sort_idx]
corr_sorted = corrs[sort_idx]
pval_sorted = pvals[sort_idx]

# Plot heatmap of 1 column (PD duration) or 2 columns (LPD, SPD)
heat_matrix = np.column_stack([corr_sorted, -corr_sorted])
c_im = ax_c.imshow(heat_matrix, aspect="auto", cmap="coolwarm", vmin=-1.0, vmax=1.0)

ax_c.set_xticks([0, 1])
ax_c.set_xticklabels(["LPD Duration", "SPD Duration"], fontsize=9.5, fontweight="bold")
ax_c.set_yticks(range(len(mod_sorted)))
ax_c.set_yticklabels(mod_sorted, fontsize=8.5, fontweight="bold")

# Annotate cells with r and (P)
for i in range(len(mod_sorted)):
    # Col 0 (LPD)
    ax_c.text(0, i, f"{corr_sorted[i]:.2f}\n({pval_sorted[i]:.2g})",
              ha="center", va="center", fontsize=7.2, color="#0F172A" if abs(corr_sorted[i]) < 0.65 else "#FFFFFF")
    # Col 1 (SPD)
    ax_c.text(1, i, f"{-corr_sorted[i]:.2f}\n({pval_sorted[i]:.2g})",
              ha="center", va="center", fontsize=7.2, color="#0F172A" if abs(corr_sorted[i]) < 0.65 else "#FFFFFF")

ax_c.set_title("Module-trait relationships", fontsize=11, fontweight="bold", pad=8)
cbar = plt.colorbar(c_im, ax=ax_c, fraction=0.046, pad=0.04)
cbar.set_label("Correlation", fontsize=8.5)
ax_c.text(-0.16, 1.05, "C", transform=ax_c.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL D: Three-way Venn Diagram (DEGs, WGCNA, Matrisome)
# ---------------------------------------------------------
# Placed on upper right or middle
gs_d = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0, 1])
ax_d = fig.add_subplot(gs_d[0])
ax_d.set_xlim(0, 10)
ax_d.set_ylim(0, 10)
ax_d.axis("off")

# Circles: DEGs (top-left, pink), WGCNA (top-right, cyan), Matrisome (bottom, gold)
c_deg = patches.Circle((3.8, 6.2), 2.7, facecolor="#F472B6", edgecolor="#DB2777", alpha=0.45, lw=1.6)
c_wgcna = patches.Circle((6.2, 6.2), 2.7, facecolor="#38BDF8", edgecolor="#0284C7", alpha=0.45, lw=1.6)
c_mat = patches.Circle((5.0, 4.0), 2.7, facecolor="#FBBF24", edgecolor="#D97706", alpha=0.45, lw=1.6)

ax_d.add_patch(c_deg)
ax_d.add_patch(c_wgcna)
ax_d.add_patch(c_mat)

# Labels
ax_d.text(2.2, 9.2, "DEGs (813)", fontsize=10, fontweight="bold", color="#BE185D", ha="center")
ax_d.text(7.8, 9.2, "WGCNA (4,549)", fontsize=10, fontweight="bold", color="#0369A1", ha="center")
ax_d.text(5.0, 0.8, "Matrisome (428)", fontsize=10, fontweight="bold", color="#B45309", ha="center")

# Overlap counts
ax_d.text(2.6, 6.6, "151", fontsize=9.5, fontweight="bold", color="#831843", ha="center") # DEG alone
ax_d.text(7.4, 6.6, "3,892", fontsize=9.5, fontweight="bold", color="#0C4A6E", ha="center") # WGCNA alone
ax_d.text(5.0, 2.5, "345", fontsize=9.5, fontweight="bold", color="#78350F", ha="center") # Mat alone
ax_d.text(5.0, 7.3, "618", fontsize=9.5, fontweight="bold", color="#1E293B", ha="center") # DEG ∩ WGCNA
ax_d.text(3.7, 4.5, "6", fontsize=9.5, fontweight="bold", color="#1E293B", ha="center") # DEG ∩ Mat
ax_d.text(6.3, 4.5, "39", fontsize=9.5, fontweight="bold", color="#1E293B", ha="center") # WGCNA ∩ Mat
ax_d.text(5.0, 5.2, "38", fontsize=11, fontweight="bold", color="#DC2626", ha="center") # Triple

ax_d.text(5.0, 9.8, "Venn Diagram of Candidate Genes", fontsize=11, fontweight="bold", ha="center")
ax_d.text(-0.06, 1.10, "D", transform=ax_d.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL E: Functional Enrichment (GO BP, CC, MF, KEGG in 2x2 grid)
# ---------------------------------------------------------
gs_e = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1, 1], hspace=0.38, wspace=0.38)

# E1: GO BP
ax_e1 = fig.add_subplot(gs_e[0, 0])
top_bp = enrich_bp.head(5).iloc[::-1]
b_bp = ax_e1.barh(range(len(top_bp)), top_bp["Count"], color="#DC2626", height=0.65)
ax_e1.set_yticks(range(len(top_bp)))
ax_e1.set_yticklabels([d[:22] + "..." if len(d) > 22 else d for d in top_bp["Description"]], fontsize=7.2)
ax_e1.set_xlabel("Count", fontsize=8)
ax_e1.set_title("GO enrichment analysis (BP)", fontsize=8.5, fontweight="bold")
ax_e1.grid(True, linestyle=":", alpha=0.4, axis="x")

# E2: GO CC
ax_e2 = fig.add_subplot(gs_e[0, 1])
cc_df = enrich_all[enrich_all["Database"] == "GO_CC"].sort_values(by="pvalue").head(5).iloc[::-1]
ax_e2.barh(range(len(cc_df)), cc_df["Count"], color="#DC2626", height=0.65)
ax_e2.set_yticks(range(len(cc_df)))
ax_e2.set_yticklabels([d[:22] + "..." if len(d) > 22 else d for d in cc_df["Description"]], fontsize=7.2)
ax_e2.set_xlabel("Count", fontsize=8)
ax_e2.set_title("GO enrichment analysis (CC)", fontsize=8.5, fontweight="bold")
ax_e2.grid(True, linestyle=":", alpha=0.4, axis="x")

# E3: GO MF
ax_e3 = fig.add_subplot(gs_e[1, 0])
mf_df = enrich_all[enrich_all["Database"] == "GO_MF"].sort_values(by="pvalue").head(5).iloc[::-1]
ax_e3.barh(range(len(mf_df)), mf_df["Count"], color="#DC2626", height=0.65)
ax_e3.set_yticks(range(len(mf_df)))
ax_e3.set_yticklabels([d[:22] + "..." if len(d) > 22 else d for d in mf_df["Description"]], fontsize=7.2)
ax_e3.set_xlabel("Count", fontsize=8)
ax_e3.set_title("GO enrichment analysis (MF)", fontsize=8.5, fontweight="bold")
ax_e3.grid(True, linestyle=":", alpha=0.4, axis="x")

# E4: KEGG
ax_e4 = fig.add_subplot(gs_e[1, 1])
kegg_sub = enrich_all[enrich_all["Database"] == "KEGG"].sort_values(by="pvalue").head(5).iloc[::-1]
ax_e4.barh(range(len(kegg_sub)), kegg_sub["Count"], color="#DC2626", height=0.65)
ax_e4.set_yticks(range(len(kegg_sub)))
ax_e4.set_yticklabels([d[:22] + "..." if len(d) > 22 else d for d in kegg_sub["Description"]], fontsize=7.2)
ax_e4.set_xlabel("Count", fontsize=8)
ax_e4.set_title("KEGG enrichment analysis", fontsize=8.5, fontweight="bold")
ax_e4.grid(True, linestyle=":", alpha=0.4, axis="x")

ax_e1.text(-0.35, 1.25, "E", transform=ax_e1.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL F: STRING PPI Network of Matrisome DEGs
# ---------------------------------------------------------
ax_f = fig.add_subplot(gs[2, 1])
if os.path.exists(string_img_path):
    img_f = Image.open(string_img_path)
    ax_f.imshow(img_f)
    ax_f.axis("off")
else:
    ax_f.text(0.5, 0.5, "STRING PPI Network", ha="center", va="center")
ax_f.set_title("PPI network of 44 matrisome-related DEGs", fontsize=11, fontweight="bold", pad=8)
ax_f.text(-0.08, 1.05, "F", transform=ax_f.transAxes, fontsize=17, fontweight="bold", va="top")

# Caption at the very bottom
fig.text(0.5, 0.012,
         "Fig. 3. Screening of matrisome-related genes in peritoneal dialysis effluent. (A) Selection of soft-threshold power and mean connectivity by WGCNA.\n"
         "(B) Gene dendrogram and co-expression module division by dynamic tree cut. (C) Module-trait relationships with PD duration. (D) Venn diagram of\n"
         "matrisome-related genes intersecting DEGs, WGCNA non-grey modules, and Naba matrisome. (E) Functional enrichment analysis (GO BP, CC, MF, and KEGG).\n"
         "(F) STRING protein-protein interaction (PPI) network of 44 matrisome-related DEGs.",
         ha="center", va="bottom", fontsize=9.2, color="#0F172A", linespacing=1.35)

plt.subplots_adjust(bottom=0.07, top=0.96)

# Save high-res PNG and PDF
png_path = os.path.join(OUT_DIR, "Fig3_matrisome_WGCNA.png")
pdf_path = os.path.join(OUT_DIR, "Fig3_matrisome_WGCNA.pdf")
fig_png = os.path.join(OUT_DIR, "Fig3.png")
fig_pdf = os.path.join(OUT_DIR, "Fig3.pdf")

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
plt.savefig(fig_png, dpi=300, bbox_inches="tight")
plt.savefig(fig_pdf, bbox_inches="tight")
plt.close()

print(f"Successfully generated Figure 3: {png_path} and {pdf_path}")
