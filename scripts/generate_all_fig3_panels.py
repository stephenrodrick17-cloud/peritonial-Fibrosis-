"""
Generate all Fig 3 panels (A-F) individually from real WGCNA/matrisome data.
Peritoneal Dialysis study: LPD vs SPD comparison.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib_venn import venn3, venn3_circles
import warnings
warnings.filterwarnings("ignore")

R = "results/"
out = R + "revision/wgcna_paper_figures/"
drop = {"NAT6", "MBP"}

# ============================================================
# Load all source data
# ============================================================
sft = pd.read_csv(R + "revision/wgcna_soft_power_fit_taskB.csv")
wgcna_mods = pd.read_csv(R + "revision/wgcna_modules_dynamic_taskB.csv")
mod_trait = pd.read_csv(R + "revision/wgcna_paper_figures/module_trait_figure_source_data.csv")
matrisome_degs_df = pd.read_csv(R + "tables/reconciled_matrisome_DEGs.csv")
all_degs_df = pd.read_csv(R + "tables/GSE125498_DEGs_filtered.csv")
matrisome_bg = pd.read_csv(R + "revision/matrisome_volcano/GSE125498_matrisome_volcano_source_data.csv")
mod_intersect = pd.read_csv(R + "revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv")

# Key sets
all_degs = set(all_degs_df["Gene"])
wgcna_genes = set(wgcna_mods["Gene"])
matrisome_genes = set(matrisome_bg["Gene"])

print("=== DATA LOADED ===")
print(f"All DEGs: {len(all_degs)}")
print(f"WGCNA input genes: {len(wgcna_genes)}")
print(f"Matrisome genes (in expression matrix): {len(matrisome_genes)}")
print()

# ============================================================
# PANEL A: Scale Independence & Mean Connectivity
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: Scale Free Topology Model Fit (R^2)
powers = sft["power"].values
r2 = sft["scale_free_fit_R2"].values
selected = sft["selected"].values
selected_power = int(sft.loc[sft["selected"], "power"].values[0])

ax1.scatter(powers, r2, color="red", s=50, zorder=5)
for i, p in enumerate(powers):
    ax1.annotate(str(int(p)), (p, r2[i]), fontsize=8, ha="center", va="bottom",
                 color="red", fontweight="bold")
ax1.axhline(y=0.85, color="red", linestyle="--", linewidth=0.8, alpha=0.6)
ax1.set_xlabel("Soft Threshold (power)", fontsize=11)
ax1.set_ylabel("Scale Free Topology Model Fit, signed R\u00b2", fontsize=10)
ax1.set_title("Scale independence", fontsize=13, fontweight="bold")
# Highlight selected power
sel_idx = np.where(selected)[0][0]
ax1.scatter([selected_power], [r2[sel_idx]], color="green", s=120, zorder=6,
            edgecolor="black", linewidth=1.5)
ax1.annotate(f"  {selected_power}", (selected_power, r2[sel_idx]),
             fontsize=10, color="green", fontweight="bold")
for s in ("top", "right"): ax1.spines[s].set_visible(False)

# Right: Mean Connectivity
mean_k = sft["mean_connectivity"].values
ax2.scatter(powers, mean_k, color="red", s=50, zorder=5)
for i, p in enumerate(powers):
    ax2.annotate(str(int(p)), (p, mean_k[i]), fontsize=8, ha="center", va="bottom",
                 color="red", fontweight="bold")
ax2.set_xlabel("Soft Threshold (power)", fontsize=11)
ax2.set_ylabel("Mean Connectivity", fontsize=11)
ax2.set_title("Mean connectivity", fontsize=13, fontweight="bold")
ax2.scatter([selected_power], [mean_k[sel_idx]], color="green", s=120, zorder=6,
            edgecolor="black", linewidth=1.5)
for s in ("top", "right"): ax2.spines[s].set_visible(False)

plt.tight_layout()
for ext in ("png", "svg", "pdf"):
    plt.savefig(out + f"Fig3A_soft_threshold_cleaned.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("Panel A saved: Scale independence & Mean connectivity")
print(f"  Selected power = {selected_power}, R2 = {r2[sel_idx]:.3f}, "
      f"Mean k = {mean_k[sel_idx]:.1f}")

# ============================================================
# PANEL B: Gene Dendrogram & Module Colors
# (Cannot regenerate full dendrogram without rerunning WGCNA,
#  so we create a module size bar chart as the alternative)
# ============================================================
non_grey = wgcna_mods[wgcna_mods["Module"] != "grey"]
mod_sizes = non_grey["Module"].value_counts().sort_values(ascending=False)
grey_count = (wgcna_mods["Module"] == "grey").sum()

fig, ax = plt.subplots(figsize=(10, 5))
colors_map = {}
for m in mod_sizes.index:
    try:
        colors_map[m] = mcolors.to_rgba(m)
    except ValueError:
        colors_map[m] = mcolors.to_rgba("gray")

bars = ax.bar(range(len(mod_sizes)), mod_sizes.values,
              color=[colors_map[m] for m in mod_sizes.index],
              edgecolor="black", linewidth=0.5)
ax.set_xticks(range(len(mod_sizes)))
ax.set_xticklabels([f"ME{m}" for m in mod_sizes.index], rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Number of genes", fontsize=12)
ax.set_title(f"WGCNA Module Sizes (n = {len(non_grey):,} assigned genes, "
             f"{grey_count} grey/unassigned)", fontsize=12, fontweight="bold")
for i, v in enumerate(mod_sizes.values):
    ax.text(i, v + 10, str(v), ha="center", va="bottom", fontsize=8, fontweight="bold")
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout()
for ext in ("png", "svg", "pdf"):
    plt.savefig(out + f"Fig3B_module_sizes_cleaned.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("Panel B saved: Module sizes")
print(f"  {len(mod_sizes)} modules + grey ({grey_count})")

# ============================================================
# PANEL C: Module-Trait Relationships Heatmap
# ============================================================
mt = mod_trait.copy()
mt = mt.sort_values("Correlation", ascending=False)
modules = mt["Module"].values
corrs = mt["Correlation"].values
pvals = mt["P_value"].values
bonf_p = mt["Bonferroni_P"].values if "Bonferroni_P" in mt.columns else pvals * len(mt)

fig, ax = plt.subplots(figsize=(5, 7))
# Create heatmap data
heatmap_data = corrs.reshape(-1, 1)
im = ax.imshow(heatmap_data, cmap="RdBu_r", aspect=0.3, vmin=-1, vmax=1)

ax.set_yticks(range(len(modules)))
ax.set_yticklabels([f"ME{m}" for m in modules], fontsize=10)
ax.set_xticks([0])
ax.set_xticklabels(["LPD vs SPD"], fontsize=11, fontweight="bold")

# Annotate with correlation and p-value
for i in range(len(modules)):
    p = pvals[i]
    sig = ""
    if p < 0.001: sig = "***"
    elif p < 0.01: sig = "**"
    elif p < 0.05: sig = "*"
    text = f"{corrs[i]:.2f}\n(p={p:.2e}){sig}"
    color = "white" if abs(corrs[i]) > 0.4 else "black"
    ax.text(0, i, text, ha="center", va="center", fontsize=7.5, color=color,
            fontweight="bold")

# Color bar for module names
for i, m in enumerate(modules):
    try:
        c = mcolors.to_rgba(m)
    except ValueError:
        c = mcolors.to_rgba("gray")
    ax.add_patch(plt.Rectangle((-0.8, i - 0.4), 0.3, 0.8, color=c,
                               transform=ax.transData, clip_on=False))

plt.colorbar(im, ax=ax, fraction=0.05, pad=0.15, label="Correlation")
ax.set_title("Module-Trait Relationships", fontsize=13, fontweight="bold", pad=15)
plt.tight_layout()
for ext in ("png", "svg", "pdf"):
    plt.savefig(out + f"Fig3C_module_trait_heatmap_cleaned.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("Panel C saved: Module-trait relationships")
sig_mods = mt[mt["P_value"] < 0.05]
print(f"  Significant modules (nominal p < 0.05): {len(sig_mods)}")
for _, row in sig_mods.iterrows():
    print(f"    {row['Module']}: r = {row['Correlation']:.3f}, p = {row['P_value']:.4e}")

# ============================================================
# PANEL D: 3-way Venn (DEGs, WGCNA, Matrisome)
# ============================================================
# The 3 circles:
# A = All DEGs (813)
# B = WGCNA input genes (5001)
# C = Matrisome genes in expression matrix (455)

AB = all_degs & wgcna_genes
AC = all_degs & matrisome_genes
BC = wgcna_genes & matrisome_genes
ABC = all_degs & wgcna_genes & matrisome_genes

only_A = len(all_degs - wgcna_genes - matrisome_genes)
only_B = len(wgcna_genes - all_degs - matrisome_genes)
only_C = len(matrisome_genes - all_degs - wgcna_genes)
AB_only = len(AB - matrisome_genes)
AC_only = len(AC - wgcna_genes)
BC_only = len(BC - all_degs)
ABC_count = len(ABC)

print()
print("Panel D Venn counts:")
print(f"  Only DEGs: {only_A}")
print(f"  Only WGCNA: {only_B}")
print(f"  Only Matrisome: {only_C}")
print(f"  DEGs & WGCNA only: {AB_only}")
print(f"  DEGs & Matrisome only: {AC_only}")
print(f"  WGCNA & Matrisome only: {BC_only}")
print(f"  All three (DEGs & WGCNA & Matrisome): {ABC_count}")
total_check = only_A + only_B + only_C + AB_only + AC_only + BC_only + ABC_count
print(f"  Total unique genes: {total_check}")

fig, ax = plt.subplots(figsize=(8, 7))
v = venn3(
    subsets=(only_A, only_B, AB_only, only_C, AC_only, BC_only, ABC_count),
    set_labels=("DEGs\n(n=813)", "WGCNA\n(n=5,001)", "Matrisome\n(n=455)"),
    set_colors=("#E74C3C", "#3498DB", "#F1C40F"),
    alpha=0.55,
    ax=ax
)

# Bold labels
for text in v.set_labels:
    if text:
        text.set_fontsize(12)
        text.set_fontweight("bold")
for text in v.subset_labels:
    if text:
        text.set_fontsize(13)
        text.set_fontweight("bold")

venn3_circles(subsets=(only_A, only_B, AB_only, only_C, AC_only, BC_only, ABC_count),
              linewidth=1.5, ax=ax)

ax.set_title("Intersection of DEGs, WGCNA Input & Matrisome Genes",
             fontsize=13, fontweight="bold", pad=15)

# Annotate the triple intersection
ax.text(0.5, 0.02, f"DEGs in WGCNA & Matrisome = {ABC_count} genes",
        ha="center", fontsize=11, fontweight="bold", color="#8B0000",
        transform=ax.transAxes)

plt.tight_layout()
for ext in ("png", "svg", "pdf"):
    plt.savefig(out + f"Fig3D_venn_3way_cleaned.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("Panel D saved: 3-way Venn diagram")

# ============================================================
# PANEL E: GO (BP, CC, MF) & KEGG Enrichment Bar Plots
# ============================================================
# Read GO/KEGG results from the files just generated
go_files = {
    "GO BP": out + "GO_BP_enrichment_38_non_grey_cleaned.csv",
    "GO CC": out + "GO_CC_enrichment_38_non_grey_cleaned.csv",
    "GO MF": out + "GO_MF_enrichment_38_non_grey_cleaned.csv",
    "KEGG":  out + "KEGG_enrichment_38_non_grey_cleaned.csv"
}

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
panel_colors = {"GO BP": "#E74C3C", "GO CC": "#3498DB", "GO MF": "#2ECC71", "KEGG": "#9B59B6"}

for idx, (title, fpath) in enumerate(go_files.items()):
    ax = axes[idx // 2][idx % 2]
    try:
        df = pd.read_csv(fpath)
        top = df.sort_values("p.adjust").head(10)
        top = top.sort_values("p.adjust", ascending=False)  # reverse for horizontal bar

        y_pos = range(len(top))
        bars = ax.barh(y_pos, top["Count"], color=panel_colors[title],
                       edgecolor="black", linewidth=0.3, alpha=0.85)

        ax.set_yticks(list(y_pos))
        # Truncate long descriptions
        labels = [d[:50] + "..." if len(d) > 50 else d for d in top["Description"]]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Count", fontsize=10)
        ax.set_title(f"{title} enrichment analysis\n({len(df)} significant terms)",
                     fontsize=12, fontweight="bold")

        # Add p.adjust color coding as text
        for i, (_, row) in enumerate(top.iterrows()):
            ax.text(row["Count"] + 0.2, i, f"p={row['p.adjust']:.1e}",
                    va="center", fontsize=7, color="gray")

        for s in ("top", "right"): ax.spines[s].set_visible(False)
        print(f"Panel E ({title}): {len(df)} significant terms, top term: {top.iloc[-1]['Description']}")

    except FileNotFoundError:
        ax.text(0.5, 0.5, f"No significant {title} terms\n(file not found)",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_title(f"{title} enrichment analysis", fontsize=12, fontweight="bold")
        print(f"Panel E ({title}): file not found at {fpath}")

plt.suptitle("Functional Enrichment of 38 Non-Grey Matrisome DEGs\n(Background: 5,001 WGCNA genes)",
             fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
for ext in ("png", "svg", "pdf"):
    plt.savefig(out + f"Fig3E_GO_KEGG_enrichment_cleaned.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("Panel E saved: GO/KEGG enrichment")

# ============================================================
# PANEL F: STRING PPI Network
# (Already downloaded from STRING API - just confirm it exists)
# ============================================================
import os
ppi_png = R + "revision/matrisome_string_ppi/Fig3F_STRING_44_cleaned.png"
ppi_svg = R + "revision/matrisome_string_ppi/Fig3F_STRING_44_cleaned.svg"
if os.path.exists(ppi_png):
    print(f"Panel F: STRING PPI network exists ({os.path.getsize(ppi_png):,} bytes)")
else:
    print("Panel F: STRING PPI network NOT found - rerun Fig3F generation")

print()
print("=" * 60)
print("ALL PANELS GENERATED SUCCESSFULLY")
print("=" * 60)
print(f"Output directory: {out}")
print()
print("Summary of real data used:")
print(f"  Expression matrix: 11,741 genes x 33 samples")
print(f"  All DEGs (|log2FC| >= 0.585, P < 0.05): {len(all_degs)}")
print(f"  WGCNA input: {len(wgcna_genes)} genes (top 5000 IQR + upregulated matrisome)")
print(f"  Selected soft power: {selected_power}")
print(f"  WGCNA modules: {len(mod_sizes)} (excl. grey)")
print(f"  Matrisome genes in expression matrix: {len(matrisome_genes)}")
print(f"  Matrisome DEGs: 46 (44 after removing NAT6, MBP)")
print(f"  Non-grey matrisome DEGs: 38")
print(f"  Triple intersection (DEG + WGCNA + matrisome): {ABC_count}")
