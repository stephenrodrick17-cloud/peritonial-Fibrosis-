import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

R = "results/"

# Load data
mod = pd.read_csv(R + "revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv")
gcol = "Gene.symbol" if "Gene.symbol" in mod.columns else "Annotated Gene"
wgcna = pd.read_csv(R + "revision/wgcna_modules_dynamic_taskB.csv")

matrisome_degs = set(mod[gcol].dropna())
wgcna_genes = set(wgcna["Gene"].dropna())

overlap = matrisome_degs & wgcna_genes
only_matrisome = matrisome_degs - wgcna_genes
only_wgcna = wgcna_genes - matrisome_degs

print(f"Matrisome DEGs total: {len(matrisome_degs)}")
print(f"WGCNA input total: {len(wgcna_genes)}")
print(f"Overlap: {len(overlap)}")
print(f"Only in matrisome DEGs: {len(only_matrisome)} -> {sorted(only_matrisome)}")
print(f"Only in WGCNA: {len(only_wgcna)}")

# --- Venn diagram ---
fig, ax = plt.subplots(figsize=(8, 6))

v = venn2(
    subsets=(len(only_matrisome), len(only_wgcna), len(overlap)),
    set_labels=("Matrisome DEGs\n(n = 46)", "WGCNA Input Genes\n(n = 5,001)"),
    set_colors=("#2196F3", "#FF9800"),
    alpha=0.55,
    ax=ax
)

# Style the labels
for text in v.set_labels:
    text.set_fontsize(13)
    text.set_fontweight("bold")

for text in v.subset_labels:
    text.set_fontsize(16)
    text.set_fontweight("bold")

# Add circles outline
venn2_circles(subsets=(len(only_matrisome), len(only_wgcna), len(overlap)),
              linewidth=1.5, ax=ax)

# Annotate the 3 non-overlapping genes
ax.text(0.15, -0.38,
        "Not in WGCNA input:\n" + ", ".join(sorted(only_matrisome)),
        ha="center", va="top", fontsize=10,
        style="italic", color="#1565C0",
        transform=ax.transAxes)

ax.set_title("Matrisome DEGs ∩ WGCNA Input Genes", fontsize=15, fontweight="bold", pad=20)
plt.tight_layout()

out = R + "revision/wgcna_paper_figures/"
for ext in ("png", "svg", "pdf"):
    plt.savefig(out + f"venn_matrisome_wgcna_intersection.{ext}", dpi=300,
                bbox_inches="tight")
print("Saved.")
