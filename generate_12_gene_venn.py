import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

# 1. Load Data
df_ecm_deg = pd.read_csv('convergent_ECM_DEGs_nominal.csv')
df_mr_sig = pd.read_csv('mr_significant_genes.csv')

# Extract unique gene symbols
genes_ecm_deg = set(df_ecm_deg['Gene Symbol'].dropna().str.strip().str.upper().unique())
mr_genes = set(df_mr_sig['GeneSymbol'].dropna().str.strip().str.upper().unique())

# Compute overlap
overlap_genes = sorted(list(genes_ecm_deg.intersection(mr_genes)))

n_ecm = len(genes_ecm_deg)    # 148 unique genes
n_mr = len(mr_genes)          # 2636 unique genes
n_overlap = len(overlap_genes) # 12 unique genes

df_deg_info = df_ecm_deg.drop_duplicates(subset=['Gene Symbol'])
up_genes = [g for g in overlap_genes if df_deg_info[df_deg_info['Gene Symbol'] == g]['Direction'].values[0] == 'UP']
down_genes = [g for g in overlap_genes if df_deg_info[df_deg_info['Gene Symbol'] == g]['Direction'].values[0] == 'DOWN']

# 2. Create High-Resolution Publication Figure
fig, ax = plt.subplots(figsize=(11, 8.5), facecolor="#F8FAFC")

# Proportional Venn diagram
v = venn2(
    subsets=(n_ecm - n_overlap, n_mr - n_overlap, n_overlap),
    set_labels=(f'Peritoneal Convergent ECM-DEGs\n(P < 0.05 & |log₂FC| ≥ 0.585)\n[N = {n_ecm}]', 
                f'Transcriptome-Wide MR Screen\n(Nominal Causal P < 0.05)\n[N = {n_mr:,}]'),
    set_colors=('#6366F1', '#10B981'),
    alpha=0.65,
    ax=ax
)

# Custom circle borders
circles = venn2_circles(
    subsets=(n_ecm - n_overlap, n_mr - n_overlap, n_overlap),
    linestyle='solid',
    linewidth=2.5,
    color='#1E293B',
    ax=ax
)

# Format label text
for text in v.set_labels:
    if text:
        text.set_fontsize(12)
        text.set_fontweight('bold')
        text.set_color('#0F172A')

if v.get_label_by_id('10'):
    v.get_label_by_id('10').set_text(str(n_ecm - n_overlap))
    v.get_label_by_id('10').set_fontsize(15)
    v.get_label_by_id('10').set_fontweight('bold')
    v.get_label_by_id('10').set_color('#1E1B4B')

if v.get_label_by_id('01'):
    v.get_label_by_id('01').set_text(f"{n_mr - n_overlap:,}")
    v.get_label_by_id('01').set_fontsize(15)
    v.get_label_by_id('01').set_fontweight('bold')
    v.get_label_by_id('01').set_color('#064E3B')

if v.get_label_by_id('11'):
    v.get_label_by_id('11').set_text(str(n_overlap))
    v.get_label_by_id('11').set_fontsize(18)
    v.get_label_by_id('11').set_fontweight('bold')
    v.get_label_by_id('11').set_color('#B91C1C')

# Title
ax.set_title("Causal Convergence: Peritoneal Fibrosis ECM-DEGs ∩ Mendelian Randomization\n(Cross-Cohort Transcriptomics ∩ eQTLGen-CKDGen Two-Sample MR)", 
             fontsize=14, fontweight='bold', color='#0F172A', pad=25)

# Annotation Box for the 12 Causal Genes
box_text = (
    "★ 12 Convergent Causal ECM Genes Identified:\n"
    "• Up-regulated Drivers (n=4):  ADAM28, ADAMTS1, P4HA2 (FDR=3.6e-7), TNC\n"
    "• Down-regulated Matrix Regulators (n=8):  AMH, BMP6 (FDR=0.013), COL4A2, CRHBP, FGL2, IGFBP3 (FDR=0.023), LTBP4, WNT11\n\n"
    "Criteria: Peritoneal DEG (P < 0.05, |log₂FC| ≥ 0.585) ∩ Human Matrisome ∩ eQTLGen-CKDGen MR (P < 0.05)"
)

ax.text(
    0.5, -0.16, box_text,
    transform=ax.transAxes,
    fontsize=10.5,
    fontweight='medium',
    ha='center',
    va='top',
    color='#0F172A',
    bbox=dict(boxstyle='round,pad=0.8', facecolor='#EFF6FF', edgecolor='#93C5FD', linewidth=1.5)
)

plt.tight_layout()

# Save figures
os.makedirs("transcriptome_wide_mr_results", exist_ok=True)
plt.savefig("venn_mr_deg_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.savefig("transcriptome_wide_mr_results/venn_mr_deg_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.savefig("transcriptome_wide_mr_results/venn_12_causal_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())

print(f"[SUCCESS] 12-gene Venn diagram successfully generated!")
