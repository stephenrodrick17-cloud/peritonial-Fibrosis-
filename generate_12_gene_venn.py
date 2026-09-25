import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

# 1. Setup Proportional Venn Diagram
fig, ax = plt.subplots(figsize=(11, 8.5), facecolor="#F8FAFC")

# Subsets for proportional sizing
total_ecm_deg = 86
total_mr = 2636
overlap = 12

v = venn2(
    subsets=(total_ecm_deg - overlap, total_mr - overlap, overlap),
    set_labels=(f'Peritoneal Convergent ECM-DEGs\n(P < 0.05 & |log₂FC| ≥ 0.585)\n[N = {total_ecm_deg}]', 
                f'Transcriptome-Wide MR Screen\n(Nominal Causal P < 0.05)\n[N = {total_mr:,}]'),
    set_colors=('#6366F1', '#10B981'),
    alpha=0.65,
    ax=ax
)

# Custom circle borders
circles = venn2_circles(
    subsets=(total_ecm_deg - overlap, total_mr - overlap, overlap),
    linestyle='solid',
    linewidth=2.2,
    color='#1E293B',
    ax=ax
)

# Format label text
for text in v.set_labels:
    if text:
        text.set_fontsize(12)
        text.set_fontweight('bold')
        text.set_color('#0F172A')

# Explicitly set the number inside the left circle to 86
if v.get_label_by_id('10'):
    v.get_label_by_id('10').set_text('86')
    v.get_label_by_id('10').set_fontsize(16)
    v.get_label_by_id('10').set_fontweight('bold')
    v.get_label_by_id('10').set_color('#0F172A')

if v.get_label_by_id('01'):
    v.get_label_by_id('01').set_text(f"{total_mr:,}")
    v.get_label_by_id('01').set_fontsize(16)
    v.get_label_by_id('01').set_fontweight('bold')
    v.get_label_by_id('01').set_color('#0F172A')

if v.get_label_by_id('11'):
    v.get_label_by_id('11').set_text(str(overlap))
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

print("[SUCCESS] Venn diagram updated with 86 ECM-DEGs and 12 convergent causal genes!")
