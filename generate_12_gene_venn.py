import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

# 1. Setup Proportional Venn Diagram
fig, ax = plt.subplots(figsize=(11, 8.5), facecolor="#F8FAFC")

# Subsets for proportional sizing
total_ecm_deg = 71
total_mr = 2636
overlap = 15
left_only = total_ecm_deg - overlap   # 56
right_only = total_mr - overlap       # 2621

v = venn2(
    subsets=(left_only, right_only, overlap),
    set_labels=(f'Peritoneal Convergent ECM-DEGs\n(Dataset ∩ Human Matrisome)\n[N = {total_ecm_deg}]', 
                f'Transcriptome-Wide MR Screen\n(Nominal Causal P < 0.05)\n[N = {total_mr:,}]'),
    set_colors=('#6366F1', '#10B981'),
    alpha=0.65,
    ax=ax
)

# Custom circle borders
circles = venn2_circles(
    subsets=(left_only, right_only, overlap),
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

# Explicitly set the number inside the left circle to 56 (71 - 15)
if v.get_label_by_id('10'):
    v.get_label_by_id('10').set_text(f'{left_only}\n(ECM-DEG only)')
    v.get_label_by_id('10').set_fontsize(14)
    v.get_label_by_id('10').set_fontweight('bold')
    v.get_label_by_id('10').set_color('#0F172A')

if v.get_label_by_id('01'):
    v.get_label_by_id('01').set_text(f'{right_only:,}\n(MR only)')
    v.get_label_by_id('01').set_fontsize(14)
    v.get_label_by_id('01').set_fontweight('bold')
    v.get_label_by_id('01').set_color('#0F172A')

if v.get_label_by_id('11'):
    v.get_label_by_id('11').set_text('15')
    v.get_label_by_id('11').set_fontsize(20)
    v.get_label_by_id('11').set_fontweight('bold')
    v.get_label_by_id('11').set_color('#B91C1C')

# Title
ax.set_title("Causal Convergence: 15 Convergent ECM-DEGs ∩ Mendelian Randomization\n(Peritoneal Fibrosis Transcriptomics ∩ eQTLGen-CKDGen TWMR)", 
             fontsize=14, fontweight='bold', color='#0F172A', pad=25)

# Annotation Box for the 15 Genes
box_text = (
    "★ 15 Convergent Causal ECM Biomarkers Identified:\n"
    "• Upregulated Matrix Drivers (n=6):\n"
    "   - P4HA2 (FDR=3.6e-7, log2FC=+1.22), ADAMTS1 (P=0.023), TNC (P=0.049), ADAM28 (P=0.031),\n"
    "   - CRISPLD2 (P=0.054, log2FC=+2.78), SPON2 (P=0.054, log2FC=+1.46)\n"
    "• Downregulated Matrix Regulators (n=9):\n"
    "   - BMP6 (FDR=0.013, log2FC=-1.48), IGFBP3 (FDR=0.023, log2FC=-1.08), WNT11 (P=0.005),\n"
    "   - CRHBP (P=0.007), AMH (P=0.019), COL4A2 (P=0.023), LTBP4 (P=0.028), FGL2 (P=0.029),\n"
    "   - FBLN5 (P=0.066, log2FC=-1.48)\n\n"
    "Criteria: Peritoneal ECM-DEGs (N=71) ∩ Transcriptome-Wide Two-Sample MR Screen (N=2,636)"
)

ax.text(
    0.5, -0.18, box_text,
    transform=ax.transAxes,
    fontsize=10.0,
    fontweight='medium',
    ha='center',
    va='top',
    color='#0F172A',
    bbox=dict(boxstyle='round,pad=0.8', facecolor='#EFF6FF', edgecolor='#93C5FD', linewidth=1.5)
)

plt.tight_layout()

# Save figures
os.makedirs("transcriptome_wide_mr_results", exist_ok=True)
plt.savefig("venn_15_causal_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.savefig("venn_mr_deg_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.savefig("transcriptome_wide_mr_results/venn_15_causal_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.savefig("transcriptome_wide_mr_results/venn_mr_deg_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())

print("[SUCCESS] Venn diagram updated with 15 convergent causal genes!")
