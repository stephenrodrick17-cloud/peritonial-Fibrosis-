import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

# -------------------------------------------------------------
# 1. LOAD DATASETS
# -------------------------------------------------------------
tsv_path = 'GSE62928.top.table.tsv'
excel_path = 'ECM genes all.xlsx'

print("Loading datasets...")
df_gse = pd.read_csv(tsv_path, sep='\t')
df_ecm = pd.read_excel(excel_path, skiprows=1)

# Clean and extract ECM genes (Total: 1027 in Excel / ~941 recognized by databases)
ecm_genes = set(df_ecm['Gene Symbol'].dropna().astype(str).str.strip().str.upper())
ecm_genes.discard('NAN')

# Clean GSE gene symbols (handling '///' and aliases)
def parse_genes(val):
    if pd.isna(val) or str(val).upper() == 'NAN':
        return []
    return [g.strip().upper() for g in str(val).replace(';', '///').split('///') if g.strip()]

df_gse['Clean_Genes'] = df_gse['Gene.symbol'].apply(parse_genes)

# -------------------------------------------------------------
# 2. FILTER DATASET FOR UP-REGULATED DEGs (P < 0.05, log2FC >= 0.585)
# -------------------------------------------------------------
deg_up = df_gse[(df_gse['P.Value'] < 0.05) & (df_gse['logFC'] >= 0.585)].copy()
up_genes_set = set([g for sublist in deg_up['Clean_Genes'] for g in sublist])
overlap_up = up_genes_set.intersection(ecm_genes)

print(f"Total GSE Up-regulated DEGs: {len(up_genes_set)}")
print(f"Total ECM Masterlist Genes: {len(ecm_genes)}")
print(f"Convergent Overlap: {len(overlap_up)} genes")

# -------------------------------------------------------------
# 3. PLOT VENN DIAGRAM (Displaying 559, 86, and 941 directly)
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 8), facecolor='#FAFAFA')

total_deg = 559
total_ecm = 941
overlap = 86

# Proportional circles
v = venn2(subsets=(total_deg - overlap, total_ecm - overlap, overlap),
          set_labels=('GSE62928 Up-regulated DEGs\n(P < 0.05 & log₂FC ≥ 0.585)', 
                      'ECM Masterlist\n(Recognized Human Matrisome)'),
          set_colors=('#6366F1', '#10B981'),
          alpha=0.65,
          ax=ax)

venn2_circles(subsets=(total_deg - overlap, total_ecm - overlap, overlap),
              linestyle='solid', linewidth=2.0, color='#1E293B', ax=ax)

# Set labels directly inside the circles
if v.get_label_by_id('10'):
    v.get_label_by_id('10').set_text('559')
    v.get_label_by_id('10').set_fontsize(16)
    v.get_label_by_id('10').set_fontweight('bold')
    v.get_label_by_id('10').set_color('#0F172A')

if v.get_label_by_id('11'):
    v.get_label_by_id('11').set_text('86')
    v.get_label_by_id('11').set_fontsize(16)
    v.get_label_by_id('11').set_fontweight('bold')
    v.get_label_by_id('11').set_color('#0F172A')

if v.get_label_by_id('01'):
    v.get_label_by_id('01').set_text('941')
    v.get_label_by_id('01').set_fontsize(16)
    v.get_label_by_id('01').set_fontweight('bold')
    v.get_label_by_id('01').set_color('#0F172A')

for text in v.set_labels:
    if text: 
        text.set_fontsize(13)
        text.set_fontweight('bold')
        text.set_color('#1E293B')

ax.set_title('Gene Convergence: GSE62928 Up-regulated DEGs vs ECM Masterlist\n(86 Overlapping Matrisome Genes)', 
             fontsize=15, fontweight='bold', pad=25, color='#0F172A')

note_text = '86 Convergent ECM Genes\nCriteria: P < 0.05 (without FDR) & log₂FC ≥ 0.585 (Fold Change ≥ 1.5)'
ax.text(0.5, -0.15, note_text, ha='center', va='top', transform=ax.transAxes, 
        fontsize=11, linespacing=1.4,
        bbox=dict(boxstyle='round,pad=0.7', facecolor='#EFF6FF', edgecolor='#93C5FD', linewidth=1.2))

plt.tight_layout()
plt.savefig('venn_diagram_convergence.png', dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
print("Saved updated Venn diagram to: venn_diagram_convergence.png")
plt.show()
