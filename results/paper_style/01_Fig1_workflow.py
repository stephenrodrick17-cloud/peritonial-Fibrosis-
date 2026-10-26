"""
Figure 1: Study Analysis Workflow Flowchart
Reproducing the workflow flow chart of Sun et al. (Electron J Biotechnol 2026)
for Peritoneal Dialysis Effluent (GSE125498: 13 LPD vs 20 SPD).
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Directory setup
OUT_DIR = os.path.join("results", "paper_style")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Source Data
source_data = pd.DataFrame([
    {"Step": 1, "Phase": "Discovery Cohort", "Dataset": "GSE125498", "Platform": "Illumina HumanHT-12 v4", 
     "Sample_N": 33, "LPD_Cases": 13, "SPD_Controls": 20, "Gene_Universe": 11741, "Status": "Reused"},
    {"Step": 2, "Phase": "Differential Expression", "Method": "limma empirical Bayes", "Threshold": "|log2FC| >= 0.585, P < 0.05",
     "Total_DEGs": 813, "Up_DEGs": 496, "Down_DEGs": 317, "FDR_05_DEGs": 143, "Status": "Reused"},
    {"Step": 3, "Phase": "Matrisome & WGCNA", "Reference": "Naba Matrisome (1,027 genes)", "Expressed_Matrisome": 428,
     "Matrisome_DEGs": 44, "Up_Matrisome_DEGs": 31, "WGCNA_Modules": 14, "WGCNA_Input": 5001, "Status": "Reused"},
    {"Step": 4, "Phase": "Machine Learning Selection", "Algorithms": "LASSO, SVM-RFE, Random Forest",
     "Pool": "38/39 ECM-DEG non-grey WGCNA genes", "LASSO_Genes": 11, "SVM_RFE_Genes": 10, "RF_Genes": 10, "Status": "Reused"},
    {"Step": 5, "Phase": "Consensus Biomarkers", "Panel_Type": "Majority Vote (>=2/3)",
     "Hub_Genes": "TNFSF15, ADAM19, FLT3LG, NAT6, ECM1 (+ LTB, EBI3)", "Panel_Size": "5-gene / 7-gene", "Status": "Reused"},
    {"Step": 6, "Phase": "Nomogram & Immune Profiling", "Model": "Five-gene ridge penalized logistic regression",
     "Validation_Tools": "Calibration quintiles, Decision curve analysis (DCA), Clinical impact",
     "Immune_Deconvolution": "ssGSEA (CellMarker 2024 / LM22) & Spearman gene-immune correlation", "Status": "Reused"},
    {"Step": 7, "Phase": "Functional Network & Pharmacology", "Network": "GeneMANIA (20 interactors, automatic weighting)",
     "Perturbation": "DSigDB candidate drug signatures", "Pathway": "Single-gene preranked KEGG GSEA", "Status": "Reused"},
    {"Step": 8, "Phase": "Molecular Subtypes", "Cohort": "LPD discovery samples (n=13)", "Method": "Consensus clustering (k=2)",
     "Metrics": "CDF, Delta area, PAC, PCA separation, Subtype immune infiltration", "Status": "Newly run (exploratory)"},
    {"Step": 9, "Phase": "Validation & Audit", "Cohorts": "GSE248762 (LV vs SV), GSE130888 (Long vs Short), GSE62928 (EPS biopsy)",
     "Verification_Rule": "Direction concordance & detectability audit; no pooling", "Status": "Audited"}
])
source_data.to_csv(os.path.join(OUT_DIR, "Fig1_source_data.csv"), index=False)

# 2. Plotting Flowchart
fig = plt.figure(figsize=(14, 11), facecolor="#FFFFFF")
ax = fig.add_subplot(111)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# Helper function to draw rounded boxes
def draw_box(ax, x, y, w, h, title, content_lines, box_color="#EBF3FA", border_color="#2B6CB0", title_color="#1A365D", text_color="#2D3748"):
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.8,rounding_size=1.5",
                                  facecolor=box_color, edgecolor=border_color, linewidth=1.8, zorder=2)
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h - 2.8, title, ha="center", va="center",
            fontsize=10.5, fontweight="bold", color=title_color, zorder=3)
    
    y_text = y + h - 5.8
    line_spacing = (h - 6.5) / max(len(content_lines), 1)
    for line in content_lines:
        ax.text(x + w / 2, y_text, line, ha="center", va="center",
                fontsize=8.5, color=text_color, zorder=3)
        y_text -= line_spacing

def draw_arrow(ax, x1, y1, x2, y2, color="#4A5568"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8, mutation_scale=14),
                zorder=4)

# Header
ax.text(50, 96.5, "Peritoneal Dialysis Effluent Diagnostic Biomarker & Subtyping Framework",
        ha="center", va="center", fontsize=15, fontweight="bold", color="#1A365D")
ax.text(50, 93.5, "Workflow replicating the matrisome-based biomarker discovery pipeline of Sun et al. (2026)",
        ha="center", va="center", fontsize=10, fontstyle="italic", color="#4A5568")

# Column 1 (Left): Data Discovery, DEGs & Matrisome WGCNA (x: 5 to 45)
draw_box(ax, 6, 76, 40, 13, "1. Discovery Cohort (GSE125498)",
         ["Platform: Illumina HumanHT-12 v4 beadchip",
          "Samples: N = 33 peritoneal effluent cell pellets",
          "Groups: 13 Long-term PD (LPD) vs 20 Short-term PD (SPD)",
          "Background universe: 11,741 unique annotated genes"],
         box_color="#EFF6FF", border_color="#3B82F6", title_color="#1E40AF")

draw_arrow(ax, 26, 76, 26, 68)

draw_box(ax, 6, 55, 40, 13, "2. Differential Expression Analysis (limma)",
         ["Linear models & empirical Bayes moderated t-statistics",
          "Thresholds: |log2FC| >= 0.585 (1.5-fold) & nominal P < 0.05",
          "Identified: 813 total DEGs (496 Up in LPD, 317 Down in LPD)",
          "Genome-wide Benjamini-Hochberg FDR < 0.05: 143 genes"],
         box_color="#F0FDF4", border_color="#22C55E", title_color="#15803D")

draw_arrow(ax, 26, 55, 26, 47)

draw_box(ax, 6, 34, 40, 13, "3. Matrisome Filtration & WGCNA",
         ["Naba reference matrisome: 1,027 genes (428 in universe)",
          "Matrisome DEGs: 44 genes (31 Up, 13 Down; 1.71x enrichment, P=0.0024)",
          "WGCNA: 5,001 variable genes clustered into 14 dynamic modules",
          "Trait correlation: Top pro-fibrotic module MEmagenta (r=0.48, P=0.0047)"],
         box_color="#FEF3C7", border_color="#F59E0B", title_color="#B45309")

draw_arrow(ax, 26, 34, 26, 26)

draw_box(ax, 6, 13, 40, 13, "4. Multi-Algorithm Machine Learning Selection",
         ["Candidate feature pool: 38/39 ECM-DEG ∩ non-grey WGCNA genes",
          "• LASSO L1-regularized logistic regression (11 selected)",
          "• SVM-RFE 10-fold recursive feature elimination (10 selected)",
          "• Random Forest mean decrease in Gini importance (10 selected)",
          "Majority-vote consensus (>=2/3 algorithms): 8 core biomarkers"],
         box_color="#FDF2F8", border_color="#EC4899", title_color="#BE185D")

# Cross Arrow from Col 1 to Col 2
draw_arrow(ax, 46, 19.5, 54, 19.5)

# Column 2 (Right): Clinical Modeling, Biology & Subtypes (x: 54 to 94)
draw_box(ax, 54, 13, 40, 13, "5. Consensus Panel & Clinical Nomogram",
         ["Diagnostic panel: TNFSF15, ADAM19, FLT3LG, NAT6, ECM1 (+ LTB, EBI3)",
          "Apparent 5-gene ridge model C-index / AUC: 0.9808 (95% CI: 0.94-1.00)",
          "Validation: Equal-count risk quintiles with exact binomial CIs",
          "Decision curve analysis (DCA) & Clinical impact curve"],
         box_color="#F5F3FF", border_color="#8B5CF6", title_color="#6D28D9")

draw_arrow(ax, 74, 26, 74, 34)

draw_box(ax, 54, 34, 40, 13, "6. Microenvironment & Immune Infiltration",
         ["ssGSEA enrichment across leukocyte and stromal cell markers",
          "Differential infiltration: LPD vs SPD with BH-FDR adjustment",
          "Immune-immune cell correlation architecture",
          "Spearman rank correlation: Diagnostic genes vs immune signatures"],
         box_color="#ECFEFF", border_color="#06B6D4", title_color="#0E7490")

draw_arrow(ax, 74, 47, 74, 55)

draw_box(ax, 54, 55, 40, 13, "7. Functional Networks & Pharmacogenomics",
         ["GeneMANIA functional association network (20 interactors)",
          "DSigDB candidate perturbational drug signatures (adjusted P < 0.05)",
          "Single-gene GSEA: KEGG pro-fibrotic pathways per biomarker",
          "Biological role: ECM remodeling, cytokine signaling, cell adhesion"],
         box_color="#FFF7ED", border_color="#F97316", title_color="#C2410C")

draw_arrow(ax, 74, 68, 74, 76)

draw_box(ax, 54, 76, 40, 13, "8. Molecular Subtyping & External Audits",
         ["Unsupervised consensus clustering on LPD cohort (k=2 / k=3)",
          "Proportion of Ambiguous Clustering (PAC) & PCA separation",
          "Subtype-specific biomarker expression & immune score divergence",
          "External audit: GSE62928 (EPS biopsy), GSE248762 & GSE130888"],
         box_color="#F8FAFC", border_color="#64748B", title_color="#334155")

# Bottom Banner / Caption
rect_cap = patches.Rectangle((4, 2.5), 92, 7.5, facecolor="#F1F5F9", edgecolor="#CBD5E1", lw=1.2)
ax.add_patch(rect_cap)
ax.text(50, 7.5, "Fig. 1. Study design and computational workflow.",
        ha="center", va="center", fontsize=10.5, fontweight="bold", color="#0F172A")
ax.text(50, 4.5, "Integrative pipeline combining transcriptomic profiling, matrisome annotation, co-expression network analysis (WGCNA),\nthree-algorithm machine learning feature selection, clinical risk nomogram modeling, and molecular subtyping.",
        ha="center", va="center", fontsize=8.5, color="#334155")

plt.tight_layout()

# Save outputs
png_path = os.path.join(OUT_DIR, "Fig1_workflow.png")
pdf_path = os.path.join(OUT_DIR, "Fig1_workflow.pdf")
fig_path = os.path.join(OUT_DIR, "Fig1.png")
fig_pdf_path = os.path.join(OUT_DIR, "Fig1.pdf")

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
plt.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.savefig(fig_pdf_path, bbox_inches="tight")
plt.close()

print(f"Successfully generated Figure 1: {png_path} and {pdf_path}")
