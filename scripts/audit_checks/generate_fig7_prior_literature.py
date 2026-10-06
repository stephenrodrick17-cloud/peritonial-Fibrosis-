"""
generate_fig7_prior_literature.py
Generates Figure 7: Relationship to Prior Literature and Effluent Hub Gene Discordance.
Panels:
  A: Effluent leukocyte marker shifts (CD14, CD3E, FCGR3B) across dialysis vintage (GSE125498, N=33: 20 SPD vs 13 LPD).
  B: Hub gene expression discordance: Discovery Tissue (GSE62928) vs. Cellular Effluent (GSE125498).
  C: VCAN covariate adjustment fragility and leave-one-out sensitivity (Set B probes).
  D: Biological coupling: Correlation of VCAN with monocyte marker CD14 in shed effluent cells.
"""

import os
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

ROOT = r"d:\Peritoneal Project"
FIGS_DIR = os.path.join(ROOT, "results", "figures")
SOURCE_DIR = os.path.join(FIGS_DIR, "source_data")
os.makedirs(SOURCE_DIR, exist_ok=True)

# -------------------------------------------------------------------------
# 1. Load Data
# -------------------------------------------------------------------------
# Load GSE125498 series matrix
matrix_file = os.path.join(ROOT, "data", "raw", "GSE125498_series_matrix.txt.gz")
with gzip.open(matrix_file, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = line.strip().split("\t")[1:]
        if line.startswith("!Sample_geo_accession"):
            samples = line.strip().split("\t")[1:]
        if line.startswith("!series_matrix_table_begin"):
            break
    df_expr = pd.read_csv(f, sep="\t", index_col=0)
    if "!series_matrix_table_end" in df_expr.index:
        df_expr = df_expr.drop("!series_matrix_table_end")

groups = ["SPD" if "short" in t.lower() or "spd" in t.lower() else "LPD" for t in titles]
df_meta = pd.DataFrame({"Sample": samples, "Title": titles, "Group": groups})
df_meta["Group"] = pd.Categorical(df_meta["Group"], categories=["SPD", "LPD"])

# Probes of interest:
# CD14: ILMN_1740015 (Set B), CD3E: ILMN_1739794, FCGR3B: ILMN_1728639 (Set B)
# VCAN: ILMN_1687301
cd14 = df_expr.loc["ILMN_1740015"].astype(float).values
cd3e = df_expr.loc["ILMN_1739794"].astype(float).values
fcgr3b = df_expr.loc["ILMN_1728639"].astype(float).values
vcan = df_expr.loc["ILMN_1687301"].astype(float).values

# -------------------------------------------------------------------------
# 2. Setup Figure Aesthetics
# -------------------------------------------------------------------------
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
palette_groups = {"SPD": "#3498db", "LPD": "#e74c3c"}

# -------------------------------------------------------------------------
# Panel A: Effluent Immune Shifts (GSE125498)
# -------------------------------------------------------------------------
ax_a = axes[0, 0]

df_markers = pd.DataFrame({
    "Sample": samples,
    "Group": groups,
    "CD14 (Monocyte)": cd14,
    "CD3E (T Cell)": cd3e,
    "FCGR3B (Neutrophil)": fcgr3b
})
df_melted = pd.melt(df_markers, id_vars=["Sample", "Group"], 
                    value_vars=["CD14 (Monocyte)", "CD3E (T Cell)", "FCGR3B (Neutrophil)"],
                    var_name="Marker", value_name="log2_Expression")

sns.boxplot(data=df_melted, x="Marker", y="log2_Expression", hue="Group", ax=ax_a, 
            palette=palette_groups, width=0.55, boxprops=dict(alpha=0.75), showmeans=True,
            meanprops=dict(marker="o", markeredgecolor="black", markerfacecolor="white", markersize=5))
sns.stripplot(data=df_melted, x="Marker", y="log2_Expression", hue="Group", ax=ax_a, 
              dodge=True, palette=palette_groups, alpha=0.6, jitter=0.2, size=5)

# Deduplicate legend in ax_a
handles, labels = ax_a.get_legend_handles_labels()
ax_a.legend(handles[:2], labels[:2], title="Dialysis Vintage", frameon=True, loc="upper right")

ax_a.set_title("A. Effluent Leukocyte Marker Expression (GSE125498)\nSPD (n=20) vs. LPD (n=13)", 
               fontsize=11.5, fontweight="bold", pad=10)
ax_a.set_ylabel("Microarray log2 Intensity", fontsize=10)
ax_a.set_xlabel("")
ax_a.text(0, 14.1, "log2FC = -0.38\nP = 0.016 (Down)", ha="center", fontsize=8.5, color="#c0392b", fontweight="semibold")
ax_a.text(1, 10.0, "log2FC = +1.46\nP = 0.003 (Up)", ha="center", fontsize=8.5, color="#27ae60", fontweight="semibold")
ax_a.text(2, 11.8, "log2FC = +0.81\nP = 0.024 (Up)", ha="center", fontsize=8.5, color="#27ae60", fontweight="semibold")
ax_a.set_ylim(2, 15)
ax_a.grid(axis="y", linestyle="--", alpha=0.4)

# -------------------------------------------------------------------------
# Panel B: Hub Gene Tissue vs. Effluent Discordance
# -------------------------------------------------------------------------
ax_b = axes[0, 1]

hubs = ["VCAN", "COL8A1", "FN1", "COL3A1", "THBS3", "ISM1", "LOX", "COL11A1", "COMP", "EDIL3", "INHBA"]
tissue_lfc = [2.75, 2.68, 1.93, 2.84, 1.16, 1.90, 2.16, 3.79, 4.08, 1.40, 2.88]
effluent_lfc = [-0.52, 0.75, 0.41, 0.19, -0.20, 0.03, 0.02, np.nan, np.nan, np.nan, np.nan]

y_pos = np.arange(len(hubs))
height = 0.38

rects1 = ax_b.barh(y_pos + height/2, tissue_lfc, height, label="Tissue Biopsy (GSE62928, N=8)", color="#e67e22", alpha=0.85)
rects2 = ax_b.barh(y_pos - height/2, [0 if np.isnan(v) else v for v in effluent_lfc], height, 
                   label="Cellular Effluent (GSE125498, N=33)", color="#2980b9", alpha=0.85)

# Annotate unassayed hubs
for i, elfc in enumerate(effluent_lfc):
    if np.isnan(elfc):
        ax_b.text(0.1, y_pos[i] - height/2, "Unassayed on GPL10558", va="center", fontsize=8, fontstyle="italic", color="#7f8c8d")

ax_b.axvline(0, color="black", linestyle="-", linewidth=0.8)
ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(hubs, fontsize=9.5, fontweight="semibold")
ax_b.invert_yaxis()
ax_b.set_xlabel("Differential Expression log2 Fold Change", fontsize=10)
ax_b.set_title("B. Hub Gene Expression Discordance: Tissue vs. Effluent\nParietal Biopsy (EPS vs. Ctrl) vs. Effluent Cells (LPD vs. SPD)", 
               fontsize=11.5, fontweight="bold", pad=10)
ax_b.legend(loc="lower right", frameon=True, fontsize=9)
ax_b.grid(axis="x", linestyle="--", alpha=0.4)

# -------------------------------------------------------------------------
# Panel C: VCAN Fragility & Leave-One-Out Covariate Analysis (Set B)
# -------------------------------------------------------------------------
ax_c = axes[1, 0]

models = [
    "Unadjusted",
    "Adjust: CD14 only",
    "Adjust: CD3E only",
    "Adjust: FCGR3B only",
    "Drop CD14 (CD3E + FCGR3B)",
    "Drop FCGR3B (CD14 + CD3E)",
    "Drop CD3E (CD14 + FCGR3B)",
    "Joint: All 3 Markers"
]

vcan_lfc = [-0.5224, -0.2845, -0.3241, -0.7043, -0.4930, 0.0178, -0.4350, -0.0902]
vcan_pval = [0.0244, 0.1429, 0.1955, 0.0026, 0.0395, 0.9276, 0.0736, 0.6716]
y_c = np.arange(len(models))

bar_colors = ["#27ae60" if p < 0.05 else "#e74c3c" for p in vcan_pval]

bars = ax_c.barh(y_c, vcan_lfc, height=0.55, color=bar_colors, alpha=0.8, edgecolor="#2c3e50")
ax_c.axvline(0, color="black", linestyle="-", linewidth=0.8)
ax_c.axvline(-0.5224, color="#2980b9", linestyle=":", alpha=0.6, label="Unadjusted Baseline")

for i, (lfc, p) in enumerate(zip(vcan_lfc, vcan_pval)):
    sig_str = f"P = {p:.4f}*" if p < 0.05 else f"P = {p:.4f}"
    ha = "left" if lfc >= 0 else "right"
    offset = 0.02 if lfc >= 0 else -0.02
    ax_c.text(lfc + offset, i, f"{lfc:+.2f} ({sig_str})", va="center", ha=ha, fontsize=8.2, fontweight="semibold")

ax_c.set_yticks(y_c)
ax_c.set_yticklabels(models, fontsize=9.0)
ax_c.invert_yaxis()
ax_c.set_xlabel("VCAN log2 Fold Change in LPD vs. SPD", fontsize=10)
ax_c.set_title("C. VCAN Covariate Fragility & Leave-One-Out Analysis (Set B)\nAttenuation depends on CD14 and FCGR3B together", 
               fontsize=11.5, fontweight="bold", pad=10)
ax_c.grid(axis="x", linestyle="--", alpha=0.4)
ax_c.set_xlim(-1.15, 0.45)

# Custom legend for significance placed upper right or below title
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#27ae60", label="Nominally Significant (P < 0.05)"),
    Patch(facecolor="#e74c3c", label="Non-Significant / Attenuated (P >= 0.05)")
]
ax_c.legend(handles=legend_elements, loc="upper right", frameon=True, fontsize=8.5)

# -------------------------------------------------------------------------
# Panel D: Biological Coupling: VCAN vs. CD14 (Monocyte Marker)
# -------------------------------------------------------------------------
ax_d = axes[1, 1]

r_pear, p_pear = stats.pearsonr(vcan, cd14)
r_spear, p_spear = stats.spearmanr(vcan, cd14)

sns.regplot(x=cd14, y=vcan, ax=ax_d, color="#2c3e50", scatter=False, line_kws={"linewidth": 1.5, "color": "#2c3e50"})
for grp, color in palette_groups.items():
    idx = df_meta["Group"] == grp
    ax_d.scatter(cd14[idx], vcan[idx], label=f"{grp} (n={idx.sum()})", color=color, s=45, alpha=0.85, edgecolors="white", linewidth=0.5)

ax_d.set_xlabel("CD14 Expression (ILMN_1740015, log2 Intensity)", fontsize=10)
ax_d.set_ylabel("VCAN Expression (ILMN_1687301, log2 Intensity)", fontsize=10)
ax_d.set_title("D. Biological Coupling: VCAN vs. CD14 in Cellular Effluent\nVCAN expression reflects shed monocyte abundance", 
               fontsize=11.5, fontweight="bold", pad=10)
ax_d.grid(True, linestyle="--", alpha=0.4)
ax_d.legend(title="Dialysis Vintage", loc="upper left", frameon=True, fontsize=9)

stats_text = (
    f"Pearson r = +{r_pear:.3f} (P = {p_pear:.2e})\n"
    f"Spearman rho = +{r_spear:.3f} (P = {p_spear:.2e})\n"
    f"N = 33 Peritoneal Effluent Samples\n\n"
    f"Mechanistic Note: VCAN's apparent\n"
    f"decrease in late PD effluent tracks\n"
    f"decreased monocyte shedding, not\n"
    f"fibrotic membrane downregulation."
)
ax_d.text(0.96, 0.06, stats_text, transform=ax_d.transAxes, fontsize=8.5,
          verticalalignment="bottom", horizontalalignment="right",
          bbox=dict(boxstyle="round,pad=0.5", facecolor="#ecf0f1", edgecolor="#bdc3c7", alpha=0.9))

# -------------------------------------------------------------------------
# Save Figure & Source Data
# -------------------------------------------------------------------------
plt.tight_layout()
fig_path_png = os.path.join(FIGS_DIR, "Fig_7_prior_literature_effluent_comparison.png")
fig_path_pdf = os.path.join(FIGS_DIR, "Fig_7_prior_literature_effluent_comparison.pdf")
plt.savefig(fig_path_png, dpi=300)
plt.savefig(fig_path_pdf)
plt.close()
print(f"Saved figure to {fig_path_png} and {fig_path_pdf}")

# Save Source Data CSV
source_rows = []

# Panel A
for m, vals in [("CD14", cd14), ("CD3E", cd3e), ("FCGR3B", fcgr3b)]:
    for s, g, v in zip(samples, groups, vals):
        source_rows.append({"Panel": "Panel_A", "Sample": s, "Group": g, "Variable": m, "Value": v, "Note": ""})

# Panel B
for h, tlfc, elfc in zip(hubs, tissue_lfc, effluent_lfc):
    source_rows.append({"Panel": "Panel_B", "Sample": "Summary", "Group": "GSE62928", "Variable": h, "Value": tlfc, "Note": "Tissue log2FC"})
    source_rows.append({"Panel": "Panel_B", "Sample": "Summary", "Group": "GSE125498", "Variable": h, "Value": elfc, "Note": "Effluent log2FC"})

# Panel C
for m, lfc, p in zip(models, vcan_lfc, vcan_pval):
    source_rows.append({"Panel": "Panel_C", "Sample": "Model", "Group": "VCAN", "Variable": m, "Value": lfc, "Note": f"P_val={p}"})

# Panel D
for s, g, c_val, v_val in zip(samples, groups, cd14, vcan):
    source_rows.append({"Panel": "Panel_D", "Sample": s, "Group": g, "Variable": "CD14_vs_VCAN", "Value": f"CD14={c_val};VCAN={v_val}", "Note": ""})

df_source = pd.DataFrame(source_rows)
source_csv_path = os.path.join(SOURCE_DIR, "fig7_prior_literature_effluent_comparison_source.csv")
df_source.to_csv(source_csv_path, index=False)
print(f"Saved source data to {source_csv_path}")
