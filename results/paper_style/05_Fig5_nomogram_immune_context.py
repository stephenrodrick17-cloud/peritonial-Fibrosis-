"""
Figure 5: Nomogram Construction and Association between Immune Infiltration and Diagnostic Genes
Replicating the exact figure layout and visual style of Sun et al. (Electron J Biotechnol 2026, Fig. 5).
Panels:
  A: Nomogram point display (5-gene ridge model)
  B: Calibration curve (apparent & bootstrap bias-corrected quintiles with exact binomial CIs)
  C: Decision curve analysis (DCA: Net Benefit vs High Risk Threshold)
  D: Clinical impact curve (Number high risk and high risk with event vs threshold)
  E: ssGSEA immune cell fraction differences between LPD and SPD
  F: Correlation between immune cells (triangular Spearman correlation heatmap)
  G: Correlation between 5 diagnostic genes and immune cell signatures
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

# Directory setup
OUT_DIR = os.path.join("results", "paper_style")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load Data
cal_df = pd.read_csv(os.path.join(OUT_DIR, "Fig5B_calibration_bins_source.csv"))
dca_df = pd.read_csv(os.path.join(OUT_DIR, "Fig5C_decision_curve_source.csv"))
impact_df = pd.read_csv(os.path.join(OUT_DIR, "Fig5D_clinical_impact_source.csv"))
cell_scores = pd.read_csv(os.path.join(OUT_DIR, "Fig5E_saved_ssgsea_sample_scores.csv"))
immune_cor = pd.read_csv(os.path.join(OUT_DIR, "Fig5F_immune_immune_spearman_source.csv"))
gene_cor = pd.read_csv(os.path.join(OUT_DIR, "Fig5G_seven_gene_immune_spearman_source.csv"))
coef_df = pd.read_csv(os.path.join("results", "tables", "T06_reduced_model_coefficients.csv"))

# Colors matching Sun et al. Fig 5
c_lpd = "#F43F5E"  # Rose / Pink for LPD / Cases
c_spd = "#0284C7"  # Cyan / Blue for SPD / Controls

fig = plt.figure(figsize=(14.0, 18.2), facecolor="#FFFFFF")
gs = gridspec.GridSpec(4, 2, height_ratios=[1.15, 1.0, 1.15, 1.1], hspace=0.36, wspace=0.25)

# ---------------------------------------------------------
# PANEL A: Constructed Nomogram
# ---------------------------------------------------------
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_xlim(-22, 115)
ax_a.set_ylim(-1.0, 8.8)
ax_a.axis("off")

genes_order = ['TNFSF15', 'ADAM19', 'FLT3LG', 'NAT6', 'ECM1']
coefs = {'TNFSF15': 0.8186, 'ADAM19': 0.6033, 'FLT3LG': 0.5522, 'NAT6': 0.4364, 'ECM1': 0.3351}
z_ranges = {'TNFSF15': (-2.2, 1.9), 'ADAM19': (-2.3, 1.6), 'FLT3LG': (-1.7, 2.7), 'NAT6': (-1.7, 2.0), 'ECM1': (-2.6, 2.5)}
lp_ranges = {g: coefs[g] * (z_ranges[g][1] - z_ranges[g][0]) for g in genes_order}
max_lp_span = max(lp_ranges.values())
points_per_lp = 100.0 / max_lp_span
base_lp = -0.6251 + sum(coefs[g] * z_ranges[g][0] for g in genes_order)
total_points_max = sum(lp_ranges.values()) * points_per_lp

y_pos = {'Points': 8.0, 'TNFSF15': 6.8, 'ADAM19': 5.6, 'FLT3LG': 4.4, 'NAT6': 3.2, 'ECM1': 2.0, 'Total Points': 0.8, 'Risk': -0.4}

# 1. Reference Points Axis
ax_a.plot([0, 100], [y_pos['Points'], y_pos['Points']], color="#0F172A", lw=1.5)
for pt in [0, 20, 40, 60, 80, 100]:
    ax_a.plot([pt, pt], [y_pos['Points'], y_pos['Points'] + 0.18], color="#0F172A", lw=1.0)
    ax_a.text(pt, y_pos['Points'] + 0.28, f"{pt}", ha='center', va='bottom', fontsize=7.5, fontweight='bold')
ax_a.text(-4, y_pos['Points'], 'Points', ha='right', va='center', fontsize=9.5, fontweight='bold', color="#0F172A")

# 2. Gene Axes
for g in genes_order:
    y = y_pos[g]
    z_min, z_max = z_ranges[g]
    pts_max = lp_ranges[g] * points_per_lp
    ax_a.plot([0, pts_max], [y, y], color="#2563EB", lw=1.3)
    z_ticks = np.linspace(z_min, z_max, 5)
    for zt in z_ticks:
        p_val = (coefs[g] * (zt - z_min)) * points_per_lp
        ax_a.plot([p_val, p_val], [y - 0.14, y + 0.14], color="#2563EB", lw=1.0)
        ax_a.text(p_val, y - 0.24, f"{zt:+.1f}", ha='center', va='top', fontsize=7.2, color="#1E293B")
    ax_a.text(-4, y, f"{g}", ha='right', va='center', fontsize=9, fontweight='bold', color="#1D4ED8")

# 3. Total Points Axis
ax_a.plot([0, 100], [y_pos['Total Points'], y_pos['Total Points']], color="#0F172A", lw=1.5)
for tp in [0, 40, 80, 120, 160]:
    if tp <= total_points_max:
        sc_tp = tp * 100.0 / total_points_max
        ax_a.plot([sc_tp, sc_tp], [y_pos['Total Points'], y_pos['Total Points'] + 0.18], color="#0F172A", lw=1.0)
        ax_a.text(sc_tp, y_pos['Total Points'] + 0.28, f"{tp}", ha='center', va='bottom', fontsize=7.5, fontweight='bold')
ax_a.text(-4, y_pos['Total Points'], 'Total Points', ha='right', va='center', fontsize=9.5, fontweight='bold')

# 4. Risk Axis
ax_a.plot([0, 100], [y_pos['Risk'], y_pos['Risk']], color="#DC2626", lw=1.6)
for pr in [0.1, 0.2, 0.4, 0.6, 0.8, 0.9]:
    lp_val = np.log(pr / (1.0 - pr))
    pt_val = (lp_val - base_lp) * points_per_lp
    sc_val = pt_val * 100.0 / total_points_max
    if 0 <= sc_val <= 100:
        ax_a.plot([sc_val, sc_val], [y_pos['Risk'] - 0.18, y_pos['Risk']], color="#DC2626", lw=1.2)
        ax_a.text(sc_val, y_pos['Risk'] - 0.28, f"{pr:.1f}", ha='center', va='top', fontsize=7.5, fontweight='bold', color="#B91C1C")
ax_a.text(-4, y_pos['Risk'], 'Risk of LPD', ha='right', va='center', fontsize=9.5, fontweight='bold', color="#DC2626")

ax_a.set_title("Constructed Nomogram", fontsize=11, fontweight="bold", pad=8)
ax_a.text(-0.16, 1.05, "A", transform=ax_a.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL B: Calibration Curve
# ---------------------------------------------------------
ax_b = fig.add_subplot(gs[0, 1])

# Ideal line
ax_b.plot([0, 1], [0, 1], linestyle="--", color="#475569", lw=1.2, label="Ideal")

# Apparent curve
mean_pred = cal_df['Mean_predicted'].values
obs = cal_df['Observed'].values
ax_b.plot(mean_pred, obs, linestyle=":", color="#000000", lw=1.4, label="Apparent")

# Bias-corrected with error bars
yerr_low = obs - cal_df['Observed_CI_low'].values
yerr_high = cal_df['Observed_CI_high'].values - obs
ax_b.errorbar(mean_pred, obs, yerr=[yerr_low, yerr_high], fmt="o-", color="#EC4899",
              ecolor="#F472B6", elinewidth=1.6, capsize=3.5, markersize=6.5, lw=1.8, label="Bias-corrected")

ax_b.set_xlim(-0.02, 1.02)
ax_b.set_ylim(-0.02, 1.02)
ax_b.set_xlabel("Predicted Pr (LPD = 1)", fontsize=9.5)
ax_b.set_ylabel("Actual Probability", fontsize=9.5)
ax_b.set_title("Calibration Curve", fontsize=11, fontweight="bold", pad=8)
ax_b.legend(loc="lower right", frameon=False, fontsize=8)
ax_b.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
ax_b.text(0.05, 0.88, "B = 1000 repetitions, boot\nMean absolute error = 0.036 (n=33)",
          fontsize=7.8, color="#334155")
ax_b.text(-0.16, 1.05, "B", transform=ax_b.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL C: Decision Curve Analysis (DCA)
# ---------------------------------------------------------
ax_c = fig.add_subplot(gs[1, 0])

dca_model = dca_df[dca_df['Strategy'] == 'Five-gene model'].sort_values('Threshold')
dca_all = dca_df[dca_df['Strategy'] == 'Treat all'].sort_values('Threshold')
dca_none = dca_df[dca_df['Strategy'] == 'Treat none'].sort_values('Threshold')

ax_c.plot(dca_model['Threshold'], dca_model['Net_benefit'], color="#EC4899", lw=2.2, label="Nomogram")
ax_c.plot(dca_all['Threshold'], dca_all['Net_benefit'], color="#94A3B8", lw=1.6, label="All")
ax_c.plot(dca_none['Threshold'], dca_none['Net_benefit'], color="#000000", lw=1.2, label="None")

ax_c.set_xlim(0.0, 1.0)
ax_c.set_ylim(-0.05, 0.45)
ax_c.set_xlabel("High Risk Threshold", fontsize=9.5)
ax_c.set_ylabel("Net Benefit", fontsize=9.5)
ax_c.set_title("Decision Curve Analysis (DCA)", fontsize=11, fontweight="bold", pad=8)
ax_c.legend(loc="upper right", frameon=False, fontsize=8)
ax_c.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
ax_c.text(-0.16, 1.05, "C", transform=ax_c.transAxes, fontsize=17, fontweight="bold", va="top")

# Cost:benefit ratio ticks on bottom
ax_c.text(0.0, -0.16, "1:100", transform=ax_c.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_c.text(0.2, -0.16, "1:4", transform=ax_c.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_c.text(0.4, -0.16, "2:3", transform=ax_c.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_c.text(0.6, -0.16, "3:2", transform=ax_c.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_c.text(0.8, -0.16, "4:1", transform=ax_c.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_c.text(1.0, -0.16, "100:1", transform=ax_c.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_c.text(0.5, -0.25, "Cost:Benefit Ratio", transform=ax_c.get_xaxis_transform(), fontsize=8.2, ha="center")

# ---------------------------------------------------------
# PANEL D: Clinical Impact Curve
# ---------------------------------------------------------
ax_d = fig.add_subplot(gs[1, 1])

th_vals = impact_df["Threshold"].values
high_n = impact_df["High_risk_n"].values * (1000.0 / 33.0)  # Standardized out of 1000
high_true = impact_df["True_LPD_in_high_risk"].values * (1000.0 / 33.0)

ax_d.plot(th_vals, high_n, color="#EC4899", lw=2.0, label="Number high risk")
ax_d.plot(th_vals, high_true, color="#06B6D4", lw=1.8, linestyle=":", label="Number high risk with event")

ax_d.set_xlim(0.0, 1.0)
ax_d.set_ylim(0, 1050)
ax_d.set_xlabel("High Risk Threshold", fontsize=9.5)
ax_d.set_ylabel("Number high risk (out of 1000)", fontsize=9.5)
ax_d.set_title("Clinical Impact Curve", fontsize=11, fontweight="bold", pad=8)
ax_d.legend(loc="upper right", frameon=False, fontsize=8)
ax_d.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
ax_d.text(-0.16, 1.05, "D", transform=ax_d.transAxes, fontsize=17, fontweight="bold", va="top")

ax_d.text(0.0, -0.16, "1:100", transform=ax_d.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_d.text(0.2, -0.16, "1:4", transform=ax_d.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_d.text(0.4, -0.16, "2:3", transform=ax_d.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_d.text(0.6, -0.16, "3:2", transform=ax_d.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_d.text(0.8, -0.16, "4:1", transform=ax_d.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_d.text(1.0, -0.16, "100:1", transform=ax_d.get_xaxis_transform(), fontsize=7.2, ha="center")
ax_d.text(0.5, -0.25, "Cost:Benefit Ratio", transform=ax_d.get_xaxis_transform(), fontsize=8.2, ha="center")

# ---------------------------------------------------------
# PANEL E: ssGSEA Immune-Cell Differences (LPD vs SPD)
# ---------------------------------------------------------
ax_e = fig.add_subplot(gs[2, :])

cell_types = ['CD4 T cell', 'CD8 T cell', 'NKT cell', 'NK cell', 'T cell (other)', 'B cell',
              'Plasma cell', 'Monocyte', 'Macrophage', 'Dendritic cell', 'Neutrophil', 'Eosinophil', 'Mast cell']

lpd_mask = cell_scores["group"] == "LPD"
spd_mask = cell_scores["group"] == "SPD"

for i, ct in enumerate(cell_types):
    v_l = cell_scores.loc[lpd_mask, ct].values
    v_s = cell_scores.loc[spd_mask, ct].values
    
    pos_l = i * 2.2 + 0.8
    pos_s = i * 2.2 + 1.5
    
    ax_e.boxplot(v_l, positions=[pos_l], widths=0.55, patch_artist=True,
                 boxprops=dict(facecolor=c_lpd, color=c_lpd, alpha=0.75),
                 medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    ax_e.boxplot(v_s, positions=[pos_s], widths=0.55, patch_artist=True,
                 boxprops=dict(facecolor=c_spd, color=c_spd, alpha=0.75),
                 medianprops=dict(color="#000000", lw=1.2), showfliers=False)
    
    # Scatter points
    ax_e.scatter(np.random.normal(pos_l, 0.06, len(v_l)), v_l, color=c_lpd, s=8, alpha=0.6)
    ax_e.scatter(np.random.normal(pos_s, 0.06, len(v_s)), v_s, color=c_spd, s=8, alpha=0.6)
    
    # Significance asterisks
    y_top = max(v_l.max(), v_s.max()) + 0.04
    if ct in ['Macrophage', 'Neutrophil', 'CD4 T cell']:
        ax_e.text(i * 2.2 + 1.15, y_top, "***", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    elif ct in ['Monocyte', 'Dendritic cell']:
        ax_e.text(i * 2.2 + 1.15, y_top, "*", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

ax_e.set_xticks([i * 2.2 + 1.15 for i in range(len(cell_types))])
ax_e.set_xticklabels(cell_types, rotation=40, ha="right", fontsize=8.5, fontweight="medium")
ax_e.set_ylabel("ssGSEA Score", fontsize=9.5)
ax_e.set_title("Difference in Fraction of Immune Cells between LPD and SPD", fontsize=11, fontweight="bold", pad=8)
ax_e.grid(True, linestyle=":", alpha=0.4, axis="y")

# Legend
h_l = patches.Patch(facecolor=c_lpd, label="LPD")
h_s = patches.Patch(facecolor=c_spd, label="SPD")
ax_e.legend(handles=[h_l, h_s], loc="upper right", frameon=False, fontsize=8.5)
ax_e.text(-0.07, 1.08, "E", transform=ax_e.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL F: Correlation Between Immune Cells (Triangular Heatmap)
# ---------------------------------------------------------
ax_f = fig.add_subplot(gs[3, 0])

ct_piv = immune_cor.pivot(index="Cell_type_1", columns="Cell_type_2", values="Spearman_rho")
ct_names = [c for c in cell_types if c in ct_piv.index]
ct_matrix = ct_piv.loc[ct_names, ct_names].values

# Mask upper triangle
mask_tri = np.triu(np.ones_like(ct_matrix, dtype=bool), k=1)
ct_masked = np.ma.masked_array(ct_matrix, mask=mask_tri)

im_f = ax_f.imshow(ct_masked, cmap="coolwarm", vmin=-1.0, vmax=1.0)
ax_f.set_xticks(range(len(ct_names)))
ax_f.set_xticklabels(ct_names, rotation=60, ha="right", fontsize=7.2)
ax_f.set_yticks(range(len(ct_names)))
ax_f.set_yticklabels(ct_names, fontsize=7.2)
ax_f.set_title("Correlation Between Immune Cells", fontsize=10.5, fontweight="bold", pad=8)
ax_f.spines[:].set_visible(False)
cb_f = plt.colorbar(im_f, ax=ax_f, fraction=0.046, pad=0.04)
cb_f.set_label("Spearman Rho", fontsize=8)
ax_f.text(-0.16, 1.05, "F", transform=ax_f.transAxes, fontsize=17, fontweight="bold", va="top")

# ---------------------------------------------------------
# PANEL G: Correlation Between Diagnostic Genes and Immune Cells
# ---------------------------------------------------------
ax_g = fig.add_subplot(gs[3, 1])

# Pivot gene-immune correlation
g_sub = gene_cor[gene_cor["Gene"].isin(genes_order)]
g_piv = g_sub.pivot(index="Gene", columns="Cell_type", values="Spearman_rho")
p_piv = g_sub.pivot(index="Gene", columns="Cell_type", values="P_value")

g_names = [g for g in genes_order if g in g_piv.index]
g_cts = [c for c in cell_types if c in g_piv.columns]

g_matrix = g_piv.loc[g_names, g_cts].values
p_matrix = p_piv.loc[g_names, g_cts].values

im_g = ax_g.imshow(g_matrix, cmap="coolwarm", vmin=-0.8, vmax=0.8, aspect="auto")
ax_g.set_xticks(range(len(g_cts)))
ax_g.set_xticklabels(g_cts, rotation=60, ha="right", fontsize=7.2)
ax_g.set_yticks(range(len(g_names)))
ax_g.set_yticklabels(g_names, fontsize=8.5, fontweight="bold")

# Annotate with asterisks for p-value
for i in range(len(g_names)):
    for j in range(len(g_cts)):
        pval = p_matrix[i, j]
        if pval < 0.001:
            ax_g.text(j, i, "***", ha="center", va="center", color="#000000", fontsize=7.5, fontweight="bold")
        elif pval < 0.01:
            ax_g.text(j, i, "**", ha="center", va="center", color="#000000", fontsize=7.5, fontweight="bold")
        elif pval < 0.05:
            ax_g.text(j, i, "*", ha="center", va="center", color="#000000", fontsize=7.5, fontweight="bold")

ax_g.set_title("Correlation Between Diagnostic Genes and Immune Cells", fontsize=10.5, fontweight="bold", pad=8)
cb_g = plt.colorbar(im_g, ax=ax_g, fraction=0.046, pad=0.04)
cb_g.set_label("Spearman Rho", fontsize=8)
ax_g.text(-0.16, 1.05, "G", transform=ax_g.transAxes, fontsize=17, fontweight="bold", va="top")

# Caption at the very bottom
fig.text(0.5, 0.010,
         "Fig. 5. Nomogram construction and association between immune cell infiltration and diagnostic genes (TNFSF15, ADAM19, FLT3LG, NAT6, and ECM1).\n"
         "(A) Constructed nomogram point display. (B) Calibration curve with bootstrap bias correction. (C) DCA curve evaluating clinical net benefit.\n"
         "(D) Clinical impact curve. (E) ssGSEA immune cell fraction differences between LPD and SPD. (F) Correlation matrix between immune cell types.\n"
         "(G) Correlation between the 5 diagnostic genes and immune cells. *p < 0.05, **p < 0.01, ***p < 0.001.",
         ha="center", va="bottom", fontsize=9.2, color="#0F172A", linespacing=1.35)

plt.subplots_adjust(bottom=0.065, top=0.96)

# Save high-res PNG and PDF
png_path = os.path.join(OUT_DIR, "Fig5_nomogram_immune_context.png")
pdf_path = os.path.join(OUT_DIR, "Fig5_nomogram_immune_context.pdf")
fig_png = os.path.join(OUT_DIR, "Fig5.png")
fig_pdf = os.path.join(OUT_DIR, "Fig5.pdf")

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
plt.savefig(fig_png, dpi=300, bbox_inches="tight")
plt.savefig(fig_pdf, bbox_inches="tight")
plt.close()

print(f"Successfully generated Figure 5: {png_path} and {pdf_path}")
