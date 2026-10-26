"""Composite figure for the five-gene ridge model, calibration, and DCA."""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Set global styling
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['text.color'] = '#0F172A'
plt.rcParams['axes.labelcolor'] = '#0F172A'
plt.rcParams['xtick.color'] = '#334155'
plt.rcParams['ytick.color'] = '#334155'

# Load verified tables
cal_df = pd.read_csv('results/tables/T06_discovery_calibration_bins.csv')
dca_df = pd.read_csv('results/tables/T06_discovery_decision_curve.csv')
coef_df = pd.read_csv('results/tables/T06_reduced_model_coefficients.csv')
auc_df = pd.read_csv('results/tables/T06_reduced_model_apparent_AUC.csv')

# Create Figure
fig = plt.figure(figsize=(16, 12.5), facecolor='#FFFFFF')
gs = gridspec.GridSpec(2, 2, width_ratios=[1.05, 1.05], height_ratios=[1.15, 1.0], hspace=0.35, wspace=0.25)

# Five-gene nomogram-style point display
ax_nomo = fig.add_subplot(gs[0, :])

genes_order = ['TNFSF15', 'ADAM19', 'FLT3LG', 'NAT6', 'ECM1']
coefs = {
    'TNFSF15': 0.818596,
    'ADAM19': 0.603281,
    'FLT3LG': 0.552221,
    'NAT6': 0.436441,
    'ECM1': 0.335072
}

# Empirical z-score limits in discovery cohort
z_ranges = {
    'TNFSF15': (-2.2, 1.9),
    'ADAM19': (-2.3, 1.6),
    'FLT3LG': (-1.7, 2.7),
    'NAT6': (-1.7, 2.0),
    'ECM1': (-2.6, 2.5)
}

lp_ranges = {g: coefs[g] * (z_ranges[g][1] - z_ranges[g][0]) for g in genes_order}
max_lp_span = max(lp_ranges.values())
points_per_lp = 100.0 / max_lp_span

# Base linear predictor at minimum predictor values
base_lp = -0.625137 + sum(coefs[g] * z_ranges[g][0] for g in genes_order)
total_points_max = sum(lp_ranges.values()) * points_per_lp

y_coords = {
    'Points': 8.2,
    'TNFSF15': 7.0,
    'ADAM19': 5.8,
    'FLT3LG': 4.6,
    'NAT6': 3.4,
    'ECM1': 2.2,
    'Total Points': 0.9,
    'Risk of LPD': -0.3
}

# 1. Reference Points Axis (0 to 100)
ax_nomo.plot([0, 100], [y_coords['Points'], y_coords['Points']], color='#0F172A', lw=1.8)
for pt in [0, 20, 40, 60, 80, 100]:
    ax_nomo.plot([pt, pt], [y_coords['Points'], y_coords['Points'] + 0.18], color='#0F172A', lw=1.3)
    ax_nomo.text(pt, y_coords['Points'] + 0.28, f"{pt}", ha='center', va='bottom', fontsize=9.5, fontweight='bold')
ax_nomo.text(-5, y_coords['Points'], 'Points', ha='right', va='center', fontsize=11, fontweight='bold', color='#0F172A')

# 2. Individual Gene Axes
for g in genes_order:
    y = y_coords[g]
    z_min, z_max = z_ranges[g]
    pts_max = lp_ranges[g] * points_per_lp
    ax_nomo.plot([0, pts_max], [y, y], color='#2563EB', lw=1.6)
    
    # 5 standard tick marks
    z_ticks = np.linspace(z_min, z_max, 5)
    for zt in z_ticks:
        p_val = (coefs[g] * (zt - z_min)) * points_per_lp
        ax_nomo.plot([p_val, p_val], [y - 0.15, y + 0.15], color='#2563EB', lw=1.2)
        ax_nomo.text(p_val, y - 0.26, f"{zt:+.1f}", ha='center', va='top', fontsize=8.5, color='#1E293B')
    ax_nomo.text(-5, y, f"{g} (z-score)", ha='right', va='center', fontsize=10.5, fontweight='bold', color='#1D4ED8')

# 3. Total Points Axis
ax_nomo.plot([0, total_points_max], [y_coords['Total Points'], y_coords['Total Points']], color='#0F172A', lw=1.8)
tp_ticks = np.arange(0, int(total_points_max) + 1, 40)
for tp in tp_ticks:
    ax_nomo.plot([tp, tp], [y_coords['Total Points'], y_coords['Total Points'] + 0.18], color='#0F172A', lw=1.3)
    ax_nomo.text(tp, y_coords['Total Points'] + 0.28, f"{tp}", ha='center', va='bottom', fontsize=9.5, fontweight='bold')
ax_nomo.text(-5, y_coords['Total Points'], 'Total Points', ha='right', va='center', fontsize=11, fontweight='bold', color='#0F172A')

# 4. Predicted Probability of LPD Axis
ax_nomo.plot([0, total_points_max], [y_coords['Risk of LPD'], y_coords['Risk of LPD']], color='#DC2626', lw=2.0)
prob_ticks = [0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.95]
for pr in prob_ticks:
    lp_val = np.log(pr / (1.0 - pr))
    pt_val = (lp_val - base_lp) * points_per_lp
    if 0 <= pt_val <= total_points_max:
        ax_nomo.plot([pt_val, pt_val], [y_coords['Risk of LPD'] - 0.18, y_coords['Risk of LPD']], color='#DC2626', lw=1.3)
        ax_nomo.text(pt_val, y_coords['Risk of LPD'] - 0.28, f"{pr:.2f}", ha='center', va='top', fontsize=9, fontweight='bold', color='#B91C1C')
ax_nomo.text(-5, y_coords['Risk of LPD'], 'Predicted LPD Risk', ha='right', va='center', fontsize=11, fontweight='bold', color='#DC2626')

ax_nomo.set_xlim(-24, total_points_max + 12)
ax_nomo.set_ylim(-0.9, 9.4)
ax_nomo.axis('off')
ax_nomo.set_title('Five-gene ridge model point display (GSE125498; long-term vs short-term PD)',
                  fontsize=13, fontweight='bold', loc='left', pad=14)

# Apparent calibration
ax_cal = fig.add_subplot(gs[1, 0])

# Perfect calibration line
ax_cal.plot([0, 1], [0, 1], linestyle='--', color='#94A3B8', lw=1.8, label='Perfect agreement')

# Observed vs Predicted error bars
mean_pred = cal_df['Mean_predicted'].values
obs = cal_df['Observed'].values
yerr_low = obs - cal_df['Observed_CI_low'].values
yerr_high = cal_df['Observed_CI_high'].values - obs

auc_result = auc_df.iloc[0]
ax_cal.errorbar(mean_pred, obs, yerr=[yerr_low, yerr_high], fmt='o-', color='#0284C7',
                ecolor='#38BDF8', elinewidth=2.0, capsize=4, capthick=1.5,
                markersize=8, markerfacecolor='#0284C7', markeredgecolor='#0C4A6E', markeredgewidth=1.5,
                lw=2.2, label=(
                    f"Apparent In-Sample Fit (AUC {auc_result['AUC']:.3f}; "
                    f"95% CI {auc_result['Bootstrap_low']:.3f}–{auc_result['Bootstrap_high']:.3f})"
                ))

# Shading for slight background contrast
ax_cal.fill_between([0, 1], [0, 1], [0, 1], color='#F0F9FF', alpha=0.5)

# Clear non-overlapping annotations for Q1 to Q5
bin_offsets = [
    (-15, 14),  # Q1
    (15, 14),   # Q2
    (-25, 14),  # Q3
    (12, -14),  # Q4
    (-32, -16)  # Q5
]
for i, row in cal_df.iterrows():
    ax_cal.annotate(f"{row['Bin']} (n={int(row['n'])})",
                    (row['Mean_predicted'], row['Observed']),
                    textcoords="offset points", xytext=bin_offsets[i],
                    fontsize=8.5, fontweight='bold', color='#0369A1',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='#FFFFFF', edgecolor='#BAE6FD', lw=0.8, alpha=0.9))

ax_cal.set_xlim(-0.02, 1.02)
ax_cal.set_ylim(-0.02, 1.02)
ax_cal.set_xlabel('Mean Predicted Probability of LPD', fontsize=11, fontweight='bold')
ax_cal.set_ylabel('Observed Case Fraction (LPD)', fontsize=11, fontweight='bold')
ax_cal.set_title('Apparent calibration by risk quintile', fontsize=13, fontweight='bold', loc='left', pad=10)
ax_cal.legend(loc='lower right', frameon=True, facecolor='#FFFFFF', edgecolor='#CBD5E1', fontsize=8.5)
ax_cal.grid(True, linestyle=':', alpha=0.5, color='#94A3B8')

# ==============================================================================
# Decision curve analysis
# ==============================================================================
ax_dca = fig.add_subplot(gs[1, 1])

dca_model = dca_df[dca_df['Strategy'] == 'Five-gene model'].sort_values('Threshold')
dca_all = dca_df[dca_df['Strategy'] == 'Treat all'].sort_values('Threshold')
dca_none = dca_df[dca_df['Strategy'] == 'Treat none'].sort_values('Threshold')

# Net benefit plot
ax_dca.plot(dca_model['Threshold'], dca_model['Net_benefit'], color='#0284C7', lw=2.6,
            label='Five-gene ridge model')
ax_dca.plot(dca_all['Threshold'], dca_all['Net_benefit'], color='#EA580C', lw=2.0, linestyle='--',
            label='Treat All (Assume All LPD)')
ax_dca.plot(dca_none['Threshold'], dca_none['Net_benefit'], color='#475569', lw=1.6, linestyle=':',
            label='Treat None (Assume No LPD)')

th_vals = dca_model['Threshold'].values
nb_model_vals = dca_model['Net_benefit'].values
nb_all_vals = dca_all['Net_benefit'].values
nb_none_vals = dca_none['Net_benefit'].values

superior_mask = (nb_model_vals > nb_all_vals) & (nb_model_vals > nb_none_vals)
max_default = np.maximum(nb_all_vals, nb_none_vals)

ax_dca.fill_between(th_vals, max_default, nb_model_vals, where=superior_mask,
                    color='#BAE6FD', alpha=0.45, label='Higher apparent net benefit')

ax_dca.set_xlim(0.0, 1.0)
ax_dca.set_ylim(-0.05, 0.45)
ax_dca.set_xlabel('Threshold Probability ($P_t$)', fontsize=11, fontweight='bold')
ax_dca.set_ylabel('Net Benefit (NB)', fontsize=11, fontweight='bold')
ax_dca.set_title('Apparent decision-curve analysis', fontsize=13, fontweight='bold', loc='left', pad=10)
ax_dca.legend(loc='upper right', frameon=True, facecolor='#FFFFFF', edgecolor='#CBD5E1', fontsize=7.8)
ax_dca.grid(True, linestyle=':', alpha=0.5, color='#94A3B8')

# Output paths
fig_png = 'results/figures/Figure_Nomogram_Calibration_DCA_Composite.png'
fig_pdf = 'results/figures/Figure_Nomogram_Calibration_DCA_Composite.pdf'
plt.savefig(fig_png, dpi=300, bbox_inches='tight')
plt.savefig(fig_pdf, bbox_inches='tight')
plt.close()
print(f"Updated scientifically audited figure: {fig_png} and {fig_pdf}")

# Copy to conversation artifact directory
artifact_dir = r'C:\Users\asus\.gemini\antigravity-ide\brain\a3bd40f9-59d0-450e-883c-c58b2c00e17e'
if os.path.exists(artifact_dir):
    shutil.copy(fig_png, os.path.join(artifact_dir, 'Figure_Nomogram_Calibration_DCA_Composite.png'))
    print(f"Updated artifact in: {artifact_dir}")
