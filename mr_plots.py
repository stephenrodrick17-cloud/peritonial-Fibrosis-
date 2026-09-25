"""
mr_plots.py
===========
Visualization module for Two-Sample Mendelian Randomization (MR).
Generates publication-quality Forest plots, Scatter plots, Funnel plots,
and Leave-One-Out sensitivity plots.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Optional


# Style presets
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
FONT_FAMILY = 'DejaVu Sans'


def plot_cross_gene_forest(summary_df: pd.DataFrame, save_path: str, title: Optional[str] = None):
    """
    Generates a publication-grade Forest plot of IVW causal estimates across all hub genes.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Columns required: gene, ivw_beta, ivw_ci_lower, ivw_ci_upper, ivw_pval, n_snps, f_stat.
    save_path : str
        Filepath to save the PNG image.
    title : str, optional
        Plot title.
    """
    if summary_df.empty or "ivw_beta" not in summary_df.columns:
        return

    # Filter rows with valid estimates
    df_plot = summary_df.dropna(subset=["ivw_beta", "ivw_ci_lower", "ivw_ci_upper"]).copy()
    if df_plot.empty:
        return

    # Sort genes for display
    df_plot = df_plot.iloc[::-1].reset_index(drop=True)
    n_genes = len(df_plot)

    fig, ax = plt.subplots(figsize=(10, max(4, n_genes * 0.8 + 2)), facecolor="white")

    y_positions = np.arange(n_genes)

    # Plot point estimates and confidence intervals
    for i, row in df_plot.iterrows():
        b = row["ivw_beta"]
        ci_l = row["ivw_ci_lower"]
        ci_u = row["ivw_ci_upper"]
        p = row["ivw_pval"]
        nsnp = row.get("n_snps", "")
        fstat = row.get("f_stat", np.nan)

        color = "#2563EB" if p >= 0.05 else "#DC2626"

        ax.errorbar(
            b, i,
            xerr=[[b - ci_l], [ci_u - b]],
            fmt='s',
            color=color,
            ecolor=color,
            elinewidth=2.2,
            capsize=4.5,
            capthick=1.8,
            markersize=7.5,
            zorder=3
        )

        # Label text on the right
        f_str = f", F={fstat:0.1f}" if pd.notna(fstat) else ""
        label_text = f"β = {b:0.3f} [{ci_l:0.3f}, {ci_u:0.3f}], p = {p:0.3e} (SNPs={nsnp}{f_str})"
        ax.text(ax.get_xlim()[1] if ax.get_xlim()[1] > ci_u else ci_u * 1.05, i,
                f"  {label_text}", va="center", fontsize=9.5, color="#1E293B")

    # Reference line at null (beta = 0)
    ax.axvline(0, color="#64748B", linestyle="--", linewidth=1.4, zorder=2)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(df_plot["gene"], fontsize=11, fontweight="bold", color="#0F172A")
    ax.set_xlabel("Causal Effect Estimate (IVW Beta & 95% CI)", fontsize=11, fontweight="bold", labelpad=10)
    
    plot_title = title or "Two-Sample Mendelian Randomization: Hub Genes vs. Kidney/Peritoneal Outcome"
    ax.set_title(plot_title, fontsize=13, fontweight="bold", pad=15, color="#0F172A")

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_scatter(
    harmonized_df: pd.DataFrame,
    gene_name: str,
    ivw_res: Dict[str, float],
    egger_res: Optional[Dict[str, float]],
    wm_res: Optional[Dict[str, float]],
    save_path: str
):
    """
    Generates a scatter plot of SNP effects on exposure vs outcome with overlaid MR regression slopes.
    """
    if harmonized_df.empty:
        return

    fig, ax = plt.subplots(figsize=(8, 6.5), facecolor="white")

    x = harmonized_df["beta_exp"].values
    y = harmonized_df["beta_out"].values
    x_err = harmonized_df["se_exp"].values
    y_err = harmonized_df["se_out"].values

    # Scatter points with 2D error bars
    ax.errorbar(
        x, y, xerr=x_err, yerr=y_err,
        fmt='o', color="#475569", ecolor="#94A3B8", elinewidth=1.2,
        capsize=2.5, markersize=6, alpha=0.85, label="Instrument SNPs"
    )

    # Reference axes
    ax.axhline(0, color="#CBD5E1", linestyle=":", linewidth=1.0)
    ax.axvline(0, color="#CBD5E1", linestyle=":", linewidth=1.0)

    # Generate line coordinates
    x_min, x_max = np.min(x) * 0.9, np.max(x) * 1.1
    if x_min > 0:
        x_min = 0.0
    x_grid = np.linspace(x_min, x_max, 100)

    # 1. IVW line (forced through origin)
    if pd.notna(ivw_res.get("beta")):
        b_ivw = ivw_res["beta"]
        p_ivw = ivw_res["pval"]
        ax.plot(x_grid, b_ivw * x_grid, color="#2563EB", linewidth=2.2,
                label=f"IVW (β={b_ivw:0.3f}, p={p_ivw:0.3e})")

    # 2. MR-Egger line
    if egger_res and pd.notna(egger_res.get("beta")) and pd.notna(egger_res.get("intercept")):
        b_egg = egger_res["beta"]
        int_egg = egger_res["intercept"]
        p_egg = egger_res["pval"]
        p_int = egger_res["intercept_pval"]
        ax.plot(x_grid, int_egg + b_egg * x_grid, color="#DC2626", linestyle="--", linewidth=2.0,
                label=f"MR-Egger (β={b_egg:0.3f}, p={p_egg:0.3e} | int_p={p_int:0.3f})")

    # 3. Weighted Median line
    if wm_res and pd.notna(wm_res.get("beta")):
        b_wm = wm_res["beta"]
        p_wm = wm_res["pval"]
        ax.plot(x_grid, b_wm * x_grid, color="#059669", linestyle=":", linewidth=2.0,
                label=f"Weighted Median (β={b_wm:0.3f}, p={p_wm:0.3e})")

    ax.set_xlabel(f"SNP Effect on Exposure ({gene_name} eQTL Beta)", fontsize=11, fontweight="bold")
    ax.set_ylabel("SNP Effect on Outcome (GWAS Beta)", fontsize=11, fontweight="bold")
    ax.set_title(f"MR Scatter Plot: {gene_name} -> Outcome", fontsize=13, fontweight="bold", pad=15)
    ax.legend(loc="best", frameon=True, framealpha=0.9, facecolor="#F8FAFC")

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_funnel(harmonized_df: pd.DataFrame, gene_name: str, ivw_res: Dict[str, float], save_path: str):
    """
    Generates a funnel plot (Precision 1/SE vs Wald ratio estimate) to visually assess asymmetry/pleiotropy.
    """
    if len(harmonized_df) < 2:
        return

    wald_ratios = harmonized_df["beta_out"] / harmonized_df["beta_exp"]
    ratio_se = harmonized_df["se_out"] / np.abs(harmonized_df["beta_exp"])
    precision = 1.0 / ratio_se

    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")

    ax.scatter(wald_ratios, precision, color="#4F46E5", edgecolor="#312E81", s=50, alpha=0.85, label="Individual SNPs")

    # Vertical line at IVW estimate
    if pd.notna(ivw_res.get("beta")):
        b_ivw = ivw_res["beta"]
        ax.axvline(b_ivw, color="#DC2626", linestyle="--", linewidth=1.8, label=f"IVW Estimate (β={b_ivw:0.3f})")

        # Funnel triangular guideline
        max_prec = np.max(precision) * 1.15
        grid_prec = np.linspace(0.1, max_prec, 100)
        ci95_l = b_ivw - 1.96 / grid_prec
        ci95_u = b_ivw + 1.96 / grid_prec
        ax.plot(ci95_l, grid_prec, color="#94A3B8", linestyle=":", linewidth=1.2)
        ax.plot(ci95_u, grid_prec, color="#94A3B8", linestyle=":", linewidth=1.2, label="95% Pseudo-CI bounds")

    ax.set_xlabel("Wald Ratio Estimate (β_out / β_exp)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Instrument Precision (1 / SE)", fontsize=11, fontweight="bold")
    ax.set_title(f"MR Funnel Plot: {gene_name}", fontsize=13, fontweight="bold", pad=15)
    ax.legend(loc="upper right", frameon=True, framealpha=0.9, facecolor="#F8FAFC")

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_leave_one_out(loo_df: pd.DataFrame, gene_name: str, full_ivw_res: Dict[str, float], save_path: str):
    """
    Generates a Forest plot for Leave-One-Out sensitivity analysis.
    """
    if loo_df.empty:
        return

    n_rows = len(loo_df) + 1
    fig, ax = plt.subplots(figsize=(9, max(4, n_rows * 0.45 + 1.5)), facecolor="white")

    y_positions = np.arange(len(loo_df))

    # Plot individual leave-one-out points
    for i, row in loo_df.iterrows():
        b = row["beta"]
        ci_l = row["ci_lower"]
        ci_u = row["ci_upper"]

        ax.errorbar(
            b, i,
            xerr=[[b - ci_l], [ci_u - b]],
            fmt='o', color="#475569", ecolor="#94A3B8",
            elinewidth=1.8, capsize=3.5, markersize=5.5
        )

    # Plot overall reference at top
    b_all = full_ivw_res["beta"]
    ci_l_all = full_ivw_res["ci_lower"]
    ci_u_all = full_ivw_res["ci_upper"]
    top_y = len(loo_df)

    ax.errorbar(
        b_all, top_y,
        xerr=[[b_all - ci_l_all], [ci_u_all - b_all]],
        fmt='s', color="#DC2626", ecolor="#DC2626",
        elinewidth=2.2, capsize=4.5, markersize=7.5, label="All SNPs included"
    )

    ax.axvline(0, color="#64748B", linestyle="--", linewidth=1.2)
    ax.axvline(b_all, color="#DC2626", linestyle=":", linewidth=1.0, alpha=0.6)

    all_labels = loo_df["SNP_removed"].tolist() + ["All SNPs (Full Estimate)"]
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(all_labels, fontsize=9.5)
    ax.set_xlabel("Leave-One-Out IVW Estimate & 95% CI", fontsize=11, fontweight="bold")
    ax.set_title(f"Leave-One-Out Sensitivity Analysis: {gene_name}", fontsize=13, fontweight="bold", pad=15)

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
