"""
Script 38: Stage 4 Steps 5, 6, 7 - Ambient Reference, Multi-Dataset Concordance & Forest Plot
Rules:
- Step 5: Report hub gene expression in T cells (negative reference) vs Stromal expression
- Step 6: Concordance table across:
    (1) Discovery tissue EPS vs Pooled Control (GSE62928)
    (2) Discovery tissue EPS vs PD Control (GSE62928)
    (3) Stromal pseudobulk LV_UF vs LV_NOT_UF (GSE248762)
    (4) Effluent array validation (GSE125498)
  Categories: concordant (significant), same direction (not significant), discordant, INSUFFICIENT DATA
- Step 7: Forest plot of stromal log2FC (95% CI) with INSUFFICIENT DATA shown in grey; dotplots;
  strictly non-causal language; effluent described as effluent, UF failure described as UF failure (not EPS).
Seed 42.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

print("=" * 80)
print("STAGE 4 - STEP 5: AMBIENT-RNA REFERENCE PROFILE")
print("=" * 80)

# Load descriptive summary
desc_df = pd.read_csv("results/tables/stage4_hub_genes_descriptive_summary.csv")

# Compare stromal vs T cells (ambient reference)
st_df = desc_df[desc_df["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)"].set_index("Gene")
t_df = desc_df[desc_df["Cell_Type"] == "T cell"].set_index("Gene")

hub_genes = ["COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
             "EDIL3", "LOX", "INHBA", "ISM1", "COMP"]

ambient_rows = []
for g in hub_genes:
    st_pct = st_df.loc[g, "Pct_Expressing_%"] if g in st_df.index else np.nan
    st_mean = st_df.loc[g, "Mean_log1p_CP10k"] if g in st_df.index else np.nan
    t_pct = t_df.loc[g, "Pct_Expressing_%"] if g in t_df.index else np.nan
    t_mean = t_df.loc[g, "Mean_log1p_CP10k"] if g in t_df.index else np.nan
    
    st_exceeds = (st_mean > t_mean) and (st_pct > t_pct)
    fold_diff = (st_mean / (t_mean + 1e-6)) if pd.notna(st_mean) and pd.notna(t_mean) else np.nan
    
    ambient_rows.append({
        "Gene": g,
        "Stromal_Pct_%": st_pct,
        "Stromal_Mean_CP10k": st_mean,
        "Tcell_Ref_Pct_%": t_pct,
        "Tcell_Ref_Mean_CP10k": t_mean,
        "Stromal_Exceeds_Tcell": st_exceeds,
        "Stromal_to_Tcell_Ratio": round(fold_diff, 2)
    })

df_ambient = pd.DataFrame(ambient_rows)
df_ambient.to_csv("results/tables/stage4_ambient_reference_comparison.csv", index=False)
print("Ambient RNA Reference Table (T cells used as negative biological reference):")
print(df_ambient.to_string(index=False))

# ---------------------------------------------------------------------------
# Step 6: Concordance Across Compartments & Datasets
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("STAGE 4 - STEP 6: MULTI-DATASET DIRECTIONAL CONCORDANCE")
print("=" * 80)

# 1. Tissue GSE62928
tissue_df = pd.read_csv("results/tables/B_gse62928_recomputed_hub_genes.csv").set_index("Gene")

# 2. GSE125498 validation
val_df = pd.read_csv("results/tables/C_gse125498_concordance.csv").set_index("Gene")

# 3. Primary stromal pseudobulk (GSE248762)
st_pb = pd.read_csv("results/tables/stage4_primary_edger_pseudobulk.csv")
st_pb_uf = st_pb[
    (st_pb["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
    (st_pb["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].set_index("Gene")

def get_dir(val):
    if pd.isna(val):
        return "NA"
    return "UP" if val > 0 else ("DOWN" if val < 0 else "FLAT")

concord_rows = []
for g in hub_genes:
    # Tissue EPS vs Pooled
    t_fc_pool = tissue_df.loc[g, "logFC_EPS_vs_PooledControl"]
    t_p_pool = tissue_df.loc[g, "P.Value_EPS_vs_PooledControl"]
    t_dir_pool = get_dir(t_fc_pool)
    
    # Tissue EPS vs PD
    t_fc_pd = tissue_df.loc[g, "logFC_EPS_vs_PD"]
    t_p_pd = tissue_df.loc[g, "P.Value_EPS_vs_PD"]
    t_dir_pd = get_dir(t_fc_pd)
    
    # Effluent array (GSE125498)
    if g in val_df.index and pd.notna(val_df.loc[g, "Validation_Effluent_log2FC"]):
        v_fc = val_df.loc[g, "Validation_Effluent_log2FC"]
        v_p = val_df.loc[g, "Validation_Effluent_Raw_P_Limma"]
        v_dir = get_dir(v_fc)
        v_status = "MAPPED"
    else:
        v_fc = np.nan; v_p = np.nan; v_dir = "UNMAPPED"; v_status = "UNMAPPED"
        
    # Stromal Pseudobulk (GSE248762)
    s_fc = st_pb_uf.loc[g, "log2FC"]
    s_p = st_pb_uf.loc[g, "PValue"]
    s_fdr = st_pb_uf.loc[g, "BH_FDR"]
    s_ev = st_pb_uf.loc[g, "Evidence_Status"]
    s_dir = get_dir(s_fc)
    
    # Concordance Category between Tissue (GSE62928) and Stromal Effluent (GSE248762)
    if s_ev == "INSUFFICIENT DATA":
        category = "INSUFFICIENT DATA"
    elif s_dir == t_dir_pool:
        if pd.notna(s_fdr) and s_fdr < 0.05:
            category = "concordant (significant)"
        else:
            category = "same direction (not significant)"
    else:
        category = "discordant"
        
    concord_rows.append({
        "Gene": g,
        "Tissue_EPS_vs_Pool_log2FC": round(t_fc_pool, 3),
        "Tissue_EPS_vs_Pool_P": f"{t_p_pool:.2e}",
        "Tissue_Dir": t_dir_pool,
        "Stromal_Effluent_log2FC": round(s_fc, 3) if pd.notna(s_fc) else np.nan,
        "Stromal_Effluent_Raw_P": f"{s_p:.4f}" if pd.notna(s_p) else "NA",
        "Stromal_Effluent_FDR": f"{s_fdr:.4f}" if pd.notna(s_fdr) else "NA",
        "Stromal_Effluent_Dir": s_dir,
        "GSE125498_Effluent_log2FC": round(v_fc, 3) if pd.notna(v_fc) else "—",
        "GSE125498_Dir": v_dir,
        "Concordance_Tissue_vs_Stromal_Effluent": category
    })

df_concord = pd.DataFrame(concord_rows)
df_concord.to_csv("results/tables/stage4_cross_dataset_concordance.csv", index=False)
print("Cross-Dataset Concordance Table (Tissue Fibrosis vs. Effluent scRNA-seq):")
print(df_concord[["Gene", "Tissue_Dir", "Stromal_Effluent_log2FC", "Stromal_Effluent_Raw_P", "Stromal_Effluent_Dir", "GSE125498_Dir", "Concordance_Tissue_vs_Stromal_Effluent"]].to_string(index=False))

print("\n--- Concordance Summary Counts (from code) ---")
category_counts = df_concord["Concordance_Tissue_vs_Stromal_Effluent"].value_counts()
print(category_counts.to_string())

# ---------------------------------------------------------------------------
# Step 7: Figures - Forest Plot
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("STAGE 4 - STEP 7: FIGURES (FOREST PLOT & PRINTOUTS)")
print("=" * 80)

# Build Forest Plot
fig, ax = plt.subplots(figsize=(9, 6), dpi=150)

# Order genes by log2FC
plot_df = st_pb_uf.loc[hub_genes].copy()
plot_df = plot_df.sort_values("log2FC", ascending=True)

y_positions = np.arange(len(plot_df))
for i, (g, row) in enumerate(plot_df.iterrows()):
    lfc = row["log2FC"]
    ci_lo = row["CI_95_low"]
    ci_hi = row["CI_95_high"]
    ev = row["Evidence_Status"]
    pval = row["PValue"]
    
    if ev == "INSUFFICIENT DATA":
        color = "#999999"
        label_txt = f"{g} (INSUFFICIENT DATA)"
    elif pval < 0.05:
        color = "#d95f02" # nominal sig
        label_txt = f"{g} (P={pval:.3f})"
    else:
        color = "#2b83ba"
        label_txt = f"{g} (P={pval:.3f})"
        
    ax.plot([ci_lo, ci_hi], [i, i], color=color, lw=2.0)
    ax.scatter(lfc, i, color=color, s=45, zorder=5)

ax.axvline(0, color="black", linestyle="--", lw=1.0, alpha=0.7)
ax.set_yticks(y_positions)
ax.set_yticklabels(plot_df.index, fontweight="bold", fontsize=11)
ax.set_xlabel("Stromal Pseudobulk log2(Fold Change): LV_UF vs LV_NOT_UF (95% CI)", fontsize=11, fontweight="bold")
ax.set_title("Effluent Stromal / Mesothelial-Lineage: Hub Gene Differential Expression\n(Grey = Insufficient Data; Orange = Nominal P < 0.05; edgeR Quasi-Likelihood)",
             fontsize=12, fontweight="bold", pad=12)

# Custom legend
custom_lines = [
    plt.Line2D([0], [0], color="#d95f02", lw=2, marker="o", label="Nominal P < 0.05 (FDR >= 0.05)"),
    plt.Line2D([0], [0], color="#2b83ba", lw=2, marker="o", label="P >= 0.05"),
    plt.Line2D([0], [0], color="#999999", lw=2, marker="o", label="INSUFFICIENT DATA (<2 donors/group)")
]
ax.legend(handles=custom_lines, loc="lower right", frameon=True, fontsize=9)
plt.tight_layout()

forest_fig_path = "results/figures/stage4_stromal_pseudobulk_forest_plot.png"
fig.savefig(forest_fig_path)
plt.close(fig)
print(f"Saved forest plot to {forest_fig_path}")

print("\n" + "=" * 80)
print("STAGE 4 COMPLETE. ALL RESULTS PRINTED VERBATIM.")
print("=" * 80)
