"""
==============================================================================
SCRIPT 06: EXTERNAL VALIDATION OF WGCNA-ECM HUB GENES IN GSE125498
Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
Dataset: GSE125498 (Human effluent-derived peritoneal cells, N = 33 patients)
Cohorts: Short-term PD (SPD, 0-24 mo, N=20, Early Stage) vs
         Long-term PD (LPD, >= 25 mo, N=13, Late Stage Progressive Fibrosis)
Panel:   11 WGCNA-ECM Consensus Hub Genes (LASSO, SVM-RFE, RF, XGBoost)
==============================================================================
"""

import os
import GEOparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression

# Set random seed
np.random.seed(42)

# Ensure directories
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

print("=" * 70)
print("STEP 1: LOADING GSE125498 DATASET & CLINICAL STAGE ANNOTATIONS")
print("=" * 70)

# Load GSE125498 soft file
gse = GEOparse.get_GEO(filepath="data/GSE125498_family.soft.gz")
piv = gse.pivot_samples("VALUE")
gpl = list(gse.gpls.values())[0]

# Probe table for GPL10558
df_gpl = gpl.table[gpl.table["ID"].isin(piv.index)].copy()
df_gpl["Sym"] = df_gpl["Symbol"].fillna("").astype(str).str.strip().str.upper()
df_gpl["ILMN"] = df_gpl["ILMN_Gene"].fillna("").astype(str).str.strip().str.upper()

# Clinical Stage Mapping
sample_stages = {}
for name, gsm in gse.gsms.items():
    title = gsm.metadata.get("title", [""])[0]
    chars = str(gsm.metadata.get("characteristics_ch1", []))
    is_long = "long-term" in title.lower() or "long-term" in chars.lower()
    sample_stages[name] = "Late_Stage_LPD" if is_long else "Early_Stage_SPD"

early_samples = [s for s, st in sample_stages.items() if st == "Early_Stage_SPD"]
late_samples = [s for s, st in sample_stages.items() if st == "Late_Stage_LPD"]
all_samples = early_samples + late_samples
y_binary = np.array([1 if sample_stages[s] == "Late_Stage_LPD" else 0 for s in all_samples])

print(f"Total Cohort Size: {len(all_samples)} patients")
print(f"  - Early Stage (Short-term PD, 0-24 mo): {len(early_samples)} patients")
print(f"  - Late Stage (Long-term PD, >= 25 mo):   {len(late_samples)} patients")

# Load WGCNA-ECM Hub Genes from Task 3
wgcna_file = "results/tables/ML_hub_genes_from_WGCNA_ECM.csv"
if os.path.exists(wgcna_file):
    df_wgcna = pd.read_csv(wgcna_file)
    target_genes = df_wgcna["Gene_Symbol"].dropna().unique().tolist()
    print(f"Loaded {len(target_genes)} WGCNA-ECM Consensus Hub Genes from {wgcna_file}")
else:
    target_genes = ["ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"]

print(f"Target Hub Genes: {', '.join(target_genes)}")

# Probe collapsing using MaxMean rule
gene_expr_dict = {}
validation_metrics = []
missing_genes = []

for gene in target_genes:
    matching = df_gpl[(df_gpl["Sym"] == gene) | (df_gpl["ILMN"] == gene)]["ID"].tolist()
    if not matching:
        missing_genes.append(gene)
        continue
    
    probe_data = piv.loc[piv.index.isin(matching)]
    if len(probe_data) == 0:
        missing_genes.append(gene)
        continue
    
    # Best probe by MaxMean
    best_probe = probe_data.mean(axis=1).idxmax()
    vals = probe_data.loc[best_probe]
    
    val_early = vals[early_samples].astype(float)
    val_late = vals[late_samples].astype(float)
    
    # Store for matrix
    gene_expr_dict[gene] = vals[all_samples].astype(float).values
    
    # Statistical tests: Mann-Whitney U
    u_stat, u_pval = mannwhitneyu(val_late, val_early, alternative="two-sided")
    
    # Log2 difference (LPD - SPD)
    mean_early = val_early.mean()
    mean_late = val_late.mean()
    median_early = val_early.median()
    median_late = val_late.median()
    log2_diff = mean_late - mean_early
    
    # AUC-ROC for individual gene
    gene_vals = vals[all_samples].astype(float).values
    direction = "UP" if log2_diff > 0 else "DOWN"
    score_for_roc = gene_vals if direction == "UP" else -gene_vals
    fpr, tpr, _ = roc_curve(y_binary, score_for_roc)
    roc_auc = auc(fpr, tpr)
    
    validation_metrics.append({
        "Gene_Symbol": gene,
        "Probe_ID": best_probe,
        "Direction": direction,
        "log2_FC_Late_vs_Early": log2_diff,
        "Mean_Early_SPD": mean_early,
        "Mean_Late_LPD": mean_late,
        "Median_Early_SPD": median_early,
        "Median_Late_LPD": median_late,
        "Mann_Whitney_U": u_stat,
        "Mann_Whitney_Pval": u_pval,
        "AUC_ROC": roc_auc
    })

df_metrics = pd.DataFrame(validation_metrics).sort_values("Mann_Whitney_Pval", ascending=True)

print("\n" + "=" * 70)
print("STEP 2: MANN-WHITNEY U VALIDATION & ROC DISCRIMINATION RESULTS")
print("=" * 70)
print(df_metrics[["Gene_Symbol", "Direction", "log2_FC_Late_vs_Early", "Mann_Whitney_Pval", "AUC_ROC"]].to_string(index=False))

if missing_genes:
    print(f"\nNote: {len(missing_genes)} hub genes lacked mapped probes on GPL10558: {', '.join(missing_genes)}")

# Export validation table
out_table = "results/tables/GSE125498_wgcna_hub_validation_metrics.csv"
df_metrics.to_csv(out_table, index=False)
print(f"\nSaved validation metrics to: {out_table}")

# ==============================================================================
# STEP 3: MULTI-GENE COMPOSITE BIOMARKER SIGNATURE
# ==============================================================================
X_val = pd.DataFrame(gene_expr_dict, index=all_samples)
available_hubs = list(gene_expr_dict.keys())

clf = LogisticRegression(random_state=42)
clf.fit(X_val[available_hubs], y_binary)
prob = clf.predict_proba(X_val[available_hubs])[:, 1]
fpr_comp, tpr_comp, _ = roc_curve(y_binary, prob)
auc_comp = auc(fpr_comp, tpr_comp)

print("\n" + "=" * 70)
print(f"COMPOSITE BIOMARKER SIGNATURE ROC PERFORMANCE (N = {len(available_hubs)} Hub Genes):")
print(f"  * Available Hubs: {', '.join(available_hubs)}")
print(f"  * Composite Signature AUC = {auc_comp:.3f}")
print("=" * 70)

# ==============================================================================
# STEP 4: PUBLICATION-GRADE VISUALIZATIONS
# ==============================================================================
print("\n" + "=" * 70)
print("STEP 4: GENERATING VALIDATION VISUALIZATIONS")
print("=" * 70)

# ------------------------------------------------------------------------------
# Plot 1: Mann-Whitney U Boxplots for Hub Biomarkers
# ------------------------------------------------------------------------------
n_plot = len(available_hubs)
n_cols = min(4, n_plot)
n_rows = int(np.ceil(n_plot / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 4.5 * n_rows), facecolor="#F8FAFC")
if n_plot == 1:
    axes = np.array([axes])
axes = axes.flatten()

for idx, gene in enumerate(available_hubs):
    ax = axes[idx]
    sub_df = pd.DataFrame({
        "Expression": X_val[gene],
        "Stage": ["Late Stage\n(LPD >=25mo)" if sample_stages[s] == "Late_Stage_LPD" else "Early Stage\n(SPD 0-24mo)" for s in all_samples]
    })
    
    # Boxplot
    sns.boxplot(
        data=sub_df, x="Stage", y="Expression",
        palette=["#3B82F6", "#EF4444"], width=0.45,
        fliersize=0, ax=ax, boxprops=dict(alpha=0.75, edgecolor="#1E293B", linewidth=1.2)
    )
    # Jitter points
    sns.stripplot(
        data=sub_df, x="Stage", y="Expression",
        color="#0F172A", alpha=0.65, size=6, jitter=0.2, ax=ax
    )
    
    # P-value & AUC
    pval = df_metrics[df_metrics["Gene_Symbol"] == gene]["Mann_Whitney_Pval"].values[0]
    auc_val = df_metrics[df_metrics["Gene_Symbol"] == gene]["AUC_ROC"].values[0]
    p_text = f"P = {pval:.3f}" if pval >= 0.001 else f"P = {pval:.2e}"
    
    ax.set_title(f"{gene}\nAUC = {auc_val:.2f}", fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel("")
    ax.set_ylabel("log2 Expression (GSE125498)", fontsize=10)
    
    # P-value banner
    y_max = sub_df["Expression"].max()
    y_min = sub_df["Expression"].min()
    h = (y_max - y_min) * 0.1
    ax.plot([0, 0, 1, 1], [y_max + h*0.2, y_max + h*0.5, y_max + h*0.5, y_max + h*0.2], color="#1E293B", linewidth=1.2)
    ax.text(0.5, y_max + h*0.6, p_text, ha="center", va="bottom", fontsize=10, fontweight="bold", color="#B91C1C")
    ax.set_ylim(y_min - h*0.5, y_max + h*1.8)
    ax.grid(axis="y", linestyle=":", alpha=0.5)

# Hide any unused subplots
for extra_idx in range(n_plot, len(axes)):
    fig.delaxes(axes[extra_idx])

plt.suptitle("Validation of WGCNA-ECM Consensus Hub Biomarkers in Independent Cohort GSE125498\nMann-Whitney U Test: Early-Stage (SPD, n=20) vs Late-Stage (LPD, n=13) Peritoneal Injury",
             fontsize=13, fontweight="bold", y=0.99)
plt.tight_layout()
plot1_file = "results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png"
plt.savefig(plot1_file, dpi=300)
plt.close()
print(f"Saved Plot 1: {plot1_file}")

# ------------------------------------------------------------------------------
# Plot 2: Multi-Curve ROC Analysis
# ------------------------------------------------------------------------------
plt.figure(figsize=(8.5, 7), facecolor="#F8FAFC")

# Plot composite model
plt.plot(fpr_comp, tpr_comp, color="#059669", linewidth=3.2, label=f"WGCNA-ECM Hub Composite (AUC = {auc_comp:.3f})")

# Plot individual genes
colors_ind = ["#D97706", "#8B5CF6", "#EC4899", "#0284C7", "#E11D48", "#4F46E5", "#0D9488"]
for idx, g in enumerate(available_hubs):
    c = colors_ind[idx % len(colors_ind)]
    dir_g = df_metrics[df_metrics["Gene_Symbol"] == g]["Direction"].values[0]
    score_g = X_val[g].values if dir_g == "UP" else -X_val[g].values
    fg, tg, _ = roc_curve(y_binary, score_g)
    ag = auc(fg, tg)
    plt.plot(fg, tg, color=c, linewidth=1.6, linestyle="--", label=f"{g} ({dir_g}, AUC = {ag:.3f})")

plt.plot([0, 1], [0, 1], color="#94A3B8", linestyle=":", linewidth=1.5, label="Random Guess (AUC = 0.500)")
plt.xlim(-0.02, 1.02)
plt.ylim(-0.02, 1.05)
plt.xlabel("False Positive Rate (1 - Specificity)", fontsize=11, fontweight="bold")
plt.ylabel("True Positive Rate (Sensitivity)", fontsize=11, fontweight="bold")
plt.title(f"Receiver Operating Characteristic (ROC) Discrimination in GSE125498\nEarly Stage (SPD, 0-24 mo) vs Late Stage (LPD, >=25 mo) Peritoneal Dialysis",
          fontsize=12, fontweight="bold", pad=15)
plt.legend(loc="lower right", frameon=True, facecolor="#EFF6FF", fontsize=9.5)
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plot2_file = "results/figures/Validation_02_roc_curves_early_vs_late.png"
plt.savefig(plot2_file, dpi=300)
plt.close()
print(f"Saved Plot 2: {plot2_file}")

# ------------------------------------------------------------------------------
# Plot 3: Stage Progression Trajectory Bar Chart
# ------------------------------------------------------------------------------
plt.figure(figsize=(10, 6), facecolor="#F8FAFC")
df_sorted = df_metrics.sort_values("log2_FC_Late_vs_Early", ascending=False)
bar_colors = ["#DC2626" if fc > 0 else "#2563EB" for fc in df_sorted["log2_FC_Late_vs_Early"]]
bars = plt.bar(df_sorted["Gene_Symbol"], df_sorted["log2_FC_Late_vs_Early"], color=bar_colors, edgecolor="#1E293B", width=0.55)
plt.axhline(0, color="#1E293B", linewidth=1.0)

for bar in bars:
    yval = bar.get_height()
    va = "bottom" if yval >= 0 else "top"
    offset = 0.02 if yval >= 0 else -0.05
    plt.text(bar.get_x() + bar.get_width()/2.0, yval + offset, f"{yval:+.2f}", ha="center", va=va, fontsize=10, fontweight="bold")

plt.title("Stage-Wise Expression Progression in GSE125498 (Late Stage LPD vs Early Stage SPD)\nWGCNA-ECM Consensus Hub Biomarkers",
          fontsize=12, fontweight="bold", pad=15)
plt.xlabel("WGCNA-ECM Consensus Hub Biomarkers", fontsize=11, fontweight="bold")
plt.ylabel("log2 Fold Difference (Late vs Early Stage)", fontsize=11, fontweight="bold")
plt.xticks(fontsize=10.5, fontweight="medium")
plt.grid(axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plot3_file = "results/figures/Validation_03_stage_progression_trajectories.png"
plt.savefig(plot3_file, dpi=300)
plt.close()
print(f"Saved Plot 3: {plot3_file}")

# ------------------------------------------------------------------------------
# Plot 4: Patient Cohort Heatmap
# ------------------------------------------------------------------------------
expr_heat = X_val[available_hubs].T
expr_heat_z = expr_heat.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

col_colors = ["#EF4444" if sample_stages[s] == "Late_Stage_LPD" else "#3B82F6" for s in all_samples]
col_colors_series = pd.Series(col_colors, index=all_samples, name="Stage")

cg = sns.clustermap(
    expr_heat_z, cmap="vlag", figsize=(11, 7),
    col_cluster=True, row_cluster=True,
    col_colors=col_colors_series,
    dendrogram_ratio=(0.18, 0.10),
    cbar_kws={"label": "Expression Z-Score"},
    linewidths=0.3, linecolor="#CBD5E1"
)
cg.figure.subplots_adjust(top=0.88)
cg.figure.suptitle("Cross-Cohort Patient Heatmap: WGCNA-ECM Hub Genes Across 33 Clinical Samples\n(Blue = Early-Stage SPD, Red = Late-Stage LPD)",
                  fontsize=12, fontweight="bold", y=0.96)
plot4_file = "results/figures/Validation_04_patient_cohort_heatmap.png"
plt.savefig(plot4_file, dpi=300)
plt.close()
print(f"Saved Plot 4: {plot4_file}")

print("\n" + "=" * 70)
print("EXTERNAL VALIDATION COMPLETED SUCCESSFULLY FOR WGCNA-ECM HUB GENES!")
print("=" * 70)
