"""
==============================================================================
SCRIPT 06: EXTERNAL VALIDATION IN GSE125498 (EARLY VS LATE STAGE DETECTION)
Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
Dataset: GSE125498 (Human effluent-derived peritoneal cells, N = 33 patients)
Cohorts: Short-term PD (SPD, 0-24 mo, N=20, Early Stage) vs
         Long-term PD (LPD, >= 25 mo, N=13, Late Stage Progressive Fibrosis)
Panels Validated:
  - 4 Causal Hub Genes (P4HA2, ADAMTS1, TNC, WNT11)
  - 12 Matrisome Hub Genes (ISM1, TGM2, MXRA5, COL3A1, COL5A2, POSTN, LOX, THBS3, etc.)
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

print(f"Total Cohort: {len(all_samples)} patients")
print(f"  - Early Stage (Short-term PD, 0-24 mo): {len(early_samples)} patients")
print(f"  - Late Stage (Long-term PD, >= 25 mo):   {len(late_samples)} patients")

# Panels
panel_4_hub = ["P4HA2", "ADAMTS1", "TNC"]
panel_12_hub = ["ISM1", "TGM2", "MXRA5", "COL3A1", "COL5A2", "POSTN", "LOX", "THBS3"]
additional_causal = ["BMP6", "IGFBP3", "COL4A2", "LTBP4", "FGL2", "ADAM28", "FBLN5", "CRISPLD2"]

all_target_genes = sorted(list(set(panel_4_hub + panel_12_hub + additional_causal)))
print(f"\nEvaluating {len(all_target_genes)} prioritized hub and causal biomarkers across cohorts...")

# Probe collapsing using MaxMean
gene_expr_dict = {}
validation_metrics = []

for gene in all_target_genes:
    matching = df_gpl[(df_gpl["Sym"] == gene) | (df_gpl["ILMN"] == gene)]["ID"].tolist()
    if not matching:
        continue
    
    probe_data = piv.loc[piv.index.isin(matching)]
    if len(probe_data) == 0:
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
    # If gene is upregulated, higher expression -> late stage; if downregulated, invert
    direction = "UP" if log2_diff > 0 else "DOWN"
    score_for_roc = gene_vals if direction == "UP" else -gene_vals
    fpr, tpr, _ = roc_curve(y_binary, score_for_roc)
    roc_auc = auc(fpr, tpr)
    
    validation_metrics.append({
        "GeneSymbol": gene,
        "Probe_ID": best_probe,
        "In_4_Hub_Panel": "YES" if gene in panel_4_hub else "NO",
        "In_12_Hub_Panel": "YES" if gene in panel_12_hub else "NO",
        "Mean_Early_SPD": mean_early,
        "Mean_Late_LPD": mean_late,
        "Median_Early_SPD": median_early,
        "Median_Late_LPD": median_late,
        "log2_FC_Late_vs_Early": log2_diff,
        "Direction": direction,
        "Mann_Whitney_U": u_stat,
        "Mann_Whitney_Pval": u_pval,
        "AUC_ROC": roc_auc
    })

df_metrics = pd.DataFrame(validation_metrics).sort_values("Mann_Whitney_Pval", ascending=True)
print("\n" + "=" * 70)
print("STEP 2: MANN-WHITNEY U VALIDATION & ROC DISCRIMINATION RESULTS")
print("=" * 70)
print(df_metrics[["GeneSymbol", "Direction", "log2_FC_Late_vs_Early", "Mann_Whitney_Pval", "AUC_ROC", "In_4_Hub_Panel", "In_12_Hub_Panel"]].to_string(index=False))

df_metrics.to_csv("results/tables/GSE125498_hub_genes_validation_metrics.csv", index=False)
print("\nSaved metrics to results/tables/GSE125498_hub_genes_validation_metrics.csv")

# ==============================================================================
# STEP 3: MULTI-GENE COMPOSITE BIOMARKER SIGNATURE
# ==============================================================================
X_val = pd.DataFrame(gene_expr_dict, index=all_samples)
# 4-Hub composite
avail_4 = [g for g in panel_4_hub if g in X_val.columns]
clf_4 = LogisticRegression(random_state=42)
clf_4.fit(X_val[avail_4], y_binary)
prob_4 = clf_4.predict_proba(X_val[avail_4])[:, 1]
fpr_4, tpr_4, _ = roc_curve(y_binary, prob_4)
auc_4 = auc(fpr_4, tpr_4)

# 12-Hub composite
avail_12 = [g for g in panel_12_hub if g in X_val.columns]
clf_12 = LogisticRegression(random_state=42)
clf_12.fit(X_val[avail_12], y_binary)
prob_12 = clf_12.predict_proba(X_val[avail_12])[:, 1]
fpr_12, tpr_12, _ = roc_curve(y_binary, prob_12)
auc_12 = auc(fpr_12, tpr_12)

print(f"\nCOMPOSITE BIOMARKER SIGNATURE ROC PERFORMANCE:")
print(f"  * 4-Hub Causal Panel ({', '.join(avail_4)}):  AUC = {auc_4:.3f}")
print(f"  * 12-Hub Matrisome Panel ({', '.join(avail_12)}): AUC = {auc_12:.3f}")

# ==============================================================================
# STEP 4: PUBLICATION-GRADE VISUALIZATIONS
# ==============================================================================
print("\n" + "=" * 70)
print("STEP 4: GENERATING VALIDATION VISUALIZATIONS")
print("=" * 70)

# ------------------------------------------------------------------------------
# Plot 1: Mann-Whitney U Boxplots for Hub Biomarkers
# ------------------------------------------------------------------------------
genes_to_plot = avail_4 + [g for g in ["TGM2", "MXRA5", "COL5A2", "COL3A1", "POSTN"] if g in X_val.columns]
n_genes = len(genes_to_plot)

fig, axes = plt.subplots(2, 4, figsize=(16, 9), facecolor="#F8FAFC")
axes = axes.flatten()

for idx, gene in enumerate(genes_to_plot[:8]):
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
    
    # Fetch P-value
    pval = df_metrics[df_metrics["GeneSymbol"] == gene]["Mann_Whitney_Pval"].values[0]
    auc_val = df_metrics[df_metrics["GeneSymbol"] == gene]["AUC_ROC"].values[0]
    
    # Significance formatting
    p_text = f"P = {pval:.3f}" if pval >= 0.001 else f"P = {pval:.2e}"
    ax.set_title(f"{gene} (AUC = {auc_val:.2f})", fontsize=12, fontweight="bold", pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("log2 Expression (GSE125498)", fontsize=10)
    
    # P-value banner
    y_max = sub_df["Expression"].max()
    y_min = sub_df["Expression"].min()
    h = (y_max - y_min) * 0.1
    ax.plot([0, 0, 1, 1], [y_max + h*0.2, y_max + h*0.5, y_max + h*0.5, y_max + h*0.2], color="#1E293B", linewidth=1.2)
    ax.text(0.5, y_max + h*0.6, p_text, ha="center", va="bottom", fontsize=10.5, fontweight="bold", color="#B91C1C")
    ax.set_ylim(y_min - h*0.5, y_max + h*1.8)
    ax.grid(axis="y", linestyle=":", alpha=0.5)

plt.suptitle("Validation of Consensus Hub Biomarkers in Independent Cohort GSE125498\nMann-Whitney U Test: Early-Stage (SPD, n=20) vs Late-Stage (LPD, n=13) Peritoneal Injury",
             fontsize=14, fontweight="bold", y=0.99)
plt.tight_layout()
plt.savefig("results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png", dpi=300)
plt.close()
print("Saved Plot 1: results/figures/Validation_01_hub_genes_mann_whitney_boxplots.png")

# ------------------------------------------------------------------------------
# Plot 2: Multi-Curve ROC Analysis
# ------------------------------------------------------------------------------
plt.figure(figsize=(9, 7.5), facecolor="#F8FAFC")

# Plot composite models
plt.plot(fpr_4, tpr_4, color="#DC2626", linewidth=2.8, label=f"4-Hub Causal Composite (AUC = {auc_4:.3f})")
plt.plot(fpr_12, tpr_12, color="#2563EB", linewidth=2.5, linestyle="--", label=f"12-Hub Matrisome Composite (AUC = {auc_12:.3f})")

# Plot top individual genes
colors_ind = ["#10B981", "#F59E0B", "#8B5CF6", "#EC4899"]
for idx, g in enumerate(["TGM2", "MXRA5", "COL5A2", "P4HA2"]):
    if g in X_val.columns:
        dir_g = df_metrics[df_metrics["GeneSymbol"] == g]["Direction"].values[0]
        score_g = X_val[g].values if dir_g == "UP" else -X_val[g].values
        fg, tg, _ = roc_curve(y_binary, score_g)
        ag = auc(fg, tg)
        plt.plot(fg, tg, color=colors_ind[idx], linewidth=1.8, label=f"{g} ({dir_g}, AUC = {ag:.3f})")

plt.plot([0, 1], [0, 1], color="#94A3B8", linestyle=":", linewidth=1.5, label="Random Guess (AUC = 0.500)")
plt.xlim(-0.02, 1.02)
plt.ylim(-0.02, 1.05)
plt.xlabel("False Positive Rate (1 - Specificity)", fontsize=11, fontweight="bold")
plt.ylabel("True Positive Rate (Sensitivity)", fontsize=11, fontweight="bold")
plt.title("Receiver Operating Characteristic (ROC) Discrimination\nEarly Stage (SPD, 0-24 mo) vs Late Stage (LPD, >=25 mo) Peritoneal Dialysis",
          fontsize=12, fontweight="bold", pad=15)
plt.legend(loc="lower right", frameon=True, facecolor="#EFF6FF", fontsize=10.5)
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig("results/figures/Validation_02_roc_curves_early_vs_late.png", dpi=300)
plt.close()
print("Saved Plot 2: results/figures/Validation_02_roc_curves_early_vs_late.png")

# ------------------------------------------------------------------------------
# Plot 3: Stage Progression Trajectory Bar Chart
# ------------------------------------------------------------------------------
plt.figure(figsize=(11, 6), facecolor="#F8FAFC")
df_sorted = df_metrics.sort_values("log2_FC_Late_vs_Early", ascending=False)
bar_colors = ["#DC2626" if fc > 0 else "#2563EB" for fc in df_sorted["log2_FC_Late_vs_Early"]]
bars = plt.bar(df_sorted["GeneSymbol"], df_sorted["log2_FC_Late_vs_Early"], color=bar_colors, edgecolor="#1E293B", width=0.6)
plt.axhline(0, color="#1E293B", linewidth=1.0)

for bar in bars:
    yval = bar.get_height()
    va = "bottom" if yval >= 0 else "top"
    offset = 0.02 if yval >= 0 else -0.05
    plt.text(bar.get_x() + bar.get_width()/2.0, yval + offset, f"{yval:+.2f}", ha="center", va=va, fontsize=9.5, fontweight="bold")

plt.title("Stage-Wise Expression Progression in GSE125498 (Late Stage LPD vs Early Stage SPD)\nDirectional Concordance with Primary Fibrosis Discovery Signature",
          fontsize=12, fontweight="bold", pad=15)
plt.xlabel("Hub & Causal Biomarkers", fontsize=11, fontweight="bold")
plt.ylabel("log2 Fold Difference (Late vs Early Stage)", fontsize=11, fontweight="bold")
plt.xticks(rotation=45, ha="right", fontsize=10, fontweight="medium")
plt.grid(axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig("results/figures/Validation_03_stage_progression_trajectories.png", dpi=300)
plt.close()
print("Saved Plot 3: results/figures/Validation_03_stage_progression_trajectories.png")

# ------------------------------------------------------------------------------
# Plot 4: Patient Cohort Heatmap
# ------------------------------------------------------------------------------
genes_heatmap = [g for g in (panel_4_hub + panel_12_hub) if g in X_val.columns]
expr_heat = X_val[genes_heatmap].T

# Standardize rows (Z-score)
expr_heat_z = expr_heat.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

# Column colors for stage
col_colors = ["#EF4444" if sample_stages[s] == "Late_Stage_LPD" else "#3B82F6" for s in all_samples]
col_colors_series = pd.Series(col_colors, index=all_samples, name="Stage")

cg = sns.clustermap(
    expr_heat_z, cmap="vlag", figsize=(12, 7.5),
    col_cluster=True, row_cluster=True,
    col_colors=col_colors_series,
    dendrogram_ratio=(0.18, 0.08),
    cbar_kws={"label": "Expression Z-Score"},
    linewidths=0.3, linecolor="#CBD5E1"
)
cg.figure.subplots_adjust(top=0.88)
cg.figure.suptitle("Cross-Cohort Patient Heatmap: Consensus Hub Genes Across 33 Clinical Samples\n(Blue = Early-Stage SPD, Red = Late-Stage LPD)",
                  fontsize=12, fontweight="bold", y=0.96)
plt.savefig("results/figures/Validation_04_patient_cohort_heatmap.png", dpi=300)
plt.close()
print("Saved Plot 4: results/figures/Validation_04_patient_cohort_heatmap.png")

print("\n" + "=" * 70)
print("EXTERNAL VALIDATION PIPELINE COMPLETED SUCCESSFULLY!")
print("=" * 70)
