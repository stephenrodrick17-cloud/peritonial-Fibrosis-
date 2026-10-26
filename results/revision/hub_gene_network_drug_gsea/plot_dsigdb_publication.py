"""
Plot top DSigDB candidate drugs associated with the diagnostic hub genes.
Shows top 10 candidates by nominal P-value and FDR, labeling overlapping genes,
and indicating the FDR < 0.05 threshold.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

out_dir = r"d:\Peritoneal Project\results\revision\hub_gene_network_drug_gsea"
art_dir = r"C:\Users\asus\.gemini\antigravity-ide\brain\a3bd40f9-59d0-450e-883c-c58b2c00e17e"

df = pd.read_csv(os.path.join(out_dir, "dsigdb_enrichment_full.csv"))
top10 = df.head(10).copy()
top10 = top10.sort_values("p_value", ascending=False)  # lowest p at the top of horizontal bar chart

fig, ax = plt.subplots(figsize=(10, 6.2), dpi=300)

y_pos = np.arange(len(top10))
neg_log10_p = -np.log10(top10["p_value"])

# Bar colors: highlight FDR < 0.05
colors = ["#C24C44" if sig else "#4A90E2" for sig in top10["significant_fdr_lt_0_05"]]

bars = ax.barh(y_pos, neg_log10_p, color=colors, height=0.68, edgecolor="none", alpha=0.9)

# Clean drug labels
clean_labels = []
for term in top10["dsigdb_term"]:
    # simplify label
    clean = term.replace(" CTD ", " (CTD: ").replace(" PC3 ", " (PC3: ").replace(" HL60 ", " (HL60: ")
    if "(" in clean:
        clean += ")"
    clean_labels.append(clean)

ax.set_yticks(y_pos)
ax.set_yticklabels(clean_labels, fontsize=10, fontweight="medium")
ax.set_xlabel(r"$-\log_{10}(P\mathrm{-value})$", fontsize=11, fontweight="bold")
ax.set_title("Candidate DSigDB Perturbagen Signatures Associated with Diagnostic Hub Genes",
             fontsize=12, fontweight="bold", pad=12)

# Add gene overlap annotations
for i, (_, row) in enumerate(top10.iterrows()):
    genes = row["overlapping_genes"]
    adj_p = row["adjusted_p_value"]
    val = neg_log10_p.iloc[i]
    sig_str = f" [Genes: {genes}; FDR = {adj_p:.3f}]" if row["significant_fdr_lt_0_05"] else f" [Genes: {genes}]"
    ax.text(val + 0.05, i, sig_str, va="center", fontsize=8.5, color="#333333")

ax.set_xlim(0, max(neg_log10_p) * 1.65)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", linestyle="--", alpha=0.4)

# Legend
custom_lines = [
    plt.Rectangle((0, 0), 1, 1, color="#C24C44"),
    plt.Rectangle((0, 0), 1, 1, color="#4A90E2")
]
ax.legend(custom_lines, ["FDR < 0.05 (Significant)", "Nominal P < 0.01 (Exploratory)"],
          loc="lower right", frameon=True, fontsize=9.5)

plt.tight_layout()

for path in [os.path.join(out_dir, "dsigdb_top_candidate_drugs.png"),
             os.path.join(art_dir, "Figure_B_DSigDB_Candidate_Drugs.png")]:
    plt.savefig(path, dpi=300, bbox_inches="tight")

for path in [os.path.join(out_dir, "dsigdb_top_candidate_drugs.pdf"),
             os.path.join(art_dir, "Figure_B_DSigDB_Candidate_Drugs.pdf")]:
    plt.savefig(path, dpi=300, bbox_inches="tight")

plt.close()
print("DSigDB publication figure updated!")
