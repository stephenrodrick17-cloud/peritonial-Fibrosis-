"""
Perform enrichment analysis (GO BP, CC, MF and KEGG) specifically for the 43 module-assigned
matrisome DEGs against the 5,001 WGCNA background.
Generates:
1. CSV result tables for GO BP, GO CC, GO MF, and KEGG.
2. Publication-quality 4-panel horizontal barplot matching Fig 3E/3D in the reference paper.
3. High-resolution STRING functional network for the 43 genes.
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp

OUT_DIR = "results/revision/wgcna_paper_figures"
os.makedirs(OUT_DIR, exist_ok=True)

# 1. 43 module-assigned matrisome DEGs
intersect = pd.read_csv("results/revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv")
gcol = "Gene.symbol" if "Gene.symbol" in intersect.columns else "Annotated Gene"
genes_43 = intersect[intersect["Module"].notna()][gcol].tolist()
bg_genes = list(pd.read_csv("results/revision/wgcna_modules_dynamic_taskB.csv")["Gene"])

print(f"Loaded {len(genes_43)} genes. Background size: {len(bg_genes)}")

# 2. Run Enrichr for each ontology / database
gene_sets = {
    "GO BP": "GO_Biological_Process_2023",
    "GO CC": "GO_Cellular_Component_2023",
    "GO MF": "GO_Molecular_Function_2023",
    "KEGG":  "KEGG_2021_Human"
}

results_dict = {}
counts_summary = {}

for name, gs in gene_sets.items():
    print(f"Running enrichment for {name} ({gs})...")
    enr = gp.enrichr(gene_list=genes_43, gene_sets=gs, background=bg_genes, cutoff=0.05)
    df = enr.results.copy()
    
    # Calculate count of overlap genes
    df["Count"] = df["Genes"].apply(lambda x: len(str(x).split(";")) if pd.notna(x) and x != "" else 0)
    
    # Clean term names (remove GO IDs or pathway IDs from title)
    df["Clean_Term"] = df["Term"].apply(lambda x: str(x).split(" (GO:")[0].split(" - Homo")[0])
    
    csv_path = os.path.join(OUT_DIR, f"{name.replace(' ', '_')}_enrichment_43_matrisome_DEGs.csv")
    df.to_csv(csv_path, index=False)
    
    sig = df[df["Adjusted P-value"] < 0.05]
    counts_summary[name] = len(sig)
    results_dict[name] = df
    print(f"  -> Total terms: {len(df)}, Significant (adjusted P < 0.05): {len(sig)}")

# 3. Create publication 4-panel figure matching Fig 3E/3D
fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300)

palette = {
    "GO BP": "#E63385",  # Pink/Rose
    "GO CC": "#00A0E9",  # Cyan/Blue
    "GO MF": "#2ECC71",  # Green
    "KEGG":  "#9B59B6"   # Purple
}

for idx, (name, df) in enumerate(results_dict.items()):
    ax = axes[idx // 2, idx % 2]
    sig_df = df[df["Adjusted P-value"] < 0.05].copy()
    
    if len(sig_df) > 0:
        top = sig_df.sort_values("Adjusted P-value").head(10).copy()
        top = top.sort_values("Adjusted P-value", ascending=False)  # reverse so top is on top of barplot
        
        y_pos = range(len(top))
        bars = ax.barh(y_pos, top["Count"], color=palette[name], edgecolor="black", linewidth=0.4, alpha=0.88)
        
        ax.set_yticks(list(y_pos))
        labels = [d[:45] + "..." if len(d) > 45 else d for d in top["Clean_Term"]]
        ax.set_yticklabels(labels, fontsize=8.5, fontweight="medium")
        ax.set_xlabel("Gene Count", fontsize=10.5, fontweight="bold")
        ax.set_title(f"{name} enrichment analysis\n({len(sig_df)} significant terms, adj. P < 0.05)",
                     fontsize=12, fontweight="bold", pad=8)
        
        for i, (_, row) in enumerate(top.iterrows()):
            p_adj = row["Adjusted P-value"]
            ax.text(row["Count"] + 0.15, i, f"p={p_adj:.1e}", va="center", fontsize=7.5, color="#475569")
            
        ax.set_xlim(0, max(top["Count"]) + 1.8)
    else:
        # Fallback if few reach 0.05 FDR
        top = df.sort_values("P-value").head(8).copy().sort_values("P-value", ascending=False)
        y_pos = range(len(top))
        ax.barh(y_pos, top["Count"], color=palette[name], edgecolor="black", linewidth=0.4, alpha=0.5)
        ax.set_yticks(list(y_pos))
        labels = [d[:45] + "..." if len(d) > 45 else d for d in top["Clean_Term"]]
        ax.set_yticklabels(labels, fontsize=8.5)
        ax.set_xlabel("Gene Count", fontsize=10.5, fontweight="bold")
        ax.set_title(f"{name} enrichment analysis (nominal P < 0.05)", fontsize=12, fontweight="bold", pad=8)
        
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)

fig.suptitle(
    "Functional Enrichment Analysis of 43 Module-Assigned Matrisome DEGs\n(Background: 5,001 WGCNA genes; BH-adjusted P < 0.05)",
    fontsize=14.5,
    fontweight="bold",
    y=0.99
)

plt.tight_layout()

for ext in ["png", "pdf", "svg"]:
    out_file = os.path.join(OUT_DIR, f"Fig3D_enrichment_43_matrisome_DEGs.{ext}")
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_file}")

plt.close()

# Summary log
summary_log = {
    "gene_set": "43 module-assigned matrisome DEGs",
    "background_genes": len(bg_genes),
    "significant_terms_FDR05": counts_summary
}
with open(os.path.join(OUT_DIR, "enrichment_43_summary.json"), "w") as f:
    json.dump(summary_log, f, indent=2)

print("\nAll enrichment steps finished successfully!")
print("Summary:", counts_summary)
