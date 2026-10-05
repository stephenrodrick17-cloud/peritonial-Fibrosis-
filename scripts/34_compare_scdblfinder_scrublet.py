"""
Script 34: Compare scDblFinder (R) vs Scrublet (fixed rule) sensitivity check
Calculates per-sample doublet calls, overlap (Jaccard index), and stromal-cluster doublet fraction.
"""
import pandas as pd
import numpy as np
import anndata as ad

def run_comparison():
    print("=" * 80)
    print("3. DOUBLETS SENSITIVITY CHECK: scDblFinder vs Scrublet")
    print("=" * 80)
    
    # Load Scrublet annotations from full QC object
    adata_qc = ad.read_h5ad("data/processed/GSE248762_hubblind_allcells_qc.h5ad")
    obs = adata_qc.obs.copy()
    obs["scrublet_doublet"] = (obs["doublet_score"] > 0.4390)
    
    # Load scDblFinder results
    df_scd = pd.read_csv("results/tables/scDblFinder_per_barcode_calls.csv")
    df_scd["is_doublet_scDblFinder"] = (df_scd["scDblFinder_class"] == "doublet")
    
    # Merge on barcode index
    obs["barcode_full"] = obs.index
    merged = pd.merge(obs, df_scd[["barcode", "scDblFinder_score", "is_doublet_scDblFinder"]],
                      left_on="barcode_full", right_on="barcode", how="inner")
    
    print(f"Total merged barcodes: {len(merged)}")
    assert len(merged) == 118895, f"Merged count {len(merged)} != 118895"
    
    # Load cell type annotations from lite h5ad
    adata_lite = ad.read_h5ad("data/processed/GSE248762_harmony_annotated_obs.h5ad")
    if "_index" in adata_lite.obs.columns:
        cell_type_map = dict(zip(adata_lite.obs["_index"], adata_lite.obs["cell_type"]))
    else:
        cell_type_map = adata_lite.obs["cell_type"].to_dict()
    merged["cell_type"] = merged["barcode_full"].map(cell_type_map)
    
    # Per-sample comparison
    sample_stats = []
    for gsm, g_df in merged.groupby("gsm", sort=False):
        donor = g_df["donor_id"].iloc[0]
        grp = g_df["group"].iloc[0]
        n = len(g_df)
        
        scrub_n = int(g_df["scrublet_doublet"].sum())
        scd_n = int(g_df["is_doublet_scDblFinder"].sum())
        
        both = int((g_df["scrublet_doublet"] & g_df["is_doublet_scDblFinder"]).sum())
        either = int((g_df["scrublet_doublet"] | g_df["is_doublet_scDblFinder"]).sum())
        jaccard = (both / either) if either > 0 else 0.0
        
        sample_stats.append({
            "GSM": gsm,
            "Sample": donor,
            "Group": grp,
            "Barcodes": n,
            "Scrublet_Calls": scrub_n,
            "Scrublet_Rate_%": round(scrub_n / n * 100, 2),
            "scDblFinder_Calls": scd_n,
            "scDblFinder_Rate_%": round(scd_n / n * 100, 2),
            "Both_Flagged": both,
            "Either_Flagged": either,
            "Jaccard_Index": round(jaccard, 3)
        })
    
    df_comp = pd.DataFrame(sample_stats)
    
    # Summary TOTAL row
    tot_barcodes = df_comp["Barcodes"].sum()
    tot_scrub = df_comp["Scrublet_Calls"].sum()
    tot_scd = df_comp["scDblFinder_Calls"].sum()
    tot_both = int((merged["scrublet_doublet"] & merged["is_doublet_scDblFinder"]).sum())
    tot_either = int((merged["scrublet_doublet"] | merged["is_doublet_scDblFinder"]).sum())
    tot_jaccard = (tot_both / tot_either) if tot_either > 0 else 0.0
    
    total_row = {
        "GSM": "—",
        "Sample": "TOTAL",
        "Group": "—",
        "Barcodes": tot_barcodes,
        "Scrublet_Calls": tot_scrub,
        "Scrublet_Rate_%": round(tot_scrub / tot_barcodes * 100, 2),
        "scDblFinder_Calls": tot_scd,
        "scDblFinder_Rate_%": round(tot_scd / tot_barcodes * 100, 2),
        "Both_Flagged": tot_both,
        "Either_Flagged": tot_either,
        "Jaccard_Index": round(tot_jaccard, 3)
    }
    
    df_comp_total = pd.concat([df_comp, pd.DataFrame([total_row])], ignore_index=True)
    
    # Assertions on TOTAL
    assert df_comp["Barcodes"].sum() == df_comp_total.loc[df_comp_total["Sample"]=="TOTAL", "Barcodes"].values[0]
    assert df_comp["Scrublet_Calls"].sum() == df_comp_total.loc[df_comp_total["Sample"]=="TOTAL", "Scrublet_Calls"].values[0]
    assert df_comp["scDblFinder_Calls"].sum() == df_comp_total.loc[df_comp_total["Sample"]=="TOTAL", "scDblFinder_Calls"].values[0]
    
    df_comp_total.to_csv("results/tables/doublet_scrublet_vs_scdblfinder_comparison.csv", index=False)
    print("\n--- Per-Sample Doublet Calls & Overlap (Jaccard) ---")
    print(df_comp_total.to_string(index=False))
    
    # -----------------------------------------------------------------------
    # Stromal cluster doublet fraction
    # -----------------------------------------------------------------------
    print("\n--- Stromal Cluster: Doublet Rate Sensitivity ---")
    st_mask = (merged["cell_type"] == "stromal / mesothelial-lineage (unresolved)")
    st_df = merged[st_mask]
    n_st = len(st_df)
    
    st_scrub = int(st_df["scrublet_doublet"].sum())
    st_scd = int(st_df["is_doublet_scDblFinder"].sum())
    st_both = int((st_df["scrublet_doublet"] & st_df["is_doublet_scDblFinder"]).sum())
    st_either = int((st_df["scrublet_doublet"] | st_df["is_doublet_scDblFinder"]).sum())
    
    st_summary = pd.DataFrame([
        {"Method": "Scrublet (fixed threshold 0.4390)", "Doublets_in_Stromal": st_scrub, "Total_Stromal": n_st, "Fraction_%": round(st_scrub / n_st * 100, 2)},
        {"Method": "scDblFinder (R independent check)", "Doublets_in_Stromal": st_scd, "Total_Stromal": n_st, "Fraction_%": round(st_scd / n_st * 100, 2)},
        {"Method": "Both flagged (Scrublet & scDblFinder)", "Doublets_in_Stromal": st_both, "Total_Stromal": n_st, "Fraction_%": round(st_both / n_st * 100, 2)},
        {"Method": "Either flagged (Scrublet OR scDblFinder)", "Doublets_in_Stromal": st_either, "Total_Stromal": n_st, "Fraction_%": round(st_either / n_st * 100, 2)}
    ])
    print(st_summary.to_string(index=False))
    st_summary.to_csv("results/tables/stromal_doublet_sensitivity.csv", index=False)

if __name__ == "__main__":
    run_comparison()
