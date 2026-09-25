"""
run_mr_pipeline.py
==================
Main execution script for Two-Sample Mendelian Randomization (MR) Pipeline.

Evaluates potential causal relationships between transcriptomic hub genes
(exposure: eQTLGen blood cis-eQTLs) and kidney/peritoneal function outcomes (e.g. CKDGen eGFR).

Supports OpenGWAS API retrieval with server-side LD clumping as well as
local CSV files (exposure_[GENE].csv and outcome.csv).
"""

import os
import sys
import pandas as pd
import numpy as np
from typing import Dict, List, Optional

# Import internal modules
from data_fetch import (
    get_instruments_opengwas,
    get_outcome_associations,
    clump_snps_opengwas,
    search_studies
)
from mr_stats import (
    harmonize_data,
    calculate_f_statistic,
    perform_ivw,
    perform_mr_egger,
    perform_weighted_median,
    calculate_cochrans_q,
    perform_leave_one_out,
    flag_outliers_mr_presso_style
)
from mr_plots import (
    plot_cross_gene_forest,
    plot_scatter,
    plot_funnel,
    plot_leave_one_out
)


# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# OpenGWAS Study IDs for Hub Genes (eQTLGen batch: 'eqtl-a-...')
# You can search using `search_studies("SOCS1", batch="eqtl-a")`
HUB_GENE_EXPOSURE_IDS: Dict[str, str] = {
    # Example format: "GENE_SYMBOL": "OPENGWAS_ID"
    # "SOCS1": "eqtl-a-ENSG00000185338",
    # "PIM2": "eqtl-a-ENSG00000102076",
    # "HSH2D": "eqtl-a-ENSG00000154864",
    # "MYO3B": "eqtl-a-ENSG00000144583",
}

# Outcome GWAS Study ID (e.g., CKDGen eGFR / Bio-PD peritoneal ultrafiltration)
# Example: "ebi-a-GCST008064" (Wuttke et al. 2019 eGFR) or custom OpenGWAS ID
OUTCOME_ID: str = "ebi-a-GCST008064"

# Significance thresholds for instrument selection
PVAL_PRIMARY: float = 5e-8      # Standard genome-wide significance
PVAL_FALLBACK: float = 1e-5     # Relaxed fallback if < 3 instruments pass primary threshold

# LD Clumping parameters (1000G EUR reference)
CLUMP_R2: float = 0.001
CLUMP_KB: int = 10000

# Directory settings
RESULTS_DIR: str = "mr_results"

# Local file fallback directory (searched if OpenGWAS IDs are empty)
LOCAL_DATA_DIR: str = "."
LOCAL_OUTCOME_FILE: str = "outcome.csv"


# =============================================================================
# PIPELINE EXECUTION ENGINE
# =============================================================================

def run_pipeline_for_gene(
    gene_symbol: str,
    exposure_id_or_file: str,
    outcome_id_or_file: str,
    out_dir: str,
    token: Optional[str] = None
) -> Optional[Dict]:
    """
    Executes the complete Mendelian Randomization workflow for a single hub gene.
    """
    print(f"\n{'='*70}")
    print(f" PROCESSING HUB GENE: {gene_symbol}")
    print(f"{'='*70}")
    
    gene_out_dir = os.path.join(out_dir, gene_symbol)
    os.makedirs(gene_out_dir, exist_ok=True)
    
    # -------------------------------------------------------------------------
    # STEP 1: Instrument Selection & Fetching
    # -------------------------------------------------------------------------
    is_local_exp = os.path.isfile(exposure_id_or_file)
    
    if is_local_exp:
        print(f"[{gene_symbol}] Loading local exposure file: {exposure_id_or_file}")
        exp_df = pd.read_csv(exposure_id_or_file)
        
        # Filter & Clump local data
        funnel_exp = {"raw_instruments": len(exp_df), "threshold_used": f"{PVAL_PRIMARY:0.1e}"}
        
        # Primary filter
        mask = exp_df["pval"] < PVAL_PRIMARY
        if mask.sum() < 3:
            print(f"[WARNING] [{gene_symbol}] Relaxing p-value threshold to {PVAL_FALLBACK:0.1e}")
            mask = exp_df["pval"] < PVAL_FALLBACK
            funnel_exp["threshold_used"] = f"{PVAL_FALLBACK:0.1e} (relaxed)"
            
        exp_filtered = exp_df[mask].copy()
        funnel_exp["pval_filtered"] = len(exp_filtered)
        
        # LD clumping
        clumped_snps = clump_snps_opengwas(
            rsids=exp_filtered["SNP"].tolist(),
            pvals=exp_filtered["pval"].tolist(),
            r2_threshold=CLUMP_R2,
            kb_window=CLUMP_KB,
            token=token
        )
        exp_instruments = exp_filtered[exp_filtered["SNP"].isin(clumped_snps)].drop_duplicates(subset=["SNP"]).copy()
        funnel_exp["clumped_independent"] = len(exp_instruments)
    else:
        print(f"[{gene_symbol}] Querying OpenGWAS API for study ID: {exposure_id_or_file}")
        exp_instruments, funnel_exp = get_instruments_opengwas(
            gwas_id=exposure_id_or_file,
            pval_primary=PVAL_PRIMARY,
            pval_fallback=PVAL_FALLBACK,
            r2=CLUMP_R2,
            kb=CLUMP_KB,
            token=token
        )
        
    if exp_instruments.empty:
        print(f"[ERROR] [{gene_symbol}] No valid instruments available. Skipping.")
        return None
        
    snp_list = exp_instruments["SNP"].tolist()
    print(f"[{gene_symbol}] {len(snp_list)} independent instrument SNPs ready.")

    # -------------------------------------------------------------------------
    # STEP 2: Outcome Data Retrieval
    # -------------------------------------------------------------------------
    is_local_out = os.path.isfile(outcome_id_or_file)
    if is_local_out:
        print(f"[{gene_symbol}] Loading local outcome file: {outcome_id_or_file}")
        out_raw = pd.read_csv(outcome_id_or_file)
        out_df = out_raw[out_raw["SNP"].isin(snp_list)].copy()
    else:
        print(f"[{gene_symbol}] Querying outcome associations from OpenGWAS ID: {outcome_id_or_file}")
        out_df = get_outcome_associations(outcome_gwas_id=outcome_id_or_file, snps=snp_list, token=token)
        
    if out_df.empty:
        print(f"[ERROR] [{gene_symbol}] None of the {len(snp_list)} instrument SNPs found in outcome GWAS. Skipping.")
        return None

    # -------------------------------------------------------------------------
    # STEP 3: Harmonization
    # -------------------------------------------------------------------------
    harm_df, harm_funnel = harmonize_data(exp_instruments, out_df, palindrome_threshold=0.08)
    
    # Print intermediate filtering funnel
    print(f"\n--- [{gene_symbol}] INSTRUMENT FILTERING FUNNEL ---")
    print(f"  1. Raw Exposure SNPs Available:      {funnel_exp.get('raw_instruments', 'N/A')}")
    print(f"  2. After P-value Cutoff ({funnel_exp.get('threshold_used', '')}): {funnel_exp.get('pval_filtered', 'N/A')}")
    print(f"  3. After LD Clumping (r² < {CLUMP_R2}):   {funnel_exp.get('clumped_independent', 'N/A')}")
    print(f"  4. Matched in Outcome GWAS:          {harm_funnel['matched_on_id']}")
    print(f"  5. Dropped (Ambiguous Palindromes):  {harm_funnel['dropped_ambiguous_palindromes']}")
    print(f"  6. Dropped (Incompatible Alleles):   {harm_funnel['dropped_incompatible_alleles']}")
    print(f"  => FINAL HARMONIZED INSTRUMENTS:     {harm_funnel['harmonized_final']}")
    print(f"----------------------------------------------------\n")
    
    if harm_df.empty or len(harm_df) < 1:
        print(f"[ERROR] [{gene_symbol}] Zero SNPs survived harmonization.")
        return None
        
    # Save harmonized data table
    harm_csv_path = os.path.join(gene_out_dir, "harmonized_data.csv")
    harm_df.to_csv(harm_csv_path, index=False)
    
    # -------------------------------------------------------------------------
    # STEP 4: Instrument Strength (F-statistic)
    # -------------------------------------------------------------------------
    f_series, mean_f, is_weak = calculate_f_statistic(harm_df)
    harm_df["F_stat"] = f_series
    print(f"[{gene_symbol}] Mean F-Statistic: {mean_f:0.2f} " + ("[WEAK INSTRUMENT WARNING (Mean F < 10)]" if is_weak else "[STRONG INSTRUMENT]"))
    
    # -------------------------------------------------------------------------
    # STEP 5: Core Mendelian Randomization Estimation
    # -------------------------------------------------------------------------
    ivw_res = perform_ivw(harm_df)
    egger_res = perform_mr_egger(harm_df) if len(harm_df) >= 3 else None
    wm_res = perform_weighted_median(harm_df, n_boot=1000) if len(harm_df) >= 3 else None
    
    # Heterogeneity
    q_res = calculate_cochrans_q(harm_df, ivw_res["beta"])
    
    # Leave-One-Out
    loo_df = perform_leave_one_out(harm_df) if len(harm_df) >= 3 else pd.DataFrame()
    if not loo_df.empty:
        loo_csv_path = os.path.join(gene_out_dir, "leave_one_out.csv")
        loo_df.to_csv(loo_csv_path, index=False)
        
    # Outlier Detection (MR-PRESSO approximation)
    clean_harm_df, outlier_snps, ivw_no_outliers = flag_outliers_mr_presso_style(harm_df, threshold_sd=3.0)
    if outlier_snps:
        print(f"[NOTE] [{gene_symbol}] Flagged {len(outlier_snps)} potential outlier SNP(s): {outlier_snps}")
        print(f"       IVW without outliers: β = {ivw_no_outliers['beta']:0.3f}, p = {ivw_no_outliers['pval']:0.3e}")

    # -------------------------------------------------------------------------
    # STEP 6: Visualization
    # -------------------------------------------------------------------------
    scatter_path = os.path.join(gene_out_dir, "scatter_plot.png")
    plot_scatter(harm_df, gene_symbol, ivw_res, egger_res, wm_res, scatter_path)
    
    if len(harm_df) >= 2:
        funnel_path = os.path.join(gene_out_dir, "funnel_plot.png")
        plot_funnel(harm_df, gene_symbol, ivw_res, funnel_path)
        
    if not loo_df.empty:
        loo_plot_path = os.path.join(gene_out_dir, "leave_one_out_plot.png")
        plot_leave_one_out(loo_df, gene_symbol, ivw_res, loo_plot_path)

    # -------------------------------------------------------------------------
    # STEP 7: Aggregate Summary Row
    # -------------------------------------------------------------------------
    summary_row = {
        "gene": gene_symbol,
        "n_snps": len(harm_df),
        "f_stat": mean_f,
        "is_weak_instrument": is_weak,
        "ivw_beta": ivw_res["beta"],
        "ivw_se": ivw_res["se"],
        "ivw_ci_lower": ivw_res["ci_lower"],
        "ivw_ci_upper": ivw_res["ci_upper"],
        "ivw_pval": ivw_res["pval"],
        "egger_beta": egger_res["beta"] if egger_res else np.nan,
        "egger_se": egger_res["se"] if egger_res else np.nan,
        "egger_pval": egger_res["pval"] if egger_res else np.nan,
        "egger_intercept": egger_res["intercept"] if egger_res else np.nan,
        "egger_intercept_pval": egger_res["intercept_pval"] if egger_res else np.nan,
        "wm_beta": wm_res["beta"] if wm_res else np.nan,
        "wm_se": wm_res["se"] if wm_res else np.nan,
        "wm_pval": wm_res["pval"] if wm_res else np.nan,
        "cochrans_q": q_res["q_stat"],
        "cochrans_q_pval": q_res["q_pval"],
        "n_outliers_flagged": len(outlier_snps),
        "ivw_beta_no_outliers": ivw_no_outliers["beta"] if outlier_snps else np.nan,
        "ivw_pval_no_outliers": ivw_no_outliers["pval"] if outlier_snps else np.nan,
    }
    
    print(f"[{gene_symbol}] SUMMARY: IVW β = {ivw_res['beta']:0.3f} (p = {ivw_res['pval']:0.3e}) | Q-pval = {q_res['q_pval']:0.3f}")
    return summary_row


def main():
    """
    Main entry point for running Two-Sample MR across all specified hub genes.
    """
    print("=" * 80)
    print(" TWO-SAMPLE MENDELIAN RANDOMIZATION PIPELINE (Python / OpenGWAS)")
    print("=" * 80)
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    token = os.environ.get("OPENGWAS_TOKEN", None)
    
    genes_to_process = {}
    
    # 1. Determine exposure inputs (API vs Local files)
    if HUB_GENE_EXPOSURE_IDS:
        print("[MODE] Using OpenGWAS API Study IDs.")
        genes_to_process = HUB_GENE_EXPOSURE_IDS
    else:
        print("[MODE] Checking for local exposure CSV files (exposure_[GENE].csv)...")
        if os.path.exists(LOCAL_DATA_DIR):
            for fname in os.listdir(LOCAL_DATA_DIR):
                if fname.startswith("exposure_") and fname.endswith(".csv"):
                    gene = fname.replace("exposure_", "").replace(".csv", "").upper()
                    genes_to_process[gene] = os.path.join(LOCAL_DATA_DIR, fname)
                    
    if not genes_to_process:
        print("\n[INFO] No hub genes configured in HUB_GENE_EXPOSURE_IDS and no local exposure_[GENE].csv found.")
        print("To run the pipeline:")
        print("1. Set your OpenGWAS token:")
        print("   set OPENGWAS_TOKEN=your_token_here (Windows) or export OPENGWAS_TOKEN=your_token_here (Linux/Mac)")
        print("2. Search for study IDs using: python -c 'from data_fetch import search_studies; print(search_studies(\"SOCS1\", batch=\"eqtl-a\"))'")
        print("3. Add the IDs to HUB_GENE_EXPOSURE_IDS in run_mr_pipeline.py and re-run.")
        sys.exit(0)
        
    # 2. Determine outcome input
    outcome_target = OUTCOME_ID
    if not OUTCOME_ID and os.path.exists(LOCAL_OUTCOME_FILE):
        outcome_target = LOCAL_OUTCOME_FILE
        
    print(f"Target Outcome: {outcome_target}")
    print(f"Hub Genes to analyze: {list(genes_to_process.keys())}\n")
    
    all_summaries = []
    
    for gene, exposure_ref in genes_to_process.items():
        try:
            summary = run_pipeline_for_gene(
                gene_symbol=gene,
                exposure_id_or_file=exposure_ref,
                outcome_id_or_file=outcome_target,
                out_dir=RESULTS_DIR,
                token=token
            )
            if summary:
                all_summaries.append(summary)
        except Exception as e:
            print(f"[ERROR] Exception occurred while processing {gene}: {e}")
            import traceback
            traceback.print_exc()

    # 3. Aggregate results and cross-gene Forest Plot
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        summary_csv_path = os.path.join(RESULTS_DIR, "mr_summary_all_genes.csv")
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"\n[SUCCESS] Master MR summary saved to: {summary_csv_path}")
        
        # Cross-gene Forest plot
        forest_plot_path = os.path.join(RESULTS_DIR, "forest_plot_all_genes.png")
        plot_cross_gene_forest(
            summary_df,
            forest_plot_path,
            title="Two-Sample MR Causal Estimates: Hub Genes -> Kidney/Peritoneal Outcome"
        )
        print(f"[SUCCESS] Cross-gene Forest Plot saved to: {forest_plot_path}")
    else:
        print("\n[WARNING] No genes successfully completed MR analysis.")


if __name__ == "__main__":
    main()
