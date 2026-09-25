"""
transcriptome_wide_mr.py
========================
High-Throughput Transcriptome-Wide Two-Sample Mendelian Randomization (MR) Screening Pipeline.

Real Confirmed Schemas:
-----------------------
1. eQTLGen:
   Columns (TAB-delimited):
   Pvalue, SNP, SNPChr, SNPPos, AssessedAllele, OtherAllele, Zscore, 
   Gene, GeneSymbol, GeneChr, GenePos, NrCohorts, NrSamples, FDR, BonferroniP

2. CKDGen eGFR (Wuttke et al. 2019, b37):
   Columns (SPACE-delimited):
   Chr, Pos_b37, RSID, Allele1, Allele2, Freq1, Effect, StdErr, P-value, n_total_sum

Key Methodological Details:
---------------------------
1. Sample Size N:
   NrSamples (NOT NrCohorts) from eQTLGen is used as sample size N.
2. Minor Allele Frequency (MAF):
   eQTLGen lacks an EAF/MAF column. The pipeline matches each SNP to CKDGen's Freq1
   as an empirical proxy MAF. For any SNP without a matched frequency, it falls back
   to an assumed MAF of 0.5.
   * Methodological limitation: Proxy MAF source is used for Z-score to beta/se scaling;
     Z-scores, p-values, and causal directions remain exact.
3. Freq1 Assumption:
   Allele1 is treated as the effect allele with Freq1 as its frequency (standard GWAS convention).
"""

import os
import sys
import gzip
import time
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
from typing import Dict, List, Optional, Tuple, Union


# =============================================================================
# CONFIGURATION
# =============================================================================

# Default file paths
EQTLGEN_FILE = os.path.join(
    "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
    "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
)
OUTCOME_GWAS_FILE = os.path.join(
    "20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt",
    "20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt"
)
DEG_ECM_CANDIDATES_FILE = "convergent_ECM_DEGs_nominal.csv"

# Output Directory
OUTPUT_DIR = "transcriptome_wide_mr_results"

# Statistical Thresholds
PVAL_EXPOSURE_THRESHOLD = 5e-8    # Primary genome-wide significance
PVAL_RELAXED_FALLBACK = 1e-5      # Relaxed fallback if gene has no genome-wide hits
PALINDROME_THRESHOLD = 0.08       # EAF tolerance around 0.5 for palindromes (0.42 - 0.58)
CLUMP_KB_WINDOW = 500             # Window in kb to select independent cis-eQTLs per gene
DEFAULT_ASSUMED_MAF = 0.5         # Fallback assumed MAF if no frequency is available


def resolve_file_path(path: str) -> str:
    """Resolves potential variations of file paths (directory nesting, .gz extension)."""
    if os.path.exists(path):
        if os.path.isdir(path):
            for fname in os.listdir(path):
                if fname.endswith(('.txt', '.gz', '.tsv', '.csv')):
                    return os.path.join(path, fname)
        return path
    
    if path.endswith('.gz') and os.path.exists(path[:-3]):
        return resolve_file_path(path[:-3])
    
    base = os.path.basename(path)
    if base.endswith('.gz'):
        base = base[:-3]
    if os.path.exists(base):
        return resolve_file_path(base)

    return path


# =============================================================================
# 1. Z-SCORE TO BETA / SE CONVERSION
# =============================================================================

def zscore_to_beta_se(
    z: np.ndarray,
    n: np.ndarray,
    maf: Optional[np.ndarray] = None,
    default_maf: float = 0.5,
    warn_fallback: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Converts eQTLGen Z-scores to effect size (beta) and standard error (se).

    Formula (Zhu et al. 2016 / eQTLGen standard):
    ---------------------------------------------
    se = 1.0 / sqrt(2 * MAF * (1 - MAF) * (N + Z^2))
    beta = Z * se
    """
    if maf is None or np.all(np.isnan(maf)):
        if warn_fallback:
            print("\n" + "!" * 80)
            print(" [METHODOLOGICAL DISCLOSURE] eQTLGen lacks an internal MAF/EAF column.")
            print(f" MAF proxy source: outcome GWAS frequency where available, else assumed {default_maf}.")
            print(" * NOTE: Z-scores, p-values, and causal directions remain exact.")
            print(" * Absolute beta/SE magnitudes are approximate under this assumption.")
            print("!" * 80 + "\n")
        maf_clean = np.full_like(z, default_maf, dtype=float)
    else:
        n_missing = np.sum(np.isnan(maf) | (maf <= 0))
        if n_missing > 0 and warn_fallback:
            print(f"[NOTE] {n_missing:,} SNPs lacked matched outcome frequency; defaulted to assumed MAF = {default_maf}.")
        maf_clean = np.where(np.isnan(maf) | (maf <= 0), default_maf, maf)
        maf_clean = np.where(maf_clean > 0.5, 1.0 - maf_clean, maf_clean)
        maf_clean = np.clip(maf_clean, 0.001, 0.5)

    denom = 2.0 * maf_clean * (1.0 - maf_clean) * (n + z**2)
    se = 1.0 / np.sqrt(denom)
    beta = z * se
    return beta, se


# =============================================================================
# 2. BULK FILE LOADERS & NORMALIZERS
# =============================================================================

def load_and_standardize_gwas(
    file_path: str,
    is_outcome: bool = True,
    nrows: Optional[int] = None
) -> pd.DataFrame:
    """
    Loads GWAS summary statistics or eQTLGen cis-eQTL file and standardizes headers.
    """
    resolved_path = resolve_file_path(file_path)
    print(f"Loading {'Outcome GWAS' if is_outcome else 'eQTLGen Exposure'} from: {resolved_path} ...")
    compression = "gzip" if resolved_path.endswith(".gz") else None
    
    with open(resolved_path, 'r', encoding='utf-8') as f:
        first_line = f.readline()
    
    sep = "\t" if "\t" in first_line else r"\s+"
    
    df = pd.read_csv(
        resolved_path,
        sep=sep,
        engine='c' if sep == "\t" else 'python',
        compression=compression,
        nrows=nrows
    )
    print(f"Loaded {len(df):,} rows with original headers: {df.columns.tolist()}")

    col_map = {
        # SNP rsID
        "snp": "SNP", "rsid": "SNP", "snpid": "SNP", "variant_id": "SNP", "markername": "SNP", "id": "SNP",
        
        # Position
        "chr": "chr", "chromosome": "chr", "snpchr": "chr", "#chrom": "chr", "genechr": "gene_chr",
        "pos": "pos", "position": "pos", "snppos": "pos", "bp": "pos", "pos_b37": "pos", "genepos": "gene_pos",
        
        # Alleles
        "effect_allele": "effect_allele", "ea": "effect_allele", "a1": "effect_allele", "assessedallele": "effect_allele", "allele1": "effect_allele",
        "other_allele": "other_allele", "nea": "other_allele", "oa": "other_allele", "a2": "other_allele", "otherallele": "other_allele", "allele2": "other_allele",
        
        # Statistics
        "beta": "beta", "b": "beta", "effect": "beta", "zscore": "zscore", "z": "zscore",
        "se": "se", "stderr": "se", "std_err": "se",
        "pval": "pval", "p": "pval", "pvalue": "pval", "p_value": "pval", "p-value": "pval", "bonferronip": "bonferroni_p",
        
        # Frequencies
        "eaf": "eaf", "frq": "eaf", "freq": "eaf", "freq1": "eaf", "maf": "maf", "eaf_ukb": "eaf", "freq_tested_allele": "eaf",
        
        # Sample size
        "samplesize": "samplesize", "n": "samplesize", "nrsamples": "samplesize", "total_n": "samplesize", 
        "totalsamplesize": "samplesize", "n_total_sum": "samplesize", "nrcohorts": "nr_cohorts",
        
        # Gene annotations
        "genesymbol": "GeneSymbol", "gene": "Gene", "hgnc_symbol": "GeneSymbol", "fdr": "eQTL_FDR"
    }

    mapped_cols = {}
    for col in df.columns:
        clean_col = col.lower().strip().replace(" ", "").replace("_", "")
        for pattern, target in col_map.items():
            if clean_col == pattern.replace("_", ""):
                mapped_cols[col] = target
                break

    df = df.rename(columns=mapped_cols)
    return df


# =============================================================================
# 3. HIGH-THROUGHPUT HARMONIZATION & VECTORIZED MR SCREENING ENGINE
# =============================================================================

def run_transcriptome_wide_mr(
    eqtl_file: str,
    outcome_file: str,
    output_dir: str = "transcriptome_wide_mr_results",
    nrows: Optional[int] = None
) -> pd.DataFrame:
    """
    Executes transcriptome-wide Mendelian Randomization across all genes in eQTLGen.
    """
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    eqtl_resolved = resolve_file_path(eqtl_file)
    outcome_resolved = resolve_file_path(outcome_file)

    # 1. Load Outcome GWAS (CKDGen)
    out_df = load_and_standardize_gwas(outcome_resolved, is_outcome=True, nrows=nrows)
    out_df["SNP"] = out_df["SNP"].astype(str).str.strip()
    out_df = out_df.dropna(subset=["SNP", "beta", "se", "pval"]).drop_duplicates(subset=["SNP"])
    print(f"Outcome GWAS indexed: {len(out_df):,} unique SNPs.")

    # 2. Load eQTLGen bulk data
    eqtl_df = load_and_standardize_gwas(eqtl_resolved, is_outcome=False, nrows=nrows)
    eqtl_df["SNP"] = eqtl_df["SNP"].astype(str).str.strip()
    
    if "GeneSymbol" not in eqtl_df.columns and "Gene" in eqtl_df.columns:
        eqtl_df["GeneSymbol"] = eqtl_df["Gene"]
        
    eqtl_df = eqtl_df.dropna(subset=["SNP", "GeneSymbol", "pval"])
    total_eqtl_genes = eqtl_df['GeneSymbol'].nunique()
    print(f"eQTLGen file indexed: {len(eqtl_df):,} associations across {total_eqtl_genes:,} unique genes.")

    # 3. Filter eQTLs by P-value threshold
    eqtl_sig = eqtl_df[eqtl_df["pval"] < PVAL_RELAXED_FALLBACK].copy()
    print(f"eQTLs passing p < {PVAL_RELAXED_FALLBACK}: {len(eqtl_sig):,} rows across {eqtl_sig['GeneSymbol'].nunique():,} genes.")

    # Suffix isolation before merging
    eqtl_renamed = eqtl_sig.rename(columns={c: f"{c}_exp" for c in eqtl_sig.columns if c not in ["SNP", "GeneSymbol"]})
    out_renamed = out_df.rename(columns={c: f"{c}_out" for c in out_df.columns if c != "SNP"})

    # 4. Fast Inner Join on SNP rsID
    print("Merging eQTLGen Exposure with Outcome GWAS on SNP ID...")
    merged = pd.merge(
        eqtl_renamed,
        out_renamed,
        on="SNP"
    )
    print(f"Matched {len(merged):,} SNP-gene pairs in Outcome GWAS.")
    if merged.empty:
        print("[ERROR] No overlapping SNPs found between eQTLGen and Outcome GWAS.")
        return pd.DataFrame()

    # 5. Convert Z-score to Beta/SE using NrSamples as N and CKDGen Freq1 as proxy MAF
    if "beta_exp" not in merged.columns or merged["beta_exp"].isna().all():
        print("Converting eQTLGen Z-scores to Beta and SE (using NrSamples as N, outcome Freq1 as proxy MAF)...")
        z_vals = merged["zscore_exp"].values if "zscore_exp" in merged.columns else (merged["pval_exp"].apply(lambda p: stats.norm.ppf(1 - p/2))).values
        n_vals = merged["samplesize_exp"].fillna(31684.0).values if "samplesize_exp" in merged.columns else np.full(len(merged), 31684.0)

        # Matched outcome EAF (Freq1) as proxy MAF
        maf_array = merged["eaf_out"].values if "eaf_out" in merged.columns else None

        b_exp, s_exp = zscore_to_beta_se(
            z=z_vals,
            n=n_vals,
            maf=maf_array,
            default_maf=DEFAULT_ASSUMED_MAF,
            warn_fallback=True
        )
        merged["beta_exp"] = b_exp
        merged["se_exp"] = s_exp

    # 6. Harmonize Effect Alleles
    ea_e = merged["effect_allele_exp"].astype(str).str.upper().values
    oa_e = merged["other_allele_exp"].astype(str).str.upper().values
    ea_o = merged["effect_allele_out"].astype(str).str.upper().values
    oa_o = merged["other_allele_out"].astype(str).str.upper().values

    direct_match = (ea_e == ea_o) & (oa_e == oa_o)
    swapped_match = (ea_e == oa_o) & (oa_e == ea_o)

    # Drop incompatible alleles
    valid_mask = direct_match | swapped_match
    merged = merged[valid_mask].copy()
    
    # Flip outcome beta for swapped alleles
    swapped_idx = merged[(merged["effect_allele_exp"].str.upper() == merged["other_allele_out"].str.upper())].index
    merged.loc[swapped_idx, "beta_out"] = -merged.loc[swapped_idx, "beta_out"]

    genes_in_harmonized = merged['GeneSymbol'].nunique()
    print(f"Surviving harmonized records: {len(merged):,} across {genes_in_harmonized:,} genes.")

    # 7. Group by Gene and Execute High-Throughput MR Screening with Progress Logging
    print("\nRunning Transcriptome-Wide MR Estimation across all genes...")
    gene_results = []
    
    grouped = merged.groupby("GeneSymbol")
    total_genes_to_process = len(grouped)
    print(f"Starting MR calculation for {total_genes_to_process:,} genes (logging every 1,000 genes)...")

    for i, (gene, group) in enumerate(grouped, 1):
        g_sorted = group.sort_values("pval_exp")
        
        # Primary vs Fallback threshold
        sig_snps = g_sorted[g_sorted["pval_exp"] < PVAL_EXPOSURE_THRESHOLD]
        if len(sig_snps) == 0:
            sig_snps = g_sorted[g_sorted["pval_exp"] < PVAL_RELAXED_FALLBACK]
            threshold_used = f"{PVAL_RELAXED_FALLBACK:0.1e} (relaxed)"
        else:
            threshold_used = f"{PVAL_EXPOSURE_THRESHOLD:0.1e}"
            
        if len(sig_snps) == 0:
            continue
            
        # Independent cis-instrument selection (greedy 500kb window pruning)
        pos_col = "pos_exp" if "pos_exp" in sig_snps.columns else ("pos_out" if "pos_out" in sig_snps.columns else None)
        if pos_col:
            sig_snps = sig_snps.sort_values("pval_exp")
            pruned_rows = []
            used_positions = []
            for _, r in sig_snps.iterrows():
                p = r[pos_col]
                if not any(abs(p - up) < (CLUMP_KB_WINDOW * 1000) for up in used_positions):
                    pruned_rows.append(r)
                    used_positions.append(p)
            instruments = pd.DataFrame(pruned_rows)
        else:
            instruments = sig_snps.drop_duplicates(subset=["SNP"]).head(10)

        n_snps = len(instruments)
        b_exp = instruments["beta_exp"].values
        s_exp = instruments["se_exp"].values
        b_out = instruments["beta_out"].values
        s_out = instruments["se_out"].values

        # Instrument Strength F-statistic
        f_stats = (b_exp / s_exp) ** 2
        mean_f = float(np.mean(f_stats))

        # MR Estimation
        if n_snps == 1:
            # Single-SNP: Wald Ratio
            b_mr = float(b_out[0] / b_exp[0])
            se_mr = float(s_out[0] / abs(b_exp[0]))
            z_mr = b_mr / se_mr if se_mr > 0 else 0.0
            p_mr = float(2.0 * (1.0 - stats.norm.cdf(abs(z_mr))))
            method = "Wald Ratio"
            q_stat, q_pval = 0.0, 1.0
        else:
            # Multi-SNP: Inverse-Variance Weighted (IVW)
            weights = (b_exp ** 2) / (s_out ** 2)
            sum_w = np.sum(weights)
            b_mr = float(np.sum((b_exp * b_out) / (s_out ** 2)) / sum_w)
            
            # Cochran's Q
            wald_r = b_out / b_exp
            q_stat = float(np.sum(weights * (wald_r - b_mr) ** 2))
            q_pval = float(1.0 - stats.chi2.cdf(q_stat, df=n_snps - 1))
            
            # Multiplicative Random Effects SE
            se_fixed = np.sqrt(1.0 / sum_w)
            if (q_stat / (n_snps - 1)) > 1.0:
                se_mr = float(se_fixed * np.sqrt(q_stat / (n_snps - 1)))
                t_crit = stats.t.ppf(0.975, df=n_snps - 1)
                p_mr = float(2.0 * (1.0 - stats.t.cdf(abs(b_mr / se_mr), df=n_snps - 1)))
            else:
                se_mr = float(se_fixed)
                t_crit = 1.96
                p_mr = float(2.0 * (1.0 - stats.norm.cdf(abs(b_mr / se_mr))))
            method = "IVW"

        ci_lower = b_mr - 1.96 * se_mr
        ci_upper = b_mr + 1.96 * se_mr

        gene_results.append({
            "GeneSymbol": gene,
            "N_SNPs": n_snps,
            "Threshold_Used": threshold_used,
            "Mean_F_stat": mean_f,
            "Is_Weak_Instrument": mean_f < 10.0,
            "MR_Method": method,
            "MR_Beta": b_mr,
            "MR_SE": se_mr,
            "MR_CI_Lower": ci_lower,
            "MR_CI_Upper": ci_upper,
            "MR_Pvalue": p_mr,
            "Cochran_Q": q_stat,
            "Cochran_Q_Pval": q_pval,
            "Top_Instrument_SNP": instruments.iloc[0]["SNP"]
        })

        if i % 1000 == 0 or i == total_genes_to_process:
            elapsed = time.time() - start_time
            print(f"  [Progress] Processed {i:,} / {total_genes_to_process:,} genes ({i/total_genes_to_process*100:.1f}%) | Elapsed: {elapsed:.1f}s")

    results_df = pd.DataFrame(gene_results)
    if results_df.empty:
        print("[WARNING] No genes successfully estimated.")
        return pd.DataFrame()

    # 8. Genome-Wide Multiple Testing Correction (Benjamini-Hochberg FDR)
    valid_p = results_df["MR_Pvalue"].dropna()
    p_sorted_idx = np.argsort(valid_p.values)
    n_tests = len(valid_p)
    fdr_values = np.zeros(n_tests)
    
    running_min = 1.0
    for rank, idx in enumerate(reversed(p_sorted_idx), 1):
        orig_rank = n_tests - rank + 1
        fdr_val = valid_p.values[idx] * n_tests / orig_rank
        running_min = min(running_min, fdr_val)
        fdr_values[idx] = min(running_min, 1.0)
        
    results_df["MR_FDR"] = fdr_values
    results_df = results_df.sort_values("MR_Pvalue").reset_index(drop=True)

    # Method breakdown
    n_wald = (results_df["MR_Method"] == "Wald Ratio").sum()
    n_ivw = (results_df["MR_Method"] == "IVW").sum()
    total_time = time.time() - start_time

    # Save complete transcriptome-wide MR screen
    all_out_path = os.path.join(output_dir, "mr_transcriptome_wide_all_results.csv")
    results_df.to_csv(all_out_path, index=False)
    results_df.to_csv("mr_transcriptome_wide_all_results.csv", index=False)
    
    # Export MR-significant subset
    mr_sig_df = results_df[results_df["MR_Pvalue"] < 0.05].copy()
    mr_sig_path = os.path.join(output_dir, "mr_significant_genes.csv")
    mr_sig_df.to_csv(mr_sig_path, index=False)
    mr_sig_df.to_csv("mr_significant_genes.csv", index=False)

    print("\n" + "=" * 80)
    print(" [FINAL SUMMARY] TRANSCRIPTOME-WIDE MR SCREEN COMPLETE")
    print("=" * 80)
    print(f"  * Total Genes in eQTLGen:                    {total_eqtl_genes:,}")
    print(f"  * Genes with Usable Harmonized Instruments:   {len(results_df):,} ({len(results_df)/total_eqtl_genes*100:.1f}%)")
    print(f"    - Single-SNP (Wald ratio):                  {n_wald:,} ({n_wald/len(results_df)*100:.1f}%)")
    print(f"    - Multi-SNP (IVW):                          {n_ivw:,} ({n_ivw/len(results_df)*100:.1f}%)")
    print(f"  * Genes with Nominal P < 0.05:               {(results_df['MR_Pvalue'] < 0.05).sum():,}")
    print(f"  * Genes with FDR < 0.05:                      {(results_df['MR_FDR'] < 0.05).sum():,}")
    print(f"  * Total Runtime:                              {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    print("=" * 80 + "\n")

    return results_df


# =============================================================================
# 4. CANDIDATE LIST INTERSECTION & DUAL-TIER CONVERGENCE
# =============================================================================

def intersect_with_candidates(
    mr_results_df: pd.DataFrame,
    candidate_file: str,
    output_dir: str = "transcriptome_wide_mr_results"
) -> pd.DataFrame:
    """
    Intersects the unbiased Transcriptome-wide MR findings with the DEG/ECM candidate list.
    Reports both FDR < 0.05 (Defensible Tier) and Nominal P < 0.05 (Exploratory Tier) side-by-side.
    """
    cand_resolved = resolve_file_path(candidate_file)
    if not os.path.exists(cand_resolved):
        print(f"[NOTE] Candidate file not found at {candidate_file}. Skipping candidate intersection.")
        return pd.DataFrame()

    print(f"\n--- INTERSECTING TRANSCRIPTOME-WIDE MR WITH CANDIDATE LIST ({cand_resolved}) ---")
    cand_df = pd.read_csv(cand_resolved)
    
    cand_col = "Gene Symbol" if "Gene Symbol" in cand_df.columns else cand_df.columns[0]
    candidate_genes = set(cand_df[cand_col].dropna().astype(str).str.strip().str.upper())
    print(f"Loaded {len(candidate_genes)} candidate genes from {cand_resolved}.")

    # MR Significant Genes (Nominal and FDR)
    mr_nominal_genes = set(mr_results_df[mr_results_df["MR_Pvalue"] < 0.05]["GeneSymbol"].str.upper())
    mr_fdr_genes = set(mr_results_df[mr_results_df["MR_FDR"] < 0.05]["GeneSymbol"].str.upper())

    overlap_nominal = candidate_genes.intersection(mr_nominal_genes)
    overlap_fdr = candidate_genes.intersection(mr_fdr_genes)

    print(f"\n==================== CONVERGENCE SUMMARY ====================")
    print(f"  * Primary Tier (Candidate List INTERSECT MR FDR < 0.05):     {len(overlap_fdr)} genes -> {sorted(list(overlap_fdr))}")
    print(f"  * Exploratory Tier (Candidate List INTERSECT MR Nom P < 0.05): {len(overlap_nominal)} genes -> {sorted(list(overlap_nominal))}")
    print(f"=============================================================\n")

    # Create merged annotated table
    intersect_df = mr_results_df[mr_results_df["GeneSymbol"].str.upper().isin(candidate_genes)].copy()
    
    merged_cand = pd.merge(
        intersect_df,
        cand_df.rename(columns={cand_col: "GeneSymbol"}),
        on="GeneSymbol",
        how="left"
    )
    
    intersect_path = os.path.join(output_dir, "mr_deg_ecm_intersection.csv")
    merged_cand.to_csv(intersect_path, index=False)
    merged_cand.to_csv("mr_deg_ecm_intersection.csv", index=False)
    print(f"[SUCCESS] Final Convergent Candidate Intersection saved to: {intersect_path}")

    # Plot Convergence Venn Diagram (Nominal and FDR overlay)
    fig, ax = plt.subplots(figsize=(9, 7), facecolor="#FAFAFA")
    v = venn2(
        subsets=(len(candidate_genes - mr_nominal_genes), len(mr_nominal_genes - candidate_genes), len(overlap_nominal)),
        set_labels=(f'DEG/ECM Candidate List\n({len(candidate_genes)} Genes)', f'Transcriptome-Wide MR Screen\n(Nominal P < 0.05)'),
        set_colors=('#6366F1', '#10B981'),
        alpha=0.65,
        ax=ax
    )
    venn2_circles(
        subsets=(len(candidate_genes - mr_nominal_genes), len(mr_nominal_genes - candidate_genes), len(overlap_nominal)),
        linestyle='solid', linewidth=2.0, color='#1E293B', ax=ax
    )
    
    for text in v.set_labels:
        if text: text.set_fontsize(12); text.set_fontweight('bold'); text.set_color('#1E293B')
    for text in v.subset_labels:
        if text: text.set_fontsize(14); text.set_fontweight('bold'); text.set_color('#0F172A')

    ax.set_title("Transcriptomic & Causal Convergence:\nDEG/ECM Candidates INTERSECT Transcriptome-Wide MR Screen", 
                 fontsize=14, fontweight='bold', pad=20, color='#0F172A')

    plt.tight_layout()
    venn_path = os.path.join(output_dir, "venn_mr_deg_convergence.png")
    plt.savefig(venn_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.savefig("venn_mr_deg_convergence.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"[SUCCESS] Venn Diagram saved to: {venn_path}")

    return merged_cand


# =============================================================================
# 5. CLI ENTRY POINT
# =============================================================================

def main():
    print("=" * 80)
    print(" TRANSCRIPTOME-WIDE TWO-SAMPLE MENDELIAN RANDOMIZATION SCREENING PIPELINE")
    print("=" * 80)
    
    if len(sys.argv) > 1 and sys.argv[1] in ["-h", "--help"]:
        print("Usage:")
        print("  python transcriptome_wide_mr.py [eqtl_file] [outcome_file] [candidates_file]")
        sys.exit(0)

    eqtl_path = sys.argv[1] if len(sys.argv) > 1 else EQTLGEN_FILE
    outcome_path = sys.argv[2] if len(sys.argv) > 2 else OUTCOME_GWAS_FILE
    cand_path = sys.argv[3] if len(sys.argv) > 3 else DEG_ECM_CANDIDATES_FILE

    eqtl_resolved = resolve_file_path(eqtl_path)
    outcome_resolved = resolve_file_path(outcome_path)
    cand_resolved = resolve_file_path(cand_path)

    if not os.path.exists(eqtl_resolved) or not os.path.exists(outcome_resolved):
        print("\n[SETUP NOTICE] Missing input bulk data files:")
        print(f"  - eQTLGen file:   '{eqtl_path}' -> '{eqtl_resolved}' (Exists: {os.path.exists(eqtl_resolved)})")
        print(f"  - Outcome GWAS:   '{outcome_path}' -> '{outcome_resolved}' (Exists: {os.path.exists(outcome_resolved)})")
        sys.exit(0)

    # Run full MR screen
    mr_results = run_transcriptome_wide_mr(eqtl_resolved, outcome_resolved, output_dir=OUTPUT_DIR)
    
    # Run candidate intersection
    if not mr_results.empty:
        intersect_with_candidates(mr_results, cand_resolved, output_dir=OUTPUT_DIR)


if __name__ == "__main__":
    main()
