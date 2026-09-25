"""
mr_stats.py
===========
Core statistical engine for Two-Sample Mendelian Randomization (MR).
Implements Harmonization, Instrument Strength (F-stat), IVW, MR-Egger,
Weighted Median, Cochran's Q, Leave-One-Out, and MR-PRESSO-style outlier checks.

All methods are documented with explicit mathematical formulas and assumptions.
"""

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from typing import Dict, Tuple, Optional, List


# -----------------------------------------------------------------------------
# 1. HARMONIZATION
# -----------------------------------------------------------------------------

def is_palindromic(a1: str, a2: str) -> bool:
    """
    Checks if a pair of alleles is palindromic (A/T or C/G) on DNA strands.
    """
    pair = {str(a1).upper(), str(a2).upper()}
    return pair in [{'A', 'T'}, {'C', 'G'}]


def harmonize_data(
    exposure_df: pd.DataFrame,
    outcome_df: pd.DataFrame,
    palindrome_threshold: float = 0.08,
    action: int = 2
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Harmonizes exposure and outcome datasets by matching SNPs, aligning effect alleles,
    and handling palindromic variants.

    Statistical Logic:
    ------------------
    - If Exposure Alleles = (A1, A2) and Outcome Alleles = (A1, A2):
        Direct match. Outcome Beta remains unchanged.
    - If Exposure Alleles = (A1, A2) and Outcome Alleles = (A2, A1):
        Alleles swapped. Outcome Beta is multiplied by -1.
    - Palindromic SNPs (A/T or C/G):
        If EAF is available and close to 0.5 (within 0.5 +/- palindrome_threshold),
        strand orientation cannot be unambiguously inferred. These SNPs are dropped.

    Parameters
    ----------
    exposure_df : pd.DataFrame
        Columns required: SNP, effect_allele, other_allele, beta, se, pval, eaf (optional).
    outcome_df : pd.DataFrame
        Columns required: SNP, effect_allele, other_allele, beta, se, pval, eaf (optional).
    palindrome_threshold : float, default=0.08
        Tolerance window around 0.5 for ambiguous palindromic allele frequencies (0.42 - 0.58).
    action : int, default=2
        Harmonization strictness level:
        1: Assume all alleles on forward strand.
        2: Try to infer forward strand from allele frequencies (standard TwoSampleMR).

    Returns
    -------
    tuple of (pd.DataFrame, dict)
        Harmonized DataFrame and funnel tracking dictionary.
    """
    funnel = {
        "exposure_snps": len(exposure_df),
        "outcome_available_snps": len(outcome_df),
        "matched_on_id": 0,
        "dropped_incompatible_alleles": 0,
        "dropped_ambiguous_palindromes": 0,
        "harmonized_final": 0
    }

    # Merge on SNP identifier
    merged = pd.merge(
        exposure_df,
        outcome_df,
        on="SNP",
        suffixes=("_exp", "_out"),
        how="inner"
    )
    funnel["matched_on_id"] = len(merged)

    if merged.empty:
        return pd.DataFrame(), funnel

    harmonized_rows = []

    for _, row in merged.iterrows():
        ea_exp = str(row["effect_allele_exp"]).upper().strip()
        oa_exp = str(row["other_allele_exp"]).upper().strip()
        ea_out = str(row["effect_allele_out"]).upper().strip()
        oa_out = str(row["other_allele_out"]).upper().strip()

        b_exp = float(row["beta_exp"])
        se_exp = float(row["se_exp"])
        b_out = float(row["beta_out"])
        se_out = float(row["se_out"])
        pval_exp = float(row["pval_exp"])
        pval_out = float(row["pval_out"])

        eaf_exp = float(row["eaf_exp"]) if ("eaf_exp" in row and pd.notna(row["eaf_exp"])) else np.nan
        eaf_out = float(row["eaf_out"]) if ("eaf_out" in row and pd.notna(row["eaf_out"])) else np.nan

        # Check palindromic status
        is_pal = is_palindromic(ea_exp, oa_exp)

        if is_pal:
            # Check if EAF allows resolving strand
            if not np.isnan(eaf_exp) and abs(eaf_exp - 0.5) < palindrome_threshold:
                funnel["dropped_ambiguous_palindromes"] += 1
                continue
            if not np.isnan(eaf_out) and abs(eaf_out - 0.5) < palindrome_threshold:
                funnel["dropped_ambiguous_palindromes"] += 1
                continue

        # Check allele alignment
        if ea_exp == ea_out and oa_exp == oa_out:
            # Direct match
            harmonized_rows.append({
                "SNP": row["SNP"],
                "effect_allele": ea_exp,
                "other_allele": oa_exp,
                "beta_exp": b_exp,
                "se_exp": se_exp,
                "pval_exp": pval_exp,
                "eaf_exp": eaf_exp,
                "beta_out": b_out,
                "se_out": se_out,
                "pval_out": pval_out,
                "eaf_out": eaf_out,
                "palindromic": is_pal,
                "flipped": False
            })
        elif ea_exp == oa_out and oa_exp == ea_out:
            # Swapped alleles: Flip outcome beta sign
            harmonized_rows.append({
                "SNP": row["SNP"],
                "effect_allele": ea_exp,
                "other_allele": oa_exp,
                "beta_exp": b_exp,
                "se_exp": se_exp,
                "pval_exp": pval_exp,
                "eaf_exp": eaf_exp,
                "beta_out": -b_out,
                "se_out": se_out,
                "pval_out": pval_out,
                "eaf_out": (1.0 - eaf_out) if not np.isnan(eaf_out) else np.nan,
                "palindromic": is_pal,
                "flipped": True
            })
        else:
            # Incompatible alleles (e.g. triallelic or strand mismatch without resolution)
            funnel["dropped_incompatible_alleles"] += 1

    harm_df = pd.DataFrame(harmonized_rows)
    funnel["harmonized_final"] = len(harm_df)
    return harm_df, funnel


# -----------------------------------------------------------------------------
# 2. INSTRUMENT STRENGTH (F-STATISTIC)
# -----------------------------------------------------------------------------

def calculate_f_statistic(df: pd.DataFrame) -> Tuple[pd.Series, float, bool]:
    """
    Calculates per-SNP F-statistic and the mean instrument F-statistic.

    Formula:
    --------
    F_j = (beta_exp,j / se_exp,j)^2
    Mean F = mean(F_j)

    Weak instrument threshold: Mean F < 10 (Stock & Yogo rule of thumb).

    Returns
    -------
    tuple of (pd.Series, float, bool)
        Per-SNP F-statistics, Mean F-statistic, and is_weak_flag (True if mean F < 10).
    """
    f_stats = (df["beta_exp"] / df["se_exp"]) ** 2
    mean_f = float(np.mean(f_stats)) if len(f_stats) > 0 else 0.0
    is_weak = mean_f < 10.0
    return f_stats, mean_f, is_weak


# -----------------------------------------------------------------------------
# 3. MR METHODS
# -----------------------------------------------------------------------------

def perform_ivw(df: pd.DataFrame) -> Dict[str, float]:
    """
    Inverse-Variance Weighted (IVW) linear regression.

    Statistical Logic:
    ------------------
    For each SNP j:
        Wald ratio estimate:  beta_j = beta_out,j / beta_exp,j
        Variance:             var_j = se_out,j^2 / beta_exp,j^2
        Weight:               w_j = 1 / var_j = beta_exp,j^2 / se_out,j^2

    Pooled IVW Estimate:
        beta_IVW = sum(w_j * beta_j) / sum(w_j)
                 = sum(beta_exp,j * beta_out,j / se_out,j^2) / sum(beta_exp,j^2 / se_out,j^2)

    Standard Error (Multiplicative Random Effects):
        SE_fixed = sqrt(1 / sum(w_j))
        Q = sum(w_j * (beta_j - beta_IVW)^2)
        SE_mre = SE_fixed * max(1.0, sqrt(Q / (K - 1)))  for K > 1 SNPs.

    Returns
    -------
    dict
        Dictionary containing beta, se, 95% CI (ci_lower, ci_upper), pval, method.
    """
    k = len(df)
    if k == 0:
        return {"beta": np.nan, "se": np.nan, "ci_lower": np.nan, "ci_upper": np.nan, "pval": np.nan, "n_snps": 0}

    beta_exp = df["beta_exp"].values
    beta_out = df["beta_out"].values
    se_out = df["se_out"].values

    if k == 1:
        # Single SNP: Wald Ratio
        b_wald = beta_out[0] / beta_exp[0]
        se_wald = se_out[0] / abs(beta_exp[0])
        z = b_wald / se_wald
        pval = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
        return {
            "beta": float(b_wald),
            "se": float(se_wald),
            "ci_lower": float(b_wald - 1.96 * se_wald),
            "ci_upper": float(b_wald + 1.96 * se_wald),
            "pval": float(pval),
            "n_snps": 1,
            "method": "Wald Ratio"
        }

    # Weights
    weights = (beta_exp ** 2) / (se_out ** 2)
    sum_w = np.sum(weights)

    # Weighted sum
    beta_ivw = np.sum((beta_exp * beta_out) / (se_out ** 2)) / sum_w

    # Cochran's Q
    wald_ratios = beta_out / beta_exp
    q_stat = np.sum(weights * (wald_ratios - beta_ivw) ** 2)

    # Multiplicative random-effects SE
    se_fixed = np.sqrt(1.0 / sum_w)
    if k > 1 and (q_stat / (k - 1)) > 1.0:
        se_ivw = se_fixed * np.sqrt(q_stat / (k - 1))
        # Use t-distribution for small N degrees of freedom
        t_crit = stats.t.ppf(0.975, df=k - 1)
        pval = 2.0 * (1.0 - stats.t.cdf(abs(beta_ivw / se_ivw), df=k - 1))
    else:
        se_ivw = se_fixed
        t_crit = 1.96
        pval = 2.0 * (1.0 - stats.norm.cdf(abs(beta_ivw / se_ivw)))

    return {
        "beta": float(beta_ivw),
        "se": float(se_ivw),
        "ci_lower": float(beta_ivw - t_crit * se_ivw),
        "ci_upper": float(beta_ivw + t_crit * se_ivw),
        "pval": float(pval),
        "n_snps": k,
        "method": "Inverse Variance Weighted"
    }


def perform_mr_egger(df: pd.DataFrame) -> Dict[str, float]:
    """
    MR-Egger Regression.

    Statistical Logic:
    ------------------
    Weighted linear regression of outcome effects on exposure effects WITH an intercept:
        (sign(beta_exp,j) * beta_out,j) = beta_0 + beta_Egger * |beta_exp,j| + epsilon_j
        Weights: w_j = 1 / se_out,j^2

    InSIDE Assumption: Instrument Strength Independent of Direct Effect.
    - Slope (beta_Egger): Causal effect estimate robust to directional pleiotropy.
    - Intercept (beta_0): Test for directional horizontal pleiotropy (H0: beta_0 = 0).

    Returns
    -------
    dict
        beta, se, ci_lower, ci_upper, pval, intercept, intercept_se, intercept_pval.
    """
    k = len(df)
    if k < 3:
        # MR-Egger requires at least 3 SNPs for estimation of slope + intercept
        return {
            "beta": np.nan, "se": np.nan, "ci_lower": np.nan, "ci_upper": np.nan, "pval": np.nan,
            "intercept": np.nan, "intercept_se": np.nan, "intercept_pval": np.nan, "n_snps": k
        }

    beta_exp = df["beta_exp"].values
    beta_out = df["beta_out"].values
    se_out = df["se_out"].values

    # Orient exposure effects to be positive (Burgess & Thompson, 2017)
    signs = np.sign(beta_exp)
    signs[signs == 0] = 1.0

    y = signs * beta_out
    x = np.abs(beta_exp)
    w = 1.0 / (se_out ** 2)

    X = sm.add_constant(x)
    model = sm.WLS(y, X, weights=w)
    results = model.fit()

    intercept = results.params[0]
    intercept_se = results.bse[0]
    intercept_pval = results.pvalues[0]

    slope = results.params[1]
    slope_se = results.bse[1]
    slope_pval = results.pvalues[1]

    t_crit = stats.t.ppf(0.975, df=k - 2)

    return {
        "beta": float(slope),
        "se": float(slope_se),
        "ci_lower": float(slope - t_crit * slope_se),
        "ci_upper": float(slope + t_crit * slope_se),
        "pval": float(slope_pval),
        "intercept": float(intercept),
        "intercept_se": float(intercept_se),
        "intercept_pval": float(intercept_pval),
        "n_snps": k,
        "method": "MR-Egger"
    }


def perform_weighted_median(df: pd.DataFrame, n_boot: int = 1000, seed: int = 42) -> Dict[str, float]:
    """
    Weighted Median Estimator (Bowden et al., 2016).

    Statistical Logic:
    ------------------
    Provides consistent causal estimates if at least 50% of the weight comes from valid instruments.
    1. Calculates Wald ratio estimates: beta_j = beta_out,j / beta_exp,j.
    2. Calculates inverse variance weights: w_j = beta_exp,j^2 / se_out,j^2.
    3. Normalizes weights: p_j = w_j / sum(w_j).
    4. Calculates the weighted median from sorted ratio estimates.
    5. Standard error is computed using parametric bootstrap resampling (n_boot iterations).

    Returns
    -------
    dict
        beta, se, ci_lower, ci_upper, pval.
    """
    k = len(df)
    if k < 3:
        return {"beta": np.nan, "se": np.nan, "ci_lower": np.nan, "ci_upper": np.nan, "pval": np.nan, "n_snps": k}

    beta_exp = df["beta_exp"].values
    beta_out = df["beta_out"].values
    se_exp = df["se_exp"].values
    se_out = df["se_out"].values

    def _calc_wm(b_e, b_o, s_o):
        ratios = b_o / b_e
        weights = (b_e ** 2) / (s_o ** 2)
        norm_w = weights / np.sum(weights)

        sort_idx = np.argsort(ratios)
        sorted_ratios = ratios[sort_idx]
        sorted_w = norm_w[sort_idx]

        cum_w = np.cumsum(sorted_w) - 0.5 * sorted_w
        # Linear interpolation for median
        return float(np.interp(0.5, cum_w, sorted_ratios))

    point_estimate = _calc_wm(beta_exp, beta_out, se_out)

    # Parametric bootstrap for standard error
    np.random.seed(seed)
    boot_estimates = []
    for _ in range(n_boot):
        b_e_boot = np.random.normal(beta_exp, se_exp)
        b_o_boot = np.random.normal(beta_out, se_out)
        wm_b = _calc_wm(b_e_boot, b_o_boot, se_out)
        boot_estimates.append(wm_b)

    se_boot = float(np.std(boot_estimates, ddof=1))
    z = point_estimate / se_boot if se_boot > 0 else 0.0
    pval = 2.0 * (1.0 - stats.norm.cdf(abs(z)))

    return {
        "beta": float(point_estimate),
        "se": float(se_boot),
        "ci_lower": float(point_estimate - 1.96 * se_boot),
        "ci_upper": float(point_estimate + 1.96 * se_boot),
        "pval": float(pval),
        "n_snps": k,
        "method": "Weighted Median"
    }


# -----------------------------------------------------------------------------
# 4. SENSITIVITY ANALYSES
# -----------------------------------------------------------------------------

def calculate_cochrans_q(df: pd.DataFrame, ivw_beta: Optional[float] = None) -> Dict[str, float]:
    """
    Cochran's Q Test for Heterogeneity across instrument SNP ratio estimates.

    Formula:
    --------
    Q = sum_j ( w_j * (beta_j - beta_IVW)^2 )
    where w_j = beta_exp,j^2 / se_out,j^2
    Degrees of freedom = K - 1.
    P-value from Chi-Square distribution: 1 - chi2.cdf(Q, df=K-1).

    Returns
    -------
    dict
        q_stat, q_df, q_pval.
    """
    k = len(df)
    if k < 2:
        return {"q_stat": 0.0, "q_df": 0, "q_pval": 1.0}

    beta_exp = df["beta_exp"].values
    beta_out = df["beta_out"].values
    se_out = df["se_out"].values

    weights = (beta_exp ** 2) / (se_out ** 2)
    wald_ratios = beta_out / beta_exp

    if ivw_beta is None:
        ivw_beta = np.sum((beta_exp * beta_out) / (se_out ** 2)) / np.sum(weights)

    q_stat = float(np.sum(weights * (wald_ratios - ivw_beta) ** 2))
    df_q = k - 1
    q_pval = float(1.0 - stats.chi2.cdf(q_stat, df=df_q))

    return {
        "q_stat": q_stat,
        "q_df": df_q,
        "q_pval": q_pval
    }


def perform_leave_one_out(df: pd.DataFrame) -> pd.DataFrame:
    """
    Leave-One-Out Sensitivity Analysis.
    Iteratively removes each SNP j and recomputes the IVW causal estimate with the remaining K - 1 SNPs.

    Returns
    -------
    pd.DataFrame
        Table containing SNP_removed, beta, se, ci_lower, ci_upper, pval.
    """
    k = len(df)
    results = []

    for i in range(k):
        sub_df = df.drop(df.index[i])
        ivw_res = perform_ivw(sub_df)
        results.append({
            "SNP_removed": df.iloc[i]["SNP"],
            "beta": ivw_res["beta"],
            "se": ivw_res["se"],
            "ci_lower": ivw_res["ci_lower"],
            "ci_upper": ivw_res["ci_upper"],
            "pval": ivw_res["pval"]
        })

    return pd.DataFrame(results)


def flag_outliers_mr_presso_style(df: pd.DataFrame, threshold_sd: float = 3.0) -> Tuple[pd.DataFrame, List[str], Dict[str, float]]:
    """
    Residual-based outlier detection (MR-PRESSO approximation).

    Statistical Logic:
    ------------------
    Standardized IVW regression residual for each SNP j:
        residual_j = (beta_out,j - beta_IVW * beta_exp,j) / se_out,j
    SNPs with |residual_j| > threshold_sd (default 3.0 standard deviations) are flagged as outliers.

    Returns
    -------
    tuple of (pd.DataFrame, list of str, dict)
        Cleaned DataFrame, list of outlier SNPs, and IVW results without outliers.
    """
    if len(df) < 3:
        return df, [], perform_ivw(df)

    ivw_full = perform_ivw(df)
    b_ivw = ivw_full["beta"]

    residuals = (df["beta_out"] - b_ivw * df["beta_exp"]) / df["se_out"]
    outlier_mask = np.abs(residuals) > threshold_sd
    outlier_snps = df.loc[outlier_mask, "SNP"].tolist()

    if outlier_snps:
        df_clean = df[~outlier_mask].copy()
        ivw_clean = perform_ivw(df_clean)
    else:
        df_clean = df.copy()
        ivw_clean = ivw_full

    return df_clean, outlier_snps, ivw_clean
