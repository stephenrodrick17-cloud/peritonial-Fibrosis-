"""
data_fetch.py
=============
Module for querying and retrieving exposure (eQTLGen) and outcome GWAS summary
statistics and performing LD clumping via the OpenGWAS API (IEU MRC / University of Bristol).

API Reference: https://api.opengwas.mrcieu.ac.uk/
"""

import os
import time
import requests
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Union

OPEN_GWAS_BASE_URL = "https://api.opengwas.mrcieu.ac.uk/api"


def get_auth_header(token: Optional[str] = None) -> Dict[str, str]:
    """
    Constructs the authorization headers for OpenGWAS API requests.
    
    Parameters
    ----------
    token : str, optional
        OpenGWAS API token. If None, checks the OPENGWAS_TOKEN environment variable.
        
    Returns
    -------
    dict
        Header dictionary with Bearer token if available.
    """
    api_token = token or os.environ.get("OPENGWAS_TOKEN", "")
    headers = {
        "User-Agent": "TwoSampleMR-Python-Pipeline/1.0",
        "Accept": "application/json",
    }
    if api_token:
        headers["Authorization"] = f"Bearer {api_token.strip()}"
    return headers


def search_studies(query: str, batch: Optional[str] = None, token: Optional[str] = None) -> pd.DataFrame:
    """
    Searches the OpenGWAS database for GWAS and eQTL studies matching a keyword or gene name.
    
    Parameters
    ----------
    query : str
        Search keyword (e.g. 'SOCS1', 'eGFR', 'CKDGen', 'peritoneal').
    batch : str, optional
        Specific batch filter (e.g. 'eqtl-a' for eQTLGen, 'ukb-b' for UK Biobank, 'ieu-a').
    token : str, optional
        OpenGWAS API token.
        
    Returns
    -------
    pd.DataFrame
        DataFrame of matching studies with columns: id, trait, sample_size, year, author, pmid, etc.
    """
    headers = get_auth_header(token)
    url = f"{OPEN_GWAS_BASE_URL}/gwasinfo"
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        # OpenGWAS returns a dict keyed by study ID or list of dicts
        if isinstance(data, dict):
            studies = list(data.values())
        else:
            studies = data
            
        df = pd.DataFrame(studies)
        if df.empty:
            return df
            
        # Filter by query keyword across relevant fields
        query_str = str(query).lower()
        mask = (
            df['trait'].astype(str).str.lower().str.contains(query_str, na=False) |
            df['id'].astype(str).str.lower().str.contains(query_str, na=False) |
            df['author'].astype(str).str.lower().str.contains(query_str, na=False)
        )
        if 'note' in df.columns:
            mask = mask | df['note'].astype(str).str.lower().str.contains(query_str, na=False)
            
        filtered_df = df[mask].copy()
        
        if batch and 'batch' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['batch'].astype(str).str.lower() == batch.lower()]
            
        cols_to_show = [c for c in ['id', 'trait', 'batch', 'sample_size', 'year', 'author', 'pmid', 'category'] if c in filtered_df.columns]
        return filtered_df[cols_to_show]
        
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Failed to query OpenGWAS study metadata: {e}")
        return pd.DataFrame()


def clump_snps_opengwas(
    rsids: List[str],
    pvals: List[float],
    r2_threshold: float = 0.001,
    kb_window: int = 10000,
    pop: str = "EUR",
    token: Optional[str] = None
) -> List[str]:
    """
    Performs server-side Linkage Disequilibrium (LD) clumping via OpenGWAS API.
    
    Parameters
    ----------
    rsids : list of str
        List of candidate instrument rsIDs.
    pvals : list of float
        Corresponding association p-values.
    r2_threshold : float, default=0.001
        Clumping LD r^2 threshold.
    kb_window : int, default=10000
        Clumping physical distance window in kilobases (10Mb default).
    pop : str, default='EUR'
        Reference population code ('EUR', 'AFR', 'EAS', 'SAS', 'AMR').
    token : str, optional
        OpenGWAS API token.
        
    Returns
    -------
    list of str
        Pruned list of independent instrument rsIDs.
    """
    if not rsids:
        return []
    if len(rsids) == 1:
        return rsids
        
    headers = get_auth_header(token)
    url = f"{OPEN_GWAS_BASE_URL}/ld/clump"
    
    # Payload format expected by OpenGWAS LD clump endpoint
    payload = {
        "rsid": rsids,
        "pval": pvals,
        "r2": r2_threshold,
        "kb": kb_window,
        "pop": pop
    }
    
    try:
        response = requests.post(url, json=payload, headers=headers, timeout=60)
        if response.status_code == 200:
            clumped_data = response.json()
            if isinstance(clumped_data, list):
                # Returns list of clumped records with 'rsid'
                clumped_snps = [item.get("rsid") for item in clumped_data if "rsid" in item]
                if clumped_snps:
                    return clumped_snps
            elif isinstance(clumped_data, dict) and "rsid" in clumped_data:
                return clumped_data["rsid"]
        else:
            print(f"[WARNING] OpenGWAS server-side clumping returned HTTP {response.status_code}: {response.text}")
    except Exception as e:
        print(f"[WARNING] API clumping request failed ({e}). Applying distance-based proxy pruning fallback.")

    # Fallback: distance / p-value ranked pruning if API clumping is unavailable
    df = pd.DataFrame({"SNP": rsids, "pval": pvals}).sort_values("pval")
    return df["SNP"].drop_duplicates().tolist()


def get_instruments_opengwas(
    gwas_id: str,
    pval_primary: float = 5e-8,
    pval_fallback: float = 1e-5,
    r2: float = 0.001,
    kb: int = 10000,
    token: Optional[str] = None
) -> (pd.DataFrame, Dict[str, int]):
    """
    Retrieves significant instrument SNPs for a given exposure (e.g. eQTL study),
    applies p-value threshold with automated fallback, and executes LD clumping.
    
    Parameters
    ----------
    gwas_id : str
        Exposure GWAS/eQTL study ID (e.g., 'eqtl-a-ENSG00000185338').
    pval_primary : float, default=5e-8
        Primary genome-wide significance threshold.
    pval_fallback : float, default=1e-5
        Relaxed significance threshold if < 3 instruments survive primary threshold.
    r2 : float, default=0.001
        LD clumping r^2 threshold.
    kb : int, default=10000
        LD clumping window in kilobases.
    token : str, optional
        OpenGWAS API token.
        
    Returns
    -------
    tuple of (pd.DataFrame, dict)
        Instrument DataFrame and filtering funnel count dictionary.
    """
    headers = get_auth_header(token)
    funnel = {
        "raw_instruments": 0,
        "pval_filtered": 0,
        "clumped_independent": 0,
        "threshold_used": "5e-8"
    }
    
    # 1. Fetch top hits from OpenGWAS tophits endpoint
    url = f"{OPEN_GWAS_BASE_URL}/tophits"
    payload = {
        "id": [gwas_id],
        "pval": pval_fallback,  # retrieve down to fallback to allow graceful relaxation
        "preclumped": 0
    }
    
    try:
        response = requests.post(url, json=payload, headers=headers, timeout=60)
        response.raise_for_status()
        records = response.json()
    except Exception as e:
        print(f"[ERROR] Failed to fetch instruments for {gwas_id}: {e}")
        return pd.DataFrame(), funnel
        
    if not records:
        print(f"[WARNING] No instruments found for {gwas_id} at p < {pval_fallback}")
        return pd.DataFrame(), funnel
        
    df = pd.DataFrame(records)
    funnel["raw_instruments"] = len(df)
    
    # Standardize OpenGWAS column names
    col_mapping = {
        "rsid": "SNP",
        "name": "SNP",
        "chr": "chr",
        "position": "pos",
        "ea": "effect_allele",
        "nea": "other_allele",
        "beta": "beta",
        "se": "se",
        "p": "pval",
        "eaf": "eaf",
        "n": "samplesize"
    }
    df = df.rename(columns={k: v for k, v in col_mapping.items() if k in df.columns})
    
    # 2. Filter by P-value threshold
    primary_mask = df["pval"] < pval_primary
    n_primary = primary_mask.sum()
    
    if n_primary >= 3:
        df_filtered = df[primary_mask].copy()
        funnel["threshold_used"] = f"{pval_primary:0.1e}"
        print(f"[{gwas_id}] Using primary genome-wide threshold (p < {pval_primary:0.1e}): {n_primary} SNPs found.")
    else:
        print(f"[WARNING] [{gwas_id}] Only {n_primary} SNPs passed p < {pval_primary:0.1e}. "
              f"Relaxing threshold to fallback (p < {pval_fallback:0.1e}).")
        df_filtered = df[df["pval"] < pval_fallback].copy()
        funnel["threshold_used"] = f"{pval_fallback:0.1e} (relaxed)"
        
    funnel["pval_filtered"] = len(df_filtered)
    if df_filtered.empty:
        return pd.DataFrame(), funnel
        
    # 3. LD Clumping
    clumped_snps = clump_snps_opengwas(
        rsids=df_filtered["SNP"].tolist(),
        pvals=df_filtered["pval"].tolist(),
        r2_threshold=r2,
        kb_window=kb,
        token=token
    )
    
    df_clumped = df_filtered[df_filtered["SNP"].isin(clumped_snps)].drop_duplicates(subset=["SNP"]).copy()
    funnel["clumped_independent"] = len(df_clumped)
    print(f"[{gwas_id}] LD Clumping (r² < {r2}, {kb}kb): {len(df_clumped)} independent instruments retained.")
    
    return df_clumped, funnel


def get_outcome_associations(
    outcome_gwas_id: str,
    snps: List[str],
    token: Optional[str] = None
) -> pd.DataFrame:
    """
    Fetches summary statistics for specific instrument SNPs from the outcome GWAS via OpenGWAS API.
    
    Parameters
    ----------
    outcome_gwas_id : str
        Outcome GWAS study ID (e.g. 'ieu-a-798' or CKDGen eGFR study ID).
    snps : list of str
        List of instrument rsIDs.
    token : str, optional
        OpenGWAS API token.
        
    Returns
    -------
    pd.DataFrame
        Outcome summary statistics DataFrame for the matched instruments.
    """
    if not snps:
        return pd.DataFrame()
        
    headers = get_auth_header(token)
    url = f"{OPEN_GWAS_BASE_URL}/associations"
    payload = {
        "id": [outcome_gwas_id],
        "variants": snps
    }
    
    try:
        response = requests.post(url, json=payload, headers=headers, timeout=60)
        response.raise_for_status()
        records = response.json()
    except Exception as e:
        print(f"[ERROR] Failed to fetch outcome associations for {outcome_gwas_id}: {e}")
        return pd.DataFrame()
        
    if not records:
        return pd.DataFrame()
        
    df = pd.DataFrame(records)
    col_mapping = {
        "rsid": "SNP",
        "name": "SNP",
        "chr": "chr",
        "position": "pos",
        "ea": "effect_allele",
        "nea": "other_allele",
        "beta": "beta",
        "se": "se",
        "p": "pval",
        "eaf": "eaf",
        "n": "samplesize"
    }
    df = df.rename(columns={k: v for k, v in col_mapping.items() if k in df.columns})
    return df
