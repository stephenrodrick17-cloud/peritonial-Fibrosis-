"""
Script 65b: Stage 5 GSE248762 Single-Cell Stromal lncRNA Detection & Credibility Audit
Implements:
  1. Parsing of GENCODE v32 lncRNA GTF (matching 10x Cell Ranger reference).
  2. Single-cell count extraction from GSE248762_hubblind_allcells_qc.h5ad for the 3,670 stromal cells.
  3. Cell-level and donor-level detection metrics:
     - pct_cells_expressing = n_cells_gt0 / 3,670 * 100%
     - n_donors_expressing = donors with >= 1 expressing stromal cell (out of 16 donors)
  4. Multi-threshold sensitivity evaluation:
     - Percentages: 0.1%, 0.5%, 1.0%, 2.0%, 5.0%, 10.0%
     - Donors: >= 1 donor, >= 2 donors, >= 3 donors
  5. Primary threshold rule: >= 1.0% of stromal cells (>= 37 cells) AND >= 2 donors.
  6. Flagging of nuclear / housekeeping / promiscuous lncRNAs as 'low-credibility ceRNA candidates':
     - MALAT1, NEAT1, XIST, GAS5, NORAD, MEG3, H19, SNHG*
     - XIST explicitly marked: EXCLUDED_FROM_INTERPRETATION_DONOR_SEX
Outputs:
  - results/tables/G_lncRNA_threshold_sensitivity.csv
  - results/tables/G_lncRNA_detected_list.csv
  - results/tables/G_lncRNA_GSE248762_detection_summary.csv
  - results/tables/G_lncRNA_provenance.csv
"""

import os
import gzip
import re
import datetime
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse

ROOT = r"d:\Peritoneal Project"
DATA_DIR = os.path.join(ROOT, "data")
OUT_DIR = os.path.join(ROOT, "results", "tables")

def main():
    print("=" * 80)
    print("SCRIPT 65b: GSE248762 STROMAL lncRNA DETECTION THRESHOLDS & CREDIBILITY AUDIT")
    print("=" * 80)

    # 1. Parse GENCODE v32 GTF
    gtf_path = os.path.join(DATA_DIR, "raw", "gencode.v32.long_noncoding_RNAs.gtf.gz")
    print(f"Loading GENCODE v32 GTF: {gtf_path}")
    lnc_dict = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) > 8 and parts[2] == "gene":
                gid_m = re.search(r'gene_id "([^"]+)"', parts[8])
                gn_m  = re.search(r'gene_name "([^"]+)"', parts[8])
                gt_m  = re.search(r'gene_type "([^"]+)"', parts[8])
                if gid_m and gn_m:
                    sym = gn_m.group(1)
                    gid = gid_m.group(1).split(".")[0]
                    gtype = gt_m.group(1) if gt_m else "lncRNA"
                    lnc_dict[sym] = {"ensembl_id": gid, "gene_type": gtype}

    print(f"Total GENCODE v32 annotated lncRNAs: {len(lnc_dict)}")

    # 2. Identify 3,670 Stromal Barcodes
    obs_annot = ad.read_h5ad(os.path.join(DATA_DIR, "processed", "GSE248762_harmony_annotated_obs.h5ad")).obs
    stromal_mask = obs_annot["cell_type"] == "stromal / mesothelial-lineage (unresolved)"
    obs_st = obs_annot[stromal_mask].copy()
    stromal_barcodes = obs_st["_index"].values
    n_stromal_cells = len(stromal_barcodes)
    print(f"Total verified stromal cells: {n_stromal_cells}")

    # 3. Load QC AnnData and Subset to Stromal Cells
    print("Loading data/processed/GSE248762_hubblind_allcells_qc.h5ad...")
    adata_qc = ad.read_h5ad(os.path.join(DATA_DIR, "processed", "GSE248762_hubblind_allcells_qc.h5ad"))
    adata_st = adata_qc[stromal_barcodes, :].copy()
    adata_st.obs["donor_id"] = obs_st["donor_id"].values

    all_genes = list(adata_st.var_names)
    print(f"Assayed genes in single-cell matrix: {len(all_genes)}")

    # Identify which assayed genes are GENCODE v32 lncRNAs
    lnc_in_matrix = [g for g in all_genes if g in lnc_dict]
    print(f"Assayed GENCODE v32 lncRNAs in matrix: {len(lnc_in_matrix)}")

    # 4. Compute Single-Cell and Donor-Level Detection per Gene
    print("\nComputing single-cell expression frequencies and donor breadth...")
    X_st = adata_st.X.tocsr()
    
    # Gene cell counts
    cell_counts_arr = np.diff(X_st.tocsc().indptr)
    gene_to_cell_count = dict(zip(all_genes, cell_counts_arr))

    # Donor counts per gene
    donors = sorted(adata_st.obs["donor_id"].unique())
    n_total_donors = len(donors)
    print(f"Total unique donors represented: {n_total_donors}")

    gene_to_donor_count = {g: 0 for g in all_genes}
    for d in donors:
        d_mask = (adata_st.obs["donor_id"] == d).values
        X_d = X_st[d_mask, :]
        d_counts = np.diff(X_d.tocsc().indptr)
        for idx, g in enumerate(all_genes):
            if d_counts[idx] > 0:
                gene_to_donor_count[g] += 1

    # 5. Evaluate Multi-Threshold Sensitivity
    pct_cutoffs = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    donor_cutoffs = [1, 2, 3]

    sens_rows = []
    for pct in pct_cutoffs:
        min_cells = int(np.ceil(pct / 100.0 * n_stromal_cells))
        for dcut in donor_cutoffs:
            pass_cnt = sum(
                (gene_to_cell_count[g] >= min_cells) and (gene_to_donor_count[g] >= dcut)
                for g in lnc_in_matrix
            )
            sens_rows.append({
                "Min_Cell_Percentage_%": pct,
                "Min_Expressing_Cells": min_cells,
                "Min_Donors": dcut,
                "Passing_lncRNAs": pass_cnt,
                "Percentage_of_Assayed_lncRNAs": f"{pass_cnt / len(lnc_in_matrix) * 100:.2f}%"
            })

    df_sens = pd.DataFrame(sens_rows)
    print("\n" + "=" * 80)
    print("lncRNA DETECTION THRESHOLD SENSITIVITY TABLE")
    print("=" * 80)
    print(df_sens.to_string(index=False))

    sens_csv = os.path.join(OUT_DIR, "G_lncRNA_threshold_sensitivity.csv")
    df_sens.to_csv(sens_csv, index=False)
    print(f"\nSaved threshold sensitivity table: {sens_csv}")

    # 6. Apply Primary Threshold: >= 1.0% of stromal cells (>=37 cells) AND >= 2 donors
    primary_pct = 1.0
    primary_min_cells = int(np.ceil(primary_pct / 100.0 * n_stromal_cells))  # 37
    primary_min_donors = 2

    # Promiscuous/housekeeping/nuclear lncRNA blacklist
    promiscuous_names = {"MALAT1", "NEAT1", "XIST", "GAS5", "NORAD", "MEG3", "H19"}

    def get_credibility_status(symbol):
        if symbol == "XIST":
            return "EXCLUDED_FROM_INTERPRETATION_DONOR_SEX", "Sex-linked (X-chromosome inactivation); variation reflects donor sex, not peritoneal biology"
        if symbol in promiscuous_names or symbol.startswith("SNHG"):
            return "LOW_CREDIBILITY_CANDIDATE", "Nuclear/housekeeping/promiscuous transcript; ubiquitous high expression non-specifically sponges miRNAs"
        return "ROBUST_CANDIDATE", "Lineage-specific or moderate expression passing strict detection thresholds"

    detected_rows = []
    for g in lnc_in_matrix:
        n_c = gene_to_cell_count[g]
        pct_c = round((n_c / n_stromal_cells) * 100.0, 3)
        n_d = gene_to_donor_count[g]
        
        pass_primary = (n_c >= primary_min_cells) and (n_d >= primary_min_donors)
        cred_status, cred_reason = get_credibility_status(g)

        detected_rows.append({
            "symbol": g,
            "ensembl_id": lnc_dict[g]["ensembl_id"],
            "gene_type": lnc_dict[g]["gene_type"],
            "stromal_cells_expressing": n_c,
            "total_stromal_cells": n_stromal_cells,
            "pct_stromal_cells_expressing_%": pct_c,
            "donors_expressing": n_d,
            "total_donors": n_total_donors,
            "passes_primary_filter_ge1pct_ge2donors": pass_primary,
            "credibility_status": cred_status if pass_primary else "BELOW_DETECTION_THRESHOLD",
            "credibility_rationale": cred_reason if pass_primary else "Detected in < 1% of stromal cells or < 2 donors"
        })

    df_detected = pd.DataFrame(detected_rows).sort_values("stromal_cells_expressing", ascending=False)
    
    # Save full detected list
    det_csv = os.path.join(OUT_DIR, "G_lncRNA_detected_list.csv")
    df_detected.to_csv(det_csv, index=False)
    print(f"Saved full detected lncRNA list: {det_csv} ({len(df_detected)} rows)")

    # 7. Summary of Passing lncRNAs and Credibility Audit
    df_pass = df_detected[df_detected["passes_primary_filter_ge1pct_ge2donors"]].copy()
    print(f"\nTotal lncRNAs passing primary filter (>= 1.0% cells & >= 2 donors): {len(df_pass)}")
    print("\nCredibility Breakdown among Passing lncRNAs:")
    print(df_pass["credibility_status"].value_counts())

    print("\nPROMISCUOUS / LOW-CREDIBILITY lncRNAs (Flagged):")
    flagged = df_pass[df_pass["credibility_status"] != "ROBUST_CANDIDATE"]
    print(flagged[["symbol", "pct_stromal_cells_expressing_%", "donors_expressing", "credibility_status", "credibility_rationale"]].to_string(index=False))

    print("\nTOP 20 ROBUST CANDIDATE lncRNAs (by stromal expressing cell percentage):")
    robust = df_pass[df_pass["credibility_status"] == "ROBUST_CANDIDATE"]
    print(robust.head(20)[["symbol", "ensembl_id", "gene_type", "pct_stromal_cells_expressing_%", "donors_expressing"]].to_string(index=False))

    # Summary table
    summary_df = pd.DataFrame([
        {"Metric": "Total GENCODE v32 lncRNAs Annotated", "Count": len(lnc_dict)},
        {"Metric": "Assayed in GSE248762 Reference Matrix", "Count": len(lnc_in_matrix)},
        {"Metric": "Expressed in >= 1 Stromal Cell", "Count": sum(gene_to_cell_count[g] >= 1 for g in lnc_in_matrix)},
        {"Metric": "Passing Primary Filter (>= 1.0% cells & >= 2 donors)", "Count": len(df_pass)},
        {"Metric": "  - Robust Candidates", "Count": sum(df_pass["credibility_status"] == "ROBUST_CANDIDATE")},
        {"Metric": "  - Low-Credibility Candidates (MALAT1, NEAT1, GAS5, SNHG*)", "Count": sum(df_pass["credibility_status"] == "LOW_CREDIBILITY_CANDIDATE")},
        {"Metric": "  - Excluded from Interpretation (XIST, donor sex-linked)", "Count": sum(df_pass["credibility_status"] == "EXCLUDED_FROM_INTERPRETATION_DONOR_SEX")}
    ])
    sum_csv = os.path.join(OUT_DIR, "G_lncRNA_GSE248762_detection_summary.csv")
    summary_df.to_csv(sum_csv, index=False)
    print(f"\nSaved detection summary: {sum_csv}")

    # Provenance
    prov_df = pd.DataFrame([{
        "dataset": "GSE248762",
        "reference_annotation": "GENCODE v32 (GRCh38.p13, Ensembl 98)",
        "gtf_file": "gencode.v32.long_noncoding_RNAs.gtf.gz",
        "stromal_cells_assayed": n_stromal_cells,
        "stromal_donors_assayed": n_total_donors,
        "primary_filter_rule": ">= 1.0% of stromal cells (>=37 cells) AND >= 2 donors",
        "passing_lncrnas": len(df_pass),
        "robust_candidates": sum(df_pass["credibility_status"] == "ROBUST_CANDIDATE"),
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat()
    }])
    prov_df.to_csv(os.path.join(OUT_DIR, "G_lncRNA_provenance.csv"), index=False)
    print("Saved provenance.")

if __name__ == "__main__":
    main()
