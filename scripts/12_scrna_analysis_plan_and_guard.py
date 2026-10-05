"""
Script 12: scRNA-seq Analysis Plan & Marker Panel Initialization (Step E0)
Timestamped Analysis Plan for GSE248762 (effluent scRNA-seq) and code-level hub-gene guard.
"""

import os
import json
import datetime

# CODE-LEVEL GUARD
BLOCKED_HUB_GENES = {
    "ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", 
    "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"
}

def check_hub_gene_guard(query_genes):
    """Raises RuntimeError if any of the 11 hub genes are queried in Stage 3."""
    if isinstance(query_genes, str):
        query_genes = [query_genes]
    violations = [g for g in query_genes if str(g).upper() in BLOCKED_HUB_GENES]
    if violations:
        raise RuntimeError(
            f"[CODE-LEVEL GUARD TRIGGERED]: Access to hub genes {violations} is strictly "
            f"prohibited during Stage 3 (hub-gene-blind QC, clustering, and cell-type annotation)!"
        )

def main():
    print("=================================================================")
    print("STEP E0: TIMESTAMPED SCRNA-SEQ ANALYSIS PLAN & MARKER PANEL")
    print("=================================================================")
    
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    os.makedirs("provenance", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)
    os.makedirs("results/tables", exist_ok=True)
    
    # 1. Marker Panel (Excluding all hub genes: NO COL3A1, NO VCAN)
    marker_panel = {
        "Mesothelial": ["WT1", "MSLN", "CALB2", "UPK3B", "KRT19"],
        "Fibroblast_Myofibroblast": ["DCN", "LUM", "PDGFRA", "COL1A1", "COL1A2", "ACTA2", "TAGLN"],
        "Macrophage": ["CD68", "CD163", "MRC1", "MARCO"],
        "Monocyte": ["CD14", "FCGR3A", "S100A8", "S100A9"],
        "Neutrophil": ["S100A8", "S100A9", "FCGR3B", "CSF3R"],
        "T_Cell": ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "TRAC"],
        "NK_Cell": ["NCAM1", "NKG7", "GNLY", "KLRD1"],
        "B_Cell": ["CD19", "MS4A1", "CD79A"],
        "Plasma_Cell": ["SDC1", "MZB1", "TNFRSF17", "IGHG1"],
        "Dendritic_Cell": ["HLA-DRA", "CD1C", "CLEC9A", "LILRA4"],
        "Mast_Cell": ["TPSAB1", "CPA3", "KIT"],
        "Endothelial": ["PECAM1", "VWF", "CDH5"]
    }
    
    # Verify guard on marker panel
    all_markers = [m for sublist in marker_panel.values() for m in sublist]
    check_hub_gene_guard(all_markers)
    print(f"Marker panel verified: {len(all_markers)} canonical lineage markers. (Zero hub genes included)")
    
    with open("provenance/marker_panel.json", "w") as f:
        json.dump(marker_panel, f, indent=2)
    print("Saved provenance/marker_panel.json")

    # 2. Complete Analysis Plan
    analysis_plan = {
        "dataset": "GSE248762 (Human peritoneal dialysis effluent scRNA-seq, 10x Genomics)",
        "cohort_size": "16 donors (6 Short Vintage, 6 Long Vintage No-UFF, 4 Long Vintage UFF)",
        "timestamp_plan_frozen_utc": timestamp,
        "hub_gene_blind_policy": "Strict code-level guard enforced during QC, normalization, clustering, and annotation",
        "qc_parameters": {
            "per_sample_adaptive_mad": "3 MAD from median on log1p(total_counts), log1p(n_genes_by_counts), and pct_counts_mt",
            "hard_floor_min_genes": 200,
            "hard_floor_min_umis": 500,
            "percent_mito_ceiling": 15.0,
            "min_cells_per_gene": 3
        },
        "doublet_detection": "Scrublet (scanpy.pp.scrublet or custom Scrublet simulation) per sample, expected doublet rate 0.06",
        "ambient_rna_policy": "Filtered matrices used; ambient RNA deconvolution not possible without raw droplet calls. Flagged that secreted ECM transcripts in effluent may reflect ambient RNA.",
        "normalization_and_hvg": {
            "target_sum": 10000,
            "transformation": "log1p",
            "n_top_genes": 2000,
            "flavor": "seurat / cell_ranger"
        },
        "dimensionality_reduction_and_integration": {
            "pca_components": 30,
            "batch_key": "donor_id / sample_id",
            "integration_algorithm": "Harmony / Scanorama across 16 donors"
        },
        "clustering_selection_rule": {
            "algorithm": "Leiden",
            "resolution_range": [0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
            "selection_criterion": "Maximum average silhouette score across tested resolutions; ties broken by lower resolution"
        },
        "annotation_rules": {
            "marker_panel_file": "provenance/marker_panel.json",
            "unassigned_rule": "Any cluster with ambiguous lineage enrichment labeled as 'Unassigned' without forced mapping"
        },
        "minimum_evidence_testability_rule": {
            "min_cells_per_donor": 20,
            "min_donors_in_target_group": 3,
            "min_donors_in_comparator_group": 2,
            "contrasts": ["LV_UF_vs_LV_NOT_UF", "LV_UF_vs_SV", "LV_NOT_UF_vs_SV"]
        }
    }
    
    with open("provenance/analysis_plan_scRNA.json", "w") as f:
        json.dump(analysis_plan, f, indent=2)
    print("Saved provenance/analysis_plan_scRNA.json (Frozen plan)")
    print("Step E0 finished successfully.")

if __name__ == "__main__":
    main()
