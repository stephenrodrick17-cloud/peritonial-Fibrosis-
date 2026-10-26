# Network and single-cell verification summary

| Claim | Status | Verified value or disposition | Supporting file |
|---|---|---|---|
| Analysis gene set is the validated seven | TRUE | 7 genes loaded from the authoritative list | `results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv` |
| Frozen panel is a subset of the validated hubs | TRUE | 5 genes in the saved panel definition | `results/revision/panel_v2/panel_v2.json` |
| Candidate axes reconciled | TRUE | 323461 total = 1260 high + 322201 low; 14950 lncRNAs = 291 high-class + 14659 low-class | `network_reconciliation.csv` |
| Primary high-confidence network | TRUE | 1260 axes are primary; 322201 low-confidence axes are background only | `network_reconciliation.csv`; `network_axes.csv` |
| Tier A hub coverage | TRUE | Tier A miRNAs by hub are listed in NETWORK_WRITEUP.md; 2 hubs have one or more Tier A miRNAs | `hub_miRNA_tier_counts.csv` |
| miRNA expression in relevant PD samples | NOT VERIFIED | No sample-level abundance matrix was verified; the six ENA runs were not quantified | `miRNA_expression_support.csv`; `data/raw/network_7hub/GSE182736_ENA_run_metadata.tsv` |
| lncRNA pseudobulk association results | TRUE | 38446 testable pairs; see full table including non-significant and untestable pairs | `lncRNA_mRNA_correlation_results.csv` |
| lncRNA pseudobulk nominal/BH summary | TRUE | 4072 nominal P < 0.05; 175 BH q < 0.05 over all testable unadjusted pairs | `lncRNA_mRNA_correlation_summary.csv` |
| Cell-count-adjusted pseudobulk BH results | TRUE | 14 global BH q < 0.05; exploratory because lineage pseudobulk sample counts are small | `lncRNA_mRNA_correlation_summary.csv` |
| multiMiR source coverage | PARTIAL | 3 of 14 listed sources were not queried; predicted results contain exactly 1000 rows, possible truncation NOT VERIFIED | `database_source_versions.csv`; `multiMiR_predicted_raw.csv` |
| ENCORI ceRNA network | PARTIAL | 3 of 7 hub queries returned parsable rows; 4 returned API errors | `ENCORI_ceRNA_query_status.csv`; `ENCORI_ceRNA_pairs.csv` |
| Current miRTarBase release results | NOT RETRIEVED | Direct current-release result could not be obtained; older multiMiR source snapshot is labeled separately | `database_query_status.csv` |
| Single-cell edgeR model refit | NOT DONE | Saved seven-gene test results were reused; plots refreshed and forest plot added | `sc7_hub_gene_lineage_results.csv`; `sc7_nominal_forest_plot_data.csv` |

Predicted and CLIP-supported edges are candidate interactions. They are not described as experimentally validated regulation unless a qualifying functional assay is recorded.

The complete itemized file-read and file-written inventories are in `FILE_INVENTORY.md`; SHA-256 hashes and package/source metadata are in `run_metadata.json`.
