# FINAL MANUSCRIPT NUMBERS

| Section | Metric | Final value | Script | Output | Status |
| :--- | :--- | ---: | :--- | :--- | :--- |
| Discovery | Gene universe | 22,049 | 01_data_preprocessing.R | results/tables/GSE62928_all_results.csv | Verified |
| Discovery | Up-DEGs (Pro-fibrotic) | 367 | 02_differential_expression.R | results/tables/GSE62928_DEGs_filtered.csv | Verified |
| Discovery | Down-DEGs | 1,263 | 02_differential_expression.R | results/tables/GSE62928_DEGs_filtered.csv | Verified |
| Matrisome | Measured Matrisome | 975 | 03_matrisome_filtering.R | results/tables/GSE62928_ECM_intersection.csv | Verified |
| Matrisome | ECM-DEGs (Up) | 71 | 03_matrisome_filtering.R | results/tables/GSE62928_ECM_DEGs_candidate_list.csv | Verified |
| Matrisome | Enrichment | 4.375x | 03_matrisome_filtering.R | Console log / Pipeline lock | Verified |
| Matrisome | Hypergeometric P | 2.5929384e-26 | 03_matrisome_filtering.R | Console log / Pipeline lock | Verified |
| WGCNA | Input | 5,038 | 02b_wgcna_analysis.R | data/WGCNA_input.RData | Verified |
| WGCNA | Modules | 14 | 02b_wgcna_analysis.R | results/tables/WGCNA_module_trait_correlations.csv | Verified |
| WGCNA | Salmon | 604 | 02b_wgcna_analysis.R | results/tables/WGCNA_module_gene_counts.csv | Verified |
| WGCNA | Correlation | 0.806 | 02b_wgcna_analysis.R | results/tables/WGCNA_module_trait_correlations.csv | Verified |
| Convergence | Genes | 40 | 04_hub_gene_selection.R | results/tables/tripartite_candidate_genes.csv | Verified |
| ML | Consensus | 11 | 05_machine_learning_selection.py | results/tables/ml_consensus_hub_genes.csv | Verified |
| External cohort | Samples | 33 | 06_external_validation_GSE125498.py | results/tables/GSE125498_clinical_metadata.csv | Verified |
| ROC | Apparent AUC | 0.869 | 08b_roc_analysis.py | results/manuscript_numbers.json | Verified (In-sample) |
| ROC | Repeated CV AUC | 0.678 | 08b_roc_analysis.py | results/manuscript_numbers.json | Verified (Primary) |
| ROC | LOOCV AUC | 0.658 | 08b_roc_analysis.py | results/manuscript_numbers.json | Verified |
| Nomogram | C-index | 0.819 | 08_nomogram_roc_analysis.py | results/manuscript_numbers.json | Verified |
| Nomogram | CV AUC | 0.550 | 08_nomogram_roc_analysis.py | results/manuscript_numbers.json | Verified |
| GSEA | EMT NES | +3.203 | 07_pathway_gsea_analysis.R | results/tables/GSEA_hallmark_summary.csv | Verified |
| GSEA | EMT FDR | <0.001 | 07_pathway_gsea_analysis.R | results/tables/GSEA_hallmark_summary.csv | Verified |
