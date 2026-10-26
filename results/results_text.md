
### T1. Batch/QC

PCA, sample-distance heatmaps, and per-sample expression summaries were generated separately for GSE125498, GSE248762, and GSE130888. The cohorts were not merged, so no batch-correction analysis was performed. GSE125498 used the saved RMA log2 matrix; effluent cohorts used the saved edgeR TMM log2-CPM matrices.

### T2. ML selection

The saved 39-gene model results were summarized without rerunning feature selection. LASSO, SVM-RFE, and random forest had 11, 10, and 10 selected genes, respectively; the saved results contained 20 unique selections and nine >=2-method consensus genes. Coefficient-path, cross-validation, and SVM-RFE error curves could not be reproduced because their intermediate model objects were unavailable.

### T3. Enrichment and matrisome composition

The saved measured-matrisome background contained 455 expressed genes; applying |log2FC| >= 0.585 and nominal P < 0.05 recovered 46 DEGs (33 upregulated and 13 downregulated). The 39-gene WGCNA intersection and 46-DEG set were tested with GO, KEGG, and Disease Ontology; no tested term had BH-adjusted q < 0.05 across the reported enrichment results. The category table also summarized all nine consensus hubs. ReactomePA was unavailable, so Reactome enrichment was omitted; no whole-genome background was substituted.

### T4. Hub expression and cross-cohort effects

Hub expression was summarized in all three cohorts, and group comparisons used Wilcoxon tests with Benjamini–Hochberg adjustment across the measured panel within each cohort. Effect directions and adjusted FDRs were displayed by cohort, with an exploratory DerSimonian–Laird random-effects summary across available estimates. The small cohorts and unverified overlap between additional cohorts limit inference from cross-cohort pooling.

### T5. Hub interaction and co-expression

A hub-centered network was extracted from the saved STRING functional-association network at a combined score threshold of 0.4, with at most 30 nodes. Because a fresh first-shell neighbor list was not available in the saved project files, this was not represented as a newly expanded STRING network; functional associations do not establish direct binding. Hub-hub Spearman correlations were summarized within each cohort.

### T6. Diagnostic modelling

The nine-gene ridge model had an apparent discovery AUC of 0.965 (bootstrap 95% CI 0.896–1.000) and a model-fitting optimism-adjusted AUC of 0.936 after 1,000 bootstrap resamples. Nested leave-one-out and repeated stratified five-fold AUCs were 0.908 and 0.914, respectively; the top five bootstrap-stability genes were TNFSF15, ADAM19, NAT6, FLT3LG, ECM1. The reduced five-gene model had an apparent AUC of 0.962 (bootstrap 95% CI 0.888–1.000); additional-cohort scores used within-cohort scaling and remained exploratory. A point-based nomogram diagram, apparent five-bin calibration plot with exact binomial intervals, and decision curve versus treat-all and treat-none strategies were generated for the reduced model (GSE125498: 20 controls, 13 long-term PD cases). Calibration and decision-curve results were apparent, not bootstrap-corrected or externally validated. The diagram visualizes the ridge model but is not an rms nomogram; panel selection on GSE125498 means that model-performance estimates may retain selection optimism.

### T7. Single-hub pathway analysis

Each measurable hub defined above- and below-median groups in GSE125498, and all expressed genes were ranked by limma t-statistics. The focused cached Hallmark analysis identified 20 hub–pathway associations at BH-adjusted P < 0.05; leading themes included TNF-alpha Signaling via NF-kB, Inflammatory Response, IL-6/JAK/STAT3 Signaling, Epithelial Mesenchymal Transition. Hallmark scores were calculated as mean within-cohort gene z-scores; GO preranked analysis exceeded the execution budget and KEGG was omitted because it requires repeated online annotation retrieval. These analyses used the discovery cohort and are post hoc; they do not constitute external pathway validation.

### T8. Immune infiltration

The available stromal score was summarized separately and was not treated as an immune-cell fraction. No two-method immune deconvolution was performed because compatible reference data and a second method were unavailable.

### T9. Drug-signature exploration

Public Enrichr queries produced no usable results in this run; no drug candidates were inferred. Drug perturbation signatures would be exploratory and would not establish therapeutic efficacy.

### T10. External-resource support

Human Protein Atlas data could not be retrieved for the consensus panel during this run. No tissue or single-cell expression values were imputed.

### T11. Summary

The workflow schematic, cohort characteristics, and integrated consensus-hub summary table were generated from saved analysis files. The final panel contained nine candidates selected by at least two of three methods; external-cohort evidence remains exploratory. Unavailable measures were left unfilled rather than inferred.


### Limitations

This analysis had small sample sizes in the additional effluent cohorts, and the discovery matrisome DEG screen used nominal P-values before filtering.
Machine-learning feature selection was performed in GSE125498; discovery-cohort ROC estimates are apparent and may be optimistic.
GSE248762 and GSE130888 are additional exploratory cohorts, but patient-level independence has not been established.
No wet-lab validation was performed; computational associations and expression-based discrimination do not establish causality or therapeutic efficacy.
