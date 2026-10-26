# Single-cell transcriptomic contextualization of the validated hub genes

## Methods

We conducted an exploratory cell-context analysis of GSE248762, a 10x
Genomics single-cell RNA-sequencing dataset from peritoneal dialysis effluent.
The analysis was restricted to the seven validated hub genes listed in
`results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`.
The five-gene panel defined in
`results/revision/panel_v2/panel_v2.json` is a subset of this validated set.
Additional genes present in older single-cell scripts and result files were
excluded from every calculation, table, and plot in this report.

Cells were retained when at least 200 genes were detected and mitochondrial
counts were no more than 20%. Broad lineages were assigned by maximum
normalized marker score, without using any of the seven hubs as annotation
markers. Eight marker-defined classes were collapsed to seven broad lineages:
lymphoid, myeloid, fibroblast/stromal, mesothelial, endothelial, epithelial,
and other/unresolved. The saved UMAP coordinates were reused; the three
coordinate entries identified as failing the recomputed mitochondrial filter
were excluded from the embedding plot.

The primary pseudobulk contrast compared long-vintage samples (LV; LV_UF and
LV_NOT_UF pooled) with short-vintage samples (SV). For a lineage to contribute
to an analysis, a sample had to contain at least 20 QC-retained cells. A
gene-lineage comparison was testable only when at least two eligible samples
in each group had at least 10 pseudobulk UMIs for the gene. Saved edgeR
quasi-likelihood statistics were used without refitting. BH q-values were
recalculated across the 23 finite P values for testable comparisons within the
validated-seven set. Positive log2 fold change denotes greater expression in
LV, based on the saved design and contrast in
`results/revision/single_cell_context/04_pseudobulk_edger.R`. The saved model
applies edgeR `filterByExpr` to matrix features and explicitly retains hubs
that pass the prespecified testability criteria.

The expression table reports the saved fraction of cells with raw count
greater than zero and mean cell-level log1p(CPM10k). It also reports
descriptive mean raw library-size CPM calculated from sample-lineage
pseudobulk count totals. These CPM values are not TMM-normalized model
estimates. Group means use samples meeting the 20-cell threshold for the
corresponding lineage.

## Results

The audit included 118,895 cells before QC and 110,353 after QC. The sample
groups were SV (6), LV_NOT_UF (6), and LV_UF (4); the primary contrast pooled
the latter two groups (LV, 10 samples) against SV (6 samples). The marker
dictionary contained none of the seven validated hubs. All seven hubs were
testable in at least one lineage.

There were seven broad lineages, giving 49 possible gene-lineage comparisons.
Twenty-three met the sample and gene-count criteria; 26 were not testable.
Six comparisons had nominal P < 0.05, and three remained significant after
BH correction across the 23 testable validated-seven comparisons. All three
BH-significant estimates were positive and arose in myeloid, mesothelial, and
epithelial lineages for the same hub. Their exact log2 fold changes, P values,
q values, eligible sample counts, and group mean CPM are in
`sc7_hub_gene_lineage_results.csv`. The other tested results, including all
nonsignificant comparisons, are in that same table. The saved statistics do
not contain confidence intervals, so none are reported.

At the BH-significance level, this one hub is the only recurrent signal across
lineages in this tested set. This does not establish a general or
independently replicated effect. The full descriptive table for the frozen
five-gene panel shows lineage-specific detected-cell fractions and mean
expression (`sc7_hub_gene_expression.csv`). In particular, the panel's
cell-level profile should not be interpreted as a T-cell-only profile:
the annotation groups T/NK and B marker classes together as lymphoid. The
separate Figure 7 analysis reported positive panel-score/T-cell
marker-enrichment correlations in GSE125498 (rho =
0.8571428571428571, BH q = 0.002352143153428128), GSE130888 (rho =
0.9393939393939393, BH q = 0.0007129268898067769), and GSE248762 (rho =
0.9558823529411766, BH q = 1.0088356800093399e-07). These are score-level
associations and do not by themselves establish cell-type specificity
(`results/revision/fig7_final_summary/FIG7_FULL_FINDINGS.md`).

Endothelial sample-lineages did not meet the 20-cell eligibility criterion.
Fibroblast/stromal eligibility was limited to four LV and five SV samples;
mesothelial eligibility was eight LV and six SV samples. Per-lineage sample
counts, eligible-cell totals, and within-sample composition summaries are
provided in `sc7_lineage_sample_cell_summary.csv` and
`sc7_lineage_composition.csv`. For context, tested pairs from the three
nonzero lineages with the lowest combined eligible-cell counts are flagged in
the results table (other/unresolved, 567 cells; fibroblast/stromal, 1,155;
mesothelial, 2,707). This relative flag does not change the prespecified
testability rules.

### Paper-ready results paragraph

In GSE248762 effluent pseudobulk data, the seven validated hub genes were
assessed in 10 LV and 6 SV samples. Of 23 testable gene-lineage comparisons,
6 were nominally significant and 3 remained significant after BH correction
within this seven-gene test family. The only hub with FDR-significant
differences had higher LV expression in myeloid (log2FC =
6.12396276741692, P = 0.0042269275581539, q =
0.0324064446125132), mesothelial (log2FC = 8.307668764101621, P =
0.0009508834036293, q = 0.011338783220178), and epithelial lineages
(log2FC = 6.60299329314642, P = 0.0009859811495807, q =
0.011338783220178). Confidence intervals were not present in the saved
results. These exploratory findings are limited by small sample sizes,
unequal lineage representation, and the absence of covariate adjustment.

## Limitations

This is an exploratory analysis of one effluent dataset, not an independent
validation of a clinical model. The primary contrast has 10 LV and 6 SV
samples, and some lineages have fewer eligible sample-level replicates; cells
were not treated as independent biological replicates. Pooling LV_UF and
LV_NOT_UF may conceal differences between those groups. The saved edgeR
design includes group only and does not adjust for clinical covariates or
batch. The saved output does not include confidence intervals or dispersion
diagnostics, so those cannot be supplied or independently reviewed here.
Broad marker-score annotation is not a fine-grained T-cell annotation, and
lineage composition may differ among samples. Finally, the analysis
contextualizes expression; it does not establish causality, mechanism, or
clinical utility.

## Figure legends

**Figure A. UMAP by broad lineage.** Saved GSE248762 UMAP coordinates are
colored by the seven broad marker-score lineages. The three cells identified
as failing the recomputed QC filter are omitted.

**Figure B. Validated-hub expression dot plot.** Rows show the seven validated
hub genes and columns show broad lineages. Dot area is proportional to the
fraction of cells with raw counts greater than zero; color represents saved
mean log1p(CPM10k).

**Figure C. Pseudobulk log2 fold change heatmap.** Tiles show saved edgeR
log2FC for the LV-versus-SV comparison. Grey tiles indicate untestable
gene-lineage pairs; a star indicates BH q < 0.05 among testable comparisons
for the validated seven.

**Figure D. Pseudobulk effect and nominal-P plot.** Each point is a testable
gene-lineage comparison, plotted by log2FC and negative log10 nominal P.
Nominal-only and BH-significant points are distinguished. The plot does not
show confidence intervals.

**Figure E. Lineage composition by original sample group.** Stacked bars show
the mean within-sample lineage percentages for SV, LV_NOT_UF, and LV_UF,
weighting samples equally.


## Refreshed output figures

The five lineage-context figures are reused from the validated-seven analysis and copied into this output folder as `fig_sc_umap_by_lineage`, `fig_sc_hub_gene_dotplot`, `fig_sc_logfc_heatmap`, `fig_sc_hub_gene_volcano`, and `fig_sc_lineage_composition` (PNG and PDF). A new forest plot of nominal comparisons, with FDR-significant points highlighted, is saved as `fig_sc_forest_nominal_and_fdr_hits` (PNG and PDF); its plotted values are in `sc7_nominal_forest_plot_data.csv`.
