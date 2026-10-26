# GSE248762 single-cell hub-gene contextualization

## Scope and source

This is an exploratory cellular-contextualization analysis, not a validation of
the diagnostic models. It uses the public GEO series [GSE248762](
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248762), a 10x Genomics
single-cell RNA-seq dataset from overnight peritoneal dialysis effluent.
The GEO sample metadata contain 16 samples: 6 `SV`, 6 `LV_NOT_UF`, and 4
`LV_UF`. These counts sum to 16. The pseudobulk contrast combines the two
long-vintage categories (`LV_NOT_UF` and `LV_UF`; n=10) and compares them with
`SV` (n=6). Participant identifiers are not reported; GSM samples are treated
as sample-level replicates, not asserted to be unique donors.

The raw GEO archive and series-matrix metadata were already cached in the
project. Their local file timestamps are 2026-10-09; an independent original
download log was not available. Analysis run date: 2026-10-13. The archive
was not re-downloaded. The starting count object is
`data/processed/GSE248762_hubblind_allcells_qc.h5ad`; it contains 118,895
cells before this analysis's filter. Reapplying the stated filter retained
110,353 cells. A prior saved run reports 110,356 cells after QC; three cells
in its UMAP coordinates exceed 20% mitochondrial counts in the current count
object and were excluded here. The three cell IDs and values are listed in
`umap_cells_excluded_by_recomputed_qc.csv`.

## Hub panels

The hub panel follows the shared analysis code's `HUBS` vector in
`scripts/_common.R`: **TNFSF15, FLT3LG, EBI3, LTB, ADAM19, NAT6, ECM1, CST7,
and CCBE1** (9 genes). Panel membership is exported in
`hub_gene_panel_membership.csv`. Five are in the five-gene model panel
(TNFSF15, ADAM19, FLT3LG, NAT6, ECM1); five are in the separate seven-gene
consensus panel (FLT3LG, TNFSF15, LTB, ADAM19, EBI3).

## QC, annotation, and embedding

Cells were retained if they had at least 200 detected genes and no more than
20% mitochondrial counts. The input count matrix is raw UMI count data.
Expression used for marker scoring and display is library-size normalized to
10,000 counts per cell and transformed as `log1p`. No sample integration or
batch correction was performed. For visualization, this analysis reuses the
saved joint UMAP coordinates generated from a group-blind set of the top 2,000
features most often selected as per-sample highly variable genes, followed by
PCA (up to 40 components), a 15-neighbor graph, and UMAP (seed 42). No Leiden
clustering was used for the present analysis; broad cell labels were assigned
by maximum mean normalized marker expression, with cells lacking a unique
positive marker-set winner labelled `Other/unresolved`. Hub genes were not
used to annotate cells.

Marker panels are exported in `lineage_marker_genes_and_sources.csv`. They
cover T/NK and B lymphoid cells; monocyte/macrophage and neutrophil myeloid
cells; fibroblast/stromal cells; mesothelial cells; and endothelial and
epithelial cells. Marker genes are literature-curated reference markers
cross-referenced to CellMarker 2.0 (Hu C, et al. *CellMarker 2.0: an updated
database of manually curated cell markers in human/mouse and web tools based
on scRNA-seq data.* *Nucleic Acids Research*. 2023;51(D1):D870-D876.
doi:10.1093/nar/gkac947; PMID:36300619). These are broad annotation markers,
not effluent-validated classifiers. The cell-level assignments and sample by
lineage counts are in `cell_type_annotation_table.csv` and
`sample_lineage_cell_counts.csv`.

## Expression and pseudobulk tests

For every hub gene and lineage, the profile table reports the fraction of
lineage cells with raw count >0 and mean `log1p` normalized expression.
Pseudobulk counts are raw UMI sums per sample and lineage. A sample-lineage
unit was eligible for testing when it contained at least 20 QC-retained cells.
A gene-lineage contrast was testable only when at least two eligible samples
in each contrast group had at least 10 pseudobulk UMIs for that gene. Both
groups needed at least two eligible samples. Untestable combinations are
retained in the result table with reasons.

edgeR 4.4.2 was used with TMM normalization and a robust quasi-likelihood
negative-binomial model. The design included the two groups only, with no
additional covariates; `LV` is the numerator, so positive log2 fold changes
indicate higher expression in long-vintage samples. Benjamini-Hochberg
correction was applied once across the complete set of testable
gene-lineage-contrast P values. The small number of samples per group limits
power; cells were not treated as independent replicates.

## Reproduction

Run from the project root in the configured environment:

```powershell
& ".venv\Scripts\python.exe" "results\revision\single_cell_context\01_audit_dataset.py"
& ".venv\Scripts\python.exe" "results\revision\single_cell_context\02_annotate_lineages.py"
& ".venv\Scripts\python.exe" "results\revision\single_cell_context\03_expression_and_pseudobulk_counts.py"
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" "results\revision\single_cell_context\04_pseudobulk_edger.R" "results\revision\single_cell_context"
& ".venv\Scripts\python.exe" "results\revision\single_cell_context\05_summarize_and_figures.py"
```

Software versions: Python 3.13.1, Scanpy 1.12.4, AnnData 0.13.4, pandas 2.3.3,
NumPy 2.4.6, SciPy 1.18.1, Matplotlib 3.11.2, seaborn 0.13.2, R 4.4.2, and
edgeR 4.4.2.

## Outputs

- `dataset_audit.json`, `sample_group_qc_audit.csv`, and
  `umap_cells_excluded_by_recomputed_qc.csv`: accession/group/cell QC audit.
- `hub_gene_panel_membership.csv`: hub-panel membership flags.
- `lineage_marker_genes_and_sources.csv`, `cell_type_annotation_table.csv`,
  `sample_lineage_cell_counts.csv`, and `umap_coordinates_annotated.csv`:
  markers, lineage labels, sample counts, and coordinates.
- `hub_gene_expression_by_lineage.csv`, `hub_gene_feature_presence.csv`,
  `hub_gene_expression_per_cell.csv`, `gene_id_symbol_map.csv`, and
  `pseudobulk_counts/`: per-lineage profiles and sample-level raw counts.
- `edger_gene_lineage_results_raw.csv`, `gene_lineage_pseudobulk_tests.csv`,
  `hub_gene_testability_by_lineage.csv`, `hub_gene_testability_overall.csv`,
  `nominal_P_lt_0_05.csv`, and
  `global_FDR_lt_0_05.csv`: pseudobulk comparisons and global correction.
- `umap_by_lineage`, `umap_by_group`, `hub_gene_lineage_dotplot`, and
  `umap_hub_gene_feature_plots` are each exported as PNG and PDF.
- `results_paragraph.txt` and `analysis_summary.json`: output-grounded result
  text and numerical summary.

Results are an exploratory contextualization of the hub genes in a single
effluent dataset. They do not validate diagnostic models and do not establish
causality or mechanism.
