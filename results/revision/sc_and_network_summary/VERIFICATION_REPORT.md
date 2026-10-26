# Verification report: single-cell context and regulatory network

## Scope and evidence rules

The single-cell analysis below uses only the seven validated hub genes loaded
from `results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`.
The frozen five-gene panel is read from `results/revision/panel_v2/panel_v2.json`;
it is a subset of those seven. Older analysis outputs containing additional
genes were not used for counts, results, or plots. Saved raw edgeR
quasi-likelihood statistics were filtered to the seven genes and the
Benjamini-Hochberg (BH) adjustment was recalculated for the validated-seven
testable pairs. The edgeR model was not refit.

“NOT VERIFIED” means that a claim could not be established from the local
underlying files. It does not mean the biological effect or interaction is
absent. File paths below are relative to the repository root.

## Part A — GSE248762 single-cell contextualization

| # | Claim | Status | Verified value or finding | Supporting file(s) |
|---|---|---|---|---|
| A1 | Cell counts and QC | TRUE | The saved audit reports 118,895 cells before QC and 110,353 after recomputed QC. The thresholds are at least 200 detected genes and at most 20% mitochondrial counts, inclusive. | `results/revision/single_cell_context/dataset_audit.json`; `results/revision/single_cell_context/01_audit_dataset.py`; `sc7_sample_qc.csv` |
| A2 | Lineage annotation is hub-blind | TRUE | Lineages are assigned from marker scores. None of the seven validated hubs occurs in the marker dictionary parsed from the annotation script (zero overlaps). The eight marker classes are collapsed into seven broad lineages, including an “other/unresolved” category. | `results/revision/single_cell_context/02_annotate_lineages.py`; `sc7_lineage_marker_overlap_audit.csv` |
| A3 | Hub set and relation to frozen panel | TRUE | Seven validated hubs are read from the cited validated-seven source. The five genes in `panel_v2.json` are a subset of this seven-gene set. | `results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`; `results/revision/panel_v2/panel_v2.json`; `sc7_hub_gene_expression.csv` |
| A4 | Number of genes testable | TRUE | All seven validated hubs were testable in at least one lineage. The previous wider single-cell result set was excluded from this count. | `sc7_hub_gene_testability.csv`; `sc7_hub_gene_lineage_results.csv` |
| A5 | Possible and actual tests | TRUE | Seven broad lineages give 7 × 7 = 49 possible gene-lineage pairs. Of these, 23 met testability criteria and 26 did not. A sample-lineage needs at least 20 QC cells; each gene needs at least 10 pseudobulk UMIs in at least two eligible samples in each group; both groups need at least two eligible samples. edgeR `filterByExpr` is applied to the matrix features, with hubs passing the explicit testability rule retained for fitting. | `results/revision/single_cell_context/04_pseudobulk_edger.R`; `sc7_hub_gene_testability.csv`; `sc7_hub_gene_lineage_results.csv` |
| A6 | Nominal and BH-significant results | TRUE | Six of 23 tested pairs had nominal P < 0.05; three were BH-significant at q < 0.05. The BH family is all 23 finite P values from testable pairs among the validated seven—not the earlier broader gene set and not all 49 possible pairs. | `sc7_hub_gene_lineage_results.csv` |
| A7 | Direction of log fold change | TRUE | The edgeR design uses separate SV and LV columns and the contrast is LV − SV. A positive log2FC therefore means higher expression in LV. | `results/revision/single_cell_context/04_pseudobulk_edger.R` |
| A8 | Directional statements | TRUE, re-derived | Each tested pair and its direction is listed in `sc7_hub_gene_lineage_results.csv`. The six nominal results and the three FDR results are summarized below. Statements involving genes outside the validated seven were not carried forward. | `sc7_hub_gene_lineage_results.csv` |
| A9 | “Most consistent signal across lineages” | SUPPORTED WITHIN THIS TEST | Only one of the seven hubs has any BH-significant pair, and that signal occurs in three lineages. This is a narrow statement about the seven-gene, seven-lineage test family, not evidence of broad biological consistency or replication. | `sc7_hub_gene_lineage_results.csv` |
| A10 | Dispersion and fitting warnings | NOT VERIFIED | Saved files include software versions and results, but no edgeR warning log or dispersion diagnostic output was found. | `results/revision/single_cell_context/pseudobulk_software_versions.csv`; `results/revision/single_cell_context/04_pseudobulk_edger.R` |

### Nominal results, including BH q-values

Positive log2FC is higher in LV. The displayed values are carried from the
saved edgeR result and the recalculated validated-seven BH family; the full
23-pair table, including nonsignificant pairs, is
`sc7_hub_gene_lineage_results.csv`.

| Hub gene | Lineage | log2FC (LV vs SV) | P | BH q (23-pair family) | Interpretation |
|---|---|---:|---:|---:|---|
| EBI3 | lymphoid | 1.45543085082006 | 0.0224587289526643 | 0.086091794318546 | Nominal only |
| TNFSF15 | myeloid | 1.96310332575269 | 0.01234165296056 | 0.070964504523218 | Nominal only |
| EBI3 | myeloid | 1.41123739720539 | 0.0186121752291994 | 0.085616006054317 | Nominal only |
| ECM1 | myeloid | 6.12396276741692 | 0.0042269275581539 | 0.0324064446125132 | BH-significant |
| ECM1 | mesothelial | 8.307668764101621 | 0.0009508834036293 | 0.011338783220178 | BH-significant |
| ECM1 | epithelial | 6.60299329314642 | 0.0009859811495807 | 0.011338783220178 | BH-significant |

The other 17 tested pairs are nonsignificant at nominal P < 0.05; untestable
pairs are also retained with their reasons in `sc7_hub_gene_lineage_results.csv`.
The results CSV includes eligible sample counts by group, numbers of samples
meeting the gene-UMI criterion, eligible-cell totals, and descriptive group
mean raw library-size CPM. CPM is not the edgeR TMM-normalized model input.

### Extra checks

| Check | Verified result | Supporting file |
|---|---|---|
| Original groups | 6 SV, 6 LV_NOT_UF, and 4 LV_UF samples; the primary contrast pools LV_NOT_UF and LV_UF (LV n=10) versus SV (n=6). | `sc7_sample_qc.csv`; `results/revision/single_cell_context/dataset_audit.json` |
| Sparse lineages | No endothelial sample-lineage passed the 20-cell threshold. Fibroblast/stromal eligibility was 4 LV and 5 SV samples; mesothelial eligibility was 8 LV and 6 SV. The per-lineage/group cell counts, minima, and medians are in `sc7_lineage_sample_cell_summary.csv`. | `sc7_lineage_sample_cell_summary.csv` |
| Lower-cell test pairs | For a descriptive relative flag, the three nonzero lineages with the fewest total cells among samples meeting the 20-cell rule are other/unresolved (567 cells), fibroblast/stromal (1,155), and mesothelial (2,707). Tested pairs from these lineages are explicitly tagged in the results CSV; this is a relative ranking, not a new exclusion threshold. | `sc7_lineage_sample_cell_summary.csv`; `sc7_hub_gene_lineage_results.csv` |
| Composition | Composition is summarized as mean within-sample percentage by original sample group, giving samples equal weight rather than pooling cells. | `sc7_lineage_composition.csv` |
| Five-gene panel expression | The five frozen panel genes are profiled in every broad lineage in `sc7_hub_gene_expression.csv`; the table gives fraction of cells with raw counts > 0, saved cell-level mean log1p(CPM10k), and LV/SV mean raw library-size CPM. Endothelial group CPM is unavailable because no sample passed the lineage eligibility filter. | `sc7_hub_gene_expression.csv`; `results/revision/single_cell_context/hub_gene_expression_by_lineage.csv` |
| Earlier T-cell enrichment association | The separate Figure 7 report records a positive panel-score/T-cell-marker-enrichment association across its cohorts. This is a score-to-signature association, not proof of T-cell-specific expression; the single-cell broad “lymphoid” category combines B and T/NK marker classes. | `results/revision/fig7_final_summary/FIG7_FULL_FINDINGS.md`; `results/revision/single_cell_context/02_annotate_lineages.py` |

The saved Figure 7 T-cell rows are:

| Cohort | n | Spearman rho | 95% CI | P | BH q |
|---|---:|---:|---|---:|---:|
| GSE125498 | 13 | 0.8571428571428571 | [0.5210537736221249, 0.9774647887323943] | 0.00018093408872524061 | 0.002352143153428128 |
| GSE130888 | 10 | 0.9393939393939393 | [0.6226415094339623, 1.0] | 5.4840529985136686e-05 | 0.0007129268898067769 |
| GSE248762 | 16 | 0.9558823529411766 | [0.8179036432188833, 0.9969653282032328] | 7.760274461610307e-09 | 1.0088356800093399e-07 |

**Source:** `results/revision/fig7_final_summary/FIG7_FULL_FINDINGS.md`.

## Part B — lncRNA–miRNA–mRNA network

| # | Claim | Status | Verified value or finding | Supporting file(s) |
|---|---|---|---|---|
| B1 | Historical 11-gene set and relation to validated seven | TRUE | The historical STRING PPI workflow specifies 11 nodes. It is a different set from the validated seven; the saved lists overlap in five genes. It is a protein interaction analysis, not an lncRNA–miRNA–mRNA network. Names are intentionally not repeated here. | `scripts/create_legacy_11_hub_ppi.py`; `results/revision/legacy_11_hub_ppi/network_nodes.csv`; `hub_set_comparison.csv` |
| B2 | miRTarBase, multiMiR, and ENCORI/starBase sources, versions, dates, filters | NOT VERIFIED | No local result/source table, version/date record, or evidence-filter record for these sources was found. No database query was performed for this task. | Local file search; `results/revision/figures_final/figure_legends.md` |
| B3 | Tier-A functional evidence | NOT VERIFIED | Tier definitions, genes meeting Tier A, and miRNA counts per gene cannot be checked without the evidence tables. | `network_claim_verification.csv`; `network_evidence_by_hub.csv` |
| B4 | Axis, miRNA, lncRNA, confidence, and arithmetic totals | NOT VERIFIED | No axis table was found. The claimed component totals and any arithmetic mismatch therefore cannot be reconciled. | `network_claim_verification.csv`; `network_axis_evidence.csv`; `network_totals_verification.csv` |
| B5 | miRNA detection and low-interaction definitions | NOT VERIFIED | No operational definitions or supporting result tables were found. The local spreadsheet `data/raw/GSE182736_All_Expressed_miRNA.xlsx` is a 1,578-row by 36-column miRNA annotation/reference workbook, not a sample-level expression matrix. | `data/raw/GSE182736_All_Expressed_miRNA.xlsx`; `network_claim_verification.csv` |
| B6 | Single-cell lncRNA–mRNA correlations and detection rates | NOT VERIFIED | The dataset, correlation method, analysis unit, adjustment model, multiplicity procedure, lncRNA detection rates, and cited counts cannot be verified from saved outputs. No corresponding correlation or detection table was found. | `network_claim_verification.csv`; `network_axis_evidence.csv` |
| B7 | “Post-translational” wording | NOT VERIFIED | No local network writeup containing the phrase was located. If describing miRNA-mediated regulation, the relevant term is post-transcriptional. | Local text search; `results/revision/figures_final/figure_legends.md` |
| B8 | Evidence-supported retained axes | NOT VERIFIED | The number of axes retained under a three-component evidence requirement cannot be determined. This must not be reported as zero biological axes; the local evidence files are missing. | `network_axis_evidence.csv` |

The old PPI node count and its five-gene overlap are set comparisons only. No
PPI edge was treated as miRNA or lncRNA evidence. The seven-hub evidence table
marks unavailable counts as `NOT VERIFIED`, not zero.

## Figure inventory

All listed figures are exported as PNG at 300 dpi and as PDF.

| Files | Description |
|---|---|
| `figA_umap_lineage.png` / `.pdf` | Saved UMAP coordinates colored by the seven broad lineages; the three cells excluded by the recomputed QC audit are omitted. |
| `figB_hub_gene_dotplot.png` / `.pdf` | Seven validated hubs by lineage; dot size is the fraction of cells with raw counts > 0, color is saved mean log1p(CPM10k). |
| `figC_logfc_heatmap.png` / `.pdf` | Log2 fold changes for all seven hubs and seven lineages; untestable pairs are grey and BH-significant pairs are starred. |
| `figD_hub_gene_volcano.png` / `.pdf` | Nominal P versus effect size for the 23 testable pairs; nominal-only and BH-significant pairs are distinguished. No confidence intervals are implied. |
| `figE_lineage_composition.png` / `.pdf` | Equal-sample-weighted lineage composition by the three original sample groups. |
| `figF_regulatory_network_status.png` / `.pdf` | Explicit evidence-status panel; no network edges are inferred. |
| `figG_supported_mirna_status.png` / `.pdf` | Explicitly marks miRNA support counts as not verified. |
| `figH_evidence_tier_status.png` / `.pdf` | Explicitly marks evidence tiers as not verified; no Sankey is inferred. |
| `figI_lncrna_correlation_status.png` / `.pdf` | Explicitly marks the unavailable lncRNA correlation outputs. |

## Files read

The complete input manifest is also recorded in `run_metadata_summary.json`.
The files read for this report were:

- `results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`
- `results/revision/panel_v2/panel_v2.json`
- `results/revision/single_cell_context/README.md`
- `results/revision/single_cell_context/analysis_summary.json`
- `results/revision/single_cell_context/dataset_audit.json`
- `results/revision/single_cell_context/01_audit_dataset.py`
- `results/revision/single_cell_context/02_annotate_lineages.py`
- `results/revision/single_cell_context/03_expression_and_pseudobulk_counts.py`
- `results/revision/single_cell_context/04_pseudobulk_edger.R`
- `results/revision/single_cell_context/edger_gene_lineage_results_raw.csv`
- `results/revision/single_cell_context/global_FDR_lt_0_05.csv`
- `results/revision/single_cell_context/nominal_P_lt_0_05.csv`
- `results/revision/single_cell_context/hub_gene_testability_by_lineage.csv`
- `results/revision/single_cell_context/hub_gene_expression_by_lineage.csv`
- `results/revision/single_cell_context/hub_gene_expression_per_cell.csv`
- `results/revision/single_cell_context/hub_gene_feature_presence.csv`
- `results/revision/single_cell_context/lineage_marker_genes_and_sources.csv`
- `results/revision/single_cell_context/sample_group_qc_audit.csv`
- `results/revision/single_cell_context/sample_lineage_cell_counts.csv`
- `results/revision/single_cell_context/umap_coordinates_annotated.csv`
- `results/revision/single_cell_context/umap_cells_excluded_by_recomputed_qc.csv`
- `results/revision/single_cell_context/pseudobulk_software_versions.csv`
- `results/revision/single_cell_context/pseudobulk_counts/counts_lymphoid.csv`
- `results/revision/single_cell_context/pseudobulk_counts/counts_myeloid.csv`
- `results/revision/single_cell_context/pseudobulk_counts/counts_fibroblast_stromal.csv`
- `results/revision/single_cell_context/pseudobulk_counts/counts_mesothelial.csv`
- `results/revision/single_cell_context/pseudobulk_counts/counts_endothelial.csv`
- `results/revision/single_cell_context/pseudobulk_counts/counts_epithelial.csv`
- `results/revision/single_cell_context/pseudobulk_counts/counts_other_unresolved.csv`
- `results/revision/legacy_11_hub_functional_500/network_nodes.csv`
- `results/revision/legacy_11_hub_ppi/network_nodes.csv`
- `scripts/create_legacy_11_hub_ppi.py`
- `data/raw/GSE182736_All_Expressed_miRNA.xlsx`
- `results/revision/figures_final/figure_legends.md`
- `results/revision/figures_final/make_figures.py`
- `results/revision/fig7_final_summary/FIG7_FULL_FINDINGS.md`

`results/revision/figures_final/figure_legends.md` and the local text search
also record that no panel-v2 miRNA/lncRNA result table was available and no
database query was made.
