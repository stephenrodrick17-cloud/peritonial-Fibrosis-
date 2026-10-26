# Seven-hub GeneMANIA, DSigDB, and KEGG analyses

All outputs from these analyses are confined to this directory. The selected
panel is the seven majority-vote diagnostic hubs in the project README and
tagged table: **FLT3LG, TNFSF15, LTB, ADAM19, SERPINA10, EBI3, and CXCL14**.
The task text also referred to five genes in its analysis steps; this run uses
the seven-gene panel explicitly selected by the user.

## Reproduction

Run from the project root with the project's Python interpreter:

```powershell
& ".venv\Scripts\python.exe" "results\revision\hub_gene_network_drug_gsea\01_genemania_network.py"
& ".venv\Scripts\python.exe" "results\revision\hub_gene_network_drug_gsea\02_dsigdb_candidates.py"
& ".venv\Scripts\python.exe" "results\revision\hub_gene_network_drug_gsea\03_kegg_single_gene_gsea.py"
```

The scripts require the already-installed `requests`, `pandas`, `numpy`,
`scipy`, `matplotlib`, `networkx`, and `gseapy` packages. Each script writes
only to this folder.

## Sources and settings

- **GeneMANIA:** Human (Homo sapiens; organism ID 4), live JSON service,
  queried 2026-10-13. Maximum of 20 additional genes; automatic/default
  network weighting. The returned interaction types and their weights are
  recorded in `genemania_interaction_type_weights.csv`; source dataset weights
  are recorded separately. The raw response is preserved.
- **Drug signatures:** Enrichr `DSigDB` library (DSigDB v1.0 as listed by
  Enrichr; library metadata and number of terms are recorded in
  `dsigdb_metadata.json`), queried 2026-10-13 for the earlier hub set. The
  validated seven-gene panel (TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, CST7)
  was queried separately on 2026-10-16; none of its 183 returned terms pass
  Benjamini-Hochberg adjusted P < 0.05 (minimum q = 0.1156). Terms are
  perturbation-signature associations, not treatment recommendations.
- **KEGG GSEA:** Human KEGG pathways obtained live from KEGG REST on
  2026-10-16. The pathway and gene database update dates and GSEA
  software/version are recorded in `validated7_kegg_gsea_metadata.json`.
  Each hub-gene analysis ranks all measured genes by Spearman correlation with
  that hub across the normalized GSE125498 discovery matrix (n=33; 13 LPD,
  20 SPD). GSEA uses 1,000 permutations, seed 20261013, gene-set sizes 15–500,
  and BH-FDR < 0.05 within each hub's pathway tests. Results are exploratory
  and hypothesis-generating.
- **Whole-cohort running-enrichment plot:** `plot_lpd_spd_kegg_running_gsea.py`
  reruns preranked GSEA using the saved unique-symbol limma ranking and the
  cached symbol-level KEGG sets. The 320 eligible pathways are tested with
  1,000 permutations (seed 20261016), gene-set sizes 15–500, and BH-FDR across
  all tested pathways. This is distinct from the earlier Entrez-based KEGG
  table and is used to plot complete pathway hit positions and running scores.

## Files

- `01_genemania_network.py`: query, full network edge and node tables, gene
  recognition, interaction-type/source-network weights, and PNG/PDF network.
- `02_dsigdb_candidates.py`: full Enrichr output, FDR-filtered candidate
  signatures, gene annotation coverage, and PNG/PDF summary plot.
- `plot_validated7_dsigdb_publication.py`: publication plot and ranked
  multi-gene nominal-term table from the frozen seven-gene DSigDB results;
  no external query is made.
- `05_seven_validated_hub_drug_gsea.py`: seven-hub DSigDB and per-gene KEGG
  GSEA workflow; produces per-gene Spearman rankings, pathway tables,
  `validated7_kegg_gsea_summary.csv`, and the seven-gene PNG/PDF GSEA panel.
- `plot_lpd_spd_kegg_running_gsea.py`: whole-cohort LPD-versus-SPD KEGG
  preranked GSEA, classic running-score figure, and full source tables.
- `plot_validated7_kegg_running_gsea.py`: classic running-enrichment plots for
  the seven hub-specific KEGG GSEA analyses, using their saved rankings and
  pathway results. Each hub displays up to three FDR-significant pathways
  selected by largest absolute NES.
- `validated7_kegg_running_gsea.png` / `.pdf`: one running-score, hit-rug, and
  hub-specific Spearman-rank panel per validated hub gene.
- `validated7_kegg_running_gsea_selected_pathways.csv`,
  `validated7_kegg_running_gsea_profiles.csv`, and
  `validated7_kegg_running_gsea_metadata.json`: plotted pathway statistics,
  per-rank curve data, and figure provenance.
- `lpd_spd_kegg_running_gsea.png` / `.pdf`: running-enrichment figure for the
  six strongest positive FDR-significant pathways.
- `lpd_spd_kegg_gsea_all_pathways.csv`,
  `lpd_spd_kegg_gsea_selected_pathways.csv`,
  `lpd_spd_kegg_gsea_ranked_genes.csv`, and
  `lpd_spd_kegg_gsea_running_profiles.csv`: complete test results, selected
  pathways, ranked genes, and per-rank curve data.
- `lpd_spd_kegg_running_gsea_metadata.json`: cohort, ranking, gene-set, and
  permutation metadata.
- `validated7_kegg_gsea_all_genes.csv`: combined pathway results and leading
  edges for all seven hub-specific analyses.
- JSON metadata files record query settings and versions for each data source.

## Limitations and deviations

- GeneMANIA interactors are database-derived functional associations, not
  differential-expression findings or diagnostic biomarkers.
- DSigDB enrichment identifies signatures associated with the input genes;
  it does not establish drug efficacy or treatment relevance to LPD.
- The discovery set is small, and the seven hubs were selected using discovery
  data. KEGG GSEA therefore has selection and small-sample limitations.
- The whole-cohort running-score plot is an exploratory analysis of the same
  discovery cohort; it is not external validation. Its symbol-based rerun and
  KEGG snapshot differ from the previously saved Entrez-based pathway table.
- The seven-hub running-score panels show within-cohort pathway associations
  with each hub's expression ranking. NES direction is not an LPD-versus-SPD
  pathway effect and should not be interpreted causally.
- Many genes have tied Spearman scores; GSEApy reports that tie ordering may
  affect preranked results. This is recorded in the run metadata.
- Empirical zero P/FDR estimates are reported at the finite resolution
  `1/(permutations + 1)` and flagged in the exported pathway tables.
- R and the requested `fgsea`/`clusterProfiler` packages were unavailable in
  the project environment. The single-gene preranked analysis uses GSEApy and
  reports this implementation deviation explicitly.
