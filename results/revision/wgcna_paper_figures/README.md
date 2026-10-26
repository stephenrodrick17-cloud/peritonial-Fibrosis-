# WGCNA paper figure panels

The figures were generated from the saved GSE125498 WGCNA workflow and checked
against the saved module assignments and module-trait statistics. Each figure
is available as a 600-dpi PNG and a vector PDF.

- `Fig3A_soft_threshold`: scale-free fit R², scale-free fit slope, and mean
  connectivity across tested powers. The selected power is marked in red;
  the dashed horizontal R² line is the prespecified 0.85 selection cutoff.
- `Fig3B_gene_dendrogram_modules`: gene dendrogram with merged module colors.
- `Fig3C_module_trait_associations`: module-eigengene correlations with PD duration.
- `Fig3D_matrisome_module_intersection`: the updated 46-to-43-to-39-to-3 gene overlap.
- `Fig3E_GO_BP_enrichment`: top 15 of the 37 BH-significant GO BP terms for the
  39 matrisome DEGs assigned to non-grey modules.

The dendrogram and module-association panels are based on the R WGCNA run,
not the separate legacy Python analysis. The 39-gene enrichment uses the full
5,001-gene WGCNA input as its background. The three-gene purple-module set
had no GO BP terms significant after BH correction. The network type was
signed, with soft-thresholding power selected as the lowest tested value
meeting R² >= 0.85; the scale-free fit slope is shown as a diagnostic, not
used as an additional selection criterion in the saved analysis.

Source data tables for the intersection and enrichment plots are in this
folder. To regenerate all panels from the repository root, run:

```powershell
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" --vanilla `
  "scripts\revision_make_wgcna_figures.R"
```
