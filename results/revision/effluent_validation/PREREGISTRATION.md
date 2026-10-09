# Preregistration: Three GEO Effluent Cohort Evaluation

## Registration boundary

This plan is fixed before downloading or reading expression data from GSE248762,
GSE130888, or GSE92455. The accession metadata and file manifests may be inspected
to establish cohort identity and file formats; expression values will not be read
until this document has been committed. Any mismatch between the stated cohort
groups and GEO metadata, unavailable/download-failed data, or incompatible file
format will be reported and the affected dataset analysis paused for user
clarification. No gene, patient, cell type, contrast, QC/detection threshold, or
decision threshold will be changed after expression results are seen.

## Fixed gene panel and hypotheses

The eight genes are `TNFSF15, FLT3LG, EBI3, LTB, ADAM19, ECM1, SERPINA10, CST7`.
The strict 3/3 subset is `TNFSF15, FLT3LG`. The expected direction for each gene
is UP in long-term PD (LT) versus short-term PD (ST), as determined independently
from the discovery results CSV for GSE125498. Discovery log2FC and nominal P will
be read from that CSV at execution time, never copied from memory or README text.

## Cohorts and fixed contrasts

- **GSE248762 (primary):** combine LV_NOT_UF and LV_UF as long-vintage (LV) and
  compare LV against short-vintage (SV). Also report LV_UF versus LV_NOT_UF and
  LV_UF versus SV as secondary descriptive contrasts. No alternate grouping or
  pooling is allowed.
- **GSE130888 (secondary):** compare long-term PD against short-term PD. Normal
  controls are descriptive only and are excluded from the primary contrast.
- **GSE92455 (supplementary):** compare epithelioid against non-epithelioid
  effluent-derived mesothelial cells using the outcome-blind probe with highest
  mean expression across all samples. Omentum samples are descriptive only.
  This is a different phenotype contrast, will not be pooled with either scRNA
  cohort, and cannot establish vintage concordance.

## Acquisition, metadata, and independence

For each accession, record the GEO supplementary/processed file names, sizes,
and stated processing state; retain real GEO source files and record their
provenance. Save sample/patient/group metadata before analysis. Derive labels
from GEO metadata and check them against the cohort descriptions above; do not
infer or repair ambiguous group labels. Compare titles, characteristics,
submitter/contributor names, and PubMed IDs across GSE125498, GSE130888, and
GSE248762. Classify independence as confirmed, unknown, or violated, documenting
the evidence. If overlap is suspected, exclude that dataset from the overall
independent-evidence verdict.

## Single-cell processing plan (GSE248762 and GSE130888)

- Use Scanpy for single-cell object handling if the GEO file format is compatible;
  record the exact Scanpy, AnnData, and Python versions. Use edgeR for raw-count
  pseudobulk differential expression and record its version. Never substitute
  normalized/log expression for raw counts in edgeR. If the available files do
  not provide usable cell-level raw counts and sample/patient assignments,
  pause and report the incompatibility.
- Apply one group-blind QC rule to every sample in both datasets: retain cells
  with at least 200 detected genes and at most 20% mitochondrial counts. Do not
  remove samples or add QC thresholds after examining outcomes. Record cell
  counts before and after QC for each sample. Use the authors' cell annotation
  when provided and its provenance is available; otherwise use the marker rules
  below, without reference to group labels.
- Fallback marker annotation (one label per cell where unambiguous): T/NK
  (`CD3D, CD3E, TRAC, NKG7, GNLY`); B (`MS4A1, CD79A, CD79B, CD74`);
  monocyte/macrophage (`LST1, TYROBP, FCER1G, CD14, LILRB1, C1QA, C1QB`);
  neutrophil (`CSF3R, FCGR3B, S100A8, S100A9, CXCR2`); mesothelial
  (`KRT8, KRT18, KRT19, MSLN, UPK3B, KRT7`); fibroblast
  (`COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA`); otherwise `other`.
  Assign a listed type only when its marker score is the unique highest score
  and at least one listed marker is detected; otherwise label `other`.
- Save the marker list and UMAP panels by cell type and metadata-derived group.
  UMAP is descriptive; it will not determine groups or cell annotations.
- For each of the fixed eight genes, report the fraction of QC-retained cells
  with nonzero raw counts and mean normalized expression by annotated cell type.
  Select each gene's most-expressing cell type by the highest mean expression
  across all samples, without group labels; ties are broken by the fixed cell
  type order above, with `other` last. Keep this cell type fixed for the
  gene-specific cell-type pseudobulk analyses.
- Aggregate integer raw counts across all QC-retained cells per patient for the
  all-cell pseudobulk. Also aggregate by patient x cell type only when at least
  20 QC-retained cells are present; save every excluded patient-cell-type
  combination and its cell count. Cells are never treated as independent
  replicates. Include patients in all-cell pseudobulk only when at least 20
  QC-retained cells are available; report exclusions.
- Use edgeR quasi-likelihood negative-binomial models for the prespecified
  pseudobulk contrasts. Report log2FC, 95% CI, nominal P, and BH-adjusted P
  across the fixed eight genes within each contrast. Report direction as
  concordant only when the gene is detected and the estimated contrast is
  positive. A not-detected gene is neither concordant nor discordant.
- For primary GSE248762 LV versus SV composition, report per-patient proportions
  of the annotated major cell types and compare each proportion between groups
  using two-sided Wilcoxon rank-sum tests; BH-adjust within each dataset across
  the tested cell types. Composition-adjusted gene models will include
  monocyte/macrophage and T/NK proportions only if each primary group has at
  least five independent patients and the design matrix is full rank. Otherwise
  state that adjustment was not feasible. Do not add other post hoc covariates.

## Microarray plan (GSE92455)

Use only the real GEO expression and platform annotation. Record whether values
are raw, log2, and/or normalized based on GEO processing metadata and inspect
file structure before choosing preprocessing. Map probes to symbols using
platform annotation. For each panel gene, choose the probe with the highest mean
expression across all samples, without group labels; report unmapped/missing
genes explicitly. Use limma for epithelioid versus non-epithelioid and report
log2FC, 95% CI, nominal P, and BH P across the eight genes. Do not treat this
as a vintage contrast, include it in the vintage gene-set tests, or count it
toward the overall verdict.

## Detection, gene-set, and resampling rules

- A gene is detected in a scRNA pseudobulk contrast only if its raw count is at
  least 10 in at least half of the independent patient pseudobulk samples in
  both compared groups. Apply this rule before testing and report it per gene.
  Genes failing the rule are reported as `not detected`; they are not assigned
  a direction and are excluded from the concordance denominator.
- For each primary contrast separately, perform a one-sided exact binomial
  sign test for detected, direction-concordant genes against probability 0.5;
  report the eight-gene panel and the two strict genes descriptively.
- Calculate each patient's per-gene log-CPM from the all-cell pseudobulk counts,
  z-score each gene across the compared patients, and take the unweighted mean
  across the eight fixed genes for the composite. Report Mann-Whitney P and
  AUC oriented so higher values indicate the long-vintage group. Do not replace
  not-detected genes or silently change the fixed eight-gene score; if any
  panel gene is not detected, report the full-panel composite as not estimable
  and separately label any detected-gene descriptive score as non-prespecified.
- Bootstrap patients within group with replacement for 10,000 replicates to
  calculate a percentile 95% CI for AUC. Use seed 42 for every bootstrap or
  permutation. Do not bootstrap individual cells.
- Draw 1,000 distinct eight-gene sets from genes meeting the same detection
  rule, matched to panel genes by decile of mean log-CPM across the compared
  patient pseudobulks. Compute the upper-tail plus-one empirical P for mean
  gene-level log2FC and composite AUC. Each draw samples one null gene per
  panel-gene decile and does not reuse genes within a set. If there are too few
  eligible genes in any required decile, report the null as not estimable; do
  not change the bins, threshold, or sampling scheme. Apply the same descriptive
  set-level analyses to the two strict genes, clearly marked secondary.
- For GSE92455 report the per-gene limma results only; do not apply the
  vintage-specific sign test, composite, or random-gene null.

## Fixed success criteria

Apply separately to the primary GSE248762 LV versus SV and GSE130888 LT versus
ST contrasts:

- **Supportive:** at least 6 of 8 detected genes concordant, composite AUC at
  least 0.80, matched-null empirical P below 0.05, and both strict genes
  concordant.
- **Partially supportive:** 5 of 8 detected genes concordant, or AUC from 0.65
  through 0.80 without significance.
- **Not supported:** otherwise.

The overall verdict is **Supportive** only if both primary contrasts are
Supportive and independence is confirmed for both. If exactly one primary
contrast is Supportive, label the combined evidence **Partially supportive**.
If neither is Supportive, report **Not supported**; suspected patient overlap
prevents an overall independent-evidence claim and will be stated explicitly.
GSE92455 is separate and never contributes to this verdict. If composition
adjustment reverses direction concordance, describe the relevant gene(s) as
“reflecting immune-cell abundance,” not as supported.

## Outputs and reporting

Write all new files under `results/revision/effluent_validation/`, with
accession-specific metadata and results in their respective subdirectories.
Include scripts, input manifests, QC and mapping tables, per-gene results,
figures, and a report generated from saved outputs. Any output number must be
computed by a script and written to a saved table or machine-readable summary
before it is included in the README or report. State actual package versions.
Discuss small patient sample sizes, suspected or unknown patient overlap,
single-cell dropout, distinct GSE92455 contrast, unavailable covariates, and
wide intervals as uncertainty rather than proof of absence.
