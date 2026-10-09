# GSE248762 analysis deviations and post-run corrections

This log records script edits made after the first full Scanpy-to-edgeR
invocation. That invocation completed QC, pseudobulk preparation, marker-based
cell annotation, and UMAP generation, but edgeR stopped while constructing
confidence intervals. It did not produce a complete result/report. Numerical
statements below distinguish corrections that affect estimates from fixes
that only enable or clarify output. The saved final outputs were regenerated
with the corrected scripts.

## Scanpy and reporting script

| Edit | Reason | Changed a number? |
|---|---|---|
| Corrected the expected-direction comparison to compare validation and discovery log2FC signs, with discovery values still read from `GSE125498_all_results.csv`. | Implement the preregistered direction-concordance definition, rather than assuming a positive discovery sign in code. | No for the current panel: all eight discovery log2FCs are positive, and the seven detected genes have positive primary log2FCs. |
| Loaded the saved cell-type profile inside the reporting function. | The first post-edgeR report attempt raised `NameError` because `profile` was local to the preprocessing function. | No effect on calculated statistics; it enabled report generation. |
| Reworked composition-adjusted result assembly to left-join the fixed panel against model results, and suppress statistics for genes failing detection. | Empty adjusted output had no `gene_id` column; after models became estimable, an undetected gene also needed to remain “not detected,” rather than receiving an adjusted estimate. | No change to edgeR fits for detected genes. It changes the reported table for undetected SERPINA10; its estimate and P are now withheld. |
| Added contrast-specific detection checks and masked undetected genes in the secondary contrasts. | Apply the frozen detection rule independently to LV_UF vs LV_NOT_UF and LV_UF vs SV. | No change to model estimates/P values for genes passing detection; it changes which entries are reported as estimates versus not detected. |
| Added the report’s per-patient audit table and ECM1 leave-one-out summary, plus clarifying limitations and software-version reporting. | Include the requested raw-count audit, influence analysis, unresolved independence, and methods explanation. | No change to primary, secondary, or adjusted model statistics. The new audit numbers are in saved patient-level tables. |

The initial draft’s unused detection helper, undefined sample-group reference,
invalid unused `logCPM` indexing, and empty-score ROC call were corrected before
the first full invocation, so they are not post-run edits and did not contribute
to the results in this report.

## edgeR script and statistical outputs

| Edit | Reason | Changed a number? |
|---|---|---|
| Replaced nonexistent `fit$var.post`/unadjusted residual degrees of freedom with the edgeR quasi-likelihood fields `fit$s2.post` and `fit$df.residual.adj` for the approximate coefficient covariance intervals. | The first full invocation stopped with a non-conformable covariance calculation: edgeR 4.4.2 exposes `s2.post` and adjusted residual degrees of freedom for this fit. The interval calculation is documented in the generated report and uses the log2 scale. | Changes confidence intervals only; it does not change edgeR log2FC or P values. |
| Retained the fixed panel genes in edgeR filtering and, for the primary contrast, retained genes meeting the preregistered raw-count detection rule in both groups. | Ensure the fixed panel and genes needed for matched-set eligibility are not removed by generic `filterByExpr` filtering before the requested analysis. | Yes. The tested feature set used for dispersion estimation is different; final edgeR estimates/P values and downstream outputs are from this corrected model run. |
| Applied edgeR TMM normalization before writing patient-level log-CPM values. | Use normalized log-CPM for the preregistered composite and expression-decile matching, rather than unnormalized CPM-derived values. | It changes saved log-CPM values and can change composite scores/AUC and matched-null deciles/P values. It does not change the fitted edgeR log2FC/P values. |
| Read composition-column names without syntactic rewriting, and kept composition-model feasibility status separate from per-cell-type detection status. | R had converted slash-containing cell-type column names, so composition covariates were not found; a reused `status` variable also risked reporting the last gene’s status as model feasibility. | Enables the previously missing adjusted model and its reported estimates/P values; it does not change the unadjusted primary model. |
| Wrote a headered empty adjusted-results table when adjustment was not feasible, and now fail explicitly if a declared-feasible adjusted fit returns no results. | Keep empty outputs readable by pandas and avoid treating an absent model as a successful adjusted analysis. | No effect on this dataset’s fitted results; its adjusted model was feasible and fitted. |

## Additional saved patient-level analyses

Two scripts were added after the initial analysis to answer the requested
influence/count audit: `gse248762_patient_expression_audit.py` and
`gse248762_ecm1_leave_one_out.R`. The first reapplies the unchanged QC and
marker rules to the real GEO count archive and saves per-patient raw
pseudobulk counts and cells with nonzero counts. The second refits the
preregistered edgeR primary model after omitting each patient in turn. These
scripts add new patient-level audit results; they do not replace or alter the
primary analysis estimates.

## Preregistered verdict interpretation

No success criterion, gene, group, contrast, or verdict was changed after
examining the results. The direction component of the supportive tier is met
(7/7 detected genes concordant), but the required full-panel composite AUC and
its matched-null P-value cannot be calculated when SERPINA10 is not detected.
The partial-support branches are exactly five concordant genes, or an
estimable AUC from 0.65 through 0.80 without null significance; neither branch
applies here. The literal fallback therefore remains **Not supported**. This
decision rule is non-monotonic: 7/7 concordant can receive a lower label than
exactly 5 concordant when the composite is unavailable. The frozen criteria
were not changed post hoc.

The verdict implementation now checks the exact-five partial-support branch
before returning a verdict for an unestimable composite. This fixes the code
path for other possible data states; it does not change this cohort's verdict
or any statistical result.
