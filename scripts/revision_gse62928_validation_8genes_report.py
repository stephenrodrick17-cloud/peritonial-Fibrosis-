import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "revision" / "gse62928_validation_8genes"


def fnum(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}f}"


summary = json.loads((OUT / "validation_summary.json").read_text(encoding="utf-8"))
per_gene = pd.read_csv(OUT / "per_gene_primary_results.csv")
sets = pd.read_csv(OUT / "gene_set_level_results.csv")
nulls = pd.read_csv(OUT / "random_gene_set_null_summary.csv")
probes = pd.read_csv(OUT / "probe_selection_and_detection.csv")
pca = pd.read_csv(OUT / "qc_pca_scores.csv")
mapping = pd.read_csv(OUT / "sample_group_mapping.csv")
sensitivity = pd.read_csv(OUT / "sensitivity_highest_iqr_results.csv")

primary = sets.loc[sets["Set"] == "Eight-gene panel"].iloc[0]
strict = sets.loc[sets["Set"] == "Strict 3/3 genes"].iloc[0]
primary_null = nulls.loc[nulls["Set"] == "Eight-gene panel"].iloc[0]
strict_null = nulls.loc[nulls["Set"] == "Strict 3/3 genes"].iloc[0]
marker = summary["marker_score"]
power = summary["power"]
qc = summary["qc"]
concordant_n = int(primary["Concordant_Genes"])
auc = float(primary["Composite_AUC"])
empirical_p = float(primary["Null_AUC_Empirical_P"])
if concordant_n >= 6 and auc >= 0.80 and empirical_p < 0.05:
    verdict = "Supportive"
elif concordant_n == 5 or (0.65 <= auc <= 0.80 and empirical_p >= 0.05):
    verdict = "Partially supportive"
else:
    verdict = "Not supported"

gene_rows = []
for _, row in per_gene.iterrows():
    gene_rows.append(
        "| {gene} | {probe} | {disc_fc} ({disc_p}) | {val_fc} [{lo}, {hi}] | "
        "{limma_p} | {mw_p} | {d} [{dlo}, {dhi}] | {bh} | {direction} |".format(
            gene=row["Gene"],
            probe=row["Probe_ID"],
            disc_fc=fnum(row["Discovery_log2FC"]),
            disc_p=f"{row['Discovery_P']:.3g}",
            val_fc=fnum(row["Validation_log2FC"]),
            lo=fnum(row["Validation_CI95_Lower"]),
            hi=fnum(row["Validation_CI95_Upper"]),
            limma_p=f"{row['Validation_Limma_P']:.4g}",
            mw_p=f"{row['Mann_Whitney_exact_two_sided_P']:.4g}",
            d=fnum(row["Cohens_d"]),
            dlo=fnum(row["Cohens_d_bootstrap_CI95_Lower"]),
            dhi=fnum(row["Cohens_d_bootstrap_CI95_Upper"]),
            bh=f"{row['Limma_BH_8_Genes']:.4g}",
            direction="Yes" if row["Direction_Concordant"] else "No",
        )
    )

mapping_rows = [
    f"| {row.sample_id} | {row.GEO_title} | {row.subgroup} | {row.subgroup_assignment_basis} |"
    for row in mapping.itertuples(index=False)
]
sensitivity_rows = []
for gene in per_gene["Gene"]:
    row = sensitivity.loc[sensitivity["Gene"] == gene].iloc[0]
    sensitivity_rows.append(
        f"| {gene} | {row['Probe_ID']} | {fnum(row['logFC'])} | {fnum(row['CI.L'])} to "
        f"{fnum(row['CI.R'])} | {fnum(row['P.Value'])} |"
    )

flagged = qc["pca_flagged_samples"]
outlier_text = ", ".join(flagged) if flagged else "None flagged by the preregistered descriptive 2-SD PC1/PC2 screen"
detected_genes = probes.loc[probes["Detected_Mean_Above_5"], "Gene"].tolist()
not_detected = probes.loc[~probes["Detected_Mean_Above_5"], "Gene"].tolist()

report = f"""# GSE62928 Preregistered Eight-Gene Validation Report

## Registration and data provenance

The fixed panel, expected direction, tests, probe rules, and decision criteria were preregistered in [PREREGISTRATION.md](PREREGISTRATION.md) before GSE62928 expression values were read. The validation matrix was downloaded from GEO to [GSE62928_series_matrix.txt.gz](GSE62928_series_matrix.txt.gz). It identifies platform GPL13158, Affymetrix HG-U133 Plus 2.0. GEO states that values were processed with Robust Multi-array Average (RMA), including background adjustment, quantile normalization, and median-polish summarization. The deposited RMA matrix is treated as log2-scale; no further normalization was applied.

### GEO sample-to-group mapping

| Sample ID | GEO title | Analysis group | Group-label basis |
|---|---|---|---|
{chr(10).join(mapping_rows)}

GEO titles and treatment characteristics identify four EPS biopsies, two PD samples, and two uremic controls. The analysis treats the EPS and PD samples as severe EPS and early PD, respectively, following the supplied cohort description; GEO's sample-level labels themselves do not specify severity or the early-PD stage. The primary comparison combines the two control subgroups (4 EPS versus 4 non-EPS). No samples were excluded.

## Probe selection and QC

For each gene, the primary probe was selected by highest mean expression across all samples, without using group labels. The IQR-selected probe and all mapped probes are recorded in [probe_selection_and_detection.csv](probe_selection_and_detection.csv), [all_panel_probe_inventory.csv](all_panel_probe_inventory.csv), and [sensitivity_all_probes_results.csv](sensitivity_all_probes_results.csv). All {len(probes)} genes mapped to GPL13158. The operational detection threshold was mean log2 expression > {qc['primary_probe_detection_cutoff_mean_log2_gt']:.1f}; {qc['primary_panel_genes_detected']} of 8 panel genes exceeded it. The genes below this reporting threshold were {", ".join(not_detected)}. This threshold is an operational cutoff, not a manufacturer negative-control background estimate or a biological detection limit.

Per-sample distribution and PCA plots are [qc_expression_boxplot.pdf](qc_expression_boxplot.pdf) and [qc_pca_by_subgroup.pdf](qc_pca_by_subgroup.pdf); PCA uses the top {qc['pca_variable_probe_count']:,} variable probes without labels. {outlier_text}. This descriptive screen did not remove samples. PCA coordinates and flags are in [qc_pca_scores.csv](qc_pca_scores.csv).

## Per-gene primary-probe results

Discovery log2FC and nominal P are taken directly from `results/tables/GSE125498_all_results.csv`. Validation log2FC, confidence intervals, and moderated P-values are from limma. Mann-Whitney P-values are exact two-sided permutation/rank probabilities for 4 versus 4. Cohen's d is pooled-SD standardized EPS-minus-control effect; its CI is a percentile bootstrap from {summary['bootstrap_replicates']:,} stratified resamples with seed 42. BH adjustment is across the eight primary-probe limma tests.

| Gene | Primary probe | Discovery log2FC (P) | Validation log2FC (95% limma CI) | limma P | Exact Mann-Whitney P | Cohen's d (bootstrap 95% CI) | BH P (8 genes) | Concordant |
|---|---|---:|---:|---:|---:|---:|---:|---|
{chr(10).join(gene_rows)}

Direction concordance was {concordant_n}/8. The one-sided exact sign-test P-value against a 0.5 concordance probability was {fnum(primary['Sign_Test_One_Sided_P'])}. The minimum attainable exact two-sided Mann-Whitney P for 4 versus 4 without ties is 2/70 = 0.0286.

The highest-IQR probe sensitivity results are in [sensitivity_highest_iqr_results.csv](sensitivity_highest_iqr_results.csv). Under that probe rule, {int((sensitivity['logFC'] > 0).sum())}/8 validation effects were positive. All-probe results are reported separately; genes with multiple probes can show probe-dependent effect directions, so the primary outcome-blind probe rule is retained for the prespecified verdict.

## Gene-set tests

| Set | Concordant | Score difference (EPS - non-EPS) | Exact Mann-Whitney P | AUC (bootstrap 95% CI) | Null mean AUC | Matched-null AUC empirical P | Mean gene log2FC | Matched-null log2FC empirical P |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Eight-gene panel | {int(primary['Concordant_Genes'])}/8 | {fnum(primary['Composite_Mean_Difference_EPS_minus_NonEPS'])} | {fnum(primary['Composite_Mann_Whitney_Exact_Two_Sided_P'])} | {fnum(primary['Composite_AUC'], 3)} [{fnum(primary['Composite_AUC_Bootstrap_CI95_Lower'], 3)}, {fnum(primary['Composite_AUC_Bootstrap_CI95_Upper'], 3)}] | {fnum(primary_null['Null_AUC_Mean'])} | {fnum(primary['Null_AUC_Empirical_P'])} | {fnum(primary['Mean_Gene_log2FC'])} | {fnum(primary['Null_Mean_Gene_log2FC_Empirical_P'])} |
| Strict genes: TNFSF15, FLT3LG | {int(strict['Concordant_Genes'])}/2 | {fnum(strict['Composite_Mean_Difference_EPS_minus_NonEPS'])} | {fnum(strict['Composite_Mann_Whitney_Exact_Two_Sided_P'])} | {fnum(strict['Composite_AUC'], 3)} [{fnum(strict['Composite_AUC_Bootstrap_CI95_Lower'], 3)}, {fnum(strict['Composite_AUC_Bootstrap_CI95_Upper'], 3)}] | {fnum(strict_null['Null_AUC_Mean'])} | {fnum(strict['Null_AUC_Empirical_P'])} | {fnum(strict['Mean_Gene_log2FC'])} | {fnum(strict['Null_Mean_Gene_log2FC_Empirical_P'])} |

For each fixed set, {int(primary_null['Null_Sets'])} random sets were sampled from genes above the operational detection threshold, matched to the panel's mean-expression decile using the detected-gene universe. Null sampling was without replacement within each set, with seed 42; empirical P-values use the plus-one correction. All eight fixed panel genes fall below the operational detection threshold and below the detected-gene expression range; they are therefore assigned to the lowest decile for this coarse matching. This materially limits how closely the detected-gene null can match their absolute abundance and is reported as a limitation, not hidden.

The score is the unweighted mean of the per-gene z-scores, where each gene was standardized across all eight samples. ROC and matched-null distributions are shown in [composite_roc_and_null.pdf](composite_roc_and_null.pdf); sample-level scores are in [composite_scores_by_sample.csv](composite_scores_by_sample.csv).

## Marker-score proxy, subgroup description, and power

The leukocyte marker-score proxy used {marker['genes_used_n']} mapped genes from the discovery Task 4 `LEUKOCYTE_MARKERS` list, excluding the fixed panel. It is a marker-expression score, not xCell/CIBERSORT or cell-fraction deconvolution. The EPS-minus-control score difference was {fnum(marker['eps_minus_control'])} (exact two-sided Mann-Whitney P = {fnum(marker['mann_whitney_exact_two_sided_p'])}); its Pearson correlation with the eight-gene composite was r = {fnum(marker['pearson_r_with_8_gene_composite'])} (P = {fnum(marker['pearson_p_with_8_gene_composite'])}). Scores are in [leukocyte_marker_score_proxy.csv](leukocyte_marker_score_proxy.csv), and genes/probes used in [leukocyte_marker_genes_used.csv](leukocyte_marker_genes_used.csv).

Expression for EPS, early PD, and uremic subgroups is displayed descriptively (no subgroup hypothesis tests; each control subgroup has n = 2) in [gene_expression_subgroups_descriptive.pdf](gene_expression_subgroups_descriptive.pdf). The minimum detectable standardized effect at 80% power is Cohen's d = {fnum(power['minimum_detectable_cohens_d'], 4)} for a two-sample equal-variance independent t-test with n = 4 per group, two-sided alpha = 0.05, and 6 degrees of freedom.

## Preregistered verdict

**{verdict}.** The fixed criteria require at least 6/8 direction-concordant genes, composite AUC >= 0.80, and matched-null AUC empirical P < 0.05 for “Supportive.” The observed panel meets the AUC and AUC-null conditions but has only {concordant_n}/8 concordant genes, so it does not meet the conjunction. It also does not meet “Partially supportive” (which requires 5/8 concordant, or AUC from 0.65 to 0.80 without null significance). The secondary strict two-gene set is not supportive.

## Limitations

- The validation comparison has only 4 EPS and 4 controls. The power calculation shows that only very large standardized effects are detectable with the stated test; wide confidence intervals indicate uncertainty, not proof of no effect.
- Discovery used dialysate effluent from long- versus short-term PD; validation uses peritoneal biopsy tissue from severe EPS versus early-PD and uremic controls. Tissue, phenotype, and comparator differences limit direct transportability.
- No age, sex, batch, or other covariate adjustment was possible in this eight-sample analysis.
- The fixed panel genes are all below the prespecified operational mean-expression threshold of 5.0. The threshold is not a platform negative-control limit, but the low means and bottom-decile-only random-gene matching constrain interpretation.
- The continuous composite AUC remains calculable from the deposited expression values, but its 0.875 estimate with a 0.500–1.000 bootstrap interval is highly uncertain. Being below the operational cutoff alone neither establishes noise-dominated signal nor provides evidence of reliable discrimination.
- The eight-gene composite has a broad bootstrap interval, and probe sensitivities show that selecting a different annotated probe can change an individual gene's effect direction.
- This is an exploratory cohort evaluation, not evidence of a validated diagnostic test or causal mechanism.

## Reproducibility

Run `Rscript scripts/revision_gse62928_validation_8genes.R <series_matrix.gz> <GPL13158_annotation> results/revision/gse62928_validation_8genes` from the repository root, then run `python scripts/revision_gse62928_validation_8genes_report.py`. Analysis code: [revision_gse62928_validation_8genes.R](../../../scripts/revision_gse62928_validation_8genes.R). Software versions and machine-readable results are in [validation_summary.json](validation_summary.json).
"""

(OUT / "gse62928_validation_report.md").write_text(report, encoding="utf-8")
print(f"Wrote {OUT / 'gse62928_validation_report.md'}")
print(f"Verdict: {verdict}")
