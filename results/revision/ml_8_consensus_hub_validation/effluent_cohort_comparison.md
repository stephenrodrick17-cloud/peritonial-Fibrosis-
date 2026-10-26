# Exploratory comparison of the revised 8-gene panel in two effluent cohorts

The panel tested was **TNFSF15, EBI3, PI3, FLT3LG, ADAM19, CLEC4F, TGM3,
CXCL14** from `results/final_analysis/ML_3_models_8_consensus_hub_genes.csv`.
The cohort pseudobulk counts, sample-level metadata, QC, and edgeR
quasi-likelihood approach were taken from each cohort's existing analysis.
The comparison is exploratory: these revised genes were not the fixed,
preregistered panel for either existing effluent-cohort analysis.

## Summary

| Dataset | Contrast | Patient pseudobulks | Detected panel genes | Same-direction among detected | Full 8-gene score |
|---|---|---:|---:|---:|---|
| GSE248762 | LV_UF + LV_NOT_UF vs SV | 10 LV, 6 SV | 7/8 | 6/7 (85.7%; sign-test P=0.0625) | Not estimable |
| GSE130888 | Long-term PD vs short-term PD | 4 LV, 6 SV | 5/8 | 4/5 (80.0%; sign-test P=0.1875) | Not estimable |

A gene was considered detected if its raw pseudobulk count was at least 10 in
at least half the patient samples in **each** group. Undetected genes were
not counted as concordant or discordant. The full-panel score was not
calculated for either cohort because one or more panel genes were undetected;
therefore there is no panel AUC or matched-null AUC to compare.

For the revised-panel re-fit, edgeR's `filterByExpr` features and genes
meeting the raw-count detection rule in both groups were retained, with all
eight target panel genes forced into the model so their detection status
could be assessed. Thus the test uses the existing cohort filtering approach
with the revised panel.

## Which cohort is the better follow-up?

**GSE248762 is the stronger exploratory check of this particular panel so
far:** it has more LV patient pseudobulks (10 vs 4), detected more genes
(7 vs 5), and had a larger direction-concordant fraction among detected genes
(6/7 vs 4/5). This is a relative comparison only, not confirmation: the
fixed eight-gene panel score could not be evaluated in either cohort.

**GSE130888 has the more explicit long-term-versus-short-term PD contrast.**
If matching the discovery study's disease-duration contrast is the priority,
it is relevant, but here it has fewer LV samples and lower panel-gene
detection. The two cohorts use effluent cells and do not directly replicate a
peritoneal-tissue biopsy cohort.

Neither cohort supports claiming that all eight genes have been validated.
The results are too small and incomplete for a definitive "best validation"
claim. A larger independent cohort with all eight genes measurable and a
prespecified composite test would be needed.

## Gene detection and direction

| Gene | GSE248762 | GSE130888 |
|---|---|---|
| TNFSF15 | Detected; concordant | Not detected |
| EBI3 | Detected; concordant | Detected; concordant |
| PI3 | Detected; concordant | Detected; concordant |
| FLT3LG | Detected; concordant | Detected; concordant |
| ADAM19 | Detected; concordant | Detected; concordant |
| CLEC4F | Detected; discordant | Detected; discordant |
| TGM3 | Not detected | Not detected |
| CXCL14 | Detected; concordant | Not detected |

In GSE248762, EBI3 and PI3 had nominal edgeR P<0.05 and remained below
BH-adjusted P=0.05 among detected panel genes; no full-panel claim should be
inferred from those individual-gene results. In GSE130888, no detected
individual gene had nominal P<0.05.

## Detailed outputs

- [GSE248762 panel summary](GSE248762/panel_validation_summary.csv)
- [GSE248762 per-gene results](GSE248762/per_gene_validation.csv)
- [GSE248762 sample groups](GSE248762/sample_group_mapping.csv)
- [GSE130888 panel summary](GSE130888/panel_validation_summary.csv)
- [GSE130888 per-gene results](GSE130888/per_gene_validation.csv)
- [GSE130888 sample groups](GSE130888/sample_group_mapping.csv)
- Analysis code: [`validate_ml_8_hubs_effluent_cohorts.R`](../../../scripts/validation/validate_ml_8_hubs_effluent_cohorts.R)
