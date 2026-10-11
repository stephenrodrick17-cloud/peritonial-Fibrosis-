# Panel detectability audit and exploratory panel v2

**Generated:** 2026-10-11

## Detectability audit

The audit reads saved cohort tables and the deposited GSE62928 series matrix; it does not rerun the preregistered analyses.

| Cohort | Panel genes detected | Detection basis |
|---|---:|---|
| GSE62928 | 0/8 | Highest-mean primary probe mean log2 expression > 5 across 8 samples |
| GSE248762 | 7/8 | Raw pseudobulk count ≥10 in at least half the samples in both LV and SV |
| GSE130888 | 6/8 | Raw pseudobulk count ≥10 in at least half the samples in both LV and SV |

### SERPINA10

- **GSE62928:** detected=False; filter=below expression cutoff; annotation present=True; mean=1.30856 (primary-probe RMA log2 expression); samples above cutoff=0/8.
- **GSE248762:** detected=False; filter=below expression cutoff; annotation present=True; mean=0.125 (raw pseudobulk counts per sample); samples above cutoff=0/16.
- **GSE130888:** detected=False; filter=zero counts; annotation present=True; mean=0 (raw pseudobulk counts per sample); samples above cutoff=0/10.

The saved effluent annotations contain SERPINA10 in both cohorts. Its GSE248762 raw counts are nonzero in the saved pseudobulk matrix but remain below the count threshold; GSE130888's saved raw pseudobulk counts are all zero. Thus it is not missing from annotation, but its undetected status has cohort-specific causes.

See `panel_detectability_audit.csv` for all genes, cohort-level means, sample counts above cutoff, and per-group counts. The biopsy cutoff is an operational reporting threshold, not a biological detection limit.
