# GSE130888 analysis deviations and corrections

The analysis used the preregistered genes, QC thresholds, long-term versus
short-term contrast, edgeR method, detection rule, seed, and cell annotations.
No sample, gene, threshold, or contrast was changed after examining expression.

| Correction | Reason | Changed fitted estimates or P values? |
|---|---|---|
| Excluded the three normal controls from inferential edgeR contrasts and retained them only in the saved descriptive expression table and plots. | The preregistration specifies long-term versus short-term PD as the primary contrast and normal controls as descriptive only. | No change to the primary contrast; prevents an unregistered control comparison. |
| Treated an empty composition-adjusted edgeR result as missing numeric estimates. | Composition adjustment is infeasible with four long-term samples under the preregistered minimum of five per group. The empty CSV initially caused a pandas type error in downstream reporting. | No. The unadjusted edgeR fit and per-gene results were unchanged. |
| Added a resume option to reuse the completed QC, pseudobulk, and edgeR files after the downstream reporting error. | Avoided rereading the 10x archive and regenerating group-blind preprocessing outputs during the correction. | No. The same saved real-count pseudobulk and edgeR outputs were used. |
| Applied the literal five-concordant-gene branch of the preregistered partial-support criterion even when the eight-gene composite was not estimable. | An initial implementation returned “Not supported” whenever the composite was unavailable, which incorrectly bypassed the separate “5 of 8 concordant” partial-support criterion. | No fitted estimate, P value, concordance, or detection status changed. The verdict was corrected from the interim “Not supported” to **“Partially supportive”** because five detected genes were concordant. |
| Added a discovery-versus-validation forest figure after confirming the discovery results CSV contains log2FC and P values but no confidence-interval columns. | Preserve the requested effect comparison without inventing or reconstructing discovery intervals that are absent from the authoritative source table. | No. Discovery points are shown without intervals; validation intervals are from edgeR. |

The final report and machine-readable summary were regenerated after these
corrections. The fixed full-panel composite and matched null remain
non-estimable because TNFSF15 and SERPINA10 fail the preregistered detection
rule; no replacement panel score was substituted.
