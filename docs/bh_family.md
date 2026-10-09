# Note on the Stage 4 Benjamini-Hochberg Multiple Testing Family (m = 98)

## 1. Why the Family Expanded from 96 to 98 Tests
In the initial implementation of Stage 4 edgeR pseudobulk modeling (`scripts/50_stage4_final.R`), the minimum evidence threshold was implemented using a global, 3-group conjunction:
```R
pass_evidence <- function(gene_name) {
  if (!gene_name %in% rownames(mat)) return(FALSE)
  row_counts <- mat[gene_name, ]
  uf_expr <- sum(row_counts[meta$group == "LV_UF"] > 0)
  not_uf_expr <- sum(row_counts[meta$group == "LV_NOT_UF"] > 0)
  sv_expr <- sum(row_counts[meta$group == "SV"] > 0)
  return(uf_expr >= 2 && not_uf_expr >= 2 && sv_expr >= 2)
}
```
Under this rule, any gene with fewer than 2 expressing donors in the `LV_UF` cohort was marked as `INSUFFICIENT DATA` across all three contrasts, including `LV_NOT_UF_vs_SV`.

In stromal cells (`stromal / mesothelial-lineage (unresolved)`), both `COL11A1` and `ISM1` have only 1 expressing donor in `LV_UF` (`LV_UF-1 = 1.0` and `3.0` counts, respectively). Consequently, both genes were marked as `INSUFFICIENT DATA` for all contrasts in `scripts/50_stage4_final.R`, yielding $m = 96$ PASS tests.

However, the contrast `LV_NOT_UF_vs_SV` evaluates differential expression exclusively between Long Vintage Non-UF (`LV_NOT_UF`, $n=6$) and Short Vintage (`SV`, $n=6$). It does not test the `LV_UF` group. When the evidence threshold is applied on a contrast-specific basis ($\ge 2$ expressing donors in each of the two groups tested in that contrast):
- `COL11A1` has 2 expressing donors in `LV_NOT_UF` (donors 1 and 2) and 4 in `SV` (donors 1, 2, 3, 4).
- `ISM1` has 3 expressing donors in `LV_NOT_UF` (donors 1, 2, 6) and 3 in `SV` (donors 1, 2, 3).

Both satisfy the per-contrast criterion of $\ge 2$ expressing donors per group. Applying this contrast-specific rule adds exactly 2 tests to the PASS family, bringing the total from $m = 96$ to $m = 98$.

## 2. Which Two Tests Were Added
The two tests added to the family are:
1. **Cell Type**: `stromal / mesothelial-lineage (unresolved)`
   **Contrast**: `LV_NOT_UF_vs_SV`
   **Gene**: `COL11A1`
   - $\log_2\text{FC} = -2.094542$
   - $\text{SE} = 0.887224$
   - $P = 0.039649$
   - $\text{BH\_FDR} = 0.258736$
   - Status: `PASS` (Previously `INSUFFICIENT DATA`)

2. **Cell Type**: `stromal / mesothelial-lineage (unresolved)`
   **Contrast**: `LV_NOT_UF_vs_SV`
   **Gene**: `ISM1`
   - $\log_2\text{FC} = 1.846394$
   - $\text{SE} = 0.778848$
   - $P = 0.127852$
   - $\text{BH\_FDR} = 0.389324$
   - Status: `PASS` (Previously `INSUFFICIENT DATA`)

Neither test achieves statistical significance under FDR control ($\text{FDR} < 0.05$ or $\text{FDR} < 0.10$).

## 3. Internal Frozen Analysis Plan Status of the Family Definition
- In `provenance/stage4_plan.json` — frozen prior to unsealing hub results in `scripts/35_stage4_preflight.py` — the multiple testing procedure was documented as an analysis plan (self-documented; committed after primary results):
  `"Benjamini-Hochberg FDR across all gene x cell type x contrast tests that pass evidence rule"`.
- The exact integer family size ($m$) was not part of an external registry; rather, it was defined within an internal frozen analysis plan, because the number of passing tests depended on data-driven filter criteria (donor-level expression).
- The transition from the 3-group global filter ($m = 96$) to the contrast-specific filter ($m = 98$) was an audit correction finalized after seeing primary results: evaluating expression in `LV_UF` for a contrast comparing only `LV_NOT_UF` and `SV` was recognized as an illogical coupling.

## 4. Whether Sensitivity Models Use the Same Family
- **No.** The sensitivity models (analyses `(b)` through `(j)` and the 16 Leave-One-Donor-Out runs, comprising 792 total rows in `stage4_sensitivity_edger_pseudobulk.csv`) do not use this $m = 98$ BH family, nor do they define a pooled 792-test family.
- All 792 sensitivity rows have `BH_FDR = NaN`.
- Per the robustness rule in the internal frozen analysis plan (`provenance/stage4_plan.json`):
  - Criterion 1 assesses primary significance ($\text{Primary FDR} < 0.05$).
  - Criteria 2 and 3 assess sign concordance ($\text{sign}(\log_2\text{FC})$ unchanged across sensitivity models and LODO iterations).
  Sensitivity analyses are thus evaluated on effect-direction robustness and nominal P-values, not on a secondary or joint FDR-adjusted family.
