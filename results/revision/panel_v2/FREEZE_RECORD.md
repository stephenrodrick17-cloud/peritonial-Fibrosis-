# Panel v2 freeze record

This record documents the panel definition frozen before the GSE92455 test.
It is not committed.

## Freeze identity

- Freeze commit: `26da2cfc081b64be841686f8ee2cf4e6e1014a82`
- Commit date: `2026-10-11T11:55:58+05:30`
- Frozen file: `results/revision/panel_v2/panel_v2.json`
- SHA-256 of the frozen file's exact Git blob bytes: `5baa81daff4c45d68adb9f73583bbbb5c9b40ba7f4a985531244fdc3074b1599`
- Audit script at the freeze: `scripts/task6_serpina10_detectability_audit.py`
- SHA-256 of the audit script's exact Git blob bytes at the freeze: `1f47d4162ab91a6f67306d0ad4080de0bd1e95a166f7811457e9b4de6a1b1c3e`
- Input hash snapshot date: `2026-10-11`

The frozen panel hash was calculated from the blob in commit `26da2cf`, not
from the working-tree file. Input hashes below are SHA-256 values of the files
currently present in the working tree on the snapshot date. For inputs marked
untracked, no committed copy is available to verify. A current hash is not
evidence that the same bytes were read when the committed panel outputs were
generated.

## Runtime inputs read by the audit script

The list covers the explicit inputs required or opened by
`task6_serpina10_detectability_audit.py`, including the cohort tables used to
derive detectability, direction rules, saved verdicts, and exploratory score
comparisons.

| Working-tree status | SHA-256 (current file) | Input path |
|---|---|---|
| **Untracked** | `2d1325c9eaa148b8c0ce21f777b30f6c523f584b493398ba37f0288dd0506ace` | `results/tables/GSE125498_all_results.csv` |
| Tracked, clean | `ae8e0b1923b118b28eae202a89e576cc970d8ad063533ff93a9b7952cbee9d10` | `results/revision/gse62928_validation_8genes/probe_selection_and_detection.csv` |
| Tracked, clean | `5528c34495f666ce30a60fe3bc9edcee4cae2b689cc92da6537373a7e7aecb94` | `results/revision/gse62928_validation_8genes/sample_group_mapping.csv` |
| Tracked, clean | `37f2378e37d06126f4c3f63ed300e90231f38f66ccba2e1412d1a9d837db7aa5` | `results/revision/gse62928_validation_8genes/GSE62928_series_matrix.txt.gz` |
| Tracked, clean | `fc4e450bec6bdaeb112bd540d809160f8fb6679313de177ca71e02aca1ea131a` | `results/revision/gse62928_validation_8genes/per_gene_primary_results.csv` |
| Tracked, clean | `09da374c32a1331da417d4596667cbcf3971a06aec56fdd49f6e987a90a18526` | `results/revision/gse62928_validation_8genes/validation_summary.json` |
| **Untracked** | `828e25004e3a542b2869631d1ed8732ac2c69627d215a60090a818c1396d5203` | `results/revision/exploratory_posthoc_20261009/detected_gene_composites_summary.json` |
| **Untracked** | `bfe06782d3d4e95b6913eb6697c225ac9cfccb3aff05eb7677737482b068a397` | `results/revision/effluent_validation/GSE248762/analysis_20261009/all_cell_pseudobulk_counts.csv` |
| **Untracked** | `55132c691f07dab4810c28dcd8d27b9d53e65a6a35aaa4fac2d08ba3cb82f38c` | `results/revision/effluent_validation/GSE248762/analysis_20261009/all_cell_pseudobulk_metadata.csv` |
| Tracked, clean | `fc0458f67a46bf7a18c5d441cf95caed828a56f37c7fb73d32851c0e6c23c821` | `results/revision/effluent_validation/GSE248762/analysis_20261009/gene_id_to_symbol.csv` |
| Tracked, clean | `8431071cc6c41b812a31eb928765c062e63e0d87af16286269e4eedbbceec404` | `results/revision/effluent_validation/GSE248762/analysis_20261009/primary_per_gene_results.csv` |
| **Untracked** | `a1b43f5e42d1386523cd7ae98b992e2368b0f647a762feb239c1b4b747ca22b8` | `results/revision/effluent_validation/GSE248762/analysis_20261009/pseudobulk_gene_detection.csv` |
| **Untracked** | `9db8febdc6b7caa073b181fe63f6c7006c7cf5f85b7f24330d86bea60ed3cceb` | `results/revision/effluent_validation/GSE248762/analysis_20261009/edger_logCPM.csv` |
| **Untracked** | `590238acd5d7ea87db4b4aa241511ac94ad7e6da6729ba71375594a52450b037` | `results/revision/effluent_validation/GSE248762/analysis_20261009/edger_all_gene_results.csv` |
| Tracked, modified | `3a7a50788afdf27cface7b1c85f020cbba98bb4e6cdfb22345fa6aa88b2429d9` | `results/revision/effluent_validation/GSE248762/analysis_20261009/gene_set_summary.json` |
| **Untracked** | `6dece1f9169a38e103a4335c06858dd8b2cd3d1cec1a85fb5bd6c468efbcc8c8` | `results/revision/effluent_validation/GSE130888/analysis_20261009/all_cell_pseudobulk_counts.csv` |
| **Untracked** | `133e6f14aad6c769da6c3cea6bef4755e331021a0b1fdc535639340e7f5e59f5` | `results/revision/effluent_validation/GSE130888/analysis_20261009/all_cell_pseudobulk_metadata.csv` |
| **Untracked** | `eb1b6338e3a1b6585ad7c5d34bacfe5adc225550c16a721a0527499f6a441e22` | `results/revision/effluent_validation/GSE130888/analysis_20261009/gene_id_to_symbol.csv` |
| **Untracked** | `3a0e6283a93134e2ae502c9dd393ca9e9eeb5c0e939b981964b749aa5baefd68` | `results/revision/effluent_validation/GSE130888/analysis_20261009/primary_per_gene_results.csv` |
| **Untracked** | `5488e38fb70a86cd1fd8a3e18aff96b9144c27d16a0c7748567d370b825a7331` | `results/revision/effluent_validation/GSE130888/analysis_20261009/pseudobulk_gene_detection.csv` |
| **Untracked** | `d305f3f3d16345a61b82d038486be146bffef8932017ac33324db19450704506` | `results/revision/effluent_validation/GSE130888/analysis_20261009/edger_logCPM.csv` |
| **Untracked** | `c0535a55cdf43c2d089f924535eb8ebd058b16178ff2b0f2f6cf44417e2a9497` | `results/revision/effluent_validation/GSE130888/analysis_20261009/edger_all_gene_results.csv` |
| **Untracked** | `c5b675ad7af1164de79d8881d2918e71eddffeda139b74751f6b8de8e006c044` | `results/revision/effluent_validation/GSE130888/analysis_20261009/gene_set_summary.json` |

The GSE248762 summary is marked modified because its verdict explanation was
edited in the working tree. The audit reads the saved verdict from this file;
the edited `verdict_rule_note` field is not used to calculate the panel
selection or score. Nevertheless, the hash above is the current modified-file
hash, and the exact summary bytes used in the original panel-v2 run cannot be
verified from this snapshot alone.
