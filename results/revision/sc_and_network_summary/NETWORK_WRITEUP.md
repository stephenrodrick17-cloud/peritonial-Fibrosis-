# lncRNA–miRNA–mRNA regulatory network: evidence audit

## Methods

The requested network was to be restricted to the seven validated hub genes
read from `results/revision/hub_gene_network_drug_gsea/validated7_kegg_gsea_summary.csv`.
The historical 11-node set was audited separately using its saved STRING PPI
source. The comparison records set sizes and overlap only; the historical set
was not substituted for the validated-seven set.

Local files under `results/`, `data/`, and `scripts/` were searched for
miRTarBase, multiMiR, ENCORI/starBase, lncRNA–miRNA–mRNA axes, evidence tiers,
miRNA expression, and single-cell lncRNA–mRNA correlation outputs. The local
miRNA workbook was inspected for structure and contents. No external database
queries or new interaction predictions were run because the required source
tables and evidence-filter definitions were not available locally.

## Results

The historical STRING PPI workflow specifies 11 nodes. Its set differs from
the validated-seven set, with an overlap of five genes
(`hub_set_comparison.csv`). A STRING protein-interaction node list is not
evidence for a miRNA-target or lncRNA-mediated regulatory axis.

No locally saved miRTarBase, multiMiR, or ENCORI/starBase result table was
found. Accordingly, evidence-source versions and retrieval dates, functional
versus non-functional validation filters, CLIP-seq stringency, minimum
support criteria, and Tier-A definitions cannot be verified. The local
workbook `data/raw/GSE182736_All_Expressed_miRNA.xlsx` contains miRNA
annotation/reference fields; it is not a sample-level expression matrix.
Thus it does not establish miRNA detection for a hub-associated network.

No underlying lncRNA–miRNA–mRNA axis table, lncRNA detection table, or
single-cell lncRNA–mRNA correlation output was found. The claimed counts of
putative axes, miRNAs, unique lncRNAs, detected miRNAs, low-interaction axes,
high-confidence axes, and retained three-component-supported axes are
therefore **NOT VERIFIED**. The draft's approximately 75 nominal positive
lncRNA–mRNA pairs and approximately 40 after donor-cell-number adjustment,
the dataset, correlation method, analysis unit, adjustment model,
multiple-testing procedure, and lncRNA detection rates are also **NOT
VERIFIED**; those draft counts are not reported as results. Without the source
tables, no arithmetic reconciliation is possible; this should not be
reported as evidence that the counts are zero. The complete count-status
table is `network_totals_verification.csv`.

The seven-row `network_evidence_by_hub.csv` and the axis-status
`network_axis_evidence.csv` preserve the requested output shape while marking
unavailable fields `NOT VERIFIED`. No predicted or experimentally supported
axis is claimed. If a future analysis identifies a miRNA-mediated mechanism,
its regulatory level should be described as **post-transcriptional**, not
post-translational; the local network writeup needed to check the draft's
wording was not found.

## Limitations

This is a local-file evidence audit, not a negative biological result. The
necessary interaction outputs, database provenance, filtering criteria,
miRNA expression data, and lncRNA analysis outputs are not present in the
audited files. The network cannot be reconstructed responsibly without those
inputs. Any future network should use the validated seven-gene list, keep
predicted database links separate from experimentally supported interactions,
define detection and support thresholds before filtering, and provide the
sample-level source data and correction family for correlation analyses.

## Figure legends

**Figure F. Regulatory network evidence status.** No local interaction table
was available from which to draw lncRNA–miRNA–mRNA edges; the graphic marks
the network as not verified and does not imply zero biological interactions.

**Figure G. miRNA support status by validated hub.** Per-hub miRNA support
counts are not verified because the cited interaction evidence tables were
not found.

**Figure H. Regulatory evidence tiers.** Tier definitions and assignments
are not verified; no Sankey or evidence-tier flow is inferred.

**Figure I. Single-cell lncRNA–mRNA association status.** Saved lncRNA
detection and correlation outputs were unavailable, so adjusted and
unadjusted distributions cannot be shown.
