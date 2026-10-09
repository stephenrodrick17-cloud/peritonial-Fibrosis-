# Task 1 Audit Report: Verification and Discrepancy Analysis

## 1. Exact Permutation P in GSE62928 Validation
In run_pipeline.py (lines 770-785), the exact permutation test fits an unpenalized logistic regression (penalty=None) with 11 features on only N=8 samples. Because the number of parameters exceeds the sample size (p=11 > n=8), the unconstrained linear model has complete separation on all C(8,4)=70 binary label combinations. Consequently, every single permutation achieves an apparent in-sample AUC of exactly 1.0 (mean AUC across permutations = 1.000). Since the apparent test statistic on the true labels was also 1.0, the condition (a >= auc_apparent) evaluates to True for all 70 permutations, yielding n_extreme_perm = 70 / 70 = 1.0000. This is a severe model-overparameterization artifact: an unpenalized model evaluated in-sample on permuted labels cannot detect a true signal. The test must be relabeled and corrected using cross-validated permutation or non-parametric rank tests on a fixed composite score.

- Combinations evaluated: 70
- Number with apparent AUC >= 1.0: 70
- Resulting reported P-value: 1.0000

## 2. Hub Derivation Verification
The three-model feature table contains 25 rows and marks 7 genes as majority hubs. However, the separate convergent_WGCNA_ECM_genes.csv table marks all seven as in_Top_WGCNA_Module=0 and Convergent=0. Therefore, the available files do not verify that these hubs were selected from a 25-gene WGCNA-convergent pool. The alternate 31-gene Up-ECM analysis is a distinct analysis and yields 8 majority hubs (TNFSF15, EBI3, FLT3LG, ADAM19, CLEC4F, CXCL14, PI3, TGM3).

- Three-model feature rows: 25; strict hubs (['FLT3LG', 'TNFSF15']); majority hubs (['FLT3LG', 'TNFSF15', 'LTB', 'ADAM19', 'SERPINA10', 'EBI3', 'CXCL14'])
- Hubs marked convergent in the WGCNA table: []
- Alternate Up-ECM pool (N=31): 2 strict hubs (['TNFSF15', 'EBI3']) and 8 majority hubs (['TNFSF15', 'EBI3', 'PI3', 'FLT3LG', 'ADAM19', 'CLEC4F', 'TGM3', 'CXCL14'])

## 3. Module Justification: MEblack vs MEmagenta
MEmagenta is the strongest module-level stage association (r = +0.4804, P = 0.00466; 37 genes) and contains zero upregulated ECM DEGs. MEblack contains 2,398 genes and 25 upregulated ECM DEGs, but its eigengene is not associated with stage (r = +0.0285, P = 0.8748). All seven majority hubs are assigned to black, yet their calculated kME_black correlations are weak and nonsignificant while their individual stage gene-significance correlations are positive. Thus, neither module is fully supported as a unified immune-matrisome hub cluster: the current files support a distinct MEmagenta module-trait association and individual hub-stage associations, but do not establish black-module co-expression or the claimed WGCNA-convergent provenance.

| Hub | kME to MEblack | kME P | Stage gene significance (r) | GS P |
|---|---:|---:|---:|---:|
| FLT3LG | -0.0957 | 0.5964 | 0.6423 | 5.576e-05 |
| TNFSF15 | 0.0393 | 0.8283 | 0.6502 | 4.214e-05 |
| LTB | -0.2996 | 0.09032 | 0.4616 | 0.006851 |
| ADAM19 | -0.0248 | 0.8911 | 0.5949 | 0.0002611 |
| SERPINA10 | -0.2215 | 0.2155 | 0.5404 | 0.001167 |
| EBI3 | -0.0773 | 0.6692 | 0.4945 | 0.003441 |
| CXCL14 | -0.0213 | 0.9063 | 0.4672 | 0.006116 |

## 4. Illumina Normalization Protocol
According to the official GEO series matrix header (!Sample_data_processing) for GSE125498, bead-level raw data were summarized using Illumina BeadStudio Software v3, imported into R, and normalized using the quantile normalization method via the Lumi package (Du et al., Bioinformatics 2008). Probes were filtered to retain only those with detectable signal intensity in at least five samples.
