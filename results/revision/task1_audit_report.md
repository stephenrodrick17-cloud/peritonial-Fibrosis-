# Task 1 Audit Report: Verification and Discrepancy Analysis

## 1. Exact Permutation P in GSE62928 Validation
In run_pipeline.py (lines 770-785), the exact permutation test fits an unpenalized logistic regression (penalty=None) with 11 features on only N=8 samples. Because the number of parameters exceeds the sample size (p=11 > n=8), the unconstrained linear model has complete separation on all C(8,4)=70 binary label combinations. Consequently, every single permutation achieves an apparent in-sample AUC of exactly 1.0 (mean AUC across permutations = 1.000). Since the apparent test statistic on the true labels was also 1.0, the condition (a >= auc_apparent) evaluates to True for all 70 permutations, yielding n_extreme_perm = 70 / 70 = 1.0000. This is a severe model-overparameterization artifact: an unpenalized model evaluated in-sample on permuted labels cannot detect a true signal. The test must be relabeled and corrected using cross-validated permutation or non-parametric rank tests on a fixed composite score.

- Combinations evaluated: 70
- Number with apparent AUC >= 1.0: 70
- Resulting reported P-value: 1.0000

## 2. Hub Derivation Verification
The primary 7 consensus hubs (FLT3LG, TNFSF15, LTB, ADAM19, SERPINA10, EBI3, CXCL14) strictly derive from the 25-gene WGCNA-convergent candidate pool. The alternate 31-gene Up-ECM analysis yields 8 majority hubs (TNFSF15, EBI3, FLT3LG, ADAM19, CLEC4F, CXCL14, PI3, TGM3). The differences stem from PI3 and TGM3 being included in the 31 Up-ECM set but excluded from the 25 WGCNA-convergent pool due to lack of positive module correlation.

- Primary WGCNA pool (N=25): 2 strict hubs (['FLT3LG', 'TNFSF15']) and 7 majority hubs (['FLT3LG', 'TNFSF15', 'LTB', 'ADAM19', 'SERPINA10', 'EBI3', 'CXCL14'])
- Alternate Up-ECM pool (N=31): 2 strict hubs (['TNFSF15', 'EBI3']) and 8 majority hubs (['TNFSF15', 'EBI3', 'PI3', 'FLT3LG', 'ADAM19', 'CLEC4F', 'TGM3', 'CXCL14'])

## 3. Module Justification: MEblack vs MEmagenta
MEmagenta has the highest overall eigengene correlation with clinical dialysis vintage (r = +0.480, P = 0.0047), representing 37 non-matrisomal genes. However, MEmagenta contains 0 ECM genes. In contrast, MEblack contains 2,398 genes and harbors 25 of the 31 pro-fibrotic ECM DEGs (including all 7 consensus hubs). Therefore, MEblack represents the true biological 'hub cluster' for ECM remodeling and immune-matrisome crosstalk, while MEmagenta represents an independent non-ECM stage-associated signature.

## 4. Illumina Normalization Protocol
According to the official GEO series matrix header (!Sample_data_processing) for GSE125498, bead-level raw data were summarized using Illumina BeadStudio Software v3, imported into R, and normalized using the quantile normalization method via the Lumi package (Du et al., Bioinformatics 2008). Probes were filtered to retain only those with detectable signal intensity in at least five samples.
