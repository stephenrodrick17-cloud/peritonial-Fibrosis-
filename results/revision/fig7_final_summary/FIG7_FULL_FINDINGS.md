# Figure 7 — Estimation of immune and stromal cell-population scores

**Evidence standard.** This report uses the saved Figure 7 outputs and their recorded source files. Values below are transcribed at the precision saved in those files; missing inference is shown as `—`, not inferred. Run 2 sensitivity analyses were recomputed from the saved expression matrices and frozen signature definition using fixed seeds. The generated results and seeds are in `run_metadata_summary.json`.

## 1. Analysis overview

The analysis asks whether a frozen five-gene panel score is associated with enrichment of immune and stromal marker signatures in peritoneal dialysis effluent. Run 1 explores whether the panel separates the 13 long-term PD samples in GSE125498 into stable groups; Run 2 uses a continuous score in the discovery cohort and two additional effluent cohorts.

The frozen `panel_v2` panel is **FLT3LG, EBI3, LTB, ADAM19, CST7**. It comes from `results/revision/panel_v2/panel_v2.json`, which records that the rule-derived panel was defined on 2026-10-11 before any *future* new test data, but also says that its scores on the existing cohorts are post hoc and exploratory, not validation claims.

All 13 frozen human CellMarker_2024 signatures were scored. ssGSEA used minimum signature size 3; these are **marker-enrichment scores, not cell fractions**.

| Signature | Source genes | GSE125498 detected | GSE130888 detected | GSE248762 detected | GSE57070 detected |
|---|---:|---:|---:|---:|---:|
| T cells | 43 | 43 | 43 | 43 | 43 |
| B cells | 47 | 47 | 47 | 47 | 46 |
| NK cells | 5 | 5 | 5 | 5 | 4 |
| Monocytes | 59 | 59 | 59 | 59 | 56 |
| M1 macrophages | 10 | 10 | 10 | 10 | 9 |
| M2 macrophages | 10 | 10 | 10 | 10 | 10 |
| Neutrophils | 22 | 22 | 22 | 22 | 21 |
| Mast cells | 22 | 22 | 21 | 22 | 21 |
| Immature dendritic cells | 22 | 22 | 22 | 22 | 21 |
| Mature dendritic cells | 4 | 4 | 4 | 4 | 4 |
| Fibroblasts | 27 | 27 | 27 | 27 | 25 |
| Myofibroblasts | 9 | 9 | 9 | 9 | 9 |
| Mesothelial cells | 9 | 9 | 9 | 9 | 9 |

**Source:** `ssgsea_signature_provenance.csv`, and `signature_detection_GSE125498.csv`, `signature_detection_GSE130888.csv`, `signature_detection_GSE248762.csv`, `signature_detection_GSE57070.csv`.

Panel/signature overlap was checked explicitly. **Only T cells share a panel gene: LTB** (one gene). B cells, NK cells, monocytes, M1/M2 macrophages, neutrophils, mast cells, immature/mature dendritic cells, fibroblasts, myofibroblasts, and mesothelial cells have zero panel-gene overlap.

**Source:** `panel_signature_overlap.csv`.

## 2. Datasets and provenance

| Accession / candidate | Tissue or sample type | Platform | Samples and group structure | Role and reason |
|---|---|---|---|---|
| GSE125498 | Peritoneal dialysis effluent cells | Illumina HumanHT-12 v4 array | 33 total: short-term PD=20; long-term PD=13. Run 1 and Run 2 discovery use the 13 LPD/long-term samples. | **Discovery.** Effluent tissue; the requested LPD subset has 13 samples. |
| GSE130888 | Peritoneal effluent cells plus normal peritoneal tissue controls | 10x Genomics scRNA-seq | 13 total: short-term PD=6; long-term PD=4; normal tissue=3. **10 effluent samples used**; normal-tissue controls excluded. | **Validation cohort in the saved analysis**, but exploratory/post hoc on existing data; pseudobulk from effluent samples. |
| GSE248762 | Overnight peritoneal dialysis effluent | 10x Genomics scRNA-seq | 16: SV=6; LV_NOT_UF=6; LV_UF=4. | **Validation cohort in the saved analysis**, but exploratory/post hoc on existing data; patient-level pseudobulk. |
| GSE57070 | Peritoneal dialysis effluent cells | Affymetrix Human Gene 1.0 ST, GPL11209 | Six arrays (AlaGln=3; control=3) from 3 patients, with paired treatment/control arrays averaged within patient; n=3 patients. | **Descriptive only.** Inference withheld at n=3; not validation. |
| GSE297829 | Peritoneal dialysis effluent cells | 10x Genomics scRNA-seq | One patient after approximately 3 years of PD. | **Not used for association.** One patient is not an independent association cohort. |
| GSE62928 | Peritoneal membrane biopsy | Affymetrix array | 8: EPS=4; PD=2; uremic controls=2. | **Not used.** Different tissue from effluent. |
| GSE92453 | Cultured and effluent-derived mesothelial cells | Expression array | 38: ex-vivo control comparisons=4; epithelioid comparisons=18; non-epithelioid comparisons=16. | **Not used.** Mixed cell culture/ex-vivo selection, not bulk effluent. |
| GSE92455 | NOT VERIFIED | Series-matrix URL returned 404 | NOT VERIFIED | **Not used.** No public GEO series could be verified from the saved audit. |
| GSE121372 | Treated human peritoneal mesothelial cells | Expression array | 4: control 6h=1; TGF-beta1 6h=1; control 24h=1; TGF-beta1 24h=1. | **Not used.** Cultured cells, not patient effluent. |
| GSE146955 | NOT VERIFIED | Series-matrix URL returned 404 | NOT VERIFIED | **Not used.** No public GEO series could be verified from the saved audit. |
| E-MTAB-5439 | Effluent-derived sorted CD14+ macrophages and CD1c+ dendritic cells | Affymetrix array | 10 arrays from 5 donors: 5 CD14+ macrophage; 5 CD1c+ dendritic-cell samples. | **Not used.** Sorted populations are not bulk effluent; saved audit found raw CEL files and a final gene list, not a processed gene-expression matrix. |
| E-MTAB-5461 | Peritoneal-dialysis effluent microRNA | miRNA array | 10 arrays: AlaGln=5; control=5; incomplete patient pairing. | **Not used.** miRNA cannot score the protein-coding mRNA panel/signatures. |
| E-MTAB-5462 | Peritoneal-dialysis effluent bulk RNA-seq | Bulk RNA-seq | 11 assay samples from 6 patients in a crossover treatment/control study. | **Not processed.** The saved audit records IDF/SDRF but no usable processed gene-level matrix; raw-read alignment/quantification was not undertaken. |

**Sources:** `dataset_audit.csv`, `sample_metadata_GSE130888.csv`, `sample_metadata_GSE248762.csv`, `sample_metadata_GSE57070.csv`, `sample_metadata_GSE57070_assays.csv`.

### SRA screening

The saved search returned **95 run records across nine BioProjects**. Each project was handled as follows:

| BioProject | Linked GEO / ArrayExpress | Saved decision |
|---|---|---|
| PRJNA541920 | GSE130888 | Analyzed through GEO raw matrices; duplicate SRA/GEO deposit. |
| PRJNA1045728 | GSE248762 | Analyzed through GEO raw matrices; duplicate SRA/GEO deposit. |
| PRJNA1266235 | GSE297829 | Excluded: one patient. |
| PRJEB19316 | E-MTAB-5462 | Not analyzed: no processed gene-level matrix; FASTQ quantification not undertaken. |
| PRJNA638799 | — | Excluded: peritoneal biopsies and small-RNA data, not effluent mRNA. |
| PRJNA672097 | GSE160145 | Excluded: heart tissue; SRA search text was a false positive for this analysis. |
| PRJNA757643 | GSE182736 | Excluded: effluent exosomal miRNA, not the panel/signature mRNA. |
| PRJNA1033066 | — | Excluded: pediatric blood TCR repertoire in sepsis; unrelated tissue/question. |
| PRJNA1212711 | — | Excluded: animal model, not human peritoneal effluent. |

**Sources:** `run_metadata.json` (search scope and run-record count), `sra_project_audit.csv`.

The download manifest records **no files downloaded during preparation of this report**. It records `GSE57070_RAW.tar` as downloaded on 2026-10-16. The listed files with an unknown original local-download date and the single dated download are:

| Source file | Saved acquisition status |
|---|---|
| `GPL10558.annot.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GPL11209.annot.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE121372_series_matrix.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE125498_non-normalized_data.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE125498_RAW.tar` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE125498_series_matrix.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE130888_RAW.tar` | Already local; original download date unknown; copied and checksum-verified 2026-10-16. |
| `GSE130888_series_matrix.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE248762_RAW.tar` | Already local; original download date unknown; copied and checksum-verified 2026-10-16. |
| `GSE248762_series_matrix.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE57070_RAW.tar` | Downloaded 2026-10-16. |
| `GSE57070_series_matrix.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE62928_series_matrix.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |
| `GSE92453_series_matrix.txt.gz` | Already local; original download date unknown; checksum-verified 2026-10-16. |

**Source:** `download_manifest.csv`. The manifest is the evidence for these dates; unknown dates remain unknown.

**Source:** `download_manifest.csv`.

## 3. Run 1 — consensus clustering, GSE125498 LPD

Run 1 used Ward/Euclidean consensus clustering, 1000 resamples, pItem=0.8, and seed 20261016; final partitions were complete-linkage cuts of distance 1-consensus. Selection used the maximum Euclidean silhouette on the gene-wise z-scored panel. The saved metrics are:

| k | CDF area | Delta area | Euclidean silhouette | Consensus-distance silhouette | Cluster sizes |
|---:|---:|---:|---:|---:|---|
| 2 | 0.25391025641025644 | — | 0.3930176034806433 | 0.7894334941569128 | 1:12; 2:1 |
| 3 | 0.4785256410256411 | 0.22461538461538466 | 0.3215720552170308 | 0.5705648489924937 | 1:11; 2:1; 3:1 |
| 4 | 0.6741025641025642 | 0.1955769230769231 | 0.2030847035440373 | 0.4834779492408137 | 1:8; 2:1; 3:3; 4:1 |
| 5 | 0.7943589743589743 | 0.1202564102564101 | 0.18054266878703057 | 0.4500770317003806 | 1:6; 2:1; 3:3; 4:2; 5:1 |

**Source:** `metrics_by_k.csv`.

Selected k=2, with sizes 12 and 1. The singleton is **GSM3575622**. Across the 20 saved seed repeats (20261017–20261036), selected k was 2 in all 20 repeats, with primary-split ARI=1.0 in each repeat. Repeatability does not fix the unusable singleton split.

**Sources:** `cluster_sizes_selected_k.csv`, `cluster_assignments_selected_k.csv`, `seed_stability.csv`, `run_metadata.json`.

Random-gene-set nulls used 1000 sets from the 43-gene pool and 1000 sets from all 11741 expressed genes; the observed values, observed percentile, and saved empirical two-sided p for each selected-k metric are:

| Universe | Metric | Observed | Observed percentile | Empirical two-sided p |
|---|---|---:|---:|---:|
| 43_gene_pool | cdf_area | 0.25391025641025644 | 1.2 | 0.025974025974025976 |
| 43_gene_pool | mean_pairwise_consensus | 0.7468578632505738 | 98.8 | 0.025974025974025976 |
| 43_gene_pool | pac_0_1 | 0.47435897435897434 | 22.8 | 0.45754245754245754 |
| 43_gene_pool | silhouette_euclidean_z | 0.3930176034806433 | 98.4 | 0.03396603396603397 |
| 43_gene_pool | silhouette_consensus_distance | 0.7894334941569128 | 62.3 | 0.7552447552447552 |
| all_expressed_genes | cdf_area | 0.25391025641025644 | 3.5000000000000004 | 0.07192807192807193 |
| all_expressed_genes | mean_pairwise_consensus | 0.7468578632505738 | 96.5 | 0.07192807192807193 |
| all_expressed_genes | pac_0_1 | 0.47435897435897434 | 24.2 | 0.4855144855144855 |
| all_expressed_genes | silhouette_euclidean_z | 0.3930176034806433 | 91.3 | 0.17582417582417584 |
| all_expressed_genes | silhouette_consensus_distance | 0.7894334941569128 | 60.3 | 0.7952047952047953 |

**Source:** `null_observed_vs_random.csv` (null samples and distributions: `null_distributions.csv`).

All 13 ssGSEA Mann–Whitney tests and leukocyte-adjusted results are shown; BH q-values are within the 13 signatures.

| Signature | MW p | MW BH q | Leukocyte-adjusted p | Leukocyte-adjusted BH q |
|---|---:|---:|---:|---:|
| T cells | 0.15384615384615385 | 0.8 | 0.2512851403111102 | 0.8816889718142173 |
| B cells | 0.9230769230769231 | 0.9230769230769231 | 0.46823168299559403 | 0.8816889718142173 |
| NK cells | 0.3076923076923077 | 0.8 | 0.956043292212275 | 0.956043292212275 |
| Monocytes | 0.3076923076923077 | 0.8 | 0.6217324241262229 | 0.8980579459600997 |
| M1 macrophages | 0.9230769230769231 | 0.9230769230769231 | 0.6980808376100366 | 0.9075050888930475 |
| M2 macrophages | 0.15384615384615385 | 0.8 | 0.19000116288817065 | 0.8816889718142173 |
| Neutrophils | 0.6153846153846154 | 0.9230769230769231 | 0.5425778288087492 | 0.8816889718142173 |
| Mast cells | 0.46153846153846156 | 0.8571428571428572 | 0.4941126660250593 | 0.8816889718142173 |
| Immature dendritic cells | 0.46153846153846156 | 0.8571428571428572 | 0.3598083840661513 | 0.8816889718142173 |
| Mature dendritic cells | 0.7692307692307693 | 0.9230769230769231 | 0.4685706557378302 | 0.8816889718142173 |
| Fibroblasts | 0.3076923076923077 | 0.8 | 0.2608310869603969 | 0.8816889718142173 |
| Myofibroblasts | 0.9230769230769231 | 0.9230769230769231 | 0.9338765344446367 | 0.956043292212275 |
| Mesothelial cells | 0.9230769230769231 | 0.9230769230769231 | 0.9338765344446367 | 0.956043292212275 |

**Sources:** `ssgsea_tests_selected_k.csv`, `ssgsea_pairwise_mann_whitney_selected_k.csv`.

**Run 1 is not a usable biological subgroup result.** k=2 isolates one sample, GSM3575622, from the other 12. For a 12-versus-1 split, the smallest possible two-sided Mann–Whitney p in the saved calculation is 0.15384615384615385; thus no signature could pass a conventional 0.05 significance threshold, even before multiple-testing correction. These are exploratory clusters, not subtypes.

## 4. Run 2 — continuous panel score and cohort correlations

Within each cohort, each panel gene was z-scored across that cohort’s analyzed samples and the panel score was the mean of those gene-wise z-scores. Each signature was scored by ssGSEA. Spearman correlations and two-sided p-values were computed; BH correction was within each cohort across 13 signatures. The 95% CIs are percentile bootstrap intervals from 2000 resamples with fixed seed 20261023. GSE57070 has only n=3, so the saved analysis reports rho only and withholds CI, p, and q.

| Cohort | n | Signature | rho | 95% CI | p | BH q |
|---|---:|---|---:|---|---:|---:|
| GSE125498 | 13 | T cells | 0.8571428571428571 | [0.5210537736221249, 0.9774647887323943] | 0.00018093408872524061 | 0.002352143153428128 |
| GSE125498 | 13 | B cells | 0.4835164835164835 | [-0.13813444401394803, 0.841852097339267] | 0.09413469259693935 | 0.40791700125340385 |
| GSE125498 | 13 | NK cells | 0.39010989010989006 | [-0.30173641739357704, 0.8651685393258428] | 0.18757411613383787 | 0.609615877434973 |
| GSE125498 | 13 | Monocytes | 0.521978021978022 | [-0.020456618401550563, 0.8324255121042828] | 0.06729171199351022 | 0.40791700125340385 |
| GSE125498 | 13 | M1 macrophages | -0.22527472527472525 | [-0.6591609143005991, 0.4190751086281808] | 0.45930514272500556 | 0.8997542512357684 |
| GSE125498 | 13 | M2 macrophages | 0.09890109890109891 | [-0.5308228297894936, 0.6667142857142854] | 0.7478683035983906 | 0.8997542512357684 |
| GSE125498 | 13 | Neutrophils | -0.11538461538461539 | [-0.6079696394686908, 0.41344867026656895] | 0.7073898735142603 | 0.8997542512357684 |
| GSE125498 | 13 | Mast cells | -0.06593406593406594 | [-0.8005839646464645, 0.6134453781512605] | 0.830542385756094 | 0.8997542512357684 |
| GSE125498 | 13 | Immature dendritic cells | 0.0934065934065934 | [-0.45918173815620994, 0.7032392444157141] | 0.7615000227827017 | 0.8997542512357684 |
| GSE125498 | 13 | Mature dendritic cells | -0.021978021978021976 | [-0.6779536771728748, 0.5687376074025109] | 0.9431862416608071 | 0.9431862416608071 |
| GSE125498 | 13 | Fibroblasts | 0.22527472527472525 | [-0.5027940851560774, 0.9047899159663864] | 0.45930514272500556 | 0.8997542512357684 |
| GSE125498 | 13 | Myofibroblasts | 0.09890109890109891 | [-0.5549919122578608, 0.8207408102905965] | 0.7478683035983906 | 0.8997542512357684 |
| GSE125498 | 13 | Mesothelial cells | 0.09890109890109891 | [-0.5595435285053929, 0.8034432120106274] | 0.7478683035983906 | 0.8997542512357684 |
| GSE130888 | 10 | T cells | 0.9393939393939393 | [0.6226415094339623, 1.0] | 5.4840529985136686e-05 | 0.0007129268898067769 |
| GSE130888 | 10 | B cells | 0.8424242424242423 | [0.37096966802006215, 1.0] | 0.002220031225916843 | 0.009620135312306318 |
| GSE130888 | 10 | NK cells | 0.709090909090909 | [0.14272321428571427, 0.9600196078431372] | 0.02166592336753244 | 0.07041425094448042 |
| GSE130888 | 10 | Monocytes | 0.5030303030303029 | [-0.20790990259740258, 0.8749999999999999] | 0.13833369839449194 | 0.25690543987548503 |
| GSE130888 | 10 | M1 macrophages | 0.5393939393939393 | [-0.20253164556962028, 0.9607906799083269] | 0.10759318778241479 | 0.23311857352856538 |
| GSE130888 | 10 | M2 macrophages | -0.13939393939393938 | [-0.8726515242558989, 0.5926511954992965] | 0.7009318849100584 | 0.9073638178128159 |
| GSE130888 | 10 | Neutrophils | 0.3696969696969697 | [-0.43756249999999997, 0.8462264150943393] | 0.29305007529282234 | 0.4762063723508363 |
| GSE130888 | 10 | Mast cells | 0.04242424242424241 | [-0.708996383363472, 0.723292127521145] | 0.9073638178128159 | 0.9073638178128159 |
| GSE130888 | 10 | Immature dendritic cells | -0.24848484848484845 | [-0.8272331154684096, 0.5586109542631272] | 0.4887763045192432 | 0.7060102176389068 |
| GSE130888 | 10 | Mature dendritic cells | 0.8787878787878788 | [0.4437086092715231, 1.0] | 0.0008138621117322101 | 0.005290103726259365 |
| GSE130888 | 10 | Fibroblasts | -0.6727272727272726 | [-0.9876543209876544, 0.03276740237690724] | 0.03304122254543772 | 0.08590717861813806 |
| GSE130888 | 10 | Myofibroblasts | -0.06666666666666665 | [-0.719751592356688, 0.65] | 0.8548130882487427 | 0.9073638178128159 |
| GSE130888 | 10 | Mesothelial cells | -0.06666666666666665 | [-0.7013588263588263, 0.6352771226415092] | 0.8548130882487427 | 0.9073638178128159 |
| GSE248762 | 16 | T cells | 0.9558823529411766 | [0.8179036432188833, 0.9969653282032328] | 7.760274461610307e-09 | 1.0088356800093399e-07 |
| GSE248762 | 16 | B cells | 0.7147058823529412 | [0.3153213751868461, 0.894910582224015] | 0.0018627410999474342 | 0.006053908574829161 |
| GSE248762 | 16 | NK cells | 0.8911764705882352 | [0.6461319112480235, 0.9821428571428572] | 3.618648907351986e-06 | 2.3521217897787908e-05 |
| GSE248762 | 16 | Monocytes | 0.3411764705882353 | [-0.24439763773388196, 0.844561445654495] | 0.19593262486788887 | 0.3183905154103194 |
| GSE248762 | 16 | M1 macrophages | 0.19411764705882353 | [-0.3889805534543621, 0.7391596365996105] | 0.47128972836210753 | 0.5569787698824907 |
| GSE248762 | 16 | M2 macrophages | -0.3411764705882353 | [-0.7347351802325975, 0.2214465204247799] | 0.19593262486788887 | 0.3183905154103194 |
| GSE248762 | 16 | Neutrophils | -0.2029411764705882 | [-0.6707367277639842, 0.3760937136361186] | 0.4509612417227097 | 0.5569787698824907 |
| GSE248762 | 16 | Mast cells | -0.3058823529411765 | [-0.712574850299401, 0.23100989709322628] | 0.24925429892891854 | 0.3600339873417712 |
| GSE248762 | 16 | Immature dendritic cells | -0.4588235294117647 | [-0.8045500113921165, 0.06983835929490506] | 0.07383553449548094 | 0.15997699140687538 |
| GSE248762 | 16 | Mature dendritic cells | 0.6941176470588235 | [0.23418194980694984, 0.9370361904249668] | 0.0028537059034569076 | 0.007419635348987959 |
| GSE248762 | 16 | Fibroblasts | -0.7470588235294118 | [-0.9245880581148123, -0.34916153179058895] | 0.0008821069929846312 | 0.003822463636266735 |
| GSE248762 | 16 | Myofibroblasts | -0.07352941176470588 | [-0.5832083958020989, 0.48889373177040324] | 0.7866792491258566 | 0.7866792491258566 |
| GSE248762 | 16 | Mesothelial cells | -0.07352941176470588 | [-0.6112780712609723, 0.5329390916961257] | 0.7866792491258566 | 0.7866792491258566 |
| GSE57070 | 3 | T cells | 1.0 | — | — | — |
| GSE57070 | 3 | B cells | 1.0 | — | — | — |
| GSE57070 | 3 | NK cells | 0.5 | — | — | — |
| GSE57070 | 3 | Monocytes | 1.0 | — | — | — |
| GSE57070 | 3 | M1 macrophages | 1.0 | — | — | — |
| GSE57070 | 3 | M2 macrophages | 1.0 | — | — | — |
| GSE57070 | 3 | Neutrophils | 0.5 | — | — | — |
| GSE57070 | 3 | Mast cells | 1.0 | — | — | — |
| GSE57070 | 3 | Immature dendritic cells | 0.5 | — | — | — |
| GSE57070 | 3 | Mature dendritic cells | 0.5 | — | — | — |
| GSE57070 | 3 | Fibroblasts | -1.0 | — | — | — |
| GSE57070 | 3 | Myofibroblasts | -0.5 | — | — | — |
| GSE57070 | 3 | Mesothelial cells | -0.5 | — | — | — |

**Source:** `correlations_all_datasets.csv`; the same machine-readable rows are copied without changes to `all_results_flat.csv`. The known values match the saved file: GSE125498 T cells rho=0.8571428571428571, q=0.002352143153428128; GSE130888 T cells rho=0.9393939393939393; GSE248762 T cells rho=0.9558823529411766.

### Findings surviving BH and direction across cohorts

- **T cells** are positively associated in all three inferential cohorts and survive BH in each: GSE125498 q=0.002352143153428128; GSE130888 q=0.0007129268898067769; GSE248762 q=1.0088356800093399e-07.
- In the two validation-design cohorts, **B cells** and **mature dendritic cells** are positive and survive BH in both; **NK cells** are positive and survive BH in GSE248762 but not GSE130888 (q=0.07041425094448042 there). These associations do not survive BH in discovery.
- **Fibroblasts are directionally inconsistent**: discovery rho=0.22527472527472525, GSE130888 rho=-0.6727272727272726, and GSE248762 rho=-0.7470588235294118. The GSE248762 result survives BH; GSE130888 does not.
- No other signature survives BH in discovery. In GSE130888, the remaining BH-significant results are absent; in GSE248762, the significant associations are T cells, B cells, NK cells, mature dendritic cells, and fibroblasts.
- GSE57070 rho values are necessarily coarse at n=3 (the saved values are exactly 1.0, 0.5, -1.0, or -0.5 here) and **must not be interpreted**.

## 5. Clinical covariates

Availability is based on the saved cohort-level sample metadata audit. “Vintage group” means a study-defined group; exact vintage in months/duration is not available.

| Cohort | Age | Sex | Vintage group | Exact vintage/duration | UF status | Transport | Peritonitis | Batch |
|---|---|---|---|---|---|---|---|---|
| GSE125498 | No | No | Yes | No | No | No | No | No |
| GSE130888 | No | No | Yes | No | No | No | No | No |
| GSE248762 | Yes | Yes | Yes | No | Yes | No | No | No |
| GSE57070 | No | No | No | No | No | No | Yes (previous peritonitis) | No |

**Source:** `clinical_covariate_availability.csv`.

All saved clinical tests were adjusted for the leukocyte proxy; effects are unstandardized OLS coefficients with 95% CIs. BH q-values are within dataset.

| Cohort | Test | n | Adjusted effect | 95% CI | p | BH q |
|---|---|---:|---:|---|---:|---:|
| GSE130888 | long_vintage_vs_short | 10 | 0.6653092384380773 | [-0.05578544307146338, 1.3864039199476181] | 0.06547696553371139 | 0.06547696553371139 |
| GSE248762 | age_years | 16 | 0.001110292454477163 | [-0.02916639476188033, 0.03138697967083466] | 0.9380609748469997 | 0.9380609748469997 |
| GSE248762 | sex_M_vs_F | 16 | -0.06309398489487406 | [-1.1093711526425212, 0.983183182852773] | 0.8983412793532496 | 0.9380609748469997 |
| GSE248762 | long_vintage_vs_short | 16 | 0.6438044835301039 | [-0.2488856385358139, 1.5364946055960216] | 0.14322486211438779 | 0.5728994484575511 |
| GSE248762 | UF_failure_vs_no_UF_in_long_vintage | 10 | 0.5144461802213838 | [-1.2059281562402744, 2.234820516683042] | 0.5023596025089293 | 0.9380609748469997 |

**Source:** `clinical_associations.csv`. **No adjusted model survived BH correction**; the smallest saved q is 0.06547696553371139.

## 6. New sensitivity checks

All sensitivity correlations below use the three inferential cohorts and the same cohort-wise gene z-scoring and Spearman test convention as Run 2. CIs use 2000 paired-sample bootstrap resamples. Fixed seeds and software versions are recorded in `run_metadata_summary.json`. The T-cell signature was rescored with LTB excluded using the frozen signature genes and saved expression matrices. GSEApy's current NES absolute scale differed from the archived score scale, but sample ranks matched exactly (rank correlations 1.0, 0.9999999999999999, and 1.0), so Spearman associations are rank-equivalent. q-values in the table below are BH within each cohort across the three LTB-exclusion contrasts.

| Cohort | LTB exclusion contrast | rho | 95% CI | p | BH q |
|---|---|---:|---|---:|---:|
| GSE125498 | Remove LTB from T-cell signature; retain 5-gene panel | 0.8571428571428571 | [0.5629945236521575, 0.9764740282077743] | 0.00018093408872524061 | 0.0005428022661757219 |
| GSE125498 | Remove LTB from panel; retain full T-cell signature | 0.7582417582417583 | [0.32579522076277895, 0.9491525423728815] | 0.0026660786121372876 | 0.0026660786121372876 |
| GSE125498 | Remove LTB from both panel and T-cell signature | 0.7582417582417583 | [0.31690960451977396, 0.9422048368879143] | 0.0026660786121372876 | 0.0026660786121372876 |
| GSE130888 | Remove LTB from T-cell signature; retain 5-gene panel | 0.9393939393939393 | [0.6226415094339623, 1.0] | 5.4840529985136686e-05 | 5.4840529985136686e-05 |
| GSE130888 | Remove LTB from panel; retain full T-cell signature | 0.9393939393939393 | [0.6226415094339623, 1.0] | 5.4840529985136686e-05 | 5.4840529985136686e-05 |
| GSE130888 | Remove LTB from both panel and T-cell signature | 0.9393939393939393 | [0.6153846153846154, 1.0] | 5.4840529985136686e-05 | 5.4840529985136686e-05 |
| GSE248762 | Remove LTB from T-cell signature; retain 5-gene panel | 0.9558823529411766 | [0.8390220686973948, 0.9969606016299911] | 7.760274461610307e-09 | 2.3280823384830922e-08 |
| GSE248762 | Remove LTB from panel; retain full T-cell signature | 0.9323529411764706 | [0.7902693410133838, 0.9822222222222222] | 1.451754439048121e-07 | 1.451754439048121e-07 |
| GSE248762 | Remove LTB from both panel and T-cell signature | 0.9323529411764706 | [0.7829947839736157, 0.9880068934142793] | 1.451754439048121e-07 | 1.451754439048121e-07 |

**Source:** `sensitivity_ltb_excluded.csv`. Panel/signature overlap contained no other shared genes/signatures, so no additional overlap-specific sensitivity was applicable (`panel_signature_overlap.csv`). The T-cell association remains positive and BH-significant across all three cohorts when the shared gene is removed from the signature, when it is removed from the panel, and when it is removed from both. In discovery, removing LTB from the panel lowers rho from 0.8571428571428571 to 0.7582417582417582; therefore the overlap contributes but does not alone account for the relationship.

Leave-one-gene-out correlations use each four-gene panel and the original full T-cell signature; only rho was requested for this comparison. The CSV also preserves the paired bootstrap intervals and nominal p-values.

| Omitted panel gene | GSE125498 rho | GSE130888 rho | GSE248762 rho |
|---|---:|---:|---:|
| FLT3LG | 0.7857142857142857 | 0.9515151515151514 | 0.9264705882352942 |
| EBI3 | 0.9175824175824175 | 0.9757575757575757 | 0.9852941176470588 |
| LTB | 0.7582417582417583 | 0.9393939393939393 | 0.9323529411764706 |
| ADAM19 | 0.8296703296703297 | 0.9393939393939393 | 0.9676470588235294 |
| CST7 | 0.8186813186813187 | 0.9393939393939393 | 0.9294117647058824 |

**Sources:** `sensitivity_leave_one_gene_out.csv`, `sensitivity_leave_one_gene_out.png`.

The leukocyte proxy is the mean of the nine non-T-cell immune ssGSEA marker-enrichment scores, avoiding direct T-cell/LTB overlap. All reported cohort correlations are:

| Cohort | n | rho | 95% CI | p | BH q across cohorts |
|---|---:|---:|---|---:|---:|
| GSE125498 | 13 | 0.32967032967032966 | [-0.27172152194211013, 0.762106315614064] | 0.27133505450523593 | 0.27133505450523593 |
| GSE130888 | 10 | 0.6121212121212121 | [-0.044886431919240466, 0.9867614173905563] | 0.05997214247725527 | 0.1799164274317658 |
| GSE248762 | 16 | 0.31764705882352945 | [-0.22623616137960398, 0.8576480200397738] | 0.2305734636218452 | 0.27133505450523593 |

**Source:** `sensitivity_leukocyte_proxy.csv`; saved original proxy results: `panel_vs_leukocyte_proxy.csv`. None survives BH across these cohorts. The estimates do not establish that the panel is independent of overall immune content.

## 7. Figures

| Figure file | What it shows |
|---|---|
| `fig7_A_correlation_heatmap.png` | Heatmap of Run 2 panel-score/Spearman-rho associations across the 13 signatures and the three inferential cohorts. |
| `fig7_B_top_associations.png` | Discovery scatterplots and smoothed trends for the three strongest absolute discovery associations shown in the saved figure: T cells, monocytes, and B cells. |
| `fig7_C_validation_replication.png` | Cohort-specific correlation estimates and bootstrap intervals across discovery and the two effluent pseudobulk cohorts; title explicitly describes cross-cohort concordance, not same-method replication. |
| `Fig7_exploratory.png` | Run 1 consensus CDF, selection metrics, k=2 consensus matrix, and PCA of the five z-scored panel genes. |
| `Fig7_gene_panel_descriptive.png` | Run 1 panel gene expression by the selected clusters; descriptive because these genes define the clusters. |
| `Fig7_immune_scores.png` | Run 1 ssGSEA marker-enrichment scores by selected cluster, with Mann–Whitney q labels. |
| `sensitivity_leave_one_gene_out.png` | New grouped bar chart of T-cell rho after omitting each panel gene, by cohort. |

GSE57070 values are based on three independent patients; its possible rho values are constrained to exactly +/-1 or +/-0.5. These points must not be interpreted as evidence or validation.

## 8. Interpretation and limitations

The continuous five-gene panel score has a positive association with the T-cell marker-enrichment score in the discovery cohort and both effluent validation-design cohorts. The T-cell direction persists when LTB is removed from the signature and when it is removed from the panel, although the discovery rho decreases when the panel excludes LTB. B-cell and mature dendritic-cell associations survive BH correction in both validation-design cohorts, while NK-cell evidence survives correction only in GSE248762. Fibroblast association is directionally inconsistent: it is slightly positive in discovery and negative in both effluent pseudobulk cohorts, with BH significance only in GSE248762. Panel-to-leukocyte-proxy correlations were positive in all three cohorts but no adjusted association survived BH; this does not rule out immune-content confounding. Discovery is bulk microarray, whereas both validation-design cohorts are single-cell-derived patient pseudobulk, so they are not same-method replications. The cohort sizes are small (13, 10, and 16), intervals are often wide, and GSE57070 has only three patients. Exact dialysis vintage and batch were unavailable, and patient-level age, sex, transport, and peritonitis information was incomplete or absent across cohorts. This work is exploratory and does not constitute clinical validation.

## 9. Paper-ready text

### Figure legend

**Figure 7. Exploratory association of a frozen five-gene panel with immune/stromal marker enrichment in peritoneal dialysis effluent.** The panel score is the mean of cohort-wise gene z-scores for FLT3LG, EBI3, LTB, ADAM19, and CST7. Spearman correlations compare that score with 13 CellMarker_2024 ssGSEA marker-enrichment scores in GSE125498 (LPD discovery), GSE130888, and GSE248762; intervals are 95% percentile-bootstrap CIs. GSE57070 is descriptive only (three patients) and is excluded from inference. Removing LTB from the T-cell signature and panel tests sensitivity to the shared marker. Enrichment scores are not cell fractions. Cohort sizes, platforms, and preparation differ, and all results are exploratory, not clinical validation.

### Results paragraph

The frozen panel score was positively associated with T-cell ssGSEA enrichment in GSE125498 (n=13; rho=0.8571428571428571; 95% CI 0.5210537736221249 to 0.9774647887323943; BH q=0.002352143153428128), GSE130888 (n=10; rho=0.9393939393939393; 95% CI 0.6226415094339623 to 1.0; q=0.0007129268898067769), and GSE248762 (n=16; rho=0.9558823529411766; 95% CI 0.8179036432188833 to 0.9969653282032328; q=1.0088356800093399e-07). With LTB removed from both panel and signature, T-cell rho remained positive in each cohort (0.7582417582417582, 0.9393939393939393, and 0.9323529411764706, respectively). B-cell and mature dendritic-cell associations passed BH correction in both validation-design cohorts; NK cells passed only in GSE248762. Fibroblast direction differed between discovery and pseudobulk cohorts. Small sample sizes, bulk-versus-pseudobulk differences, incomplete clinical covariates, and the exploratory panel definition preclude clinical validation.

### Plain-language takeaway

The panel tracks T-cell marker enrichment across the three larger effluent cohorts even without LTB, but the small and methodologically different datasets do not establish a clinically validated cell-population or fibrosis measure.

## 10. What is still not done

- GSE57070 CEL reprocessing was not done; the analysis used the GEO processed log2 series matrix, and its three-patient correlations remain descriptive only.
- E-MTAB-5462 raw-read alignment and quantification were not done; the saved audit found no processed gene-level matrix.
- Patient-level exact dialysis vintage and batch were not available in the analyzed GEO sample metadata; obtaining them requires non-public or separately curated clinical records.
- Independent prospective or genuinely new external validation remains outstanding; the panel’s results on the existing datasets are explicitly post hoc/exploratory.
- Clinical validation and adequately powered adjusted models remain outstanding; no adjusted clinical association passed BH correction.

## Files read

The following are the files whose contents were read to prepare this report, calculate the requested checks, or verify the outputs. The separately listed `data/raw/fig7_real_data` archives were inventory-listed and referenced by the manifest but their large raw contents were not re-read for this report.

**Run 1 (`results/revision/fig7_lpd_current/`):**

- `run_fig7_lpd.py`
- `fig7_results_and_legend.md`
- `run_metadata.json`
- `metrics_by_k.csv`
- `ssgsea_tests_selected_k.csv`
- `null_observed_vs_random.csv`
- `null_distributions.csv`
- `ssgsea_signature_provenance.csv`
- `seed_stability.csv`
- `cluster_assignments_selected_k.csv`
- `cluster_sizes_selected_k.csv`
- `clinical_covariate_availability.csv`
- `ssgsea_pairwise_mann_whitney_selected_k.csv`
- `ssgsea_scores_selected_k.csv`
- `Fig7_exploratory.png`
- `Fig7_gene_panel_descriptive.png`
- `Fig7_immune_scores.png`

**Run 2 (`results/revision/fig7_real_data/`):**

- `run_fig7_real_data.py`
- `fig7_results_and_legend.md`
- `run_metadata.json`
- `correlations_all_datasets.csv`
- `panel_signature_overlap.csv`
- `panel_gene_detection.csv`
- `signature_detection_GSE125498.csv`
- `signature_detection_GSE130888.csv`
- `signature_detection_GSE248762.csv`
- `signature_detection_GSE57070.csv`
- `dataset_audit.csv`
- `download_manifest.csv`
- `sra_project_audit.csv`
- `ltb_removed_t_cell_sensitivity.csv`
- `clinical_covariate_availability.csv`
- `clinical_associations.csv`
- `panel_vs_leukocyte_proxy.csv`
- `sample_scores_GSE125498.csv`
- `sample_scores_GSE130888.csv`
- `sample_scores_GSE248762.csv`
- `sample_scores_GSE57070.csv`
- `ssgsea_scores_GSE125498.csv`
- `ssgsea_scores_GSE130888.csv`
- `ssgsea_scores_GSE248762.csv`
- `ssgsea_scores_GSE125498_T_cells_no_LTB.csv`
- `ssgsea_scores_GSE130888_T_cells_no_LTB.csv`
- `ssgsea_scores_GSE248762_T_cells_no_LTB.csv`
- `array_preprocessing_qc.csv`
- `array_preprocessing_qc_GSE57070.csv`
- `sample_metadata_GSE130888.csv`
- `sample_metadata_GSE248762.csv`
- `sample_metadata_GSE57070.csv`
- `sample_metadata_GSE57070_assays.csv`
- `pseudobulk_expression_GSE130888.csv`
- `pseudobulk_expression_GSE248762.csv`
- `gene_expression_GSE125498_LPD.csv`
- `null_distributions_not_run.txt`
- `fig7_A_correlation_heatmap.png`
- `fig7_B_top_associations.png`
- `fig7_C_validation_replication.png`

**Frozen references and newly written outputs:**

- `results/revision/panel_v2/panel_v2.json`
- `results/revision/microenvironment_deconvolution/final_signature_table_both_variants.csv`
- `results/revision/fig7_final_summary/all_results_flat.csv`
- `results/revision/fig7_final_summary/FIG7_FULL_FINDINGS.md` (read for output validation)
- `results/revision/fig7_final_summary/sensitivity_ltb_excluded.csv`
- `results/revision/fig7_final_summary/sensitivity_leave_one_gene_out.csv`
- `results/revision/fig7_final_summary/sensitivity_leukocyte_proxy.csv`
- `results/revision/fig7_final_summary/sensitivity_leave_one_gene_out.png`
- `results/revision/fig7_final_summary/run_metadata_summary.json`
