args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("ggplot2", "dplyr", "tidyr"), "T11")

bundles <- load_all_cohorts()
cohort_table <- do.call(rbind, lapply(bundles, function(b) data.frame(
  GSE = b$cohort,
  Platform_or_assay = if (b$cohort == "GSE125498") "Illumina microarray, RMA log2 expression" else "Effluent single-cell RNA-seq, patient pseudobulk",
  Control_group = if (b$cohort == "GSE125498") "Short-term PD (SPD)" else "Short-vintage/short-term PD (SV)",
  Case_group = if (b$cohort == "GSE125498") "Long-term PD (LPD)" else "Long-vintage/long-term PD (LV)",
  n_control = sum(b$meta$Group == "Control"),
  n_case = sum(b$meta$Group == "Case"),
  n_total = nrow(b$meta),
  Analysis_scale = b$scale
)))
write.csv(cohort_table, file.path(TABLES, "T11_Table1_cohort_characteristics.csv"), row.names = FALSE)

ml <- read_csv_safe(file.path(RESULTS, "revision", "ml_39_matrisome_wgcna", "ml_39_gene_results.csv"))
panel <- ml[ml$Consensus_At_Least_2_of_3 == 1, , drop = FALSE]
cats <- read_csv_safe(file.path(RESULTS, "revision", "ml_39_matrisome_wgcna", "wgcna_matrisome_39_gene_input.csv"))
panel <- merge(panel, cats[, c("Gene", "Matrisome_Category")], by = "Gene", all.x = TRUE)
effects_path <- file.path(TABLES, "T04_hub_effects_by_cohort.csv")
if (file.exists(effects_path)) {
  effects <- read_csv_safe(effects_path)
  for (cohort in c("GSE125498", "GSE248762", "GSE130888")) {
    e <- effects[effects$Cohort == cohort, c("Gene", "log2FC", "FDR_BH_across_hub_tests", "Detected"), drop = FALSE]
    names(e)[-1] <- paste0(cohort, c("_log2FC", "_FDR_hub_BH", "_Detected"))
    panel <- merge(panel, e, by = "Gene", all.x = TRUE)
  }
  concordance <- effects[effects$Cohort != "GSE125498", c("Gene", "Cohort", "log2FC")]
  disc_dir <- setNames(sign(panel$Discovery_log2FC), panel$Gene)
  concordance$Direction_concordant <- sign(concordance$log2FC) == disc_dir[concordance$Gene]
  for (cohort in c("GSE248762", "GSE130888")) {
    e <- concordance[concordance$Cohort == cohort, c("Gene", "Direction_concordant"), drop = FALSE]
    names(e)[2] <- paste0(cohort, "_direction_concordant")
    panel <- merge(panel, e, by = "Gene", all.x = TRUE)
  }
}
for (cohort in c("GSE248762", "GSE130888")) {
  f <- file.path(RESULTS, "revision", "ml_39_matrisome_wgcna", "validation_auc", "all_detectable_consensus_hub_gene_ROC_auc.csv")
  if (file.exists(f)) {
    a <- read_csv_safe(f)
    a <- a[a$cohort == cohort, c("gene", "auc", "ci_low", "ci_high"), drop = FALSE]
    names(a) <- c("Gene", paste0(cohort, c("_AUC", "_AUC_CI_low", "_AUC_CI_high")))
    panel <- merge(panel, a, by = "Gene", all.x = TRUE)
  }
}
disc <- load_discovery()
discovery_auc <- do.call(rbind, lapply(HUBS, function(gene) {
  if (!gene %in% rownames(disc$expr)) return(data.frame(Gene = gene, GSE125498_apparent_AUC = NA_real_,
    GSE125498_AUC_CI_low = NA_real_, GSE125498_AUC_CI_high = NA_real_))
  y <- as.integer(disc$meta$stage_binary)
  score <- as.numeric(disc$expr[gene, ])
  ci <- bootstrap_auc(y, score)
  data.frame(Gene = gene, GSE125498_apparent_AUC = auc_rank(y, score),
    GSE125498_AUC_CI_low = ci[1], GSE125498_AUC_CI_high = ci[2])
}))
panel <- merge(panel, discovery_auc, by = "Gene", all.x = TRUE)
meta_path <- file.path(TABLES, "T04_random_effects_meta_analysis.csv")
if (file.exists(meta_path)) panel <- merge(panel, read_csv_safe(meta_path), by = "Gene", all.x = TRUE)
panel$Feature_selection_methods <- vapply(seq_len(nrow(panel)), function(i) {
  flags <- c(LASSO = panel$LASSO_selected[i], SVM_RFE = panel$SVM_RFE_selected[i], RF = panel$RF_selected[i])
  paste(names(flags)[flags == 1], collapse = ";")
}, character(1))
keep <- unique(c("Gene", "Matrisome_Category", "Matrisome_DEG_Status", "Discovery_log2FC", "Discovery_nominal_P",
  "Discovery_BH_adjusted_P", "Feature_selection_methods",
  grep("^(GSE125498|GSE248762|GSE130888)_(log2FC|FDR|Detected|AUC|AUC_CI|direction_concordant)", names(panel), value = TRUE),
  "log2FC_random", "CI_low", "CI_high", "I2_percent", "tau2_DL"))
panel_table <- panel[, intersect(keep, names(panel)), drop = FALSE]
write.csv(panel_table, file.path(TABLES, "T11_Table2_consensus_hub_summary.csv"), row.names = FALSE)

workflow <- data.frame(step = factor(c("Discovery cohort\nGSE125498", "Matrisome DEGs\n46 genes", "WGCNA intersection\n39 genes", "Three feature selectors\nLASSO / SVM-RFE / RF", "Consensus panel\n9 candidates", "Additional cohorts\nGSE248762 + GSE130888"),
                                    levels = c("Discovery cohort\nGSE125498", "Matrisome DEGs\n46 genes", "WGCNA intersection\n39 genes", "Three feature selectors\nLASSO / SVM-RFE / RF", "Consensus panel\n9 candidates", "Additional cohorts\nGSE248762 + GSE130888")),
                      x = 1:6, y = 1)
pw <- ggplot2::ggplot(workflow, ggplot2::aes(x, y)) + ggplot2::geom_label(ggplot2::aes(label = step),
  fill = "#EAF2F8", colour = "#17324D", size = 3.3, label.size = 0.4, label.padding = grid::unit(0.23, "lines")) +
  ggplot2::geom_segment(data = data.frame(x = 1:5, xend = 2:6, y = 1, yend = 1),
    ggplot2::aes(x = x + 0.2, xend = xend - 0.2, y = y, yend = yend),
    inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.16, "inches")), colour = "#0072B2") +
  ggplot2::scale_x_continuous(limits = c(0.5, 6.5), breaks = NULL) + ggplot2::scale_y_continuous(limits = c(0.7, 1.3), breaks = NULL) +
  ggplot2::labs(title = "Analysis workflow") + theme_project() +
  ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank())
save_plot(pw, "T11_workflow_schematic", 14, 4)

task_log("T11_summary_tables", c("Generated cohort characteristics and consensus-hub summary by joining saved selection, effect-size, meta-analysis, and ROC outputs.",
  "Discovery AUCs are apparent single-gene score AUCs; external cohorts are additional/exploratory, not formally independent.",
  "Table 2 records missing or unavailable comparisons as NA rather than imputing values."))
append_results_text("T11. Summary", paste("The workflow schematic, cohort characteristics, and integrated consensus-hub summary table were generated from saved analysis files.",
  "The final panel contained nine candidates selected by at least two of three methods; external-cohort evidence remains exploratory.",
  "Unavailable measures were left unfilled rather than inferred."))

text <- readLines(file.path(RESULTS, "results_text.md"), warn = FALSE)
legacy_heading <- match("### T08_immune_infiltration", text)
if (!is.na(legacy_heading)) {
  next_heading <- which(seq_along(text) > legacy_heading & grepl("^### ", text))
  end <- if (length(next_heading)) next_heading[[1]] - 1L else length(text)
  text <- text[-seq.int(legacy_heading, end)]
}
limitation_start <- match("### Limitations", text)
if (!is.na(limitation_start)) text <- if (limitation_start > 1L) text[seq_len(limitation_start - 1L)] else character()
text <- c(text, "", "### Limitations", "",
  "This analysis had small sample sizes in the additional effluent cohorts, and the discovery matrisome DEG screen used nominal P-values before filtering.",
  "Machine-learning feature selection was performed in GSE125498; discovery-cohort ROC estimates are apparent and may be optimistic.",
  "GSE248762 and GSE130888 are additional exploratory cohorts, but patient-level independence has not been established.",
  "No wet-lab validation was performed; computational associations and expression-based discrimination do not establish causality or therapeutic efficacy.")
writeLines(text, file.path(RESULTS, "results_text.md"))

files <- list.files(c(FIGURES, TABLES), pattern = "^(T0[1-9]|T10|T11)_", recursive = FALSE, full.names = FALSE)
scripts <- c("01_batch_qc.R", "02_ml_selection.R", "03_enrichment_matrisome.R", "04_hub_expression_meta.R",
             "05_interaction_regulation.R", "06_diagnostic_models.R", "07_single_gene_gsea.R",
             "08_immune_infiltration.R", "09_drug_prediction.R", "10_external_resources.R", "11_summary_tables.R")
readme <- c("# Reproducible peritoneal membrane matrisome analysis", "",
  "Run `Rscript scripts/run_all.R` from the repository root. The scripts use R >= 4.3; the verified environment was R 4.4.2.",
  "Random seed: 2026. All results are written under `results/figures`, `results/tables`, and `results/logs`.",
  "", "## Figure and table inventory", "",
  vapply(files, function(f) {
    task <- sub("^((T0[1-9]|T10|T11)).*", "\\1", f)
    creator <- scripts[match(task, paste0("T", sprintf("%02d", 1:11)))]
    if (!length(creator) || is.na(creator)) creator <- "run_all.R"
    paste0("- `", f, "` — created by `scripts/", creator, "`.")
  }, character(1)),
  "", "## Conclusion", "",
  "Three feature-selection methods applied to the 39-gene WGCNA–MatrisomeAnalyzeR intersection yielded nine consensus candidates; the additional-cohort evidence is exploratory, and CCBE1 showed an inverse expression-based AUC in GSE248762.",
  "", "## Conservative omissions", "",
  "See `results/logs/` for unavailable ReactomePA and rms-based calibration, immune deconvolution, external-resource-query, and other method-specific notes. A coefficient-based ridge nomogram visualization and apparent calibration/DCA plots were generated; these are not rms bootstrap-calibrated outputs.")
writeLines(readme, file.path(RESULTS, "README.md"))
