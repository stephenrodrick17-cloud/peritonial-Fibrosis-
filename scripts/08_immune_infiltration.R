args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("ggplot2", "dplyr"), "T08")

old <- file.path(TABLES, "immune_infiltration_scores.csv")
if (!file.exists(old)) {
  reason <- "No deconvolution package or compatible cell-signature reference was available."
  write.csv(data.frame(status = "not_run", reason = reason),
    file.path(TABLES, "T08_immune_deconvolution_status.csv"), row.names = FALSE)
  task_log("T08_immune_infiltration", reason, "skipped")
  append_results_text("T8. Immune infiltration", "No compatible deconvolution package or reference matrix was available, so immune-cell fractions were not inferred. The saved myofibroblast score was not treated as immune deconvolution.")
} else {
  score <- read_csv_safe(old)
  if (!all(c("sample_id", "Myofibroblast_Score") %in% names(score))) {
    reason <- "The saved score table contained no compatible immune cell-type fractions; only myofibroblast scores were available."
    write.csv(data.frame(status = "not_run", reason = reason), file.path(TABLES, "T08_immune_deconvolution_status.csv"), row.names = FALSE)
    task_log("T08_immune_infiltration", reason, "skipped")
    append_results_text("T8. Immune infiltration", "No compatible immune-cell fractions or two-method deconvolution outputs were available. The saved myofibroblast score was kept distinct and was not interpreted as immune infiltration.")
  } else {
    # This output is explicitly a stromal score, not immune deconvolution.
    meta <- load_discovery()$meta
    d <- merge(score, meta[, c("sample_id", "Group")], by = "sample_id")
    d$Group <- factor(d$Group, levels = c("Control", "Case"))
    test <- stats::wilcox.test(Myofibroblast_Score ~ Group, d, exact = FALSE)
    write.csv(data.frame(Measure = "Myofibroblast_Score", n_control = sum(d$Group == "Control"), n_case = sum(d$Group == "Case"),
      P = test$p.value, FDR = p.adjust(test$p.value, "BH"), Interpretation = "stromal score, not immune deconvolution"),
      file.path(TABLES, "T08_myofibroblast_score_only.csv"), row.names = FALSE)
    p <- ggplot2::ggplot(d, ggplot2::aes(Group, Myofibroblast_Score, fill = Group)) +
      ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.12) +
      ggplot2::scale_fill_manual(values = GROUP_COLORS) +
      ggplot2::labs(title = "Saved myofibroblast score (not immune deconvolution)", x = NULL, y = "Score") + theme_project()
    save_plot(p, "T08_myofibroblast_score_only", 6, 5)
    reason <- "A myofibroblast score was available, but neither CIBERSORT-compatible references nor a second deconvolution method were available."
    task_log("T08_immune_infiltration", reason, "skipped")
    append_results_text("T8. Immune infiltration", "The available stromal score was summarized separately and was not treated as an immune-cell fraction. No two-method immune deconvolution was performed because compatible reference data and a second method were unavailable.")
  }
}
