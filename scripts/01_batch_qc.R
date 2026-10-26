args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("ggplot2", "tidyr"), "T01")

bundles <- load_all_cohorts()
log <- character()
for (b in bundles) {
  expr <- b$expr
  meta <- b$meta
  vars <- apply(expr, 1, stats::var, na.rm = TRUE)
  selected <- names(sort(vars, decreasing = TRUE))[seq_len(min(1000L, sum(is.finite(vars))))]
  pca <- stats::prcomp(t(expr[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
  score <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  score$sample_id <- rownames(score)
  score <- merge(score, meta[, c("sample_id", "Group")], by = "sample_id")
  p <- ggplot2::ggplot(score, ggplot2::aes(PC1, PC2, colour = Group, label = sample_id)) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(check_overlap = TRUE, vjust = -0.7, size = 2.7, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = GROUP_COLORS) +
    ggplot2::labs(title = paste(b$cohort, "PCA"), subtitle = paste(b$scale, "| top variable genes:", length(selected)),
                  x = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
                  y = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]), colour = "Group") +
    theme_project()
  save_plot(p, paste0("T01_", b$cohort, "_PCA"), 7, 5)

  distance <- as.matrix(stats::dist(t(expr[selected, , drop = FALSE])))
  dist_df <- as.data.frame(as.table(distance))
  names(dist_df) <- c("Sample1", "Sample2", "Distance")
  p_dist <- ggplot2::ggplot(dist_df, ggplot2::aes(Sample1, Sample2, fill = Distance)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::labs(title = paste(b$cohort, "sample distance"), x = NULL, y = NULL, fill = "Euclidean\n distance") +
    theme_project(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  save_plot(p_dist, paste0("T01_", b$cohort, "_sample_distance"), 8, 7)

  q <- t(apply(expr, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE))
  box <- data.frame(sample_id = rownames(q), p05 = q[, 1], q1 = q[, 2], median = q[, 3], q3 = q[, 4], p95 = q[, 5])
  box <- merge(box, meta[, c("sample_id", "Group")], by = "sample_id")
  p_box <- ggplot2::ggplot(box, ggplot2::aes(sample_id, median, colour = Group)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = p05, ymax = p95), linewidth = 0.6) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = q1, ymax = q3), linewidth = 2.1) +
    ggplot2::geom_point(size = 1.8) + ggplot2::scale_colour_manual(values = GROUP_COLORS) +
    ggplot2::labs(title = paste(b$cohort, "normalized expression by sample"), subtitle = b$scale,
                  x = "Sample", y = "Expression (5th–95th percentile; thick bar = IQR)", colour = "Group") +
    theme_project() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  save_plot(p_box, paste0("T01_", b$cohort, "_expression_QC"), 10, 5)
  write.csv(box, file.path(TABLES, paste0("T01_", b$cohort, "_expression_QC.csv")), row.names = FALSE)
  log <- c(log, paste0(b$cohort, ": ", ncol(expr), " samples; ", nrow(expr), " genes; expression source ", b$scale,
                       "; no cohorts were merged, so no batch correction was performed."))
}
task_log("T01_batch_qc", log)
append_results_text("T1. Batch/QC", paste("PCA, sample-distance heatmaps, and per-sample expression summaries were generated separately for GSE125498, GSE248762, and GSE130888.",
  "The cohorts were not merged, so no batch-correction analysis was performed.",
  "GSE125498 used the saved RMA log2 matrix; effluent cohorts used the saved edgeR TMM log2-CPM matrices."))
