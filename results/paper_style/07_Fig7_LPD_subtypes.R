source(file.path("results", "paper_style", "00_helpers.R"))

set.seed(20261016)
expr_df <- read_project_csv(file.path(TABLES, "GSE125498_full_expression_matrix.csv"))
meta <- read_project_csv(file.path(TABLES, "GSE125498_sample_metadata.csv"))
panel_df <- read_project_csv(file.path(TABLES, "ML_3_models_hub_genes.csv"))
genes <- panel_df$Gene[panel_df$Is_Consensus_Hub_2of3 == 1]
meta_lpd <- meta[meta$group == "LPD", ]
expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- expr_df$Gene
if (length(genes) != 7L || !all(genes %in% rownames(expr))) {
  stop("The saved seven-gene consensus panel is incomplete in the expression matrix.")
}
X <- t(expr[genes, meta_lpd$sample_id, drop = FALSE])
X <- scale(X)
n <- nrow(X)
k_values <- 2:6
repetitions <- 1000L
sample_fraction <- 0.80
consensus <- setNames(vector("list", length(k_values)), as.character(k_values))
metrics <- list()
for (k in k_values) {
  co_cluster <- matrix(0, n, n, dimnames = list(meta_lpd$sample_id, meta_lpd$sample_id))
  co_sampled <- matrix(0, n, n, dimnames = list(meta_lpd$sample_id, meta_lpd$sample_id))
  for (b in seq_len(repetitions)) {
    keep <- sort(sample(seq_len(n), size = max(k + 1L, floor(sample_fraction * n))))
    km <- stats::kmeans(X[keep, , drop = FALSE], centers = k, nstart = 20, iter.max = 100)
    co_sampled[keep, keep] <- co_sampled[keep, keep] + 1
    same <- outer(km$cluster, km$cluster, `==`)
    co_cluster[keep, keep] <- co_cluster[keep, keep] + same
  }
  cm <- co_cluster / pmax(co_sampled, 1)
  diag(cm) <- 1
  consensus[[as.character(k)]] <- cm
  offdiag <- cm[upper.tri(cm)]
  cdf <- data.frame(k = k, consensus_value = seq(0, 1, by = 0.01))
  cdf$cumulative_fraction <- vapply(cdf$consensus_value, function(v) mean(offdiag <= v), numeric(1))
  cdf$PAC_0_1 <- mean(offdiag > 0.1 & offdiag < 0.9)
  cdf$CDF_area <- mean(cdf$cumulative_fraction)
  metrics[[as.character(k)]] <- cdf
}
cdf_table <- do.call(rbind, metrics)
metric_summary <- unique(cdf_table[, c("k", "PAC_0_1", "CDF_area")])
metric_summary <- metric_summary[order(metric_summary$k), ]
metric_summary$Delta_area <- c(NA, diff(metric_summary$CDF_area))
best_k <- metric_summary$k[which.max(metric_summary$Delta_area)]
best_matrix <- consensus[[as.character(best_k)]]
best_cluster <- stats::cutree(stats::hclust(stats::as.dist(1 - best_matrix), method = "complete"),
                              k = best_k)
cluster_table <- data.frame(sample_id = names(best_cluster), cluster = as.integer(best_cluster))
cluster_table <- merge(cluster_table, meta_lpd[, c("sample_id", "group")], by = "sample_id", sort = FALSE)
cluster_sizes <- as.data.frame(table(cluster = cluster_table$cluster))
names(cluster_sizes)[2] <- "sample_count"
cluster_sizes$selected_k <- best_k
save_csv(cdf_table, "Fig7A_consensus_CDF_source.csv")
save_csv(metric_summary, "Fig7A_consensus_metrics.csv")
save_csv(as.data.frame(as.table(best_matrix)), "Fig7A_best_k_consensus_matrix.csv")
save_csv(cluster_table, "Fig7_cluster_assignments.csv")
save_csv(cluster_sizes, "Fig7_cluster_sizes.csv")

p_cdf <- ggplot(cdf_table, aes(consensus_value, cumulative_fraction, colour = factor(k))) +
  geom_line(linewidth = 0.65) +
  labs(title = "Consensus-clustering CDF", x = "Consensus index",
       y = "Cumulative fraction", colour = "k") + theme_paper()
save_plot_pair(add_letter(p_cdf, "A"), "Fig7A_consensus_CDF", 7, 5)

p_delta <- ggplot(metric_summary[!is.na(metric_summary$Delta_area), ], aes(k, Delta_area)) +
  geom_col(fill = "#2C7FB8") +
  geom_line(aes(group = 1), colour = "#16324F") +
  geom_point(colour = "#16324F") +
  labs(title = "Consensus CDF area change", x = "Number of clusters (k)",
       y = "Change in CDF area") + theme_paper()
save_plot_pair(p_delta, "Fig7A_consensus_delta_area", 6, 5)

p_pac <- ggplot(metric_summary, aes(k, PAC_0_1)) +
  geom_line(colour = "#16324F") + geom_point(size = 2, colour = "#2C7FB8") +
  labs(title = "Proportion of ambiguously clustered pairs", x = "Number of clusters (k)",
       y = "PAC (consensus 0.1–0.9)") + theme_paper()
save_csv(metric_summary[, c("k", "PAC_0_1")], "Fig7B_PAC_source.csv")
save_plot_pair(add_letter(p_pac, "B"), "Fig7B_ambiguous_clustering", 6, 5)

matrix_long <- as.data.frame(as.table(best_matrix))
names(matrix_long) <- c("sample_1", "sample_2", "consensus")
save_csv(matrix_long, "Fig7A_best_k_consensus_heatmap_source.csv")
p_matrix <- ggplot(matrix_long, aes(sample_1, sample_2, fill = consensus)) +
  geom_tile() + scale_fill_viridis_c(limits = c(0, 1), name = "Consensus") +
  labs(title = paste("Consensus matrix at selected k =", best_k), x = NULL, y = NULL) +
  theme_paper(8) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                         panel.grid = element_blank())
save_plot_pair(p_matrix, "Fig7A_best_k_consensus_matrix", 7, 6)

pca <- stats::prcomp(X, center = FALSE, scale. = FALSE)
pca_df <- data.frame(sample_id = rownames(X), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df <- merge(pca_df, cluster_table, by = "sample_id", sort = FALSE)
pca_df$cluster <- factor(pca_df$cluster)
save_csv(pca_df, "Fig7B_LPD_panel_PCA_source.csv")
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = cluster, label = sample_id)) +
  geom_point(size = 2.8) + geom_text(nudge_y = 0.15, size = 2.2, show.legend = FALSE) +
  labs(title = "PCA of the seven-gene panel in LPD samples",
       subtitle = paste("Cluster assignments from consensus matrix; selected k =", best_k),
       colour = "Cluster") + theme_paper()
save_plot_pair(add_letter(p_pca, "B"), "Fig7B_LPD_panel_PCA", 7, 6)

gene_long <- do.call(rbind, lapply(genes, function(g) {
  data.frame(sample_id = rownames(X), Gene = g, standardized_expression = X[, g])
}))
gene_long <- merge(gene_long, cluster_table[, c("sample_id", "cluster")], by = "sample_id")
save_csv(gene_long, "Fig7C_panel_expression_by_cluster_source.csv")
p_gene <- ggplot(gene_long, aes(factor(cluster), standardized_expression, fill = factor(cluster))) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.1, size = 1.1) + facet_wrap(~Gene, ncol = 4) +
  labs(title = "Consensus-panel expression by LPD cluster", x = "Cluster",
       y = "Gene-wise standardized expression", fill = "Cluster") + theme_paper(8)
save_plot_pair(add_letter(p_gene, "C"), "Fig7C_panel_expression_by_cluster", 10, 7)

cell_scores <- read_project_csv(file.path(REV, "panel_v2_context", "immune_celltype_ssgsea_scores.csv"))
score_cols <- setdiff(names(cell_scores), c("sample_id", "analysis_label", "title", "stage_binary",
  "stage_label", "group_binary_trait", "binary_numeric", "group"))
score_lpd <- merge(cell_scores[cell_scores$sample_id %in% meta_lpd$sample_id, ],
                   cluster_table[, c("sample_id", "cluster")], by = "sample_id")
cell_tests <- do.call(rbind, lapply(score_cols, function(ct) {
  x <- score_lpd[[ct]]
  groups <- split(x, score_lpd$cluster)
  if (length(groups) < 2L || any(lengths(groups) < 2L)) {
    return(data.frame(cell_type = ct, P_value = NA_real_, direction = "not testable",
                      cluster_medians = NA_character_))
  }
  medians <- vapply(groups, median, numeric(1))
  test <- if (length(groups) == 2L) {
    stats::wilcox.test(groups[[1]], groups[[2]], exact = FALSE)
  } else stats::kruskal.test(x ~ factor(score_lpd$cluster))
  data.frame(cell_type = ct, P_value = test$p.value,
    direction = paste0("highest median in cluster ", names(which.max(medians))),
    cluster_medians = paste(names(medians), sprintf("%.6g", medians), sep = ":", collapse = ";"))
}))
cell_tests$BH_FDR <- p.adjust(cell_tests$P_value, method = "BH")
save_csv(cell_tests, "Fig7D_ssgsea_cluster_comparisons.csv")
long_scores <- tidyr::pivot_longer(score_lpd, all_of(score_cols), names_to = "cell_type", values_to = "ssgsea_score")
save_csv(long_scores, "Fig7D_ssgsea_cluster_scores_source.csv")
p_cell <- ggplot(long_scores, aes(factor(cluster), ssgsea_score, fill = factor(cluster))) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.1, size = 0.7) + facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
  labs(title = "Saved relative ssGSEA scores by LPD cluster",
       subtitle = "13 available scores; exploratory comparisons with BH correction",
       x = "Cluster", y = "Relative ssGSEA score", fill = "Cluster") + theme_paper(8)
save_plot_pair(add_letter(p_cell, "D"), "Fig7D_ssgsea_by_cluster", 11, 8)

write_status("Fig7", "complete", paste(
  "Custom consensus resampling on", n, "LPD samples; 1,000 random 80% sample resamples per k = 2–6.",
  "Selected k =", best_k, "by the largest incremental CDF area; cluster sizes:",
  paste(cluster_sizes$sample_count, collapse = ", "),
  ". Clustering is exploratory; package ConsensusClusterPlus was unavailable, so base-R implementation was used."))
