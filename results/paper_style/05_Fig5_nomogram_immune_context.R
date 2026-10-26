source(file.path("results", "paper_style", "00_helpers.R"))

expr_df <- read_project_csv(file.path(TABLES, "GSE125498_full_expression_matrix.csv"))
meta <- read_project_csv(file.path(TABLES, "GSE125498_sample_metadata.csv"))
coef_df <- read_project_csv(file.path(TABLES, "T06_reduced_model_coefficients.csv"))
auc <- read_project_csv(file.path(TABLES, "T06_reduced_model_apparent_AUC.csv"))
cal <- read_project_csv(file.path(TABLES, "T06_discovery_calibration_bins.csv"))
dca <- read_project_csv(file.path(TABLES, "T06_discovery_decision_curve.csv"))

expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- expr_df$Gene
genes <- setdiff(coef_df$Term, "(Intercept)")
X <- t(expr[genes, meta$sample_id, drop = FALSE])
Xz <- scale(X)
beta <- setNames(coef_df$Coefficient, coef_df$Term)

point_rows <- do.call(rbind, lapply(genes, function(g) {
  z <- Xz[, g]
  b <- beta[[g]]
  data.frame(Gene = g, z_value = c(min(z), 0, max(z)),
    Linear_predictor_contribution = b * c(min(z), 0, max(z)),
    Points = 100 * b * c(min(z), 0, max(z)) / max(abs(beta[genes])))
}))
save_csv(point_rows, "Fig5A_five_gene_ridge_point_display.csv")
p_a <- ggplot(point_rows, aes(z_value, Points, colour = Gene)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  facet_wrap(~Gene, scales = "free_y") +
  labs(title = "Point display of the existing five-gene ridge model",
       subtitle = "Illustrative points scaled from the saved ridge coefficients; not a refit",
       x = "Gene-wise standardized expression", y = "Scaled model points", colour = NULL) +
  theme_paper(8)
save_plot_pair(add_letter(p_a, "A"), "Fig5A_ridge_point_display", 10, 6)

save_csv(cal, "Fig5B_calibration_bins_source.csv")
p_b <- ggplot(cal, aes(Mean_predicted, Observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = Observed_CI_low, ymax = Observed_CI_high), width = 0.015,
                colour = "#4C72B0") +
  geom_point(aes(size = n), colour = "#C44E52") +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Apparent calibration in equal-count risk bins",
       subtitle = "Observed event-rate intervals are exact binomial intervals",
       x = "Mean predicted risk", y = "Observed LPD fraction", size = "Samples") +
  theme_paper()
save_plot_pair(add_letter(p_b, "B"), "Fig5B_calibration", 6, 6)

save_csv(dca, "Fig5C_decision_curve_source.csv")
p_c <- ggplot(dca, aes(Threshold, Net_benefit, colour = Strategy)) +
  geom_line(linewidth = 0.75) +
  scale_colour_manual(values = c("Five-gene model" = "#D95F02",
    "Treat all" = "#4C72B0", "Treat none" = "#555555")) +
  labs(title = "Decision-curve analysis using apparent predictions",
       subtitle = "Treat-all crosses zero at the sample prevalence (13/33)",
       x = "Threshold probability", y = "Net benefit", colour = NULL) + theme_paper()
save_plot_pair(add_letter(p_c, "C"), "Fig5C_decision_curve", 7, 5)

prob <- plogis(as.numeric(beta["(Intercept)"] + Xz[, genes, drop = FALSE] %*% beta[genes]))
risk_thresholds <- seq(0.05, 0.95, by = 0.01)
impact <- do.call(rbind, lapply(risk_thresholds, function(th) {
  high <- prob >= th
  data.frame(Threshold = th, High_risk_n = sum(high),
             True_LPD_in_high_risk = sum(high & meta$group == "LPD"),
             Total_LPD = sum(meta$group == "LPD"))
}))
save_csv(impact, "Fig5D_clinical_impact_source.csv")
impact_long <- tidyr::pivot_longer(impact, c("High_risk_n", "True_LPD_in_high_risk"),
                                   names_to = "Measure", values_to = "Patients")
p_d <- ggplot(impact_long, aes(Threshold, Patients, colour = Measure)) +
  geom_line(linewidth = 0.75) +
  scale_colour_manual(values = c("High_risk_n" = "#4C72B0",
    "True_LPD_in_high_risk" = "#C44E52"),
    labels = c("High-risk predictions", "LPD cases in high-risk predictions")) +
  labs(title = "Clinical impact curve from apparent model predictions",
       x = "Threshold probability", y = "Patients", colour = NULL) + theme_paper()
save_plot_pair(add_letter(p_d, "D"), "Fig5D_clinical_impact", 7, 5)

micro_dir <- file.path(REV, "panel_v2_context")
cell_scores <- read_project_csv(file.path(micro_dir, "immune_celltype_ssgsea_scores.csv"))
cell_diff <- read_project_csv(file.path(micro_dir, "immune_celltype_group_comparison.csv"))
save_csv(cell_scores, "Fig5E_saved_ssgsea_sample_scores.csv")
save_csv(cell_diff, "Fig5E_saved_ssgsea_group_comparisons.csv")
cell_diff$Cell_type <- factor(cell_diff$cell_type, levels = cell_diff$cell_type[order(cell_diff$LPD_minus_SPD_mean)])
p_e <- ggplot(cell_diff, aes(LPD_minus_SPD_mean, Cell_type, colour = BH_q_across_cell_types)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_point(size = 2.2) + scale_colour_viridis_c(trans = "log10", name = "BH FDR") +
  labs(title = "Saved ssGSEA score differences by cell type",
       subtitle = "Available saved analysis contains 13 CellMarker_2024 scores, not 28 cell types",
       x = "Mean ssGSEA score difference (LPD - SPD)", y = NULL) + theme_paper()
save_plot_pair(add_letter(p_e, "E"), "Fig5E_ssgsea_celltype_differences", 8, 6)

score_cols <- setdiff(names(cell_scores), c("sample_id", "analysis_label", "title",
  "stage_binary", "stage_label", "group_binary_trait", "binary_numeric", "group"))
score_matrix <- as.matrix(cell_scores[, score_cols, drop = FALSE])
storage.mode(score_matrix) <- "double"
immune_cor <- stats::cor(score_matrix, method = "spearman", use = "pairwise.complete.obs")
immune_cor_df <- as.data.frame(as.table(immune_cor), stringsAsFactors = FALSE)
names(immune_cor_df) <- c("Cell_type_1", "Cell_type_2", "Spearman_rho")
save_csv(immune_cor_df, "Fig5F_immune_immune_spearman_source.csv")
p_f <- ggplot(immune_cor_df, aes(Cell_type_1, Cell_type_2, fill = Spearman_rho)) +
  geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                     midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Immune-score correlation", x = NULL, y = NULL, fill = "Spearman\nrho") +
  theme_paper(8) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                         panel.grid = element_blank())
save_plot_pair(add_letter(p_f, "F"), "Fig5F_immune_immune_correlation", 9, 8)

hub_panel <- read_project_csv(file.path(TABLES, "ML_3_models_hub_genes.csv"))
hub_genes <- hub_panel$Gene[hub_panel$Is_Consensus_Hub_2of3 == 1]
gene_matrix <- t(expr[hub_genes, meta$sample_id, drop = FALSE])
immune_wide <- cell_scores[match(meta$sample_id, cell_scores$sample_id), score_cols, drop = FALSE]
cor_rows <- lapply(hub_genes, function(g) {
  do.call(rbind, lapply(score_cols, function(ct) {
    test <- suppressWarnings(cor.test(gene_matrix[, g], immune_wide[[ct]], method = "spearman", exact = FALSE))
    data.frame(Gene = g, Cell_type = ct, Spearman_rho = unname(test$estimate),
               P_value = test$p.value)
  }))
})
gene_immune <- do.call(rbind, cor_rows)
gene_immune$BH_FDR <- p.adjust(gene_immune$P_value, method = "BH")
save_csv(gene_immune, "Fig5G_seven_gene_immune_spearman_source.csv")
p_g <- ggplot(gene_immune, aes(Cell_type, Gene, fill = Spearman_rho)) +
  geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                     midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Seven-gene consensus panel versus ssGSEA scores",
       x = NULL, y = NULL, fill = "Spearman\nrho") +
  theme_paper(8) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                         panel.grid = element_blank())
save_plot_pair(add_letter(p_g, "G"), "Fig5G_gene_immune_spearman", 9, 5)

write_status("Fig5", "partial", paste(
  "Reused apparent five-gene ridge model calibration and DCA. Immune source table contains",
  length(score_cols), "CellMarker_2024 ssGSEA scores; requested 28-cell ssGSEA source was unavailable."))
