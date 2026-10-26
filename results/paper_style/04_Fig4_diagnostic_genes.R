source(file.path("results", "paper_style", "00_helpers.R"))

ml39 <- read_project_csv(file.path(REV, "ml_39_matrisome_wgcna", "ml_39_gene_results.csv"))
patterns <- read_project_csv(file.path(TABLES, "T02_ML_selection_intersections.csv"))
final_votes <- read_project_csv(file.path(TABLES, "ML_3_models_hub_genes.csv"))
final_panel <- final_votes$Gene[final_votes$Is_Consensus_Hub_2of3 == 1]
five_panel <- c("TNFSF15", "ADAM19", "FLT3LG", "NAT6", "ECM1")
stopifnot(length(final_panel) == 7L)

method_counts <- data.frame(
  Method = c("LASSO", "SVM-RFE", "Random forest"),
  Selected_genes = c(sum(ml39$LASSO_selected == 1),
                     sum(ml39$SVM_RFE_selected == 1),
                     sum(ml39$RF_selected == 1))
)
save_csv(method_counts, "Fig4A-C_saved_selector_counts.csv")

panel_provenance <- data.frame(
  Gene = sort(unique(c(final_panel, five_panel))),
  In_seven_gene_consensus = sort(unique(c(final_panel, five_panel))) %in% final_panel,
  In_five_gene_ridge_model = sort(unique(c(final_panel, five_panel))) %in% five_panel,
  Selection_source = "ML_3_models_hub_genes.csv; majority vote in saved broader-pool screen"
)
save_csv(panel_provenance, "Fig4_final_panel_provenance.csv")
save_csv(final_votes[final_votes$Is_Consensus_Hub_2of3 == 1, ], "Fig4_seven_gene_consensus_panel.csv")
save_csv(ml39, "Fig4_39_gene_pool_selector_results_reused.csv")
save_csv(patterns, "Fig4_39_gene_pool_intersections_reused.csv")

counts_plot <- ggplot(method_counts, aes(Method, Selected_genes, fill = Method)) +
  geom_col(show.legend = FALSE) + geom_text(aes(label = Selected_genes), vjust = -0.35) +
  labs(title = "Saved feature-selector counts in the 39-gene pool",
       subtitle = "Selection curves and fold-level fits were not retained",
       x = NULL, y = "Genes selected") + theme_paper()
save_plot_pair(counts_plot, "Fig4_selector_counts_supplement", 6, 5)

p_lasso <- ggplot(ml39[ml39$LASSO_selected == 1, ],
                  aes(reorder(Gene, LASSO_abs_coefficient), LASSO_abs_coefficient)) +
  geom_col(fill = "#2C7FB8") + coord_flip() +
  labs(title = "Saved LASSO-selected coefficients", x = NULL,
       y = "Absolute coefficient (saved fit)") + theme_paper()
save_plot_pair(add_letter(p_lasso, "A"), "Fig4A_lasso_selected_coefficients", 7, 6)

p_svm <- ggplot(ml39[ml39$SVM_RFE_selected == 1, ],
                aes(reorder(Gene, -SVM_RFE_rank), SVM_RFE_rank)) +
  geom_col(fill = "#41AB5D") + coord_flip() +
  labs(title = "Saved SVM-RFE selected-feature ranks", x = NULL,
       y = "Rank (smaller is earlier)") + theme_paper()
save_plot_pair(add_letter(p_svm, "B"), "Fig4B_svm_rfe_selected_ranks", 7, 6)

p_rf <- ggplot(ml39[ml39$RF_selected == 1, ],
               aes(reorder(Gene, RF_importance), RF_importance)) +
  geom_col(fill = "#DD8452") + coord_flip() +
  labs(title = "Saved random-forest feature importance",
       subtitle = "The original error curve was not retained",
       x = NULL, y = "Mean decrease in Gini") + theme_paper()
save_plot_pair(add_letter(p_rf, "C"), "Fig4C_random_forest_importance", 7, 6)

pattern_plot <- patterns[patterns$Gene_count > 0, ]
p_venn <- ggplot(pattern_plot,
                 aes(reorder(Selection_pattern_LASSO_SVM_RF, Gene_count), Gene_count)) +
  geom_col(fill = "#7A5195") + geom_text(aes(label = Gene_count), hjust = -0.2) +
  coord_flip() + labs(title = "Three-selector intersections in the saved 39-gene pool",
    subtitle = "Pattern order: LASSO / SVM-RFE / random forest",
    x = "Membership pattern", y = "Genes") + theme_paper()
save_plot_pair(add_letter(p_venn, "D"), "Fig4D_selector_intersections", 8, 6)
selector_sets <- list(
LASSO = ml39$Gene[ml39$LASSO_selected == 1],
SVM_RFE = ml39$Gene[ml39$SVM_RFE_selected == 1],
Random_forest = ml39$Gene[ml39$RF_selected == 1]
)
selector_membership <- data.frame(Gene = sort(unique(unlist(selector_sets))))
for (method in names(selector_sets)) {
selector_membership[[method]] <- selector_membership$Gene %in% selector_sets[[method]]
}
selector_membership$Membership_code <- apply(
selector_membership[, names(selector_sets), drop = FALSE],
1, function(x) paste(as.integer(x), collapse = "")
)
save_csv(selector_membership, "Fig4D_selector_venn_source.csv")
venn_grob <- VennDiagram::venn.diagram(
x = selector_sets, filename = NULL, fill = c("#2C7FB8", "#41AB5D", "#DD8452"),
alpha = 0.42, cat.col = c("#1E5A84", "#2C783E", "#A64B1C"),
cex = 1.1, cat.cex = 0.9, margin = 0.08
)
venn_png <- file.path(OUT, "Fig4D_selector_venn.png")
venn_pdf <- file.path(OUT, "Fig4D_selector_venn.pdf")
new_path(venn_png); new_path(venn_pdf)
grDevices::png(venn_png, width = 8, height = 6, units = "in", res = 300, bg = "white")
grid::grid.newpage(); grid::grid.draw(venn_grob)
grid::grid.text("D", x = 0.02, y = 0.98, just = c("left", "top"),
              gp = grid::gpar(fontface = "bold", fontsize = 14, col = "#16324F"))
grDevices::dev.off()
grDevices::cairo_pdf(venn_pdf, width = 8, height = 6, bg = "white")
grid::grid.newpage(); grid::grid.draw(venn_grob)
grid::grid.text("D", x = 0.02, y = 0.98, just = c("left", "top"),
              gp = grid::gpar(fontface = "bold", fontsize = 14, col = "#16324F"))
grDevices::dev.off()

expr_df <- read_project_csv(file.path(TABLES, "GSE125498_full_expression_matrix.csv"))
meta <- read_project_csv(file.path(TABLES, "GSE125498_sample_metadata.csv"))
coef_df <- read_project_csv(file.path(TABLES, "T06_reduced_model_coefficients.csv"))
auc_saved <- read_project_csv(file.path(TABLES, "T06_reduced_model_apparent_AUC.csv"))
expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- expr_df$Gene
X <- t(expr[intersect(five_panel, rownames(expr)), meta$sample_id, drop = FALSE])
X <- scale(X)
beta <- setNames(coef_df$Coefficient, coef_df$Term)
needed <- setdiff(names(beta), "(Intercept)")
if (!all(needed %in% colnames(X))) stop("Reduced model predictors do not match the saved expression matrix.")
lp <- as.numeric(beta["(Intercept)"] + X[, needed, drop = FALSE] %*% beta[needed])
prob <- plogis(lp)
roc_data <- data.frame(sample_id = meta$sample_id, group = meta$group,
                       case = as.integer(meta$group == "LPD"), predicted_probability = prob)
thresholds <- sort(unique(c(Inf, prob, -Inf)), decreasing = TRUE)
roc_points <- do.call(rbind, lapply(thresholds, function(th) {
  call <- prob >= th
  data.frame(threshold = th, TPR = sum(call & roc_data$case == 1) / sum(roc_data$case == 1),
             FPR = sum(call & roc_data$case == 0) / sum(roc_data$case == 0))
}))
roc_points <- roc_points[order(roc_points$FPR, roc_points$TPR), ]
roc_points$Cohort <- "GSE125498 apparent five-gene ridge"
save_csv(roc_data, "Fig4E_training_model_predictions.csv")
save_csv(roc_points, "Fig4E_training_ROC_points.csv")
p_roc <- ggplot(roc_points, aes(FPR, TPR)) +
  geom_step(colour = "#D95F02", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  coord_equal() +
  labs(title = "Apparent training ROC: existing five-gene ridge model",
       subtitle = sprintf("AUC %.3f (bootstrap 95%% CI %.3f–%.3f)",
         auc_saved$AUC[1], auc_saved$Bootstrap_low[1], auc_saved$Bootstrap_high[1]),
       x = "False-positive rate", y = "True-positive rate") + theme_paper()
save_plot_pair(add_letter(p_roc, "E"), "Fig4E_training_ROC", 6, 6)

box_source <- do.call(rbind, lapply(intersect(final_panel, rownames(expr)), function(g) {
  data.frame(Gene = g, sample_id = meta$sample_id, expression = as.numeric(expr[g, meta$sample_id]),
             group = meta$group)
}))
save_csv(box_source, "Fig4G_discovery_expression_source.csv")
p_box <- ggplot(box_source, aes(group, expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.7) +
  facet_wrap(~Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(LPD = "#C44E52", SPD = "#4C72B0")) +
  labs(title = "Seven-gene consensus panel in the discovery cohort",
       x = NULL, y = "Normalized expression") + theme_paper(8)
save_plot_pair(add_letter(p_box, "G"), "Fig4G_discovery_expression", 10, 7)

validation <- read_project_csv(file.path(TABLES, "T06_additional_cohort_model_performance.csv"))
save_csv(validation, "Fig4F_saved_additional_cohort_AUC_summaries.csv")
write_status("Fig4", "partial", paste(
  "Saved selector/intersection tables reused. Curves for LASSO CV, SVM-RFE CV and RF error were not saved.",
  "Training ROC was regenerated from the saved five-gene ridge coefficients and expression matrix.",
  "Validation ROC coordinates and matching seven-gene validation expression matrices were unavailable;",
  nrow(validation), "additional-cohort AUC summaries use a different nine-gene model and were not presented as this panel's validation."))
