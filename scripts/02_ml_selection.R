args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("ggplot2", "dplyr", "tidyr"), "T02")

path <- file.path(RESULTS, "revision", "ml_39_matrisome_wgcna", "ml_39_gene_results.csv")
if (!file.exists(path)) stop("Saved 39-feature ML selection table is missing: ", path)
x <- read_csv_safe(path)
flags <- c("LASSO_selected", "SVM_RFE_selected", "RF_selected")
if (!all(flags %in% names(x))) stop("The ML result table lacks expected selection flags.")
sets <- c(LASSO = sum(x[[flags[1]]] == 1), SVM_RFE = sum(x[[flags[2]]] == 1), Random_Forest = sum(x[[flags[3]]] == 1))
pattern <- apply(x[, flags], 1, paste0, collapse = "")
upset <- as.data.frame(table(pattern), stringsAsFactors = FALSE)
names(upset) <- c("Selection_pattern_LASSO_SVM_RF", "Gene_count")
upset <- upset[upset$Gene_count > 0, ]
upset$Genes <- vapply(upset$Selection_pattern_LASSO_SVM_RF, function(pat) paste(x$Gene[pattern == pat], collapse = ";"), character(1))
write.csv(upset, file.path(TABLES, "T02_ML_selection_intersections.csv"), row.names = FALSE)
write.csv(x, file.path(TABLES, "T02_ML_selection_reproduced_from_saved_results.csv"), row.names = FALSE)
p <- ggplot2::ggplot(upset, ggplot2::aes(reorder(Selection_pattern_LASSO_SVM_RF, Gene_count), Gene_count)) +
  ggplot2::geom_col(fill = "#0072B2") + ggplot2::coord_flip() +
  ggplot2::labs(title = "UpSet-style ML selection intersections", x = "LASSO / SVM-RFE / RF (1 = selected)", y = "Gene count") + theme_project()
save_plot(p, "T02_ML_selection_intersections", 8, 5)

coef <- x[x$LASSO_selected == 1, c("Gene", "LASSO_abs_coefficient"), drop = FALSE]
p_lasso <- ggplot2::ggplot(coef, ggplot2::aes(reorder(Gene, LASSO_abs_coefficient), LASSO_abs_coefficient)) +
  ggplot2::geom_col(fill = "#0072B2") + ggplot2::coord_flip() +
  ggplot2::labs(title = "Saved LASSO selected coefficient magnitudes", subtitle = "Coefficient path and CV curves were not saved; not reconstructed from selected coefficients.",
                x = NULL, y = "Absolute selected coefficient") + theme_project()
save_plot(p_lasso, "T02_LASSO_selected_coefficients", 7, 5)
svm <- x[x$SVM_RFE_selected == 1, c("Gene", "SVM_RFE_rank"), drop = FALSE]
p_svm <- ggplot2::ggplot(svm, ggplot2::aes(reorder(Gene, -SVM_RFE_rank), SVM_RFE_rank)) +
  ggplot2::geom_segment(ggplot2::aes(xend = Gene, y = 0, yend = SVM_RFE_rank), colour = "#56B4E9") +
  ggplot2::geom_point(size = 2, colour = "#0072B2") +
  ggplot2::labs(title = "Saved SVM-RFE feature ranks", subtitle = "Error-versus-feature-count curve unavailable in saved results.",
                x = NULL, y = "Saved SVM-RFE rank") + theme_project()
save_plot(p_svm, "T02_SVM_RFE_saved_ranks", 7, 5)
rf <- x[x$RF_selected == 1, c("Gene", "RF_importance"), drop = FALSE]
p_rf <- ggplot2::ggplot(rf, ggplot2::aes(reorder(Gene, RF_importance), RF_importance)) +
  ggplot2::geom_segment(ggplot2::aes(xend = Gene, y = 0, yend = RF_importance), colour = "#CC79A7") +
  ggplot2::geom_point(size = 2, colour = "#882255") + ggplot2::coord_flip() +
  ggplot2::labs(title = "Saved random-forest selected-feature importance", x = NULL, y = "Saved RF importance") + theme_project()
save_plot(p_rf, "T02_RF_importance", 7, 5)
task_log("T02_ml_selection", c(paste("Loaded saved results; did not repeat feature selection.", paste(names(sets), sets, sep = "=", collapse = "; ")),
  paste("Selection union:", sum(rowSums(x[, flags]) > 0), "; consensus >=2 methods:", sum(x$Consensus_At_Least_2_of_3 == 1), "; strict 3/3:", sum(rowSums(x[, flags]) == 3)),
  "LASSO coefficient path/CV diagnostics and SVM-RFE error curve were not present in saved files and were not fabricated."))
append_results_text("T2. ML selection", paste("The saved 39-gene model results were summarized without rerunning feature selection.",
  "LASSO, SVM-RFE, and random forest had 11, 10, and 10 selected genes, respectively; the saved results contained 20 unique selections and nine >=2-method consensus genes.",
  "Coefficient-path, cross-validation, and SVM-RFE error curves could not be reproduced because their intermediate model objects were unavailable."))
