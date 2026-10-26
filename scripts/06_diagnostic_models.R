args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("glmnet", "ggplot2", "dplyr", "tidyr"), "T06")

B <- 1000L
set.seed(2026)
disc <- load_discovery()
genes <- intersect(HUBS, rownames(disc$expr))
X0 <- t(disc$expr[genes, , drop = FALSE])
y <- as.integer(disc$meta$stage_binary)
zscore <- function(x) {
  mu <- colMeans(x)
  sig <- apply(x, 2, stats::sd)
  keep <- is.finite(sig) & sig > 0
  list(x = scale(x[, keep, drop = FALSE], center = mu[keep], scale = sig[keep]), center = mu[keep], scale = sig[keep])
}
genes <- colnames(X0)
scale_train_test <- function(train_x, test_x) {
  center <- colMeans(train_x)
  spread <- apply(train_x, 2, stats::sd)
  keep <- is.finite(spread) & spread > 0
  if (!any(keep)) stop("No non-constant predictors in the training split.")
  train_scaled <- scale(train_x[, keep, drop = FALSE], center = center[keep], scale = spread[keep])
  test_scaled <- scale(test_x[, keep, drop = FALSE], center = center[keep], scale = spread[keep])
  list(train = unclass(train_scaled), test = unclass(test_scaled), genes = names(center[keep]))
}
make_foldid <- function(labels, k) {
  foldid <- integer(length(labels))
  for (cl in sort(unique(labels))) {
    ii <- which(labels == cl)
    foldid[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  }
  foldid
}
fit_tuned_ridge <- function(train_x, train_y, test_x) {
  k <- min(5L, min(table(train_y)))
  if (k < 3L) return(list(train = rep(NA_real_, length(train_y)), test = rep(NA_real_, nrow(test_x))))
  scaled <- scale_train_test(train_x, test_x)
  foldid <- make_foldid(train_y, k)
  cv_fit <- glmnet::cv.glmnet(scaled$train, train_y, family = "binomial", alpha = 0,
    foldid = foldid, type.measure = "deviance", standardize = FALSE)
  fit <- glmnet::glmnet(scaled$train, train_y, family = "binomial", alpha = 0,
    lambda = cv_fit$lambda.min, standardize = FALSE)
  list(train = as.numeric(stats::predict(fit, scaled$train, type = "response")),
    test = as.numeric(stats::predict(fit, scaled$test, type = "response")))
}
full_scaled <- scale_train_test(X0, X0)
X <- full_scaled$train
cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 0, foldid = make_foldid(y, 5L),
                        type.measure = "deviance", standardize = FALSE)
lambda <- cv$lambda.min
fit <- glmnet::glmnet(X, y, family = "binomial", alpha = 0, lambda = lambda, standardize = FALSE)
pred <- as.numeric(stats::predict(fit, X, type = "response"))

delong_ci <- function(labels, scores) {
  pos <- scores[labels == 1]; neg <- scores[labels == 0]
  z <- (outer(pos, neg, ">") + 0.5 * outer(pos, neg, "=="))
  a <- mean(z)
  se <- sqrt(stats::var(rowMeans(z)) / length(pos) + stats::var(colMeans(z)) / length(neg))
  c(max(0, a - 1.96 * se), min(1, a + 1.96 * se))
}
auc_ci_row <- function(name, labels, scores, B = 2000L) {
  a <- auc_rank(labels, scores)
  ci <- bootstrap_auc(labels, scores, B)
  dci <- delong_ci(labels, scores)
  data.frame(Model = name, AUC = a, Bootstrap_low = ci[1], Bootstrap_high = ci[2],
             DeLong_low = dci[1], DeLong_high = dci[2], n_case = sum(labels == 1), n_control = sum(labels == 0))
}
apparent <- auc_ci_row("GSE125498 apparent nine-gene ridge", y, pred)

optimism <- numeric(B)
selected <- matrix(FALSE, nrow = B, ncol = length(genes), dimnames = list(NULL, genes))
for (b in seq_len(B)) {
  idx <- sample(seq_along(y), replace = TRUE)
  if (length(unique(y[idx])) < 2 || min(table(y[idx])) < 3L) { optimism[b] <- NA_real_; next }
  boot_fit <- fit_tuned_ridge(X0[idx, , drop = FALSE], y[idx], rbind(X0[idx, , drop = FALSE], X0))
  n_boot <- length(idx)
  optimism[b] <- auc_rank(y[idx], boot_fit$train) -
    auc_rank(y, boot_fit$test[(n_boot + 1L):(n_boot + length(y))])
  boot_scaled <- scale_train_test(X0[idx, , drop = FALSE], X0[idx, , drop = FALSE])
  lasso_cv <- glmnet::cv.glmnet(boot_scaled$train, y[idx], family = "binomial", alpha = 1,
    foldid = make_foldid(y[idx], min(5L, min(table(y[idx])))), type.measure = "deviance", standardize = FALSE)
  lasso <- glmnet::glmnet(boot_scaled$train, y[idx], family = "binomial", alpha = 1,
    lambda = lasso_cv$lambda.1se, standardize = FALSE)
  co <- as.matrix(stats::coef(lasso))[-1, 1]
  selected[b, names(co)] <- abs(co) > 1e-8
}
optimism_corrected <- apparent$AUC - mean(optimism, na.rm = TRUE)
apparent$Optimism_corrected_AUC <- optimism_corrected
apparent$Mean_bootstrap_optimism <- mean(optimism, na.rm = TRUE)
apparent$Correction_scope <- "Model-fitting/tuning conditional on the preselected nine-gene panel; panel-selection optimism was not corrected"
write.csv(apparent, file.path(TABLES, "T06_discovery_apparent_optimism_corrected_AUC.csv"), row.names = FALSE)
write.csv(data.frame(replicate = seq_along(optimism), optimism = optimism), file.path(TABLES, "T06_bootstrap_optimism_replicates.csv"), row.names = FALSE)

freq <- colMeans(selected)
stability <- data.frame(Gene = names(freq), Selection_frequency = as.numeric(freq))
stability <- stability[order(-stability$Selection_frequency, stability$Gene), ]
write.csv(stability, file.path(TABLES, "T06_lasso_bootstrap_stability.csv"), row.names = FALSE)
top_genes <- head(stability$Gene, min(5L, length(genes)))
reduced <- glmnet::cv.glmnet(X[, top_genes, drop = FALSE], y, family = "binomial", alpha = 0,
                             nfolds = 5, type.measure = "deviance", standardize = FALSE)
reduced_fit <- glmnet::glmnet(X[, top_genes, drop = FALSE], y, family = "binomial",
                              alpha = 0, lambda = reduced$lambda.min, standardize = FALSE)
reduced_pred <- as.numeric(stats::predict(reduced_fit, X[, top_genes, drop = FALSE], type = "response"))
reduced_apparent <- auc_ci_row("GSE125498 apparent reduced ridge", y, reduced_pred)
reduced_apparent$Genes_used <- paste(top_genes, collapse = ";")
write.csv(reduced_apparent, file.path(TABLES, "T06_reduced_model_apparent_AUC.csv"), row.names = FALSE)
reduced_coef <- as.matrix(stats::coef(reduced_fit))[, 1]
write.csv(data.frame(Term = names(reduced_coef), Coefficient = unname(reduced_coef)),
  file.path(TABLES, "T06_reduced_model_coefficients.csv"), row.names = FALSE)

predictor_range <- t(vapply(top_genes, function(gene) {
  z <- X[, gene]
  c(min = min(z), max = max(z), beta = reduced_coef[[gene]])
}, numeric(3)))
lp_ranges <- predictor_range[, "beta"] * (predictor_range[, "max"] - predictor_range[, "min"])
max_lp_range <- max(abs(lp_ranges))
points_per_lp <- 100 / max_lp_range
predictor_rows <- lapply(seq_along(top_genes), function(i) {
  gene <- top_genes[[i]]
  z <- X[, gene]
  beta <- reduced_coef[[gene]]
  values <- seq(min(z), max(z), length.out = 5)
  data.frame(Gene = gene, z_value = values,
    Points = (beta * values - min(beta * c(min(z), max(z)))) * points_per_lp,
    Row = length(top_genes) - i + 2L)
})
nomogram_ticks <- do.call(rbind, predictor_rows)
total_points_max <- sum(abs(lp_ranges)) * points_per_lp
base_lp <- unname(reduced_coef[["(Intercept)"]]) +
  sum(pmin(predictor_range[, "beta"] * predictor_range[, "min"],
           predictor_range[, "beta"] * predictor_range[, "max"]))
score_ticks <- seq(0, total_points_max, length.out = 6)
score_rows <- data.frame(Points = score_ticks, Row = 1, Label = sprintf("%.0f", score_ticks))
points_axis <- data.frame(Points = seq(0, 100, by = 25), Row = 1, Label = as.character(seq(0, 100, by = 25)))
prob_ticks <- seq(0.05, 0.95, by = 0.15)
prob_points <- (stats::qlogis(prob_ticks) - base_lp) * points_per_lp
prob_keep <- prob_points >= 0 & prob_points <= total_points_max
prob_rows <- data.frame(Points = prob_points[prob_keep], Row = 0,
  Label = sprintf("%.2f", prob_ticks[prob_keep]))
axis_segments <- rbind(
  data.frame(x = 0, xend = 100, y = 1, yend = 1),
  data.frame(x = nomogram_ticks$Points, xend = nomogram_ticks$Points,
    y = nomogram_ticks$Row, yend = nomogram_ticks$Row),
  data.frame(x = 0, xend = total_points_max, y = 0, yend = 0),
  data.frame(x = 0, xend = total_points_max, y = -1, yend = -1)
)
nomogram_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = axis_segments, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45, colour = "#17324D") +
  ggplot2::geom_segment(data = points_axis, ggplot2::aes(x = Points, xend = Points,
    y = Row - 0.08, yend = Row + 0.08), linewidth = 0.35) +
  ggplot2::geom_text(data = points_axis, ggplot2::aes(Points, Row + 0.17, label = Label),
    size = 2.8, vjust = 0) +
  ggplot2::geom_segment(data = nomogram_ticks, ggplot2::aes(x = Points, xend = Points,
    y = Row - 0.08, yend = Row + 0.08), linewidth = 0.35) +
  ggplot2::geom_text(data = nomogram_ticks, ggplot2::aes(Points, Row - 0.18, label = sprintf("%.1f", z_value)),
    size = 2.8, vjust = 1) +
  ggplot2::geom_segment(data = score_rows, ggplot2::aes(x = Points, xend = Points, y = 0.08, yend = -0.08),
    linewidth = 0.35) +
  ggplot2::geom_text(data = score_rows, ggplot2::aes(Points, 0.16, label = Label), size = 2.8, vjust = 0) +
  ggplot2::geom_segment(data = prob_rows, ggplot2::aes(x = Points, xend = Points, y = -0.92, yend = -1.08),
    linewidth = 0.35) +
  ggplot2::geom_text(data = prob_rows, ggplot2::aes(Points, -1.18, label = Label), size = 2.8, vjust = 1) +
  ggplot2::scale_x_continuous(limits = c(0, total_points_max), expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
  ggplot2::scale_y_continuous(breaks = c(-1, 0, 1, seq(2, length(top_genes) + 1L)),
    labels = c("Predicted LPD risk", "Total points", "Points", rev(top_genes)),
    limits = c(-1.45, length(top_genes) + 1.45)) +
  ggplot2::labs(title = "Five-gene ridge nomogram (GSE125498)",
    subtitle = "Predictors are discovery-cohort standardized expression values; apparent risk scale",
    x = NULL, y = NULL) + theme_project(10) +
  ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(colour = "#17324D"),
    panel.grid = ggplot2::element_blank())
save_plot(nomogram_plot, "T06_reduced_five_gene_nomogram", 10, 5.5)

calibration <- data.frame(Observed = y, Predicted = reduced_pred)
calibration$bin <- factor(pmin(5L, ceiling(rank(calibration$Predicted, ties.method = "first") /
  length(y) * 5L)), levels = 1:5, labels = paste0("Q", 1:5))
cal_summary <- do.call(rbind, lapply(split(calibration, calibration$bin), function(d) {
  ci <- stats::binom.test(sum(d$Observed), nrow(d))$conf.int
  data.frame(Bin = as.character(d$bin[[1]]), n = nrow(d), Observed = mean(d$Observed),
    Observed_CI_low = ci[[1]], Observed_CI_high = ci[[2]], Mean_predicted = mean(d$Predicted))
}))
write.csv(cal_summary, file.path(TABLES, "T06_discovery_calibration_bins.csv"), row.names = FALSE)
calibration_plot <- ggplot2::ggplot(cal_summary, ggplot2::aes(Mean_predicted, Observed)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey50") +
  ggplot2::geom_line(colour = "#0072B2", linewidth = 0.7) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Observed_CI_low, ymax = Observed_CI_high),
    width = 0.02, colour = "#0072B2") +
  ggplot2::geom_point(size = 2.5, colour = "#0072B2") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(title = "Apparent calibration of the five-gene ridge model",
    subtitle = "GSE125498; equal-count bins with exact binomial intervals",
    x = "Mean predicted probability", y = "Observed case fraction") + theme_project()
save_plot(calibration_plot, "T06_discovery_calibration", 6.5, 6)

thresholds <- seq(0.05, 0.95, by = 0.025)
prevalence <- mean(y == 1L)
dca <- do.call(rbind, lapply(thresholds, function(th) {
  classified <- reduced_pred >= th
  tp <- sum(classified & y == 1L); fp <- sum(classified & y == 0L)
  model_nb <- tp / length(y) - fp / length(y) * th / (1 - th)
  treat_all_nb <- prevalence - (1 - prevalence) * th / (1 - th)
  data.frame(Threshold = th, Strategy = c("Five-gene model", "Treat all", "Treat none"),
    Net_benefit = c(model_nb, treat_all_nb, 0))
}))
write.csv(dca, file.path(TABLES, "T06_discovery_decision_curve.csv"), row.names = FALSE)
dca_plot <- ggplot2::ggplot(dca, ggplot2::aes(Threshold, Net_benefit, colour = Strategy, linetype = Strategy)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::scale_colour_manual(values = c("Five-gene model" = "#0072B2", "Treat all" = "#D55E00",
    "Treat none" = "#333333")) +
  ggplot2::scale_linetype_manual(values = c("Five-gene model" = "solid", "Treat all" = "dashed",
    "Treat none" = "dotted")) +
  ggplot2::labs(title = "Apparent decision-curve analysis",
    subtitle = "Reduced five-gene ridge model; GSE125498",
    x = "Threshold probability", y = "Net benefit", colour = NULL, linetype = NULL) + theme_project()
save_plot(dca_plot, "T06_discovery_decision_curve", 7, 5)
dca_model <- dca[dca$Strategy == "Five-gene model", ]
dca_all <- dca[dca$Strategy == "Treat all", ]
dca_none <- dca[dca$Strategy == "Treat none", ]
dca_better <- dca_model$Net_benefit > dca_all$Net_benefit & dca_model$Net_benefit > dca_none$Net_benefit
write.csv(data.frame(
  Metric = c("n_control", "n_case", "Better-than-both thresholds", "Mean apparent net benefit over thresholds"),
  Value = c(sum(y == 0L), sum(y == 1L),
    if (any(dca_better)) paste(range(dca_model$Threshold[dca_better]), collapse = " to ") else "none",
    mean(dca_model$Net_benefit)),
  stringsAsFactors = FALSE
), file.path(TABLES, "T06_calibration_dca_summary.csv"), row.names = FALSE)

stability_plot <- ggplot2::ggplot(stability, ggplot2::aes(reorder(Gene, Selection_frequency), Selection_frequency)) +
  ggplot2::geom_col(fill = "#0072B2") + ggplot2::coord_flip() +
  ggplot2::labs(title = "Bootstrap LASSO selection stability", subtitle = paste(B, "bootstrap resamples; fixed full-data lambda.1se"),
                x = NULL, y = "Selection frequency") + theme_project()
save_plot(stability_plot, "T06_lasso_bootstrap_stability", 7, 5)

roc_sets <- list(GSE125498 = list(y = y, score = pred))
external_rows <- list()
for (cohort in c("GSE248762", "GSE130888")) {
  ext <- load_effluent(cohort)
  available <- intersect(genes, rownames(ext$expr))
  if (!length(available)) next
  # A missing hub is handled by refitting on discovery samples using the measurable subset.
  train_subset <- X0[, available, drop = FALSE]
  train_subset_z <- zscore(train_subset)
  mcv <- glmnet::cv.glmnet(train_subset_z$x, y, family = "binomial", alpha = 0,
                           nfolds = 5, type.measure = "deviance", standardize = FALSE)
  trained <- glmnet::glmnet(train_subset_z$x, y, family = "binomial", alpha = 0,
                            lambda = mcv$lambda.min, standardize = FALSE)
  eX <- t(ext$expr[available, , drop = FALSE])
  eXz <- zscore(eX)$x
  py <- as.integer(ext$meta$Group == "Case")
  pp <- as.numeric(stats::predict(trained, eXz, type = "response"))
  row <- auc_ci_row(paste0(cohort, " additional-cohort ridge"), py, pp)
  row$Genes_used <- paste(available, collapse = ";")
  row$Refit_subset_due_to_missing_genes <- length(available) < length(genes)
  row$Youden_cutoff <- {
    rp <- roc_points(py, pp)
    cut <- unique(pp)
    vals <- vapply(cut, function(th) {
      sens <- mean(pp[py == 1] >= th)
      spec <- mean(pp[py == 0] < th)
      sens + spec - 1
    }, numeric(1))
    cut[which.max(vals)]
  }
  row$Sensitivity_at_Youden <- mean(pp[py == 1] >= row$Youden_cutoff)
  row$Specificity_at_Youden <- mean(pp[py == 0] < row$Youden_cutoff)
  external_rows[[cohort]] <- row
  roc_sets[[cohort]] <- list(y = py, score = pp)
}
external <- if (length(external_rows)) do.call(rbind, external_rows) else data.frame()
if (nrow(external)) write.csv(external, file.path(TABLES, "T06_additional_cohort_model_performance.csv"), row.names = FALSE)

cv_rows <- list()
oof_loo <- rep(NA_real_, length(y))
for (i in seq_along(y)) {
  tr <- setdiff(seq_along(y), i)
  oof_loo[[i]] <- fit_tuned_ridge(X0[tr, , drop = FALSE], y[tr], X0[i, , drop = FALSE])$test
}
cv_rows[[1]] <- data.frame(Method = "LOO ridge", AUC = auc_rank(y, oof_loo))
repeat_oof <- matrix(NA_real_, nrow = length(y), ncol = 50L)
for (r in seq_len(ncol(repeat_oof))) {
  fold <- integer(length(y))
  for (cl in 0:1) {
    ii <- which(y == cl)
    fold[ii] <- sample(rep(seq_len(5), length.out = length(ii)))
  }
  for (k in 1:5) {
    te <- which(fold == k); tr <- setdiff(seq_along(y), te)
    repeat_oof[te, r] <- fit_tuned_ridge(X0[tr, , drop = FALSE], y[tr], X0[te, , drop = FALSE])$test
  }
}
cv_rows[[2]] <- data.frame(Method = "Repeated stratified 5-fold ridge (50 repeats)", AUC = mean(apply(repeat_oof, 2, function(z) auc_rank(y, z)), na.rm = TRUE))
write.csv(do.call(rbind, cv_rows), file.path(TABLES, "T06_discovery_cross_validation_AUC.csv"), row.names = FALSE)

roc_data <- do.call(rbind, lapply(names(roc_sets), function(nm) {
  d <- roc_points(roc_sets[[nm]]$y, roc_sets[[nm]]$score); d$Cohort <- nm; d
}))
pr <- ggplot2::ggplot(roc_data, ggplot2::aes(FPR, TPR, colour = Cohort)) +
  ggplot2::geom_step(linewidth = 0.8) + ggplot2::geom_abline(linetype = 2, colour = "grey50") +
  ggplot2::coord_equal() + ggplot2::labs(title = "Hub-model ROC curves", subtitle = "Discovery curve is apparent; additional cohorts are exploratory",
    x = "False-positive rate", y = "True-positive rate") + theme_project()
save_plot(pr, "T06_combined_ROC", 7, 6)

score_rows <- list()
for (b in bundles <- load_all_cohorts()) {
  present <- intersect(HUBS, rownames(b$expr))
  z <- scale(t(b$expr[present, , drop = FALSE]))
  score <- rowMeans(z, na.rm = TRUE)
  yy <- as.integer(b$meta$Group == "Case")
  ci <- bootstrap_auc(yy, score)
  score_rows[[length(score_rows) + 1L]] <- data.frame(Cohort = b$cohort, n_genes = length(present),
    genes_used = paste(present, collapse = ";"), n_case = sum(yy == 1), n_control = sum(yy == 0),
    AUC = auc_rank(yy, score), Bootstrap_low = ci[1], Bootstrap_high = ci[2])
  score_df <- data.frame(sample_id = b$meta$sample_id, Group = b$meta$Group, Hub_score = score)
  write.csv(score_df, file.path(TABLES, paste0("T06_matrisome_hub_score_", b$cohort, ".csv")), row.names = FALSE)
}
score_table <- do.call(rbind, score_rows)
write.csv(score_table, file.path(TABLES, "T06_matrisome_hub_score_AUC.csv"), row.names = FALSE)
score_plot_data <- do.call(rbind, lapply(bundles, function(b) {
  d <- read_csv_safe(file.path(TABLES, paste0("T06_matrisome_hub_score_", b$cohort, ".csv")))
  d$Cohort <- b$cohort; d
}))
ps <- ggplot2::ggplot(score_plot_data, ggplot2::aes(Group, Hub_score, fill = Group)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.12, size = 1.3) +
  ggplot2::facet_wrap(~Cohort, scales = "free_y") + ggplot2::scale_fill_manual(values = GROUP_COLORS) +
  ggplot2::labs(title = "Exploratory Matrisome Hub Score", x = NULL, y = "Mean within-cohort gene z-score") + theme_project()
save_plot(ps, "T06_matrisome_hub_score", 9, 5)

status <- c("TMM log2-CPM applied to external cohorts; missing-gene model refit on the measurable subset in discovery data.",
  paste("The 9-gene discovery model is apparent; bootstrap optimism estimated with", B, "resamples, conditional on the preselected panel."),
  "The consensus genes were selected on GSE125498; panel-selection optimism was not included in the bootstrap correction or fold-based estimates.",
  "A point-based visualization was generated from the fitted five-gene ridge model; it is not an rms nomogram.",
  "Calibration is apparent, grouped into five equal-count bins with exact binomial intervals; it is not bootstrap-corrected.",
  "Decision-curve net benefit is apparent and compared with treat-all and treat-none strategies; it is not externally validated.")
warning_messages <- warnings()
if (length(warning_messages)) {
  warning_lines <- capture.output(print(warning_messages))
  status <- c(status, paste("R warnings (first entries):", paste(head(warning_lines, 8), collapse = " ")))
}
task_log("T06_diagnostic_models", status, "partial")
dca_better_range <- if (any(dca_better)) paste(range(dca_model$Threshold[dca_better]), collapse = "–") else "none"
calibration_description <- paste(vapply(seq_len(nrow(cal_summary)), function(i) {
  sprintf("%s (n=%d): predicted %.2f, observed %.2f",
    cal_summary$Bin[[i]], cal_summary$n[[i]], cal_summary$Mean_predicted[[i]], cal_summary$Observed[[i]])
}, character(1)), collapse = "; ")
append_results_text("T6. Diagnostic modelling", paste(
  sprintf("The nine-gene ridge model had an apparent discovery AUC of %.3f (bootstrap 95%% CI %.3f–%.3f) and a model-fitting optimism-adjusted AUC of %.3f after 1,000 bootstrap resamples.",
    apparent$AUC, apparent$Bootstrap_low, apparent$Bootstrap_high, optimism_corrected),
  sprintf("Nested leave-one-out and repeated stratified five-fold AUCs were %.3f and %.3f, respectively; the top five bootstrap-stability genes were %s.",
    cv_rows[[1]]$AUC, cv_rows[[2]]$AUC, paste(top_genes, collapse = ", ")),
  sprintf("The reduced five-gene model had an apparent AUC of %.3f (bootstrap 95%% CI %.3f–%.3f); additional-cohort scores used within-cohort scaling and remained exploratory.",
    reduced_apparent$AUC, reduced_apparent$Bootstrap_low, reduced_apparent$Bootstrap_high),
  sprintf("A point-based nomogram diagram, apparent five-bin calibration plot with exact binomial intervals, and decision curve versus treat-all and treat-none strategies were generated for the reduced model (GSE125498: %d controls, %d long-term PD cases).",
    sum(y == 0L), sum(y == 1L)),
  paste0("The calibration-bin means were ", calibration_description, ". Apparent model net benefit exceeded both treat-all and treat-none strategies across thresholds ", dca_better_range, "."),
  "Calibration and decision-curve results were apparent, not bootstrap-corrected or externally validated. The diagram visualizes the ridge model but is not an rms nomogram; panel selection on GSE125498 means that model-performance estimates may retain selection optimism."))
