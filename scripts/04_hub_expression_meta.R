args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("ggplot2", "dplyr", "tidyr"), "T04")

bundles <- load_all_cohorts()
long <- do.call(rbind, lapply(bundles, hub_expr_long))
if (!nrow(long)) stop("No consensus hub genes were measurable in loaded expression matrices.")
long$Cohort <- factor(long$cohort, levels = c("GSE125498", "GSE248762", "GSE130888"))
tests <- do.call(rbind, lapply(split(long, interaction(long$Cohort, long$Gene, drop = TRUE)), function(d) {
  p <- tryCatch(stats::wilcox.test(Expression ~ Group, d, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.frame(Cohort = as.character(d$Cohort[[1]]), Gene = d$Gene[[1]], n_control = sum(d$Group == "Control"),
             n_case = sum(d$Group == "Case"), median_control = median(d$Expression[d$Group == "Control"]),
             median_case = median(d$Expression[d$Group == "Case"]), wilcoxon_P = p)
}))
tests$BH_FDR_9_genes_within_cohort <- ave(tests$wilcoxon_P, tests$Cohort, FUN = function(p) p.adjust(p, method = "BH"))
write.csv(tests, file.path(TABLES, "T04_hub_expression_wilcoxon.csv"), row.names = FALSE)
p <- ggplot2::ggplot(long, ggplot2::aes(Group, Expression, fill = Group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.62) +
  ggplot2::geom_jitter(width = 0.12, size = 1.3, alpha = 0.8) +
  ggplot2::facet_grid(Gene ~ Cohort, scales = "free_y") +
  ggplot2::scale_fill_manual(values = GROUP_COLORS) +
  ggplot2::labs(title = "Consensus hub expression across cohorts", x = NULL, y = "Normalized expression (within-cohort scale)", fill = "Group") +
  theme_project(9)
save_plot(p, "T04_hub_expression_all_cohorts", 12, 14)

effect_rows <- list()
disc <- read_csv_safe(file.path(TABLES, "GSE125498_gene_level_toptable.csv"))
for (gene in HUBS) {
  r <- disc[disc$Gene == gene, , drop = FALSE]
  if (nrow(r)) {
    se <- if (is.finite(r$t[[1]]) && r$t[[1]] != 0) abs(r$logFC[[1]] / r$t[[1]]) else NA_real_
    effect_rows[[length(effect_rows) + 1L]] <- data.frame(Gene = gene, Cohort = "GSE125498", log2FC = r$logFC[[1]],
      SE = se, P_raw = r$P.Value[[1]], FDR_original = r$adj.P.Val[[1]], Detected = TRUE)
  }
}
for (cohort in c("GSE248762", "GSE130888")) {
  f <- read_csv_safe(file.path(RESULTS, "revision", "effluent_validation", cohort, "analysis_20261009", "edger_all_gene_results.csv"))
  d <- file.path(RESULTS, "revision", "effluent_validation", cohort, "analysis_20261009")
  panel_detection <- read_csv_safe(file.path(d, "primary_per_gene_results.csv"))
  detected <- setNames(as.logical(panel_detection$detected), panel_detection$gene)
  f <- f[f$contrast == "LV_vs_SV" & f$gene_symbol %in% HUBS & vapply(f$gene_symbol, function(g) isTRUE(detected[g]), logical(1)), , drop = FALSE]
  for (i in seq_len(nrow(f))) {
    se <- (f$ci_upper[[i]] - f$ci_lower[[i]]) / (2 * 1.96)
    effect_rows[[length(effect_rows) + 1L]] <- data.frame(Gene = f$gene_symbol[[i]], Cohort = cohort, log2FC = f$logFC[[i]],
      SE = se, P_raw = f$PValue[[i]], FDR_original = f$FDR_all_genes[[i]], Detected = TRUE)
  }
}
effects <- do.call(rbind, effect_rows)
effects$FDR_BH_across_hub_tests <- ave(effects$P_raw, effects$Cohort, FUN = function(p) p.adjust(p, method = "BH"))
grid <- expand.grid(Gene = HUBS, Cohort = c("GSE125498", "GSE248762", "GSE130888"), stringsAsFactors = FALSE)
effects <- merge(grid, effects, by = c("Gene", "Cohort"), all.x = TRUE, sort = FALSE)
effects$Detected[is.na(effects$Detected)] <- FALSE
write.csv(effects, file.path(TABLES, "T04_hub_effects_by_cohort.csv"), row.names = FALSE)
heat <- effects
heat$label <- ifelse(heat$FDR_BH_across_hub_tests < 0.05, "*", "")
pheat <- ggplot2::ggplot(heat, ggplot2::aes(Cohort, Gene, fill = log2FC)) +
  ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = paste0(sprintf("%.2f", log2FC), label)), size = 3) +
  ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, na.value = "grey80") +
  ggplot2::labs(title = "Consensus-hub effect directions across cohorts", subtitle = "* BH-adjusted FDR < 0.05 within cohort",
                x = NULL, y = NULL, fill = "log2 fold change") + theme_project()
save_plot(pheat, "T04_cross_cohort_log2FC_heatmap", 8, 5.5)

meta_rows <- list()
for (gene in HUBS) {
  d <- effects[effects$Gene == gene & is.finite(effects$log2FC) & is.finite(effects$SE) & effects$SE > 0, ]
  if (nrow(d) < 2) next
  w <- 1 / d$SE^2
  fixed <- sum(w * d$log2FC) / sum(w)
  Q <- sum(w * (d$log2FC - fixed)^2)
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (nrow(d) - 1)) / C)
  wr <- 1 / (d$SE^2 + tau2)
  est <- sum(wr * d$log2FC) / sum(wr)
  se <- sqrt(1 / sum(wr))
  meta_rows[[length(meta_rows) + 1L]] <- data.frame(Gene = gene, k = nrow(d), log2FC_random = est,
    CI_low = est - 1.96 * se, CI_high = est + 1.96 * se, tau2_DL = tau2,
    I2_percent = if (Q > 0) max(0, (Q - (nrow(d) - 1)) / Q) * 100 else 0,
    Q = Q, Q_P = stats::pchisq(Q, df = nrow(d) - 1, lower.tail = FALSE),
    Direction_consistent_across_available_cohorts = length(unique(sign(d$log2FC))) == 1L)
}
meta <- if (length(meta_rows)) do.call(rbind, meta_rows) else data.frame()
write.csv(meta, file.path(TABLES, "T04_random_effects_meta_analysis.csv"), row.names = FALSE)
if (nrow(meta)) {
  pf <- ggplot2::ggplot(meta, ggplot2::aes(log2FC_random, reorder(Gene, log2FC_random))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = 0.18) +
    ggplot2::geom_point(size = 2.3, colour = "#0072B2") +
    ggplot2::labs(title = "Exploratory random-effects meta-analysis", subtitle = "DerSimonian–Laird; cohorts may not be independent",
                  x = "Pooled log2 fold change (95% CI)", y = NULL) + theme_project()
  save_plot(pf, "T04_random_effects_forest", 8, 5)
}
task_log("T04_hub_expression_meta", c(paste("Wilcoxon group comparisons:", nrow(tests), "; BH adjustment applied within each cohort across measured hubs."),
  "Effect-size tables include raw P and BH adjustment across the hub tests within cohort; meta-analysis used inverse-variance DerSimonian-Laird random effects.",
  "Standard errors came from discovery t-statistics and external edgeR confidence intervals; undetected genes were excluded from pooling.",
  "External cohorts were described as additional; independence is unverified. The meta-analysis is exploratory."))
append_results_text("T4. Hub expression and cross-cohort effects", paste("Hub expression was summarized in all three cohorts, and group comparisons used Wilcoxon tests with Benjamini–Hochberg adjustment across the measured panel within each cohort.",
  "Effect directions and adjusted FDRs were displayed by cohort, with an exploratory DerSimonian–Laird random-effects summary across available estimates.",
  "The small cohorts and unverified overlap between additional cohorts limit inference from cross-cohort pooling."))
