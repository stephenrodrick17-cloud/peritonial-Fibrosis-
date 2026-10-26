args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("limma", "clusterProfiler", "org.Hs.eg.db", "ggplot2", "jsonlite", "dplyr"), "T07")

b <- load_discovery()
expr <- b$expr
hallmark_path <- file.path(RESULTS, "revision", "panel_v2_context", "hallmark_gene_sets.json")
if (!file.exists(hallmark_path)) {
  write_skip("T07_single_gene_gsea", "saved Hallmark gene-set JSON was missing.")
} else {
  hm <- jsonlite::fromJSON(hallmark_path, simplifyVector = FALSE)[["gene_sets"]]
  themes <- c("Epithelial Mesenchymal Transition", "TGF-beta Signaling", "Inflammatory Response",
              "IL-6/JAK/STAT3 Signaling", "TNF-alpha Signaling via NF-kB", "Angiogenesis", "Hypoxia")
  hm <- hm[intersect(themes, names(hm))]
  t2g <- do.call(rbind, lapply(names(hm), function(term) data.frame(term = term, gene = unlist(hm[[term]], use.names = FALSE))))
  all_gsea <- list()
  status <- character()
  for (hub in HUBS) {
    if (!hub %in% rownames(expr)) {
      status <- c(status, paste(hub, "not present in discovery matrix"))
      next
    }
    high <- expr[hub, ] > stats::median(expr[hub, ], na.rm = TRUE)
    if (length(unique(high)) < 2) {
      status <- c(status, paste(hub, "median split produced one group; skipped"))
      next
    }
    design <- stats::model.matrix(~ factor(high, levels = c(FALSE, TRUE)))
    fit <- limma::eBayes(limma::lmFit(expr, design))
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    tt$Gene <- rownames(tt)
    tt <- tt[!duplicated(tt$Gene), ]
    ranks <- tt$t
    names(ranks) <- tt$Gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
    hm_res <- tryCatch(clusterProfiler::GSEA(ranks, TERM2GENE = t2g, pvalueCutoff = 1,
      pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, verbose = FALSE),
      error = function(e) { status <<- c(status, paste(hub, "Hallmark GSEA failed:", conditionMessage(e))); NULL })
    if (!is.null(hm_res) && nrow(as.data.frame(hm_res))) {
      df <- as.data.frame(hm_res); df$Hub <- hub; df$Collection <- "Saved Hallmark subset"; all_gsea[[length(all_gsea) + 1L]] <- df
    }
    status <- c(status, paste(hub, "ranked all expressed discovery genes by limma t-statistic for above- vs below-median hub expression."))
  }
  status <- c(status, "GO and KEGG GSEA were omitted: GO preranked analysis exceeded the execution budget and KEGG requires repeated online annotation retrieval. The cached, focused Hallmark subset was analyzed instead.")
  gsea <- if (length(all_gsea)) do.call(rbind, all_gsea) else data.frame()
  write.csv(gsea, file.path(TABLES, "T07_single_hub_GSEA.csv"), row.names = FALSE)
  if (nrow(gsea)) {
    top <- gsea[gsea$Collection == "Saved Hallmark subset", , drop = FALSE]
    top <- do.call(rbind, lapply(split(top, top$Hub), function(d) head(d[order(d$p.adjust), ], 5)))
    write.csv(top, file.path(TABLES, "T07_top_Hallmark_terms_per_hub.csv"), row.names = FALSE)
    if (nrow(top)) {
      pg <- ggplot2::ggplot(top, ggplot2::aes(NES, reorder(Description, NES), colour = p.adjust, size = setSize)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~Hub, scales = "free_y") +
        ggplot2::scale_colour_viridis_c(direction = -1) +
        ggplot2::labs(title = "Top saved Hallmark enrichments per hub", x = "Normalized enrichment score", y = NULL, colour = "BH FDR") +
        theme_project(8)
      save_plot(pg, "T07_hallmark_GSEA_top_terms", 13, 10)
    }
  }
  # The cached MyGeneSet collection is a focused subset, not the full MSigDB Hallmark release.
  sig_names <- names(hm)
  if (length(sig_names)) {
    signature_scores <- list()
    for (cohort in load_all_cohorts()) {
      scores <- sapply(hm, function(gs) {
        common <- intersect(unlist(gs, use.names = FALSE), rownames(cohort$expr))
        if (length(common) < 5) return(rep(NA_real_, ncol(cohort$expr)))
        rowMeans(scale(t(cohort$expr[common, , drop = FALSE])), na.rm = TRUE)
      })
      if (is.null(dim(scores))) scores <- matrix(scores, ncol = 1)
      rownames(scores) <- colnames(cohort$expr)
      hubs <- intersect(HUBS, rownames(cohort$expr))
      for (hub in hubs) for (sig in colnames(scores)) {
        valid <- is.finite(scores[, sig])
        if (sum(valid) < 5) next
        test <- suppressWarnings(stats::cor.test(cohort$expr[hub, valid], scores[valid, sig], method = "spearman", exact = FALSE))
        signature_scores[[length(signature_scores) + 1L]] <- data.frame(Cohort = cohort$cohort, Hub = hub, Signature = sig,
          Spearman_rho = unname(test$estimate), P = test$p.value, n = sum(valid))
      }
    }
    corr <- if (length(signature_scores)) do.call(rbind, signature_scores) else data.frame()
    if (nrow(corr)) {
      corr$BH_FDR_within_cohort <- ave(corr$P, corr$Cohort, FUN = function(p) p.adjust(p, method = "BH"))
      write.csv(corr, file.path(TABLES, "T07_hub_hallmark_signature_correlations.csv"), row.names = FALSE)
      ph <- ggplot2::ggplot(corr, ggplot2::aes(Signature, Hub, fill = Spearman_rho)) +
        ggplot2::geom_tile() + ggplot2::facet_wrap(~Cohort) +
        ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", limits = c(-1, 1)) +
        ggplot2::labs(title = "Hub expression vs saved Hallmark signature scores", x = NULL, y = NULL, fill = "Spearman rho") +
        theme_project(8) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      save_plot(ph, "T07_hub_hallmark_signature_correlations", 12, 7)
    }
  }
  task_log("T07_single_gene_gsea", c(status,
    paste("Used cached, project-local Hallmark gene sets for these themes:", paste(sig_names, collapse = "; ")),
    "High/low is an exploratory median split in the same discovery cohort; results are not independent validation.",
    "A fresh msigdbr download was not performed; msigdbr is not installed. Signature scores are mean within-cohort z-scores, not GSVA/ssGSEA."), "partial")
  sig_hallmark <- gsea[gsea$Collection == "Saved Hallmark subset" & !is.na(gsea$p.adjust) & gsea$p.adjust < 0.05, , drop = FALSE]
  top_signals <- unique(sig_hallmark$Description)
  append_results_text("T7. Single-hub pathway analysis", paste(
    sprintf("Each measurable hub defined above- and below-median groups in GSE125498, and all expressed genes were ranked by limma t-statistics. The focused cached Hallmark analysis identified %d hub–pathway associations at BH-adjusted P < 0.05; leading themes included %s.",
      nrow(sig_hallmark), if (length(top_signals)) paste(head(top_signals, 4), collapse = ", ") else "none"),
    "Hallmark scores were calculated as mean within-cohort gene z-scores; GO preranked analysis exceeded the execution budget and KEGG was omitted because it requires repeated online annotation retrieval.",
    "These analyses used the discovery cohort and are post hoc; they do not constitute external pathway validation."))
}
