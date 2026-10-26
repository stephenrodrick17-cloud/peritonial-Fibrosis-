source(file.path("results", "paper_style", "00_helpers.R"))

limma <- read_project_csv(file.path(TABLES, "GSE125498_all_results.csv"))
expr_df <- read_project_csv(file.path(TABLES, "GSE125498_full_expression_matrix.csv"))
meta <- read_project_csv(file.path(TABLES, "GSE125498_sample_metadata.csv"))
stopifnot(all(c("Gene.symbol", "logFC", "t", "P.Value", "adj.P.Val") %in% names(limma)))

limma <- limma[!is.na(limma$Gene.symbol) & nzchar(limma$Gene.symbol), ]
limma <- limma[order(-abs(limma$t), limma$Gene.symbol), ]
limma <- limma[!duplicated(limma$Gene.symbol), ]
limma$minus_log10_FDR <- -log10(pmax(limma$adj.P.Val, .Machine$double.xmin))
limma$Direction <- ifelse(limma$logFC > 0, "Higher in LPD", "Higher in SPD")
limma$Call <- ifelse(limma$adj.P.Val < 0.05 & abs(limma$logFC) >= 0.585,
                     "BH FDR < 0.05 and |log2FC| >= 0.585", "Other")
save_csv(limma, "Fig2_limma_gene_results_unique_symbols.csv")

volcano <- ggplot(limma, aes(logFC, minus_log10_FDR, colour = Call)) +
  geom_point(size = 0.8, alpha = 0.65) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = 2, colour = "grey45") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey45") +
  scale_colour_manual(values = c("Other" = "#B7C1CB",
    "BH FDR < 0.05 and |log2FC| >= 0.585" = "#D95F02")) +
  labs(title = "Differential expression: LPD versus SPD",
       subtitle = "BH-adjusted P < 0.05; absolute log2 fold change >= 0.585",
       x = "log2 fold change (LPD relative to SPD)", y = "-log10(BH-adjusted P)",
       colour = NULL) + theme_paper()
save_plot_pair(add_letter(volcano, "B"), "Fig2B_volcano", 7, 5)

expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- expr_df$Gene
if (!all(meta$sample_id %in% colnames(expr))) stop("Metadata samples do not match expression columns.")
expr <- expr[, meta$sample_id, drop = FALSE]
top <- head(limma[order(limma$adj.P.Val, -abs(limma$logFC)), "Gene.symbol"], 50)
top <- intersect(top, rownames(expr))
z <- t(scale(t(expr[top, , drop = FALSE])))
heat_long <- as.data.frame(as.table(z), stringsAsFactors = FALSE)
names(heat_long) <- c("Gene", "sample_id", "Row_Z")
heat_long$Gene <- factor(heat_long$Gene, levels = rev(top))
heat_long$sample_id <- factor(heat_long$sample_id, levels = meta$sample_id)
heat_long <- merge(heat_long, meta[, c("sample_id", "group")], by = "sample_id", sort = FALSE)
save_csv(heat_long, "Fig2C_top50_heatmap_source.csv")
save_csv(data.frame(Gene = top, rank = seq_along(top)), "Fig2C_top50_genes.csv")
heat <- ggplot(heat_long, aes(sample_id, Gene, fill = Row_Z)) +
  geom_tile() + facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Gene-wise\nz-score") +
  labs(title = "Top 50 genes by BH-adjusted P value", x = NULL, y = NULL) +
  theme_paper(8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                         panel.grid = element_blank())
save_plot_pair(add_letter(heat, "C"), "Fig2C_top50_heatmap", 11, 9)

gsea_status <- "not_run"
gsea_detail <- "clusterProfiler GSEA was not attempted."
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  mapped <- tryCatch(
    clusterProfiler::bitr(limma$Gene.symbol, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE),
    error = function(e) e
  )
  if (inherits(mapped, "error")) {
    gsea_status <- "failed"
    gsea_detail <- paste("SYMBOL-to-Entrez mapping failed:", conditionMessage(mapped))
  } else {
    ranks <- limma[, c("Gene.symbol", "t")]
    names(ranks) <- c("SYMBOL", "score")
    ranks <- merge(ranks, mapped, by = "SYMBOL")
    ranks <- ranks[order(-abs(ranks$score)), ]
    ranks <- ranks[!duplicated(ranks$ENTREZID), ]
    gene_list <- setNames(ranks$score, ranks$ENTREZID)
    gene_list <- sort(gene_list, decreasing = TRUE)
    save_csv(data.frame(ENTREZID = names(gene_list), moderated_t = unname(gene_list)),
             "Fig2_DEG_moderated_t_ranking_Entrez.csv")
    gsea <- tryCatch(
      clusterProfiler::gseKEGG(geneList = gene_list, organism = "hsa",
        minGSSize = 15, maxGSSize = 500, pAdjustMethod = "BH",
        pvalueCutoff = 1, verbose = FALSE),
      error = function(e) e
    )
    if (inherits(gsea, "error")) {
      gsea_status <- "failed"
      gsea_detail <- paste("KEGG preranked GSEA failed:", conditionMessage(gsea))
    } else {
      gtab <- as.data.frame(gsea)
      save_csv(gtab, "Fig2_DEG_preranked_KEGG_all.csv")
      sig <- if (nrow(gtab) && "p.adjust" %in% names(gtab)) {
        gtab[!is.na(gtab$p.adjust) & gtab$p.adjust < 0.05, , drop = FALSE]
      } else gtab
      up <- sig[sig$NES > 0, , drop = FALSE]
      down <- sig[sig$NES < 0, , drop = FALSE]
      save_csv(up, "Fig2D_KEGG_up_LPD_FDR05.csv")
      save_csv(down, "Fig2E_KEGG_up_SPD_FDR05.csv")
      gsea_status <- "complete"
      gsea_detail <- paste(nrow(gtab), "KEGG sets tested;", nrow(up),
                           "FDR-significant positive-NES sets and", nrow(down),
                           "FDR-significant negative-NES sets.")
      if (nrow(up)) {
        p_up <- enrichplot::ridgeplot(gsea, showCategory = up$Description) +
          labs(title = "Up-regulated KEGG pathways (positive NES = higher in LPD)")
        save_plot_pair(add_letter(p_up, "D"), "Fig2D_KEGG_ridge_up", 9, 6)
      }
      if (nrow(down)) {
        p_down <- enrichplot::ridgeplot(gsea, showCategory = down$Description) +
          labs(title = "Down-regulated KEGG pathways (negative NES = higher in SPD)")
        save_plot_pair(add_letter(p_down, "E"), "Fig2E_KEGG_ridge_down", 9, 6)
      }
      if (!nrow(up) || !nrow(down)) {
        gsea_detail <- paste(gsea_detail, "No ridge plot was emitted for any direction with zero significant sets.")
      }
    }
  }
}
save_csv(data.frame(status = gsea_status, detail = gsea_detail), "Fig2_GSEA_run_status.csv")
write_status("Fig2", if (gsea_status == "failed") "partial" else "complete",
             paste("Single discovery cohort; PCA/batch correction not applicable.", gsea_detail))
