#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve this script's path.")
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(out_dir, "..", "..", ".."), winslash = "/")
result_path <- file.path(out_dir, "hallmark_gsea_all.csv")
ranking_path <- file.path(out_dir, "hallmark_ranked_genes.csv")
gene_sets_path <- file.path(root, "results", "revision", "panel_v2_context", "hallmark_gene_sets.json")
for (path in c(result_path, ranking_path, gene_sets_path)) {
  if (!file.exists(path)) stop("Required figure input missing: ", path)
}

res <- fread(result_path)
ranked <- fread(ranking_path)
stats <- setNames(ranked$t, ranked$Gene)
stats <- sort(stats, decreasing = TRUE)
cached <- jsonlite::fromJSON(gene_sets_path, simplifyVector = FALSE)
pathways <- lapply(cached$gene_sets, function(genes) intersect(unique(toupper(unlist(genes, use.names = FALSE))), names(stats)))
pathways <- pathways[lengths(pathways) >= 15L & lengths(pathways) <= 500L]
sig <- res[!is.na(FDR) & FDR < 0.05]

plot_nes <- function() {
  if (!nrow(sig)) {
    plot.new()
    title("No Hallmark sets passed FDR < 0.05")
    return(invisible(NULL))
  }
  ordered <- copy(sig[order(NES)])
  colors <- ifelse(ordered$NES > 0, "#C84B31", "#2878A5")
  par(mar = c(5, 18, 3, 2))
  plot(
    ordered$NES, seq_len(nrow(ordered)), type = "n", yaxt = "n",
    xlab = "Normalized enrichment score (NES)", ylab = "",
    main = "Significant MSigDB Hallmark sets; FDR < 0.05"
  )
  abline(v = 0, col = "grey60", lty = 2)
  segments(0, seq_len(nrow(ordered)), ordered$NES, seq_len(nrow(ordered)),
           lwd = 5, col = colors)
  axis(2, at = seq_len(nrow(ordered)), labels = ordered$pathway, las = 2, cex.axis = 0.72)
  legend("bottomright", legend = c("Positive NES (LPD)", "Negative NES (SPD)"),
         col = c("#C84B31", "#2878A5"), lwd = 4, bty = "n", cex = 0.8)
}

save_base_plot <- function(stem, draw, width = 10, height = 8) {
  grDevices::png(file.path(out_dir, paste0(stem, ".png")),
                 width = width, height = height, units = "in", res = 300)
  draw()
  grDevices::dev.off()
  grDevices::pdf(file.path(out_dir, paste0(stem, ".pdf")),
                 width = width, height = height, onefile = TRUE)
  draw()
  grDevices::dev.off()
}

save_base_plot("hallmark_significant_NES", plot_nes, width = 12, height = max(6, 0.28 * nrow(sig) + 2))

save_enrichment <- function(direction, stem) {
  candidates <- sig[if (direction == "up") NES > 0 else NES < 0]
  if (!nrow(candidates)) {
    writeLines(
      paste("No significant", direction, "Hallmark set was available for an enrichment plot."),
      file.path(out_dir, paste0(stem, "_not_available.txt"))
    )
    return(invisible(NULL))
  }
  set <- candidates[which.max(abs(NES)), pathway]
  plot_fun <- function() {
    fgsea::plotEnrichment(pathways[[set]], stats) +
      ggplot2::labs(
        title = set,
        subtitle = if (direction == "up") "Positive NES; enrichment toward LPD" else "Negative NES; enrichment toward SPD",
        x = "Rank in moderated-t gene list", y = "Running enrichment score"
      ) +
      ggplot2::theme_classic(base_size = 11)
  }
  plot <- plot_fun()
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = 8, height = 6, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = 8, height = 6, device = grDevices::cairo_pdf, bg = "white")
  source <- data.table(direction = direction, pathway = set, NES = candidates[pathway == set, NES],
                       FDR = candidates[pathway == set, FDR])
  fwrite(source, file.path(out_dir, paste0(stem, "_source.csv")))
}

save_enrichment("up", "hallmark_top_up_enrichment")
save_enrichment("down", "hallmark_top_down_enrichment")
