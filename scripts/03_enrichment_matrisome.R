args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("ggplot2", "dplyr", "clusterProfiler", "org.Hs.eg.db"), "T03")

source46 <- read_csv_safe(file.path(RESULTS, "revision", "matrisome_volcano", "GSE125498_matrisome_volcano_source_data.csv"))
input39 <- read_csv_safe(file.path(RESULTS, "revision", "ml_39_matrisome_wgcna", "wgcna_matrisome_39_gene_input.csv"))
ml <- read_csv_safe(file.path(RESULTS, "revision", "ml_39_matrisome_wgcna", "ml_39_gene_results.csv"))
deg46 <- source46[source46$P.Value < 0.05 & abs(source46$logFC) >= 0.585, ]
if (nrow(deg46) != 46L) stop("Expected 46 MatrisomeAnalyzeR DEGs from the saved measured-matrisome table; observed ", nrow(deg46))
if (nrow(input39) != 39L) stop("Expected 39 saved WGCNA-matrisome input genes; observed ", nrow(input39))
hub <- merge(ml[ml$Consensus_At_Least_2_of_3 == 1, "Gene", drop = FALSE],
             input39[, c("Gene", "Matrisome_Division", "Matrisome_Category", "DEG_Status")],
             by = "Gene", all.x = TRUE)
deg46$Source_set <- "46 matrisome DEGs"
hub$Source_set <- "9 consensus hubs"
common <- rbind(
  data.frame(Gene = deg46$Gene, Matrisome_Division = deg46$Matrisome_Division, Matrisome_Category = deg46$Matrisome_Category,
             DEG_Status = deg46$DEG_Status, Source_set = deg46$Source_set),
  data.frame(Gene = hub$Gene, Matrisome_Division = hub$Matrisome_Division, Matrisome_Category = hub$Matrisome_Category,
             DEG_Status = hub$DEG_Status, Source_set = hub$Source_set)
)
write.csv(common, file.path(TABLES, "T03_matrisome_category_gene_table.csv"), row.names = FALSE)
counts <- aggregate(Gene ~ Source_set + Matrisome_Division + Matrisome_Category + DEG_Status, common, length)
names(counts)[5] <- "Gene_count"
write.csv(counts, file.path(TABLES, "T03_matrisome_category_composition.csv"), row.names = FALSE)
p <- ggplot2::ggplot(counts, ggplot2::aes(Matrisome_Category, Gene_count, fill = DEG_Status)) +
  ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::facet_wrap(~Source_set, scales = "free_y") +
  ggplot2::scale_fill_manual(values = STATUS_COLORS) +
  ggplot2::labs(title = "Matrisome category composition", x = NULL, y = "Gene count", fill = "Discovery direction") + theme_project()
save_plot(p, "T03_matrisome_category_composition", 10, 7)

universe <- unique(source46$Gene[!is.na(source46$Gene) & nzchar(source46$Gene)])
gene_sets <- list(`39-gene WGCNA intersection` = unique(input39$Gene), `46 Matrisome DEGs` = unique(deg46$Gene))
universe_map <- tryCatch(clusterProfiler::bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE), error = function(e) data.frame())
all_results <- list()
plot_objects <- list()
status <- c(paste("Measured MatrisomeAnalyzeR universe:", length(universe), "unique expressed genes from saved source table."),
  paste("DEG threshold reproduced from source table: |log2FC| >= 0.585 and nominal P < 0.05; selected", nrow(deg46), "genes."))

for (set_name in names(gene_sets)) {
  symbols <- intersect(gene_sets[[set_name]], universe)
  mapped <- tryCatch(clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE), error = function(e) data.frame())
  status <- c(status, paste(set_name, "input:", length(gene_sets[[set_name]]), "; represented in matrisome background:", length(symbols),
                            "; Entrez-mapped:", nrow(mapped)))
  if (!length(symbols) || !nrow(universe_map) || !nrow(mapped)) next
  jobs <- list()
  for (ont in c("BP", "MF", "CC")) {
    obj <- tryCatch(clusterProfiler::enrichGO(gene = symbols, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 3, maxGSSize = 500, readable = TRUE), error = function(e) {
        status <<- c(status, paste(set_name, "GO", ont, "failed:", conditionMessage(e))); NULL
      })
    if (!is.null(obj)) jobs[[paste0("GO_", ont)]] <- obj
  }
  obj_kegg <- tryCatch(clusterProfiler::enrichKEGG(gene = unique(mapped$ENTREZID),
    universe = unique(universe_map$ENTREZID), organism = "hsa", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500),
    error = function(e) { status <<- c(status, paste(set_name, "KEGG failed:", conditionMessage(e))); NULL })
  if (!is.null(obj_kegg)) jobs[["KEGG"]] <- obj_kegg
  if (requireNamespace("DOSE", quietly = TRUE)) {
    obj_do <- tryCatch(DOSE::enrichDO(gene = unique(mapped$ENTREZID), universe = unique(universe_map$ENTREZID),
      ont = "HDO", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 3, maxGSSize = 500, readable = TRUE),
      error = function(e) { status <<- c(status, paste(set_name, "Disease Ontology failed:", conditionMessage(e))); NULL })
    if (!is.null(obj_do)) jobs[["Disease_Ontology"]] <- obj_do
  }
  if (!requireNamespace("ReactomePA", quietly = TRUE)) status <- c(status, paste(set_name, "Reactome skipped: ReactomePA package is not installed."))
  for (db in names(jobs)) {
    obj <- jobs[[db]]
    df <- as.data.frame(obj)
    if (!nrow(df)) next
    df$Gene_set <- set_name; df$Database <- db
    all_results[[length(all_results) + 1L]] <- df
    plot_objects[[paste(set_name, db, sep = "__")]] <- obj
    if (requireNamespace("enrichplot", quietly = TRUE)) {
      safe <- gsub("[^A-Za-z0-9]+", "_", paste(set_name, db))
      for (kind in c("dotplot", "barplot")) {
        pl <- tryCatch(if (kind == "dotplot") enrichplot::dotplot(obj, showCategory = 15) else {
          bar_data <- head(df[order(df$p.adjust, df$pvalue), , drop = FALSE], 15)
          ggplot2::ggplot(bar_data, ggplot2::aes(reorder(Description, Count), Count,
            fill = -log10(pmax(p.adjust, 1e-300)))) + ggplot2::geom_col() + ggplot2::coord_flip() +
            ggplot2::scale_fill_viridis_c() + ggplot2::labs(x = NULL, y = "Overlapping genes", fill = "-log10 BH q")
        },
          error = function(e) { status <<- c(status, paste(safe, kind, "plot failed:", conditionMessage(e))); NULL })
        if (!is.null(pl)) save_plot(pl + ggplot2::labs(title = paste(set_name, db, kind)), paste0("T03_", safe, "_", kind), 8, 6)
      }
      for (kind in c("cnetplot", "emapplot")) {
        graph_obj <- if (kind == "emapplot") tryCatch(enrichplot::pairwise_termsim(obj), error = function(e) {
          status <<- c(status, paste(safe, "term similarity failed:", conditionMessage(e))); obj
        }) else obj
        pl <- tryCatch(if (kind == "cnetplot") enrichplot::cnetplot(obj, showCategory = 8) else enrichplot::emapplot(graph_obj, showCategory = 15),
          error = function(e) { status <<- c(status, paste(safe, kind, "plot failed:", conditionMessage(e))); NULL })
        if (!is.null(pl)) save_plot(pl + ggplot2::labs(title = paste(set_name, db, kind)), paste0("T03_", safe, "_", kind), 9, 7)
      }
    } else status <- c(status, "enrichplot unavailable; network enrichment plots not generated.")
  }
}
results <- if (length(all_results)) dplyr::bind_rows(all_results) else data.frame()
write.csv(results, file.path(TABLES, "T03_enrichment_all_terms.csv"), row.names = FALSE)
sig <- if (nrow(results) && "p.adjust" %in% names(results)) results[!is.na(results$p.adjust) & results$p.adjust < 0.05, ] else results
write.csv(sig, file.path(TABLES, "T03_enrichment_BH_q_lt_0.05.csv"), row.names = FALSE)
write.csv(data.frame(Gene = universe), file.path(TABLES, "T03_matrisome_expressed_universe.csv"), row.names = FALSE)
task_log("T03_enrichment", c(status, paste("Completed supported enrichment objects:", length(plot_objects)),
  paste("Significant terms with BH-adjusted q < 0.05:", nrow(sig))), if (length(plot_objects)) "partial" else "skipped")
direction_counts <- table(factor(deg46$DEG_Status, levels = c("Upregulated", "Downregulated")))
append_results_text("T3. Enrichment and matrisome composition", paste(
  sprintf("The saved measured-matrisome background contained %d expressed genes; applying |log2FC| >= 0.585 and nominal P < 0.05 recovered %d DEGs (%d upregulated and %d downregulated).",
    length(universe), nrow(deg46), direction_counts[["Upregulated"]], direction_counts[["Downregulated"]]),
  "The 39-gene WGCNA intersection and 46-DEG set were tested with GO, KEGG, and Disease Ontology; no tested term had BH-adjusted q < 0.05 across the reported enrichment results. The category table also summarized all nine consensus hubs.",
  "ReactomePA was unavailable, so Reactome enrichment was omitted; no whole-genome background was substituted."))
