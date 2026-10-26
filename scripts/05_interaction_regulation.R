args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("ggplot2", "igraph", "dplyr", "tidyr"), "T05")

string_dir <- file.path(RESULTS, "revision", "matrisome_string_ppi", "functional")
edge_file <- file.path(string_dir, "string_functional_interactions.tsv")
if (!file.exists(edge_file)) {
  write_skip("T05_interaction", "saved STRING functional edge table was missing; no new network query was attempted.")
} else {
  edges <- read.delim(edge_file, check.names = FALSE)
  edges <- edges[edges$score >= 0.4, , drop = FALSE]
  candidate <- unique(c(edges$input_gene_A[edges$input_gene_A %in% HUBS],
                        edges$input_gene_B[edges$input_gene_B %in% HUBS]))
  incident <- edges[edges$input_gene_A %in% HUBS | edges$input_gene_B %in% HUBS, , drop = FALSE]
  neighbours <- sort(unique(c(incident$input_gene_A, incident$input_gene_B)))
  nodes <- unique(c(HUBS[HUBS %in% neighbours], setdiff(neighbours, HUBS)))
  nodes <- head(nodes, 30L)
  incident <- incident[incident$input_gene_A %in% nodes & incident$input_gene_B %in% nodes, , drop = FALSE]
  nodes_df <- read_csv_safe(file.path(string_dir, "matrisome_deg_source_table.csv"))
  node_fc <- setNames(nodes_df$logFC, nodes_df$Gene)
  graph <- igraph::graph_from_data_frame(incident[, c("input_gene_A", "input_gene_B", "score")], directed = FALSE,
                                          vertices = data.frame(name = nodes))
  fc <- node_fc[igraph::V(graph)$name]
  pal <- grDevices::colorRampPalette(c("#0072B2", "white", "#D55E00"))(101)
  finite_fc <- fc[is.finite(fc)]
  range_fc <- if (length(finite_fc)) range(finite_fc) else c(-1, 1)
  scaled <- if (diff(range_fc) == 0) rep(51, length(fc)) else pmax(1, pmin(101, round((fc - range_fc[1]) / diff(range_fc) * 100 + 1)))
  node_colors <- pal[scaled]
  node_colors[!is.finite(fc)] <- "grey80"
  grDevices::png(file.path(FIGURES, "T05_hub_STRING_network.png"), width = 2400, height = 1800, res = 300)
  plot(graph, vertex.label = igraph::V(graph)$name, vertex.color = node_colors,
       vertex.size = ifelse(igraph::V(graph)$name %in% HUBS, 26, 18),
       vertex.frame.color = ifelse(igraph::V(graph)$name %in% HUBS, "#111111", "grey60"),
       vertex.label.cex = 0.8, edge.width = 1 + 3 * incident$score,
       edge.color = grDevices::adjustcolor("#555555", alpha.f = 0.6),
       main = "STRING functional associations (score >= 0.4)")
  legend("topleft", legend = c("Negative log2FC", "Positive log2FC", "Consensus hub"),
         pch = c(21, 21, 21), pt.bg = c("#0072B2", "#D55E00", "white"),
         pt.cex = c(1.4, 1.4, 1.7), bty = "n")
  grDevices::dev.off()
  grDevices::cairo_pdf(file.path(FIGURES, "T05_hub_STRING_network.pdf"), width = 8, height = 6)
  plot(graph, vertex.label = igraph::V(graph)$name, vertex.color = node_colors,
       vertex.size = ifelse(igraph::V(graph)$name %in% HUBS, 26, 18),
       vertex.frame.color = ifelse(igraph::V(graph)$name %in% HUBS, "#111111", "grey60"),
       vertex.label.cex = 0.8, edge.width = 1 + 3 * incident$score,
       edge.color = grDevices::adjustcolor("#555555", alpha.f = 0.6),
       main = "STRING functional associations (score >= 0.4)")
  dev.off()
  write.csv(data.frame(gene = nodes, log2FC_discovery = unname(node_fc[nodes]), is_hub = nodes %in% HUBS),
            file.path(TABLES, "T05_STRING_nodes.csv"), row.names = FALSE)
  write.csv(incident, file.path(TABLES, "T05_STRING_edges.csv"), row.names = FALSE)

  correlations <- list()
  for (b in load_all_cohorts()) {
    present <- intersect(HUBS, rownames(b$expr))
    cmat <- stats::cor(t(b$expr[present, , drop = FALSE]), method = "spearman", use = "pairwise.complete.obs")
    d <- as.data.frame(as.table(cmat))
    names(d) <- c("Gene1", "Gene2", "Spearman_rho")
    d$Cohort <- b$cohort
    correlations[[b$cohort]] <- d
    pcor <- ggplot2::ggplot(d, ggplot2::aes(Gene1, Gene2, fill = Spearman_rho)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", limits = c(-1, 1)) +
      ggplot2::labs(title = paste(b$cohort, "hub-hub Spearman correlation"), x = NULL, y = NULL, fill = "rho") +
      theme_project(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_plot(pcor, paste0("T05_", b$cohort, "_hub_correlation"), 7, 6)
  }
  write.csv(do.call(rbind, correlations), file.path(TABLES, "T05_hub_correlation_long.csv"), row.names = FALSE)
  task_log("T05_interaction", c("Used the saved STRING functional-association edges for the 46-DEG network; applied score >= 0.4 and retained at most 30 connected nodes.",
    "No new first-shell expansion was available from saved sources. Functional associations are not necessarily direct physical binding.",
    paste("Retained network:", igraph::vcount(graph), "nodes,", igraph::ecount(graph), "edges; hub-hub correlations generated in three cohorts.")), "partial")
  append_results_text("T5. Hub interaction and co-expression", paste("A hub-centered network was extracted from the saved STRING functional-association network at a combined score threshold of 0.4, with at most 30 nodes.",
    "Because a fresh first-shell neighbor list was not available in the saved project files, this was not represented as a newly expanded STRING network; functional associations do not establish direct binding.",
    "Hub-hub Spearman correlations were summarized within each cohort."))
}
