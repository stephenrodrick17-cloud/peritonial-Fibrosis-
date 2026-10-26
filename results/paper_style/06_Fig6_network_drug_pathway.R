source(file.path("results", "paper_style", "00_helpers.R"))

network_dir <- file.path(REV, "hub_gene_network_drug_gsea")
nodes <- read_project_csv(file.path(network_dir, "genemania_network_nodes.csv"))
edges <- read_project_csv(file.path(network_dir, "genemania_network_edges.csv"))
nodes$query_gene <- tolower(as.character(nodes$query_gene)) == "true"
stopifnot(sum(nodes$query_gene) == 7L)
edges$pair <- paste(pmin(edges$gene_a, edges$gene_b), pmax(edges$gene_a, edges$gene_b), sep = "||")
edge_summary <- aggregate(combined_weight_contribution ~ pair, edges, sum)
split_pair <- strsplit(edge_summary$pair, "\\|\\|")
edge_summary$gene_a <- vapply(split_pair, `[`, character(1), 1)
edge_summary$gene_b <- vapply(split_pair, `[`, character(1), 2)
node_names <- nodes$gene
angles <- seq(0, 2 * pi, length.out = length(node_names) + 1L)[-(length(node_names) + 1L)]
coords <- data.frame(gene = node_names, x = cos(angles), y = sin(angles))
edge_plot <- merge(edge_summary, coords, by.x = "gene_a", by.y = "gene", all = FALSE)
names(edge_plot)[names(edge_plot) == "x"] <- "x"
names(edge_plot)[names(edge_plot) == "y"] <- "y"
edge_plot <- merge(edge_plot, coords, by.x = "gene_b", by.y = "gene", all = FALSE,
                   suffixes = c("_a", "_b"))
save_csv(nodes, "Fig6A_GeneMANIA_nodes_source.csv")
save_csv(edges, "Fig6A_GeneMANIA_edges_source.csv")
save_csv(edge_summary, "Fig6A_GeneMANIA_summarized_edges.csv")
p_a <- ggplot() +
  geom_segment(data = edge_plot, aes(x = x_a, y = y_a, xend = x_b, yend = y_b,
                                     linewidth = combined_weight_contribution),
               colour = "#9BA8B4", alpha = 0.35) +
  geom_point(data = merge(coords, nodes, by.x = "gene", by.y = "gene"),
             aes(x, y, colour = query_gene), size = 3) +
  geom_text(data = merge(coords, nodes, by.x = "gene", by.y = "gene"),
            aes(x, y, label = gene), size = 2.5, nudge_y = 0.055) +
  scale_colour_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#2C7FB8"),
                      labels = c("Seven-gene input", "Additional gene")) +
  scale_linewidth(range = c(0.15, 1.1), guide = "none") +
  coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  labs(title = "GeneMANIA functional-association network",
       subtitle = "Seven-gene consensus input; up to 20 additional genes; automatic weighting",
       x = NULL, y = NULL, colour = NULL) + theme_void() +
  theme(plot.title = element_text(face = "bold", colour = "#16324F"),
        plot.subtitle = element_text(colour = "#16324F"), legend.position = "bottom")
save_plot_pair(add_letter(p_a, "A"), "Fig6A_GeneMANIA_network", 9, 8)

dsig <- read_project_csv(file.path(network_dir, "dsigdb_enrichment_full.csv"))
dsig <- dsig[order(dsig$adjusted_p_value, dsig$p_value), ]
dsig_top <- head(dsig, 12)
save_csv(dsig, "Fig6B_DSigDB_all_terms_source.csv")
save_csv(dsig_top, "Fig6B_DSigDB_top_terms.csv")
dsig_top$Signature <- factor(dsig_top$dsigdb_term,
  levels = rev(dsig_top$dsigdb_term))
p_b <- ggplot(dsig_top, aes(-log10(pmax(adjusted_p_value, 1e-300)), Signature,
                            fill = adjusted_p_value < 0.05)) +
  geom_col() + scale_fill_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#A6B4C2"),
                                 labels = c("Not FDR-significant", "BH FDR < 0.05")) +
  labs(title = "DSigDB candidate drug-induced signatures",
       subtitle = "Signature enrichment is not evidence of treatment efficacy",
       x = "-log10(BH-adjusted P)", y = NULL, fill = NULL) + theme_paper(8)
save_plot_pair(add_letter(p_b, "B"), "Fig6B_DSigDB_signatures", 9, 7)

kegg <- read_project_csv(file.path(network_dir, "kegg_gsea_figure_data.csv"))
save_csv(kegg, "Fig6C_single_gene_KEGG_GSEA_all_source.csv")
if (all(c("gene", "pathway", "NES", "FDR") %in% names(kegg))) {
  kegg$Direction <- ifelse(kegg$NES > 0, "Higher in LPD", "Higher in SPD")
  kegg_top <- do.call(rbind, lapply(split(kegg, kegg$gene), function(x) {
    x <- x[order(x$FDR, -abs(x$NES)), ]
    head(x, 5)
  }))
  save_csv(kegg_top, "Fig6C_single_gene_KEGG_GSEA_top_terms.csv")
  kegg_top$pathway <- factor(kegg_top$pathway,
    levels = rev(unique(kegg_top$pathway[order(kegg_top$NES)])))
  p_c <- ggplot(kegg_top, aes(NES, pathway, colour = FDR, size = gene_set_size)) +
    geom_point() + facet_wrap(~gene, scales = "free_y", ncol = 2) +
    scale_colour_viridis_c(trans = "log10", name = "FDR") +
    labs(title = "Single-gene preranked KEGG GSEA",
         subtitle = "Positive NES indicates higher association with LPD; non-significant terms remain labeled by FDR",
         x = "Normalized enrichment score", y = NULL, size = "Set size") +
    theme_paper(8)
  save_plot_pair(add_letter(p_c, "C"), "Fig6C_single_gene_KEGG_GSEA", 12, 10)
} else {
  save_csv(data.frame(status = "not available",
    detail = "The saved KEGG GSEA figure data did not include the expected gene, pathway, NES and FDR columns."),
    "Fig6C_KEGG_status.csv")
}
write_status("Fig6", "complete",
  "Reused seven-gene GeneMANIA, DSigDB and KEGG GSEA exports; candidate signatures are not treatments.")
