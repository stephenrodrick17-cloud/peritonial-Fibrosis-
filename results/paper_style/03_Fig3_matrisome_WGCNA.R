source(file.path("results", "paper_style", "00_helpers.R"))

revision <- file.path(ROOT, "results", "revision")
soft_path <- file.path(revision, "wgcna_soft_power_fit_taskB.csv")
module_path <- file.path(revision, "wgcna_modules_dynamic_taskB.csv")
soft <- read_project_csv(soft_path)
modules <- read_project_csv(module_path)
trait <- read_project_csv(file.path(REV, "wgcna_paper_figures", "module_trait_figure_source_data.csv"))
intersection <- read_project_csv(file.path(REV, "wgcna_paper_figures", "matrisome_deg_module_intersection.csv"))

save_csv(soft, "Fig3A_soft_threshold_source.csv")
p_a_data <- tidyr::pivot_longer(soft,
  cols = c("scale_free_fit_R2", "mean_connectivity"),
  names_to = "Metric", values_to = "Value")
save_csv(p_a_data, "Fig3A_soft_threshold_plot_source.csv")
p_a <- ggplot(p_a_data, aes(power, Value)) +
  geom_line(colour = "#2C7FB8") + geom_point(colour = "#2C7FB8") +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "WGCNA soft-threshold selection", x = "Soft-thresholding power",
       y = "Saved metric") + theme_paper()
save_plot_pair(add_letter(p_a, "A"), "Fig3A_soft_threshold", 7, 5)

save_csv(modules, "Fig3B_module_assignments_source.csv")
module_counts <- as.data.frame(table(Module = modules$Module))
names(module_counts)[2] <- "Gene_count"
save_csv(module_counts, "Fig3B_module_counts.csv")
p_b_supp <- ggplot(module_counts, aes(reorder(Module, Gene_count), Gene_count, fill = Module)) +
  geom_col(show.legend = FALSE) + coord_flip() +
  scale_fill_hue() +
  labs(title = "WGCNA module membership", x = "Module", y = "Genes") + theme_paper()
save_plot_pair(p_b_supp, "Fig3B_module_counts_supplement", 6, 6)
saved_dendrogram_png <- file.path(REV, "wgcna_paper_figures", "Fig3B_gene_dendrogram_modules.png")
saved_dendrogram_pdf <- file.path(REV, "wgcna_paper_figures", "Fig3B_gene_dendrogram_modules.pdf")
if (!file.exists(saved_dendrogram_png) || !file.exists(saved_dendrogram_pdf)) {
  stop("The saved WGCNA dendrogram image is missing.")
}
dendro_png <- file.path(OUT, "Fig3B_gene_dendrogram.png")
dendro_pdf <- file.path(OUT, "Fig3B_gene_dendrogram.pdf")
new_path(dendro_png); new_path(dendro_pdf)
if (!file.copy(saved_dendrogram_png, dendro_png) || !file.copy(saved_dendrogram_pdf, dendro_pdf)) {
  stop("Could not copy the saved WGCNA dendrogram into the new output folder.")
}

save_csv(trait, "Fig3C_module_trait_source.csv")
trait$Trait <- "PD duration"
p_c <- ggplot(trait, aes(Trait, reorder(Module, Correlation), fill = Correlation)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  geom_text(aes(label = sprintf("r = %.2f\nP = %.3g", Correlation, P_value)), size = 2.3) +
  labs(title = "Module association with PD duration", x = "WGCNA module",
       y = NULL, fill = "Correlation") +
  theme_paper() + theme(axis.text.y = element_text(colour = "#16324F"),
                        panel.grid = element_blank())
save_plot_pair(add_letter(p_c, "C"), "Fig3C_module_trait", 8, 6)

deg <- read_project_csv(file.path(TABLES, "GSE125498_all_results.csv"))
ecm <- read_project_csv(file.path(TABLES, "GSE125498_ECM_intersection.csv"))
expr <- read_project_csv(file.path(TABLES, "GSE125498_full_expression_matrix.csv"))
deg_set <- unique(deg$Gene.symbol[!is.na(deg$Gene.symbol) & deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 0.585])
module_set <- unique(modules$Gene[modules$Module != "grey"])
matrisome_set <- unique(ecm$Gene.symbol[!is.na(ecm$Gene.symbol) & nzchar(ecm$Gene.symbol)])
universe <- unique(expr$Gene)
gene_sets <- list(DEG = deg_set, WGCNA_non_grey = module_set, Matrisome_measured = matrisome_set)
all_genes <- sort(unique(unlist(gene_sets)))
membership <- data.frame(Gene = all_genes)
for (nm in names(gene_sets)) membership[[nm]] <- membership$Gene %in% gene_sets[[nm]]
membership$Intersection_code <- apply(membership[, names(gene_sets), drop = FALSE], 1,
  function(x) paste(as.integer(x), collapse = ""))
save_csv(membership, "Fig3D_three_way_membership.csv")
counts <- as.data.frame(table(membership$Intersection_code), stringsAsFactors = FALSE)
names(counts) <- c("Intersection_code", "Gene_count")
counts$DEG <- substring(counts$Intersection_code, 1, 1) == "1"
counts$WGCNA_non_grey <- substring(counts$Intersection_code, 2, 2) == "1"
counts$Matrisome_measured <- substring(counts$Intersection_code, 3, 3) == "1"
save_csv(counts, "Fig3D_three_way_intersection_counts.csv")
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  stop("The VennDiagram package is required for the requested three-set Venn panel.")
}
venn_grob <- VennDiagram::venn.diagram(
  x = gene_sets, filename = NULL, fill = c("#4C78A8", "#F58518", "#54A24B"),
  alpha = 0.42, cat.col = c("#2E4C6D", "#A34E00", "#2E6930"),
  cex = 1.1, cat.cex = 0.9, margin = 0.08
)
venn_png <- file.path(OUT, "Fig3D_three_way_venn.png")
venn_pdf <- file.path(OUT, "Fig3D_three_way_venn.pdf")
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
p_d <- ggplot(counts[counts$Gene_count > 0, ], aes(reorder(Intersection_code, Gene_count), Gene_count)) +
  geom_col(fill = "#3B82A0") + geom_text(aes(label = Gene_count), hjust = -0.15, size = 3) +
  coord_flip() + labs(title = "DEG, WGCNA-module and matrisome overlap",
    subtitle = "Intersection codes follow DEG / non-grey module / measured matrisome order",
    x = "Set-membership code", y = "Genes") + theme_paper()
save_plot_pair(add_letter(p_d, "D"), "Fig3D_three_way_overlap", 8, 6)

triple_genes <- membership$Gene[membership$Intersection_code == "111"]
enrich_full <- data.frame()
enrich_error <- character()
if (length(triple_genes) >= 3 &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  go_obj <- tryCatch(clusterProfiler::enrichGO(
    gene = triple_genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1,
    qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500, readable = TRUE
  ), error = function(e) { enrich_error <<- c(enrich_error, paste("GO BP:", conditionMessage(e))); NULL })
  if (!is.null(go_obj)) {
    go_df <- as.data.frame(go_obj)
    if (nrow(go_df)) { go_df$Database <- "GO_BP"; enrich_full <- dplyr::bind_rows(enrich_full, go_df) }
  }
  universe_entrez <- tryCatch(clusterProfiler::bitr(
    universe, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE
  ), error = function(e) { enrich_error <<- c(enrich_error, paste("Background mapping:", conditionMessage(e))); data.frame() })
  genes_entrez <- tryCatch(clusterProfiler::bitr(
    triple_genes, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE
  ), error = function(e) { enrich_error <<- c(enrich_error, paste("Intersection mapping:", conditionMessage(e))); data.frame() })
  if (nrow(universe_entrez) && nrow(genes_entrez)) {
    kegg_obj <- tryCatch(clusterProfiler::enrichKEGG(
      gene = unique(genes_entrez$ENTREZID), universe = unique(universe_entrez$ENTREZID),
      organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 3, maxGSSize = 500
    ), error = function(e) { enrich_error <<- c(enrich_error, paste("KEGG:", conditionMessage(e))); NULL })
    if (!is.null(kegg_obj)) {
      kegg_df <- as.data.frame(kegg_obj)
      if (nrow(kegg_df)) { kegg_df$Database <- "KEGG"; enrich_full <- dplyr::bind_rows(enrich_full, kegg_df) }
    }
  }
}
if (nrow(enrich_full)) {
  enrich_full$FDR <- enrich_full$p.adjust
  save_csv(enrich_full, "Fig3E_intersection_enrichment_source.csv")
  enrich_sig <- enrich_full[!is.na(enrich_full$FDR) & enrich_full$FDR < 0.05, ]
  save_csv(enrich_sig, "Fig3E_intersection_enrichment_FDR05.csv")
  if (nrow(enrich_sig)) {
    enrich_top <- do.call(rbind, lapply(split(enrich_sig, enrich_sig$Database), function(x) {
      x <- x[order(x$FDR, x$pvalue), ]
      head(x, 10)
    }))
    save_csv(enrich_top, "Fig3E_top_enrichment_terms.csv")
    enrich_top$Term <- enrich_top$Description
    p_e <- ggplot(enrich_top, aes(reorder(Term, -log10(pmax(FDR, 1e-300))),
                                  Count, fill = -log10(pmax(FDR, 1e-300)))) +
      geom_col() + coord_flip() + facet_wrap(~Database, scales = "free_y") +
      scale_fill_viridis_c(name = "-log10(BH FDR)") +
      labs(title = "GO BP and KEGG enrichment of the three-way intersection",
           x = NULL, y = "Overlapping genes") + theme_paper(8)
    save_plot_pair(add_letter(p_e, "E"), "Fig3E_intersection_enrichment", 10, 8)
  } else {
    p_e <- ggplot() + annotate("text", x = 0, y = 0,
      label = "No GO BP or KEGG terms passed BH FDR < 0.05") +
      xlim(-1, 1) + ylim(-1, 1) + labs(title = "Enrichment of three-way intersection") +
      theme_void()
    save_plot_pair(add_letter(p_e, "E"), "Fig3E_intersection_enrichment", 8, 4)
    save_csv(data.frame(status = "no terms pass BH FDR < 0.05",
      detail = paste(enrich_error, collapse = "; ")), "Fig3E_enrichment_status.csv")
  }
} else {
  p_e <- ggplot() + annotate("text", x = 0, y = 0,
    label = "No GO BP or KEGG enrichment rows were returned") +
    xlim(-1, 1) + ylim(-1, 1) + labs(title = "Enrichment of three-way intersection") +
    theme_void()
  save_plot_pair(add_letter(p_e, "E"), "Fig3E_intersection_enrichment", 8, 4)
  save_csv(data.frame(status = "no enrichment rows",
    detail = paste(enrich_error, collapse = "; ")), "Fig3E_enrichment_status.csv")
}

string_edges <- read_project_csv(file.path(TABLES, "T05_STRING_edges.csv"))
string_nodes <- read_project_csv(file.path(TABLES, "T05_STRING_nodes.csv"))
string_edges <- string_edges[string_edges$score >= 0.4, ]
save_csv(string_edges, "Fig3F_saved_STRING_edges_score_ge_0_4.csv")
save_csv(string_nodes, "Fig3F_saved_STRING_nodes.csv")
if (nrow(string_edges)) {
  nodes <- sort(unique(c(string_edges$input_gene_A, string_edges$input_gene_B)))
  adj <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  for (i in seq_len(nrow(string_edges))) {
    adj[string_edges$input_gene_A[i], string_edges$input_gene_B[i]] <- 1L
    adj[string_edges$input_gene_B[i], string_edges$input_gene_A[i]] <- 1L
  }
  factorial_score <- function(n) if (n <= 1L) 1 else prod(seq_len(n - 1L))
  graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  cliques <- if (requireNamespace("igraph", quietly = TRUE)) {
    igraph::max_cliques(graph)
  } else list()
  mcc <- setNames(rep(0, length(nodes)), nodes)
  for (cl in cliques) {
    members <- igraph::V(graph)$name[as.integer(cl)]
    if (length(members) > 1L) mcc[members] <- mcc[members] + factorial_score(length(members))
  }
  mnc <- vapply(nodes, function(v) {
    nn <- names(which(adj[v, ] > 0))
    if (!length(nn)) return(0)
    if (requireNamespace("igraph", quietly = TRUE)) {
      max(igraph::components(igraph::induced_subgraph(
        igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE),
        vids = nn))$csize)
    } else length(nn)
  }, numeric(1))
  set.seed(20261016)
  epc <- setNames(rep(0, length(nodes)), nodes)
  weighted <- igraph::E(graph)
  ends <- igraph::ends(graph, weighted)
  pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
  edge_prob <- setNames(string_edges$score, pair_key(string_edges$input_gene_A, string_edges$input_gene_B))
  probabilities <- unname(edge_prob[pair_key(ends[, 1], ends[, 2])])
  if (anyNA(probabilities)) stop("A saved STRING edge could not be matched to its percolation weight.")
  for (b in seq_len(1000L)) {
    g <- graph
    if (length(weighted)) {
      keep <- stats::runif(length(weighted)) < probabilities
      g <- igraph::delete_edges(g, weighted[!keep])
    }
    comp <- igraph::components(g)
    epc <- epc + tabulate(comp$membership, nbins = max(comp$membership))[comp$membership]
  }
  epc <- epc / 1000
  centrality <- data.frame(Gene = nodes, Degree = rowSums(adj), MCC = unname(mcc[nodes]),
                           MNC = unname(mnc[nodes]), EPC = unname(epc[nodes]))
  save_csv(centrality, "Fig3F_STRING_centrality_source.csv")
  p_f <- ggplot(centrality, aes(reorder(Gene, Degree), Degree, colour = EPC, size = MCC)) +
    geom_point() + coord_flip() + scale_colour_viridis_c() +
    labs(title = "Saved STRING subnetwork centrality",
         subtitle = "MCC, EPC (1,000 edge-percolation draws), MNC and degree; subnetwork scope is limited to saved edges",
         x = NULL, y = "Degree", colour = "EPC", size = "MCC") + theme_paper()
  save_plot_pair(add_letter(p_f, "F"), "Fig3F_STRING_network_centrality", 8, 6)
} else {
  save_csv(data.frame(status = "not available", detail = "Saved STRING edge table has no score >= 0.4."),
           "Fig3F_centrality_status.csv")
}

write_status("Fig3", "partial",
  paste("Reused saved WGCNA and STRING tables; reran GO BP/KEGG enrichment on the six-gene three-way intersection;",
        nrow(string_edges), "saved STRING edges passed score >= 0.4; no new full matrisome network query was made."))
