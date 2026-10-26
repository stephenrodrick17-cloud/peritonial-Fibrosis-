#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

input_path <- file.path("results", "final_analysis", "ML_3_models_hub_genes_all31_up_ecm.csv")
output_dir <- file.path("results", "revision", "ml_3_model_hub_overlap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Missing three-model candidate table: ", input_path)
}

results <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Gene", "LASSO_Vote", "SVM_RFE_Vote", "RF_Vote", "Total_Votes_out_of_3")
if (!all(required %in% names(results))) {
  stop("Candidate table is missing gene or model-vote columns.")
}
if (anyDuplicated(results$Gene)) {
  stop("Candidate table contains duplicated gene symbols.")
}

vote_columns <- c("LASSO_Vote", "SVM_RFE_Vote", "RF_Vote")
if (any(!as.matrix(results[, vote_columns]) %in% c(0, 1))) {
  stop("Model votes must be binary 0/1 values.")
}
results$Vote_Count <- rowSums(results[, vote_columns])
if (any(results$Vote_Count != results$Total_Votes_out_of_3)) {
  stop("Stored vote totals do not match the three model-vote columns.")
}

lasso <- results$Gene[results$LASSO_Vote == 1]
svm_rfe <- results$Gene[results$SVM_RFE_Vote == 1]
random_forest <- results$Gene[results$RF_Vote == 1]
union_genes <- sort(unique(c(lasso, svm_rfe, random_forest)))
consensus_genes <- sort(results$Gene[results$Vote_Count >= 2])
strict_genes <- sort(results$Gene[results$Vote_Count == 3])

if (length(lasso) != 10L ||
    length(svm_rfe) != 10L ||
    length(random_forest) != 10L ||
    length(union_genes) != 20L ||
    length(consensus_genes) != 8L ||
    !setequal(strict_genes, c("EBI3", "TNFSF15"))) {
  stop("The algorithm sets do not match the expected 10/10/10 set sizes, 20-gene union, and 8-gene consensus.")
}

regions <- list(
  LASSO_only = setdiff(lasso, union(svm_rfe, random_forest)),
  SVM_RFE_only = setdiff(svm_rfe, union(lasso, random_forest)),
  Random_Forest_only = setdiff(random_forest, union(lasso, svm_rfe)),
  LASSO_SVM_RFE_only = setdiff(intersect(lasso, svm_rfe), random_forest),
  LASSO_RF_only = setdiff(intersect(lasso, random_forest), svm_rfe),
  SVM_RFE_RF_only = setdiff(intersect(svm_rfe, random_forest), lasso),
  All_three = Reduce(intersect, list(lasso, svm_rfe, random_forest))
)
if (sum(lengths(regions)) != length(union_genes)) {
  stop("Venn intersection regions do not partition the union.")
}

membership <- data.frame(
  Gene = results$Gene,
  LASSO = results$LASSO_Vote == 1,
  SVM_RFE = results$SVM_RFE_Vote == 1,
  Random_Forest = results$RF_Vote == 1,
  Vote_Count = results$Vote_Count,
  In_Union = results$Vote_Count >= 1,
  Consensus_At_Least_2_of_3 = results$Vote_Count >= 2,
  Strict_3_of_3 = results$Vote_Count == 3,
  stringsAsFactors = FALSE
)
write.csv(membership, file.path(output_dir, "gene_membership_and_votes.csv"), row.names = FALSE)

region_table <- data.frame(
  Region = names(regions),
  Gene_Count = lengths(regions),
  Genes = vapply(regions, function(x) paste(sort(x), collapse = ", "), character(1)),
  stringsAsFactors = FALSE
)
write.csv(region_table, file.path(output_dir, "venn_intersection_regions.csv"), row.names = FALSE)

circle_path <- function(cx, cy, radius, label) {
  angle <- seq(0, 2 * pi, length.out = 500)
  data.frame(
    x = cx + radius * cos(angle),
    y = cy + radius * sin(angle),
    Set = label,
    stringsAsFactors = FALSE
  )
}

circles <- rbind(
  circle_path(0, 0, 1.58, "LASSO"),
  circle_path(2.10, 0, 1.58, "SVM-RFE"),
  circle_path(1.05, 1.72, 1.58, "Random Forest")
)

format_region <- function(heading, genes) {
  paste0(
    heading, " (n=", length(genes), ")",
    if (length(genes)) paste0("\n", paste(sort(genes), collapse = "\n")) else "\n—"
  )
}

region_labels <- data.frame(
  x = c(-0.66, 2.76, 1.05, 1.05, 0.18, 1.92, 1.05),
  y = c(0.12, 0.12, 2.79, -0.98, 1.22, 1.22, 0.34),
  label = c(
    format_region("LASSO only", regions$LASSO_only),
    format_region("SVM-RFE only", regions$SVM_RFE_only),
    format_region("RF only", regions$Random_Forest_only),
    format_region("LASSO + SVM-RFE", regions$LASSO_SVM_RFE_only),
    format_region("LASSO + RF", regions$LASSO_RF_only),
    format_region("SVM-RFE + RF", regions$SVM_RFE_RF_only),
    format_region("All three", regions$All_three)
  ),
  Region_Type = c(
    "Single method", "Single method", "Single method",
    "Pairwise overlap", "Pairwise overlap", "Pairwise overlap",
    "Consensus"
  ),
  stringsAsFactors = FALSE
)

plot <- ggplot2::ggplot() +
  ggplot2::geom_path(
    data = circles,
    ggplot2::aes(x = x, y = y, group = Set, colour = Set),
    linewidth = 1.1
  ) +
  ggplot2::geom_text(
    data = region_labels,
    ggplot2::aes(x = x, y = y, label = label, fontface = ifelse(Region_Type == "Consensus", "bold", "plain")),
    size = 3.15,
    lineheight = 0.93,
    colour = "#1E293B"
  ) +
  ggplot2::scale_colour_manual(
    values = c("LASSO" = "#2F6B9A", "SVM-RFE" = "#D17A22", "Random Forest" = "#4B8B62"),
    breaks = c("LASSO", "SVM-RFE", "Random Forest")
  ) +
  ggplot2::coord_equal(xlim = c(-1.82, 3.92), ylim = c(-1.72, 3.58), clip = "off") +
  ggplot2::labs(
    title = "Three-model feature-selection overlap",
    subtitle = paste0(
      "31 upregulated ECM candidates | Union: ", length(union_genes),
      " genes | Consensus (≥2/3): ", length(consensus_genes),
      " | Strict 3/3: ", length(strict_genes)
    ),
    colour = "Selection method",
    caption = "Overlap regions list the genes in each exact combination; circles are schematic and not area-proportional."
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
    plot.caption = ggplot2::element_text(size = 8.5, colour = "grey35", hjust = 0.5),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(14, 20, 12, 20)
  )

ggplot2::ggsave(
  file.path(output_dir, "three_model_gene_overlap.png"),
  plot = plot,
  width = 12,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggplot2::ggsave(
  file.path(output_dir, "three_model_gene_overlap.pdf"),
  plot = plot,
  width = 12,
  height = 9,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat("Input candidates: ", nrow(results), "\n", sep = "")
cat("Set sizes: LASSO=", length(lasso), "; SVM-RFE=", length(svm_rfe),
    "; RF=", length(random_forest), "\n", sep = "")
cat("Union of genes selected by >=1 method: ", length(union_genes), "\n", sep = "")
cat("Consensus selected by >=2 of 3 methods: ", length(consensus_genes), "\n", sep = "")
cat("Strict consensus selected by all 3 methods: ", length(strict_genes), "\n", sep = "")
cat("Figure and source tables saved in: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
