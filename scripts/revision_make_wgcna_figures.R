#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
WGCNA::disableWGCNAThreads()
set.seed(42)

revision_dir <- file.path("results", "revision")
figure_dir <- file.path(revision_dir, "wgcna_paper_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path("results", "tables", "GSE125498_full_expression_matrix.csv")
metadata_path <- file.path("results", "tables", "GSE125498_sample_metadata.csv")
deg_path <- file.path("results", "tables", "GSE125498_all_results.csv")
network_ecm_path <- file.path("results", "final_analysis", "FINAL_DEG_Matrisomal_Genes.csv")
module_path <- file.path(revision_dir, "wgcna_modules_dynamic_taskB.csv")
soft_power_path <- file.path(revision_dir, "wgcna_soft_power_fit_taskB.csv")
module_trait_path <- file.path(revision_dir, "wgcna_module_trait_dynamic_taskB.csv")

required_files <- c(expr_path, metadata_path, deg_path, network_ecm_path, module_path, soft_power_path, module_trait_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required WGCNA input/result files are missing: ", paste(missing_files, collapse = ", "))
}

expr_df <- read.csv(expr_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)
deg_table <- read.csv(deg_path, check.names = FALSE)
network_ecm <- read.csv(network_ecm_path, check.names = FALSE)
saved_modules <- read.csv(module_path, check.names = FALSE)
soft_power_table <- read.csv(soft_power_path, check.names = FALSE)
saved_module_trait <- read.csv(module_trait_path, check.names = FALSE)

if (!"Gene" %in% names(expr_df) ||
    !all(c("sample_id", "stage_binary") %in% names(metadata)) ||
    !all(c("Gene", "DEG_Status") %in% names(deg_table)) ||
    !all(c("Gene", "DEG_Status") %in% names(network_ecm)) ||
    !all(c("Gene", "Module") %in% names(saved_modules))) {
  stop("An input file does not contain the expected gene, module, or phenotype columns.")
}

expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- as.character(expr_df$Gene)
sample_ids <- as.character(metadata$sample_id)
if (!all(sample_ids %in% colnames(expr))) {
  stop("Some metadata samples are missing from the expression matrix.")
}
expr <- expr[, sample_ids, drop = FALSE]

gene_iqr <- apply(expr, 1L, IQR, na.rm = TRUE)
top_variable <- names(sort(gene_iqr, decreasing = TRUE))[seq_len(min(5000L, length(gene_iqr)))]
up_matrisome <- unique(as.character(network_ecm$Gene[network_ecm$DEG_Status == "Upregulated"]))
network_genes <- sort(unique(c(top_variable, intersect(up_matrisome, rownames(expr)))))
dat_expr <- t(expr[network_genes, , drop = FALSE])
dat_expr <- dat_expr[, apply(dat_expr, 2L, function(x) all(is.finite(x)) && stats::sd(x) > 0), drop = FALSE]

if (ncol(dat_expr) != 5001L || nrow(dat_expr) != nrow(metadata)) {
  stop("Reconstructed WGCNA input dimensions differ from the saved analysis.")
}

power_row <- which(soft_power_table$selected)
if (length(power_row) != 1L) {
  stop("Expected exactly one selected soft-threshold power in the saved table.")
}
soft_power <- as.numeric(soft_power_table$power[power_row])
if (soft_power != 14) {
  stop("The saved WGCNA analysis does not use the expected soft-threshold power of 14.")
}

net <- WGCNA::blockwiseModules(
  dat_expr,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 100,
  deepSplit = 2,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(dat_expr) + 1L,
  randomSeed = 42,
  nThreads = 1,
  verbose = 2
)

module_table <- data.frame(
  Gene = colnames(dat_expr),
  Module = as.character(net$colors),
  stringsAsFactors = FALSE
)
saved_match <- match(module_table$Gene, saved_modules$Gene)
if (anyNA(saved_match) || any(module_table$Module != as.character(saved_modules$Module[saved_match]))) {
  stop("Reconstructed module assignments do not match the saved WGCNA assignments; figures were not written.")
}
if (sum(module_table$Module == "grey") != 452L ||
    length(setdiff(unique(module_table$Module), "grey")) != 12L) {
  stop("Reconstructed module counts do not match the saved WGCNA summary.")
}

suppressPackageStartupMessages(library(MatrisomeAnalyzeR))

module_sizes <- as.data.frame(table(Module = module_table$Module), stringsAsFactors = FALSE)
names(module_sizes)[2] <- "Gene_Count"
mes <- WGCNA::orderMEs(net$MEs)
trait <- as.numeric(metadata$stage_binary)
module_trait <- do.call(rbind, lapply(colnames(mes), function(me) {
  test <- stats::cor.test(mes[, me], trait, method = "pearson")
  data.frame(
    Module = sub("^ME", "", me),
    Correlation = unname(test$estimate),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
module_trait <- merge(module_sizes, module_trait, by = "Module", all.x = TRUE, sort = FALSE)
module_trait$Bonferroni_P <- pmin(1, module_trait$P_value * nrow(module_trait))
module_trait <- module_trait[order(module_trait$P_value), ]
saved_trait_match <- match(module_trait$Module, saved_module_trait$Module)
if (anyNA(saved_trait_match) ||
    any(abs(module_trait$Correlation - saved_module_trait$Correlation[saved_trait_match]) > 1e-8) ||
    any(abs(module_trait$P_value - saved_module_trait$P_value[saved_trait_match]) > 1e-8)) {
  stop("Recomputed module-trait statistics do not match the saved analysis.")
}
write.csv(module_trait, file.path(figure_dir, "module_trait_figure_source_data.csv"), row.names = FALSE)

suppressPackageStartupMessages({
  matrisome_annotated <- MatrisomeAnalyzeR::matriannotate(
  data = deg_table,
  gene.column = "Gene",
  species = "human"
  )
})
matrisome_degs <- matrisome_annotated[
  matrisome_annotated[["Annotated Matrisome Division"]] != "Non-matrisome" &
    matrisome_annotated[["DEG_Status"]] != "Not Significant",
  ,
  drop = FALSE
]
matrisome_degs <- matrisome_degs[!duplicated(matrisome_degs[["Annotated Gene"]]), , drop = FALSE]
if (nrow(matrisome_degs) != 46L) {
  stop("Expected 46 unique MatrisomeAnalyzeR DEGs; found ", nrow(matrisome_degs), ".")
}

matched <- merge(
  matrisome_degs,
  module_table,
  by.x = "Annotated Gene",
  by.y = "Gene",
  all.x = TRUE,
  sort = FALSE
)
if (sum(!is.na(matched$Module)) != 43L ||
    sum(!is.na(matched$Module) & matched$Module != "grey") != 39L ||
    sum(matched$Module == "purple", na.rm = TRUE) != 3L) {
  stop("Matrisome DEG/WGCNA intersection counts do not match the expected 43/39/3.")
}
purple_genes <- sort(unique(trimws(as.character(matched[["Annotated Gene"]][matched$Module == "purple"]))))
if (!identical(purple_genes, c("CCBE1", "SEMA3E", "TNFSF15"))) {
  stop(
    "The purple-module intersection differs from the expected three genes: ",
    paste(sort(purple_genes), collapse = ", ")
  )
}
write.csv(matched, file.path(figure_dir, "matrisome_deg_module_intersection.csv"), row.names = FALSE)

network_universe <- unique(module_table$Gene)
non_grey_genes <- unique(matched[["Annotated Gene"]][!is.na(matched$Module) & matched$Module != "grey"])
go_result <- clusterProfiler::enrichGO(
  gene = non_grey_genes,
  universe = network_universe,
  OrgDb = getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"),
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)
go_table <- as.data.frame(go_result)
if (!nrow(go_table)) {
  stop("GO BP enrichment returned no terms for the 39-gene non-grey intersection.")
}
go_significant <- go_table[!is.na(go_table$p.adjust) & go_table$p.adjust < 0.05, , drop = FALSE]
if (nrow(go_significant) != 37L) {
  stop("Expected 37 BH-significant GO BP terms; found ", nrow(go_significant), ".")
}
write.csv(go_table, file.path(figure_dir, "GO_BP_enrichment_39_non_grey_intersection.csv"), row.names = FALSE)
write.csv(go_significant, file.path(figure_dir, "GO_BP_enrichment_39_non_grey_intersection_FDR05.csv"), row.names = FALSE)

save_ggplot <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

theme_paper <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(face = "bold")
  )

fit_data <- soft_power_table
fit_data$selected <- as.logical(fit_data$selected)
threshold_plot <- ggplot2::ggplot(
  fit_data,
  ggplot2::aes(x = power, y = scale_free_fit_R2)
) +
  ggplot2::geom_hline(yintercept = 0.85, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_line(colour = "#38598C", linewidth = 0.65) +
  ggplot2::geom_point(ggplot2::aes(colour = selected), size = 2.3) +
  ggplot2::geom_vline(xintercept = soft_power, linetype = "dotted", colour = "#C43C39") +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "#38598C", `TRUE` = "#C43C39"), guide = "none") +
  ggplot2::scale_x_continuous(breaks = fit_data$power) +
  ggplot2::labs(title = "Scale-free topology fit", x = "Soft-thresholding power", y = expression(R^2)) +
  theme_paper +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

connectivity_plot <- ggplot2::ggplot(
  fit_data,
  ggplot2::aes(x = power, y = mean_connectivity)
) +
  ggplot2::geom_line(colour = "#38598C", linewidth = 0.65) +
  ggplot2::geom_point(ggplot2::aes(colour = selected), size = 2.3) +
  ggplot2::geom_vline(xintercept = soft_power, linetype = "dotted", colour = "#C43C39") +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "#38598C", `TRUE` = "#C43C39"), guide = "none") +
  ggplot2::scale_x_continuous(breaks = fit_data$power) +
  ggplot2::labs(title = "Mean network connectivity", x = "Soft-thresholding power", y = "Mean connectivity") +
  theme_paper +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

slope_plot <- ggplot2::ggplot(
  fit_data,
  ggplot2::aes(x = power, y = scale_free_fit_slope)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_line(colour = "#38598C", linewidth = 0.65) +
  ggplot2::geom_point(ggplot2::aes(colour = selected), size = 2.3) +
  ggplot2::geom_vline(xintercept = soft_power, linetype = "dotted", colour = "#C43C39") +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "#38598C", `TRUE` = "#C43C39"), guide = "none") +
  ggplot2::scale_x_continuous(breaks = fit_data$power) +
  ggplot2::labs(title = "Scale-free fit slope", x = "Soft-thresholding power", y = "Linear-fit slope") +
  theme_paper +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

plot_a <- threshold_plot + slope_plot + connectivity_plot +
  patchwork::plot_layout(ncol = 3)
save_ggplot(plot_a, "Fig3A_soft_threshold", 14.2, 4.8)

draw_dendrogram <- function() {
  WGCNA::plotDendroAndColors(
    dendro = net$dendrograms[[1]],
    colors = as.character(net$colors[net$blockGenes[[1]]]),
    groupLabels = "Module",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste0(
      "GSE125498 gene dendrogram and module assignments (n = ",
      ncol(dat_expr), " genes)"
    ),
    cex.colorLabels = 0.8
  )
}
grDevices::png(
  file.path(figure_dir, "Fig3B_gene_dendrogram_modules.png"),
  width = 12,
  height = 6,
  units = "in",
  res = 600,
  type = "cairo"
)
draw_dendrogram()
grDevices::dev.off()
grDevices::cairo_pdf(file.path(figure_dir, "Fig3B_gene_dendrogram_modules.pdf"), width = 12, height = 6)
draw_dendrogram()
grDevices::dev.off()

trait_plot_data <- module_trait
trait_plot_data$Module_Label <- paste0(trait_plot_data$Module, " (n=", trait_plot_data$Gene_Count, ")")
trait_plot_data$Module_Label <- factor(
  trait_plot_data$Module_Label,
  levels = trait_plot_data$Module_Label[order(trait_plot_data$Correlation)]
)
trait_plot_data$Association <- ifelse(trait_plot_data$P_value < 0.05, "Nominal P < 0.05", "Not nominally significant")
trait_plot_data$Module_Label <- as.character(trait_plot_data$Module_Label)
trait_plot_data$Module_Label[trait_plot_data$Module == "purple"] <- paste0(
  "purple (n=", trait_plot_data$Gene_Count[trait_plot_data$Module == "purple"], "; strongest)"
)
trait_plot_data$Module_Label <- factor(
  trait_plot_data$Module_Label,
  levels = trait_plot_data$Module_Label[order(trait_plot_data$Correlation)]
)
plot_c <- ggplot2::ggplot(
  trait_plot_data,
  ggplot2::aes(x = Correlation, y = Module_Label, fill = Association)
) +
  ggplot2::geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  ggplot2::scale_fill_manual(values = c("Nominal P < 0.05" = "#B33A3A", "Not nominally significant" = "#9AA8B8")) +
  ggplot2::labs(
    title = "Module–phenotype correlations",
    subtitle = "Long-term PD versus short-term PD; nominal Pearson tests",
    x = "Pearson correlation (r)",
    y = NULL,
    fill = NULL
  ) +
  theme_paper +
  ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 9))
save_ggplot(plot_c, "Fig3C_module_trait_associations", 8.2, 6.2)

intersection_steps <- data.frame(
  Step = factor(
    c("Matrisome DEGs", "Present in WGCNA input", "Assigned to non-grey module", "Assigned to purple module"),
    levels = c("Matrisome DEGs", "Present in WGCNA input", "Assigned to non-grey module", "Assigned to purple module")
  ),
  Count = c(46L, 43L, 39L, 3L),
  Detail = c("46 unique genes", "3 not in network input", "4 assigned grey", "CCBE1, SEMA3E, TNFSF15")
)
intersection_steps$Step <- factor(intersection_steps$Step, levels = rev(levels(intersection_steps$Step)))
plot_d <- ggplot2::ggplot(intersection_steps, ggplot2::aes(x = Count, y = Step, fill = Step)) +
  ggplot2::geom_col(width = 0.62, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(Count, "  |  ", Detail)),
    hjust = -0.04,
    size = 3.5
  ) +
  ggplot2::scale_fill_manual(values = c(
    "Matrisome DEGs" = "#4976A8",
    "Present in WGCNA input" = "#6E9C83",
    "Assigned to non-grey module" = "#D49B45",
    "Assigned to purple module" = "#8667A8"
  )) +
  ggplot2::scale_x_continuous(limits = c(0, 58), breaks = seq(0, 50, 10), expand = c(0, 0)) +
  ggplot2::labs(
    title = "Matrisome DEG intersection with WGCNA modules",
    subtitle = "46 unique MatrisomeAnalyzeR DEGs; grey genes are unassigned",
    x = "Number of genes",
    y = NULL
  ) +
  theme_paper +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
save_ggplot(plot_d, "Fig3D_matrisome_module_intersection", 9, 4.4)

go_plot_data <- go_significant[order(go_significant$p.adjust), , drop = FALSE]
go_plot_data <- utils::head(go_plot_data, 15L)
go_plot_data$Description <- factor(go_plot_data$Description, levels = rev(go_plot_data$Description))
go_plot_data$GeneRatioValue <- vapply(strsplit(go_plot_data$GeneRatio, "/", fixed = TRUE), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
}, numeric(1))
plot_e <- ggplot2::ggplot(
  go_plot_data,
  ggplot2::aes(x = GeneRatioValue, y = Description, size = Count, colour = -log10(p.adjust))
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_colour_viridis_c(option = "C", name = expression(-log[10]("BH-adjusted P"))) +
  ggplot2::scale_size_continuous(range = c(2.5, 7), name = "Genes") +
  ggplot2::labs(
    title = "GO Biological Process enrichment",
    subtitle = "39 matrisome DEGs assigned to non-grey WGCNA modules; top 15 of 37 significant terms",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_paper +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8.5), legend.position = "right")
save_ggplot(plot_e, "Fig3E_GO_BP_enrichment", 9.2, 6.4)

cat("Figures written to: ", normalizePath(figure_dir, winslash = "/"), "\n", sep = "")
cat("Verified WGCNA network genes: ", ncol(dat_expr), "\n", sep = "")
cat("Verified non-grey modules / grey genes: 12 / 452\n")
cat("Verified MatrisomeAnalyzeR DEGs / network overlap / non-grey / purple: 46 / 43 / 39 / 3\n")
cat("Verified purple genes: ", paste(sort(purple_genes), collapse = ", "), "\n", sep = "")
cat("Verified BH-significant GO BP terms for 39-gene intersection: ", nrow(go_significant), "\n", sep = "")
