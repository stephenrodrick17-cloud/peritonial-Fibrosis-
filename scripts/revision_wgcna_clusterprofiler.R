#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
WGCNA::disableWGCNAThreads()
set.seed(42)

revision_dir <- file.path("results", "revision")
dir.create(revision_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path("results", "tables", "GSE125498_full_expression_matrix.csv")
metadata_path <- file.path("results", "tables", "GSE125498_sample_metadata.csv")
ecm_path <- file.path("results", "final_analysis", "FINAL_DEG_Matrisomal_Genes.csv")

expr_df <- read.csv(expr_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)
ecm <- read.csv(ecm_path, check.names = FALSE)

if (!"Gene" %in% names(expr_df) || !all(c("sample_id", "stage_binary") %in% names(metadata))) {
  stop("Expected Gene expression identifiers and sample_id/stage_binary metadata columns.")
}
if (!all(c("Gene", "DEG_Status") %in% names(ecm))) {
  stop("Expected Gene and DEG_Status columns in the final matrisomal DEG table.")
}

expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- as.character(expr_df$Gene)

sample_ids <- as.character(metadata$sample_id)
if (!all(sample_ids %in% colnames(expr))) {
  stop("Some metadata samples are missing from the expression matrix.")
}
expr <- expr[, sample_ids, drop = FALSE]

gene_iqr <- apply(expr, 1, IQR, na.rm = TRUE)
top_variable <- names(sort(gene_iqr, decreasing = TRUE))[seq_len(min(5000L, length(gene_iqr)))]
up_ecm <- unique(as.character(ecm$Gene[ecm$DEG_Status == "Upregulated"]))
network_genes <- sort(unique(c(top_variable, intersect(up_ecm, rownames(expr)))))
dat_expr <- t(expr[network_genes, , drop = FALSE])
dat_expr <- dat_expr[, apply(dat_expr, 2, function(x) all(is.finite(x)) && sd(x) > 0), drop = FALSE]

if (ncol(dat_expr) < 1000L || nrow(dat_expr) != nrow(metadata)) {
  stop("Unexpected WGCNA input dimensions after filtering.")
}

power_vector <- c(1:10, seq(12, 30, by = 2))
sft <- WGCNA::pickSoftThreshold(
  dat_expr,
  powerVector = power_vector,
  networkType = "signed",
  RsquaredCut = 0.85,
  verbose = 5
)
fit <- as.data.frame(sft$fitIndices)
names(fit)[names(fit) == "Power"] <- "power"
names(fit)[names(fit) == "SFT.R.sq"] <- "scale_free_fit_R2"
names(fit)[names(fit) == "slope"] <- "scale_free_fit_slope"
names(fit)[names(fit) == "mean.k."] <- "mean_connectivity"
names(fit)[names(fit) == "median.k."] <- "median_connectivity"
names(fit)[names(fit) == "max.k."] <- "maximum_connectivity"
fit$selected <- FALSE

acceptable <- which(is.finite(fit$scale_free_fit_R2) & fit$scale_free_fit_R2 >= 0.85)
if (length(acceptable)) {
  selected_row <- acceptable[1]
  power_rule <- "lowest tested power with scale-free fit R2 >= 0.85"
} else {
  valid <- which(is.finite(fit$scale_free_fit_R2))
  if (!length(valid)) stop("pickSoftThreshold returned no finite scale-free fit values.")
  selected_row <- valid[which.max(fit$scale_free_fit_R2[valid])]
  power_rule <- "highest scale-free fit R2 among tested powers; none reached 0.85"
}
soft_power <- as.numeric(fit$power[selected_row])
fit$selected[selected_row] <- TRUE
write.csv(fit, file.path(revision_dir, "wgcna_soft_power_fit_taskB.csv"), row.names = FALSE)

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
  verbose = 3
)

module_labels <- as.character(net$colors)
names(module_labels) <- colnames(dat_expr)
module_table <- data.frame(
  Gene = names(module_labels),
  Module = module_labels,
  stringsAsFactors = FALSE
)
write.csv(module_table, file.path(revision_dir, "wgcna_modules_dynamic_taskB.csv"), row.names = FALSE)

module_sizes <- as.data.frame(table(Module = module_table$Module), stringsAsFactors = FALSE)
module_sizes$Gene_Count <- as.integer(module_sizes$Freq)
module_sizes$Freq <- NULL

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
module_trait <- module_trait[order(module_trait$P_value), ]
write.csv(module_trait, file.path(revision_dir, "wgcna_module_trait_dynamic_taskB.csv"), row.names = FALSE)

universe_symbols <- unique(colnames(dat_expr))
enrichment_all <- list()
modules_for_enrichment <- setdiff(unique(module_table$Module), "grey")
for (module in modules_for_enrichment) {
  genes <- module_table$Gene[module_table$Module == module]
  result <- clusterProfiler::enrichGO(
    gene = genes,
    universe = universe_symbols,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
    result_df <- as.data.frame(result)
    result_df$Module <- module
    enrichment_all[[module]] <- result_df
  }
}

if (length(enrichment_all)) {
  enrichment_df <- do.call(rbind, enrichment_all)
  enrichment_df <- enrichment_df[, c("Module", setdiff(names(enrichment_df), "Module"))]
} else {
  enrichment_df <- data.frame()
}
significant_enrichment <- if (nrow(enrichment_df)) {
  enrichment_df[!is.na(enrichment_df$p.adjust) & enrichment_df$p.adjust < 0.05, , drop = FALSE]
} else {
  enrichment_df
}
write.csv(enrichment_df, file.path(revision_dir, "clusterprofiler_GO_BP_all_taskB.csv"), row.names = FALSE)
write.csv(significant_enrichment, file.path(revision_dir, "clusterprofiler_GO_BP_FDR05_taskB.csv"), row.names = FALSE)

module_gene_count_sig_go <- if (nrow(significant_enrichment)) {
  sum(module_sizes$Gene_Count[module_sizes$Module %in% unique(significant_enrichment$Module)])
} else {
  0L
}
summary <- list(
  dataset = "GSE125498",
  expression_input = normalizePath(expr_path, winslash = "/", mustWork = TRUE),
  samples = nrow(dat_expr),
  genes_tested = ncol(dat_expr),
  gene_filter = "Top 5,000 genes by IQR plus measured Upregulated genes from FINAL_DEG_Matrisomal_Genes.csv",
  network_type = "signed",
  power_vector = power_vector,
  selected_power = soft_power,
  selected_power_rule = power_rule,
  selected_power_scale_free_fit_R2 = as.numeric(fit$scale_free_fit_R2[selected_row]),
  selected_power_mean_connectivity = as.numeric(fit$mean_connectivity[selected_row]),
  dynamic_tree_cut = list(minModuleSize = 100, deepSplit = 2),
  module_merge_cut_height = 0.25,
  total_modules_excluding_grey = length(setdiff(unique(module_labels), "grey")),
  grey_unassigned_genes = sum(module_labels == "grey"),
  module_sizes = stats::setNames(module_sizes$Gene_Count, module_sizes$Module),
  stage_associated_modules_nominal_p_lt_0_05 = sum(module_trait$P_value < 0.05, na.rm = TRUE),
  enrichment = list(
    package = "clusterProfiler",
    version = as.character(utils::packageVersion("clusterProfiler")),
    ontology = "GO Biological Process",
    adjustment = "Benjamini-Hochberg",
    significant_terms_adjusted_p_lt_0_05 = nrow(significant_enrichment),
    modules_with_significant_terms = length(unique(significant_enrichment$Module)),
    genes_in_modules_with_significant_terms = module_gene_count_sig_go
  ),
  WGCNA_version = as.character(utils::packageVersion("WGCNA")),
  R_version = as.character(getRversion()),
  random_seed = 42
)
jsonlite::write_json(
  summary,
  file.path(revision_dir, "wgcna_clusterprofiler_summary_taskB.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

cat("WGCNA/clusterProfiler run complete.\n")
cat("Samples:", nrow(dat_expr), "Genes:", ncol(dat_expr), "\n")
cat("Selected power:", soft_power, "Rule:", power_rule, "\n")
cat("Modules (excluding grey):", summary$total_modules_excluding_grey, "\n")
cat("Genes in modules:", sum(module_labels != "grey"), "Unassigned/grey:", sum(module_labels == "grey"), "\n")
cat("GO BP terms with BH-adjusted P < 0.05:", nrow(significant_enrichment), "\n")
