# Post hoc / exploratory sensitivity fits use raw patient-level pseudobulk counts.
# Filtering uses edgeR::filterByExpr, then explicitly retains ECM1 so the
# requested target remains testable even if it would otherwise be filtered out.
# Normalization and quasi-likelihood testing follow the primary edgeR workflow
# (TMM, robust dispersion, glmQLFit/glmQLFTest); the approximate CI is obtained
# by inverting the one-degree-of-freedom quasi-likelihood F statistic.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: ecm1_sensitivity_edgeR.R <analysis-output-directory> <output-directory>")
}
analysis_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages(library(edgeR))

read_matrix <- function(path) {
  table <- read.csv(path, row.names = 1L, check.names = FALSE)
  matrix <- as.matrix(table)
  storage.mode(matrix) <- "integer"
  if (anyNA(matrix) || any(matrix < 0L)) {
    stop(sprintf("Invalid raw pseudobulk counts in %s", path))
  }
  matrix
}

all_counts <- read_matrix(file.path(analysis_dir, "all_cell_pseudobulk_counts.csv"))
all_meta <- read.csv(
  file.path(analysis_dir, "all_cell_pseudobulk_metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
celltype_counts <- read_matrix(file.path(analysis_dir, "celltype_pseudobulk_counts.csv"))
celltype_meta <- read.csv(
  file.path(analysis_dir, "celltype_pseudobulk_metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
gene_map <- read.csv(file.path(analysis_dir, "gene_id_to_symbol.csv"))
ecm1_ids <- gene_map$gene_id[toupper(gene_map$gene_symbol) == "ECM1"]
if (length(ecm1_ids) != 1L || !ecm1_ids %in% rownames(all_counts)) {
  stop("Expected exactly one mapped ECM1 feature in the pseudobulk matrix")
}

run_ecm1 <- function(counts, sample_meta, contrast_name) {
  if (!all(colnames(counts) %in% sample_meta$sample_id)) {
    stop(sprintf("%s: count columns do not map to sample metadata", contrast_name))
  }
  sample_meta <- sample_meta[
    match(colnames(counts), sample_meta$sample_id),
    ,
    drop = FALSE
  ]
  if (anyNA(sample_meta$sample_id)) {
    stop(sprintf("%s: missing sample metadata", contrast_name))
  }
  groups <- factor(sample_meta$primary_group, levels = c("SV", "LV"))
  if (any(table(groups) < 2L)) {
    stop(sprintf("%s: fewer than two patient pseudobulks in a group", contrast_name))
  }
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- c("groupSV", "groupLV")
  if (qr(design)$rank != ncol(design)) {
    stop(sprintf("%s: rank-deficient design", contrast_name))
  }
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)
  keep[ecm1_ids] <- TRUE
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  contrast <- c(groupSV = -1, groupLV = 1)
  test <- glmQLFTest(fit, contrast = contrast)
  tab <- topTags(test, n = Inf, sort.by = "none")$table
  gene_index <- match(ecm1_ids, rownames(tab))
  if (is.na(gene_index)) {
    stop(sprintf("%s: ECM1 absent from fitted results", contrast_name))
  }
  i <- gene_index
  df <- fit$df.residual.adj[i]
  statistic <- tab$F[i]
  estimate <- tab$logFC[i]
  if (!is.finite(statistic) || statistic <= 0 || !is.finite(df) || df <= 0) {
    stop(sprintf("%s: invalid ECM1 QL statistic or degrees of freedom", contrast_name))
  }
  # For this one-degree-of-freedom contrast, sqrt(F) is the QL t statistic.
  # Invert that same test statistic so the interval and P-value agree.
  se <- abs(estimate) / sqrt(statistic)
  data.frame(
    analysis_label = "post hoc / exploratory",
    contrast = contrast_name,
    gene = "ECM1",
    n_SV = sum(groups == "SV"),
    n_LV = sum(groups == "LV"),
    log2FC = estimate,
    ci_lower_95 = estimate - qt(0.975, df = df) * se,
    ci_upper_95 = estimate + qt(0.975, df = df) * se,
    P_value = tab$PValue[i],
    edgeR_version = as.character(packageVersion("edgeR")),
    R_version = as.character(getRversion()),
    stringsAsFactors = FALSE
  )
}

keep_without_outlier <- all_meta$sample_id != "GSM7919584"
if (!anyNA(keep_without_outlier) &&
    sum(keep_without_outlier) == ncol(all_counts)) {
  stop("GSM7919584 was not found in all-cell metadata")
}
without_outlier_meta <- all_meta[keep_without_outlier, , drop = FALSE]
without_outlier_counts <- all_counts[
  ,
  without_outlier_meta$sample_id,
  drop = FALSE
]

mesothelial_columns <- grep("\\|mesothelial$", colnames(celltype_counts), value = TRUE)
mesothelial_samples <- sub("\\|mesothelial$", "", mesothelial_columns)
mesothelial_meta <- celltype_meta[
  celltype_meta$cell_type == "mesothelial" &
    celltype_meta$sample_id %in% mesothelial_samples &
    celltype_meta$n_qc_cells >= 20L,
  ,
  drop = FALSE
]
mesothelial_columns <- paste0(mesothelial_meta$sample_id, "|mesothelial")
if (!all(mesothelial_columns %in% colnames(celltype_counts))) {
  stop("Mesothelial pseudobulk metadata and matrix columns differ")
}
mesothelial_counts <- celltype_counts[
  ,
  mesothelial_columns,
  drop = FALSE
]
mesothelial_sample_ids <- sub(
  "\\|mesothelial$",
  "",
  colnames(mesothelial_counts)
)
mesothelial_meta <- mesothelial_meta[
  match(mesothelial_sample_ids, mesothelial_meta$sample_id),
  ,
  drop = FALSE
]
colnames(mesothelial_counts) <- mesothelial_sample_ids
if (anyNA(mesothelial_meta$sample_id)) {
  stop("Could not align mesothelial pseudobulk columns to patient metadata")
}

results <- rbind(
  run_ecm1(
    without_outlier_counts,
    without_outlier_meta,
    "LV_vs_SV_excluding_GSM7919584"
  ),
  run_ecm1(
    mesothelial_counts,
    mesothelial_meta,
    "mesothelial_only_LV_vs_SV"
  )
)
write.csv(
  results,
  file.path(output_dir, "ecm1_posthoc_edgeR_sensitivities.csv"),
  row.names = FALSE,
  na = ""
)
cat("Post hoc / exploratory ECM1 edgeR sensitivities\n")
print(results, row.names = FALSE)
