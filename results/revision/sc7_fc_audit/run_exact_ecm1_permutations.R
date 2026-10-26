args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_exact_ecm1_permutations.R <repo-root> <output-dir> <lineage>")
}

root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
lineage <- args[[3]]
suppressPackageStartupMessages(library(edgeR))

single_cell_dir <- file.path(root, "results", "revision", "single_cell_context")
network_dir <- file.path(root, "results", "revision", "network_7hub")
hub_table <- read.csv(
  file.path(root, "results", "revision", "hub_gene_network_drug_gsea",
            "validated7_kegg_gsea_summary.csv"),
  stringsAsFactors = FALSE
)
hub_genes <- toupper(as.character(hub_table$gene))
saved_results <- read.csv(
  file.path(network_dir, "sc7_hub_gene_lineage_results.csv"),
  stringsAsFactors = FALSE
)
bh_results <- saved_results[
  saved_results$status == "testable" &
    is.finite(saved_results$BH_FDR_7hub_family) &
    saved_results$BH_FDR_7hub_family < 0.05,
]
target_genes <- unique(toupper(as.character(bh_results$gene)))
if (length(target_genes) != 1L || !target_genes[[1]] %in% hub_genes) {
  stop("Could not resolve the single BH-significant hub gene from saved results.")
}
target_gene <- target_genes[[1]]
pair <- bh_results[toupper(bh_results$gene) == target_gene &
                     bh_results$lineage == lineage, ]
if (nrow(pair) != 1L) {
  stop(sprintf("Expected one saved BH result for %s / %s.", target_gene, lineage))
}

lineage_path <- file.path(
  single_cell_dir, "pseudobulk_counts",
  paste0("counts_", tolower(gsub("[ /-]", "_", lineage)), ".csv")
)
table <- read.csv(lineage_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_columns <- setdiff(names(table), c("feature_id", "gene_symbol"))
count_matrix <- as.matrix(table[, sample_columns, drop = FALSE])
storage.mode(count_matrix) <- "integer"
symbols <- toupper(as.character(table$gene_symbol))
empty <- is.na(symbols) | symbols == ""
symbols[empty] <- as.character(table$feature_id)[empty]
counts <- rowsum(count_matrix, group = symbols, reorder = FALSE)

cell_counts <- read.csv(
  file.path(single_cell_dir, "sample_lineage_cell_counts.csv"),
  stringsAsFactors = FALSE
)
info <- cell_counts[cell_counts$lineage == lineage, , drop = FALSE]
info <- info[match(sample_columns, info$sample_id), , drop = FALSE]
if (anyNA(info$sample_id)) {
  stop(sprintf("%s: sample metadata do not match the count matrix.", lineage))
}
eligible <- info$n_qc_cells >= 20L
fit_info <- info[eligible, , drop = FALSE]
fit_counts <- counts[, fit_info$sample_id, drop = FALSE]
if (!target_gene %in% rownames(fit_counts)) {
  stop(sprintf("%s is absent from the %s count matrix.", target_gene, lineage))
}

tested_mask <- vapply(hub_genes, function(gene) {
  if (!gene %in% rownames(fit_counts)) return(FALSE)
  lv <- fit_info$primary_group == "LV"
  sv <- fit_info$primary_group == "SV"
  sum(fit_counts[gene, lv] >= 10L) >= 2L &&
    sum(fit_counts[gene, sv] >= 10L) >= 2L
}, logical(1))
tested_hubs <- hub_genes[tested_mask]
if (!target_gene %in% tested_hubs) {
  stop(sprintf("%s / %s fails the saved model's testability rule.", target_gene, lineage))
}

fit_two_group <- function(group) {
  group <- factor(as.character(group), levels = c("SV", "LV"))
  dge <- DGEList(counts = fit_counts)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  keep <- filterByExpr(dge, group = group)
  keep[rownames(dge) %in% tested_hubs] <- TRUE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)
  contrast <- rep(0, ncol(design))
  names(contrast) <- colnames(design)
  contrast["LV"] <- 1
  contrast["SV"] <- -1
  qlf <- glmQLFTest(fit, contrast = contrast)
  if (!target_gene %in% rownames(qlf$table)) {
    stop(sprintf("%s was omitted from the fitted table.", target_gene))
  }
  qlf$table[target_gene, , drop = FALSE]
}

original_group <- as.character(fit_info$primary_group)
observed <- fit_two_group(original_group)
observed_logfc <- unname(observed$logFC[[1]])
observed_p <- unname(observed$PValue[[1]])
lv_count <- sum(original_group == "LV")
allocations <- combn(seq_along(original_group), lv_count)
n_allocations <- ncol(allocations)

result_path <- file.path(
  out_dir, paste0("sc7_exact_permutation_results_", tolower(lineage), ".csv")
)
summary_path <- file.path(
  out_dir, paste0("sc7_exact_permutation_summary_", tolower(lineage), ".csv")
)
completed <- integer()
if (file.exists(result_path)) {
  prior <- read.csv(result_path, stringsAsFactors = FALSE)
  completed <- as.integer(prior$allocation_index)
}
if (length(completed) > 0L && max(completed) >= n_allocations) {
  stop("Existing results contain an allocation index outside the current test space.")
}

chunk_size <- 25L
chunk <- vector("list", chunk_size)
chunk_n <- 0L
fits_succeeded <- length(completed)
started <- Sys.time()
for (allocation_index in seq_len(n_allocations)) {
  if (allocation_index %in% completed) next
  group <- rep("SV", length(original_group))
  selected <- allocations[, allocation_index]
  group[selected] <- "LV"
  fitted <- fit_two_group(group)
  chunk_n <- chunk_n + 1L
  fits_succeeded <- fits_succeeded + 1L
  chunk[[chunk_n]] <- data.frame(
    allocation_index = allocation_index,
    LV_sample_ids = paste(fit_info$sample_id[selected], collapse = ";"),
    log2FC = unname(fitted$logFC[[1]]),
    edgeR_P = unname(fitted$PValue[[1]]),
    observed_log2FC = observed_logfc,
    observed_edgeR_P = observed_p,
    stringsAsFactors = FALSE
  )
  if (chunk_n == chunk_size || allocation_index == n_allocations) {
    chunk_table <- do.call(rbind, chunk[seq_len(chunk_n)])
    write.table(
      chunk_table, result_path, sep = ",", row.names = FALSE,
      col.names = !file.exists(result_path), append = file.exists(result_path),
      na = ""
    )
    chunk <- vector("list", chunk_size)
    chunk_n <- 0L
    cat(sprintf(
      "%s: completed %d / %d allocations at %s\n",
      lineage, fits_succeeded, n_allocations, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
    flush.console()
  }
}

results <- read.csv(result_path, stringsAsFactors = FALSE)
if (nrow(results) != n_allocations ||
    anyDuplicated(results$allocation_index) ||
    any(!seq_len(n_allocations) %in% results$allocation_index)) {
  stop(sprintf("%s: exact allocation output is incomplete or duplicated.", lineage))
}
extreme_count <- sum(abs(results$log2FC) >= abs(observed_logfc))
elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
summary <- data.frame(
  gene = target_gene,
  lineage = lineage,
  n_LV = lv_count,
  n_SV = length(original_group) - lv_count,
  n_allocations = n_allocations,
  n_fits_succeeded = nrow(results),
  observed_log2FC = observed_logfc,
  observed_edgeR_P = observed_p,
  n_abs_permuted_log2FC_ge_observed = extreme_count,
  exact_two_sided_permutation_P = extreme_count / n_allocations,
  statistic = "absolute edgeR robust quasi-likelihood log2FC; refit each allocation",
  elapsed_seconds_this_invocation = elapsed_seconds,
  R_version = as.character(getRversion()),
  edgeR_version = as.character(packageVersion("edgeR")),
  stringsAsFactors = FALSE
)
write.csv(summary, summary_path, row.names = FALSE, na = "")
print(summary)
