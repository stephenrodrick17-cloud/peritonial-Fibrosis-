args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_exact_ecm1_permutations_parallel.R <repo-root> <output-dir> <workers>")
}

root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
workers <- as.integer(args[[3]])
if (!is.finite(workers) || workers < 1L) stop("workers must be a positive integer.")

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
  stop("Could not resolve the single BH-significant hub from saved results.")
}
target_gene <- target_genes[[1]]
lineages <- unique(as.character(bh_results$lineage))
cell_counts <- read.csv(
  file.path(single_cell_dir, "sample_lineage_cell_counts.csv"),
  stringsAsFactors = FALSE
)

collapse_counts <- function(path) {
  table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  sample_columns <- setdiff(names(table), c("feature_id", "gene_symbol"))
  count_matrix <- as.matrix(table[, sample_columns, drop = FALSE])
  storage.mode(count_matrix) <- "integer"
  symbols <- toupper(as.character(table$gene_symbol))
  empty <- is.na(symbols) | symbols == ""
  symbols[empty] <- as.character(table$feature_id)[empty]
  list(counts = rowsum(count_matrix, group = symbols, reorder = FALSE),
       samples = sample_columns)
}

fit_exact_one <- function(allocation_index) {
  selected <- allocations[, allocation_index]
  group <- rep("SV", length(sample_ids))
  group[selected] <- "LV"
  group <- factor(group, levels = c("SV", "LV"))
  dge <- edgeR::DGEList(counts = fit_counts)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  keep <- edgeR::filterByExpr(dge, group = group)
  keep[rownames(dge) %in% tested_hubs] <- TRUE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  contrast <- rep(0, ncol(design))
  names(contrast) <- colnames(design)
  contrast["LV"] <- 1
  contrast["SV"] <- -1
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  if (!target_gene %in% rownames(qlf$table)) {
    stop(sprintf("%s was omitted from an exact allocation fit.", target_gene))
  }
  data.frame(
    allocation_index = allocation_index,
    LV_sample_ids = paste(sample_ids[selected], collapse = ";"),
    log2FC = unname(qlf$table[target_gene, "logFC"]),
    edgeR_P = unname(qlf$table[target_gene, "PValue"]),
    observed_log2FC = observed_logfc,
    observed_edgeR_P = observed_p,
    stringsAsFactors = FALSE
  )
}

started <- Sys.time()
cluster <- parallel::makeCluster(workers, type = "PSOCK")
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterEvalQ(cluster, suppressPackageStartupMessages(library(edgeR)))

for (lineage in lineages) {
  matrix_data <- collapse_counts(file.path(
    single_cell_dir, "pseudobulk_counts",
    paste0("counts_", tolower(gsub("[ /-]", "_", lineage)), ".csv")
  ))
  info <- cell_counts[cell_counts$lineage == lineage, , drop = FALSE]
  info <- info[match(matrix_data$samples, info$sample_id), , drop = FALSE]
  if (anyNA(info$sample_id)) {
    stop(sprintf("%s: sample metadata do not match the count matrix.", lineage))
  }
  eligible <- info$n_qc_cells >= 20L
  fit_info <- info[eligible, , drop = FALSE]
  sample_ids <- as.character(fit_info$sample_id)
  fit_counts <- matrix_data$counts[, sample_ids, drop = FALSE]
  original_group <- as.character(fit_info$primary_group)
  lv_count <- sum(original_group == "LV")
  if (lv_count < 2L || sum(original_group == "SV") < 2L) {
    stop(sprintf("%s: insufficient eligible samples for exact enumeration.", lineage))
  }

  tested_mask <- vapply(hub_genes, function(gene) {
    if (!gene %in% rownames(fit_counts)) return(FALSE)
    lv <- original_group == "LV"
    sv <- original_group == "SV"
    sum(fit_counts[gene, lv] >= 10L) >= 2L &&
      sum(fit_counts[gene, sv] >= 10L) >= 2L
  }, logical(1))
  tested_hubs <- hub_genes[tested_mask]
  if (!target_gene %in% tested_hubs) {
    stop(sprintf("%s / %s fails the saved model's testability rule.", target_gene, lineage))
  }

  group_factor <- factor(original_group, levels = c("SV", "LV"))
  dge <- DGEList(counts = fit_counts)
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- sub("^group_factor", "", colnames(design))
  keep <- filterByExpr(dge, group = group_factor)
  keep[rownames(dge) %in% tested_hubs] <- TRUE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)
  contrast <- rep(0, ncol(design))
  names(contrast) <- colnames(design)
  contrast["LV"] <- 1
  contrast["SV"] <- -1
  baseline_table <- glmQLFTest(fit, contrast = contrast)$table
  if (!target_gene %in% rownames(baseline_table)) {
    stop(sprintf("%s / %s baseline fit is missing the target.", target_gene, lineage))
  }
  observed_logfc <- unname(baseline_table[target_gene, "logFC"])
  observed_p <- unname(baseline_table[target_gene, "PValue"])
  allocations <- combn(seq_along(sample_ids), lv_count)
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
    stop(sprintf("%s: existing results exceed current allocation space.", lineage))
  }
  remaining <- setdiff(seq_len(n_allocations), completed)
  parallel::clusterExport(
    cluster,
    c("fit_exact_one", "fit_counts", "sample_ids", "tested_hubs",
      "target_gene", "observed_logfc", "observed_p", "allocations"),
    envir = environment()
  )
  if (length(remaining) > 0L) {
    chunk_size <- max(1L, workers * 4L)
    for (start in seq.int(1L, length(remaining), by = chunk_size)) {
      chunk_indices <- remaining[start:min(start + chunk_size - 1L, length(remaining))]
      chunk_table <- do.call(
        rbind,
        parallel::parLapplyLB(cluster, chunk_indices, fit_exact_one)
      )
      write.table(
        chunk_table, result_path, sep = ",", row.names = FALSE,
        col.names = !file.exists(result_path), append = file.exists(result_path),
        na = ""
      )
      cat(sprintf(
        "%s: completed %d / %d allocations at %s\n",
        lineage, length(completed) + min(start + length(chunk_indices) - 1L, length(remaining)),
        n_allocations, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ))
      flush.console()
    }
  }

  results <- read.csv(result_path, stringsAsFactors = FALSE)
  if (nrow(results) != n_allocations ||
      anyDuplicated(results$allocation_index) ||
      any(!seq_len(n_allocations) %in% results$allocation_index)) {
    stop(sprintf("%s: exact output is incomplete or duplicated.", lineage))
  }
  extreme_count <- sum(abs(results$log2FC) >= abs(observed_logfc))
  summary <- data.frame(
    gene = target_gene,
    lineage = lineage,
    n_LV = lv_count,
    n_SV = length(sample_ids) - lv_count,
    n_allocations = n_allocations,
    n_fits_succeeded = nrow(results),
    observed_log2FC = observed_logfc,
    observed_edgeR_P = observed_p,
    n_abs_permuted_log2FC_ge_observed = extreme_count,
    exact_two_sided_permutation_P = extreme_count / n_allocations,
    statistic = "absolute edgeR robust quasi-likelihood log2FC; refit each allocation",
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    n_workers = workers,
    R_version = as.character(getRversion()),
    edgeR_version = as.character(packageVersion("edgeR")),
    stringsAsFactors = FALSE
  )
  write.csv(summary, summary_path, row.names = FALSE, na = "")
  print(summary)
}
