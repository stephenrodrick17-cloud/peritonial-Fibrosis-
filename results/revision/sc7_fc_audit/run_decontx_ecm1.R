args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_decontx_ecm1.R <repo-root> <output-dir> <seed>")
}

root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
seed <- as.integer(args[[3]])
suppressPackageStartupMessages({
  library(celda)
  library(edgeR)
  library(Matrix)
})

single_cell_dir <- file.path(root, "results", "revision", "single_cell_context")
raw_dir <- file.path(root, "data", "raw", "GSE248762_extracted")
pseudobulk_dir <- file.path(out_dir, "decontx_pseudobulk")
dir.create(pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)

annotations <- read.csv(
  file.path(single_cell_dir, "cell_type_annotation_table.csv"),
  stringsAsFactors = FALSE
)
sample_metadata <- read.csv(
  file.path(single_cell_dir, "sample_group_qc_audit.csv"),
  stringsAsFactors = FALSE
)
sample_ids <- as.character(sample_metadata$sample_id)
target_table <- read.csv(
  file.path(root, "results", "revision", "hub_gene_network_drug_gsea",
            "validated7_kegg_gsea_summary.csv"),
  stringsAsFactors = FALSE
)
hub_genes <- toupper(as.character(target_table$gene))
saved_results <- read.csv(
  file.path(root, "results", "revision", "network_7hub",
            "sc7_hub_gene_lineage_results.csv"),
  stringsAsFactors = FALSE
)
bh_genes <- unique(toupper(as.character(saved_results$gene[
  saved_results$status == "testable" &
    is.finite(saved_results$BH_FDR_7hub_family) &
    saved_results$BH_FDR_7hub_family < 0.05
])))
if (length(bh_genes) != 1L || !bh_genes[[1]] %in% hub_genes) {
  stop("Could not resolve the single BH-significant hub from saved results.")
}
target_gene <- bh_genes[[1]]
lineages <- unique(as.character(saved_results$lineage[
  toupper(saved_results$gene) == target_gene &
    saved_results$status == "testable" &
    is.finite(saved_results$BH_FDR_7hub_family) &
    saved_results$BH_FDR_7hub_family < 0.05
]))
tested_hubs_by_lineage <- list()
corrected_by_lineage <- list()
contamination_rows <- list()
sample_summary_rows <- list()
three_group_saved <- read.csv(
  file.path(out_dir, "sc7_three_group_descriptive.csv"),
  stringsAsFactors = FALSE
)

find_file <- function(sample_id, suffix) {
  files <- list.files(raw_dir, full.names = TRUE)
  matches <- files[grepl(
    paste0("^", sample_id, "_.*", suffix, "$"),
    basename(files)
  )]
  if (length(matches) != 1L) {
    stop(sprintf("Expected one %s input for %s; found %d.", suffix, sample_id, length(matches)))
  }
  matches[[1]]
}

for (sample_id in sample_ids) {
  sample_annotations <- annotations[annotations$sample_id == sample_id, , drop = FALSE]
  if (nrow(sample_annotations) == 0L) {
    stop(sprintf("%s has no QC-retained cell annotations.", sample_id))
  }
  sample_group <- as.character(
    sample_metadata$group_from_metadata[match(sample_id, sample_metadata$sample_id)]
  )
  sample_ids_in_cell <- sub(paste0("^", sample_id, "_"), "", sample_annotations$cell_id)
  if (any(sample_ids_in_cell == sample_annotations$cell_id)) {
    stop(sprintf("%s: cell IDs do not have the expected sample prefix.", sample_id))
  }

  matrix_path <- find_file(sample_id, "\\.matrix\\.mtx\\.gz")
  barcode_path <- find_file(sample_id, "\\.barcodes\\.tsv\\.gz")
  feature_path <- find_file(sample_id, "\\.features\\.tsv\\.gz")
  raw_counts <- Matrix::readMM(gzfile(matrix_path))
  barcodes <- read.delim(
    gzfile(barcode_path), header = FALSE, stringsAsFactors = FALSE
  )[[1]]
  features <- read.delim(
    gzfile(feature_path), header = FALSE, stringsAsFactors = FALSE
  )
  if (nrow(raw_counts) != nrow(features) || ncol(raw_counts) != length(barcodes)) {
    stop(sprintf("%s: matrix dimensions do not match feature/barcode tables.", sample_id))
  }
  if (anyDuplicated(barcodes)) stop(sprintf("%s has duplicate barcodes.", sample_id))
  cell_indices <- match(sample_ids_in_cell, barcodes)
  if (anyNA(cell_indices) || anyDuplicated(cell_indices)) {
    stop(sprintf("%s: QC-retained cell barcodes do not match the raw matrix.", sample_id))
  }
  raw_counts <- Matrix::Matrix(raw_counts[, cell_indices, drop = FALSE], sparse = TRUE)
  colnames(raw_counts) <- sample_annotations$cell_id
  feature_ids <- as.character(features[[1]])
  symbols <- if (ncol(features) >= 2L) toupper(as.character(features[[2]])) else feature_ids
  missing_symbol <- is.na(symbols) | symbols == ""
  symbols[missing_symbol] <- feature_ids[missing_symbol]
  rownames(raw_counts) <- feature_ids
  cell_types <- as.factor(sample_annotations$cell_type)

  contamination_path <- file.path(
    out_dir, paste0("decontx_contamination_", sample_id, ".csv")
  )
  counts_path <- file.path(
    pseudobulk_dir, paste0("decontx_sample_counts_", sample_id, ".csv")
  )
  if (file.exists(contamination_path) && file.exists(counts_path)) {
    saved_contamination <- read.csv(contamination_path, stringsAsFactors = FALSE)
    saved_counts <- read.csv(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(saved_contamination) != nrow(sample_annotations)) {
      stop(sprintf("%s: saved DecontX outputs do not match current cells.", sample_id))
    }
    contamination_rows[[sample_id]] <- saved_contamination
    sample_summary_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      original_group = sample_group,
      n_qc_cells = nrow(sample_annotations),
      mean_contamination = mean(saved_contamination$estimated_contamination_raw),
      median_contamination = median(saved_contamination$estimated_contamination_raw),
      stringsAsFactors = FALSE
    )
    next
  }

  log_path <- file.path(out_dir, paste0("decontx_", sample_id, ".log"))
  corrected <- celda::decontX(
    raw_counts,
    z = cell_types,
    seed = seed,
    logfile = log_path,
    verbose = TRUE
  )
  corrected_counts <- corrected$decontXcounts
  if (is.null(corrected_counts) ||
      nrow(corrected_counts) != nrow(raw_counts) ||
      ncol(corrected_counts) != ncol(raw_counts)) {
    stop(sprintf("%s: DecontX returned corrected counts with unexpected dimensions.", sample_id))
  }
  estimated_contamination <- as.numeric(corrected$contamination)
  if (length(estimated_contamination) != ncol(raw_counts) ||
      any(!is.finite(estimated_contamination))) {
    stop(sprintf("%s: DecontX returned invalid contamination estimates.", sample_id))
  }

  corrected_counts <- Matrix::Matrix(round(corrected_counts), sparse = TRUE)
  symbol_levels <- unique(symbols)
  collapse_map <- Matrix::sparseMatrix(
    i = match(symbols, symbol_levels),
    j = seq_along(symbols),
    x = 1,
    dims = c(length(symbol_levels), length(symbols))
  )
  collapsed <- collapse_map %*% corrected_counts
  rownames(collapsed) <- symbol_levels
  colnames(collapsed) <- sample_annotations$cell_id

  count_frame <- data.frame(gene_symbol = symbol_levels, stringsAsFactors = FALSE)
  for (lineage in lineages) {
    local_cells <- which(sample_annotations$lineage == lineage)
    if (length(local_cells) == 0L) next
    count_frame[[lineage]] <- as.numeric(
      Matrix::rowSums(collapsed[, local_cells, drop = FALSE])
    )
  }

  contamination_frame <- data.frame(
    cell_id = sample_annotations$cell_id,
    sample_id = sample_id,
    original_group = sample_group,
    cell_type = sample_annotations$cell_type,
    lineage = sample_annotations$lineage,
    estimated_contamination_raw = estimated_contamination,
    stringsAsFactors = FALSE
  )
  write.csv(contamination_frame, contamination_path, row.names = FALSE)
  write.csv(count_frame, counts_path, row.names = FALSE, na = "")
  contamination_rows[[sample_id]] <- contamination_frame
  sample_summary_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    original_group = sample_group,
    n_qc_cells = nrow(sample_annotations),
    mean_contamination = mean(estimated_contamination),
    median_contamination = median(estimated_contamination),
    stringsAsFactors = FALSE
  )
  cat(sprintf(
    "Completed DecontX for %s (%d QC cells); mean estimate %.6g.\n",
    sample_id, nrow(sample_annotations), mean(estimated_contamination)
  ))
  flush.console()
  rm(raw_counts, corrected, corrected_counts, collapsed, collapse_map)
  gc(verbose = FALSE)
}

contamination_table <- do.call(rbind, contamination_rows)
write.csv(
  contamination_table,
  file.path(out_dir, "sc7_decontx_contamination_per_cell.csv"),
  row.names = FALSE
)
sample_summary <- do.call(rbind, sample_summary_rows)
write.csv(
  sample_summary,
  file.path(out_dir, "sc7_decontx_contamination_by_sample.csv"),
  row.names = FALSE
)

for (lineage in lineages) {
  lineage_counts <- list()
  all_symbols <- character()
  for (sample_id in sample_ids) {
    counts_file <- file.path(
      pseudobulk_dir, paste0("decontx_sample_counts_", sample_id, ".csv")
    )
    sample_counts <- read.csv(counts_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (!lineage %in% names(sample_counts)) next
    lineage_counts[[sample_id]] <- setNames(sample_counts[[lineage]], sample_counts$gene_symbol)
    all_symbols <- union(all_symbols, sample_counts$gene_symbol)
  }
  count_matrix <- matrix(
    0, nrow = length(all_symbols), ncol = length(sample_ids),
    dimnames = list(all_symbols, sample_ids)
  )
  for (sample_id in names(lineage_counts)) {
    values <- lineage_counts[[sample_id]]
    count_matrix[names(values), sample_id] <- values
  }
  write.csv(
    data.frame(gene_symbol = rownames(count_matrix), count_matrix, check.names = FALSE),
    file.path(out_dir, paste0("sc7_decontx_counts_", tolower(lineage), ".csv")),
    row.names = FALSE,
    na = ""
  )

  original_counts <- read.csv(
    file.path(single_cell_dir, "pseudobulk_counts",
              paste0("counts_", tolower(gsub("[ /-]", "_", lineage)), ".csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  original_symbols <- toupper(as.character(original_counts$gene_symbol))
  original_symbols[is.na(original_symbols) | original_symbols == ""] <-
    as.character(original_counts$feature_id)[
      is.na(original_counts$gene_symbol) | original_counts$gene_symbol == ""
    ]
  original_matrix <- as.matrix(original_counts[, sample_ids, drop = FALSE])
  storage.mode(original_matrix) <- "integer"
  original_collapsed <- rowsum(original_matrix, group = original_symbols, reorder = FALSE)
  n_cells <- read.csv(
    file.path(single_cell_dir, "sample_lineage_cell_counts.csv"),
    stringsAsFactors = FALSE
  )
  eligible_info <- n_cells[n_cells$lineage == lineage, , drop = FALSE]
  eligible_info <- eligible_info[match(sample_ids, eligible_info$sample_id), , drop = FALSE]
  eligible <- eligible_info$n_qc_cells >= 20L
  fit_info <- eligible_info[eligible, , drop = FALSE]
  corrected_fit <- count_matrix[, fit_info$sample_id, drop = FALSE]
  corrected_fit <- rowsum(corrected_fit, group = rownames(corrected_fit), reorder = FALSE)
  if (!target_gene %in% rownames(corrected_fit)) {
    stop(sprintf("%s absent from corrected %s counts.", target_gene, lineage))
  }

  original_fit <- original_collapsed[, fit_info$sample_id, drop = FALSE]
  tested_hubs <- vapply(hub_genes, function(gene) {
    if (!gene %in% rownames(original_fit)) return(FALSE)
    lv <- fit_info$primary_group == "LV"
    sv <- fit_info$primary_group == "SV"
    sum(original_fit[gene, lv] >= 10L) >= 2L &&
      sum(original_fit[gene, sv] >= 10L) >= 2L
  }, logical(1))
  tested_hubs <- hub_genes[tested_hubs]

  fit_model <- function(counts, groups, levels, contrast_name) {
    groups <- factor(as.character(groups), levels = levels)
    dge <- DGEList(counts = counts)
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- sub("^groups", "", colnames(design))
    keep <- filterByExpr(dge, group = groups)
    keep[rownames(dge) %in% tested_hubs] <- TRUE
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design, robust = TRUE)
    fit <- glmQLFit(dge, design, robust = TRUE)
    contrast <- rep(0, ncol(design))
    names(contrast) <- colnames(design)
    contrast[contrast_name] <- 1
    contrast["SV"] <- -1
    qlf <- glmQLFTest(fit, contrast = contrast)
    qlf$table[target_gene, , drop = FALSE]
  }

  metadata_groups <- as.character(
    sample_metadata$group_from_metadata[
      match(fit_info$sample_id, sample_metadata$sample_id)
    ]
  )
  comparisons <- list(
    "LV vs SV" = list(
      groups = as.character(fit_info$primary_group),
      levels = c("SV", "LV"),
      contrast_name = "LV",
      saved = saved_results[
        toupper(saved_results$gene) == target_gene &
          saved_results$lineage == lineage, ,
        drop = FALSE
      ]
    ),
    "LV_UF vs SV" = list(
      groups = metadata_groups,
      levels = c("SV", "LV_NOT_UF", "LV_UF"),
      contrast_name = "LV_UF",
      saved = three_group_saved[
        toupper(three_group_saved$gene) == target_gene &
          three_group_saved$lineage == lineage &
          three_group_saved$contrast == "LV_UF vs SV", ,
        drop = FALSE
      ]
    ),
    "LV_NOT_UF vs SV" = list(
      groups = metadata_groups,
      levels = c("SV", "LV_NOT_UF", "LV_UF"),
      contrast_name = "LV_NOT_UF",
      saved = three_group_saved[
        toupper(three_group_saved$gene) == target_gene &
          three_group_saved$lineage == lineage &
          three_group_saved$contrast == "LV_NOT_UF vs SV", ,
        drop = FALSE
      ]
    )
  )

  result_rows <- list()
  for (comparison in names(comparisons)) {
    details <- comparisons[[comparison]]
    saved_row <- details$saved
    if (nrow(saved_row) != 1L) {
      stop(sprintf("Expected one saved before-correction result for %s / %s.", lineage, comparison))
    }
    groups <- details$groups
    test_group <- details$contrast_name
    sufficient <- sum(corrected_fit[target_gene, groups == test_group] >= 10L) >= 2L &&
      sum(corrected_fit[target_gene, groups == "SV"] >= 10L) >= 2L
    if (sufficient) {
      after <- fit_model(
        corrected_fit, groups, details$levels, details$contrast_name
      )
      after_logfc <- unname(after$logFC[[1]])
      after_p <- unname(after$PValue[[1]])
      status <- "tested after DecontX"
    } else {
      after_logfc <- NA_real_
      after_p <- NA_real_
      status <- "not testable after DecontX under >=2 samples/group with >=10 UMIs"
    }
    result_rows[[comparison]] <- data.frame(
      gene = target_gene,
      lineage = lineage,
      contrast = comparison,
      n_contrast_samples = sum(groups == test_group),
      n_SV_samples = sum(groups == "SV"),
      before_log2FC = if ("log2FC" %in% names(saved_row)) saved_row$log2FC[[1]] else
        saved_row$log2FC_LV_vs_SV[[1]],
      before_P = if ("P_value" %in% names(saved_row)) saved_row$P_value[[1]] else
        saved_row$P_value[[1]],
      after_log2FC = after_logfc,
      after_P = after_p,
      after_status = status,
      stringsAsFactors = FALSE
    )
  }
  write.csv(
    do.call(rbind, result_rows),
    file.path(out_dir, paste0("sc7_decontx_contrasts_", tolower(lineage), ".csv")),
    row.names = FALSE,
    na = ""
  )
}

all_contrasts <- do.call(rbind, lapply(lineages, function(lineage) {
  read.csv(
    file.path(out_dir, paste0("sc7_decontx_contrasts_", tolower(lineage), ".csv")),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  all_contrasts,
  file.path(out_dir, "sc7_decontx_contrast_results.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    software = c("R", "edgeR", "celda", "Matrix"),
    version = c(
      as.character(getRversion()),
      as.character(packageVersion("edgeR")),
      as.character(packageVersion("celda")),
      as.character(packageVersion("Matrix"))
    ),
    seed = c(seed, seed, seed, seed),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "sc7_decontx_software_versions.csv"),
  row.names = FALSE
)
cat(sprintf(
  "Saved DecontX summaries and %d ECM1 contrast rows. celda %s; edgeR %s.\n",
  nrow(all_contrasts), as.character(packageVersion("celda")),
  as.character(packageVersion("edgeR"))
))
