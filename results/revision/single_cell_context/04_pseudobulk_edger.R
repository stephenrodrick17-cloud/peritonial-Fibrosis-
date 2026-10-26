args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 04_pseudobulk_edger.R <single-cell-context-output-directory>")
}
out_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages(library(edgeR))

hub_genes <- c(
  "TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19",
  "NAT6", "ECM1", "CST7", "CCBE1"
)
min_cells_per_sample <- 20L
min_pseudobulk_umis <- 10L
min_detected_samples_per_group <- 2L

cell_counts <- read.csv(
  file.path(out_dir, "sample_lineage_cell_counts.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
lineages <- unique(cell_counts$lineage)
results <- list()

for (lineage in lineages) {
  safe_lineage <- tolower(gsub("[ /-]", "_", lineage))
  input_path <- file.path(
    out_dir, "pseudobulk_counts", paste0("counts_", safe_lineage, ".csv")
  )
  table <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("feature_id", "gene_symbol") %in% names(table))) {
    stop(sprintf("%s: pseudobulk matrix has missing feature columns", lineage))
  }
  sample_columns <- setdiff(names(table), c("feature_id", "gene_symbol"))
  count_matrix <- as.matrix(table[, sample_columns, drop = FALSE])
  storage.mode(count_matrix) <- "integer"
  symbols <- toupper(as.character(table$gene_symbol))
  symbols[is.na(symbols) | symbols == ""] <- as.character(table$feature_id)[
    is.na(symbols) | symbols == ""
  ]
  collapsed <- rowsum(count_matrix, group = symbols, reorder = FALSE)
  ids_by_symbol <- split(as.character(table$feature_id), symbols)

  sample_info <- cell_counts[cell_counts$lineage == lineage, , drop = FALSE]
  eligible <- sample_info[
    sample_info$n_qc_cells >= min_cells_per_sample,
    ,
    drop = FALSE
  ]
  eligible <- eligible[match(sample_columns, eligible$sample_id), , drop = FALSE]
  eligible <- eligible[!is.na(eligible$sample_id), , drop = FALSE]
  group <- factor(eligible$primary_group, levels = c("SV", "LV"))
  lv_samples <- eligible$sample_id[group == "LV"]
  sv_samples <- eligible$sample_id[group == "SV"]
  enough_replicates <- length(lv_samples) >= min_detected_samples_per_group &&
    length(sv_samples) >= min_detected_samples_per_group

  gene_counts <- if (enough_replicates) {
    collapsed[, eligible$sample_id, drop = FALSE]
  } else {
    matrix(integer(), nrow = nrow(collapsed), ncol = 0L)
  }
  testability <- lapply(hub_genes, function(gene) {
    feature_ids <- if (gene %in% names(ids_by_symbol)) {
      paste(ids_by_symbol[[gene]], collapse = ";")
    } else {
      ""
    }
    if (!gene %in% rownames(collapsed)) {
      count_lv <- count_sv <- 0L
      status <- "not_testable"
      reason <- "gene symbol absent from feature matrix"
    } else if (!enough_replicates) {
      count_lv <- count_sv <- 0L
      status <- "not_testable"
      reason <- sprintf(
        "fewer than %d eligible samples in one or both groups",
        min_detected_samples_per_group
      )
    } else {
      counts <- gene_counts[gene, ]
      count_lv <- sum(counts[eligible$primary_group == "LV"] >= min_pseudobulk_umis)
      count_sv <- sum(counts[eligible$primary_group == "SV"] >= min_pseudobulk_umis)
      passes <- count_lv >= min_detected_samples_per_group &&
        count_sv >= min_detected_samples_per_group
      status <- if (passes) "testable" else "not_testable"
      reason <- if (passes) {
        ""
      } else {
        sprintf(
          "fewer than %d samples per group had >=%d pseudobulk UMIs",
          min_detected_samples_per_group,
          min_pseudobulk_umis
        )
      }
    }
    data.frame(
      gene = gene,
      feature_ids = feature_ids,
      lineage = lineage,
      status = status,
      not_testable_reason = reason,
      n_LV_samples_with_ge_10_umis = count_lv,
      n_SV_samples_with_ge_10_umis = count_sv,
      n_LV_samples_with_ge_20_cells = length(lv_samples),
      n_SV_samples_with_ge_20_cells = length(sv_samples),
      n_LV_cells_in_eligible_samples = sum(
        eligible$n_qc_cells[eligible$primary_group == "LV"]
      ),
      n_SV_cells_in_eligible_samples = sum(
        eligible$n_qc_cells[eligible$primary_group == "SV"]
      ),
      log2FC_LV_vs_SV = NA_real_,
      P_value = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  lineage_result <- do.call(rbind, testability)

  tested_genes <- lineage_result$gene[lineage_result$status == "testable"]
  if (length(tested_genes) > 0L) {
    dge <- DGEList(counts = gene_counts)
    design <- model.matrix(~ 0 + group)
    colnames(design) <- sub("^group", "", colnames(design))
    keep <- filterByExpr(dge, group = group)
    keep[rownames(dge) %in% tested_genes] <- TRUE
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design, robust = TRUE)
    fit <- glmQLFit(dge, design, robust = TRUE)
    contrast <- rep(0, ncol(design))
    names(contrast) <- colnames(design)
    contrast["LV"] <- 1
    contrast["SV"] <- -1
    qlf <- glmQLFTest(fit, contrast = contrast)
    fitted <- qlf$table
    for (gene in tested_genes) {
      if (!gene %in% rownames(fitted)) {
        stop(sprintf("%s / %s: testable gene was removed before fitting", gene, lineage))
      }
      target <- lineage_result$gene == gene
      lineage_result$log2FC_LV_vs_SV[target] <- fitted[gene, "logFC"]
      lineage_result$P_value[target] <- fitted[gene, "PValue"]
    }
  }
  results[[lineage]] <- lineage_result
}

result_table <- do.call(rbind, results)
if (any(result_table$status == "testable" & !is.finite(result_table$P_value))) {
  stop("One or more testable gene-lineage contrasts lacks a finite P value.")
}
write.csv(
  result_table,
  file.path(out_dir, "edger_gene_lineage_results_raw.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    software = c("R", "edgeR"),
    version = c(as.character(getRversion()), as.character(packageVersion("edgeR")))
  ),
  file.path(out_dir, "pseudobulk_software_versions.csv"),
  row.names = FALSE
)
cat(
  sprintf(
    "Saved %d gene-lineage rows; %d testable contrasts across %d lineages.\n",
    nrow(result_table),
    sum(result_table$status == "testable"),
    length(lineages)
  )
)
