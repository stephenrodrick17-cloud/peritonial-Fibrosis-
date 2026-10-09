args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: gse130888_step2_4_edger.R <analysis-output-directory>")
}
out_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(edgeR))

counts_path <- file.path(out_dir, "all_cell_pseudobulk_counts.csv")
metadata_path <- file.path(out_dir, "all_cell_pseudobulk_metadata.csv")
symbols_path <- file.path(out_dir, "gene_id_to_symbol.csv")
counts_table <- read.csv(counts_path, row.names = 1L, check.names = FALSE)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
gene_map <- read.csv(symbols_path, stringsAsFactors = FALSE)
counts <- as.matrix(counts_table)
storage.mode(counts) <- "integer"

if (!all(colnames(counts) %in% metadata$sample_id)) {
  stop("Pseudobulk columns do not all map to GEO patient metadata")
}
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id)) {
  stop("Could not map all pseudobulk columns to patients")
}
if (any(counts < 0L)) {
  stop("Pseudobulk matrix contains negative counts")
}
if (!all(c("gene_id", "gene_symbol") %in% names(gene_map))) {
  stop("Gene symbol map is missing required columns")
}
symbol_by_id <- setNames(gene_map$gene_symbol, gene_map$gene_id)
panel_symbols <- c(
  "TNFSF15", "FLT3LG", "EBI3", "LTB",
  "ADAM19", "ECM1", "SERPINA10", "CST7"
)
panel_gene_ids <- gene_map$gene_id[
  toupper(gene_map$gene_symbol) %in% panel_symbols
]
if (length(panel_gene_ids) != length(panel_symbols) ||
    anyDuplicated(panel_gene_ids)) {
  stop("Gene-ID map does not contain exactly one feature for each panel gene")
}

primary_group <- ifelse(
  metadata$original_group == "SV",
  "SV",
  ifelse(metadata$original_group == "LV", "LV", NA_character_)
)
metadata$primary_group <- factor(primary_group, levels = c("SV", "LV"))
primary_keep <- !is.na(metadata$primary_group)
primary_meta <- metadata[primary_keep, , drop = FALSE]
primary_counts <- counts[, primary_keep, drop = FALSE]

run_ql <- function(count_matrix, sample_meta, group_factor, contrast_vector,
                   contrast_name, extra_covariates = NULL) {
  if (any(table(group_factor) < 2L)) {
    return(NULL)
  }
  sample_meta$group <- factor(group_factor)
  if (is.null(extra_covariates)) {
    design <- model.matrix(~ 0 + group, data = sample_meta)
  } else {
    design_data <- cbind(
      data.frame(group = factor(group_factor)),
      extra_covariates
    )
    design <- model.matrix(
      ~ 0 + group + monocyte_fraction + tnk_fraction,
      data = design_data
    )
  }
  if (qr(design)$rank < ncol(design)) {
    return(NULL)
  }
  contrast <- rep(0, ncol(design))
  names(contrast) <- colnames(design)
  if (!all(names(contrast_vector) %in% names(contrast))) {
    stop(
      sprintf(
        "%s: contrast columns absent from design (%s)",
        contrast_name,
        paste(names(contrast_vector), collapse = ", ")
      )
    )
  }
  contrast[names(contrast_vector)] <- contrast_vector

  y <- DGEList(counts = count_matrix)
  keep <- filterByExpr(y, design = design)
  keep[rownames(y) %in% panel_gene_ids] <- TRUE
  if (identical(contrast_name, "LV_vs_SV")) {
    detected_by_group <- vapply(levels(group_factor), function(group) {
      samples <- group_factor == group
      if (!any(samples)) {
        return(rep(FALSE, nrow(count_matrix)))
      }
      rowSums(count_matrix[, samples, drop = FALSE] >= 10L) >=
        ceiling(sum(samples) / 2)
    }, logical(nrow(count_matrix)))
    keep <- keep | apply(detected_by_group, 1L, all)
  }
  if (!any(keep)) {
    return(NULL)
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, contrast = contrast)
  table <- topTags(test, n = Inf, sort.by = "none")$table
  gene_ids <- rownames(table)
  mean_log_cpm <- rowMeans(cpm(y, log = TRUE, prior.count = 2))
  fdr_all <- p.adjust(table$PValue, method = "BH")
  design_matrix <- fit$design
  ql_variance <- fit$s2.post
  nb_dispersion <- fit$dispersion
  residual_df <- fit$df.residual.adj
  if (length(ql_variance) == 1L) {
    ql_variance <- rep(ql_variance, nrow(table))
  }
  if (length(nb_dispersion) == 1L) {
    nb_dispersion <- rep(nb_dispersion, nrow(table))
  }
  if (length(residual_df) == 1L) {
    residual_df <- rep(residual_df, nrow(table))
  }
  ci_lower <- ci_upper <- rep(NA_real_, nrow(table))
  contrast_vector <- contrast[names(contrast)]
  for (i in seq_len(nrow(table))) {
    mu <- fit$fitted.values[i, ]
    weights <- mu / (1 + nb_dispersion[i] * mu)
    information <- crossprod(design_matrix, design_matrix * weights)
    covariance <- tryCatch(
      solve(information) * ql_variance[i],
      error = function(error) NULL
    )
    if (is.null(covariance)) {
      next
    }
    # The GLM covariance is on the natural-log scale; report log2 intervals.
    variance <- drop(
      t(contrast_vector) %*% covariance %*% contrast_vector
    ) / (log(2)^2)
    if (!is.finite(variance) || variance < 0) {
      next
    }
    standard_error <- sqrt(variance)
    df <- residual_df[i]
    if (!is.finite(df) || df <= 0) {
      next
    }
    critical <- qt(0.975, df = df)
    ci_lower[i] <- table$logFC[i] - critical * standard_error
    ci_upper[i] <- table$logFC[i] + critical * standard_error
  }
  data.frame(
    contrast = contrast_name,
    gene_id = gene_ids,
    gene_symbol = unname(symbol_by_id[gene_ids]),
    logFC = table$logFC,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    PValue = table$PValue,
    FDR_all_genes = fdr_all,
    F = table$F,
    mean_logCPM = unname(mean_log_cpm[gene_ids]),
    stringsAsFactors = FALSE
  )
}

all_results <- list()
primary_results <- run_ql(
  primary_counts,
  primary_meta,
  primary_meta$primary_group,
  c(groupLV = 1, groupSV = -1),
  "LV_vs_SV"
)
if (is.null(primary_results)) {
  stop("edgeR primary LV vs SV model was not estimable")
}
all_results[["LV_vs_SV"]] <- primary_results

write.csv(
  do.call(rbind, all_results),
  file.path(out_dir, "edger_all_gene_results.csv"),
  row.names = FALSE,
  na = ""
)

logcpm_dge <- calcNormFactors(DGEList(counts = counts))
write.csv(
  cpm(logcpm_dge, log = TRUE, prior.count = 2),
  file.path(out_dir, "edger_logCPM.csv"),
  quote = FALSE
)

composition_path <- file.path(out_dir, "patient_celltype_composition.csv")
composition <- read.csv(
  composition_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
composition <- composition[
  composition$sample_id %in% primary_meta$sample_id &
    composition$n_qc_cells >= 20,
  ,
  drop = FALSE
]
composition <- composition[
  match(primary_meta$sample_id, composition$sample_id),
  ,
  drop = FALSE
]
composition_status <- "not feasible: required composition covariates unavailable"
if (all(c("fraction_monocyte/macrophage", "fraction_T/NK") %in%
        names(composition))) {
  adjusted_meta <- primary_meta
  adjusted_meta$monocyte_fraction <- composition$`fraction_monocyte/macrophage`
  adjusted_meta$tnk_fraction <- composition$`fraction_T/NK`
  n_lv <- sum(adjusted_meta$primary_group == "LV")
  n_sv <- sum(adjusted_meta$primary_group == "SV")
  adjusted_group <- adjusted_meta$primary_group
  adjusted_design_data <- data.frame(
    group = adjusted_group,
    monocyte_fraction = adjusted_meta$monocyte_fraction,
    tnk_fraction = adjusted_meta$tnk_fraction
  )
  adjusted_design <- model.matrix(
    ~ 0 + group + monocyte_fraction + tnk_fraction,
    data = adjusted_design_data
  )
  if (n_lv >= 5L && n_sv >= 5L && qr(adjusted_design)$rank == ncol(adjusted_design)) {
    adjusted <- run_ql(
      primary_counts,
      adjusted_meta,
      adjusted_meta$primary_group,
      c(groupLV = 1, groupSV = -1),
      "LV_vs_SV_composition_adjusted",
      extra_covariates = adjusted_meta[
        , c("monocyte_fraction", "tnk_fraction"), drop = FALSE
      ]
    )
    if (is.null(adjusted)) {
      stop("Composition-adjusted edgeR model returned no result")
    }
    write.csv(
      adjusted,
      file.path(out_dir, "edger_composition_adjusted_results.csv"),
      row.names = FALSE,
      na = ""
    )
    composition_status <-
      "feasible: both groups have at least five patients; design full rank"
  } else {
    composition_status <- sprintf(
      "not feasible: n_LV=%d, n_SV=%d, rank=%d/%d",
      n_lv, n_sv, qr(adjusted_design)$rank, ncol(adjusted_design)
    )
    write.csv(
      data.frame(
        contrast = character(),
        gene_id = character(),
        gene_symbol = character(),
        logFC = numeric(),
        ci_lower = numeric(),
        ci_upper = numeric(),
        PValue = numeric()
      ),
      file.path(out_dir, "edger_composition_adjusted_results.csv"),
      row.names = FALSE
    )
  }
} else {
  write.csv(
    data.frame(
      contrast = character(),
      gene_id = character(),
      gene_symbol = character(),
      logFC = numeric(),
      ci_lower = numeric(),
      ci_upper = numeric(),
      PValue = numeric()
    ),
    file.path(out_dir, "edger_composition_adjusted_results.csv"),
    row.names = FALSE
  )
}

cell_counts_path <- file.path(out_dir, "celltype_pseudobulk_counts.csv")
cell_metadata_path <- file.path(out_dir, "celltype_pseudobulk_metadata.csv")
cell_counts_frame <- read.csv(
  cell_counts_path, row.names = 1L, check.names = FALSE
)
cell_metadata <- read.csv(
  cell_metadata_path, stringsAsFactors = FALSE, check.names = FALSE
)
gene_expression <- read.csv(
  file.path(out_dir, "gene_expression_by_cell_type.csv"),
  stringsAsFactors = FALSE
)
top_cell_type <- unique(
  gene_expression[, c("gene", "gene_id", "gene_top_cell_type")]
)
celltype_results <- vector("list", nrow(top_cell_type))
for (i in seq_len(nrow(top_cell_type))) {
  row <- top_cell_type[i, ]
  selected <- cell_metadata$cell_type == row$gene_top_cell_type &
    cell_metadata$primary_group %in% c("SV", "LV")
  meta <- cell_metadata[selected, , drop = FALSE]
  if (nrow(meta) == 0L) {
    celltype_results[[i]] <- data.frame(
      gene = row$gene,
      gene_id = row$gene_id,
      top_cell_type = row$gene_top_cell_type,
      status = "no patient-cell-type pseudobulks",
      n_LV = 0L,
      n_SV = 0L,
      logFC = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      PValue = NA_real_
    )
    next
  }
  sample_cols <- paste(meta$sample_id, meta$cell_type, sep = "|")
  count_matrix <- as.matrix(cell_counts_frame[, sample_cols, drop = FALSE])
  storage.mode(count_matrix) <- "integer"
  detection <- vapply(c("SV", "LV"), function(group) {
    group_counts <- count_matrix[
      row$gene_id,
      meta$primary_group == group,
      drop = TRUE
    ]
    sum(group_counts >= 10L)
  }, integer(1))
  group_n <- table(factor(meta$primary_group, levels = c("SV", "LV")))
  detected <- detection[["SV"]] >= ceiling(group_n[["SV"]] / 2) &&
    detection[["LV"]] >= ceiling(group_n[["LV"]] / 2)
  cell_status <- if (!detected) "not detected in top cell type" else "tested"
  result <- NULL
  if (detected && all(group_n >= 2L)) {
    result <- run_ql(
      count_matrix,
      meta,
      factor(meta$primary_group, levels = c("SV", "LV")),
      c(groupLV = 1, groupSV = -1),
      paste0("top_cell_type_", row$gene)
    )
  }
  if (is.null(result) || !row$gene_id %in% result$gene_id) {
    celltype_results[[i]] <- data.frame(
      gene = row$gene,
      gene_id = row$gene_id,
      top_cell_type = row$gene_top_cell_type,
      status = if (detected) "not estimable by edgeR" else cell_status,
      n_LV = as.integer(group_n[["LV"]]),
      n_SV = as.integer(group_n[["SV"]]),
      cells_min_count_LV = detection[["LV"]],
      cells_min_count_SV = detection[["SV"]],
      logFC = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      PValue = NA_real_
    )
  } else {
    tested <- result[result$gene_id == row$gene_id, , drop = FALSE]
    celltype_results[[i]] <- data.frame(
      gene = row$gene,
      gene_id = row$gene_id,
      top_cell_type = row$gene_top_cell_type,
      status = cell_status,
      n_LV = as.integer(group_n[["LV"]]),
      n_SV = as.integer(group_n[["SV"]]),
      cells_min_count_LV = detection[["LV"]],
      cells_min_count_SV = detection[["SV"]],
      logFC = tested$logFC,
      ci_lower = tested$ci_lower,
      ci_upper = tested$ci_upper,
      PValue = tested$PValue
    )
  }
}
write.csv(
  do.call(rbind, celltype_results),
  file.path(out_dir, "edger_celltype_gene_results.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  data.frame(
    R_version = R.version.string,
    edgeR_version = as.character(packageVersion("edgeR")),
    composition_adjustment = composition_status,
    primary_patients_LV = sum(primary_meta$primary_group == "LV"),
    primary_patients_SV = sum(primary_meta$primary_group == "SV"),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "edger_analysis_config.csv"),
  row.names = FALSE
)
