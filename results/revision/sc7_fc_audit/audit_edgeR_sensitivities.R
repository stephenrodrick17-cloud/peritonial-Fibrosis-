args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: audit_edgeR_sensitivities.R <repo-root> <output-dir> <seed> <permutations>")
}

root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
seed <- as.integer(args[[3]])
n_permutations <- as.integer(args[[4]])
set.seed(seed)
suppressPackageStartupMessages(library(edgeR))

single_cell_dir <- file.path(root, "results", "revision", "single_cell_context")
network_dir <- file.path(root, "results", "revision", "network_7hub")
hub_table <- read.csv(
  file.path(root, "results", "revision", "hub_gene_network_drug_gsea",
            "validated7_kegg_gsea_summary.csv"),
  stringsAsFactors = FALSE
)
hub_genes <- toupper(as.character(hub_table$gene))
saved_hub_results <- read.csv(
  file.path(network_dir, "sc7_hub_gene_lineage_results.csv"),
  stringsAsFactors = FALSE
)
cell_counts <- read.csv(
  file.path(single_cell_dir, "sample_lineage_cell_counts.csv"),
  stringsAsFactors = FALSE
)
sample_audit <- read.csv(
  file.path(single_cell_dir, "sample_group_qc_audit.csv"),
  stringsAsFactors = FALSE
)

nominal <- saved_hub_results[
  saved_hub_results$status == "testable" &
    is.finite(saved_hub_results$P_value) &
    saved_hub_results$P_value < 0.05,
]
if (nrow(nominal) != 6L) {
  stop(sprintf("Expected six saved nominal pairs; found %d", nrow(nominal)))
}
bh_pairs <- nominal[nominal$BH_FDR_7hub_family < 0.05, ]
if (nrow(bh_pairs) != 3L) {
  stop(sprintf("Expected three saved BH pairs; found %d", nrow(bh_pairs)))
}

lineage_path <- function(lineage) {
  safe <- tolower(gsub("[ /-]", "_", lineage))
  file.path(single_cell_dir, "pseudobulk_counts",
            paste0("counts_", safe, ".csv"))
}

collapse_counts <- function(path) {
  table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  sample_columns <- setdiff(names(table), c("feature_id", "gene_symbol"))
  count_matrix <- as.matrix(table[, sample_columns, drop = FALSE])
  storage.mode(count_matrix) <- "integer"
  symbols <- toupper(as.character(table$gene_symbol))
  empty <- is.na(symbols) | symbols == ""
  symbols[empty] <- as.character(table$feature_id)[empty]
  collapsed <- rowsum(count_matrix, group = symbols, reorder = FALSE)
  list(table = table, counts = collapsed, samples = sample_columns)
}

pair_testability <- function(counts, eligible_info, genes = hub_genes) {
  groups <- as.character(eligible_info$primary_group)
  vapply(genes, function(gene) {
    if (!gene %in% rownames(counts)) return(FALSE)
    lv <- groups == "LV"
    sv <- groups == "SV"
    sum(counts[gene, lv] >= 10L) >= 2L &&
      sum(counts[gene, sv] >= 10L) >= 2L
  }, logical(1))
}

fit_two_group <- function(counts, group, tested_genes, target_gene) {
  group <- factor(as.character(group), levels = c("SV", "LV"))
  dge <- DGEList(counts = counts)
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
  table <- qlf$table
  list(
    table = table,
    fit = fit,
    dge = dge,
    result = if (target_gene %in% rownames(table)) table[target_gene, , drop = FALSE] else NULL
  )
}

fit_three_group <- function(counts, group, tested_genes, target_gene) {
  group_values <- as.character(group)
  group_levels <- c("SV", "LV_NOT_UF", "LV_UF")
  group <- factor(group_values, levels = group_levels)
  dge <- DGEList(counts = counts)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  keep <- filterByExpr(dge, group = group)
  keep[rownames(dge) %in% tested_genes] <- TRUE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)
  results <- list()
  for (comparison in c("LV_NOT_UF", "LV_UF")) {
    contrast <- rep(0, ncol(design))
    names(contrast) <- colnames(design)
    contrast[comparison] <- 1
    contrast["SV"] <- -1
    qlf <- glmQLFTest(fit, contrast = contrast)
    tbl <- qlf$table
    results[[comparison]] <- if (target_gene %in% rownames(tbl)) {
      tbl[target_gene, , drop = FALSE]
    } else {
      NULL
    }
  }
  results
}

sample_rows <- list()
group_rows <- list()
fc_rows <- list()
loo_rows <- list()
three_group_rows <- list()
confounding_rows <- list()
expression_rows <- list()
negative_control_rows <- list()
permutation_rows <- list()

for (lineage in unique(nominal$lineage)) {
  matrix_data <- collapse_counts(lineage_path(lineage))
  counts_all <- matrix_data$counts
  sample_ids <- matrix_data$samples
  info <- cell_counts[cell_counts$lineage == lineage, , drop = FALSE]
  info <- info[match(sample_ids, info$sample_id), , drop = FALSE]
  if (anyNA(info$sample_id)) stop(sprintf("%s: sample metadata do not match matrix", lineage))
  eligible <- info$n_qc_cells >= 20L
  fit_info <- info[eligible, , drop = FALSE]
  fit_counts <- counts_all[, fit_info$sample_id, drop = FALSE]
  tested_mask <- pair_testability(fit_counts, fit_info)
  tested_hubs <- names(tested_mask)[tested_mask]
  if (!all(nominal$gene[nominal$lineage == lineage] %in% tested_hubs)) {
    stop(sprintf("%s: a saved nominal hub fails reconstructed testability", lineage))
  }

  baseline <- fit_two_group(
    fit_counts, fit_info$primary_group, tested_hubs,
    target_gene = nominal$gene[nominal$lineage == lineage][[1]]
  )
  all_gene_table <- baseline$table
  all_gene_table$gene <- rownames(all_gene_table)
  all_gene_table$lineage <- lineage
  ql_dispersion <- baseline$fit$s2.post
  names(ql_dispersion) <- rownames(baseline$dge)
  all_gene_table$ql_dispersion <- ql_dispersion[all_gene_table$gene]
  all_gene_table$ql_dispersion[is.na(all_gene_table$ql_dispersion)] <- NA_real_
  negative_control_rows[[lineage]] <- all_gene_table

  lineage_nominal <- nominal[nominal$lineage == lineage, , drop = FALSE]
  for (pair_index in seq_len(nrow(lineage_nominal))) {
    pair <- lineage_nominal[pair_index, , drop = FALSE]
    gene <- as.character(pair$gene[[1]])
    sample_counts <- counts_all[gene, ]
    library_sizes <- colSums(counts_all)
    cpm <- sample_counts / library_sizes * 1e6
    original_group <- as.character(info$original_group)
    contrast_group <- as.character(info$primary_group)
    pair_sample <- data.frame(
      gene = gene,
      lineage = lineage,
      sample_id = sample_ids,
      original_group = original_group,
      contrast_group = contrast_group,
      n_qc_cells_in_lineage = info$n_qc_cells,
      raw_gene_count = as.numeric(sample_counts),
      total_library_size = as.numeric(library_sizes),
      CPM = as.numeric(cpm),
      log2_CPM_plus_1 = log2(cpm + 1),
      excluded_by_20_cell_rule = !eligible,
      stringsAsFactors = FALSE
    )
    sample_rows[[paste(gene, lineage, sep = "|")]] <- pair_sample

    eligible_rows <- pair_sample[!pair_sample$excluded_by_20_cell_rule, ]
    for (group_name in c("SV", "LV")) {
      grp <- eligible_rows[eligible_rows$contrast_group == group_name, ]
      group_rows[[paste(gene, lineage, group_name, sep = "|")]] <- data.frame(
        gene = gene,
        lineage = lineage,
        contrast_group = group_name,
        n_samples_eligible = nrow(grp),
        n_zero = sum(grp$raw_gene_count == 0),
        n_count_1_to_5 = sum(grp$raw_gene_count >= 1 & grp$raw_gene_count <= 5),
        n_count_gt_5 = sum(grp$raw_gene_count > 5),
        mean_raw_count = mean(grp$raw_gene_count),
        median_raw_count = median(grp$raw_gene_count),
        mean_CPM = mean(grp$CPM),
        median_CPM = median(grp$CPM),
        min_library_size = min(grp$total_library_size),
        max_library_size = max(grp$total_library_size),
        max_sample_id_by_raw_count = grp$sample_id[which.max(grp$raw_gene_count)],
        max_sample_count = max(grp$raw_gene_count),
        max_sample_CPM = max(grp$CPM),
        max_sample_CPM_fraction_of_group_sum = if (sum(grp$CPM) > 0) max(grp$CPM) / sum(grp$CPM) else NA_real_,
        group_mostly_zero = sum(grp$raw_gene_count == 0) > nrow(grp) / 2,
        one_sample_supplies_majority_CPM = if (sum(grp$CPM) > 0) max(grp$CPM) / sum(grp$CPM) > 0.5 else NA,
        stringsAsFactors = FALSE
      )
    }

    group <- as.character(fit_info$primary_group)
    group_lv <- fit_counts[gene, group == "LV"]
    group_sv <- fit_counts[gene, group == "SV"]
    lib_lv <- colSums(fit_counts)[group == "LV"]
    lib_sv <- colSums(fit_counts)[group == "SV"]
    mean_cpm_lv <- mean(group_lv / lib_lv * 1e6)
    mean_cpm_sv <- mean(group_sv / lib_sv * 1e6)
    prior_cpm_log <- log2((mean_cpm_lv + 1) / (mean_cpm_sv + 1))
    prior_logcpm <- cpm(baseline$dge, log = TRUE, prior.count = 1)
    prior_count_logfc <- mean(prior_logcpm[gene, group == "LV"]) -
      mean(prior_logcpm[gene, group == "SV"])
    saved <- saved_hub_results[
      saved_hub_results$gene == gene & saved_hub_results$lineage == lineage,
    ]
    fit_row <- baseline$table[gene, , drop = FALSE]
    dispersion <- if (gene %in% names(ql_dispersion)) {
      unname(ql_dispersion[gene])
    } else {
      NA_real_
    }
    fc_rows[[paste(gene, lineage, sep = "|")]] <- data.frame(
      gene = gene,
      lineage = lineage,
      saved_log2FC = saved$log2FC_LV_vs_SV[[1]],
      recomputed_edgeR_log2FC = fit_row$logFC[[1]],
      recomputed_edgeR_P = fit_row$PValue[[1]],
      edgeR_AveLogCPM = fit_row[["logCPM"]][[1]],
      edgeR_QL_dispersion = dispersion,
      mean_CPM_LV = mean_cpm_lv,
      mean_CPM_SV = mean_cpm_sv,
      simple_CPM_ratio_log2FC_pseudocount_1 = prior_cpm_log,
      edgeR_logCPM_prior_count_1_log2FC = prior_count_logfc,
      absolute_saved_vs_prior_count_gap = abs(saved$log2FC_LV_vs_SV[[1]] - prior_count_logfc),
      stringsAsFactors = FALSE
    )

    high_idx <- which(!pair_sample$excluded_by_20_cell_rule)
    highest_idx <- high_idx[which.max(pair_sample$raw_gene_count[high_idx])]
    loo_info <- fit_info[-match(pair_sample$sample_id[highest_idx], fit_info$sample_id), , drop = FALSE]
    loo_counts <- fit_counts[, loo_info$sample_id, drop = FALSE]
    loo_testable <- pair_testability(loo_counts, loo_info)
    loo_tested_hubs <- names(loo_testable)[loo_testable]
    if (gene %in% loo_tested_hubs) {
      loo_fit <- fit_two_group(loo_counts, loo_info$primary_group, loo_tested_hubs, gene)
      loo_row <- loo_fit$result
      loo_rows[[paste(gene, lineage, sep = "|")]] <- data.frame(
        gene = gene,
        lineage = lineage,
        omitted_sample_id = pair_sample$sample_id[highest_idx],
        omitted_original_group = pair_sample$original_group[highest_idx],
        omitted_raw_count = pair_sample$raw_gene_count[highest_idx],
        omitted_CPM = pair_sample$CPM[highest_idx],
        remaining_LV_n = sum(loo_info$primary_group == "LV"),
        remaining_SV_n = sum(loo_info$primary_group == "SV"),
        leave_one_out_log2FC = loo_row$logFC[[1]],
        leave_one_out_P = loo_row$PValue[[1]],
        leave_one_out_nominal_P_lt_0_05 = loo_row$PValue[[1]] < 0.05,
        status = "tested",
        stringsAsFactors = FALSE
      )
    } else {
      loo_rows[[paste(gene, lineage, sep = "|")]] <- data.frame(
        gene = gene,
        lineage = lineage,
        omitted_sample_id = pair_sample$sample_id[highest_idx],
        omitted_original_group = pair_sample$original_group[highest_idx],
        omitted_raw_count = pair_sample$raw_gene_count[highest_idx],
        omitted_CPM = pair_sample$CPM[highest_idx],
        remaining_LV_n = sum(loo_info$primary_group == "LV"),
        remaining_SV_n = sum(loo_info$primary_group == "SV"),
        leave_one_out_log2FC = NA_real_,
        leave_one_out_P = NA_real_,
        leave_one_out_nominal_P_lt_0_05 = NA,
        status = "not testable after omission",
        stringsAsFactors = FALSE
      )
    }

    eligible_samples <- eligible_rows$sample_id
    eligible_cpm <- pair_sample[pair_sample$sample_id %in% eligible_samples, ]
    sample_group_numeric <- as.numeric(eligible_cpm$contrast_group == "LV")
    confounding_rows[[paste(gene, lineage, sep = "|")]] <- data.frame(
      gene = gene,
      lineage = lineage,
      n_samples = nrow(eligible_cpm),
      Pearson_r_group_vs_lineage_cell_count = cor(
        sample_group_numeric, eligible_cpm$n_qc_cells_in_lineage, method = "pearson"
      ),
      Pearson_r_group_vs_library_size = cor(
        sample_group_numeric, eligible_cpm$total_library_size, method = "pearson"
      ),
      stringsAsFactors = FALSE
    )

    three_info <- fit_info
    three_info$three_group <- ifelse(
      three_info$original_group == "SV", "SV",
      ifelse(three_info$original_group == "LV_NOT_UF", "LV_NOT_UF", "LV_UF")
    )
    for (comparison in c("LV_NOT_UF", "LV_UF")) {
      three_comparison_info <- three_info
      three_comparison_counts <- fit_counts[, three_comparison_info$sample_id, drop = FALSE]
      count_testable <- vapply(c("SV", comparison), function(g) {
        sum(three_comparison_counts[gene, three_comparison_info$three_group == g] >= 10L) >= 2L
      }, logical(1))
      if (all(count_testable)) {
        three_fit <- fit_three_group(
          three_comparison_counts,
          three_comparison_info$three_group,
          gene,
          gene
        )
        three_row <- three_fit[[comparison]]
        three_group_rows[[paste(gene, lineage, comparison, sep = "|")]] <- data.frame(
          gene = gene, lineage = lineage, contrast = paste(comparison, "vs SV"),
          n_comparison_group = sum(three_comparison_info$three_group == comparison),
          n_SV = sum(three_comparison_info$three_group == "SV"),
          log2FC = three_row$logFC[[1]], P_value = three_row$PValue[[1]],
          status = "descriptive; underpowered", stringsAsFactors = FALSE
        )
      } else {
        three_group_rows[[paste(gene, lineage, comparison, sep = "|")]] <- data.frame(
          gene = gene, lineage = lineage, contrast = paste(comparison, "vs SV"),
          n_comparison_group = sum(three_comparison_info$three_group == comparison),
          n_SV = sum(three_comparison_info$three_group == "SV"),
          log2FC = NA_real_, P_value = NA_real_,
          status = "not testable by >=2 samples/group with >=10 UMIs",
          stringsAsFactors = FALSE
        )
      }
    }

    expression_rows[[paste(gene, lineage, sep = "|")]] <- data.frame(
      sample_id = eligible_cpm$sample_id,
      original_group = eligible_cpm$original_group,
      lineage = lineage,
      gene = gene,
      raw_count = eligible_cpm$raw_gene_count,
      library_size = eligible_cpm$total_library_size,
      CPM = eligible_cpm$CPM,
      stringsAsFactors = FALSE
    )
  }

  lineage_bh <- bh_pairs[bh_pairs$lineage == lineage, , drop = FALSE]
  for (pair_index in seq_len(nrow(lineage_bh))) {
    pair <- lineage_bh[pair_index, , drop = FALSE]
    gene <- as.character(pair$gene[[1]])
    observed <- baseline$table[gene, "logFC"]
    observed_abs <- abs(observed)
    perm_group <- as.character(fit_info$primary_group)
    valid_perm <- 0L
    for (iteration in seq_len(n_permutations)) {
      shuffled <- sample(perm_group, length(perm_group), replace = FALSE)
      permutation_fit <- tryCatch(
        fit_two_group(fit_counts, shuffled, tested_hubs, gene),
        error = function(e) NULL
      )
      if (is.null(permutation_fit) || is.null(permutation_fit$result)) next
      valid_perm <- valid_perm + 1L
      permutation_rows[[length(permutation_rows) + 1L]] <- data.frame(
        gene = gene,
        lineage = lineage,
        permutation = iteration,
        permuted_log2FC = permutation_fit$result$logFC[[1]],
        permuted_P = permutation_fit$result$PValue[[1]],
        observed_log2FC = observed,
        stringsAsFactors = FALSE
      )
    }
    if (valid_perm == 0L) stop(sprintf("No successful permutations for %s / %s", gene, lineage))
  }
}

write_table <- function(rows, filename) {
  if (length(rows) == 0L) {
    write.csv(data.frame(), file.path(out_dir, filename), row.names = FALSE)
  } else {
    write.csv(do.call(rbind, rows), file.path(out_dir, filename),
              row.names = FALSE, na = "")
  }
}

write_table(sample_rows, "sc7_nominal_pair_sample_counts.csv")
write_table(group_rows, "sc7_pair_group_diagnostics.csv")
write_table(fc_rows, "sc7_fold_change_comparison.csv")
write_table(loo_rows, "sc7_leave_one_out.csv")
write_table(three_group_rows, "sc7_three_group_descriptive.csv")
write_table(confounding_rows, "sc7_confounding_correlations.csv")
write_table(expression_rows, "sc7_sample_lineage_expression.csv")
if (length(negative_control_rows) > 0L) {
  write.csv(do.call(rbind, negative_control_rows),
            file.path(out_dir, "sc7_negative_control_all_genes.csv"),
            row.names = FALSE, na = "")
}
if (length(permutation_rows) > 0L) {
  permutation_table <- do.call(rbind, permutation_rows)
  write.csv(permutation_table, file.path(out_dir, "sc7_permutation_results.csv"),
            row.names = FALSE, na = "")
}
write.csv(
  data.frame(
    software = c("R", "edgeR"),
    version = c(as.character(getRversion()), as.character(packageVersion("edgeR"))),
    seed = c(seed, seed),
    permutation_iterations_requested = c(n_permutations, n_permutations)
  ),
  file.path(out_dir, "edgeR_software_versions.csv"),
  row.names = FALSE
)
