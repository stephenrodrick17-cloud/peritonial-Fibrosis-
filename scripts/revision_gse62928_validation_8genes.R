#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(limma))
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript revision_gse62928_validation_8genes.R <series_matrix.gz> <GPL13158.soft> <output_dir>")
}
matrix_path <- args[[1]]
platform_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel8 <- c("TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "ECM1", "SERPINA10", "CST7")
strict2 <- c("TNFSF15", "FLT3LG")
expected_discovery_direction <- setNames(rep("UP", length(panel8)), panel8)
bootstrap_reps <- 10000L
null_reps <- 1000L
detection_cutoff <- 5.0

read_gse_preamble <- function(path) {
  con <- gzfile(path, open = "rt")
  lines <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("GEO series matrix ended before its expression table.")
    if (startsWith(line, "!series_matrix_table_begin")) break
    lines <- c(lines, line)
  }
  list(lines = lines, connection = con)
}

parse_geo_vector <- function(lines, prefix) {
  line <- lines[startsWith(lines, prefix)]
  if (length(line) != 1L) stop("Expected exactly one GEO metadata row for ", prefix)
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  trimws(gsub('^"|"$', "", fields[-1]))
}

geo <- read_gse_preamble(matrix_path)
platform_line <- geo$lines[startsWith(geo$lines, "!Series_platform_id")]
if (length(platform_line) != 1L || !grepl("GPL13158", platform_line, fixed = TRUE)) {
  stop("Downloaded matrix does not identify GPL13158.")
}
sample_ids <- parse_geo_vector(geo$lines, "!Sample_geo_accession")
sample_titles <- parse_geo_vector(geo$lines, "!Sample_title")
sample_desc <- parse_geo_vector(geo$lines, "!Sample_description")
char_lines <- geo$lines[startsWith(geo$lines, "!Sample_characteristics_ch1")]
if (!length(char_lines)) stop("GEO sample group characteristics are missing.")
char_mat <- lapply(char_lines, function(line) {
  trimws(gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]][-1]))
})
if (any(lengths(char_mat) != length(sample_ids))) stop("GEO sample characteristics are misaligned.")
char_text <- vapply(seq_along(sample_ids), function(i) {
  paste(vapply(char_mat, `[[`, character(1), i), collapse = "; ")
}, character(1))

group3 <- ifelse(grepl("^EPS", sample_titles, ignore.case = TRUE), "EPS",
                 ifelse(grepl("^PD", sample_titles, ignore.case = TRUE), "Early PD",
                        ifelse(grepl("^UREMIC", sample_titles, ignore.case = TRUE), "Uremic", NA_character_)))
if (anyNA(group3)) stop("A sample title does not map to EPS, early PD, or uremic.")
if (any(group3 == "EPS" & !grepl("treatment: EPS", char_text, ignore.case = TRUE)) ||
    any(group3 == "Early PD" & !grepl("treatment: PD", char_text, ignore.case = TRUE)) ||
    any(group3 == "Uremic" & !grepl("UREMIC", char_text, ignore.case = TRUE))) {
  stop("GEO title and treatment characteristics disagree for at least one sample.")
}
group_primary <- ifelse(group3 == "EPS", "EPS", "Non-EPS")
subgroup_basis <- ifelse(
  group3 == "Early PD",
  "GEO title/treatment says PD; early-PD status follows supplied cohort description",
  ifelse(group3 == "EPS",
         "GEO identifies EPS; severe-EPS status follows supplied cohort description",
         "GEO title/treatment metadata")
)
mapping <- data.frame(
  sample_id = sample_ids,
  GEO_title = sample_titles,
  GEO_description = sample_desc,
  GEO_characteristics = char_text,
  subgroup = group3,
  subgroup_assignment_basis = subgroup_basis,
  primary_group = group_primary
)
write.csv(mapping, file.path(out_dir, "sample_group_mapping.csv"), row.names = FALSE)
observed_counts <- table(factor(group3, levels = c("EPS", "Early PD", "Uremic")))
if (!identical(as.integer(observed_counts), c(4L, 2L, 2L))) {
  stop(
    "Sample groups differ from the preregistered cohort description: ",
    paste(names(observed_counts), as.integer(observed_counts), collapse = "; ")
  )
}
if (length(sample_ids) != 8L || anyDuplicated(sample_ids)) stop("Expected eight unique validation samples.")

expr_df <- read.delim(
  geo$connection,
  header = TRUE,
  check.names = FALSE,
  quote = "\"",
  comment.char = "",
  stringsAsFactors = FALSE
)
close(geo$connection)
expr_df <- expr_df[!grepl("^!", as.character(expr_df[[1]])), , drop = FALSE]
probe_ids <- as.character(expr_df[[1]])
expr <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- probe_ids
if (!all(sample_ids %in% colnames(expr))) stop("Expression matrix sample IDs do not match GEO metadata.")
expr <- expr[, sample_ids, drop = FALSE]
if (any(!is.finite(expr))) stop("Expression matrix contains non-finite values.")

processing <- geo$lines[startsWith(geo$lines, "!Sample_data_processing")]
processing_text <- paste(processing, collapse = " ")
if (!grepl("RMA", processing_text, ignore.case = TRUE)) {
  stop("GEO metadata does not confirm RMA normalization; stop and review preprocessing.")
}

annotation_lines <- readLines(platform_path, warn = FALSE)
begin <- which(annotation_lines == "!platform_table_begin")
end <- which(annotation_lines == "!platform_table_end")
if (length(begin) != 1L || length(end) != 1L || end <= begin) {
  stop("Could not identify GPL13158 annotation table.")
}
annotation <- read.delim(
  text = paste(annotation_lines[(begin + 1L):(end - 1L)], collapse = "\n"),
  check.names = FALSE,
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)
if (!all(c("ID", "Gene Symbol") %in% names(annotation))) {
  stop("GPL13158 annotation lacks ID and Gene Symbol columns.")
}
annotation <- annotation[, c("ID", "Gene Symbol")]
annotation$ID <- as.character(annotation$ID)
annotation$Gene_Symbol <- as.character(annotation[["Gene Symbol"]])
annotation <- annotation[annotation$ID %in% probe_ids, , drop = FALSE]
annotation$Gene.Symbols <- lapply(annotation$Gene_Symbol, function(x) {
  if (is.na(x) || !nzchar(trimws(x)) || x %in% c("---", "NA")) return(character())
  symbols <- trimws(strsplit(x, "///", fixed = TRUE)[[1]])
  symbols[nzchar(symbols) & !symbols %in% c("---", "NA")]
})

expanded <- do.call(rbind, lapply(seq_len(nrow(annotation)), function(i) {
  symbols <- annotation$Gene.Symbols[[i]]
  if (!length(symbols)) return(NULL)
  data.frame(Probe_ID = annotation$ID[i], Gene = symbols, stringsAsFactors = FALSE)
}))
if (is.null(expanded)) stop("GPL13158 annotation produced no gene mappings.")
expanded <- unique(expanded)
panel_probes <- expanded[expanded$Gene %in% panel8, , drop = FALSE]
missing_genes <- setdiff(panel8, unique(panel_probes$Gene))
if (length(missing_genes)) stop("Fixed gene(s) have no GPL13158 probe mapping: ", paste(missing_genes, collapse = ", "))

probe_stats <- data.frame(
  Probe_ID = probe_ids,
  Mean_Log2_Expression = rowMeans(expr),
  IQR_Log2_Expression = apply(expr, 1, IQR),
  stringsAsFactors = FALSE
)
probe_stats <- merge(probe_stats, expanded, by = "Probe_ID", all = FALSE, sort = FALSE)
panel_stats <- probe_stats[probe_stats$Gene %in% panel8, , drop = FALSE]

select_probe <- function(gene, metric) {
  sub <- panel_stats[panel_stats$Gene == gene, , drop = FALSE]
  if (!nrow(sub)) stop("No annotated probe for ", gene)
  sub <- sub[order(-sub[[metric]], sub$Probe_ID), , drop = FALSE]
  sub$Probe_ID[[1]]
}

primary_probe <- setNames(vapply(panel8, select_probe, character(1), metric = "Mean_Log2_Expression"), panel8)
iqr_probe <- setNames(vapply(panel8, select_probe, character(1), metric = "IQR_Log2_Expression"), panel8)
probe_inventory <- do.call(rbind, lapply(panel8, function(gene) {
  sub <- panel_stats[panel_stats$Gene == gene, , drop = FALSE]
  data.frame(
    Gene = gene,
    All_Probe_IDs = paste(sort(unique(sub$Probe_ID)), collapse = ";"),
    Probe_Count = length(unique(sub$Probe_ID)),
    Primary_Highest_Mean_Probe = primary_probe[[gene]],
    Primary_Mean_Log2_Expression = probe_stats$Mean_Log2_Expression[match(primary_probe[[gene]], probe_stats$Probe_ID)],
    Primary_IQR = probe_stats$IQR_Log2_Expression[match(primary_probe[[gene]], probe_stats$Probe_ID)],
    Sensitivity_Highest_IQR_Probe = iqr_probe[[gene]],
    Sensitivity_IQR_Mean_Log2_Expression = probe_stats$Mean_Log2_Expression[match(iqr_probe[[gene]], probe_stats$Probe_ID)],
    Sensitivity_IQR = probe_stats$IQR_Log2_Expression[match(iqr_probe[[gene]], probe_stats$Probe_ID)],
    Detected_Mean_Above_5 = probe_stats$Mean_Log2_Expression[match(primary_probe[[gene]], probe_stats$Probe_ID)] > detection_cutoff,
    Detection_Cutoff_Mean_Log2 = detection_cutoff,
    stringsAsFactors = FALSE
  )
}))
write.csv(probe_inventory, file.path(out_dir, "probe_selection_and_detection.csv"), row.names = FALSE)
all_probe_inventory <- merge(
  panel_stats[, c("Probe_ID", "Gene", "Mean_Log2_Expression", "IQR_Log2_Expression"), drop = FALSE],
  data.frame(
    Gene = panel8,
    Primary_Highest_Mean_Probe = unname(primary_probe),
    Sensitivity_Highest_IQR_Probe = unname(iqr_probe),
    stringsAsFactors = FALSE
  ),
  by = "Gene",
  all.x = TRUE,
  sort = FALSE
)
all_probe_inventory$Detected_Mean_Above_5 <- all_probe_inventory$Mean_Log2_Expression > detection_cutoff
write.csv(all_probe_inventory, file.path(out_dir, "all_panel_probe_inventory.csv"), row.names = FALSE)

group <- factor(group_primary, levels = c("Non-EPS", "EPS"))
design <- model.matrix(~0 + group)
colnames(design) <- c("Control", "EPS")
contrast <- limma::makeContrasts(EPS - Control, levels = design)
fit <- limma::lmFit(expr, design)
fit <- limma::contrasts.fit(fit, contrast)
fit <- limma::eBayes(fit)
limma_all <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none", confint = TRUE)
limma_all$Probe_ID <- rownames(limma_all)
write.csv(limma_all, file.path(out_dir, "limma_all_probe_results.csv"), row.names = FALSE)

discovery_path <- file.path("results", "tables", "GSE125498_all_results.csv")
discovery <- read.csv(discovery_path, check.names = FALSE)
if (!all(c("Gene", "logFC", "P.Value") %in% names(discovery))) {
  stop("Discovery results CSV lacks Gene, logFC, or P.Value.")
}
if (anyDuplicated(discovery$Gene)) stop("Discovery results CSV contains duplicate gene symbols.")
if (!all(panel8 %in% discovery$Gene)) stop("A fixed validation gene is absent from discovery results CSV.")
disc8 <- discovery[match(panel8, discovery$Gene), c("Gene", "logFC", "P.Value")]
if (any(!is.finite(disc8$logFC)) || any(disc8$logFC <= 0)) {
  stop("Discovery results CSV does not show all preregistered genes upregulated.")
}
names(disc8)[names(disc8) == "logFC"] <- "Discovery_log2FC"
names(disc8)[names(disc8) == "P.Value"] <- "Discovery_P"

exact_wilcox_p <- function(x, y) {
  if (length(x) != 4L || length(y) != 4L) stop("Exact rank test expects groups of four.")
  observed <- sum(outer(x, y, function(a, b) as.numeric(a > b) + 0.5 * as.numeric(a == b)))
  combos <- utils::combn(seq_len(8L), 4L)
  pooled <- c(x, y)
  null_u <- apply(combos, 2, function(idx) {
    a <- pooled[idx]
    b <- pooled[-idx]
    sum(outer(a, b, function(v, w) as.numeric(v > w) + 0.5 * as.numeric(v == w)))
  })
  mean(abs(null_u - 8) >= abs(observed - 8) - 1e-12)
}

cohen_d <- function(x, y) {
  pooled <- sqrt(((length(x) - 1) * stats::var(x) + (length(y) - 1) * stats::var(y)) /
                   (length(x) + length(y) - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (mean(x) - mean(y)) / pooled
}

bootstrap_d_ci <- function(x, y, reps = bootstrap_reps) {
  set.seed(42)
  vals <- replicate(reps, {
    bx <- sample(x, length(x), replace = TRUE)
    by <- sample(y, length(y), replace = TRUE)
    cohen_d(bx, by)
  })
  unname(stats::quantile(vals[is.finite(vals)], c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

primary_rows <- match(primary_probe, limma_all$Probe_ID)
if (anyNA(primary_rows)) stop("Primary selected probe missing from limma output.")
primary_results <- limma_all[primary_rows, , drop = FALSE]
primary_results$Gene <- panel8
names(primary_results)[names(primary_results) == "logFC"] <- "Validation_log2FC"
names(primary_results)[names(primary_results) == "CI.L"] <- "Validation_CI95_Lower"
names(primary_results)[names(primary_results) == "CI.R"] <- "Validation_CI95_Upper"
names(primary_results)[names(primary_results) == "P.Value"] <- "Validation_Limma_P"
names(primary_results)[names(primary_results) == "adj.P.Val"] <- "Validation_BH_P"
primary_results$Probe_ID <- unname(primary_probe)
primary_results <- merge(primary_results, disc8, by = "Gene", sort = FALSE)
primary_results <- primary_results[match(panel8, primary_results$Gene), , drop = FALSE]
primary_results$Discovery_expected_direction <- unname(expected_discovery_direction[primary_results$Gene])
primary_results$Direction_Concordant <- primary_results$Validation_log2FC > 0
primary_results$Mann_Whitney_exact_two_sided_P <- vapply(seq_along(panel8), function(i) {
  vals <- expr[primary_probe[[panel8[i]]], ]
  exact_wilcox_p(vals[group_primary == "EPS"], vals[group_primary == "Non-EPS"])
}, numeric(1))
primary_results$Cohens_d <- NA_real_
primary_results$Cohens_d_bootstrap_CI95_Lower <- NA_real_
primary_results$Cohens_d_bootstrap_CI95_Upper <- NA_real_
for (i in seq_along(panel8)) {
  vals <- expr[primary_probe[[panel8[i]]], ]
  eps <- vals[group_primary == "EPS"]
  ctrl <- vals[group_primary == "Non-EPS"]
  primary_results$Cohens_d[i] <- cohen_d(eps, ctrl)
  ci <- bootstrap_d_ci(eps, ctrl)
  primary_results$Cohens_d_bootstrap_CI95_Lower[i] <- ci[1]
  primary_results$Cohens_d_bootstrap_CI95_Upper[i] <- ci[2]
}
primary_results$Detection_Mean_Above_5 <- probe_inventory$Detected_Mean_Above_5[match(primary_results$Gene, probe_inventory$Gene)]
primary_results$Limma_BH_8_Genes <- p.adjust(primary_results$Validation_Limma_P, method = "BH")
write.csv(primary_results, file.path(out_dir, "per_gene_primary_results.csv"), row.names = FALSE)

selected_probe_analysis <- function(genes, probes, name) {
  rows <- match(probes[genes], limma_all$Probe_ID)
  out <- limma_all[rows, , drop = FALSE]
  out$Gene <- genes
  out$Probe_ID <- unname(probes[genes])
  stats <- lapply(seq_along(genes), function(i) {
    vals <- expr[probes[[genes[i]]], ]
    data.frame(
      Gene = genes[i],
      Mann_Whitney_exact_two_sided_P = exact_wilcox_p(vals[group_primary == "EPS"], vals[group_primary == "Non-EPS"]),
      stringsAsFactors = FALSE
    )
  })
  out$Mann_Whitney_exact_two_sided_P <- vapply(stats, `[[`, numeric(1), "Mann_Whitney_exact_two_sided_P")
  out$Probe_Selection = name
  out
}
sensitivity_iqr <- selected_probe_analysis(panel8, iqr_probe, "Highest_IQR")
write.csv(sensitivity_iqr, file.path(out_dir, "sensitivity_highest_iqr_results.csv"), row.names = FALSE)
all_probe_results <- merge(
  all_probe_inventory,
  limma_all,
  by.x = "Probe_ID",
  by.y = "Probe_ID",
  all.x = TRUE,
  sort = FALSE
)
all_probe_results$Mann_Whitney_exact_two_sided_P <- vapply(seq_len(nrow(all_probe_results)), function(i) {
  vals <- expr[all_probe_results$Probe_ID[i], ]
  exact_wilcox_p(vals[group_primary == "EPS"], vals[group_primary == "Non-EPS"])
}, numeric(1))
write.csv(all_probe_results, file.path(out_dir, "sensitivity_all_probes_results.csv"), row.names = FALSE)

gene_values <- function(genes, probes) {
  mat <- do.call(rbind, lapply(genes, function(gene) expr[probes[[gene]], sample_ids]))
  rownames(mat) <- genes
  mat
}

auc_rank <- function(scores, positive) {
  ranks <- rank(scores, ties.method = "average")
  n1 <- sum(positive)
  n0 <- length(positive) - n1
  (sum(ranks[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

composite_score <- function(genes, probes) {
  mat <- gene_values(genes, probes)
  z <- t(scale(t(mat), center = TRUE, scale = TRUE))
  list(matrix = mat, z = z, score = colMeans(z))
}

bootstrap_auc_ci <- function(scores, reps = bootstrap_reps) {
  positive <- which(group_primary == "EPS")
  negative <- which(group_primary == "Non-EPS")
  set.seed(42)
  vals <- replicate(reps, {
    idx <- c(sample(positive, length(positive), replace = TRUE),
             sample(negative, length(negative), replace = TRUE))
    auc_rank(scores[idx], group_primary[idx] == "EPS")
  })
  unname(stats::quantile(vals, c(0.025, 0.975), names = FALSE))
}

mean_score_difference <- function(scores) mean(scores[group_primary == "EPS"]) - mean(scores[group_primary == "Non-EPS"])
summarize_set <- function(genes, label) {
  z <- composite_score(genes, primary_probe)
  score <- z$score
  positive <- group_primary == "EPS"
  concordant <- sum(primary_results$Direction_Concordant[match(genes, primary_results$Gene)])
  auc <- auc_rank(score, positive)
  ci <- bootstrap_auc_ci(score)
  mw <- exact_wilcox_p(score[positive], score[!positive])
  binom <- stats::binom.test(concordant, length(genes), p = 0.5, alternative = "greater")
  list(
    genes = genes,
    label = label,
    score = score,
    matrix = z$matrix,
    z = z$z,
    concordant = concordant,
    n = length(genes),
    sign_test_one_sided_p = binom$p.value,
    score_difference = mean_score_difference(score),
    mann_whitney_exact_two_sided_p = mw,
    auc = auc,
    auc_ci = ci
  )
}

set8 <- summarize_set(panel8, "Eight-gene panel")
set2 <- summarize_set(strict2, "Strict 3/3 genes")
write.csv(
  data.frame(
    sample_id = sample_ids,
    subgroup = group3,
    primary_group = group_primary,
    score_8_gene = set8$score,
    score_2_strict = set2$score,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "composite_scores_by_sample.csv"),
  row.names = FALSE
)

gene_primary_all <- do.call(rbind, lapply(split(probe_stats, probe_stats$Gene), function(sub) {
  sub <- sub[order(-sub$Mean_Log2_Expression, sub$Probe_ID), , drop = FALSE]
  sub[1, c("Gene", "Probe_ID", "Mean_Log2_Expression"), drop = FALSE]
}))
gene_primary_all$Expression_Decile <- pmin(
  10L,
  ceiling(rank(gene_primary_all$Mean_Log2_Expression, ties.method = "average") /
            nrow(gene_primary_all) * 10)
)
detected_genes <- gene_primary_all[gene_primary_all$Mean_Log2_Expression > detection_cutoff, , drop = FALSE]
if (!all(panel8 %in% detected_genes$Gene)) {
  warning("At least one fixed panel gene is below the preregistered operational detection cutoff.")
}
detected_genes$Expression_Decile <- pmin(
  10L,
  ceiling(rank(detected_genes$Mean_Log2_Expression, ties.method = "average") /
            nrow(detected_genes) * 10)
)
matched_null <- function(genes, observed_auc, observed_mean_fc, label) {
  detected_ecdf <- stats::ecdf(detected_genes$Mean_Log2_Expression)
  hub_expression <- gene_primary_all$Mean_Log2_Expression[match(genes, gene_primary_all$Gene)]
  hub_bins <- pmax(1L, pmin(10L, ceiling(detected_ecdf(hub_expression) * 10)))
  pool <- detected_genes[!detected_genes$Gene %in% c(panel8, strict2), , drop = FALSE]
  need <- table(factor(hub_bins, levels = 1:10))
  available <- table(factor(pool$Expression_Decile, levels = 1:10))
  if (any(available < need)) {
    stop("Insufficient detected-gene candidates in exact expression decile(s) for ", label, " null.")
  }
  set.seed(42)
  aucs <- numeric(null_reps)
  mean_fcs <- numeric(null_reps)
  for (b in seq_len(null_reps)) {
    draws <- unlist(lapply(1:10, function(decile) {
      n <- as.integer(need[decile])
      if (!n) return(character())
      candidates <- pool$Gene[pool$Expression_Decile == decile]
      candidates[sample.int(length(candidates), n, replace = FALSE)]
    }), use.names = FALSE)
    selected <- gene_primary_all[match(draws, gene_primary_all$Gene), ]
    mat <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
      expr[selected$Probe_ID[i], sample_ids]
    }))
    z <- t(scale(t(mat), center = TRUE, scale = TRUE))
    score <- colMeans(z)
    aucs[b] <- auc_rank(score, group_primary == "EPS")
    mean_fcs[b] <- mean(rowMeans(mat[, group_primary == "EPS", drop = FALSE]) -
                          rowMeans(mat[, group_primary == "Non-EPS", drop = FALSE]))
  }
  p_auc <- (1 + sum(aucs >= observed_auc)) / (null_reps + 1)
  p_fc <- (1 + sum(mean_fcs >= observed_mean_fc)) / (null_reps + 1)
  write.csv(
    data.frame(iteration = seq_len(null_reps), Composite_AUC = aucs, Mean_Gene_log2FC = mean_fcs),
    file.path(out_dir, paste0("random_gene_set_null_", if (length(genes) == 8L) "8gene" else "2strict", ".csv")),
    row.names = FALSE
  )
  list(
    label = label,
    n_sets = null_reps,
    mean_auc = mean(aucs),
    sd_auc = stats::sd(aucs),
    auc_empirical_p = p_auc,
    n_null_auc_ge_observed = sum(aucs >= observed_auc),
    mean_fc = mean(mean_fcs),
    sd_fc = stats::sd(mean_fcs),
    mean_fc_empirical_p = p_fc,
    n_null_fc_ge_observed = sum(mean_fcs >= observed_mean_fc),
    observed_auc = observed_auc,
    observed_mean_gene_log2fc = observed_mean_fc,
    hub_expression_deciles = as.list(stats::setNames(hub_bins, genes)),
    hub_mean_expression = as.list(stats::setNames(hub_expression, genes)),
    minimum_detected_gene_mean = min(detected_genes$Mean_Log2_Expression),
    maximum_detected_gene_mean = max(detected_genes$Mean_Log2_Expression),
    decile_reference = "Empirical expression distribution among detected genes; below-range hubs assigned to the lowest decile"
  )
}

observed_fc8 <- mean(primary_results$Validation_log2FC)
observed_fc2 <- mean(primary_results$Validation_log2FC[match(strict2, primary_results$Gene)])
null8 <- matched_null(panel8, set8$auc, observed_fc8, "eight-gene panel")
null2 <- matched_null(strict2, set2$auc, observed_fc2, "strict two-gene panel")

marker_source <- readLines(file.path("scripts", "task4_marker_sets.py"), warn = FALSE)
marker_text <- paste(marker_source, collapse = "\n")
marker_match <- regexec("LEUKOCYTE_MARKERS\\s*=\\s*\\(([^)]*)\\)", marker_text)
marker_parts <- regmatches(marker_text, marker_match)[[1]]
if (length(marker_parts) < 2L) stop("Could not find LEUKOCYTE_MARKERS in Task 4 source.")
markers <- unique(trimws(gsub('"', "", unlist(strsplit(marker_parts[2], "[,\n]")))))
markers <- markers[nzchar(markers) & !markers %in% panel8]
mapped_markers <- intersect(markers, gene_primary_all$Gene)
score_markers <- mapped_markers
if (length(score_markers) < 2L) stop("Fewer than two leukocyte markers map to the validation platform.")
marker_probes <- setNames(
  gene_primary_all$Probe_ID[match(score_markers, gene_primary_all$Gene)],
  score_markers
)
marker_values <- gene_values(score_markers, marker_probes)
marker_z <- t(scale(t(marker_values), center = TRUE, scale = TRUE))
marker_score <- colMeans(marker_z)
marker_test <- exact_wilcox_p(marker_score[group_primary == "EPS"], marker_score[group_primary == "Non-EPS"])
marker_correlation <- stats::cor.test(marker_score, set8$score, method = "pearson")
write.csv(
  data.frame(sample_id = sample_ids, subgroup = group3, marker_score_proxy = marker_score,
             hub_composite_score = set8$score),
  file.path(out_dir, "leukocyte_marker_score_proxy.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Gene = score_markers, Probe_ID = unname(marker_probes), stringsAsFactors = FALSE),
  file.path(out_dir, "leukocyte_marker_genes_used.csv"),
  row.names = FALSE
)

subgroup_long <- do.call(rbind, lapply(panel8, function(gene) {
  probe <- primary_probe[[gene]]
  data.frame(
    Gene = gene,
    Probe_ID = probe,
    sample_id = sample_ids,
    subgroup = group3,
    Expression = as.numeric(expr[probe, ]),
    stringsAsFactors = FALSE
  )
}))
write.csv(subgroup_long, file.path(out_dir, "subgroup_descriptive_expression.csv"), row.names = FALSE)

power_search <- function() {
  ncp_fun <- function(d) {
    df <- 6
    ncp <- d * sqrt(4 * 4 / (4 + 4))
    critical <- stats::qt(0.975, df)
    stats::pt(-critical, df, ncp = ncp) + stats::pt(critical, df, ncp = ncp, lower.tail = FALSE) - 0.8
  }
  uniroot(ncp_fun, c(0, 10))$root
}
minimum_detectable_d <- power_search()

discovery_expr_path <- file.path("results", "tables", "GSE125498_full_expression_matrix.csv")
discovery_meta_path <- file.path("results", "tables", "GSE125498_sample_metadata.csv")
disc_expr_df <- read.csv(discovery_expr_path, check.names = FALSE)
disc_meta <- read.csv(discovery_meta_path, check.names = FALSE)
disc_expr <- as.matrix(disc_expr_df[, setdiff(names(disc_expr_df), "Gene"), drop = FALSE])
storage.mode(disc_expr) <- "double"
rownames(disc_expr) <- disc_expr_df$Gene
disc_ids <- as.character(disc_meta$sample_id)
if (!all(panel8 %in% rownames(disc_expr)) || !all(disc_ids %in% colnames(disc_expr))) {
  stop("Discovery expression table is incomplete for the fixed genes or samples.")
}
disc_group <- as.character(disc_meta$stage_label)
if (!all(c("Late_Stage_LPD", "Early_Stage_SPD") %in% disc_group)) {
  stop("Discovery metadata lacks the preregistered LPD/SPD groups.")
}
disc_fc_ci <- do.call(rbind, lapply(panel8, function(gene) {
  vals <- disc_expr[gene, disc_ids]
  a <- vals[disc_group == "Late_Stage_LPD"]
  b <- vals[disc_group == "Early_Stage_SPD"]
  test <- stats::t.test(a, b, var.equal = FALSE)
  data.frame(
    Gene = gene,
    Discovery_Welch_CI95_Lower = unname(test$conf.int[1]),
    Discovery_Welch_CI95_Upper = unname(test$conf.int[2]),
    stringsAsFactors = FALSE
  )
}))
forest <- merge(primary_results, disc_fc_ci, by = "Gene", sort = FALSE)
forest <- forest[match(panel8, forest$Gene), ]
write.csv(forest, file.path(out_dir, "forest_plot_data.csv"), row.names = FALSE)

score_summary <- data.frame(
  Set = c("Eight-gene panel", "Strict 3/3 genes"),
  Gene_Count = c(length(panel8), length(strict2)),
  Concordant_Genes = c(set8$concordant, set2$concordant),
  Sign_Test_One_Sided_P = c(set8$sign_test_one_sided_p, set2$sign_test_one_sided_p),
  Composite_Mean_Difference_EPS_minus_NonEPS = c(set8$score_difference, set2$score_difference),
  Composite_Mann_Whitney_Exact_Two_Sided_P = c(set8$mann_whitney_exact_two_sided_p, set2$mann_whitney_exact_two_sided_p),
  Composite_AUC = c(set8$auc, set2$auc),
  Composite_AUC_Bootstrap_CI95_Lower = c(set8$auc_ci[1], set2$auc_ci[1]),
  Composite_AUC_Bootstrap_CI95_Upper = c(set8$auc_ci[2], set2$auc_ci[2]),
  Mean_Gene_log2FC = c(observed_fc8, observed_fc2),
  Null_AUC_Empirical_P = c(null8$auc_empirical_p, null2$auc_empirical_p),
  Null_Mean_Gene_log2FC_Empirical_P = c(null8$mean_fc_empirical_p, null2$mean_fc_empirical_p),
  stringsAsFactors = FALSE
)
write.csv(score_summary, file.path(out_dir, "gene_set_level_results.csv"), row.names = FALSE)

null_summary <- data.frame(
  Set = c("Eight-gene panel", "Strict 3/3 genes"),
  Null_Sets = c(null8$n_sets, null2$n_sets),
  Observed_AUC = c(null8$observed_auc, null2$observed_auc),
  Null_AUC_Mean = c(null8$mean_auc, null2$mean_auc),
  Null_AUC_SD = c(null8$sd_auc, null2$sd_auc),
  Null_AUC_Count_Greater_Or_Equal = c(null8$n_null_auc_ge_observed, null2$n_null_auc_ge_observed),
  Null_AUC_Empirical_P = c(null8$auc_empirical_p, null2$auc_empirical_p),
  Observed_Mean_Gene_log2FC = c(null8$observed_mean_gene_log2fc, null2$observed_mean_gene_log2fc),
  Null_Mean_Gene_log2FC_Mean = c(null8$mean_fc, null2$mean_fc),
  Null_Mean_Gene_log2FC_SD = c(null8$sd_fc, null2$sd_fc),
  Null_log2FC_Count_Greater_Or_Equal = c(null8$n_null_fc_ge_observed, null2$n_null_fc_ge_observed),
  Null_log2FC_Empirical_P = c(null8$mean_fc_empirical_p, null2$mean_fc_empirical_p),
  stringsAsFactors = FALSE
)
write.csv(null_summary, file.path(out_dir, "random_gene_set_null_summary.csv"), row.names = FALSE)

# QC: boxplots use all RMA-normalized probe sets; PCA uses the 1,000 most variable
# probes, selected without labels. No sample is removed by the QC procedure.
grDevices::pdf(file.path(out_dir, "qc_expression_boxplot.pdf"), width = 10, height = 6)
graphics::boxplot(as.data.frame(expr), las = 2, outline = FALSE, ylab = "RMA log2 expression",
                  main = "GSE62928 per-sample expression distribution")
grDevices::dev.off()
grDevices::png(file.path(out_dir, "qc_expression_boxplot.png"), width = 1400, height = 850, res = 140)
graphics::boxplot(as.data.frame(expr), las = 2, outline = FALSE, ylab = "RMA log2 expression",
                  main = "GSE62928 per-sample expression distribution")
grDevices::dev.off()

probe_iqr <- apply(expr, 1, IQR)
top_pca_probes <- names(sort(probe_iqr, decreasing = TRUE))[seq_len(min(1000L, length(probe_iqr)))]
pca <- stats::prcomp(t(expr[top_pca_probes, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
pca_scores$sample_id <- rownames(pca_scores)
pca_scores$subgroup <- group3[match(pca_scores$sample_id, sample_ids)]
pca_sd <- apply(pca_scores[, c("PC1", "PC2")], 2, stats::sd)
pca_scores$PC1_abs_gt_2SD <- abs(pca_scores$PC1) > 2 * pca_sd["PC1"]
pca_scores$PC2_abs_gt_2SD <- abs(pca_scores$PC2) > 2 * pca_sd["PC2"]
write.csv(pca_scores, file.path(out_dir, "qc_pca_scores.csv"), row.names = FALSE)
variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)
grDevices::pdf(file.path(out_dir, "qc_pca_by_subgroup.pdf"), width = 8, height = 6)
graphics::plot(pca_scores$PC1, pca_scores$PC2, pch = 19,
               col = as.integer(factor(pca_scores$subgroup)),
               xlab = sprintf("PC1 (%.1f%%)", 100 * variance_explained[1]),
               ylab = sprintf("PC2 (%.1f%%)", 100 * variance_explained[2]),
               main = "GSE62928 PCA (top 1,000 variable probes)")
graphics::text(pca_scores$PC1, pca_scores$PC2, labels = pca_scores$sample_id, pos = 3, cex = 0.65)
graphics::legend("topright", legend = levels(factor(pca_scores$subgroup)),
                 col = seq_along(levels(factor(pca_scores$subgroup))), pch = 19, bty = "n")
grDevices::dev.off()

grDevices::pdf(file.path(out_dir, "gene_expression_primary_probe.pdf"), width = 12, height = 9)
graphics::par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
for (gene in panel8) {
  vals <- expr[primary_probe[[gene]], ]
  graphics::boxplot(vals ~ factor(group_primary, levels = c("Non-EPS", "EPS")),
                    outline = FALSE, ylab = "RMA log2 expression", xlab = "",
                    main = gene, col = c("#8DB4E2", "#D88080"))
  graphics::stripchart(vals ~ factor(group_primary, levels = c("Non-EPS", "EPS")),
                       vertical = TRUE, method = "jitter", add = TRUE, pch = 19,
                       col = "#333333", cex = 0.8)
}
grDevices::dev.off()

grDevices::pdf(file.path(out_dir, "gene_expression_subgroups_descriptive.pdf"), width = 12, height = 9)
graphics::par(mfrow = c(2, 4), mar = c(5, 4, 2.5, 1))
subgroup_levels <- c("EPS", "Early PD", "Uremic")
for (gene in panel8) {
  vals <- expr[primary_probe[[gene]], ]
  graphics::boxplot(vals ~ factor(group3, levels = subgroup_levels), outline = FALSE,
                    ylab = "RMA log2 expression", xlab = "", main = gene,
                    col = c("#D88080", "#8DB4E2", "#A7C7A1"))
  graphics::stripchart(vals ~ factor(group3, levels = subgroup_levels),
                       vertical = TRUE, method = "jitter", add = TRUE, pch = 19,
                       col = "#333333", cex = 0.8)
}
grDevices::dev.off()

grDevices::pdf(file.path(out_dir, "forest_discovery_validation.pdf"), width = 10, height = 7)
graphics::par(mar = c(5, 8, 2, 2))
ypos <- rev(seq_along(panel8))
xrange <- range(c(forest$Discovery_Welch_CI95_Lower, forest$Discovery_Welch_CI95_Upper,
                  forest$Validation_CI95_Lower, forest$Validation_CI95_Upper), finite = TRUE)
graphics::plot(NA, xlim = xrange, ylim = c(0.5, length(panel8) + 0.5),
               yaxt = "n", ylab = "", xlab = "log2 fold change (higher in EPS / LPD)",
               main = "Discovery and external validation effects")
graphics::axis(2, at = ypos, labels = panel8, las = 1)
graphics::abline(v = 0, lty = 2, col = "grey50")
for (i in seq_along(panel8)) {
  row <- forest[forest$Gene == panel8[i], ]
  y <- ypos[i]
  graphics::segments(row$Discovery_Welch_CI95_Lower, y + 0.10, row$Discovery_Welch_CI95_Upper, y + 0.10,
                     col = "#1F77B4", lwd = 2)
  graphics::points(row$Discovery_log2FC, y + 0.10, pch = 16, col = "#1F77B4")
  graphics::segments(row$Validation_CI95_Lower, y - 0.10, row$Validation_CI95_Upper, y - 0.10,
                     col = "#D62728", lwd = 2)
  graphics::points(row$Validation_log2FC, y - 0.10, pch = 17, col = "#D62728")
}
graphics::legend("bottomright", legend = c("Discovery estimate (Welch 95% CI)", "Validation limma estimate (moderated 95% CI)"),
                 col = c("#1F77B4", "#D62728"), pch = c(16, 17), lty = 1, bty = "n", cex = 0.8)
grDevices::dev.off()

roc_points <- function(scores, positive) {
  thresholds <- sort(unique(c(-Inf, scores, Inf)), decreasing = TRUE)
  tpr <- vapply(thresholds, function(t) mean(scores[positive] >= t), numeric(1))
  fpr <- vapply(thresholds, function(t) mean(scores[!positive] >= t), numeric(1))
  ord <- order(fpr, tpr)
  data.frame(FPR = fpr[ord], TPR = tpr[ord])
}
grDevices::pdf(file.path(out_dir, "composite_roc_and_null.pdf"), width = 11, height = 5.5)
graphics::par(mfrow = c(1, 2))
roc <- roc_points(set8$score, group_primary == "EPS")
graphics::plot(roc$FPR, roc$TPR, type = "l", lwd = 2, col = "#7C3AED",
               xlim = c(0, 1), ylim = c(0, 1), xlab = "False positive rate",
               ylab = "True positive rate",
               main = sprintf("Eight-gene composite ROC (AUC %.3f)", set8$auc))
graphics::abline(0, 1, lty = 2, col = "grey50")
null_auc8 <- read.csv(file.path(out_dir, "random_gene_set_null_8gene.csv"))$Composite_AUC
graphics::hist(null_auc8, breaks = seq(0, 1, by = 0.05), col = "grey80", border = "white",
               main = sprintf("Matched 8-gene null (empirical P %.4f)", null8$auc_empirical_p),
               xlab = "Null composite AUC")
graphics::abline(v = set8$auc, col = "#D62728", lwd = 2)
grDevices::dev.off()

qc_summary <- list(
  dataset = "GSE62928",
  platform = "GPL13158 Affymetrix HG-U133 Plus 2.0",
  sample_count = length(sample_ids),
  sample_groups = as.list(observed_counts),
  data_processing = processing_text,
  matrix_logged_rma = TRUE,
  normalization_applied_again = FALSE,
  probe_sets_in_matrix = nrow(expr),
  mapped_unique_gene_symbols = length(unique(expanded$Gene)),
  primary_probe_detection_cutoff_mean_log2_gt = detection_cutoff,
  primary_panel_genes_detected = sum(probe_inventory$Detected_Mean_Above_5),
  panel_genes_without_probe = missing_genes,
  pca_variable_probe_count = length(top_pca_probes),
  pca_descriptive_outlier_rule = "Flag if absolute PC1 or PC2 score exceeds 2 sample SD; no samples removed.",
  pca_flagged_samples = pca_scores$sample_id[pca_scores$PC1_abs_gt_2SD | pca_scores$PC2_abs_gt_2SD],
  sample_ids_removed = character(),
  random_seed = 42
)
summary <- list(
  qc = qc_summary,
  discovery_values_source = normalizePath(discovery_path, winslash = "/", mustWork = TRUE),
  preregistered_genes = panel8,
  strict_3_of_3_genes = strict2,
  per_gene_concordant_count = sum(primary_results$Direction_Concordant),
  per_gene_total = length(panel8),
  gene_set_results = score_summary,
  null_results = null_summary,
  marker_score = list(
    label = "Leukocyte marker-score proxy, not cell-fraction deconvolution",
    source = "scripts/task4_marker_sets.py: LEUKOCYTE_MARKERS",
    genes_requested_after_hub_removal = markers,
    genes_used = score_markers,
    genes_used_n = length(score_markers),
    eps_mean = mean(marker_score[group_primary == "EPS"]),
    control_mean = mean(marker_score[group_primary == "Non-EPS"]),
    eps_minus_control = mean(marker_score[group_primary == "EPS"]) - mean(marker_score[group_primary == "Non-EPS"]),
    mann_whitney_exact_two_sided_p = marker_test,
    pearson_r_with_8_gene_composite = unname(marker_correlation$estimate),
    pearson_p_with_8_gene_composite = marker_correlation$p.value
  ),
  power = list(
    model = "Two-sample independent t-test, equal group sizes, two-sided alpha=0.05, 80% power, df=6",
    n_per_group = 4,
    alpha = 0.05,
    target_power = 0.80,
    minimum_detectable_cohens_d = minimum_detectable_d
  ),
  interpretation_rule = list(
    supportive = "At least 6/8 concordant AND composite AUC >= 0.80 AND empirical matched-null AUC P < 0.05",
    partially_supportive = "5/8 concordant OR AUC 0.65-0.80 with empirical matched-null AUC P >= 0.05",
    not_supported = "Otherwise"
  ),
  software = list(R = as.character(getRversion()), limma = as.character(utils::packageVersion("limma"))),
  bootstrap_replicates = bootstrap_reps,
  random_gene_set_replicates = null_reps,
  random_seed = 42
)
summary_path <- file.path(out_dir, "validation_summary.json")
jsonlite::write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, na = "null")

cat("GSE62928 real-data analysis complete.\n")
cat("Samples:", length(sample_ids), "; groups EPS/early PD/uremic:",
    paste(as.integer(observed_counts), collapse = "/"), "\n")
cat("Probe sets:", nrow(expr), "; unique mapped genes:", length(unique(expanded$Gene)), "\n")
cat("Panel mapped:", length(panel8) - length(missing_genes), "/8; detected above mean log2", detection_cutoff, ":",
    sum(probe_inventory$Detected_Mean_Above_5), "\n")
cat("Eight-gene concordance:", set8$concordant, "/", set8$n,
    "; AUC:", signif(set8$auc, 5), "; null empirical P:", signif(null8$auc_empirical_p, 5), "\n")
cat("Strict two-gene AUC:", signif(set2$auc, 5), "; null empirical P:", signif(null2$auc_empirical_p, 5), "\n")
cat("Minimum detectable Cohen's d:", signif(minimum_detectable_d, 5), "\n")
