#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(limma))
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript scripts/validation/validate_ml_8_consensus_hubs_gse62928.R ",
    "<series_matrix.gz> <GPL13158.soft.gz> <hub_genes.csv> <output_dir>"
  )
}

matrix_path <- args[[1]]
platform_path <- args[[2]]
panel_path <- args[[3]]
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_source <- read.csv(panel_path, check.names = FALSE)
if (!all(c("Gene", "Discovery_logFC", "Discovery_P_Value") %in% names(panel_source))) {
  stop("Hub CSV must contain Gene, Discovery_logFC, and Discovery_P_Value.")
}
panel <- unique(as.character(panel_source$Gene))
if (length(panel) != 8L || anyNA(panel) || any(!nzchar(panel)) ||
    anyDuplicated(panel)) {
  stop("Expected exactly eight unique, non-empty genes in the hub CSV.")
}
if (!all(panel_source$Discovery_logFC > 0)) {
  stop("This analysis expects all panel genes to be upregulated in discovery.")
}
discovery <- panel_source[match(panel, panel_source$Gene), , drop = FALSE]

read_geo_preamble <- function(path) {
  con <- gzfile(path, open = "rt")
  lines <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("GEO matrix ended before its expression table.")
    if (startsWith(line, "!series_matrix_table_begin")) break
    lines <- c(lines, line)
  }
  list(lines = lines, connection = con)
}

parse_geo_vector <- function(lines, prefix) {
  line <- lines[startsWith(lines, prefix)]
  if (length(line) != 1L) stop("Expected one GEO metadata row for ", prefix)
  trimws(gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]][-1]))
}

geo <- read_geo_preamble(matrix_path)
if (length(geo$lines[startsWith(geo$lines, "!Series_platform_id")]) != 1L ||
    !grepl("GPL13158", paste(geo$lines, collapse = "\n"), fixed = TRUE)) {
  stop("GEO matrix does not identify GPL13158.")
}
sample_ids <- parse_geo_vector(geo$lines, "!Sample_geo_accession")
sample_titles <- parse_geo_vector(geo$lines, "!Sample_title")
sample_desc <- parse_geo_vector(geo$lines, "!Sample_description")
char_lines <- geo$lines[startsWith(geo$lines, "!Sample_characteristics_ch1")]
if (!length(char_lines)) stop("GEO sample characteristics are missing.")
char_mat <- lapply(char_lines, function(line) {
  trimws(gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]][-1]))
})
if (any(lengths(char_mat) != length(sample_ids))) stop("Sample characteristics misaligned.")
char_text <- vapply(seq_along(sample_ids), function(i) {
  paste(vapply(char_mat, `[[`, character(1), i), collapse = "; ")
}, character(1))
subgroup <- ifelse(
  grepl("^EPS", sample_titles, ignore.case = TRUE), "EPS",
  ifelse(grepl("^PD", sample_titles, ignore.case = TRUE), "Early PD",
         ifelse(grepl("^UREMIC", sample_titles, ignore.case = TRUE), "Uremic", NA_character_))
)
if (anyNA(subgroup) || length(sample_ids) != 8L || anyDuplicated(sample_ids)) {
  stop("Expected eight uniquely identified EPS, PD, and uremic samples.")
}
if (any(subgroup == "EPS" & !grepl("treatment: EPS", char_text, ignore.case = TRUE)) ||
    any(subgroup == "Early PD" & !grepl("treatment: PD", char_text, ignore.case = TRUE)) ||
    any(subgroup == "Uremic" & !grepl("UREMIC", char_text, ignore.case = TRUE))) {
  stop("GEO sample title and characteristics disagree.")
}
if (!identical(as.integer(table(factor(subgroup, levels = c("EPS", "Early PD", "Uremic")))),
               c(4L, 2L, 2L))) {
  stop("Expected four EPS, two PD, and two uremic samples.")
}
primary_group <- ifelse(subgroup == "EPS", "EPS", "Non-EPS")
processing <- paste(geo$lines[startsWith(geo$lines, "!Sample_data_processing")], collapse = " ")
if (!grepl("RMA", processing, ignore.case = TRUE)) {
  stop("GEO metadata does not confirm RMA normalization.")
}

expr_df <- read.delim(
  geo$connection, header = TRUE, check.names = FALSE, quote = "\"",
  comment.char = "", stringsAsFactors = FALSE
)
close(geo$connection)
expr_df <- expr_df[!grepl("^!", as.character(expr_df[[1]])), , drop = FALSE]
probe_ids <- as.character(expr_df[[1]])
expr <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- probe_ids
if (!all(sample_ids %in% colnames(expr)) || any(!is.finite(expr))) {
  stop("Expression data do not match GEO metadata or contain non-finite values.")
}
expr <- expr[, sample_ids, drop = FALSE]

soft_con <- gzfile(platform_path, open = "rt")
annotation_lines <- readLines(soft_con, warn = FALSE)
close(soft_con)
begin <- which(annotation_lines == "!platform_table_begin")
end <- which(annotation_lines == "!platform_table_end")
if (length(begin) != 1L || length(end) != 1L || end <= begin) {
  stop("Could not locate GPL13158 annotation table.")
}
annotation <- read.delim(
  text = paste(annotation_lines[(begin + 1L):(end - 1L)], collapse = "\n"),
  check.names = FALSE, quote = "", fill = TRUE, stringsAsFactors = FALSE
)
if (!all(c("ID", "Gene Symbol") %in% names(annotation))) {
  stop("GPL13158 annotation lacks ID and Gene Symbol.")
}
annotation <- annotation[annotation$ID %in% probe_ids, c("ID", "Gene Symbol"), drop = FALSE]
expanded <- do.call(rbind, lapply(seq_len(nrow(annotation)), function(i) {
  symbol <- annotation[["Gene Symbol"]][i]
  if (is.na(symbol) || !nzchar(trimws(symbol)) || symbol %in% c("---", "NA")) return(NULL)
  symbols <- trimws(strsplit(symbol, "///", fixed = TRUE)[[1]])
  symbols <- symbols[nzchar(symbols) & !symbols %in% c("---", "NA")]
  if (!length(symbols)) return(NULL)
  data.frame(Probe_ID = annotation$ID[i], Gene = symbols, stringsAsFactors = FALSE)
}))
if (is.null(expanded)) stop("No mapped genes in GPL13158 annotation.")
expanded <- unique(expanded)
panel_probes <- expanded[expanded$Gene %in% panel, , drop = FALSE]
if (!setequal(panel, unique(panel_probes$Gene))) {
  stop("Missing GPL13158 probe mapping for: ",
       paste(setdiff(panel, unique(panel_probes$Gene)), collapse = ", "))
}

probe_stats <- data.frame(
  Probe_ID = probe_ids,
  Mean_Log2_Expression = rowMeans(expr),
  stringsAsFactors = FALSE
)
probe_stats <- merge(probe_stats, expanded, by = "Probe_ID", all = FALSE, sort = FALSE)
select_probe <- function(gene) {
  rows <- probe_stats[probe_stats$Gene == gene, , drop = FALSE]
  rows <- rows[order(-rows$Mean_Log2_Expression, rows$Probe_ID), , drop = FALSE]
  if (!nrow(rows)) stop("No probe available for ", gene)
  rows[1, c("Probe_ID", "Mean_Log2_Expression"), drop = FALSE]
}
selected <- do.call(rbind, lapply(panel, select_probe))
selected$Gene <- panel
selected <- selected[, c("Gene", "Probe_ID", "Mean_Log2_Expression")]
selected$Detected_Mean_Above_5 <- selected$Mean_Log2_Expression > 5
write.csv(selected, file.path(out_dir, "probe_selection.csv"), row.names = FALSE)
primary_probe <- setNames(selected$Probe_ID, selected$Gene)

group <- factor(primary_group, levels = c("Non-EPS", "EPS"))
design <- model.matrix(~0 + group)
colnames(design) <- c("Control", "EPS")
fit <- limma::lmFit(expr, design)
fit <- limma::contrasts.fit(fit, limma::makeContrasts(EPS - Control, levels = design))
fit <- limma::eBayes(fit)
all_results <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none", confint = TRUE)
all_results$Probe_ID <- rownames(all_results)

exact_mw_p <- function(x, y) {
  observed <- sum(outer(x, y, function(a, b) as.numeric(a > b) + 0.5 * as.numeric(a == b)))
  combos <- utils::combn(seq_len(8L), 4L)
  pooled <- c(x, y)
  null_u <- apply(combos, 2, function(idx) {
    sum(outer(pooled[idx], pooled[-idx],
              function(a, b) as.numeric(a > b) + 0.5 * as.numeric(a == b)))
  })
  mean(abs(null_u - 8) >= abs(observed - 8) - 1e-12)
}

per_gene <- do.call(rbind, lapply(seq_along(panel), function(i) {
  gene <- panel[i]
  probe <- primary_probe[[gene]]
  result <- all_results[match(probe, all_results$Probe_ID), , drop = FALSE]
  values <- expr[probe, ]
  data.frame(
    Gene = gene,
    Discovery_logFC = discovery$Discovery_logFC[i],
    Discovery_P_Value = discovery$Discovery_P_Value[i],
    Probe_ID = probe,
    Mean_Log2_Expression = selected$Mean_Log2_Expression[i],
    Detection_Mean_Above_5 = selected$Detected_Mean_Above_5[i],
    Validation_logFC = result$logFC,
    Validation_CI95_Lower = result$CI.L,
    Validation_CI95_Upper = result$CI.R,
    Validation_Limma_P = result$P.Value,
    Validation_BH_P = NA_real_,
    Mann_Whitney_exact_two_sided_P = exact_mw_p(
      values[primary_group == "EPS"], values[primary_group == "Non-EPS"]
    ),
    Direction_Concordant = result$logFC > 0,
    stringsAsFactors = FALSE
  )
}))
per_gene$Validation_BH_P <- p.adjust(per_gene$Validation_Limma_P, method = "BH")
write.csv(per_gene, file.path(out_dir, "per_gene_validation.csv"), row.names = FALSE)

mat <- do.call(rbind, lapply(panel, function(gene) expr[primary_probe[[gene]], ]))
rownames(mat) <- panel
gene_z <- t(scale(t(mat), center = TRUE, scale = TRUE))
scores <- colMeans(gene_z)
auc_rank <- function(values, positive) {
  ranks <- rank(values, ties.method = "average")
  n1 <- sum(positive)
  n0 <- length(values) - n1
  (sum(ranks[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
positive <- primary_group == "EPS"
observed_auc <- auc_rank(scores, positive)
set.seed(42)
boot_auc <- replicate(10000L, {
  ids <- c(sample(which(positive), 4L, replace = TRUE),
           sample(which(!positive), 4L, replace = TRUE))
  auc_rank(scores[ids], positive[ids])
})
observed_fc <- mean(per_gene$Validation_logFC)

gene_primary <- do.call(rbind, lapply(split(probe_stats, probe_stats$Gene), function(rows) {
  rows <- rows[order(-rows$Mean_Log2_Expression, rows$Probe_ID), , drop = FALSE]
  rows[1, c("Gene", "Probe_ID", "Mean_Log2_Expression"), drop = FALSE]
}))
detected <- gene_primary[gene_primary$Mean_Log2_Expression > 5, , drop = FALSE]
if (!nrow(detected)) stop("No genes above the operational detection cutoff for null matching.")
detected$Expression_Decile <- pmin(
  10L, ceiling(rank(detected$Mean_Log2_Expression, ties.method = "average") / nrow(detected) * 10)
)
gene_primary$Expression_Decile <- NA_integer_
gene_primary$Expression_Decile[match(detected$Gene, gene_primary$Gene)] <- detected$Expression_Decile
detected_ecdf <- stats::ecdf(detected$Mean_Log2_Expression)
hub_expression <- gene_primary$Mean_Log2_Expression[match(panel, gene_primary$Gene)]
hub_bins <- pmax(1L, pmin(10L, ceiling(detected_ecdf(hub_expression) * 10)))
pool <- detected[!detected$Gene %in% panel, , drop = FALSE]
need <- table(factor(hub_bins, levels = 1:10))
available <- table(factor(pool$Expression_Decile, levels = 1:10))
if (any(available < need)) stop("Insufficient genes for expression-decile matched null.")

set.seed(42)
null_auc <- numeric(1000L)
null_fc <- numeric(1000L)
for (b in seq_len(1000L)) {
  draws <- unlist(lapply(1:10, function(decile) {
    n <- as.integer(need[decile])
    if (!n) return(character())
    candidates <- pool$Gene[pool$Expression_Decile == decile]
    candidates[sample.int(length(candidates), n, replace = FALSE)]
  }), use.names = FALSE)
  chosen <- gene_primary[match(draws, gene_primary$Gene), , drop = FALSE]
  null_mat <- do.call(rbind, lapply(seq_len(nrow(chosen)), function(i) {
    expr[chosen$Probe_ID[i], ]
  }))
  null_z <- t(scale(t(null_mat), center = TRUE, scale = TRUE))
  null_score <- colMeans(null_z)
  null_auc[b] <- auc_rank(null_score, positive)
  null_fc[b] <- mean(rowMeans(null_mat[, positive, drop = FALSE]) -
                       rowMeans(null_mat[, !positive, drop = FALSE]))
}
write.csv(
  data.frame(Iteration = seq_len(1000L), Null_AUC = null_auc, Null_Mean_Gene_logFC = null_fc),
  file.path(out_dir, "matched_gene_set_null.csv"), row.names = FALSE
)
panel_summary <- data.frame(
  Gene_Count = length(panel),
  Concordant_Genes = sum(per_gene$Direction_Concordant),
  Composite_Mean_Difference_EPS_minus_NonEPS = mean(scores[positive]) - mean(scores[!positive]),
  Composite_Mann_Whitney_exact_two_sided_P = exact_mw_p(scores[positive], scores[!positive]),
  Composite_AUC = observed_auc,
  Composite_AUC_Bootstrap_CI95_Lower = unname(quantile(boot_auc, 0.025)),
  Composite_AUC_Bootstrap_CI95_Upper = unname(quantile(boot_auc, 0.975)),
  Mean_Gene_logFC = observed_fc,
  Matched_Null_Mean_AUC = mean(null_auc),
  Matched_Null_AUC_Empirical_P = (1 + sum(null_auc >= observed_auc)) / 1001,
  Matched_Null_Mean_Gene_logFC_Empirical_P = (1 + sum(null_fc >= observed_fc)) / 1001,
  Detected_Panel_Genes = sum(selected$Detected_Mean_Above_5),
  stringsAsFactors = FALSE
)
write.csv(panel_summary, file.path(out_dir, "panel_validation_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(
    Sample_ID = sample_ids,
    GEO_Title = sample_titles,
    Subgroup = subgroup,
    Analysis_Group = primary_group,
    Composite_Z_Score = as.numeric(scores)
  ),
  file.path(out_dir, "sample_composite_scores.csv"),
  row.names = FALSE
)

message("Validated ", length(panel), " genes against GSE62928 (4 EPS vs 4 non-EPS).")
message("Direction-concordant genes: ", sum(per_gene$Direction_Concordant), "/", length(panel))
message("Composite AUC: ", signif(observed_auc, 4))
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
