options(stringsAsFactors = FALSE, scipen = 999)
set.seed(2026)

project_root <- function() {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg)) return(normalizePath(file.path(dirname(sub("^--file=", "", arg[[1]])), ".."), winslash = "/", mustWork = TRUE))
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- project_root()
RESULTS <- file.path(ROOT, "results")
FIGURES <- file.path(RESULTS, "figures")
TABLES <- file.path(RESULTS, "tables")
LOGS <- file.path(RESULTS, "logs")
for (path in c(FIGURES, TABLES, LOGS)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

HUBS <- c("TNFSF15", "FLT3LG", "EBI3", "LTB", "ADAM19", "NAT6", "ECM1", "CST7", "CCBE1")
GROUP_COLORS <- c(Control = "#0072B2", Case = "#D55E00")
STATUS_COLORS <- c(Upregulated = "#D55E00", Downregulated = "#0072B2")
GENE_COLORS <- c(
  TNFSF15 = "#0072B2", EBI3 = "#E69F00", FLT3LG = "#009E73",
  LTB = "#D55E00", ADAM19 = "#CC79A7", ECM1 = "#56B4E9",
  CST7 = "#000000", NAT6 = "#F0E442", CCBE1 = "#882255"
)

require_pkgs <- function(pkgs, task) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop(task, ": missing R packages: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

task_log <- function(task, lines, status = "completed") {
  path <- file.path(LOGS, paste0(task, ".log"))
  writeLines(c(paste0("Task: ", task), paste0("Status: ", status), paste0("Seed: 2026"), "", lines), path)
  invisible(path)
}

append_results_text <- function(task, paragraph) {
  path <- file.path(RESULTS, "results_text.md")
  heading <- paste0("### ", task)
  text <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  starts <- grep("^### ", text)
  prefix <- if (length(starts) && starts[[1]] > 1L) text[seq_len(starts[[1]] - 1L)] else character()
  blocks <- lapply(seq_along(starts), function(i) {
    end <- if (i < length(starts)) starts[[i + 1L]] - 1L else length(text)
    text[starts[[i]]:end]
  })
  titles <- vapply(blocks, `[[`, character(1), 1L)
  blocks <- blocks[titles != heading]
  blocks[[length(blocks) + 1L]] <- c(heading, "", paragraph, "")
  titles <- vapply(blocks, `[[`, character(1), 1L)
  order_key <- suppressWarnings(as.integer(sub("^### T0*([0-9]+).*", "\\1", titles)))
  order_key[is.na(order_key)] <- Inf
  blocks <- blocks[order(order_key, titles)]
  writeLines(c(prefix, unlist(blocks, use.names = FALSE)), path)
}

save_plot <- function(plot, stem, width = 9, height = 6) {
  ggplot2::ggsave(file.path(FIGURES, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(FIGURES, paste0(stem, ".pdf")), plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
}

theme_project <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_rect(fill = "grey95", colour = NA))
}

read_csv_safe <- function(path, ...) {
  if (!file.exists(path)) stop("Required input file does not exist: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

load_discovery <- function() {
  expr_path <- file.path(TABLES, "GSE125498_full_expression_matrix.csv")
  meta_path <- file.path(TABLES, "GSE125498_sample_metadata.csv")
  expr_df <- read_csv_safe(expr_path)
  meta <- read_csv_safe(meta_path)
  if (!"Gene" %in% names(expr_df)) stop("Discovery expression matrix has no Gene column.")
  expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- as.character(expr_df$Gene)
  if (!all(meta$sample_id %in% colnames(expr))) stop("Discovery sample IDs are not all represented in the expression matrix.")
  expr <- expr[, meta$sample_id, drop = FALSE]
  meta$Group <- ifelse(meta$stage_binary == 1, "Case", "Control")
  meta$Group <- factor(meta$Group, levels = c("Control", "Case"))
  list(expr = expr, meta = meta, cohort = "GSE125498", scale = "RMA log2 microarray")
}

load_effluent <- function(cohort) {
  d <- file.path(RESULTS, "revision", "effluent_validation", cohort, "analysis_20261009")
  meta <- read_csv_safe(file.path(d, "all_cell_pseudobulk_metadata.csv"))
  meta <- meta[meta$primary_group %in% c("LV", "SV"), , drop = FALSE]
  meta <- meta[!duplicated(meta$sample_id), , drop = FALSE]
  meta$Group <- ifelse(meta$primary_group == "LV", "Case", "Control")
  meta$Group <- factor(meta$Group, levels = c("Control", "Case"))
  raw <- read_csv_safe(file.path(d, "edger_logCPM.csv"), row.names = 1)
  expr <- as.matrix(raw)
  storage.mode(expr) <- "double"
  rownames(expr) <- as.character(rownames(raw))
  mapping <- read_csv_safe(file.path(d, "gene_id_to_symbol.csv"))
  mapping <- mapping[!duplicated(mapping$gene_id), , drop = FALSE]
  symbols <- mapping$gene_symbol[match(rownames(expr), mapping$gene_id)]
  keep <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  expr <- rowsum(expr, group = symbols, reorder = FALSE) /
    as.numeric(table(factor(symbols, levels = unique(symbols))))
  if (!all(meta$sample_id %in% colnames(expr))) stop(cohort, " metadata samples are missing from logCPM.")
  expr <- expr[, meta$sample_id, drop = FALSE]
  list(expr = expr, meta = meta, cohort = cohort, scale = "edgeR TMM log2-CPM")
}

load_all_cohorts <- function() {
  list(load_discovery(), load_effluent("GSE248762"), load_effluent("GSE130888"))
}

auc_rank <- function(labels, scores) {
  ok <- is.finite(scores) & !is.na(labels)
  labels <- as.integer(labels[ok])
  scores <- scores[ok]
  n1 <- sum(labels == 1L)
  n0 <- sum(labels == 0L)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(scores, ties.method = "average")[labels == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

bootstrap_auc <- function(labels, scores, B = 2000L) {
  pos <- scores[labels == 1L]
  neg <- scores[labels == 0L]
  if (!length(pos) || !length(neg)) return(c(NA_real_, NA_real_))
  draws <- replicate(B, {
    p <- sample(pos, length(pos), replace = TRUE)
    n <- sample(neg, length(neg), replace = TRUE)
    auc_rank(c(rep(1L, length(p)), rep(0L, length(n))), c(p, n))
  })
  unname(stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

roc_points <- function(labels, scores) {
  ord <- order(scores, decreasing = TRUE)
  y <- as.integer(labels[ord])
  s <- scores[ord]
  thresholds <- unique(s)
  tpr <- fpr <- numeric(length(thresholds) + 1L)
  for (i in seq_along(thresholds)) {
    pred <- s >= thresholds[[i]]
    tpr[[i + 1L]] <- sum(pred & y == 1L) / sum(y == 1L)
    fpr[[i + 1L]] <- sum(pred & y == 0L) / sum(y == 0L)
  }
  data.frame(FPR = fpr, TPR = tpr)
}

hub_expr_long <- function(bundle) {
  genes <- intersect(HUBS, rownames(bundle$expr))
  if (!length(genes)) return(data.frame())
  x <- as.data.frame(t(bundle$expr[genes, , drop = FALSE]))
  x$sample_id <- rownames(x)
  x <- merge(x, bundle$meta[, c("sample_id", "Group")], by = "sample_id", sort = FALSE)
  out <- tidyr::pivot_longer(x, cols = all_of(genes), names_to = "Gene", values_to = "Expression")
  out$cohort <- bundle$cohort
  out
}

write_skip <- function(task, reason) {
  msg <- paste0("Skipped conservatively: ", reason)
  task_log(task, msg, status = "skipped")
  append_results_text(task, paste0(msg, " No result was inferred without the required input or method."))
  message(task, ": ", msg)
}
