#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript scripts/validation/validate_ml_8_hubs_effluent_cohorts.R ",
    "<GSE248762|GSE130888> <analysis_input_dir> <hub_genes.csv> <output_dir>"
  )
}

accession <- args[[1]]
input_dir <- args[[2]]
panel_path <- args[[3]]
out_dir <- args[[4]]
if (!accession %in% c("GSE248762", "GSE130888")) stop("Unsupported accession: ", accession)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_source <- read.csv(panel_path, check.names = FALSE)
required_panel_columns <- c("Gene", "Discovery_logFC", "Discovery_P_Value")
if (!all(required_panel_columns %in% names(panel_source))) {
  stop("Hub CSV must contain: ", paste(required_panel_columns, collapse = ", "))
}
panel <- as.character(panel_source$Gene)
if (length(panel) != 8L || anyNA(panel) || any(!nzchar(panel)) || anyDuplicated(panel)) {
  stop("Expected exactly eight unique, non-empty panel genes.")
}
panel_source <- panel_source[match(panel, panel_source$Gene), , drop = FALSE]
if (any(!is.finite(panel_source$Discovery_logFC))) stop("Discovery logFCs must be finite.")

counts_table <- read.csv(
  file.path(input_dir, "all_cell_pseudobulk_counts.csv"),
  row.names = 1L, check.names = FALSE
)
counts <- as.matrix(counts_table)
storage.mode(counts) <- "integer"
metadata <- read.csv(
  file.path(input_dir, "all_cell_pseudobulk_metadata.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
gene_map <- read.csv(
  file.path(input_dir, "gene_id_to_symbol.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
if (!all(colnames(counts) %in% metadata$sample_id)) {
  stop("Some pseudobulk columns have no matching sample metadata.")
}
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id) || any(counts < 0L)) {
  stop("Sample metadata mismatch or negative pseudobulk count.")
}
if (!all(c("gene_id", "gene_symbol") %in% names(gene_map))) {
  stop("Gene mapping must contain gene_id and gene_symbol.")
}

gene_map$gene_symbol <- toupper(trimws(gene_map$gene_symbol))
panel_map <- gene_map[gene_map$gene_symbol %in% panel, , drop = FALSE]
if (!setequal(panel, panel_map$gene_symbol)) {
  stop("Missing panel gene mapping(s): ",
       paste(setdiff(panel, panel_map$gene_symbol), collapse = ", "))
}
if (anyDuplicated(panel_map$gene_symbol) || anyDuplicated(panel_map$gene_id)) {
  stop("Panel symbols must map one-to-one to unique gene IDs.")
}
gene_ids <- setNames(panel_map$gene_id, panel_map$gene_symbol)[panel]
if (!all(gene_ids %in% rownames(counts))) stop("A panel gene ID is missing from count matrix.")

if (!"original_group" %in% names(metadata)) stop("Missing original_group metadata.")
metadata$validation_group <- ifelse(
  grepl("^LV", metadata$original_group), "LV",
  ifelse(metadata$original_group == "SV", "SV", NA_character_)
)
keep_samples <- !is.na(metadata$validation_group)
metadata <- metadata[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples, drop = FALSE]
metadata$validation_group <- factor(metadata$validation_group, levels = c("SV", "LV"))
group_sizes <- table(metadata$validation_group)
if (any(group_sizes < 2L)) stop("Both LV and SV must have at least two patient pseudobulks.")
if (anyDuplicated(metadata$sample_id)) stop("Patient pseudobulk sample IDs must be unique.")

detection_by_group <- lapply(levels(metadata$validation_group), function(level) {
  selected_samples <- metadata$validation_group == level
  rowSums(counts[gene_ids, selected_samples, drop = FALSE] >= 10L)
})
names(detection_by_group) <- levels(metadata$validation_group)
detected <- vapply(seq_along(panel), function(i) {
  all(vapply(names(detection_by_group), function(level) {
    detection_by_group[[level]][i] >= ceiling(sum(metadata$validation_group == level) / 2)
  }, logical(1)))
}, logical(1))

y <- DGEList(counts = counts)
design <- model.matrix(~0 + validation_group, data = metadata)
colnames(design) <- c("SV", "LV")
keep <- filterByExpr(y, design = design)
detected_for_fit <- vapply(levels(metadata$validation_group), function(level) {
  selected_samples <- metadata$validation_group == level
  rowSums(counts[, selected_samples, drop = FALSE] >= 10L) >=
    ceiling(sum(selected_samples) / 2)
}, logical(nrow(counts)))
keep <- keep | apply(detected_for_fit, 1L, all)
keep[gene_ids] <- TRUE
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
contrast <- makeContrasts(LV - SV, levels = design)
ql_result <- glmQLFTest(fit, contrast = contrast)
ql_table <- topTags(ql_result, n = Inf, sort.by = "none")$table

ql_variance <- fit$s2.post
dispersion <- fit$dispersion
residual_df <- fit$df.residual.adj
if (length(ql_variance) == 1L) ql_variance <- rep(ql_variance, nrow(ql_table))
if (length(dispersion) == 1L) dispersion <- rep(dispersion, nrow(ql_table))
if (length(residual_df) == 1L) residual_df <- rep(residual_df, nrow(ql_table))
design_matrix <- fit$design
contrast_vector <- drop(contrast)
ci_lower <- ci_upper <- rep(NA_real_, nrow(ql_table))
for (i in seq_len(nrow(ql_table))) {
  mu <- fit$fitted.values[i, ]
  weights <- mu / (1 + dispersion[i] * mu)
  covariance <- tryCatch(
    solve(crossprod(design_matrix, design_matrix * weights)) * ql_variance[i],
    error = function(error) NULL
  )
  if (is.null(covariance)) next
  variance <- drop(t(contrast_vector) %*% covariance %*% contrast_vector) / log(2)^2
  if (!is.finite(variance) || variance < 0 || !is.finite(residual_df[i]) ||
      residual_df[i] <= 0) next
  margin <- qt(0.975, df = residual_df[i]) * sqrt(variance)
  ci_lower[i] <- ql_table$logFC[i] - margin
  ci_upper[i] <- ql_table$logFC[i] + margin
}

panel_result <- do.call(rbind, lapply(seq_along(panel), function(i) {
  gene <- panel[i]
  gene_id <- gene_ids[[gene]]
  row_index <- match(gene_id, rownames(ql_table))
  log_cpm <- cpm(y, log = TRUE, prior.count = 2)[gene_id, ]
  is_detected <- detected[i]
  data.frame(
    Gene = gene,
    Discovery_logFC = panel_source$Discovery_logFC[i],
    Discovery_P_Value = panel_source$Discovery_P_Value[i],
    Gene_ID = gene_id,
    Detected = is_detected,
    LV_Patients_ge10 = detection_by_group$LV[i],
    LV_Patients_total = sum(metadata$validation_group == "LV"),
    SV_Patients_ge10 = detection_by_group$SV[i],
    SV_Patients_total = sum(metadata$validation_group == "SV"),
    Mean_logCPM = mean(log_cpm),
    Validation_logFC = if (is_detected) ql_table$logFC[row_index] else NA_real_,
    CI95_Lower = if (is_detected) ci_lower[row_index] else NA_real_,
    CI95_Upper = if (is_detected) ci_upper[row_index] else NA_real_,
    P_Value = if (is_detected) ql_table$PValue[row_index] else NA_real_,
    BH_Panel = NA_real_,
    Direction_Concordant = if (is_detected) {
      sign(ql_table$logFC[row_index]) == sign(panel_source$Discovery_logFC[i])
    } else {
      NA
    },
    stringsAsFactors = FALSE
  )
}))
detected_rows <- which(panel_result$Detected)
panel_result$BH_Panel[detected_rows] <- p.adjust(panel_result$P_Value[detected_rows], method = "BH")
write.csv(panel_result, file.path(out_dir, "per_gene_validation.csv"), row.names = FALSE, na = "")
write.csv(
  data.frame(
    Dataset = accession,
    Sample_ID = metadata$sample_id,
    Group = as.character(metadata$validation_group),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "sample_group_mapping.csv"),
  row.names = FALSE
)

concordant_count <- sum(panel_result$Direction_Concordant %in% TRUE)
all_detected <- all(panel_result$Detected)
composite_auc <- composite_p <- auc_ci_lower <- auc_ci_upper <- NA_real_
score_difference <- NA_real_
if (all_detected) {
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  expression <- logcpm[gene_ids, metadata$sample_id, drop = FALSE]
  gene_sd <- apply(expression, 1, sd)
  if (any(!is.finite(gene_sd) | gene_sd == 0)) {
    stop("Cannot score the fixed panel: at least one gene has zero variance.")
  }
  z <- t(scale(t(expression), center = TRUE, scale = TRUE))
  scores <- colMeans(z)
  lv <- metadata$validation_group == "LV"
  auc_rank <- function(values, positive) {
    ranks <- rank(values, ties.method = "average")
    n1 <- sum(positive)
    n0 <- length(values) - n1
    (sum(ranks[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  composite_auc <- auc_rank(scores, lv)
  score_difference <- mean(scores[lv]) - mean(scores[!lv])
  composite_p <- wilcox.test(scores[lv], scores[!lv], exact = FALSE)$p.value
  set.seed(42)
  bootstrap_auc <- replicate(10000L, {
    sampled <- c(
      sample(which(lv), sum(lv), replace = TRUE),
      sample(which(!lv), sum(!lv), replace = TRUE)
    )
    auc_rank(scores[sampled], c(rep(TRUE, sum(lv)), rep(FALSE, sum(!lv))))
  })
  ci <- quantile(bootstrap_auc, c(0.025, 0.975), names = FALSE)
  auc_ci_lower <- ci[1]
  auc_ci_upper <- ci[2]
  write.csv(
    data.frame(
      Sample_ID = metadata$sample_id,
      Group = as.character(metadata$validation_group),
      Composite_Z_Score = as.numeric(scores)
    ),
    file.path(out_dir, "sample_composite_scores.csv"),
    row.names = FALSE
  )
}

sign_p <- if (sum(panel_result$Detected) > 0L) {
  binom.test(
    concordant_count, sum(panel_result$Detected),
    p = 0.5, alternative = "greater"
  )$p.value
} else {
  NA_real_
}
contrast_label <- if (accession == "GSE248762") {
  "LV (LV_UF + LV_NOT_UF) vs SV"
} else {
  "Long-term PD (LV) vs short-term PD (SV)"
}
summary <- data.frame(
  Dataset = accession,
  Contrast = contrast_label,
  LV_Patients = sum(metadata$validation_group == "LV"),
  SV_Patients = sum(metadata$validation_group == "SV"),
  Panel_Genes = length(panel),
  Detected_Genes = sum(panel_result$Detected),
  Direction_Concordant_Detected_Genes = concordant_count,
  Direction_Concordance_Fraction = if (sum(panel_result$Detected) > 0L) {
    concordant_count / sum(panel_result$Detected)
  } else {
    NA_real_
  },
  One_Sided_Sign_Test_P = sign_p,
  Full_Panel_Composite_Estimable = all_detected,
  Composite_AUC = composite_auc,
  Composite_AUC_CI95_Lower = auc_ci_lower,
  Composite_AUC_CI95_Upper = auc_ci_upper,
  Composite_Wilcoxon_P = composite_p,
  Composite_Mean_Score_Difference_LV_minus_SV = score_difference,
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_dir, "panel_validation_summary.csv"), row.names = FALSE, na = "")

message(
  accession, ": detected ", sum(panel_result$Detected), "/", length(panel),
  " genes; concordant ", concordant_count, "/",
  sum(panel_result$Detected),
  if (all_detected) paste0("; full-panel AUC=", round(composite_auc, 3)) else
    "; fixed full-panel composite not estimable"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
