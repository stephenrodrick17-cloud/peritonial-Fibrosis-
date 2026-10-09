args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: gse248762_ecm1_leave_one_out.R <analysis-output-directory>")
}
out_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages(library(edgeR))

counts <- as.matrix(read.csv(
  file.path(out_dir, "all_cell_pseudobulk_counts.csv"),
  row.names = 1L,
  check.names = FALSE
))
storage.mode(counts) <- "integer"
metadata <- read.csv(
  file.path(out_dir, "all_cell_pseudobulk_metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!identical(colnames(counts), metadata$sample_id)) {
  stop("Pseudobulk counts and metadata columns are not aligned")
}
ecm1_id <- "ENSG00000143369"
if (!ecm1_id %in% rownames(counts)) {
  stop("ECM1 feature ID missing from pseudobulk counts")
}

estimate_ecm1 <- function(selected_samples, omitted_sample = NA_character_) {
  selected <- metadata$sample_id %in% selected_samples
  meta <- metadata[selected, , drop = FALSE]
  matrix <- counts[, selected, drop = FALSE]
  group <- factor(meta$primary_group, levels = c("SV", "LV"))
  if (any(table(group) < 2L)) {
    return(data.frame(
      omitted_sample = omitted_sample,
      omitted_group = ifelse(is.na(omitted_sample), NA_character_,
        metadata$original_group[match(omitted_sample, metadata$sample_id)]),
      n_SV = sum(group == "SV"),
      n_LV = sum(group == "LV"),
      log2FC = NA_real_,
      PValue = NA_real_,
      status = "fewer than two patients in a group"
    ))
  }
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  y <- DGEList(counts = matrix)
  keep <- filterByExpr(y, design = design)
  detected_by_group <- vapply(levels(group), function(level) {
    samples_in_group <- group == level
    rowSums(matrix[, samples_in_group, drop = FALSE] >= 10L) >=
      ceiling(sum(samples_in_group) / 2)
  }, logical(nrow(matrix)))
  keep <- keep | apply(detected_by_group, 1L, all)
  keep[rownames(y) == ecm1_id] <- TRUE
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, contrast = c(SV = -1, LV = 1))
  table <- topTags(test, n = Inf, sort.by = "none")$table
  if (!ecm1_id %in% rownames(table)) {
    stop("ECM1 was unexpectedly filtered from leave-one-out analysis")
  }
  data.frame(
    omitted_sample = omitted_sample,
    omitted_group = ifelse(is.na(omitted_sample), NA_character_,
      metadata$original_group[match(omitted_sample, metadata$sample_id)]),
    n_SV = sum(group == "SV"),
    n_LV = sum(group == "LV"),
    log2FC = table[ecm1_id, "logFC"],
    PValue = table[ecm1_id, "PValue"],
    status = "estimated"
  )
}

full <- estimate_ecm1(metadata$sample_id)
full$change_from_full <- NA_real_
full$absolute_change_from_full <- NA_real_
leave_one_out <- do.call(
  rbind,
  lapply(metadata$sample_id, function(sample_id) {
    estimate_ecm1(
      setdiff(metadata$sample_id, sample_id),
      omitted_sample = sample_id
    )
  })
)
leave_one_out$change_from_full <- leave_one_out$log2FC - full$log2FC
leave_one_out$absolute_change_from_full <- abs(leave_one_out$change_from_full)
write.csv(
  rbind(
    transform(full, analysis = "full"),
    transform(leave_one_out, analysis = "leave_one_out")
  ),
  file.path(out_dir, "ecm1_leave_one_patient_out.csv"),
  row.names = FALSE,
  na = ""
)
most_influential <- leave_one_out[
  which.max(leave_one_out$absolute_change_from_full),
  ,
  drop = FALSE
]
write.csv(
  data.frame(
    full_log2FC = full$log2FC,
    most_influential_patient = most_influential$omitted_sample,
    patient_group = most_influential$omitted_group,
    leave_one_out_log2FC = most_influential$log2FC,
    absolute_change = most_influential$absolute_change_from_full,
    all_leave_one_out_log2FC_positive = all(leave_one_out$log2FC > 0),
    minimum_leave_one_out_log2FC = min(leave_one_out$log2FC, na.rm = TRUE),
    maximum_leave_one_out_log2FC = max(leave_one_out$log2FC, na.rm = TRUE)
  ),
  file.path(out_dir, "ecm1_influence_summary.csv"),
  row.names = FALSE,
  na = ""
)
print(full)
print(most_influential)
