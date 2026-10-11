#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
script_path <- normalizePath(script_arg, winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

input_dir <- file.path(repo_root, "results", "revision", "effluent_validation", "GSE92455")
panel_path <- file.path(repo_root, "results", "revision", "panel_v2", "panel_v2.json")
output_dir <- file.path(repo_root, "results", "revision", "panel_v2", "GSE92455_20261011")
matrix_path <- file.path(input_dir, "GSE92455-GPL6480_series_matrix.txt.gz")
annotation_path <- file.path(input_dir, "GPL6480.annot.gz")
metadata_path <- file.path(input_dir, "sample_metadata.csv")

required <- c(matrix_path, annotation_path, metadata_path, panel_path)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Required input file(s) missing: ", paste(missing, collapse = ", "))
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = FALSE))) {
  stop("Refusing to overwrite non-empty output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_geo_row <- function(lines, key) {
  matches <- lines[startsWith(lines, paste0(key, "\t"))]
  if (!length(matches)) {
    stop("GEO metadata row not found: ", key)
  }
  parsed <- read.delim(
    text = matches[[1]],
    header = FALSE,
    sep = "\t",
    quote = "\"",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  as.character(unlist(parsed[1, -1, drop = FALSE], use.names = FALSE))
}

matrix_lines <- readLines(gzfile(matrix_path), warn = FALSE)
table_start <- which(matrix_lines == "!series_matrix_table_begin")
table_end <- which(matrix_lines == "!series_matrix_table_end")
if (length(table_start) != 1L || length(table_end) != 1L || table_end <= table_start) {
  stop("Could not identify one valid GEO series matrix table.")
}
matrix_text <- paste(matrix_lines[(table_start + 1L):(table_end - 1L)], collapse = "\n")
series_matrix <- read.delim(
  text = matrix_text,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!"ID_REF" %in% names(series_matrix)) {
  stop("GSE92455 series matrix has no ID_REF column.")
}

sample_ids <- names(series_matrix)[-1]
source_ch1 <- read_geo_row(matrix_lines, "!Sample_source_name_ch1")
source_ch2 <- read_geo_row(matrix_lines, "!Sample_source_name_ch2")
label_ch1 <- read_geo_row(matrix_lines, "!Sample_label_ch1")
label_ch2 <- read_geo_row(matrix_lines, "!Sample_label_ch2")
if (any(lengths(list(source_ch1, source_ch2, label_ch1, label_ch2)) != length(sample_ids))) {
  stop("GEO channel metadata length does not match the series matrix columns.")
}

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata[metadata$platform == "GPL6480", , drop = FALSE]
if (!setequal(metadata$sample_id, sample_ids)) {
  stop("GPL6480 metadata sample IDs do not match the series matrix columns.")
}
metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
metadata$source_ch1 <- source_ch1
metadata$source_ch2 <- source_ch2
metadata$label_ch1 <- label_ch1
metadata$label_ch2 <- label_ch2
metadata$matrix_column <- seq_along(sample_ids)

clinical_groups <- c("Epithelioid effluent HPMC", "Non-epithelioid effluent HPMC")
metadata <- metadata[metadata$group_from_metadata %in% clinical_groups, , drop = FALSE]
metadata$group <- ifelse(
  metadata$group_from_metadata == "Epithelioid effluent HPMC",
  "Epithelioid",
  "NonEpithelioid"
)
metadata$group <- factor(metadata$group, levels = c("NonEpithelioid", "Epithelioid"))

is_patient_ch1 <- grepl("Effluent-derived HPMC", metadata$source_ch1, fixed = TRUE)
is_reference_ch2 <- grepl("Pooled HPMCs derived from 4 omentum", metadata$source_ch2, fixed = TRUE)
selected_samples <- metadata[is_patient_ch1 & is_reference_ch2, , drop = FALSE]
if (nrow(selected_samples) != 17L ||
    anyDuplicated(selected_samples$patient_id_reported) ||
    anyNA(selected_samples$patient_id_reported)) {
  stop("Expected exactly one consistently oriented GPL6480 array per each of 17 patient IDs.")
}
group_counts <- table(selected_samples$group)
if (!identical(as.integer(group_counts[c("Epithelioid", "NonEpithelioid")]), c(9L, 8L))) {
  stop("GSE92455 patient-level group counts do not match the metadata-defined 9/8 contrast.")
}
if (!all(selected_samples$label_ch1 == "Cy3") || !all(selected_samples$label_ch2 == "Cy5")) {
  stop("Selected arrays do not share the expected channel-label orientation.")
}

all_values <- as.matrix(series_matrix[, selected_samples$sample_id, drop = FALSE])
storage.mode(all_values) <- "numeric"
rownames(all_values) <- make.unique(as.character(series_matrix$ID_REF))
if (any(is.infinite(all_values))) {
  stop("Selected GEO processed matrix contains infinite values; inspect before analysis.")
}

annotation <- read.delim(
  gzfile(annotation_path),
  skip = 28,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(c("ID", "Gene symbol") %in% names(annotation))) {
  stop("GPL6480 annotation lacks expected ID and Gene symbol columns.")
}
annotation <- annotation[, c("ID", "Gene symbol"), drop = FALSE]
annotation$ID <- as.character(annotation$ID)
annotation$Gene.symbol.tokens <- strsplit(
  ifelse(is.na(annotation$`Gene symbol`), "", annotation$`Gene symbol`),
  "\\s*///\\s*"
)
row_annotation <- annotation[match(rownames(all_values), annotation$ID), , drop = FALSE]

panel <- fromJSON(panel_path, simplifyVector = TRUE)
genes <- as.character(panel$genes_passing_selection_rules)
if (!length(genes) || anyDuplicated(genes)) {
  stop("Committed panel_v2.json has an empty or duplicated gene list.")
}

probe_rows <- lapply(genes, function(gene) {
  rows <- which(vapply(
    row_annotation$Gene.symbol.tokens,
    function(tokens) gene %in% tokens,
    logical(1)
  ))
  if (!length(rows)) {
    return(data.frame(
      gene = gene,
      probe_id = NA_character_,
      mean_GEO_value = NA_real_,
      n_annotated_probes = 0L,
      status = "no annotated probe",
      stringsAsFactors = FALSE
    ))
  }
  means <- rowMeans(all_values[rows, , drop = FALSE], na.rm = TRUE)
  valid <- is.finite(means)
  if (!any(valid)) {
    return(data.frame(
      gene = gene,
      probe_id = NA_character_,
      mean_GEO_value = NA_real_,
      n_annotated_probes = length(rows),
      status = "annotated probes have no finite values",
      stringsAsFactors = FALSE
    ))
  }
  winner_index <- which(valid)[which.max(means[valid])]
  winner <- rows[winner_index]
  data.frame(
    gene = gene,
    probe_id = rownames(all_values)[winner],
    mean_GEO_value = means[winner_index],
    n_annotated_probes = length(rows),
    status = "mapped",
    stringsAsFactors = FALSE
  )
})
probe_selection <- do.call(rbind, probe_rows)
write.csv(probe_selection, file.path(output_dir, "probe_selection.csv"), row.names = FALSE, na = "")
write.csv(selected_samples, file.path(output_dir, "patient_level_sample_metadata.csv"), row.names = FALSE, na = "")

mapped <- probe_selection$status == "mapped"
per_gene <- NULL
if (any(mapped)) {
  design <- model.matrix(~ group, data = selected_samples)
  fit <- eBayes(lmFit(all_values, design))
  coefficient <- which(colnames(design) == "groupEpithelioid")
  if (length(coefficient) != 1L) {
    stop("Could not identify the Epithelioid-versus-NonEpithelioid limma coefficient.")
  }
  top <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none", adjust.method = "none")
  top$probe_id <- rownames(top)
  selected_rows <- match(probe_selection$probe_id[mapped], top$probe_id)
  if (anyNA(selected_rows)) {
    stop("Selected probe IDs were not present in the limma results.")
  }
  fit_rows <- match(probe_selection$probe_id[mapped], rownames(fit$coefficients))
  moderated_se <- fit$stdev.unscaled[fit_rows, coefficient] * sqrt(fit$s2.post[fit_rows])
  total_df <- fit$df.total[fit_rows]
  result <- probe_selection[mapped, c("gene", "probe_id", "mean_GEO_value", "n_annotated_probes"), drop = FALSE]
  result$logFC_GEO_value_scale <- top$logFC[selected_rows]
  result$CI95_lower <- result$logFC_GEO_value_scale - qt(0.975, total_df) * moderated_se
  result$CI95_upper <- result$logFC_GEO_value_scale + qt(0.975, total_df) * moderated_se
  result$P_value <- top$P.Value[selected_rows]
  result$BH_q_within_mapped_panel_genes <- p.adjust(result$P_value, method = "BH")
  result$direction_higher_in_epithelioid <- ifelse(
    result$logFC_GEO_value_scale > 0,
    "Epithelioid",
    ifelse(result$logFC_GEO_value_scale < 0, "NonEpithelioid", "tie")
  )
  result$status <- "limma tested"
  missing_genes <- probe_selection[!mapped, , drop = FALSE]
  if (nrow(missing_genes)) {
    missing_genes$logFC_GEO_value_scale <- NA_real_
    missing_genes$CI95_lower <- NA_real_
    missing_genes$CI95_upper <- NA_real_
    missing_genes$P_value <- NA_real_
    missing_genes$BH_q_within_mapped_panel_genes <- NA_real_
    missing_genes$direction_higher_in_epithelioid <- NA_character_
    missing_genes$status <- "not tested: no annotated probe"
  }
  per_gene <- rbind(result, missing_genes)
  per_gene <- per_gene[match(genes, per_gene$gene), , drop = FALSE]
} else {
  per_gene <- probe_selection
  per_gene$logFC_GEO_value_scale <- NA_real_
  per_gene$CI95_lower <- NA_real_
  per_gene$CI95_upper <- NA_real_
  per_gene$P_value <- NA_real_
  per_gene$BH_q_within_mapped_panel_genes <- NA_real_
  per_gene$direction_higher_in_epithelioid <- NA_character_
  per_gene$status <- "not tested: no annotated probe"
}
write.csv(per_gene, file.path(output_dir, "panel_v2_per_gene_limma.csv"), row.names = FALSE, na = "")

score_summary <- list(
  accession = "GSE92455",
  analysis_label = "post-freeze panel v2 exploratory phenotype test; not preregistered",
  panel_freeze_commit = system("git rev-parse HEAD", intern = TRUE),
  contrast = "Epithelioid versus Non-epithelioid effluent-derived HPMCs",
  panel_genes = genes,
  n_epithelioid_patients = as.integer(group_counts[["Epithelioid"]]),
  n_non_epithelioid_patients = as.integer(group_counts[["NonEpithelioid"]]),
  score_estimable = all(mapped),
  auc_higher_score_epithelioid = NA_real_,
  auc_bootstrap_ci95 = c(NA_real_, NA_real_),
  mann_whitney_p_two_sided = NA_real_,
  mann_whitney_method = NA_character_,
  bootstrap_resamples = 10000L,
  seed = 42L,
  matched_null_p = "not computed; GSE92455 is a distinct phenotype contrast and the frozen vintage null rule does not apply",
  interpretation = paste(
    "The contrast is not a long-vintage versus short-vintage test. The panel-v2",
    "composite was added only after the dated panel-freeze commit; it does not",
    "alter the original GSE92455 preregistered per-gene-only analysis."
  )
)

sample_scores <- data.frame(
  sample_id = selected_samples$sample_id,
  patient_id = selected_samples$patient_id_reported,
  group = as.character(selected_samples$group),
  panel_v2_composite = NA_real_,
  analysis_label = score_summary$analysis_label,
  stringsAsFactors = FALSE
)
if (all(mapped)) {
  selected_probe_ids <- probe_selection$probe_id[match(genes, probe_selection$gene)]
  gene_values <- all_values[selected_probe_ids, , drop = FALSE]
  rownames(gene_values) <- genes
  gene_sd <- apply(gene_values, 1, sd)
  if (all(is.finite(gene_sd)) && all(gene_sd > 0)) {
    gene_z <- t(scale(t(gene_values)))
    composite <- colMeans(gene_z)
    sample_scores$panel_v2_composite <- composite[match(sample_scores$sample_id, names(composite))]
    epi <- sample_scores$panel_v2_composite[sample_scores$group == "Epithelioid"]
    non_epi <- sample_scores$panel_v2_composite[sample_scores$group == "NonEpithelioid"]
    pairwise <- outer(epi, non_epi, "-")
    auc <- mean(pairwise > 0) + 0.5 * mean(pairwise == 0)
    has_ties <- anyDuplicated(c(epi, non_epi)) > 0L
    wilcox_exact <- !has_ties
    wilcox_p <- wilcox.test(
      epi,
      non_epi,
      exact = wilcox_exact,
      correct = TRUE,
      alternative = "two.sided"
    )$p.value
    set.seed(42)
    bootstrap_auc <- replicate(10000, {
      epi_draw <- sample(epi, length(epi), replace = TRUE)
      non_epi_draw <- sample(non_epi, length(non_epi), replace = TRUE)
      differences <- outer(epi_draw, non_epi_draw, "-")
      mean(differences > 0) + 0.5 * mean(differences == 0)
    })
    score_summary$auc_higher_score_epithelioid <- auc
    score_summary$auc_bootstrap_ci95 <- as.numeric(quantile(
      bootstrap_auc,
      probs = c(0.025, 0.975),
      names = FALSE,
      type = 7
    ))
    score_summary$mann_whitney_p_two_sided <- wilcox_p
    score_summary$mann_whitney_method <- if (wilcox_exact) {
      "exact two-sided Wilcoxon rank-sum (no ties)"
    } else {
      "asymptotic two-sided Wilcoxon rank-sum with continuity correction (ties)"
    }
    score_summary$direction_of_auc_was_predeclared <- FALSE
  } else {
    score_summary$score_estimable <- FALSE
    score_summary$interpretation <- paste(
      score_summary$interpretation,
      "The fixed panel score was not estimable because at least one mapped gene had zero variance."
    )
  }
} else {
  score_summary$missing_panel_genes <- probe_selection$gene[!mapped]
  score_summary$interpretation <- paste(
    score_summary$interpretation,
    "The fixed panel score was not computed because one or more panel genes had no annotated probe."
  )
}
write.csv(sample_scores, file.path(output_dir, "panel_v2_sample_scores.csv"), row.names = FALSE, na = "")
write_json(score_summary, file.path(output_dir, "panel_v2_score_summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

hashes <- tools::md5sum(required)
manifest <- list(
  analysis_date = as.character(Sys.Date()),
  panel_freeze_commit = system("git rev-parse HEAD", intern = TRUE),
  panel_freeze_commit_date = system("git show -s --format=%ci HEAD", intern = TRUE),
  script = normalizePath(script_path, winslash = "/"),
  inputs = as.list(setNames(as.character(hashes), basename(names(hashes)))),
  R_version = R.version.string,
  limma_version = as.character(packageVersion("limma")),
  jsonlite_version = as.character(packageVersion("jsonlite")),
  session_info = capture.output(sessionInfo())
)
write_json(manifest, file.path(output_dir, "analysis_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

format_num <- function(x, digits = 4) {
  if (length(x) == 0L || is.na(x)) "not estimable" else formatC(x, digits = digits, format = "f")
}
report <- c(
  "# GSE92455: first panel-v2 test",
  "",
  paste0("**Analysis date:** ", as.character(Sys.Date())),
  paste0("**Panel freeze commit:** `", score_summary$panel_freeze_commit, "`."),
  "**Label:** post-freeze, exploratory panel-v2 test; not preregistered.",
  "",
  "## Cohort and processing",
  "",
  "The tested contrast is the supplementary phenotype comparison of epithelioid versus non-epithelioid effluent-derived HPMCs. It is not a vintage contrast and cannot establish long-vintage/short-vintage concordance. GPL6480 ex-vivo data were used; GPL6848 in-vitro treatment arrays and pooled-omentum reference samples were excluded.",
  "",
  paste0(
    "GEO metadata identified ", nrow(selected_samples), " distinct patient labels (",
    score_summary$n_epithelioid_patients, " epithelioid and ",
    score_summary$n_non_epithelioid_patients, " non-epithelioid). One array per patient was used, selecting the consistently oriented layout with the effluent sample in Cy3/channel 1 and pooled reference in Cy5/channel 2; reciprocal technical arrays were excluded so technical replicates were not treated as independent patients."
  ),
  "",
  "The series-matrix processing metadata states background subtraction and LOESS/LOWESS plus quantile normalization. The GEO record does not explicitly identify the logarithm base for the submitted matrix values. No additional normalization or channel sign reversal was applied; all analyzed arrays have the same channel layout. Limma estimates are reported on the submitted GEO value scale.",
  "",
  "## Fixed panel results",
  "",
  paste0(
    "The panel was read from the committed `panel_v2.json` before expression analysis. The fixed five-gene composite was ",
    ifelse(score_summary$score_estimable, "estimable.", "not estimable.")
  ),
  "",
  paste0(
    "AUC (higher score = epithelioid): ",
    format_num(score_summary$auc_higher_score_epithelioid),
    " (bootstrap 95% CI ",
    format_num(score_summary$auc_bootstrap_ci95[1]),
    "–",
    format_num(score_summary$auc_bootstrap_ci95[2]),
    "); ",
    score_summary$mann_whitney_method,
    ", P=",
    format_num(score_summary$mann_whitney_p_two_sided),
    ". Bootstrap resampling was stratified by patient group (10,000 resamples; seed 42). The phenotype direction was not predeclared, so the AUC was not flipped to maximize performance."
  ),
  "",
  "| Gene | Selected probe | Limma estimate (Epithelioid − Non-epithelioid) | 95% CI | P | BH q (panel) |",
  "|---|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(per_gene))) {
  row <- per_gene[i, ]
  report <- c(report, paste0(
    "| ", row$gene, " | ", ifelse(is.na(row$probe_id), "not mapped", row$probe_id), " | ",
    format_num(row$logFC_GEO_value_scale), " | ",
    format_num(row$CI95_lower), "–", format_num(row$CI95_upper), " | ",
    format_num(row$P_value), " | ", format_num(row$BH_q_within_mapped_panel_genes), " |"
  ))
}
report <- c(
  report,
  "",
  "No random-gene matched-null P was calculated. This distinct phenotype contrast is not covered by the frozen vintage-specific matched-null rule. The original preregistered GSE92455 analysis remains per-gene only; this additional composite is a new exploratory panel-v2 analysis and does not replace or alter it.",
  "",
  "## Reproducibility",
  "",
  "Saved outputs in this directory include the panel-v2 sample scores, patient-level metadata, outcome-blind probe selection, per-gene limma results, machine-readable score summary, and input/package manifest. The script refuses to overwrite existing non-empty outputs."
)
writeLines(report, file.path(output_dir, "GSE92455_panel_v2_report.md"), useBytes = TRUE)

cat("Panel v2 test completed for GSE92455.\n")
cat("Patients: ", nrow(selected_samples), " (", group_counts[["Epithelioid"]], " Epi, ",
    group_counts[["NonEpithelioid"]], " Non-Epi)\n", sep = "")
cat("Panel score estimable: ", score_summary$score_estimable, "\n", sep = "")
cat("AUC: ", format_num(score_summary$auc_higher_score_epithelioid),
    "; MW P: ", format_num(score_summary$mann_whitney_p_two_sided), "\n", sep = "")
cat("Outputs: ", output_dir, "\n", sep = "")
