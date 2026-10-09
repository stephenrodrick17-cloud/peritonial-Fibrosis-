args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: compartment_edger_gsea_input.R <input-dir> <output-rankings.csv> <output-models.csv>")
}

input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
rankings_path <- args[[2]]
models_path <- args[[3]]
suppressPackageStartupMessages(library(edgeR))

counts_table <- read.csv(
  file.path(input_dir, "celltype_pseudobulk_counts.csv"),
  row.names = 1L,
  check.names = FALSE
)
metadata <- read.csv(
  file.path(input_dir, "celltype_pseudobulk_metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
gene_map <- read.csv(
  file.path(input_dir, "gene_id_to_symbol.csv"),
  stringsAsFactors = FALSE
)
counts <- as.matrix(counts_table)
storage.mode(counts) <- "integer"
if (any(counts < 0L)) {
  stop("Cell-type pseudobulk input contains negative counts")
}
symbol_by_id <- setNames(gene_map$gene_symbol, gene_map$gene_id)
groups <- c("LV_NOT_UF", "LV_UF")
compartments <- c(
  "monocyte/macrophage",
  "T/NK",
  "fibroblast",
  "mesothelial"
)
ranking_rows <- list()
model_rows <- list()

for (cell_type in compartments) {
  selected <- metadata$cell_type == cell_type &
    metadata$original_group %in% groups
  meta <- metadata[selected, , drop = FALSE]
  sample_columns <- paste(meta$sample_id, cell_type, sep = "|")
  matched_columns <- match(sample_columns, colnames(counts))
  if (anyNA(matched_columns)) {
    stop(sprintf("%s: patient-to-pseudobulk mapping failed", cell_type))
  }
  meta <- meta[order(matched_columns), , drop = FALSE]
  sample_columns <- sample_columns[order(matched_columns)]
  counts_ct <- counts[, sample_columns, drop = FALSE]
  group <- factor(meta$original_group, levels = groups)
  if (any(table(group) < 2L)) {
    stop(sprintf(
      "%s: fewer than two eligible patient pseudobulks in a contrast group",
      cell_type
    ))
  }

  design <- model.matrix(~ group)
  y <- DGEList(counts = counts_ct)
  keep <- filterByExpr(y, design = design)
  if (!any(keep)) {
    stop(sprintf("%s: edgeR filterByExpr retained no genes", cell_type))
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "groupLV_UF")
  table <- topTags(test, n = Inf, sort.by = "none")$table
  ids <- rownames(table)
  symbols <- unname(symbol_by_id[ids])
  symbols[is.na(symbols) | symbols == ""] <- ids[
    is.na(symbols) | symbols == ""
  ]
  ranking_rows[[cell_type]] <- data.frame(
    compartment = cell_type,
    gene_id = ids,
    gene_symbol = symbols,
    logFC = table$logFC,
    F = table$F,
    PValue = table$PValue,
    FDR = p.adjust(table$PValue, method = "BH"),
    rank_score = sign(table$logFC) * sqrt(table$F),
    stringsAsFactors = FALSE
  )
  model_rows[[cell_type]] <- data.frame(
    compartment = cell_type,
    n_LV_UF = sum(group == "LV_UF"),
    n_LV_NOT_UF = sum(group == "LV_NOT_UF"),
    edgeR_version = as.character(packageVersion("edgeR")),
    R_version = R.version.string,
    design = "~ group; LV_UF relative to LV_NOT_UF",
    filtering = "edgeR filterByExpr",
    normalization = "TMM",
    test = "robust quasi-likelihood F-test",
    stringsAsFactors = FALSE
  )
}

write.csv(do.call(rbind, ranking_rows), rankings_path, row.names = FALSE)
write.csv(do.call(rbind, model_rows), models_path, row.names = FALSE)
