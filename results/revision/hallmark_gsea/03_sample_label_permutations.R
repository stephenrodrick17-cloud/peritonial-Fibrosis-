#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(limma)
})

SEED <- 20261014L
LABEL_PERMUTATIONS <- 1000L
GSEA_PERMUTATIONS_PER_LABEL <- 500L
MIN_SIZE <- 15L
MAX_SIZE <- 500L
GSEA_EXPONENT <- 1

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve this script's path.")
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(out_dir, "..", "..", ".."), winslash = "/")
expr_path <- file.path(root, "results", "tables", "GSE125498_full_expression_matrix.csv")
meta_path <- file.path(root, "results", "tables", "GSE125498_sample_metadata.csv")
observed_path <- file.path(out_dir, "hallmark_gsea_all.csv")
gene_sets_path <- file.path(root, "results", "revision", "panel_v2_context", "hallmark_gene_sets.json")
for (path in c(expr_path, meta_path, observed_path, gene_sets_path)) {
  if (!file.exists(path)) stop("Required input missing: ", path)
}

expr_df <- fread(expr_path)
if (!"Gene" %in% names(expr_df) || anyDuplicated(expr_df$Gene)) stop("Expression matrix must have unique Gene symbols.")
expr <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene"), with = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- toupper(expr_df$Gene)
meta <- fread(meta_path)
if (!all(meta$sample_id %in% colnames(expr))) stop("Expression matrix lacks metadata samples.")
meta <- meta[match(colnames(expr), sample_id)]
if (!all(meta$group %in% c("SPD", "LPD")) || sum(meta$group == "LPD") != 13L || sum(meta$group == "SPD") != 20L) {
  stop("Expected fixed-size labels of 13 LPD and 20 SPD.")
}

cached <- jsonlite::fromJSON(gene_sets_path, simplifyVector = FALSE)
pathways <- lapply(cached$gene_sets, function(genes) intersect(unique(toupper(unlist(genes, use.names = FALSE))), rownames(expr)))
pathways <- pathways[lengths(pathways) >= MIN_SIZE & lengths(pathways) <= MAX_SIZE]
observed <- fread(observed_path)
observed <- observed[pathway %in% names(pathways)]
if (nrow(observed) != length(pathways)) stop("Observed GSEA pathways and permutation pathways do not match.")
permuted_path <- file.path(out_dir, "hallmark_sample_label_permutation_NES.csv")

fit_ranked_gsea <- function(labels) {
  group <- factor(labels, levels = c("SPD", "LPD"))
  design <- model.matrix(~ group)
  fit <- limma::eBayes(limma::lmFit(expr, design))
  tstats <- fit$t[, "groupLPD"]
  names(tstats) <- rownames(expr)
  tstats <- sort(tstats[is.finite(tstats)], decreasing = TRUE)
  gsea <- fgsea::fgseaSimple(
    pathways = pathways,
    stats = tstats,
    nperm = GSEA_PERMUTATIONS_PER_LABEL,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    gseaParam = GSEA_EXPONENT,
    nproc = 1
  )
  as.data.table(gsea)[, .(pathway, ES, NES)]
}

set.seed(SEED)
if (file.exists(permuted_path)) {
  permuted <- fread(permuted_path)
  complete <- nrow(permuted) == LABEL_PERMUTATIONS * length(pathways) &&
    uniqueN(permuted$permutation) == LABEL_PERMUTATIONS &&
    setequal(permuted$pathway, names(pathways))
  if (!complete) stop("Existing sample-label permutation table is incomplete or incompatible.")
  message("Reusing the completed ", LABEL_PERMUTATIONS, "-permutation NES table.")
} else {
  permutation_rows <- vector("list", LABEL_PERMUTATIONS)
  for (i in seq_len(LABEL_PERMUTATIONS)) {
    perm_labels <- sample(meta$group, length(meta$group), replace = FALSE)
    result <- fit_ranked_gsea(perm_labels)
    result[, permutation := i]
    permutation_rows[[i]] <- result
    if (i %% 25L == 0L) message("Completed ", i, " / ", LABEL_PERMUTATIONS, " random label permutations.")
  }
  permuted <- rbindlist(permutation_rows)
  setcolorder(permuted, c("permutation", "pathway", "ES", "NES"))
  fwrite(permuted, permuted_path)
}

setkey(permuted, pathway)
empirical <- observed[, .(pathway, observed_ES = ES, observed_NES = NES, observed_FDR = FDR)]
empirical[, empirical_P := vapply(seq_len(.N), function(j) {
  pathway_name <- pathway[[j]]
  obs_nes <- observed_NES[[j]]
  (sum(abs(permuted[pathway == pathway_name, NES]) >= abs(obs_nes)) + 1) /
    (LABEL_PERMUTATIONS + 1)
}, numeric(1))]
empirical[, abs_observed_NES := abs(observed_NES)]
setorder(empirical, empirical_P, -abs_observed_NES)
empirical[, abs_observed_NES := NULL]
fwrite(empirical, file.path(out_dir, "hallmark_sample_label_empirical_P.csv"))

observed_top <- observed[which.max(abs(NES))]
permuted_top_abs <- permuted[, max(abs(NES)), by = permutation]$V1
top_pathways <- permuted[, pathway[which.max(abs(NES))], by = permutation]$V1
top_rank <- data.table(
  observed_top_pathway = observed_top$pathway,
  observed_top_NES = observed_top$NES,
  rank_among_permuted_top_sets = 1L + sum(permuted_top_abs >= abs(observed_top$NES)),
  permutations_with_same_top_pathway = sum(top_pathways == observed_top$pathway),
  random_label_permutations = LABEL_PERMUTATIONS
)
fwrite(top_rank, file.path(out_dir, "hallmark_sample_label_top_set_rank.csv"))

settings <- data.table(
  item = c("run_date", "R_version", "limma_version", "fgsea_version",
           "random_label_permutations", "GSEA_simple_permutations_per_label", "seed",
           "group_sizes_preserved", "limma_model", "contrast", "covariates",
           "empirical_P_definition", "top_set_rank_definition"),
  value = c(
    as.character(Sys.Date()), as.character(getRversion()), as.character(packageVersion("limma")),
    as.character(packageVersion("fgsea")), LABEL_PERMUTATIONS, GSEA_PERMUTATIONS_PER_LABEL, SEED,
    "yes: 13 LPD and 20 SPD in every permutation",
    "limma lmFit followed by eBayes for each random label assignment",
    "LPD versus SPD; positive moderated t-statistics indicate LPD",
    "none",
    "(count of permuted |NES| >= observed |NES| + 1) / (random label permutations + 1)",
    "1 + number of permuted maximum-|NES| values at least as large as observed maximum |NES|"
  )
)
fwrite(settings, file.path(out_dir, "hallmark_sample_label_permutation_metadata.csv"))
message("Completed all random sample-label permutations.")
