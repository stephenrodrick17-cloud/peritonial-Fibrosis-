#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(jsonlite)
})

SEED <- 20261013L
N_PERM <- 10000L
MIN_SIZE <- 15L
MAX_SIZE <- 500L
GSEA_EXPONENT <- 1
LARGE_NES_CUTOFF <- 1.5

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve this script's path.")
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(out_dir, "..", "..", ".."), winslash = "/")
ranking_path <- file.path(out_dir, "hallmark_ranked_genes.csv")
library_path <- file.path(root, "results", "revision", "panel_v2_context", "hallmark_gene_sets.json")
if (!file.exists(ranking_path)) stop("Run 01_rank_limma.R before this step.")
if (!file.exists(library_path)) stop("Cached MSigDB Hallmark gene-set snapshot is missing.")

ranked <- fread(ranking_path)
if (anyDuplicated(ranked$Gene) || any(!is.finite(ranked$t))) stop("Ranking genes must be unique with finite t statistics.")
stats <- setNames(ranked$t, ranked$Gene)
stats <- sort(stats, decreasing = TRUE)

cached <- jsonlite::fromJSON(library_path, simplifyVector = FALSE)
library_name <- cached$source_library
if (!identical(library_name, "MSigDB_Hallmark_2020")) {
  stop("Unexpected cached Hallmark source: ", library_name)
}
raw_sets <- cached$gene_sets
pathways <- lapply(raw_sets, function(genes) intersect(unique(toupper(unlist(genes, use.names = FALSE))), names(stats)))
pathway_sizes <- lengths(pathways)
eligible <- pathway_sizes >= MIN_SIZE & pathway_sizes <= MAX_SIZE
excluded <- data.table(
  pathway = names(pathways)[!eligible],
  genes_in_ranked_universe = pathway_sizes[!eligible],
  reason = ifelse(pathway_sizes[!eligible] < MIN_SIZE, "below minimum set size", "above maximum set size")
)
fwrite(excluded, file.path(out_dir, "hallmark_sets_excluded_by_size.csv"))
pathways <- pathways[eligible]
if (!length(pathways)) stop("No Hallmark sets met the configured size filters.")

set.seed(SEED)
res <- fgsea::fgseaSimple(
  pathways = pathways,
  stats = stats,
  nperm = N_PERM,
  minSize = MIN_SIZE,
  maxSize = MAX_SIZE,
  gseaParam = GSEA_EXPONENT,
  nproc = 1
)
res <- as.data.table(res)
res[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
setnames(res, c("pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge"),
         c("pathway", "size", "ES", "NES", "P_value", "FDR", "leading_edge_genes"))
setorder(res, FDR, P_value)
res[, passes_FDR_0_05 := !is.na(FDR) & FDR < 0.05]
fwrite(res, file.path(out_dir, "hallmark_gsea_all.csv"))
fwrite(res[NES > 0 & passes_FDR_0_05], file.path(out_dir, "hallmark_gsea_upregulated.csv"))
fwrite(res[NES < 0 & passes_FDR_0_05], file.path(out_dir, "hallmark_gsea_downregulated.csv"))
high_nonsig <- res[abs(NES) >= LARGE_NES_CUTOFF & (is.na(FDR) | FDR >= 0.05)]
high_nonsig[, significance_label := "Non-significant"]
fwrite(high_nonsig, file.path(out_dir, "hallmark_gsea_large_NES_non_significant.csv"))

meta <- data.table(
  item = c(
    "run_date", "cohort", "MSigDB_collection_identifier", "gene_set_source",
    "total_cached_Hallmark_sets", "sets_tested", "sets_excluded_by_size",
    "ranked_gene_count", "limma_source", "contrast", "covariates",
    "GSEA_software", "R_version", "limma_version", "fgsea_version",
    "data_table_version", "jsonlite_version", "gene_set_permutations",
    "min_gene_set_size", "max_gene_set_size", "scoring_exponent",
    "FDR_cutoff", "seed", "NES_direction", "non_significant_large_NES_rule"
  ),
  value = c(
    as.character(Sys.Date()), "GSE125498", library_name,
    "Project-cached Enrichr MSigDB_Hallmark_2020 human gene sets; retrieved 2026-10-11",
    length(raw_sets), length(pathways), nrow(excluded), length(stats),
    "Existing Validation/GSE125498.top.table.tsv", "LPD versus SPD; LPD numerator",
    "none", "fgseaSimple", as.character(getRversion()),
    as.character(packageVersion("limma")), as.character(packageVersion("fgsea")),
    as.character(packageVersion("data.table")), as.character(packageVersion("jsonlite")), N_PERM,
    MIN_SIZE, MAX_SIZE, GSEA_EXPONENT, 0.05, SEED,
    "positive NES = LPD; negative NES = SPD",
    paste0("abs(NES) >= ", LARGE_NES_CUTOFF, " and FDR >= 0.05")
  )
)
fwrite(meta, file.path(out_dir, "hallmark_gsea_metadata.csv"))
message("Tested ", length(pathways), " Hallmark sets; ", sum(res$passes_FDR_0_05), " passed FDR < 0.05.")
