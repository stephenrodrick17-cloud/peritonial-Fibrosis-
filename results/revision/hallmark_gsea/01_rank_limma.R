#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve this script's path.")
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(out_dir, "..", "..", ".."), winslash = "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

probe_table_path <- file.path(root, "Validation", "GSE125498.top.table.tsv")
metadata_path <- file.path(root, "results", "tables", "GSE125498_sample_metadata.csv")
if (!file.exists(probe_table_path)) stop("Saved probe-level limma table is missing: ", probe_table_path)
if (!file.exists(metadata_path)) stop("Discovery metadata is missing: ", metadata_path)

meta <- fread(metadata_path)
if (!all(c("sample_id", "group") %in% names(meta))) stop("Sample metadata lacks sample_id/group.")
if (nrow(meta) != 33L || sum(meta$group == "LPD") != 13L || sum(meta$group == "SPD") != 20L) {
  stop("Expected GSE125498 metadata with 13 LPD and 20 SPD samples.")
}

probe <- fread(probe_table_path)
required <- c("ID", "t", "logFC", "P.Value", "adj.P.Val", "Gene.symbol")
if (!all(required %in% names(probe))) stop("Saved limma table columns do not match expected schema.")
probe[, Gene := trimws(as.character(Gene.symbol))]
probe[Gene %in% c("", "NA", "NaN", "nan", "None"), Gene := NA_character_]
probe[, t := as.numeric(t)]
probe[, logFC := as.numeric(logFC)]
probe[, P.Value := as.numeric(P.Value)]
probe[, adj.P.Val := as.numeric(adj.P.Val)]
probe <- probe[!is.na(Gene) & is.finite(t)]
if (!nrow(probe)) stop("No annotated finite limma statistics were found.")

probe[, abs_t := abs(t)]
setorder(probe, Gene, -abs_t, ID)
collapsed <- probe[, .SD[1L], by = Gene]
setorder(collapsed, -t, Gene)
collapsed[, rank := seq_len(.N)]
collapsed <- collapsed[, .(rank, Gene, probe_id = ID, t, logFC, P.Value, adj.P.Val)]
fwrite(collapsed, file.path(out_dir, "hallmark_ranked_genes.csv"))

metadata <- data.table(
  item = c(
    "source_limma_table", "cohort", "samples", "LPD_samples", "SPD_samples",
    "contrast", "covariates", "probe_rows_with_gene_and_finite_t",
    "unique_genes_ranked", "duplicate_gene_collapse",
    "ranking_statistic", "direction"
  ),
  value = c(
    "Validation/GSE125498.top.table.tsv", "GSE125498", nrow(meta), sum(meta$group == "LPD"),
    sum(meta$group == "SPD"), "LPD versus SPD; LPD numerator", "none",
    nrow(probe), nrow(collapsed), "largest absolute moderated t-statistic per gene",
    "limma empirical Bayes moderated t-statistic, descending",
    "positive t-statistic indicates higher expression in LPD"
  )
)
fwrite(metadata, file.path(out_dir, "ranking_metadata.csv"))
message("Reused saved limma table; ranked ", nrow(collapsed), " unique genes.")
