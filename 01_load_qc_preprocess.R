# ==============================================================================
# SCRIPT 01: DATA LOADING, QC, PREPROCESSING & PROBE COLLAPSING
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Dataset: GSE62928 (Platform: GPL10558 Illumina HumanHT-12 V4.0)
# Comparison: Peritoneal fibrosis (PF) vs Control
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(ggrepel)
})

# Define output directories
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

cat("=================================================================\n")
cat("STEP 1: LOADING DATASET GSE62928\n")
cat("=================================================================\n")

# GSE62928 already exists as a pre-processed top-table TSV
# Load the limma top-table directly
gse_file <- file.path("GSE62928.top.table.tsv")
if (!file.exists(gse_file)) {
  stop("GSE62928.top.table.tsv not found in project directory!")
}

cat("Loading GSE62928 top-table from:", gse_file, "\n")
gse_toptable <- read.table(gse_file, sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE, quote = "")

cat("Loaded top-table dimensions (probes x columns):", nrow(gse_toptable), "x", ncol(gse_toptable), "\n")
cat("Columns:", paste(colnames(gse_toptable), collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# STEP 2: CLEAN & COLLAPSE PROBES TO UNIQUE GENE SYMBOLS
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 2: PROBE ANNOTATION & COLLAPSING (Best P-Value per Gene)\n")
cat("=================================================================\n")

# Remove probes without gene annotation
gse_annotated <- gse_toptable[!is.na(gse_toptable$Gene.symbol) &
                                 gse_toptable$Gene.symbol != "" &
                                 gse_toptable$Gene.symbol != "---", ]
cat("Probes with gene annotation:", nrow(gse_annotated), "/", nrow(gse_toptable), "\n")

# Expand multi-gene probes (e.g., "GENE1 /// GENE2")
gse_annotated$Gene.symbol <- trimws(gse_annotated$Gene.symbol)
has_multi <- grepl("///", gse_annotated$Gene.symbol)
cat("Multi-gene probes:", sum(has_multi), "\n")

# Expand by splitting on /// — vectorised (fast)
split_syms <- strsplit(gse_annotated$Gene.symbol, " ?/// ?")
n_syms     <- lengths(split_syms)
gse_expanded <- gse_annotated[rep(seq_len(nrow(gse_annotated)), n_syms), ]
gse_expanded$Gene.symbol <- trimws(unlist(split_syms))
cat("Rows after expansion:", nrow(gse_expanded), "\n")

# Keep best probe per gene (lowest P.Value)
gse_best <- gse_expanded[order(gse_expanded$P.Value), ]
gse_best <- gse_best[!duplicated(gse_best$Gene.symbol), ]
rownames(gse_best) <- gse_best$Gene.symbol

cat("Unique gene symbols (best probe per gene):", nrow(gse_best), "\n")

# Summary statistics of the collapsed dataset
cat("\nExpression statistics (logFC distribution):\n")
cat("  Min logFC:    ", round(min(gse_best$logFC, na.rm = TRUE), 3), "\n")
cat("  Max logFC:    ", round(max(gse_best$logFC, na.rm = TRUE), 3), "\n")
cat("  Mean |logFC|: ", round(mean(abs(gse_best$logFC), na.rm = TRUE), 3), "\n")
cat("  Genes P < 0.05:", sum(gse_best$P.Value < 0.05, na.rm = TRUE), "\n")
cat("  Genes FDR < 0.05:", sum(gse_best$adj.P.Val < 0.05, na.rm = TRUE), "\n")

# ------------------------------------------------------------------------------
# STEP 3: QUALITY CONTROL PLOTS
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 3: QUALITY CONTROL PLOTS\n")
cat("=================================================================\n")

# 1. LogFC Distribution
p_hist <- ggplot(gse_best, aes(x = logFC)) +
  geom_histogram(bins = 80, fill = "#4575b4", color = "white", alpha = 0.85) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed",
             color = "#d73027", linewidth = 0.8) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE62928: log2 Fold Change Distribution",
       subtitle = "Dashed lines at |log2FC| = 0.585 (1.5-fold change threshold)",
       x = "log2 Fold Change", y = "Number of Genes")

ggsave("results/figures/QC_01_logFC_distribution.pdf", p_hist, width = 9, height = 6)
ggsave("results/figures/QC_01_logFC_distribution.png", p_hist, width = 9, height = 6, dpi = 300)
cat("Saved logFC distribution to results/figures/QC_01_logFC_distribution.png\n")

# 2. P-value Distribution
p_pval <- ggplot(gse_best, aes(x = P.Value)) +
  geom_histogram(bins = 50, fill = "#74c476", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#d73027", linewidth = 0.8) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE62928: P-Value Distribution",
       subtitle = "Dashed line at P = 0.05",
       x = "Nominal P-Value", y = "Number of Genes")

ggsave("results/figures/QC_02_pvalue_distribution.pdf", p_pval, width = 8, height = 6)
ggsave("results/figures/QC_02_pvalue_distribution.png", p_pval, width = 8, height = 6, dpi = 300)
cat("Saved P-value distribution to results/figures/QC_02_pvalue_distribution.png\n")

# Save preprocessed gene-level data
all_genes <- gse_best
save(all_genes, gse_toptable, file = "data/GSE62928_preprocessed.RData")
write.csv(gse_best,
          "results/tables/GSE62928_gene_level_toptable.csv", row.names = FALSE)

cat("\nSuccessfully saved preprocessed data to data/GSE62928_preprocessed.RData\n")
cat("Script 01 finished successfully.\n")
