# ==============================================================================
# SCRIPT 01: DATA LOADING, QC, PREPROCESSING & PROBE COLLAPSING
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Dataset: GSE125498 (Platform: GPL10558 Illumina HumanHT-12 V4.0)
# Comparison: Long-term PD (LPD, n=13, >=25 mo) vs Short-term PD (SPD, n=20, 0-24 mo)
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
cat("STEP 1: LOADING DATASET GSE125498\n")
cat("=================================================================\n")

matrix_file <- file.path("data", "GSE125498_series_matrix.txt.gz")
if (!file.exists(matrix_file)) {
  cat("Downloading GSE125498 series matrix...\n")
  gse_list <- getGEO("GSE125498", GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data")
  gse <- gse_list[[1]]
} else {
  cat("Loading local Series Matrix file:", matrix_file, "\n")
  gse_res <- GEOquery:::parseGSEMatrix(matrix_file, destdir = "data", AnnotGPL = FALSE, getGPL = FALSE)
  gse <- gse_res$eset
}

# ------------------------------------------------------------------------------
# Metadata extraction & Group labeling
# ------------------------------------------------------------------------------
cat("\nParsing clinical metadata and phenotype annotations...\n")
pheno <- pData(gse)

# Classify samples into SPD (0-24 months) and LPD (>=25 months) based on sample titles
group_vec <- rep(NA, nrow(pheno))
group_vec[grepl("long", pheno$title, ignore.case = TRUE)] <- "LPD"
group_vec[grepl("short", pheno$title, ignore.case = TRUE)] <- "SPD"

if (any(is.na(group_vec))) {
  char_cols <- grep("characteristics|source_name", colnames(pheno), value = TRUE)
  raw_chars <- apply(pheno[, char_cols, drop=FALSE], 1, paste, collapse=" ")
  group_vec[is.na(group_vec) & grepl("long|LPD|>= ?25", raw_chars, ignore.case = TRUE)] <- "LPD"
  group_vec[is.na(group_vec) & grepl("short|SPD|0-24|< ?25", raw_chars, ignore.case = TRUE)] <- "SPD"
}

pheno$Group <- factor(group_vec, levels = c("SPD", "LPD"))
pData(gse)$Group <- pheno$Group

cat("\nSample distribution:\n")
print(table(Group = pheno$Group))
cat("Total samples:", nrow(pheno), "\n")

# ------------------------------------------------------------------------------
# STEP 2: QUALITY CONTROL (Boxplots, Density, PCA, Outlier Analysis)
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 2: QUALITY CONTROL & DISTRIBUTION CHECKS\n")
cat("=================================================================\n")

raw_expr <- exprs(gse)
cat("Expression matrix dimensions (Probes x Samples):", nrow(raw_expr), "x", ncol(raw_expr), "\n")

# Check if log2 transformation is required
qx <- as.numeric(quantile(raw_expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  cat("Applying log2 transformation to expression matrix (max value was:", round(qx[6], 2), ")...\n")
  raw_expr[raw_expr <= 0] <- NaN
  expr_mat <- log2(raw_expr)
} else {
  cat("Data is already on log2 scale (max value:", round(qx[6], 2), ").\n")
  expr_mat <- raw_expr
}

# Handle any missing values via row median imputation
if (any(is.na(expr_mat))) {
  n_na <- sum(is.na(expr_mat))
  cat("Found", n_na, "missing values (", round(n_na / length(expr_mat) * 100, 3), "%). Imputing with row medians...\n")
  row_meds <- apply(expr_mat, 1, median, na.rm = TRUE)
  for (i in which(rowSums(is.na(expr_mat)) > 0)) {
    expr_mat[i, is.na(expr_mat[i, ])] <- row_meds[i]
  }
}

# 1. Boxplot of Expression Distributions
set.seed(42)
sampled_rows <- sample(seq_len(nrow(expr_mat)), min(5000, nrow(expr_mat)))
df_long <- data.frame(
  Sample = rep(colnames(expr_mat), each = length(sampled_rows)),
  Group = rep(pheno$Group, each = length(sampled_rows)),
  Expression = as.vector(expr_mat[sampled_rows, ])
)

p_box <- ggplot(df_long, aes(x = Sample, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("SPD" = "#2b8cbe", "LPD" = "#de2d26")) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        legend.position = "top") +
  labs(title = "GSE125498 Sample Expression Distributions",
       x = "Sample ID (GEO GSM)", y = "Log2 Signal Intensity")

ggsave("results/figures/QC_01_expression_boxplots.pdf", p_box, width = 11, height = 6)
ggsave("results/figures/QC_01_expression_boxplots.png", p_box, width = 11, height = 6, dpi = 300)
cat("Saved QC boxplot to results/figures/QC_01_expression_boxplots.png\n")

# 2. PCA Plot
pca_res <- prcomp(t(expr_mat), scale. = TRUE)
pca_df <- data.frame(
  Sample = colnames(expr_mat),
  Group = pheno$Group,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(aes(shape = Group), size = 4, alpha = 0.9) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.90) +
  scale_color_manual(values = c("SPD" = "#2b8cbe", "LPD" = "#de2d26")) +
  scale_fill_manual(values = c("SPD" = "#2b8cbe", "LPD" = "#de2d26")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE125498 Unsupervised PCA: LPD vs SPD",
       x = paste0("PC1 (", var_explained[1], "% variance)"),
       y = paste0("PC2 (", var_explained[2], "% variance)"))

ggsave("results/figures/QC_02_PCA_plot.pdf", p_pca, width = 8, height = 6)
ggsave("results/figures/QC_02_PCA_plot.png", p_pca, width = 8, height = 6, dpi = 300)
cat("Saved PCA plot to results/figures/QC_02_PCA_plot.png\n")

# 3. Hierarchical Clustering Dendrogram
sample_dists <- dist(t(expr_mat))
sample_clust <- hclust(sample_dists, method = "average")
png("results/figures/QC_03_sample_hierarchical_clustering.png", width = 2400, height = 1600, res = 300)
plot(sample_clust, main = "GSE125498 Sample Hierarchical Clustering (Average Linkage)",
     xlab = "Samples", sub = "", col = ifelse(pheno$Group == "LPD", "#de2d26", "#2b8cbe"),
     lwd = 1.5, cex = 0.8)
legend("topright", legend = c("LPD", "SPD"), fill = c("#de2d26", "#2b8cbe"))
dev.off()
cat("Saved Hierarchical Clustering to results/figures/QC_03_sample_hierarchical_clustering.png\n")

# ------------------------------------------------------------------------------
# STEP 3: PROBE-TO-GENE ANNOTATION & COLLAPSING (MaxMean Strategy)
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 3: PROBE-TO-GENE ANNOTATION & COLLAPSING (MaxMean)\n")
cat("=================================================================\n")

annot_file <- "data/GPL10558.annot.gz"
probe_gene_map <- NULL

if (file.exists(annot_file)) {
  cat("Parsing platform annotation table from", annot_file, "...\n")
  gpl_con <- gzfile(annot_file, "rt")
  lines <- readLines(gpl_con, n = 100)
  close(gpl_con)
  skip_lines <- grep("!platform_table_begin", lines)
  if (length(skip_lines) > 0) {
    gpl_table <- read.table(gzfile(annot_file), skip = skip_lines, header = TRUE, 
                            sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE,
                            fill = TRUE)
    probe_col <- grep("^ID$", colnames(gpl_table), value = TRUE)
    sym_col <- grep("Gene.symbol|Gene symbol", colnames(gpl_table), value = TRUE)
    if (length(sym_col) > 0 && length(probe_col) > 0) {
      probe_gene_map <- data.frame(
        ID = gpl_table[[probe_col[1]]],
        Gene = gpl_table[[sym_col[1]]],
        stringsAsFactors = FALSE
      )
    }
  }
}

if (is.null(probe_gene_map)) {
  feature_data <- fData(gse)
  sym_col <- grep("symbol|gene_symbol|Gene Symbol|ILMN_Gene", colnames(feature_data), value = TRUE, ignore.case = TRUE)
  if (length(sym_col) > 0) {
    probe_gene_map <- data.frame(
      ID = rownames(feature_data),
      Gene = as.character(feature_data[[sym_col[1]]]),
      stringsAsFactors = FALSE
    )
  }
}

probe_gene_map <- probe_gene_map %>% filter(ID %in% rownames(expr_mat))
probe_gene_map$Gene <- sub(" ///.*", "", probe_gene_map$Gene)
probe_gene_map$Gene <- trimws(probe_gene_map$Gene)
probe_gene_map$Gene[probe_gene_map$Gene == "" | probe_gene_map$Gene == "---" | probe_gene_map$Gene == "NA"] <- NA

cat("Total Probes in Expression Matrix:", nrow(expr_mat), "\n")
cat("Probes mapped to valid Gene Symbols:", sum(!is.na(probe_gene_map$Gene)), "\n")
cat("Unique Gene Symbols mapped:", length(unique(na.omit(probe_gene_map$Gene))), "\n")

valid_map <- probe_gene_map %>% filter(!is.na(Gene))
expr_annot <- expr_mat[valid_map$ID, ]
genes_annot <- valid_map$Gene

cat("\nCollapsing multiple probes per gene using Maximum Mean Expression (MaxMean)...\n")
row_means <- rowMeans(expr_annot)
split_idx <- split(seq_len(nrow(expr_annot)), genes_annot)
selected_probe_indices <- sapply(split_idx, function(idx) {
  idx[which.max(row_means[idx])]
})

expr_collapsed <- expr_annot[selected_probe_indices, ]
rownames(expr_collapsed) <- names(split_idx)

cat("Dimensions of gene-level expression matrix:", nrow(expr_collapsed), "genes x", ncol(expr_collapsed), "samples\n")

vars <- apply(expr_collapsed, 1, var)
keep_var <- vars > quantile(vars, 0.10)
expr_clean <- expr_collapsed[keep_var, ]
cat("Retained", nrow(expr_clean), "genes after low-variance filtering.\n")

save(gse, pheno, expr_mat, expr_clean, file = "data/GSE125498_preprocessed.RData")
write.csv(data.frame(Gene = rownames(expr_clean), expr_clean), 
          "results/tables/GSE125498_normalized_gene_expression.csv", row.names = FALSE)

cat("Successfully saved preprocessed data to data/GSE125498_preprocessed.RData\n")
cat("Script 01 finished successfully.\n")
