# ==============================================================================
# SCRIPT 02: DIFFERENTIAL EXPRESSION ANALYSIS
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Dataset: GSE62928 (Pre-computed limma top-table)
# Comparison: Peritoneal Fibrosis vs Control
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(ggrepel)
})

cat("=================================================================\n")
cat("STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS (GSE62928 Top-Table)\n")
cat("=================================================================\n")

if (!file.exists("data/GSE62928_preprocessed.RData")) {
  stop("Preprocessed data not found! Please run 01_load_qc_preprocess.R first.")
}
load("data/GSE62928_preprocessed.RData")

# all_genes is the gene-level best-probe collapsed data frame
cat("Loaded", nrow(all_genes), "unique genes from GSE62928\n")
cat("Available columns:", paste(colnames(all_genes), collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# 1. DEG Classification using Significance Thresholds
# Thresholds: |log2FC| >= 0.585 (1.5-fold) and P.Value < 0.05 (nominal)
# ------------------------------------------------------------------------------
logfc_cutoff  <- 0.585
pval_cutoff   <- 0.05
adjp_cutoff   <- 0.05

all_results <- all_genes %>%
  mutate(Gene = Gene.symbol)

all_results$DEG_Status <- "Not Significant"
all_results$DEG_Status[all_results$logFC >  logfc_cutoff & all_results$P.Value < pval_cutoff] <- "Upregulated"
all_results$DEG_Status[all_results$logFC < -logfc_cutoff & all_results$P.Value < pval_cutoff] <- "Downregulated"
all_results$DEG_Status <- factor(all_results$DEG_Status,
                                 levels = c("Upregulated", "Downregulated", "Not Significant"))

deg_table <- all_results %>% filter(DEG_Status != "Not Significant")

cat("\n=================================================================\n")
cat("DEG FILTERING RESULTS (|log2FC| >=", logfc_cutoff, ", P.Value <", pval_cutoff, "):\n")
cat("Total genes tested:    ", nrow(all_results), "\n")
cat("Total DEGs identified: ", nrow(deg_table), "\n")
cat("  Upregulated:         ", sum(all_results$DEG_Status == "Upregulated"), "\n")
cat("  Downregulated:       ", sum(all_results$DEG_Status == "Downregulated"), "\n")
cat("Genes with FDR < 0.05:", sum(all_results$adj.P.Val < adjp_cutoff, na.rm = TRUE), "\n")
cat("=================================================================\n")

write.csv(all_results, "results/tables/GSE62928_all_results.csv", row.names = FALSE)
write.csv(deg_table,   "results/tables/GSE62928_DEGs_filtered.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. Volcano Plot
# ------------------------------------------------------------------------------
top_up   <- deg_table %>% filter(DEG_Status == "Upregulated")   %>% arrange(P.Value) %>% head(15)
top_down <- deg_table %>% filter(DEG_Status == "Downregulated") %>% arrange(P.Value) %>% head(15)
labeled_genes <- rbind(top_up, top_down)

p_volcano <- ggplot(all_results, aes(x = logFC, y = -log10(P.Value), color = DEG_Status)) +
  geom_point(alpha = 0.55, size = 1.6) +
  scale_color_manual(values = c("Upregulated"     = "#d73027",
                                "Downregulated"   = "#4575b4",
                                "Not Significant" = "#bdbdbd")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_text_repel(data = labeled_genes, aes(label = Gene),
                  size = 3.5, fontface = "bold", max.overlaps = 25,
                  box.padding = 0.4, point.padding = 0.3, show.legend = FALSE) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE62928 Volcano Plot: Peritoneal Fibrosis vs Control",
       subtitle = paste0("Cutoffs: |log2FC| >= ", logfc_cutoff, ", P < ", pval_cutoff,
                         "  (Up: ", sum(all_results$DEG_Status == "Upregulated"),
                         ", Down: ", sum(all_results$DEG_Status == "Downregulated"), ")"),
       x = "log2 Fold Change",
       y = "-log10(P-Value)")

ggsave("results/figures/DEG_01_volcano_plot.pdf", p_volcano, width = 9, height = 7)
ggsave("results/figures/DEG_01_volcano_plot.png", p_volcano, width = 9, height = 7, dpi = 300)
cat("Saved volcano plot to results/figures/DEG_01_volcano_plot.png\n")

# ------------------------------------------------------------------------------
# 3. B-statistic vs logFC Plot (B = log-odds of DE; proxy for expression intensity)
# ------------------------------------------------------------------------------
p_ma <- ggplot(all_results, aes(x = B, y = logFC, color = DEG_Status)) +
  geom_point(alpha = 0.45, size = 1.4) +
  scale_color_manual(values = c("Upregulated"     = "#d73027",
                                "Downregulated"   = "#4575b4",
                                "Not Significant" = "#d9d9d9")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-logfc_cutoff, logfc_cutoff),
             linetype = "dashed", color = "#636363", linewidth = 0.5) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE62928: log-Odds (B) vs log2FC",
       subtitle = "B = log-odds of differential expression (limma)",
       x = "B statistic (log-odds of DE)", y = "log2 Fold Change")

ggsave("results/figures/DEG_02_MA_plot.pdf", p_ma, width = 9, height = 6)
ggsave("results/figures/DEG_02_MA_plot.png", p_ma, width = 9, height = 6, dpi = 300)
cat("Saved MA plot to results/figures/DEG_02_MA_plot.png\n")

# Save results
save(all_results, deg_table, file = "data/GSE62928_limma_results.RData")
cat("Script 02 finished successfully.\n")
