# ==============================================================================
# SCRIPT 02: DIFFERENTIAL EXPRESSION ANALYSIS (LIMMA)
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Dataset: GSE125498 (Platform: GPL10558 Illumina HumanHT-12 V4.0)
# Comparison: LPD (Long-term, n=13) vs SPD (Short-term, n=20)
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(ggrepel)
})

cat("=================================================================\n")
cat("STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS (LIMMA)\n")
cat("=================================================================\n")

if (!file.exists("data/GSE125498_preprocessed.RData")) {
  stop("Preprocessed data not found! Please run 01_load_qc_preprocess.R first.")
}
load("data/GSE125498_preprocessed.RData")

group <- pheno$Group
cat("Sample groups for linear modeling:\n")
print(table(Group = group))

# ------------------------------------------------------------------------------
# 1. Limma Linear Modeling & Empirical Bayes Moderation
# ------------------------------------------------------------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- c("SPD", "LPD")
rownames(design) <- colnames(expr_clean)

fit <- lmFit(expr_clean, design)
contrast_mat <- makeContrasts(LPD_vs_SPD = LPD - SPD, levels = design)
fit_contrasts <- contrasts.fit(fit, contrast_mat)
fit_bayes <- eBayes(fit_contrasts)

all_results <- topTable(fit_bayes, coef = "LPD_vs_SPD", number = Inf, adjust.method = "BH", sort.by = "P")
all_results$Gene <- rownames(all_results)

cat("\nDifferential expression summary across all", nrow(all_results), "genes:\n")
cat("Genes with unadjusted P < 0.05:", sum(all_results$P.Value < 0.05), "\n")
cat("Genes with FDR (adj.P.Val) < 0.05:", sum(all_results$adj.P.Val < 0.05), "\n")

# ------------------------------------------------------------------------------
# 2. DEG Filtering using Significance Thresholds
# Thresholds: |log2FC| > 1.0 and adj.P.Val < 0.05
# ------------------------------------------------------------------------------
logfc_cutoff <- 1.0
adjp_cutoff <- 0.05

all_results$DEG_Status <- "Not Significant"
all_results$DEG_Status[all_results$logFC > logfc_cutoff & all_results$adj.P.Val < adjp_cutoff] <- "Upregulated in LPD"
all_results$DEG_Status[all_results$logFC < -logfc_cutoff & all_results$adj.P.Val < adjp_cutoff] <- "Downregulated in LPD"
all_results$DEG_Status <- factor(all_results$DEG_Status, 
                                 levels = c("Upregulated in LPD", "Downregulated in LPD", "Not Significant"))

deg_table <- all_results %>% filter(DEG_Status != "Not Significant")

cat("\n=================================================================\n")
cat("DEG FILTERING RESULTS (|log2FC| >", logfc_cutoff, ", adj.P.Val <", adjp_cutoff, "):\n")
cat("Total DEGs identified:", nrow(deg_table), "\n")
cat("Upregulated in LPD:", sum(all_results$DEG_Status == "Upregulated in LPD"), "\n")
cat("Downregulated in LPD:", sum(all_results$DEG_Status == "Downregulated in LPD"), "\n")
cat("=================================================================\n")

write.csv(all_results, "results/tables/GSE125498_all_limma_results.csv", row.names = FALSE)
write.csv(deg_table, "results/tables/GSE125498_DEGs_filtered.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Volcano Plot
# ------------------------------------------------------------------------------
top_up <- deg_table %>% filter(DEG_Status == "Upregulated in LPD") %>% arrange(adj.P.Val) %>% head(15)
top_down <- deg_table %>% filter(DEG_Status == "Downregulated in LPD") %>% arrange(adj.P.Val) %>% head(15)
labeled_genes <- rbind(top_up, top_down)

p_volcano <- ggplot(all_results, aes(x = logFC, y = -log10(adj.P.Val), color = DEG_Status)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Upregulated in LPD" = "#d73027", 
                                "Downregulated in LPD" = "#4575b4", 
                                "Not Significant" = "#bdbdbd")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = -log10(adjp_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_text_repel(data = labeled_genes, aes(label = Gene), 
                  size = 3.5, fontface = "bold", max.overlaps = 25, 
                  box.padding = 0.4, point.padding = 0.3, show.legend = FALSE) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE125498 Volcano Plot: Long-term PD (LPD) vs Short-term PD (SPD)",
       subtitle = paste0("Cutoffs: |log2FC| > ", logfc_cutoff, ", FDR < ", adjp_cutoff, 
                         " (Up: ", sum(all_results$DEG_Status == "Upregulated in LPD"), 
                         ", Down: ", sum(all_results$DEG_Status == "Downregulated in LPD"), ")"),
       x = "log2 Fold Change (LPD vs SPD)",
       y = "-log10(Adjusted P-Value)")

ggsave("results/figures/DEG_01_volcano_plot.pdf", p_volcano, width = 9, height = 7)
ggsave("results/figures/DEG_01_volcano_plot.png", p_volcano, width = 9, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 4. Expression Heatmap of Top DEGs
# ------------------------------------------------------------------------------
if (nrow(deg_table) > 0) {
  top_deg_genes <- deg_table %>% arrange(adj.P.Val) %>% head(min(50, nrow(deg_table))) %>% pull(Gene)
  heat_mat <- expr_clean[top_deg_genes, ]
  heat_mat_scaled <- t(scale(t(heat_mat)))
  
  annotation_col <- data.frame(Group = pheno$Group)
  rownames(annotation_col) <- colnames(heat_mat_scaled)
  
  ann_colors <- list(Group = c("SPD" = "#2b8cbe", "LPD" = "#de2d26"))
  
  png("results/figures/DEG_02_top_degs_heatmap.png", width = 3000, height = 3600, res = 300)
  pheatmap(heat_mat_scaled,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("#4575b4", "#f7f7f7", "#d73027"))(100),
           show_colnames = TRUE,
           show_rownames = TRUE,
           fontsize_row = 8,
           fontsize_col = 8,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           main = "Top Differentially Expressed Genes (GSE125498: LPD vs SPD)")
  dev.off()
}

save(all_results, deg_table, file = "data/GSE125498_limma_results.RData")
cat("Script 02 finished successfully.\n")
