# ==============================================================================
# SCRIPT 03: ECM & MATRISOME FILTERING — INTERSECTION
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Dataset: GSE62928 DEGs ∩ ECM genes all.xlsx (Naba et al. Human Matrisome)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readxl)
})

cat("=================================================================\n")
cat("STEP 5: ECM GENE FILTERING — GSE62928 ∩ ECM MASTERLIST\n")
cat("=================================================================\n")

if (!file.exists("data/GSE62928_limma_results.RData")) {
  stop("Limma results not found! Please run 02_differential_expression.R first.")
}
load("data/GSE62928_limma_results.RData")

# ------------------------------------------------------------------------------
# 1. Load ECM Masterlist from "ECM genes all.xlsx"
# ------------------------------------------------------------------------------
ecm_file <- "ECM genes all.xlsx"
if (!file.exists(ecm_file)) {
  stop("ECM genes all.xlsx not found! Please ensure the file is in the project directory.")
}

cat("Loading ECM masterlist from:", ecm_file, "\n")
ecm_raw <- read_excel(ecm_file, skip = 1, col_names = TRUE)
colnames(ecm_raw)[1:4] <- c("Matrisome_Division", "Matrisome_Category",
                              "Gene_Symbol", "Gene_Name")

# Clean ECM gene symbols
ecm_clean <- ecm_raw %>%
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "Gene Symbol") %>%
  mutate(Gene_Symbol = trimws(as.character(Gene_Symbol))) %>%
  filter(Gene_Symbol != "" & Gene_Symbol != "NA")

cat("Total ECM masterlist genes:", nrow(ecm_clean), "\n")
cat("ECM Categories:\n")
print(table(ecm_clean$Matrisome_Division))

# Build a reference matrisome data frame
matrisome_df <- ecm_clean %>%
  select(Gene = Gene_Symbol,
         Category = Matrisome_Division,
         Subcategory = Matrisome_Category) %>%
  distinct(Gene, .keep_all = TRUE)

ecm_ref_genes <- unique(matrisome_df$Gene)
cat("\nUnique ECM reference genes:", length(ecm_ref_genes), "\n")

# ------------------------------------------------------------------------------
# 2. Find All GSE62928 Genes ∩ ECM Masterlist
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("COMPUTING INTERSECTION: GSE62928 ∩ ECM Masterlist\n")
cat("=================================================================\n")

gse_genes <- all_results$Gene

# All intersection (regardless of DEG status)
intersection_all <- all_results %>%
  filter(Gene %in% ecm_ref_genes) %>%
  left_join(matrisome_df, by = "Gene") %>%
  arrange(P.Value)

cat("GSE62928 total genes tested:    ", nrow(all_results), "\n")
cat("ECM masterlist genes:           ", length(ecm_ref_genes), "\n")
cat("Intersection (all):             ", nrow(intersection_all), "\n")

# DEG intersection only (|logFC| >= 0.585 & P < 0.05)
ecm_degs_table <- intersection_all %>%
  filter(DEG_Status != "Not Significant") %>%
  arrange(P.Value)

cat("ECM-DEGs (|logFC|>=0.585, P<0.05):", nrow(ecm_degs_table), "\n")
cat("  Upregulated ECM-DEGs:   ", sum(ecm_degs_table$DEG_Status == "Upregulated"), "\n")
cat("  Downregulated ECM-DEGs: ", sum(ecm_degs_table$DEG_Status == "Downregulated"), "\n")
cat("=================================================================\n")

# Borderline ECM genes (present in ECM but not reaching DEG threshold)
borderline_ecm <- intersection_all %>%
  filter(DEG_Status == "Not Significant" & P.Value < 0.05) %>%
  arrange(P.Value)
cat("Borderline ECM genes (P<0.05, non-DEG):", nrow(borderline_ecm), "\n")

# ------------------------------------------------------------------------------
# 3. Save Intersection Results
# ------------------------------------------------------------------------------
write.csv(intersection_all,   "results/tables/GSE62928_ECM_intersection.csv",       row.names = FALSE)
write.csv(ecm_degs_table,     "results/tables/GSE62928_ECM_DEGs_candidate_list.csv", row.names = FALSE)
write.csv(borderline_ecm,     "results/tables/GSE62928_ECM_Borderline.csv",          row.names = FALSE)

cat("\nSaved intersection tables to results/tables/\n")

# Venn summary
venn_counts <- data.frame(
  Set   = c("GSE62928 DEGs Only", "ECM-DEGs (Intersection)", "ECM Masterlist Only"),
  Count = c(nrow(deg_table) - nrow(ecm_degs_table),
            nrow(ecm_degs_table),
            length(ecm_ref_genes) - nrow(ecm_degs_table))
)
write.csv(venn_counts, "results/tables/GSE62928_Matrisome_Overlap_Counts.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. ECM Volcano Plot (Highlighting Intersection)
# ------------------------------------------------------------------------------
all_results$Is_ECM <- all_results$Gene %in% ecm_ref_genes
all_results$ECM_Highlight <- "Other Genes"
all_results$ECM_Highlight[all_results$Is_ECM & all_results$DEG_Status != "Not Significant"] <- "ECM-DEG"
all_results$ECM_Highlight[all_results$Is_ECM & all_results$DEG_Status == "Not Significant"] <- "Non-DEG Matrisome"
all_results$ECM_Highlight <- factor(all_results$ECM_Highlight,
                                    levels = c("ECM-DEG", "Non-DEG Matrisome", "Other Genes"))

ecm_labeled <- all_results %>% filter(ECM_Highlight == "ECM-DEG")

p_ecm_volcano <- ggplot(all_results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = all_results %>% filter(ECM_Highlight == "Other Genes"),
             color = "#d9d9d9", alpha = 0.35, size = 1.4) +
  geom_point(data = all_results %>% filter(ECM_Highlight == "Non-DEG Matrisome"),
             color = "#4575b4", alpha = 0.55, size = 1.8) +
  geom_point(data = all_results %>% filter(ECM_Highlight == "ECM-DEG"),
             color = "#d73027", alpha = 0.9, size = 3.2) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  geom_text_repel(data = ecm_labeled, aes(label = Gene),
                  size = 3.8, fontface = "bold", color = "#800026",
                  max.overlaps = 40, box.padding = 0.5, point.padding = 0.3) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE62928: ECM & Matrisome Gene Distribution",
       subtitle = paste0("ECM-DEGs (red): n=", nrow(ecm_labeled),
                         " | Non-DEG Matrisome (blue): n=",
                         sum(all_results$ECM_Highlight == "Non-DEG Matrisome")),
       x = "log2 Fold Change",
       y = "-log10(P-Value)")

ggsave("results/figures/ECM_01_volcano_highlight.pdf", p_ecm_volcano, width = 10, height = 7)
ggsave("results/figures/ECM_01_volcano_highlight.png", p_ecm_volcano, width = 10, height = 7, dpi = 300)
cat("Saved ECM volcano plot to results/figures/ECM_01_volcano_highlight.png\n")

# ------------------------------------------------------------------------------
# 5. Top ECM-DEGs Bar Chart
# ------------------------------------------------------------------------------
if (nrow(ecm_degs_table) > 0) {
  top_n_show <- min(30, nrow(ecm_degs_table))
  top_ecm <- ecm_degs_table %>%
    arrange(P.Value) %>%
    head(top_n_show) %>%
    mutate(Gene = factor(Gene, levels = rev(Gene)),
           Direction = ifelse(logFC > 0, "Upregulated", "Downregulated"))

  p_bar <- ggplot(top_ecm, aes(x = Gene, y = logFC, fill = Direction)) +
    geom_col(alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = c("Upregulated" = "#d73027", "Downregulated" = "#4575b4")) +
    coord_flip() +
    theme_classic(base_size = 12) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = paste0("Top ", top_n_show, " ECM-DEGs: GSE62928 ∩ ECM Masterlist"),
         subtitle = "Ordered by statistical significance (P-Value)",
         x = "Gene Symbol", y = "log2 Fold Change")

  ggsave("results/figures/ECM_02_top_ecm_degs_bar.pdf", p_bar, width = 9, height = 8)
  ggsave("results/figures/ECM_02_top_ecm_degs_bar.png", p_bar, width = 9, height = 8, dpi = 300)
  cat("Saved ECM top-DEGs bar chart to results/figures/ECM_02_top_ecm_degs_bar.png\n")
}

# Save RData
save(ecm_degs_table, intersection_all, borderline_ecm, matrisome_df,
     file = "data/GSE62928_ECM_DEGs.RData")
cat("\nScript 03 finished successfully.\n")
cat("Key output: results/tables/GSE62928_ECM_intersection.csv\n")
cat("  - Total ECM intersection genes:", nrow(intersection_all), "\n")
cat("  - ECM-DEGs (significant):      ", nrow(ecm_degs_table), "\n")
