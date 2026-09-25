# ==============================================================================
# SCRIPT 04: FUNCTIONAL ENRICHMENT ANALYSIS (CLUSTERPROFILER)
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Dataset: GSE62928 ECM-DEGs (GSE62928 ∩ ECM genes all.xlsx)
# Functional Profiling: GO (BP, CC, MF) and KEGG Pathways
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

cat("=================================================================\n")
cat("STEP 6: FUNCTIONAL ENRICHMENT OF GSE62928 ECM-DEGs\n")
cat("=================================================================\n")

if (!file.exists("data/GSE62928_ECM_DEGs.RData")) {
  stop("ECM-DEG results not found! Please run 03_matrisome_filtering.R first.")
}
load("data/GSE62928_ECM_DEGs.RData")
load("data/GSE62928_limma_results.RData")

target_genes  <- unique(ecm_degs_table$Gene)
all_deg_genes <- unique(deg_table$Gene)

cat("ECM-DEGs for enrichment:", length(target_genes), "\n")
cat("All DEGs for background enrichment:", length(all_deg_genes), "\n")

# Entrez ID mapping
entrez_mapping <- tryCatch({
  bitr(target_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) return(NULL))

all_deg_entrez <- tryCatch({
  bitr(all_deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) return(NULL))

# ------------------------------------------------------------------------------
# GO Enrichment Helper
# ------------------------------------------------------------------------------
run_go_enrichment <- function(genes, ont_type, title_prefix = "ECM-DEGs") {
  cat("\nRunning GO (", ont_type, ") enrichment for", title_prefix, "...\n")
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "SYMBOL",
                  ont           = ont_type,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.20,
                  readable      = TRUE)
  return(ego)
}

# ------------------------------------------------------------------------------
# GO Enrichment: BP, CC, MF on ECM-DEGs
# ------------------------------------------------------------------------------
go_bp <- run_go_enrichment(target_genes, "BP", "GSE62928 ECM-DEGs")
go_cc <- run_go_enrichment(target_genes, "CC", "GSE62928 ECM-DEGs")
go_mf <- run_go_enrichment(target_genes, "MF", "GSE62928 ECM-DEGs")

if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  df_bp <- as.data.frame(go_bp)
  write.csv(df_bp, "results/tables/GO_BP_enrichment_ECM_DEGs.csv", row.names = FALSE)
  p_bp <- dotplot(go_bp, showCategory = min(15, nrow(df_bp)),
                  title = "GO Biological Process (GSE62928 ECM-DEGs)") +
    theme_bw(base_size = 12) +
    scale_color_gradient(low = "#d73027", high = "#4575b4")
  ggsave("results/figures/Enrichment_01_GO_BP_dotplot.pdf", p_bp, width = 8, height = 6)
  ggsave("results/figures/Enrichment_01_GO_BP_dotplot.png", p_bp, width = 8, height = 6, dpi = 300)
  cat("GO BP: Found", nrow(df_bp), "enriched terms.\n")
}

if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
  df_cc <- as.data.frame(go_cc)
  write.csv(df_cc, "results/tables/GO_CC_enrichment_ECM_DEGs.csv", row.names = FALSE)
  p_cc <- dotplot(go_cc, showCategory = min(15, nrow(df_cc)),
                  title = "GO Cellular Component (GSE62928 ECM-DEGs)") +
    theme_bw(base_size = 12) +
    scale_color_gradient(low = "#d73027", high = "#4575b4")
  ggsave("results/figures/Enrichment_02_GO_CC_dotplot.pdf", p_cc, width = 8, height = 6)
  ggsave("results/figures/Enrichment_02_GO_CC_dotplot.png", p_cc, width = 8, height = 6, dpi = 300)
  cat("GO CC: Found", nrow(df_cc), "enriched terms.\n")
}

if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
  df_mf <- as.data.frame(go_mf)
  write.csv(df_mf, "results/tables/GO_MF_enrichment_ECM_DEGs.csv", row.names = FALSE)
  p_mf <- dotplot(go_mf, showCategory = min(15, nrow(df_mf)),
                  title = "GO Molecular Function (GSE62928 ECM-DEGs)") +
    theme_bw(base_size = 12) +
    scale_color_gradient(low = "#d73027", high = "#4575b4")
  ggsave("results/figures/Enrichment_03_GO_MF_dotplot.pdf", p_mf, width = 8, height = 6)
  ggsave("results/figures/Enrichment_03_GO_MF_dotplot.png", p_mf, width = 8, height = 6, dpi = 300)
  cat("GO MF: Found", nrow(df_mf), "enriched terms.\n")
}

# ------------------------------------------------------------------------------
# KEGG Pathway Enrichment
# ------------------------------------------------------------------------------
if (!is.null(entrez_mapping) && nrow(entrez_mapping) > 0) {
  kegg_res <- tryCatch({
    enrichKEGG(gene          = entrez_mapping$ENTREZID,
               organism      = "hsa",
               pvalueCutoff  = 0.10,
               pAdjustMethod = "BH")
  }, error = function(e) return(NULL))

  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    df_kegg  <- as.data.frame(kegg_res)
    write.csv(df_kegg, "results/tables/KEGG_enrichment_ECM_DEGs.csv", row.names = FALSE)
    p_kegg <- dotplot(kegg_res, showCategory = min(15, nrow(df_kegg)),
                      title = "KEGG Pathway Enrichment (GSE62928 ECM-DEGs)") +
      theme_bw(base_size = 12) +
      scale_color_gradient(low = "#d73027", high = "#4575b4")
    ggsave("results/figures/Enrichment_04_KEGG_dotplot.pdf", p_kegg, width = 8, height = 6)
    ggsave("results/figures/Enrichment_04_KEGG_dotplot.png", p_kegg, width = 8, height = 6, dpi = 300)
    cat("KEGG: Found", nrow(df_kegg), "enriched pathways.\n")
  }
}

# ------------------------------------------------------------------------------
# Global DEG GO BP enrichment (all DEGs background)
# ------------------------------------------------------------------------------
go_all_bp <- tryCatch({
  enrichGO(gene          = all_deg_genes,
           OrgDb         = org.Hs.eg.db,
           keyType       = "SYMBOL",
           ont           = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           readable      = TRUE)
}, error = function(e) return(NULL))

if (!is.null(go_all_bp) && nrow(as.data.frame(go_all_bp)) > 0) {
  df_all_bp <- as.data.frame(go_all_bp)
  write.csv(df_all_bp, "results/tables/GO_BP_enrichment_All_DEGs.csv", row.names = FALSE)
  p_all_bp <- dotplot(go_all_bp, showCategory = 15,
                      title = "GO Biological Process (All DEGs: GSE62928)") +
    theme_bw(base_size = 12) +
    scale_color_gradient(low = "#d73027", high = "#4575b4")
  ggsave("results/figures/Enrichment_05_Global_DEGs_GO_BP.png", p_all_bp, width = 9, height = 7, dpi = 300)
}

save(go_bp, go_cc, go_mf, file = "data/GSE62928_Enrichment_results.RData")
cat("Script 04 finished successfully.\n")
