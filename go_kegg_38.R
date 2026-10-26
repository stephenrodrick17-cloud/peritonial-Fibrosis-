#!/usr/bin/env Rscript
# GO/KEGG enrichment for 38 non-grey matrisome DEGs (cleaned: NAT6 & MBP removed)
# Background = 5,001 WGCNA input genes
# Test set = 38 matrisome DEGs in non-grey modules

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# --- Load data ---
mod <- read.csv("results/revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv",
                check.names = FALSE)
gcol <- if ("Gene.symbol" %in% names(mod)) "Gene.symbol" else "Annotated Gene"

# Remove NAT6 and MBP
drop <- c("NAT6", "MBP")
m <- mod[!mod[[gcol]] %in% drop, ]

# 38 non-grey matrisome DEGs
non_grey <- m[!is.na(m$Module) & m$Module != "grey", ]
test_genes <- unique(non_grey[[gcol]])
cat("Test genes:", length(test_genes), "\n")
cat("Test gene list:", paste(test_genes, collapse=", "), "\n\n")

# Background = 5,001 WGCNA input genes
wgcna_mods <- read.csv("results/revision/wgcna_modules_dynamic_taskB.csv",
                        check.names = FALSE)
bg_genes <- unique(wgcna_mods$Gene)
cat("Background genes:", length(bg_genes), "\n\n")

# Convert symbols to Entrez IDs
test_entrez <- bitr(test_genes, fromType="SYMBOL", toType="ENTREZID",
                    OrgDb=org.Hs.eg.db)
bg_entrez   <- bitr(bg_genes,   fromType="SYMBOL", toType="ENTREZID",
                    OrgDb=org.Hs.eg.db)

cat("Mapped test genes (Entrez):", nrow(test_entrez), "\n")
cat("Mapped background genes (Entrez):", nrow(bg_entrez), "\n\n")

# --- GO enrichment (BP, CC, MF) ---
for (ont in c("BP", "CC", "MF")) {
  cat("============================\n")
  cat("GO", ont, "enrichment\n")
  cat("============================\n")
  ego <- enrichGO(
    gene         = test_entrez$ENTREZID,
    universe     = bg_entrez$ENTREZID,
    OrgDb        = org.Hs.eg.db,
    ont          = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 10,
    maxGSSize     = 500,
    readable      = TRUE
  )
  res <- as.data.frame(ego)
  cat("Significant terms (BH p.adjust < 0.05):", nrow(res), "\n")
  if (nrow(res) > 0) {
    cat("\nTop 5 terms:\n")
    top5 <- head(res[order(res$p.adjust), c("ID","Description","GeneRatio",
                                             "BgRatio","pvalue","p.adjust","Count")], 5)
    print(top5, row.names = FALSE)
    # Save full results
    write.csv(res,
              paste0("results/revision/wgcna_paper_figures/GO_", ont,
                     "_enrichment_38_non_grey_cleaned.csv"),
              row.names = FALSE)
  }
  cat("\n")
}

# --- KEGG enrichment ---
cat("============================\n")
cat("KEGG pathway enrichment\n")
cat("============================\n")
ekegg <- enrichKEGG(
  gene         = test_entrez$ENTREZID,
  universe     = bg_entrez$ENTREZID,
  organism     = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 10,
  maxGSSize     = 500
)
kegg_res <- as.data.frame(ekegg)
cat("Significant KEGG pathways (BH p.adjust < 0.05):", nrow(kegg_res), "\n")
if (nrow(kegg_res) > 0) {
  cat("\nTop 5 pathways:\n")
  top5k <- head(kegg_res[order(kegg_res$p.adjust),
                         c("ID","Description","GeneRatio","BgRatio",
                           "pvalue","p.adjust","Count")], 5)
  print(top5k, row.names = FALSE)
  write.csv(kegg_res,
            "results/revision/wgcna_paper_figures/KEGG_enrichment_38_non_grey_cleaned.csv",
            row.names = FALSE)
}
cat("\n")

cat("=== SUMMARY ===\n")
cat("Background genes:", length(bg_genes), "\n")
cat("Test genes:", length(test_genes), "\n")
cat("Done.\n")
