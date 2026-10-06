
suppressPackageStartupMessages({
  library(edgeR)
})

hub_genes <- c("COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
               "EDIL3", "LOX", "INHBA", "ISM1", "COMP")

run_edger <- function(counts_file, meta_file, add_covariate=FALSE) {
  dt <- data.table::fread(counts_file, data.table = FALSE)
  rownames(dt) <- dt[[1]]
  counts <- dt[, -1, drop = FALSE]
  meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
  
  common <- intersect(rownames(meta), rownames(counts))
  meta <- meta[common, ]
  counts <- counts[common, ]
  
  mat <- t(as.matrix(counts))
  groups <- factor(meta$group, levels = c("LV_NOT_UF", "LV_UF", "SV"))
  
  if (add_covariate) {
    cov_val <- meta$meso_frac
    design <- model.matrix(~ 0 + groups + cov_val)
    colnames(design) <- c("LV_NOT_UF", "LV_UF", "SV", "meso_frac")
  } else {
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- levels(groups)
  }
  
  dge <- DGEList(counts = mat, group = groups)
  keep <- filterByExpr(dge, design)
  hub_in_mat <- intersect(hub_genes, rownames(mat))
  keep[hub_in_mat] <- TRUE
  dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
  dge_filt <- calcNormFactors(dge_filt, method = "TMM")
  dge_filt <- estimateDisp(dge_filt, design, robust = TRUE)
  fit <- glmQLFit(dge_filt, design, robust = TRUE)
  
  contrast_uf <- makeContrasts(LV_UF - LV_NOT_UF, levels = design)
  qlf <- glmQLFTest(fit, contrast = contrast_uf)
  tab <- topTags(qlf, n = Inf)$table
  
  res <- list()
  for (g in hub_genes) {
    if (g %in% rownames(tab)) {
      res[[g]] <- data.frame(Gene = g, log2FC = tab[g, "logFC"], PValue = tab[g, "PValue"])
    } else {
      res[[g]] <- data.frame(Gene = g, log2FC = NA, PValue = NA)
    }
  }
  do.call(rbind, res)
}

cat("=== 1. Stromal with Mesothelial Fraction Covariate ===\n")
res_cov <- run_edger("data/processed/pseudobulk/pb_audit_stromal_all_counts.csv", 
                     "data/processed/pseudobulk/pb_audit_stromal_all_meta.csv", add_covariate=TRUE)
print(res_cov)

cat("\n=== 2. Mesothelial-Scored Cells Subset (>=20 cells/donor) ===\n")
res_meso <- run_edger("data/processed/pseudobulk/pb_audit_stromal_meso_counts.csv", 
                      "data/processed/pseudobulk/pb_audit_stromal_meso_meta.csv", add_covariate=FALSE)
print(res_meso)

cat("\n=== 3. Fibroblast-Scored Cells Subset (>=20 cells/donor) ===\n")
res_fibro <- run_edger("data/processed/pseudobulk/pb_audit_stromal_fibro_counts.csv", 
                       "data/processed/pseudobulk/pb_audit_stromal_fibro_meta.csv", add_covariate=FALSE)
print(res_fibro)
