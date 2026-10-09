
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})

set.seed(42)

pb_dir <- "data/processed/pseudobulk"
hub_genes <- c("COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
               "EDIL3", "LOX", "INHBA", "ISM1", "COMP")

run_de <- function(cfile, mfile, ana_name) {
  dt <- data.table::fread(cfile, data.table = FALSE)
  rownames(dt) <- dt[[1]]
  counts <- dt[, -1, drop = FALSE]
  meta <- read.csv(mfile, row.names = 1, stringsAsFactors = FALSE)
  
  common_donors <- intersect(rownames(meta), rownames(counts))
  meta <- meta[common_donors, ]
  counts <- counts[common_donors, ]
  
  mat <- t(as.matrix(counts))
  groups <- factor(meta$group, levels = c("LV_NOT_UF", "LV_UF", "SV"))
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  dge <- DGEList(counts = mat, group = groups)
  keep <- filterByExpr(dge, design)
  hub_in_mat <- intersect(hub_genes, rownames(mat))
  keep[hub_in_mat] <- TRUE
  
  dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
  dge_filt <- calcNormFactors(dge_filt, method = "TMM")
  dge_filt <- estimateDisp(dge_filt, design, robust = TRUE)
  fit <- glmQLFit(dge_filt, design, robust = TRUE)
  
  contrasts_list <- list(
    LV_UF_vs_LV_NOT_UF = makeContrasts(LV_UF - LV_NOT_UF, levels = design),
    LV_UF_vs_SV = makeContrasts(LV_UF - SV, levels = design)
  )
  
  res_list <- list()
  for (cname in names(contrasts_list)) {
    qlf <- glmQLFTest(fit, contrast = contrasts_list[[cname]])
    tab <- topTags(qlf, n = Inf)$table
    df_denom <- qlf$df.total[1]
    t_crit <- qt(0.975, df = df_denom)
    
    n_uf <- sum(meta$group == "LV_UF")
    n_not_uf <- sum(meta$group == "LV_NOT_UF")
    n_sv <- sum(meta$group == "SV")
    
    for (g in hub_genes) {
      if (g %in% rownames(tab)) {
        row <- tab[g, ]
        lfc <- row$logFC
        pval <- row$PValue
        f_stat <- row$F
        se <- if (f_stat > 0) abs(lfc) / sqrt(f_stat) else NA
        ci_lo <- lfc - t_crit * se
        ci_hi <- lfc + t_crit * se
      } else {
        lfc <- NA; pval <- NA; se <- NA; ci_lo <- NA; ci_hi <- NA
      }
      res_list[[paste(cname, g, sep = "__")]] <- data.frame(
        Analysis = ana_name,
        Cell_Type = "stromal / mesothelial-lineage (unresolved)",
        Contrast = cname,
        Gene = g,
        log2FC = round(lfc, 4),
        SE = round(se, 4),
        CI_95_low = round(ci_lo, 4),
        CI_95_high = round(ci_hi, 4),
        PValue = pval,
        QL_Denom_df = round(df_denom, 2),
        Donors_LV_UF = n_uf,
        Donors_LV_NOT_UF = n_not_uf,
        Donors_SV = n_sv,
        Cells_total = sum(meta$n_cells),
        stringsAsFactors = FALSE
      )
    }
  }
  res_df <- do.call(rbind, res_list)
  res_df$BH_FDR <- p.adjust(res_df$PValue, method = "BH")
  res_df
}

res_i <- run_de("data/processed/pseudobulk/pb_sens_i_no_upper_ceilings_counts.csv", "data/processed/pseudobulk/pb_sens_i_no_upper_ceilings_metadata.csv", "sens_i_no_upper_ceilings")
res_j <- run_de("data/processed/pseudobulk/pb_sens_j_pooled_fixed_ceiling_counts.csv", "data/processed/pseudobulk/pb_sens_j_pooled_fixed_ceiling_metadata.csv", "sens_j_pooled_fixed_ceiling")

res_ij <- rbind(res_i, res_j)
write.csv(res_ij, "results/tables/stage4_sensitivity_i_j_edger_pseudobulk.csv", row.names = FALSE)
cat("Sensitivity (i) and (j) edgeR differential testing complete.\n")
