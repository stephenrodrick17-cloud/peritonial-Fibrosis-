# Script 42: edgeR Differential Testing for Sensitivities (g) and (h)
# Rules:
# - Run edgeR quasi-likelihood across ALL 37,487 genes (filterByExpr + TMM)
# - Contrasts: LV_UF_vs_LV_NOT_UF, LV_UF_vs_SV, LV_NOT_UF_vs_SV
# - Compute log2FC, SE, 95% CI, raw P, BH-FDR
# - Compare signs against PRIMARY stromal results
# - Seed 42, no causal language.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})

set.seed(42)

pb_dir <- "data/processed/pseudobulk"
hub_genes <- c("COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
               "EDIL3", "LOX", "INHBA", "ISM1", "COMP")

cat("================================================================================\n")
cat("STAGE 4: edgeR PSEUDOBULK FOR SENSITIVITY (g) & (h)\n")
cat("================================================================================\n")

run_edger_de_gh <- function(counts_file, meta_file, analysis_name) {
  dt <- data.table::fread(counts_file, data.table = FALSE)
  rownames(dt) <- dt[[1]]
  counts <- dt[, -1, drop = FALSE]
  meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
  
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
    LV_UF_vs_SV = makeContrasts(LV_UF - SV, levels = design),
    LV_NOT_UF_vs_SV = makeContrasts(LV_NOT_UF - SV, levels = design)
  )
  
  results_list <- list()
  for (cname in names(contrasts_list)) {
    qlf <- glmQLFTest(fit, contrast = contrasts_list[[cname]])
    tab <- topTags(qlf, n = Inf)$table
    
    n_uf <- sum(meta$group == "LV_UF")
    n_not_uf <- sum(meta$group == "LV_NOT_UF")
    n_sv <- sum(meta$group == "SV")
    
    for (g in hub_genes) {
      if (g %in% rownames(tab)) {
        row <- tab[g, ]
        lfc <- row$logFC
        p_val <- row$PValue
        f_stat <- row$F
        se <- if (f_stat > 0) abs(lfc) / sqrt(f_stat) else NA
        ci_lo <- lfc - 1.96 * se
        ci_hi <- lfc + 1.96 * se
      } else {
        lfc <- NA; p_val <- NA; se <- NA; ci_lo <- NA; ci_hi <- NA
      }
      
      results_list[[paste(cname, g, sep = "__")]] <- data.frame(
        Analysis = analysis_name,
        Cell_Type = "stromal / mesothelial-lineage (unresolved)",
        Contrast = cname,
        Gene = g,
        log2FC = round(lfc, 4),
        SE = round(se, 4),
        CI_95_low = round(ci_lo, 4),
        CI_95_high = round(ci_hi, 4),
        PValue = p_val,
        Donors_LV_UF = n_uf,
        Donors_LV_NOT_UF = n_not_uf,
        Donors_SV = n_sv,
        Cells_total = sum(meta$n_cells),
        stringsAsFactors = FALSE
      )
    }
  }
  res_df <- do.call(rbind, results_list)
  # BH-FDR within the analysis across tested hub genes
  res_df$BH_FDR <- p.adjust(res_df$PValue, method = "BH")
  res_df
}

# Run Sensitivity (g)
cat("Running edgeR for Sensitivity (g): Ceiling-Rescued Stromal...\n")
res_g <- run_edger_de_gh(
  file.path(pb_dir, "pb_sens_g_ceiling_rescued_stromal_counts.csv"),
  file.path(pb_dir, "pb_sens_g_ceiling_rescued_stromal_metadata.csv"),
  "sens_g_ceiling_rescued"
)

# Run Sensitivity (h)
cat("Running edgeR for Sensitivity (h): High-Confidence Donors (>=50 cells & med genes >=1000)...\n")
res_h <- run_edger_de_gh(
  file.path(pb_dir, "pb_sens_h_ge50cells_med1000genes_counts.csv"),
  file.path(pb_dir, "pb_sens_h_ge50cells_med1000genes_metadata.csv"),
  "sens_h_ge50cells_med1000genes"
)

df_gh <- rbind(res_g, res_h)
write.csv(df_gh, "results/tables/stage4_sensitivity_g_h_edger_pseudobulk.csv", row.names = FALSE)
cat("Saved sensitivity (g) and (h) results to results/tables/stage4_sensitivity_g_h_edger_pseudobulk.csv\n\n")

# Load Primary Stromal Results for Comparison
df_primary <- read.csv("results/tables/stage4_primary_edger_pseudobulk.csv", stringsAsFactors = FALSE)
df_prim_st <- df_primary[df_primary$Cell_Type == "stromal / mesothelial-lineage (unresolved)" & df_primary$Contrast == "LV_UF_vs_LV_NOT_UF", ]

cat("================================================================================\n")
cat(">>> SENSITIVITY (g) RESULTS: Ceiling-Rescued Stromal (All 3 Contrasts)\n")
cat("================================================================================\n")
print(res_g[, c("Contrast", "Gene", "log2FC", "SE", "PValue", "BH_FDR", "Donors_LV_UF", "Donors_LV_NOT_UF", "Donors_SV", "Cells_total")], row.names = FALSE)

cat("\n================================================================================\n")
cat(">>> SENSITIVITY (h) RESULTS: High-Confidence Donors (All 3 Contrasts)\n")
cat("================================================================================\n")
print(res_h[, c("Contrast", "Gene", "log2FC", "SE", "PValue", "BH_FDR", "Donors_LV_UF", "Donors_LV_NOT_UF", "Donors_SV", "Cells_total")], row.names = FALSE)

cat("\n================================================================================\n")
cat(">>> SIGN COMPARISON: PRIMARY vs SENSITIVITIES (g) and (h) [LV_UF vs LV_NOT_UF]\n")
cat("================================================================================\n")

comp_list <- list()
for (g in hub_genes) {
  p_row <- df_prim_st[df_prim_st$Gene == g, ]
  g_row <- res_g[res_g$Contrast == "LV_UF_vs_LV_NOT_UF" & res_g$Gene == g, ]
  h_row <- res_h[res_h$Contrast == "LV_UF_vs_LV_NOT_UF" & res_h$Gene == g, ]
  
  p_lfc <- if (nrow(p_row) > 0) p_row$log2FC else NA
  p_pval <- if (nrow(p_row) > 0) p_row$PValue else NA
  p_fdr <- if (nrow(p_row) > 0) p_row$BH_FDR else NA
  p_ev <- if (nrow(p_row) > 0) p_row$Evidence_Status else NA
  
  g_lfc <- if (nrow(g_row) > 0) g_row$log2FC else NA
  g_pval <- if (nrow(g_row) > 0) g_row$PValue else NA
  g_fdr <- if (nrow(g_row) > 0) g_row$BH_FDR else NA
  
  h_lfc <- if (nrow(h_row) > 0) h_row$log2FC else NA
  h_pval <- if (nrow(h_row) > 0) h_row$PValue else NA
  h_fdr <- if (nrow(h_row) > 0) h_row$BH_FDR else NA
  
  sign_p <- ifelse(is.na(p_lfc), "NA", ifelse(p_lfc > 0, "+", "-"))
  sign_g <- ifelse(is.na(g_lfc), "NA", ifelse(g_lfc > 0, "+", "-"))
  sign_h <- ifelse(is.na(h_lfc), "NA", ifelse(h_lfc > 0, "+", "-"))
  
  match_g <- ifelse(sign_p == sign_g, "MATCH", "OPPOSITE")
  match_h <- ifelse(sign_p == sign_h, "MATCH", "OPPOSITE")
  
  comp_list[[g]] <- data.frame(
    Gene = g,
    Evidence = p_ev,
    PRIMARY_log2FC = p_lfc,
    PRIMARY_P = p_pval,
    PRIMARY_FDR = p_fdr,
    SENS_g_log2FC = g_lfc,
    SENS_g_P = g_pval,
    SENS_g_FDR = g_fdr,
    Sign_g_vs_Primary = match_g,
    SENS_h_log2FC = h_lfc,
    SENS_h_P = h_pval,
    SENS_h_FDR = h_fdr,
    Sign_h_vs_Primary = match_h,
    stringsAsFactors = FALSE
  )
}

df_comp <- do.call(rbind, comp_list)
write.csv(df_comp, "results/tables/stage4_sensitivity_g_h_comparison_with_primary.csv", row.names = FALSE)
print(df_comp, row.names = FALSE)

cat("\nAnalysis of sensitivities (g) and (h) completed successfully.\n")
cat("================================================================================\n")
