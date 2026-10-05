# Script 44: Stage 4 Re-run of Sensitivities (b)-(f) and Leave-One-Donor-Out (LODO)
# Rules:
# - Re-run edgeR quasi-likelihood across ALL 37,487 genes (filterByExpr + TMM)
# - Compare TMM normalization factors for (f) vs primary to explain identical fold changes
# - Save to results/tables/stage4_sensitivity_edger_pseudobulk.csv
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
cat("STAGE 4: RE-RUNNING SENSITIVITY ANALYSES (b)-(f) AND LODO IN edgeR\n")
cat("================================================================================\n")

run_edger_de_detailed <- function(counts_file, meta_file, cell_type_label, analysis_name="sens", top_n_hvg=NULL) {
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
  
  # Check if top-N HVG normalization is requested
  if (!is.null(top_n_hvg)) {
    log_cpm <- cpm(dge_filt, log = TRUE)
    vars <- apply(log_cpm, 1, var)
    top_genes <- names(sort(vars, decreasing = TRUE))[1:min(top_n_hvg, length(vars))]
    top_genes <- unique(c(top_genes, hub_in_mat))
    
    cat(sprintf("  [Analysis %s] Normalization subset: %d top HVGs selected out of %d expressed genes.\n", 
                analysis_name, length(top_genes), nrow(dge_filt)))
    
    dge_tmm <- dge_filt[top_genes, , keep.lib.sizes = FALSE]
    dge_tmm <- calcNormFactors(dge_tmm, method = "TMM")
    
    # Standard TMM for comparison
    dge_std <- calcNormFactors(dge_filt, method = "TMM")
    cat("  [Analysis Comparison] TMM norm factors (top 2000 HVG vs All Genes):\n")
    comp_nf <- data.frame(
      Donor = rownames(dge_filt$samples),
      NormFactor_Top2000HVG = round(dge_tmm$samples$norm.factors, 6),
      NormFactor_AllGenes = round(dge_std$samples$norm.factors, 6),
      AbsDiff = round(abs(dge_tmm$samples$norm.factors - dge_std$samples$norm.factors), 6)
    )
    print(comp_nf)
    
    dge_filt$samples$norm.factors <- dge_tmm$samples$norm.factors
  } else {
    dge_filt <- calcNormFactors(dge_filt, method = "TMM")
  }
  
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
    
    # Extract denominator df for QL F-test
    df_denom <- qlf$df.total
    if (length(df_denom) > 1) {
      df_denom_scalar <- df_denom[1]
    } else {
      df_denom_scalar <- df_denom
    }
    
    for (g in hub_genes) {
      if (g %in% rownames(tab)) {
        row <- tab[g, ]
        lfc <- row$logFC
        p_val <- row$PValue
        f_stat <- row$F
        se <- if (f_stat > 0) abs(lfc) / sqrt(f_stat) else NA
        # Student's t critical value
        t_crit <- qt(0.975, df = df_denom_scalar)
        ci_lo <- lfc - t_crit * se
        ci_hi <- lfc + t_crit * se
      } else {
        lfc <- NA; p_val <- NA; se <- NA; ci_lo <- NA; ci_hi <- NA
      }
      
      results_list[[paste(cname, g, sep = "__")]] <- data.frame(
        Analysis = analysis_name,
        Cell_Type = cell_type_label,
        Contrast = cname,
        Gene = g,
        log2FC = round(lfc, 4),
        SE = round(se, 4),
        CI_95_low = round(ci_lo, 4),
        CI_95_high = round(ci_hi, 4),
        PValue = p_val,
        QL_Denom_df = round(df_denom_scalar, 2),
        Donors_LV_UF = n_uf,
        Donors_LV_NOT_UF = n_not_uf,
        Donors_SV = n_sv,
        Cells_total = sum(meta$n_cells),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, results_list)
}

sens_results_list <- list()

# (b) stromal without LV_UF-3
cat("\nRunning Sensitivity (b): Stromal without LV_UF-3...\n")
res_b <- run_edger_de_detailed(file.path(pb_dir, "pb_sens_b_no_LV_UF3_counts.csv"),
                               file.path(pb_dir, "pb_sens_b_no_LV_UF3_metadata.csv"),
                               "stromal / mesothelial-lineage (unresolved)",
                               analysis_name = "sens_b_no_LV_UF3")
sens_results_list[["sens_b"]] <- res_b

# (c) after removing scDblFinder doublets
cat("\nRunning Sensitivity (c): Stromal after removing scDblFinder doublets...\n")
res_c <- run_edger_de_detailed(file.path(pb_dir, "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_counts.csv"),
                               file.path(pb_dir, "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_metadata.csv"),
                               "stromal / mesothelial-lineage (unresolved)",
                               analysis_name = "sens_c_no_scDblFinder")
sens_results_list[["sens_c"]] <- res_c

# (d) stromal using only cells with >= 700 genes
cat("\nRunning Sensitivity (d): Stromal with genes >= 700...\n")
res_d <- run_edger_de_detailed(file.path(pb_dir, "pb_sens_d_ge700genes_counts.csv"),
                               file.path(pb_dir, "pb_sens_d_ge700genes_metadata.csv"),
                               "stromal / mesothelial-lineage (unresolved)",
                               analysis_name = "sens_d_ge700genes")
sens_results_list[["sens_d"]] <- res_d

# (e) stromal pseudobulk only for donors with >= 50 stromal cells
cat("\nRunning Sensitivity (e): Stromal donors with >= 50 cells...\n")
res_e <- run_edger_de_detailed(file.path(pb_dir, "pb_sens_e_ge50cells_counts.csv"),
                               file.path(pb_dir, "pb_sens_e_ge50cells_metadata.csv"),
                               "stromal / mesothelial-lineage (unresolved)",
                               analysis_name = "sens_e_ge50cells")
sens_results_list[["sens_e"]] <- res_e

# (f) top 2000 most variable genes for TMM
cat("\nRunning Sensitivity (f): Stromal with top 2000 HVGs for TMM...\n")
res_f <- run_edger_de_detailed(file.path(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv"),
                               file.path(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv"),
                               "stromal / mesothelial-lineage (unresolved)",
                               analysis_name = "sens_f_top2000_hvg_tmm",
                               top_n_hvg = 2000)
sens_results_list[["sens_f"]] <- res_f

# (a) Leave-One-Donor-Out (LODO) for all 16 donors
all_donors <- c("LV_NOT_UF-1", "LV_NOT_UF-2", "LV_NOT_UF-3", "LV_NOT_UF-4", "LV_NOT_UF-5", "LV_NOT_UF-6",
                "LV_UF-1", "LV_UF-2", "LV_UF-3", "LV_UF-4",
                "SV-1", "SV-2", "SV-3", "SV-4", "SV-5", "SV-6")

cat("\nRunning Leave-One-Donor-Out (16 iterations)...\n")
for (d_out in all_donors) {
  cfile <- file.path(pb_dir, paste0("pb_sens_lodo_st_", d_out, "_counts.csv"))
  mfile <- file.path(pb_dir, paste0("pb_sens_lodo_st_", d_out, "_metadata.csv"))
  res_lodo <- run_edger_de_detailed(cfile, mfile,
                                    "stromal / mesothelial-lineage (unresolved)",
                                    analysis_name = paste0("lodo_out_", d_out))
  sens_results_list[[paste0("lodo_", d_out)]] <- res_lodo
}

df_all_sens <- do.call(rbind, sens_results_list)
out_csv <- "results/tables/stage4_sensitivity_edger_pseudobulk.csv"
write.csv(df_all_sens, out_csv, row.names = FALSE)
cat(sprintf("\nAll sensitivity analyses completed and saved to: %s\n", out_csv))
cat(sprintf("Total rows written: %d\n", nrow(df_all_sens)))
cat("================================================================================\n")
