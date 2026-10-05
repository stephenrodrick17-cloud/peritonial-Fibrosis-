# ==============================================================================
# Script 50: Stage 4 Final Unified edgeR Pipeline
# Strictly builds all primary and sensitivity models from pseudobulk inputs
# Moves legacy tables to results/legacy_stage4/
# Outputs verified edgeR results with full Student's t distribution metrics
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})

set.seed(42)

cat("================================================================================\n")
cat("SCRIPT 50: STAGE 4 UNIFIED edgeR PSEUDOBULK RECONSTRUCTION\n")
cat("================================================================================\n\n")

# ------------------------------------------------------------------------------
# 1. Archive legacy Stage 4 results
# ------------------------------------------------------------------------------
legacy_dir <- "results/legacy_stage4"
dir.create(legacy_dir, showWarnings = FALSE, recursive = TRUE)

old_stage4_files <- list.files("results/tables", pattern = "^stage4_.*\\.csv$", full.names = TRUE)
if (length(old_stage4_files) > 0) {
  cat(sprintf("Archiving %d legacy stage4 files to %s/ ...\n", length(old_stage4_files), legacy_dir))
  for (f in old_stage4_files) {
    dest <- file.path(legacy_dir, basename(f))
    file.rename(f, dest)
  }
} else {
  cat("No legacy stage4 files found to archive.\n")
}

# ------------------------------------------------------------------------------
# 2. Setup Hub Genes and Directories
# ------------------------------------------------------------------------------
hub_genes <- c("COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
               "EDIL3", "LOX", "INHBA", "ISM1", "COMP")

pb_dir <- "data/processed/pseudobulk"
out_dir <- "results/tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Target cell types
cell_types <- list(
  "stromal" = list(
    label = "stromal / mesothelial-lineage (unresolved)",
    counts = file.path(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv"),
    meta = file.path(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv")
  ),
  "monocyte" = list(
    label = "Monocyte / macrophage",
    counts = file.path(pb_dir, "pb_primary_Monocyte_macrophage_counts.csv"),
    meta = file.path(pb_dir, "pb_primary_Monocyte_macrophage_metadata.csv")
  ),
  "cdc" = list(
    label = "cDC",
    counts = file.path(pb_dir, "pb_primary_cDC_counts.csv"),
    meta = file.path(pb_dir, "pb_primary_cDC_metadata.csv")
  ),
  "t_cell" = list(
    label = "T cell",
    counts = file.path(pb_dir, "pb_primary_T_cell_counts.csv"),
    meta = file.path(pb_dir, "pb_primary_T_cell_metadata.csv")
  ),
  "nk_cell" = list(
    label = "NK cell",
    counts = file.path(pb_dir, "pb_primary_NK_cell_counts.csv"),
    meta = file.path(pb_dir, "pb_primary_NK_cell_metadata.csv")
  )
)

# ------------------------------------------------------------------------------
# 3. Core edgeR Runner Function
# ------------------------------------------------------------------------------
run_edger_de <- function(counts_file, meta_file, cell_type_label, analysis_name="primary", 
                         top_n_hvg=NULL, is_not_estimable=FALSE) {
  
  meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
  n_uf <- sum(meta$group == "LV_UF")
  n_not_uf <- sum(meta$group == "LV_NOT_UF")
  n_sv <- sum(meta$group == "SV")
  tot_cells <- sum(meta$n_cells)
  
  contrasts_to_run <- c("LV_UF_vs_LV_NOT_UF", "LV_UF_vs_SV", "LV_NOT_UF_vs_SV")
  
  # If explicitly not estimable (or insufficient donors in any primary group)
  if (is_not_estimable || n_uf < 2 || n_not_uf < 2 || n_sv < 2) {
    rows <- list()
    for (cname in contrasts_to_run) {
      for (g in hub_genes) {
        rows[[paste(cname, g, sep="__")]] <- data.frame(
          Analysis = analysis_name,
          Cell_Type = cell_type_label,
          Contrast = cname,
          Gene = g,
          log2FC = NA_real_,
          SE = NA_real_,
          CI_95_low = NA_real_,
          CI_95_high = NA_real_,
          t_stat = NA_real_,
          PValue = NA_real_,
          BH_FDR = NA_real_,
          QL_Denom_df = NA_real_,
          Donors_LV_UF = n_uf,
          Donors_LV_NOT_UF = n_not_uf,
          Donors_SV = n_sv,
          Cells_total = tot_cells,
          Evidence_Status = "NOT ESTIMABLE",
          stringsAsFactors = FALSE
        )
      }
    }
    return(do.call(rbind, rows))
  }
  
  dt <- data.table::fread(counts_file, data.table = FALSE)
  rownames(dt) <- dt[[1]]
  counts <- dt[, -1, drop = FALSE]
  
  common_donors <- intersect(rownames(meta), rownames(counts))
  meta <- meta[common_donors, ]
  counts <- counts[common_donors, ]
  
  mat <- t(as.matrix(counts))
  
  groups <- factor(meta$group, levels = c("LV_NOT_UF", "LV_UF", "SV"))
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  dge <- DGEList(counts = mat, group = groups)
  
  # Filter low expression
  keep <- filterByExpr(dge, design)
  hub_in_mat <- intersect(hub_genes, rownames(mat))
  keep[hub_in_mat] <- TRUE
  dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalization: top-2000 HVGs vs full TMM
  if (!is.null(top_n_hvg)) {
    log_cpm <- cpm(dge_filt, log = TRUE)
    vars <- apply(log_cpm, 1, var)
    top_genes <- names(sort(vars, decreasing = TRUE))[1:min(top_n_hvg, length(vars))]
    top_genes <- unique(c(top_genes, hub_in_mat))
    dge_tmm <- dge_filt[top_genes, , keep.lib.sizes = FALSE]
    dge_tmm <- calcNormFactors(dge_tmm, method = "TMM")
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
  
  # Descriptive evidence check: >=2 donors in each group with >= 1 count
  pass_evidence <- function(gene_name) {
    if (!gene_name %in% rownames(mat)) return(FALSE)
    row_counts <- mat[gene_name, ]
    uf_expr <- sum(row_counts[meta$group == "LV_UF"] > 0)
    not_uf_expr <- sum(row_counts[meta$group == "LV_NOT_UF"] > 0)
    sv_expr <- sum(row_counts[meta$group == "SV"] > 0)
    return(uf_expr >= 2 && not_uf_expr >= 2 && sv_expr >= 2)
  }
  
  results_list <- list()
  for (cname in names(contrasts_list)) {
    qlf <- glmQLFTest(fit, contrast = contrasts_list[[cname]])
    tab <- topTags(qlf, n = Inf)$table
    
    for (g in hub_genes) {
      ev_status <- if (pass_evidence(g)) "PASS" else "INSUFFICIENT DATA"
      
      if (g %in% rownames(tab)) {
        row <- tab[g, ]
        lfc <- as.numeric(row$logFC)
        p_val <- as.numeric(row$PValue)
        f_stat <- as.numeric(row$F)
        
        # Denominator degrees of freedom from empirical Bayes shrinkage in qlf$df.total
        gene_idx <- which(rownames(fit$counts) == g)
        df_denom <- if (length(qlf$df.total) == 1) as.numeric(qlf$df.total) else as.numeric(qlf$df.total[gene_idx])
        
        # Standard error and t-statistic
        se <- if (f_stat > 0) abs(lfc) / sqrt(f_stat) else NA_real_
        t_stat <- if (!is.na(se) && se > 0) lfc / se else NA_real_
        
        # Student's t critical value
        t_crit <- qt(0.975, df = df_denom)
        ci_lo <- lfc - t_crit * se
        ci_hi <- lfc + t_crit * se
      } else {
        lfc <- NA_real_; p_val <- NA_real_; se <- NA_real_
        ci_lo <- NA_real_; ci_hi <- NA_real_; t_stat <- NA_real_
        df_denom <- NA_real_
      }
      
      results_list[[paste(cname, g, sep = "__")]] <- data.frame(
        Analysis = analysis_name,
        Cell_Type = cell_type_label,
        Contrast = cname,
        Gene = g,
        log2FC = lfc,
        SE = se,
        CI_95_low = ci_lo,
        CI_95_high = ci_hi,
        t_stat = t_stat,
        PValue = p_val,
        BH_FDR = NA_real_, # assigned globally across 96 tests
        QL_Denom_df = df_denom,
        Donors_LV_UF = n_uf,
        Donors_LV_NOT_UF = n_not_uf,
        Donors_SV = n_sv,
        Cells_total = tot_cells,
        Evidence_Status = ev_status,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, results_list)
}

# ------------------------------------------------------------------------------
# 4. Run Primary Models (5 Cell Types x 3 Contrasts x 11 Genes)
# ------------------------------------------------------------------------------
cat("--- Running Primary edgeR Models across 5 Cell Types ---\n")
primary_list <- list()
for (ct_key in names(cell_types)) {
  ct_info <- cell_types[[ct_key]]
  cat(sprintf("  Processing %s ...\n", ct_info$label))
  res <- run_edger_de(ct_info$counts, ct_info$meta, ct_info$label, analysis_name = "primary")
  primary_list[[ct_key]] <- res
}
df_primary <- do.call(rbind, primary_list)
rownames(df_primary) <- NULL

# Compute BH-FDR across all m=96 PASS tests in primary analysis
pass_mask <- df_primary$Evidence_Status == "PASS" & !is.na(df_primary$PValue)
m_family <- sum(pass_mask)
cat(sprintf("Total tests in primary analysis: %d | PASS evidence rule (m): %d\n", nrow(df_primary), m_family))

df_primary$BH_FDR[pass_mask] <- p.adjust(df_primary$PValue[pass_mask], method = "BH")

primary_csv <- file.path(out_dir, "stage4_primary_edger_pseudobulk.csv")
write.csv(df_primary, primary_csv, row.names = FALSE)
cat(sprintf("Saved Primary Results to: %s (%d rows)\n\n", primary_csv, nrow(df_primary)))

# ------------------------------------------------------------------------------
# 5. Run Sensitivity Models for Stromal Cells
# ------------------------------------------------------------------------------
cat("--- Running Stromal Sensitivity Models (b)-(j) + LODO ---\n")
st_label <- "stromal / mesothelial-lineage (unresolved)"

sens_configs <- list(
  # (b) drop LV_UF-3
  "sens_b_no_LV_UF3" = list(
    counts = file.path(pb_dir, "pb_sens_b_no_LV_UF3_counts.csv"),
    meta = file.path(pb_dir, "pb_sens_b_no_LV_UF3_metadata.csv"),
    hvg = NULL, not_est = FALSE
  ),
  # (c) drop cells flagged by Scrublet or scDblFinder
  "sens_c_no_scDblFinder" = list(
    counts = file.path(pb_dir, "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_counts.csv"),
    meta = file.path(pb_dir, "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_metadata.csv"),
    hvg = NULL, not_est = FALSE
  ),
  # (d) stromal cells with n_genes >= 700
  "sens_d_ge700genes" = list(
    counts = file.path(pb_dir, "pb_sens_d_ge700genes_counts.csv"),
    meta = file.path(pb_dir, "pb_sens_d_ge700genes_metadata.csv"),
    hvg = NULL, not_est = FALSE
  ),
  # (e) donors with >=50 stromal cells (3 LV_UF / 2 LV_NOT_UF / 5 SV)
  "sens_e_ge50cells" = list(
    counts = file.path(pb_dir, "pb_sens_e_ge50cells_counts.csv"),
    meta = file.path(pb_dir, "pb_sens_e_ge50cells_metadata.csv"),
    hvg = NULL, not_est = FALSE
  ),
  # (f) TMM computed on the top-2000 variable genes
  "sens_f_top2000_hvg_tmm" = list(
    counts = cell_types$stromal$counts,
    meta = cell_types$stromal$meta,
    hvg = 2000, not_est = FALSE
  ),
  # (h) donors with >=50 cells AND median genes >=1000 (one LV_UF donor -> NOT ESTIMABLE)
  "sens_h_ge50cells_med1000genes" = list(
    counts = file.path(pb_dir, "pb_sens_h_ge50cells_med1000genes_counts.csv"),
    meta = file.path(pb_dir, "pb_sens_h_ge50cells_med1000genes_metadata.csv"),
    hvg = NULL, not_est = TRUE
  ),
  # (i) no upper ceilings, remove Scrublet+scDblFinder calls
  "sens_i_no_upper_ceilings" = list(
    counts = file.path(pb_dir, "pb_sens_i_no_upper_ceilings_counts.csv"),
    meta = file.path(pb_dir, "pb_sens_i_no_upper_ceilings_metadata.csv"),
    hvg = NULL, not_est = FALSE
  ),
  # (j) one pooled 99.5th-percentile ceiling
  "sens_j_pooled_fixed_ceiling" = list(
    counts = file.path(pb_dir, "pb_sens_j_pooled_fixed_ceiling_counts.csv"),
    meta = file.path(pb_dir, "pb_sens_j_pooled_fixed_ceiling_metadata.csv"),
    hvg = NULL, not_est = FALSE
  )
)

sens_results_list <- list()

for (s_name in names(sens_configs)) {
  cfg <- sens_configs[[s_name]]
  cat(sprintf("  Running %s ...\n", s_name))
  res <- run_edger_de(cfg$counts, cfg$meta, st_label, analysis_name = s_name, 
                      top_n_hvg = cfg$hvg, is_not_estimable = cfg$not_est)
  sens_results_list[[s_name]] <- res
}

# ------------------------------------------------------------------------------
# 6. Run Leave-One-Donor-Out (LODO) across all 16 Donors
# ------------------------------------------------------------------------------
cat("--- Running 16 Leave-One-Donor-Out (LODO) Iterations ---\n")
lodo_files <- list.files(pb_dir, pattern = "^pb_sens_lodo_st_.*_metadata\\.csv$", full.names = TRUE)
donors_dropped <- gsub("^pb_sens_lodo_st_|_metadata\\.csv$", "", basename(lodo_files))

for (donor_out in sort(donors_dropped)) {
  c_file <- file.path(pb_dir, sprintf("pb_sens_lodo_st_%s_counts.csv", donor_out))
  m_file <- file.path(pb_dir, sprintf("pb_sens_lodo_st_%s_metadata.csv", donor_out))
  analysis_label <- sprintf("lodo_out_%s", donor_out)
  cat(sprintf("  LODO dropping %s ...\n", donor_out))
  res_lodo <- run_edger_de(c_file, m_file, st_label, analysis_name = analysis_label)
  sens_results_list[[analysis_label]] <- res_lodo
}

df_sens <- do.call(rbind, sens_results_list)
rownames(df_sens) <- NULL

sens_csv <- file.path(out_dir, "stage4_sensitivity_edger_pseudobulk.csv")
write.csv(df_sens, sens_csv, row.names = FALSE)
cat(sprintf("Saved Sensitivity Results to: %s (%d rows)\n\n", sens_csv, nrow(df_sens)))

# ------------------------------------------------------------------------------
# 7. Print Session Info and Call Signatures
# ------------------------------------------------------------------------------
cat("================================================================================\n")
cat("R SESSION INFORMATION & edgeR CODE CALLS\n")
cat("================================================================================\n")
print(sessionInfo())

cat("\nEXACT edgeR CALLS USED IN PIPELINE:\n")
cat("  1. dge <- DGEList(counts = mat, group = groups)\n")
cat("  2. keep <- filterByExpr(dge, design); keep[hub_genes] <- TRUE; dge_filt <- dge[keep, , keep.lib.sizes=FALSE]\n")
cat("  3. dge_filt <- calcNormFactors(dge_filt, method = 'TMM')\n")
cat("  4. dge_filt <- estimateDisp(dge_filt, design, robust = TRUE)\n")
cat("  5. fit <- glmQLFit(dge_filt, design, robust = TRUE)\n")
cat("  6. qlf <- glmQLFTest(fit, contrast = makeContrasts(contrast_str, levels = design))\n")
cat("  7. t_stat <- sign(log2FC) * sqrt(F_stat); SE <- abs(log2FC) / sqrt(F_stat)\n")
cat("  8. CI_95 <- log2FC +/- qt(0.975, df = fit$df.total) * SE\n")
cat("================================================================================\n")
cat("Script 50 execution completed successfully.\n")
