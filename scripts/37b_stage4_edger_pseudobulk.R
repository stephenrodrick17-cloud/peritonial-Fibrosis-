# Script 37b: Stage 4 Step 3 (Primary edgeR Pseudobulk) & Step 4 (Sensitivity Analyses)
# Uses edgeR quasi-likelihood (filterByExpr, TMM over ALL 37,487 genes)
# Seed 42, no causal language.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

set.seed(42)

pb_dir <- "data/processed/pseudobulk"
hub_genes <- c("COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
               "EDIL3", "LOX", "INHBA", "ISM1", "COMP")

# Evidence rule: pairs with < 2 donors in ANY group with counts are INSUFFICIENT DATA
# Load descriptive summary to get evidence status
desc_df <- read.csv("results/tables/stage4_hub_genes_descriptive_summary.csv", stringsAsFactors = FALSE)

cat("================================================================================\n")
cat("STAGE 4 - STEP 3: PRIMARY PSEUDOBULK (edgeR Quasi-Likelihood over ALL genes)\n")
cat("================================================================================\n")

run_edger_de <- function(counts_file, meta_file, cell_type_label, analysis_name="primary", top_n_hvg=NULL) {
  dt <- data.table::fread(counts_file, data.table = FALSE)
  rownames(dt) <- dt[[1]]
  counts <- dt[, -1, drop = FALSE]
  meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
  
  # Ensure matching donors
  common_donors <- intersect(rownames(meta), rownames(counts))
  meta <- meta[common_donors, ]
  counts <- counts[common_donors, ]
  
  # Transpose: genes x donors
  mat <- t(as.matrix(counts))
  
  # Group factor
  groups <- factor(meta$group, levels = c("LV_NOT_UF", "LV_UF", "SV"))
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  dge <- DGEList(counts = mat, group = groups)
  
  # Filter by expression across all genes
  keep <- filterByExpr(dge, design)
  # Ensure all hub genes are retained for inspection even if lowly expressed
  hub_in_mat <- intersect(hub_genes, rownames(mat))
  keep[hub_in_mat] <- TRUE
  
  dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalization: top HVG check (f) or standard TMM over all genes
  if (!is.null(top_n_hvg)) {
    # Compute gene variances on log-CPM
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
  
  results_list <- list()
  for (cname in names(contrasts_list)) {
    qlf <- glmQLFTest(fit, contrast = contrasts_list[[cname]])
    tab <- topTags(qlf, n = Inf)$table
    
    # Donors per group
    n_uf <- sum(meta$group == "LV_UF")
    n_not_uf <- sum(meta$group == "LV_NOT_UF")
    n_sv <- sum(meta$group == "SV")
    
    for (g in hub_genes) {
      if (g %in% rownames(tab)) {
        row <- tab[g, ]
        lfc <- row$logFC
        p_val <- row$PValue
        f_stat <- row$F
        # Standard error: SE = abs(logFC) / sqrt(F)
        se <- if (f_stat > 0) abs(lfc) / sqrt(f_stat) else NA
        ci_lo <- lfc - 1.96 * se
        ci_hi <- lfc + 1.96 * se
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

# 1. PRIMARY ANALYSIS
primary_cts <- list(
  "stromal / mesothelial-lineage (unresolved)" = "stromal_mesothelial-lineage_unresolved",
  "Monocyte / macrophage" = "Monocyte_macrophage",
  "cDC" = "cDC",
  "T cell" = "T_cell",
  "NK cell" = "NK_cell"
)

primary_results_list <- list()
for (ct in names(primary_cts)) {
  clean_name <- primary_cts[[ct]]
  cfile <- file.path(pb_dir, paste0("pb_primary_", clean_name, "_counts.csv"))
  mfile <- file.path(pb_dir, paste0("pb_primary_", clean_name, "_metadata.csv"))
  res <- run_edger_de(cfile, mfile, ct, analysis_name = "primary")
  primary_results_list[[ct]] <- res
}
df_primary <- do.call(rbind, primary_results_list)

# Match Evidence Status from Descriptive Profile
df_primary$Evidence_Status <- "PASS"
for (i in seq_len(nrow(df_primary))) {
  g <- df_primary$Gene[i]
  ct <- df_primary$Cell_Type[i]
  match_ev <- desc_df$Evidence_Status[desc_df$Gene == g & desc_df$Cell_Type == ct]
  if (length(match_ev) > 0 && match_ev[1] == "INSUFFICIENT DATA") {
    df_primary$Evidence_Status[i] <- "INSUFFICIENT DATA"
  }
}

# Calculate BH-FDR across all gene x cell type x contrast tests that pass evidence rule
pass_evidence_idx <- which(df_primary$Evidence_Status == "PASS" & !is.na(df_primary$PValue))
df_primary$BH_FDR <- NA
df_primary$BH_FDR[pass_evidence_idx] <- p.adjust(df_primary$PValue[pass_evidence_idx], method = "BH")

write.csv(df_primary, "results/tables/stage4_primary_edger_pseudobulk.csv", row.names = FALSE)
cat("Primary edgeR pseudobulk analysis complete! Saved to results/tables/stage4_primary_edger_pseudobulk.csv\n\n")

# Print primary stromal results verbatim
cat("--- PRIMARY ANALYSIS RESULTS: Stromal / Mesothelial-Lineage (Primary contrast: LV_UF vs LV_NOT_UF) ---\n")
st_prim <- df_primary[df_primary$Cell_Type == "stromal / mesothelial-lineage (unresolved)" & df_primary$Contrast == "LV_UF_vs_LV_NOT_UF", ]
print(st_prim[, c("Gene", "log2FC", "SE", "CI_95_low", "CI_95_high", "PValue", "BH_FDR", "Evidence_Status")], row.names = FALSE)

cat("\n================================================================================\n")
cat("STAGE 4 - STEP 4: SENSITIVITY ANALYSES (b, c, d, e, f, and Leave-One-Donor-Out)\n")
cat("================================================================================\n")

sens_results_list <- list()

# (b) stromal without LV_UF-3
cat("Running Sensitivity (b): Stromal without LV_UF-3...\n")
res_b <- run_edger_de(file.path(pb_dir, "pb_sens_b_no_LV_UF3_counts.csv"),
                      file.path(pb_dir, "pb_sens_b_no_LV_UF3_metadata.csv"),
                      "stromal / mesothelial-lineage (unresolved)",
                      analysis_name = "sens_b_no_LV_UF3")
sens_results_list[["sens_b"]] <- res_b

# (c) after removing scDblFinder doublets
cat("Running Sensitivity (c): Stromal after removing scDblFinder doublets...\n")
res_c <- run_edger_de(file.path(pb_dir, "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_counts.csv"),
                      file.path(pb_dir, "pb_sens_c_no_scDblFinder_stromal_mesothelial-lineage_unresolved_metadata.csv"),
                      "stromal / mesothelial-lineage (unresolved)",
                      analysis_name = "sens_c_no_scDblFinder")
sens_results_list[["sens_c"]] <- res_c

# (d) stromal using only cells with >= 700 genes
cat("Running Sensitivity (d): Stromal with genes >= 700...\n")
res_d <- run_edger_de(file.path(pb_dir, "pb_sens_d_ge700genes_counts.csv"),
                      file.path(pb_dir, "pb_sens_d_ge700genes_metadata.csv"),
                      "stromal / mesothelial-lineage (unresolved)",
                      analysis_name = "sens_d_ge700genes")
sens_results_list[["sens_d"]] <- res_d

# (e) stromal pseudobulk only for donors with >= 50 stromal cells
cat("Running Sensitivity (e): Stromal donors with >= 50 cells...\n")
res_e <- run_edger_de(file.path(pb_dir, "pb_sens_e_ge50cells_counts.csv"),
                      file.path(pb_dir, "pb_sens_e_ge50cells_metadata.csv"),
                      "stromal / mesothelial-lineage (unresolved)",
                      analysis_name = "sens_e_ge50cells")
sens_results_list[["sens_e"]] <- res_e

# (f) top-N most variable genes for TMM (normalization check: top 2000 HVG)
cat("Running Sensitivity (f): TMM on top 2000 HVG...\n")
res_f <- run_edger_de(file.path(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv"),
                      file.path(pb_dir, "pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv"),
                      "stromal / mesothelial-lineage (unresolved)",
                      analysis_name = "sens_f_top2000_hvg_tmm",
                      top_n_hvg = 2000)
sens_results_list[["sens_f"]] <- res_f

# (a) Leave-one-donor-out (16 runs for stromal)
cat("Running Sensitivity (a): Leave-One-Donor-Out across all 16 donors...\n")
lodo_files <- list.files(pb_dir, pattern = "^pb_sens_lodo_st_.*_counts\\.csv$", full.names = TRUE)
for (lf in lodo_files) {
  d_out <- sub("_counts\\.csv$", "", sub("^.*pb_sens_lodo_st_", "", lf))
  mf <- file.path(pb_dir, paste0("pb_sens_lodo_st_", d_out, "_metadata.csv"))
  res_lodo <- run_edger_de(lf, mf, "stromal / mesothelial-lineage (unresolved)",
                           analysis_name = paste0("lodo_out_", d_out))
  sens_results_list[[paste0("lodo_", d_out)]] <- res_lodo
}

all_sens <- do.call(rbind, sens_results_list)
write.csv(all_sens, "results/tables/stage4_sensitivity_edger_pseudobulk.csv", row.names = FALSE)
cat("All sensitivity edgeR pseudobulk runs complete! Saved to results/tables/stage4_sensitivity_edger_pseudobulk.csv\n\n")

# ---------------------------------------------------------------------------
# Robustness Rule Assessment for Primary Contrast (LV_UF vs LV_NOT_UF)
# ---------------------------------------------------------------------------
cat("================================================================================\n")
cat("ROBUSTNESS RULE VERDICT (Stromal: Primary contrast LV_UF vs LV_NOT_UF)\n")
cat("================================================================================\n")

# Filter primary stromal
st_prim_uf <- df_primary[df_primary$Cell_Type == "stromal / mesothelial-lineage (unresolved)" & df_primary$Contrast == "LV_UF_vs_LV_NOT_UF", ]

robust_rows <- list()
for (g in hub_genes) {
  p_row <- st_prim_uf[st_prim_uf$Gene == g, ]
  if (nrow(p_row) == 0 || is.na(p_row$log2FC)) {
    next
  }
  
  prim_lfc <- p_row$log2FC
  prim_p <- p_row$PValue
  prim_fdr <- p_row$BH_FDR
  ev_status <- p_row$Evidence_Status
  
  if (ev_status == "INSUFFICIENT DATA") {
    robust_rows[[g]] <- data.frame(
      Gene = g,
      Evidence = "INSUFFICIENT DATA",
      Primary_log2FC = prim_lfc,
      Primary_P = prim_p,
      Primary_FDR = prim_fdr,
      Sign_b = "—", Sign_c = "—", Sign_d = "—", Sign_e = "—", Sign_f = "—",
      LODO_sign_consistent = "—",
      Verdict = "INSUFFICIENT DATA",
      stringsAsFactors = FALSE
    )
    next
  }
  
  # Sensitivity signs
  get_sign <- function(df, aname) {
    sub <- df[df$Analysis == aname & df$Gene == g & df$Contrast == "LV_UF_vs_LV_NOT_UF", ]
    if (nrow(sub) > 0 && !is.na(sub$log2FC)) {
      if (sub$log2FC > 0) "+" else if (sub$log2FC < 0) "-" else "0"
    } else "NA"
  }
  
  prim_sign <- if (prim_lfc > 0) "+" else if (prim_lfc < 0) "-" else "0"
  s_b <- get_sign(all_sens, "sens_b_no_LV_UF3")
  s_c <- get_sign(all_sens, "sens_c_no_scDblFinder")
  s_d <- get_sign(all_sens, "sens_d_ge700genes")
  s_e <- get_sign(all_sens, "sens_e_ge50cells")
  s_f <- get_sign(all_sens, "sens_f_top2000_hvg_tmm")
  
  # LODO consistency: sign unchanged in EVERY leave-one-donor-out run
  lodo_subs <- all_sens[grepl("^lodo_out_", all_sens$Analysis) & all_sens$Gene == g & all_sens$Contrast == "LV_UF_vs_LV_NOT_UF", ]
  lodo_signs <- ifelse(lodo_subs$log2FC > 0, "+", ifelse(lodo_subs$log2FC < 0, "-", "0"))
  lodo_all_match <- (length(lodo_signs) == 16) && all(lodo_signs == prim_sign)
  n_lodo_match <- sum(lodo_signs == prim_sign)
  
  # Robustness rule:
  # (1) primary FDR < 0.05
  # (2) sign unchanged in (b), (c), (d), (e)
  # (3) sign unchanged in every LODO run (16/16)
  c1 <- (!is.na(prim_fdr) && prim_fdr < 0.05)
  c2 <- (s_b == prim_sign && s_c == prim_sign && s_d == prim_sign && s_e == prim_sign)
  c3 <- lodo_all_match
  
  verdict <- if (c1 && c2 && c3) "supported" else "not robust"
  
  robust_rows[[g]] <- data.frame(
    Gene = g,
    Evidence = "PASS",
    Primary_log2FC = prim_lfc,
    Primary_P = signif(prim_p, 4),
    Primary_FDR = signif(prim_fdr, 4),
    Sign_Primary = prim_sign,
    Sign_b_noUF3 = s_b,
    Sign_c_noDbl = s_c,
    Sign_d_ge700g = s_d,
    Sign_e_ge50cells = s_e,
    Sign_f_HVG = s_f,
    LODO_Consistent = sprintf("%d/16", n_lodo_match),
    Verdict = verdict,
    stringsAsFactors = FALSE
  )
}

df_robust <- do.call(rbind, robust_rows)
write.csv(df_robust, "results/tables/stage4_stromal_robustness_verdicts.csv", row.names = FALSE)
print(df_robust, row.names = FALSE)
cat("\nStep 3 & Step 4 completed successfully.\n")
