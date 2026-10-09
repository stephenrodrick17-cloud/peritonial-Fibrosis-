# ==============================================================================
# SCRIPT 02: POSITIVE CONTROLS & SENSITIVITY CONTRASTS WITH UNCERTAINTY (GSE62928)
# Fixed Seed: 42
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(dplyr)
})

cat("=================================================================\n")
cat("1. FACTOR LEVEL ORDER & CONTRAST MATRIX VERIFICATION\n")
cat("=================================================================\n")

matrix_file <- "results/tables/GSE62928_full_expression_matrix.csv"
meta_file   <- "results/tables/GSE62928_sample_metadata.csv"

if (!file.exists(matrix_file) || !file.exists(meta_file)) {
  stop("GSE62928 full expression matrix or metadata not found!")
}

expr_mat <- read.csv(matrix_file, row.names = 1, check.names = FALSE)
meta_df  <- read.csv(meta_file, stringsAsFactors = FALSE)

# Explicit Factor Ordering: Uremic (control) -> PD -> EPS
meta_df$group_3class <- factor(meta_df$group_3class, levels = c("Uremic", "PD", "EPS"))
cat("Sample IDs, Titles, and Clinical Classifications:\n")
print(meta_df[, c("sample_id", "sample_title", "group_3class")])

cat("\nFactor Levels explicitly defined as: ", paste(levels(meta_df$group_3class), collapse = " -> "), "\n")

# Design Matrix without intercept (cell means model)
design_3grp <- model.matrix(~ 0 + group_3class, data = meta_df)
colnames(design_3grp) <- levels(meta_df$group_3class)
cat("\nDesign Matrix:\n")
print(design_3grp)

# Contrast Matrix:
# EPS_vs_PD: EPS - PD (Expect positive for pro-fibrotic genes)
# EPS_vs_Uremic: EPS - Uremic (Expect positive for pro-fibrotic genes)
# PD_vs_Uremic: PD - Uremic
# EPS_vs_PooledControl: EPS - (PD + Uremic)/2
contrast_matrix <- makeContrasts(
  EPS_vs_PD = EPS - PD,
  EPS_vs_Uremic = EPS - Uremic,
  PD_vs_Uremic = PD - Uremic,
  EPS_vs_PooledControl = EPS - (PD + Uremic)/2,
  levels = design_3grp
)

cat("\nContrast Matrix:\n")
print(contrast_matrix)

# Fit Limma model on full MaxMean expression matrix
fit <- lmFit(expr_mat, design_3grp)
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrasts)

cat("\n=================================================================\n")
cat("2. POSITIVE-CONTROL SIGN CHECK (GSE62928 - Reimold et al. 2013)\n")
cat("=================================================================\n")

pos_controls <- c("FN1", "THBS1", "COL1A1", "ACTA2", "SULF1", "LEP", "RBP4")
pos_ctrl_rows <- list()

for (gene in pos_controls) {
  if (!gene %in% rownames(expr_mat)) {
    cat(sprintf("WARNING: Positive control %s not in matrix!\n", gene))
    next
  }
  for (c_name in colnames(contrast_matrix)) {
    lfc  <- fit_ebayes$coefficients[gene, c_name]
    tval <- fit_ebayes$t[gene, c_name]
    se   <- fit_ebayes$stdev.unscaled[gene, c_name] * sqrt(fit_ebayes$s2.post[gene])
    pval <- fit_ebayes$p.value[gene, c_name]
    df_t <- fit_ebayes$df.total[gene]
    crit_t <- qt(0.975, df = df_t)
    ci_low <- lfc - crit_t * se
    ci_high <- lfc + crit_t * se
    
    all_p <- fit_ebayes$p.value[, c_name]
    adjp  <- p.adjust(all_p, method = "BH")[gene]
    
    expected_sign <- if (gene %in% c("FN1", "THBS1", "COL1A1", "ACTA2", "SULF1")) {
      if (c_name %in% c("EPS_vs_PD", "EPS_vs_Uremic", "EPS_vs_PooledControl")) "UP (+)" else "N/A"
    } else {
      if (c_name %in% c("EPS_vs_Uremic", "PD_vs_Uremic", "EPS_vs_PooledControl")) "DOWN (-)" else "N/A"
    }
    
    observed_sign <- if (lfc > 0) "UP (+)" else "DOWN (-)"
    sign_check <- if (expected_sign == "N/A") {
      "DESCRIPTIVE"
    } else if (expected_sign == observed_sign) {
      "PASSED"
    } else {
      "FAILED_REVERSED"
    }
    
    pos_ctrl_rows[[length(pos_ctrl_rows) + 1]] <- data.frame(
      Gene = gene,
      Contrast = c_name,
      log2FC = lfc,
      SE = se,
      CI_95_Lower = ci_low,
      CI_95_Upper = ci_high,
      t_statistic = tval,
      Raw_P = pval,
      Adj_P_BH = adjp,
      Expected_Sign = expected_sign,
      Observed_Sign = observed_sign,
      Sign_Check = sign_check,
      stringsAsFactors = FALSE
    )
  }
}

pos_ctrl_df <- do.call(rbind, pos_ctrl_rows)
write.csv(pos_ctrl_df, "results/tables/B_gse62928_positive_controls.csv", row.names = FALSE)
cat("Saved positive controls table to results/tables/B_gse62928_positive_controls.csv\n\n")

print(pos_ctrl_df[, c("Gene", "Contrast", "log2FC", "SE", "Raw_P", "Adj_P_BH", "Expected_Sign", "Sign_Check")])

# Verify whether any critical sign check failed
failed_checks <- pos_ctrl_df %>% filter(Sign_Check == "FAILED_REVERSED")
if (nrow(failed_checks) > 0) {
  cat("\nCRITICAL FAILURE: Positive control sign check failed!\n")
  print(failed_checks)
  stop("HALTING: Positive control direction is reversed. Diagnose contrast coding before proceeding.")
} else {
  cat("\nALL POSITIVE CONTROL SIGN CHECKS PASSED SUCCESSFULLY!\n")
}

cat("\n=================================================================\n")
cat("3. SENSITIVITY TABLE WITH UNCERTAINTY & SIDE-BY-SIDE COLLAPSE METHODS\n")
cat("=================================================================\n")

# Hub genes list
hub_genes <- c("ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", 
               "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX")

# A. MaxMean sensitivity across all 4 contrasts
hub_sens_rows <- list()
for (gene in hub_genes) {
  if (!gene %in% rownames(expr_mat)) next
  for (c_name in colnames(contrast_matrix)) {
    lfc  <- fit_ebayes$coefficients[gene, c_name]
    tval <- fit_ebayes$t[gene, c_name]
    se   <- fit_ebayes$stdev.unscaled[gene, c_name] * sqrt(fit_ebayes$s2.post[gene])
    pval <- fit_ebayes$p.value[gene, c_name]
    df_t <- fit_ebayes$df.total[gene]
    crit_t <- qt(0.975, df = df_t)
    ci_low <- lfc - crit_t * se
    ci_high <- lfc + crit_t * se
    all_p <- fit_ebayes$p.value[, c_name]
    adjp  <- p.adjust(all_p, method = "BH")[gene]
    
    hub_sens_rows[[length(hub_sens_rows) + 1]] <- data.frame(
      Gene = gene,
      Contrast = c_name,
      log2FC_MaxMean = lfc,
      SE_MaxMean = se,
      CI_95_Lower_MaxMean = ci_low,
      CI_95_Upper_MaxMean = ci_high,
      t_MaxMean = tval,
      Raw_P_MaxMean = pval,
      Adj_P_BH_MaxMean = adjp,
      stringsAsFactors = FALSE
    )
  }
}
hub_sens_df <- do.call(rbind, hub_sens_rows)

# B. Best-P Probe Collapse from toptable
toptable_orig_file <- "GSE62928.top.table.tsv"
top_raw <- read.table(toptable_orig_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
top_annot <- top_raw[!is.na(top_raw$Gene.symbol) & top_raw$Gene.symbol != "" & top_raw$Gene.symbol != "---", ]
top_annot$Gene.symbol <- trimws(top_annot$Gene.symbol)
split_syms <- strsplit(top_annot$Gene.symbol, " ?/// ?")
n_syms <- lengths(split_syms)
top_exp <- top_annot[rep(seq_len(nrow(top_annot)), n_syms), ]
top_exp$Gene.symbol <- trimws(unlist(split_syms))
top_best_p <- top_exp[order(top_exp$P.Value), ]
top_best_p <- top_best_p[!duplicated(top_best_p$Gene.symbol), ]
rownames(top_best_p) <- top_best_p$Gene.symbol

# Side-by-side Table Construction
side_by_side_rows <- list()
for (gene in hub_genes) {
  for (c_name in colnames(contrast_matrix)) {
    mm_sub <- hub_sens_df %>% filter(Gene == gene, Contrast == c_name)
    
    # Best-P values from discovery toptable for EPS_vs_PooledControl
    if (c_name == "EPS_vs_PooledControl" && gene %in% rownames(top_best_p)) {
      lfc_bp  <- top_best_p[gene, "logFC"]
      t_bp    <- top_best_p[gene, "t"]
      pval_bp <- top_best_p[gene, "P.Value"]
      adjp_bp <- top_best_p[gene, "adj.P.Val"]
      se_bp   <- abs(lfc_bp / t_bp)
      ci_low_bp  <- lfc_bp - 1.96 * se_bp
      ci_high_bp <- lfc_bp + 1.96 * se_bp
      ci_str_bp  <- sprintf("[%+.3f, %+.3f]", ci_low_bp, ci_high_bp)
    } else {
      lfc_bp  <- NA
      se_bp   <- NA
      ci_str_bp <- "N/A"
      pval_bp <- NA
      adjp_bp <- NA
    }
    
    perm_note <- if (c_name %in% c("EPS_vs_PD", "EPS_vs_Uremic")) {
      "Descriptive (4 vs 2: 15 splits, min 1-sided P = 1/15 = 0.067, min 2-sided P = 2/15 = 0.133)"
    } else if (c_name == "PD_vs_Uremic") {
      "Descriptive (2 vs 2: 6 splits, min 1-sided P = 1/6 = 0.167, min 2-sided P = 2/6 = 0.333)"
    } else {
      "Discovery (4 vs 4: 70 splits, min 1-sided P = 1/70 = 0.0143, min 2-sided P = 2/70 = 0.0286)"
    }
    
    side_by_side_rows[[length(side_by_side_rows) + 1]] <- data.frame(
      Gene = gene,
      Contrast = c_name,
      log2FC_BestProbe = lfc_bp,
      SE_BestProbe = se_bp,
      CI_95_BestProbe = ci_str_bp,
      Raw_P_BestProbe = pval_bp,
      Adj_P_BestProbe = adjp_bp,
      log2FC_MaxMean = mm_sub$log2FC_MaxMean,
      SE_MaxMean = mm_sub$SE_MaxMean,
      CI_95_MaxMean = sprintf("[%+.3f, %+.3f]", mm_sub$CI_95_Lower_MaxMean, mm_sub$CI_95_Upper_MaxMean),
      Raw_P_MaxMean = mm_sub$Raw_P_MaxMean,
      Adj_P_BH_MaxMean = mm_sub$Adj_P_BH_MaxMean,
      Statistical_Limitation_Note = perm_note,
      stringsAsFactors = FALSE
    )
  }
}

side_by_side_df <- do.call(rbind, side_by_side_rows)
write.csv(side_by_side_df, "results/tables/B_gse62928_sensitivity_contrasts_uncertainty.csv", row.names = FALSE)
cat("Saved detailed sensitivity contrasts with uncertainty to results/tables/B_gse62928_sensitivity_contrasts_uncertainty.csv\n")

cat("\nScript 02 completed successfully!\n")
