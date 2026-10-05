# ==============================================================================
# SCRIPT 04: HOUSEKEEPING FIXES - ALL CONTRASTS FOR ALL 7 POSITIVE CONTROLS
# Fixed Seed: 42
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(dplyr)
})

cat("=================================================================\n")
cat("HOUSEKEEPING H1: COMPLETE POSITIVE CONTROL CONTRASTS (GSE62928)\n")
cat("=================================================================\n")

matrix_file <- "results/tables/GSE62928_full_expression_matrix.csv"
meta_file   <- "results/tables/GSE62928_sample_metadata.csv"

expr_mat <- read.csv(matrix_file, row.names = 1, check.names = FALSE)
meta_df  <- read.csv(meta_file, stringsAsFactors = FALSE)

meta_df$group_3class <- factor(meta_df$group_3class, levels = c("Uremic", "PD", "EPS"))
design_3grp <- model.matrix(~ 0 + group_3class, data = meta_df)
colnames(design_3grp) <- levels(meta_df$group_3class)

contrast_matrix <- makeContrasts(
  EPS_vs_PD = EPS - PD,
  EPS_vs_Uremic = EPS - Uremic,
  PD_vs_Uremic = PD - Uremic,
  EPS_vs_PooledControl = EPS - (PD + Uremic)/2,
  levels = design_3grp
)

fit <- lmFit(expr_mat, design_3grp)
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrasts)

pos_controls <- c("FN1", "THBS1", "COL1A1", "ACTA2", "SULF1", "LEP", "RBP4")
pos_ctrl_rows <- list()

for (gene in pos_controls) {
  if (!gene %in% rownames(expr_mat)) {
    cat(sprintf("WARNING: %s not in matrix\n", gene))
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
      if (c_name %in% c("EPS_vs_PD", "EPS_vs_Uremic", "EPS_vs_PooledControl")) {
        "UP (+)"
      } else if (c_name == "PD_vs_Uremic") {
        "UP (+) [Peritoneal Dialysis effect]"
      } else {
        "DESCRIPTIVE"
      }
    } else {
      # LEP, RBP4
      if (c_name %in% c("EPS_vs_Uremic", "PD_vs_Uremic", "EPS_vs_PooledControl")) {
        "DOWN (-)"
      } else if (c_name == "EPS_vs_PD") {
        "DESCRIPTIVE"
      } else {
        "DESCRIPTIVE"
      }
    }
    
    observed_sign <- if (lfc > 0) "UP (+)" else "DOWN (-)"
    
    sign_check <- if (grepl("DESCRIPTIVE", expected_sign)) {
      "DESCRIPTIVE"
    } else if (grepl("UP", expected_sign) && observed_sign == "UP (+)") {
      "PASSED"
    } else if (grepl("DOWN", expected_sign) && observed_sign == "DOWN (-)") {
      "PASSED"
    } else {
      "FAILED_REVERSED"
    }
    
    perm_note <- if (c_name %in% c("EPS_vs_PD", "EPS_vs_Uremic")) {
      "4 vs 2 (15 splits: min 1-sided P = 1/15 = 0.067, min 2-sided P = 2/15 = 0.133)"
    } else if (c_name == "PD_vs_Uremic") {
      "2 vs 2 (6 splits: min 1-sided P = 1/6 = 0.167, min 2-sided P = 2/6 = 0.333)"
    } else {
      "4 vs 4 (70 splits: min 1-sided P = 1/70 = 0.0143, min 2-sided P = 2/70 = 0.0286)"
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
      Permutation_Note = perm_note,
      stringsAsFactors = FALSE
    )
  }
}

pos_ctrl_df <- do.call(rbind, pos_ctrl_rows)
write.csv(pos_ctrl_df, "results/tables/B_gse62928_positive_controls.csv", row.names = FALSE)
cat(sprintf("Saved updated complete positive controls table with all 28 rows (7 genes x 4 contrasts) to results/tables/B_gse62928_positive_controls.csv\n"))
