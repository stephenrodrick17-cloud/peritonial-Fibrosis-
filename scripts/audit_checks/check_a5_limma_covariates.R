# scripts/audit_checks/check_a5_limma_covariates.R
suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(GEOquery)
})

matrix_file <- "data/raw/GSE125498_series_matrix.txt.gz"
gse <- getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
expr_mat <- exprs(gse)
pdata <- pData(gse)

# Classify samples into Short-term PD (SPD) vs Long-term PD (LPD)
group <- ifelse(grepl("SPD|short", paste(pdata$title, pdata$characteristics_ch1), ignore.case = TRUE), "SPD", "LPD")
group <- factor(group, levels = c("SPD", "LPD"))
cat("Sample classification breakdown:\n")
print(table(group))

# Probes of interest:
# VCAN: ILMN_1687301
# FN1: ILMN_2366463
# CD3E: ILMN_1739794
# CD14: ILMN_1740015, ILMN_2396444
# FCGR3B: ILMN_1728639, ILMN_2134453

cat("\n--- Covariate Probe Details ---\n")
cov_probes <- c("ILMN_1739794", "ILMN_1740015", "ILMN_2396444", "ILMN_1728639", "ILMN_2134453")
for (p in cov_probes) {
  cat(sprintf("%s: mean = %.3f, sd = %.3f\n", p, mean(expr_mat[p, ]), sd(expr_mat[p, ])))
}

# 1. Unadjusted Model (Primary Limma as in Script 10)
cat("\n=== 1. UNADJUSTED MODEL: ~ group ===\n")
design_unadj <- model.matrix(~ group)
fit_unadj <- lmFit(expr_mat, design_unadj)
fit_unadj <- eBayes(fit_unadj)

for (gene in c("VCAN", "FN1")) {
  p_id <- ifelse(gene == "VCAN", "ILMN_1687301", "ILMN_2366463")
  lfc <- fit_unadj$coefficients[p_id, "groupLPD"]
  pval <- fit_unadj$p.value[p_id, "groupLPD"]
  tval <- fit_unadj$t[p_id, "groupLPD"]
  cat(sprintf("%s (%s): log2FC = %.4f, t = %.4f, P = %.5e\n", gene, p_id, lfc, tval, pval))
}

# 2. Adjusted Model using primary probes for CD14, CD3E, FCGR3B
# CD14: ILMN_1740015
# CD3E: ILMN_1739794
# FCGR3B: ILMN_1728639
cd14_expr   <- expr_mat["ILMN_1740015", ]
cd3e_expr   <- expr_mat["ILMN_1739794", ]
fcgr3b_expr <- expr_mat["ILMN_1728639", ]

cat("\n=== 2. ADJUSTED MODEL: ~ group + CD14 (ILMN_1740015) + CD3E (ILMN_1739794) + FCGR3B (ILMN_1728639) ===\n")
design_adj <- model.matrix(~ group + cd14_expr + cd3e_expr + fcgr3b_expr)
fit_adj <- lmFit(expr_mat, design_adj)
fit_adj <- eBayes(fit_adj)

for (gene in c("VCAN", "FN1")) {
  p_id <- ifelse(gene == "VCAN", "ILMN_1687301", "ILMN_2366463")
  lfc <- fit_adj$coefficients[p_id, "groupLPD"]
  pval <- fit_adj$p.value[p_id, "groupLPD"]
  tval <- fit_adj$t[p_id, "groupLPD"]
  cat(sprintf("%s (%s): log2FC = %.4f, t = %.4f, P = %.5e\n", gene, p_id, lfc, tval, pval))
}

# 3. Adjusted Model with alternative FCGR3B probe ILMN_2134453
fcgr3b_expr2 <- expr_mat["ILMN_2134453", ]
cat("\n=== 3. ADJUSTED MODEL with FCGR3B (ILMN_2134453) ===\n")
design_adj2 <- model.matrix(~ group + cd14_expr + cd3e_expr + fcgr3b_expr2)
fit_adj2 <- lmFit(expr_mat, design_adj2)
fit_adj2 <- eBayes(fit_adj2)

for (gene in c("VCAN", "FN1")) {
  p_id <- ifelse(gene == "VCAN", "ILMN_1687301", "ILMN_2366463")
  lfc <- fit_adj2$coefficients[p_id, "groupLPD"]
  pval <- fit_adj2$p.value[p_id, "groupLPD"]
  tval <- fit_adj2$t[p_id, "groupLPD"]
  cat(sprintf("%s (%s): log2FC = %.4f, t = %.4f, P = %.5e\n", gene, p_id, lfc, tval, pval))
}

# 4. Single-covariate adjustments
cat("\n=== 4. SINGLE-COVARIATE ADJUSTMENTS ===\n")
for (cov_name in c("CD14", "CD3E", "FCGR3B")) {
  if (cov_name == "CD14") c_val <- cd14_expr
  if (cov_name == "CD3E") c_val <- cd3e_expr
  if (cov_name == "FCGR3B") c_val <- fcgr3b_expr
  d_single <- model.matrix(~ group + c_val)
  f_single <- eBayes(lmFit(expr_mat, d_single))
  cat(sprintf("Adjusting for %s only:\n", cov_name))
  for (gene in c("VCAN", "FN1")) {
    p_id <- ifelse(gene == "VCAN", "ILMN_1687301", "ILMN_2366463")
    lfc <- f_single$coefficients[p_id, "groupLPD"]
    pval <- f_single$p.value[p_id, "groupLPD"]
    cat(sprintf("  %s: log2FC = %.4f, P = %.5e\n", gene, lfc, pval))
  }
}
