# scripts/audit_checks/check_directive5_limma.R
suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(GEOquery)
})

matrix_file <- "data/raw/GSE125498_series_matrix.txt.gz"
gse <- getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
expr_mat <- exprs(gse)
pdata <- pData(gse)

group <- ifelse(grepl("SPD|short", paste(pdata$title, pdata$characteristics_ch1), ignore.case = TRUE), "SPD", "LPD")
group <- factor(group, levels = c("SPD", "LPD"))

# Top-table verified probes:
# CD3E:   ILMN_1739794 (logFC = +1.4609, P = 0.00305)
# CD14:   ILMN_2396444 (logFC = -0.2435, P = 0.0445)
# FCGR3B: ILMN_2134453 (logFC = +0.4872, P = 0.0281)
# Alternative probes:
# CD14_alt:   ILMN_1740015 (logFC = -0.3791, P = 0.0928)
# FCGR3B_alt: ILMN_1728639 (logFC = +0.9670, P = 0.0532)
# FN1 probes:
# FN1_top:     ILMN_1778237 (logFC = +0.4068, P = 0.2690)
# FN1_maxmean: ILMN_2366463 (logFC = +0.3171, P = 0.3188)
# VCAN:        ILMN_1687301 (logFC = -0.5224, P = 0.0244)

cat("=== TOP-TABLE VERIFICATION (UNADJUSTED) ===\n")
fit0 <- eBayes(lmFit(expr_mat, model.matrix(~ group)))
test_probes <- c("ILMN_1739794", "ILMN_2396444", "ILMN_1740015", "ILMN_2134453", "ILMN_1728639", "ILMN_1778237", "ILMN_2366463", "ILMN_1687301")
probe_names <- c("CD3E", "CD14 (ILMN_2396444)", "CD14 (ILMN_1740015)", "FCGR3B (ILMN_2134453)", "FCGR3B (ILMN_1728639)", "FN1 (ILMN_1778237, top-table)", "FN1 (ILMN_2366463, MaxMean)", "VCAN (ILMN_1687301)")

for (i in seq_along(test_probes)) {
  p <- test_probes[i]
  cat(sprintf("%-32s (%s): log2FC = %+.4f, t = %+.4f, P = %.5f\n", 
              probe_names[i], p, fit0$coefficients[p, "groupLPD"], fit0$t[p, "groupLPD"], fit0$p.value[p, "groupLPD"]))
}

# Run models across both probe choices:
# Set A: Top-table probes (CD14=ILMN_2396444, FCGR3B=ILMN_2134453)
# Set B: Alternate probes (CD14=ILMN_1740015, FCGR3B=ILMN_1728639)

run_covariate_suite <- run_models <- function(cd14_p, fcgr3b_p, label) {
  cat(sprintf("\n======================================================================\n"))
  cat(sprintf("COVARIATE SUITE: %s\n", label))
  cat(sprintf("CD3E: ILMN_1739794 | CD14: %s | FCGR3B: %s\n", cd14_p, fcgr3b_p))
  cat(sprintf("======================================================================\n"))
  
  c_cd3e   <- expr_mat["ILMN_1739794", ]
  c_cd14   <- expr_mat[cd14_p, ]
  c_fcgr3b <- expr_mat[fcgr3b_p, ]
  
  # Single covariates
  models <- list(
    "Unadjusted"          = model.matrix(~ group),
    "Adjust: CD14 only"   = model.matrix(~ group + c_cd14),
    "Adjust: CD3E only"   = model.matrix(~ group + c_cd3e),
    "Adjust: FCGR3B only" = model.matrix(~ group + c_fcgr3b),
    "Drop CD14 (CD3E+FCGR3B)" = model.matrix(~ group + c_cd3e + c_fcgr3b),
    "Drop FCGR3B (CD14+CD3E)" = model.matrix(~ group + c_cd14 + c_cd3e),
    "Drop CD3E (CD14+FCGR3B)" = model.matrix(~ group + c_cd14 + c_fcgr3b),
    "Joint: All 3 covs"   = model.matrix(~ group + c_cd14 + c_cd3e + c_fcgr3b)
  )
  
  for (mname in names(models)) {
    fit_m <- eBayes(lmFit(expr_mat, models[[mname]]))
    cat(sprintf("\n--- %s ---\n", mname))
    for (target in c("VCAN (ILMN_1687301)", "FN1 (ILMN_1778237, top)", "FN1 (ILMN_2366463, maxmean)")) {
      p_id <- if (grepl("VCAN", target)) "ILMN_1687301" else if (grepl("1778237", target)) "ILMN_1778237" else "ILMN_2366463"
      lfc <- fit_m$coefficients[p_id, "groupLPD"]
      pval <- fit_m$p.value[p_id, "groupLPD"]
      tval <- fit_m$t[p_id, "groupLPD"]
      cat(sprintf("  %-30s: log2FC = %+.4f, t = %+.4f, P = %.5e\n", target, lfc, tval, pval))
    }
  }
}

run_covariate_suite("ILMN_2396444", "ILMN_2134453", "Set A: Top-Table Probe Set")
run_covariate_suite("ILMN_1740015", "ILMN_1728639", "Set B: Probe Set B")
