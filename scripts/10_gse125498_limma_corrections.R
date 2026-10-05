# ==============================================================================
# SCRIPT 10: GSE125498 LIMMA DIFFERENTIAL EXPRESSION & DETECTION RECOMPUTATION
# Fixed Seed: 42
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(dplyr)
})

cat("=================================================================\n")
cat("1. LIMMA DIFFERENTIAL EXPRESSION IN GSE125498 (LPD vs SPD)\n")
cat("=================================================================\n")

matrix_file <- "data/raw/GSE125498_series_matrix.txt.gz"
annot_file  <- "data/raw/GPL10558.annot.gz"
raw_supp_file <- "data/raw/GSE125498_non-normalized_data.txt.gz"

# 1. Load Series Matrix
gse <- GEOquery::getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
expr_mat <- exprs(gse)
pdata <- pData(gse)

cat(sprintf("Loaded GSE125498 Expression Matrix: %d probes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# Classify samples into Short-term PD (SPD) vs Long-term PD (LPD)
group <- ifelse(grepl("SPD|short", paste(pdata$title, pdata$characteristics_ch1), ignore.case = TRUE), "SPD", "LPD")
group <- factor(group, levels = c("SPD", "LPD"))
cat("Sample classification breakdown:\n")
print(table(group))

# Setup Limma design (cell means model: SPD as reference)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_mat, design)

contrast_mat <- makeContrasts(LPD_vs_SPD = LPD - SPD, levels = design)
fit_contr <- contrasts.fit(fit, contrast_mat)
fit_eb <- eBayes(fit_contr)

# Full array results
tt_all <- topTable(fit_eb, coef = "LPD_vs_SPD", number = Inf, sort.by = "none")
tt_all$Probe_ID <- rownames(tt_all)
tt_all$Adj_P_BH_Array <- p.adjust(tt_all$P.Value, method = "BH")

# Primary MaxMean mapped probes for the 11 hub genes
hub_probes <- data.frame(
  Gene = c("ISM1", "FN1", "VCAN", "COL3A1", "COL8A1", "THBS3", "LOX"),
  Probe_ID = c("ILMN_3239288", "ILMN_2366463", "ILMN_1687301", "ILMN_1773079", "ILMN_1685433", "ILMN_1804663", "ILMN_1695880"),
  stringsAsFactors = FALSE
)

# Extract moderated statistics for mapped hubs
hub_limma_rows <- list()
for (i in 1:nrow(hub_probes)) {
  g <- hub_probes$Gene[i]
  p <- hub_probes$Probe_ID[i]
  
  lfc  <- fit_eb$coefficients[p, "LPD_vs_SPD"]
  tval <- fit_eb$t[p, "LPD_vs_SPD"]
  se   <- fit_eb$stdev.unscaled[p, "LPD_vs_SPD"] * sqrt(fit_eb$s2.post[p])
  pval <- fit_eb$p.value[p, "LPD_vs_SPD"]
  df_t <- fit_eb$df.total[p]
  crit_t <- qt(0.975, df = df_t)
  ci_low <- lfc - crit_t * se
  ci_high <- lfc + crit_t * se
  adjp_array <- tt_all$Adj_P_BH_Array[tt_all$Probe_ID == p]
  
  hub_limma_rows[[i]] <- data.frame(
    Gene = g,
    Probe_ID = p,
    log2FC = lfc,
    SE_Limma = se,
    CI_95_Lower_Limma = ci_low,
    CI_95_Upper_Limma = ci_high,
    Moderated_t = tval,
    Moderated_df = df_t,
    Raw_P_Limma = pval,
    Adj_P_BH_Array = adjp_array,
    N_SPD = sum(group == "SPD"),
    N_LPD = sum(group == "LPD"),
    Compartment = "Peritoneal Effluent Cells (PD-Duration Proxy, NOT Tissue, NOT EPS)",
    stringsAsFactors = FALSE
  )
}

df_hub_limma <- do.call(rbind, hub_limma_rows)
# FDR across mapped hub genes
df_hub_limma$Adj_P_BH_Hubs <- p.adjust(df_hub_limma$Raw_P_Limma, method = "BH")

write.csv(df_hub_limma, "results/tables/C_gse125498_hub_limma.csv", row.names = FALSE)
cat("Saved C_gse125498_hub_limma.csv with moderated limma statistics.\n")
print(df_hub_limma[, c("Gene", "Probe_ID", "log2FC", "SE_Limma", "Raw_P_Limma", "Adj_P_BH_Hubs", "Adj_P_BH_Array")])

cat("\n=================================================================\n")
cat("2. REGENERATE C_gse125498_concordance.csv WITH LIMMA MODERATED RESULTS\n")
cat("=================================================================\n")

# Load Discovery sensitivity results
sens_file <- "results/tables/B_gse62928_sensitivity_contrasts_uncertainty.csv"
df_sens <- read.csv(sens_file, stringsAsFactors = FALSE)

all_11_hubs <- c("ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", 
                 "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX")

conc_rows <- list()
for (g in all_11_hubs) {
  disc_sub <- df_sens %>% filter(Gene == g, Contrast == "EPS_vs_PooledControl")
  disc_lfc <- disc_sub$log2FC_MaxMean[1]
  disc_p   <- disc_sub$Raw_P_MaxMean[1]
  
  if (g %in% df_hub_limma$Gene) {
    val_sub <- df_hub_limma %>% filter(Gene == g)
    val_lfc <- val_sub$log2FC[1]
    val_se  <- val_sub$SE_Limma[1]
    val_p   <- val_sub$Raw_P_Limma[1]
    ci_str  <- sprintf("[%+.3f, %+.3f]", val_sub$CI_95_Lower_Limma[1], val_sub$CI_95_Upper_Limma[1])
    
    # Concordance rule:
    # Concordant (significant, P < 0.05): same sign & P < 0.05
    # Same direction (not significant, P >= 0.05): same sign & P >= 0.05
    # Discordant: opposite sign
    # Flat (|log2FC| < 0.1)
    if (abs(val_lfc) < 0.10) {
      cat_label <- "Flat (|log2FC| < 0.1)"
    } else if ((val_lfc > 0 && disc_lfc > 0) || (val_lfc < 0 && disc_lfc < 0)) {
      if (val_p < 0.05) {
        cat_label <- "Concordant (significant, P < 0.05)"
      } else {
        cat_label <- "Same direction (not significant, P >= 0.05)"
      }
    } else {
      cat_label <- "Discordant"
    }
    
    conc_rows[[length(conc_rows) + 1]] <- data.frame(
      Gene = g,
      Status = "MAPPED",
      Discovery_Tissue_EPS_vs_Pooled_log2FC = disc_lfc,
      Discovery_Tissue_Raw_P = disc_p,
      Validation_Effluent_log2FC = val_lfc,
      Validation_Effluent_SE_Limma = val_se,
      Validation_Effluent_95_CI = ci_str,
      Validation_Effluent_Raw_P_Limma = val_p,
      Concordance_Category = cat_label,
      Compartment_Note = "Effluent cells (PD duration proxy) vs Tissue (EPS)",
      stringsAsFactors = FALSE
    )
  } else {
    conc_rows[[length(conc_rows) + 1]] <- data.frame(
      Gene = g,
      Status = "NOT MAPPED ON GPL10558",
      Discovery_Tissue_EPS_vs_Pooled_log2FC = disc_lfc,
      Discovery_Tissue_Raw_P = disc_p,
      Validation_Effluent_log2FC = NA,
      Validation_Effluent_SE_Limma = NA,
      Validation_Effluent_95_CI = "N/A",
      Validation_Effluent_Raw_P_Limma = NA,
      Concordance_Category = "Not mapped (No probe on array)",
      Compartment_Note = "Unmapped on Illumina HT-12 v4 beadchip",
      stringsAsFactors = FALSE
    )
  }
}

df_conc <- do.call(rbind, conc_rows)
write.csv(df_conc, "results/tables/C_gse125498_concordance.csv", row.names = FALSE)
cat("Saved C_gse125498_concordance.csv.\n")
print(df_conc[, c("Gene", "Status", "Validation_Effluent_log2FC", "Validation_Effluent_Raw_P_Limma", "Concordance_Category")])

cat("\nScript 10 completed successfully!\n")
