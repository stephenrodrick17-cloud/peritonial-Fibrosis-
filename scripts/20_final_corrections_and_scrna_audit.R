# Script 20: Limma Degrees of Freedom, exact CI, and Table Regeneration for GSE125498
suppressPackageStartupMessages({
  library(limma)
})

cat("=================================================================\n")
cat("SCRIPT 20: GSE125498 LIMMA DEGREES OF FREEDOM & CI REGENERATION\n")
cat("=================================================================\n")

# Load raw expression matrix and metadata
expr_file <- "results/tables/GSE125498_expression_matrix_hubs_maxmean_primary.csv"
meta_file <- "results/tables/GSE125498_sample_metadata.csv"

expr_df <- read.csv(expr_file, row.names = 1, check.names = FALSE)
meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)

# Ensure matching sample order
common_samples <- intersect(colnames(expr_df), meta_df$sample_id)
expr_mat <- as.matrix(expr_df[, common_samples])
meta_df <- meta_df[match(common_samples, meta_df$sample_id), ]

group <- factor(meta_df$stage_binary, levels = c("0", "1"), labels = c("SPD", "LPD"))
design <- model.matrix(~ 0 + group)
colnames(design) <- c("SPD", "LPD")

# Fit linear model
fit <- lmFit(expr_mat, design)
contrast_mat <- makeContrasts(LPD_vs_SPD = LPD - SPD, levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit_eb <- eBayes(fit2, trend = FALSE, robust = FALSE)

# Extract Degrees of Freedom
df_residual <- fit_eb$df.residual[1]
df_prior <- fit_eb$df.prior
df_total <- fit_eb$df.total[1]

cat(sprintf("fit$df.residual: %d\n", df_residual))
cat(sprintf("fit$df.prior:    %.4f\n", df_prior))
cat(sprintf("fit$df.total:    %.4f (per gene: %.4f)\n", df_total, fit_eb$df.total[1]))

# Compute moderated SE and exact t-based 95% CI per gene
# In limma, standard error of the contrast is stdev.unscaled * sqrt(s2.post)
se_moderated <- fit_eb$stdev.unscaled[, 1] * sqrt(fit_eb$s2.post)
coefs <- fit_eb$coefficients[, 1]
t_stats <- fit_eb$t[, 1]
p_raw <- fit_eb$p.value[, 1]

# qt(0.975, df.total)
t_crit <- qt(0.975, df = fit_eb$df.total)
ci_lower <- coefs - t_crit * se_moderated
ci_upper <- coefs + t_crit * se_moderated

# Probes and genes
gene_names <- rownames(expr_mat)
probe_mapping <- data.frame(
  Gene = c("ISM1", "FN1", "VCAN", "COL3A1", "COL8A1", "THBS3", "LOX"),
  Probe_ID = c("ILMN_3239288", "ILMN_2366463", "ILMN_1687301", "ILMN_1773079", "ILMN_1685433", "ILMN_1804663", "ILMN_1695880"),
  stringsAsFactors = FALSE
)
rownames(probe_mapping) <- probe_mapping$Gene

old_limma <- read.csv("results/tables/C_gse125498_hub_limma.csv", stringsAsFactors = FALSE)
rownames(old_limma) <- old_limma$Gene

hub_fdr <- p.adjust(p_raw, method = "BH")

# Legacy Welch P values from audit
legacy_welch <- c(
  ISM1 = 0.903080, FN1 = 0.268752, VCAN = 0.031535,
  COL3A1 = 0.636603, COL8A1 = 0.392095, THBS3 = 0.424177, LOX = 0.975440
)

# Full array BH-FDR from topTable
array_fdr <- c(
  ISM1 = 0.993875, FN1 = 0.861748, VCAN = 0.305679,
  COL3A1 = 0.949791, COL8A1 = 0.899602, THBS3 = 0.909989, LOX = 0.998495
)

res_df <- data.frame(
  Gene = gene_names,
  Probe_ID = probe_mapping[gene_names, "Probe_ID"],
  Coefficient_log2FC = round(coefs, 6),
  SE_Limma = round(se_moderated, 6),
  Moderated_df = round(fit_eb$df.total, 4),
  CI_95_Lower = round(ci_lower, 6),
  CI_95_Upper = round(ci_upper, 6),
  Moderated_t = round(t_stats, 6),
  Raw_P_Limma = round(p_raw, 6),
  BH_FDR_Across_11_Hubs = round(hub_fdr, 6),
  BH_FDR_Across_Array = round(array_fdr[gene_names], 6),
  N_SPD = sum(group == "SPD"),
  N_LPD = sum(group == "LPD"),
  Tissue_Context = "Peritoneal Effluent Cells (PD-Duration Proxy, NOT Tissue, NOT EPS)",
  Raw_P_Welch_Legacy = round(legacy_welch[gene_names], 6),
  stringsAsFactors = FALSE
)

# Write regenerated limma table
write.csv(res_df, "results/tables/C_gse125498_hub_limma.csv", row.names = FALSE)
cat("\nSaved regenerated results/tables/C_gse125498_hub_limma.csv\n")
print(res_df)

# Regenerate concordance table
discovery_tissue <- data.frame(
  Gene = c("ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"),
  Discovery_Tissue_EPS_vs_Pooled_log2FC = c(1.725877, 4.066665, 2.176920, 2.470463, 2.277435, 3.840111, 2.446207, 1.057567, 3.844884, 2.739727, 2.391670),
  Discovery_Tissue_Raw_P = c(0.000011, 0.000018, 0.003524, 0.000426, 0.000029, 0.000114, 0.000593, 0.002315, 0.000282, 0.000393, 0.002473),
  stringsAsFactors = FALSE
)

conc_records <- list()
for (i in 1:nrow(discovery_tissue)) {
  g <- discovery_tissue$Gene[i]
  d_fc <- discovery_tissue$Discovery_Tissue_EPS_vs_Pooled_log2FC[i]
  d_p <- discovery_tissue$Discovery_Tissue_Raw_P[i]
  
  if (g %in% res_df$Gene) {
    sub_r <- res_df[res_df$Gene == g, ]
    v_fc <- sub_r$Coefficient_log2FC
    v_se <- sub_r$SE_Limma
    v_ci <- sprintf("[%+.3f, %+.3f]", sub_r$CI_95_Lower, sub_r$CI_95_Upper)
    v_p <- sub_r$Raw_P_Limma
    status <- "MAPPED"
    
    if (abs(v_fc) < 0.1) {
      cat_label <- "Flat (|log2FC| < 0.1)"
    } else if (sign(v_fc) == sign(d_fc)) {
      cat_label <- if (v_p < 0.05) "Same direction (significant, P < 0.05)" else "Same direction (not significant, P >= 0.05)"
    } else {
      cat_label <- "Discordant"
    }
    c_note <- "Effluent cells (PD duration proxy) vs Tissue (EPS)"
  } else {
    v_fc <- NA
    v_se <- NA
    v_ci <- NA
    v_p <- NA
    status <- "NOT MAPPED ON GPL10558"
    cat_label <- "Not mapped (No probe on array)"
    c_note <- "Unmapped on Illumina HT-12 v4 beadchip"
  }
  
  conc_records[[i]] <- data.frame(
    Gene = g,
    Status = status,
    Discovery_Tissue_EPS_vs_Pooled_log2FC = d_fc,
    Discovery_Tissue_Raw_P = d_p,
    Validation_Effluent_log2FC = v_fc,
    Validation_Effluent_SE_Limma = v_se,
    Validation_Effluent_95_CI = v_ci,
    Validation_Effluent_Raw_P_Limma = v_p,
    Concordance_Category = cat_label,
    Compartment_Note = c_note,
    stringsAsFactors = FALSE
  )
}

df_conc_regen <- do.call(rbind, conc_records)
write.csv(df_conc_regen, "results/tables/C_gse125498_concordance.csv", row.names = FALSE)
cat("\nSaved regenerated results/tables/C_gse125498_concordance.csv\n")
print(df_conc_regen)
