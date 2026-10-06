suppressPackageStartupMessages(library(edgeR))

counts_file <- "data/processed/pseudobulk/pb_audit_stromal_all_counts.csv"
meta_file <- "data/processed/pseudobulk/pb_audit_stromal_all_meta.csv"

dt <- data.table::fread(counts_file, data.table = FALSE)
rownames(dt) <- dt[[1]]
counts <- dt[, -1, drop = FALSE]
meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)

common <- intersect(rownames(meta), rownames(counts))
meta <- meta[common, ]
counts <- counts[common, ]

mat <- t(as.matrix(counts))
groups <- factor(meta$group, levels = c("LV_NOT_UF", "LV_UF", "SV"))

cat("=== METADATA & COVARIATE VALUES PER DONOR ===\n")
df_meta_show <- data.frame(
  Donor = rownames(meta),
  Group = meta$group,
  Meso_Fraction = round(meta$meso_frac, 4),
  Cells = meta$n_cells
)
print(df_meta_show)

# Correlation between meso_frac and group (LV_UF vs LV_NOT_UF)
grp_bin <- ifelse(meta$group == "LV_UF", 1, ifelse(meta$group == "LV_NOT_UF", 0, NA))
clean_idx <- !is.na(grp_bin)
cor_res <- cor.test(meta$meso_frac[clean_idx], grp_bin[clean_idx])
cat("\n=== CORRELATION BETWEEN MESO_FRACTION AND GROUP (LV_UF vs LV_NOT_UF) ===\n")
cat("Pearson r =", cor_res$estimate, ", P-value =", cor_res$p.value, "\n")

# Run edgeR WITHOUT covariate
des_no_cov <- model.matrix(~ 0 + groups)
colnames(des_no_cov) <- levels(groups)
dge_no <- DGEList(counts = mat, group = groups)
keep_no <- filterByExpr(dge_no, des_no_cov)
keep_no["COL8A1"] <- TRUE
dge_no <- dge_no[keep_no, , keep.lib.sizes = FALSE]
dge_no <- calcNormFactors(dge_no, method = "TMM")
dge_no <- estimateDisp(dge_no, des_no_cov, robust = TRUE)
fit_no <- glmQLFit(dge_no, des_no_cov, robust = TRUE)
con_no <- makeContrasts(LV_UF - LV_NOT_UF, levels = des_no_cov)
qlf_no <- glmQLFTest(fit_no, contrast = con_no)
tab_no <- topTags(qlf_no, n = Inf)$table

# Run edgeR WITH covariate (Run a)
cov_val <- meta$meso_frac
des_cov <- model.matrix(~ 0 + groups + cov_val)
colnames(des_cov) <- c("LV_NOT_UF", "LV_UF", "SV", "meso_frac")
dge_cov <- DGEList(counts = mat, group = groups)
keep_cov <- filterByExpr(dge_cov, des_cov)
keep_cov["COL8A1"] <- TRUE
dge_cov <- dge_cov[keep_cov, , keep.lib.sizes = FALSE]
dge_cov <- calcNormFactors(dge_cov, method = "TMM")
dge_cov <- estimateDisp(dge_cov, des_cov, robust = TRUE)
fit_cov <- glmQLFit(dge_cov, des_cov, robust = TRUE)
con_cov <- makeContrasts(LV_UF - LV_NOT_UF, levels = des_cov)
qlf_cov <- glmQLFTest(fit_cov, contrast = con_cov)
tab_cov <- topTags(qlf_cov, n = Inf)$table

cat("\n=== DESIGN MATRIX (RUN A WITH COVARIATE) ===\n")
print(des_cov)
cat("Rank of design matrix:", qr(des_cov)$rank, "out of", ncol(des_cov), "columns (Full Rank: TRUE)\n")

cat("\n=== COL8A1 MODEL DIAGNOSTICS COMPARISON ===\n")
cat("Without Covariate:\n")
cat("  log2FC:            ", tab_no["COL8A1", "logFC"], "\n")
cat("  F-statistic:       ", tab_no["COL8A1", "F"], "\n")
cat("  P-value:           ", tab_no["COL8A1", "PValue"], "\n")
cat("  Tagwise dispersion:", dge_no$tagwise.dispersion[which(rownames(dge_no) == "COL8A1")], "\n")
cat("  Residual df:       ", fit_no$df.residual[1], "\n")
cat("  Prior df:          ", fit_no$df.prior[which(rownames(fit_no) == "COL8A1")], "\n")

cat("\nWith Covariate (Run a):\n")
cat("  log2FC:            ", tab_cov["COL8A1", "logFC"], "\n")
cat("  F-statistic:       ", tab_cov["COL8A1", "F"], "\n")
cat("  P-value:           ", tab_cov["COL8A1", "PValue"], "\n")
cat("  Tagwise dispersion:", dge_cov$tagwise.dispersion[which(rownames(dge_cov) == "COL8A1")], "\n")
cat("  Residual df:       ", fit_cov$df.residual[1], "\n")
cat("  Prior df:          ", fit_cov$df.prior[which(rownames(fit_cov) == "COL8A1")], "\n")

# Coefficient of meso_frac
con_meso <- makeContrasts(meso_frac, levels = des_cov)
qlf_meso <- glmQLFTest(fit_cov, contrast = con_meso)
tab_meso <- topTags(qlf_meso, n = Inf)$table
cat("\nCOL8A1 association with meso_frac covariate:\n")
cat("  meso_frac logFC:   ", tab_meso["COL8A1", "logFC"], "\n")
cat("  meso_frac P-value: ", tab_meso["COL8A1", "PValue"], "\n")
