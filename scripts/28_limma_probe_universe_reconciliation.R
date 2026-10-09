
library(limma)
library(GEOquery)

# Load raw series matrix for GSE125498
gse_file <- 'data/raw/GSE125498_series_matrix.txt.gz'
gse <- getGEO(filename=gse_file, GSEMatrix=TRUE, AnnotGPL=FALSE)
mat_full <- exprs(gse) # 47323 x 33
pdata <- pData(gse)

# Classify groups
spd_mask <- grepl("short", pdata$characteristics_ch1.1, ignore.case=TRUE) | grepl("SPD", pdata$title) | grepl("short", pdata$title)
group <- factor(ifelse(spd_mask, "SPD", "LPD"), levels=c("SPD", "LPD"))
design <- model.matrix(~ group)

# Full array eBayes
fit_full <- lmFit(mat_full, design)
eb_full <- eBayes(fit_full)

# Filtered 19,164 universe (from Validation/GSE125498.top.table.tsv)
top_tab <- read.table('Validation/GSE125498.top.table.tsv', sep='\t', header=TRUE)
valid_probes <- intersect(rownames(mat_full), top_tab$ID)
mat_19k <- mat_full[valid_probes, ]

fit_19k <- lmFit(mat_19k, design)
eb_19k <- eBayes(fit_19k)

# Extract VCAN (ILMN_1687301)
vcan_id <- "ILMN_1687301"

se_full <- fit_full$stdev.unscaled[vcan_id, 2] * sqrt(eb_full$s2.post[vcan_id])
se_19k <- fit_19k$stdev.unscaled[vcan_id, 2] * sqrt(eb_19k$s2.post[vcan_id])

df_res <- data.frame(
  Universe = c("Full Array (47,323 probes)", "Filtered Top Table (19,164 probes)"),
  N_Probes = c(nrow(mat_full), nrow(mat_19k)),
  df_residual = c(fit_full$df.residual[1], fit_19k$df.residual[1]),
  df_prior = c(eb_full$df.prior, eb_19k$df.prior),
  s2_prior = c(eb_full$s2.prior, eb_19k$s2.prior),
  VCAN_log2FC = c(fit_full$coefficients[vcan_id, 2], fit_19k$coefficients[vcan_id, 2]),
  VCAN_SE = c(se_full, se_19k),
  VCAN_t = c(eb_full$t[vcan_id, 2], eb_19k$t[vcan_id, 2]),
  VCAN_P = c(eb_full$p.value[vcan_id, 2], eb_19k$p.value[vcan_id, 2])
)

write.csv(df_res, 'results/tables/C2_limma_probe_universe_comparison.csv', index=FALSE)
print(df_res)
