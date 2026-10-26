suppressPackageStartupMessages({
  library(limma)
})

set.seed(42)

expr_file <- file.path("results", "tables", "GSE125498_full_expression_matrix.csv")
meta_file <- file.path("results", "tables", "GSE125498_sample_metadata.csv")
matri_file <- file.path("results", "revision", "matrisome_volcano", "GSE125498_matrisome_volcano_source_data.csv")

expr <- read.csv(expr_file, row.names=1, check.names=FALSE)
meta <- read.csv(meta_file, stringsAsFactors=FALSE)
matri <- read.csv(matri_file, stringsAsFactors=FALSE)

sample_ids <- meta$sample_id
group_labels <- meta$group
X <- as.matrix(expr[, sample_ids])

matri_genes <- intersect(matri$Gene, rownames(X))
cat("Total genes:", nrow(X), "\n")
cat("Matrisome genes:", length(matri_genes), "\n")

# Actual DEG counts in observed data
# Re-running limma once on true labels to get observed cutoffs:
design_true <- model.matrix(~ factor(group_labels, levels=c("SPD", "LPD")))
fit_true <- eBayes(lmFit(X, design_true))
tt_true <- topTable(fit_true, coef=2, number=Inf, sort.by="none")
tt_matri_true <- tt_true[matri_genes, ]

obs_all <- sum(abs(tt_matri_true$logFC) >= 0.585 & tt_matri_true$P.Value < 0.05)
obs_up <- sum(tt_matri_true$logFC >= 0.585 & tt_matri_true$P.Value < 0.05)
obs_down <- sum(tt_matri_true$logFC <= -0.585 & tt_matri_true$P.Value < 0.05)

cat(sprintf("Observed Matrisome DEGs in GSE125498: All=%d, Up=%d, Down=%d\n", obs_all, obs_up, obs_down))

# Run 10,000 permutations
N_PERM <- 10000
count_extreme_all <- 0
count_extreme_up <- 0
count_extreme_down <- 0

# To make it fast, notice that in lmFit for a 2-group design:
# logFC is simply the difference between group means!
# eBayes only shrinks the residual variances across all 11,741 genes.
# We will run full lmFit + eBayes each time.
cat("Starting 10,000 limma permutations (seed=42)...\n")
t0 <- Sys.time()

for (i in 1:N_PERM) {
  shuff_grp <- sample(group_labels)
  design <- model.matrix(~ factor(shuff_grp, levels=c("SPD", "LPD")))
  fit <- eBayes(lmFit(X, design))
  
  # We only need p-values and logFC for matri_genes
  # In limma, coefficients are in fit$coefficients[, 2]
  # p-values are in fit$p.value[, 2]
  lfc <- fit$coefficients[matri_genes, 2]
  pval <- fit$p.value[matri_genes, 2]
  
  n_all <- sum(abs(lfc) >= 0.585 & pval < 0.05)
  n_up <- sum(lfc >= 0.585 & pval < 0.05)
  n_down <- sum(lfc <= -0.585 & pval < 0.05)
  
  if (n_all >= obs_all) count_extreme_all <- count_extreme_all + 1
  if (n_up >= obs_up) count_extreme_up <- count_extreme_up + 1
  if (n_down >= obs_down) count_extreme_down <- count_extreme_down + 1
  
  if (i %% 1000 == 0) {
    t_curr <- Sys.time()
    elapsed <- as.numeric(t_curr - t0, units="secs")
    cat(sprintf("Completed %d / %d perms (elapsed: %.1fs, ~%.1fs remaining)\n", 
                i, N_PERM, elapsed, elapsed / i * (N_PERM - i)))
  }
}

p_emp_all <- (count_extreme_all + 1) / (N_PERM + 1)
p_emp_up <- (count_extreme_up + 1) / (N_PERM + 1)
p_emp_down <- (count_extreme_down + 1) / (N_PERM + 1)

res_df <- data.frame(
  Test = c("All Matrisome DEGs", "Upregulated Matrisome DEGs", "Downregulated Matrisome DEGs"),
  Observed_Count = c(obs_all, obs_up, obs_down),
  Permutations_ge_Observed = c(count_extreme_all, count_extreme_up, count_extreme_down),
  Total_Permutations = N_PERM,
  Empirical_Permutation_P = c(p_emp_all, p_emp_up, p_emp_down),
  Random_Seed = 42
)

write.csv(res_df, "results/tables/matrisome_limma_10k_permutation_results.csv", row.names=FALSE)
cat("\n=== Final 10,000 Limma Permutation Results ===\n")
print(res_df)
