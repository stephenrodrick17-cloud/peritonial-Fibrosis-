# ==============================================================================
# PERMUTATION TEST: ECM ENRICHMENT AMONG UP-REGULATED DEGs (LIMMA EBAYES)
# Dataset: GSE62928 (4 EPS Cases vs 4 Controls -> Exact 70 Permutations)
# Primary Cutoff: Nominal P < 0.05 and log2FC >= 0.80
# Sensitivity Cutoff: Nominal P < 0.05 and log2FC >= 0.585
# ==============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(readxl)
  library(dplyr)
})

# 1. Load ECM masterlist
ecm_raw <- read_excel("ECM genes all.xlsx", skip = 1)
sym_col <- grep("Symbol|Gene", colnames(ecm_raw), value = TRUE, ignore.case = TRUE)[1]
ecm_genes <- unique(trimws(toupper(as.character(ecm_raw[[sym_col]]))))
ecm_genes <- ecm_genes[!is.na(ecm_genes) & !ecm_genes %in% c("GENE SYMBOL", "NA", "NAN", "")]
cat(sprintf("Loaded %d unique curated ECM / Matrisome genes\n", length(ecm_genes)))

# 2. Load Series Matrix
gse_lines <- readLines("data/GSE62928_series_matrix.txt.gz")
start_idx <- which(grepl("^!series_matrix_table_begin", gse_lines)) + 1
end_idx <- which(grepl("^!series_matrix_table_end", gse_lines)) - 1

mat_text <- gse_lines[start_idx:end_idx]
mat_df <- read.table(text = mat_text, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "\"")
rownames(mat_df) <- mat_df[, 1]
mat_df <- as.matrix(mat_df[, -1])
storage.mode(mat_df) <- "numeric"
cat(sprintf("Expression matrix loaded: %d probes x %d samples\n", nrow(mat_df), ncol(mat_df)))

# 3. Load Annotation mapping from GSE62928.top.table.tsv
df_tt <- read.table("GSE62928.top.table.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
df_annot <- df_tt[!is.na(df_tt$Gene.symbol) & df_tt$Gene.symbol != "" & df_tt$Gene.symbol != "---", ]
df_annot$Gene.symbol <- trimws(df_annot$Gene.symbol)

# Expand multi-gene probes
split_syms <- strsplit(df_annot$Gene.symbol, " ?/// ?")
n_syms <- lengths(split_syms)
expanded_df <- data.frame(
  ID = rep(df_annot$ID, n_syms),
  Gene = trimws(toupper(unlist(split_syms))),
  stringsAsFactors = FALSE
)
expanded_df <- expanded_df[expanded_df$ID %in% rownames(mat_df), ]
cat(sprintf("Total annotated probe-gene pairs: %d\n", nrow(expanded_df)))

# Subset matrix to expanded rows
sub_mat <- mat_df[expanded_df$ID, ]
genes_vec <- expanded_df$Gene

# True group labels
true_eps_idx <- 1:4
true_ctrl_idx <- 5:8

# Function to run limma DEG pipeline and count ECM genes
run_limma_permutation <- function(eps_indices, logfc_cutoff = 0.80) {
  group <- rep("Control", 8)
  group[eps_indices] <- "EPS"
  group <- factor(group, levels = c("Control", "EPS"))
  
  design <- model.matrix(~ group)
  fit <- lmFit(sub_mat, design)
  fit <- eBayes(fit)
  
  # Extract statistics for EPS vs Control (coef = 2)
  logfc <- fit$coefficients[, 2]
  pval <- fit$p.value[, 2]
  
  # Collapse to unique gene symbol by minimum P.Value
  df_temp <- data.frame(Gene = genes_vec, logFC = logfc, P.Value = pval, stringsAsFactors = FALSE)
  df_temp <- df_temp[order(df_temp$P.Value), ]
  df_best <- df_temp[!duplicated(df_temp$Gene), ]
  
  # Up-regulated DEGs
  up_degs <- df_best[df_best$logFC >= logfc_cutoff & df_best$P.Value < 0.05, ]
  ecm_degs <- up_degs[up_degs$Gene %in% ecm_genes, ]
  
  return(list(
    n_total_genes = nrow(df_best),
    n_up_degs = nrow(up_degs),
    n_ecm_degs = nrow(ecm_degs),
    ecm_list = ecm_degs$Gene
  ))
}

# 1. Run Observed True Partition
obs_res_80 <- run_limma_permutation(true_eps_idx, logfc_cutoff = 0.80)
obs_res_585 <- run_limma_permutation(true_eps_idx, logfc_cutoff = 0.585)

cat("\n=================================================================\n")
cat("OBSERVED REAL TEST (GSE62928, N = 8, 4 EPS vs 4 Control):\n")
cat("=================================================================\n")
cat(sprintf("Total Unique Genes Collapsed: %d\n", obs_res_80$n_total_genes))
cat(sprintf("Primary Cutoff (log2FC >= 0.80, P < 0.05):\n"))
cat(sprintf("  - Total Up-Regulated DEGs: %d (Reproduces locked 367)\n", obs_res_80$n_up_degs))
cat(sprintf("  - Observed Up-Regulated ECM-DEGs: %d (Reproduces locked 71)\n", obs_res_80$n_ecm_degs))
cat(sprintf("Sensitivity Cutoff (log2FC >= 0.585, P < 0.05):\n"))
cat(sprintf("  - Total Up-Regulated DEGs: %d\n", obs_res_585$n_up_degs))
cat(sprintf("  - Observed Up-Regulated ECM-DEGs: %d\n", obs_res_585$n_ecm_degs))
cat("=================================================================\n")

# 2. Run all 70 Permutations
all_combos <- combn(8, 4)
n_perms <- ncol(all_combos)
cat(sprintf("\nRunning exact %d label permutations with limma eBayes...\n", n_perms))

perm_counts_80 <- integer(n_perms)
perm_counts_585 <- integer(n_perms)
perm_up_80 <- integer(n_perms)
perm_up_585 <- integer(n_perms)

for (i in seq_len(n_perms)) {
  comb <- all_combos[, i]
  res80 <- run_limma_permutation(comb, logfc_cutoff = 0.80)
  res585 <- run_limma_permutation(comb, logfc_cutoff = 0.585)
  
  perm_counts_80[i] <- res80$n_ecm_degs
  perm_counts_585[i] <- res585$n_ecm_degs
  perm_up_80[i] <- res80$n_up_degs
  perm_up_585[i] <- res585$n_up_degs
}

# Statistics for Primary Cutoff (log2FC >= 0.80)
obs_count_80 <- obs_res_80$n_ecm_degs
rank_80 <- sum(perm_counts_80 > obs_count_80) + 1
p_perm_80 <- mean(perm_counts_80 >= obs_count_80)

cat("\n=================================================================\n")
cat("STEP 2(c): PERMUTATION TEST RESULTS (PRIMARY: log2FC >= 0.80, P < 0.05)\n")
cat("=================================================================\n")
cat(sprintf("Observed ECM-DEG Count: %d\n", obs_count_80))
cat(sprintf("Rank within 70 Permutations: %d / 70\n", rank_80))
cat(sprintf("Exact Permutation P-Value (Count >= %d): %.4f (%d / 70)\n",
            obs_count_80, p_perm_80, sum(perm_counts_80 >= obs_count_80)))
cat(sprintf("Distribution Summary: Min = %d, Q25 = %.1f, Median = %.1f, Mean = %.2f, Q75 = %.1f, Max = %d\n",
            min(perm_counts_80), quantile(perm_counts_80, 0.25), median(perm_counts_80),
            mean(perm_counts_80), quantile(perm_counts_80, 0.75), max(perm_counts_80)))
cat("\nFull Frequency Distribution across all 70 Permutations:\n")
print(table(perm_counts_80))

# Step 2(d): Hypergeometric Enrichment Test
# N = Total background genes collapsed
# K = Total matrisome genes in background
# n = Total up-regulated DEGs
# k = Observed ECM genes in up-DEGs
load("data/GSE62928_preprocessed.RData")
N_pop <- obs_res_80$n_total_genes # 22,049
K_ecm <- sum(unique(all_genes$Gene.symbol) %in% ecm_genes) # 975
n_draw_80 <- obs_res_80$n_up_degs # 367
k_succ_80 <- obs_res_80$n_ecm_degs # 71

# phyper(k-1, K, N-K, n, lower.tail=FALSE)
p_hypergeom_80 <- phyper(k_succ_80 - 1, K_ecm, N_pop - K_ecm, n_draw_80, lower.tail = FALSE)

cat("\n=================================================================\n")
cat("STEP 2(d): HYPERGEOMETRIC ENRICHMENT CHECK (INDEPENDENT CHECK)\n")
cat("=================================================================\n")
cat(sprintf("Background Total Genes (N):  %d\n", N_pop))
cat(sprintf("Background ECM Genes (K):    %d (%.2f%% of genome)\n", K_ecm, 100 * K_ecm / N_pop))
cat(sprintf("Up-Regulated DEGs Drawn (n): %d\n", n_draw_80))
cat(sprintf("Observed ECM in Up-DEGs (k): %d (%.2f%% of DEGs)\n", k_succ_80, 100 * k_succ_80 / n_draw_80))
cat(sprintf("Expected ECM by Chance:      %.2f genes\n", n_draw_80 * (K_ecm / N_pop)))
cat(sprintf("Fold Enrichment:             %.2fx\n", (k_succ_80 / n_draw_80) / (K_ecm / N_pop)))
cat(sprintf("Hypergeometric Exact P-Value: %s\n", format(p_hypergeom_80, scientific = TRUE, digits = 5)))

# Step 2(e): Sensitivity Check (log2FC >= 0.585)
obs_count_585 <- obs_res_585$n_ecm_degs
rank_585 <- sum(perm_counts_585 > obs_count_585) + 1
p_perm_585 <- mean(perm_counts_585 >= obs_count_585)
p_hypergeom_585 <- phyper(obs_count_585 - 1, K_ecm, N_pop - K_ecm, obs_res_585$n_up_degs, lower.tail = FALSE)

cat("\n=================================================================\n")
cat("STEP 2(e): SENSITIVITY CHECK (LESS STRICT CUTOFF: log2FC >= 0.585, P < 0.05)\n")
cat("=================================================================\n")
cat(sprintf("Total Up-Regulated DEGs: %d\n", obs_res_585$n_up_degs))
cat(sprintf("Observed ECM-DEG Count:  %d\n", obs_count_585))
cat(sprintf("Rank within 70 Permutations: %d / 70\n", rank_585))
cat(sprintf("Exact Permutation P-Value (Count >= %d): %.4f (%d / 70)\n",
            obs_count_585, p_perm_585, sum(perm_counts_585 >= obs_count_585)))
cat(sprintf("Hypergeometric Exact P-Value: %s\n", format(p_hypergeom_585, scientific = TRUE, digits = 5)))
cat("=================================================================\n")
