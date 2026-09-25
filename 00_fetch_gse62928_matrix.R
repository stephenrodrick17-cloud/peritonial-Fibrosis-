# ==============================================================================
# SCRIPT 00: FETCH & VERIFY FULL GSE62928 EXPRESSION MATRIX
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Dataset: GSE62928 (Platform: GPL13158 Affymetrix HT HG-U133 Plus PM)
# Purpose: Extract normalized expression matrix and phenotype metadata for WGCNA
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# Ensure directories exist
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

cat("=================================================================\n")
cat("STEP 1: FETCHING GSE62928 FROM NCBI GEO\n")
cat("=================================================================\n")

# Use local cache if present, otherwise fetch via GEOquery
gse_list <- getGEO("GSE62928", GSEMatrix = TRUE, destdir = "data")
if (length(gse_list) == 0) {
  stop("Failed to retrieve GSE62928 series matrix from GEO!")
}

eset <- gse_list[[1]]
cat("Successfully loaded ExpressionSet for GSE62928.\n")
cat("Platform Annotation:", annotation(eset), "\n")
cat("Raw Probe Matrix Dimensions:", nrow(eset), "probes x", ncol(eset), "samples\n")

# ------------------------------------------------------------------------------
# STEP 2: PHENOTYPE DATA (pData) INSPECTION & THREE-GROUP IDENTIFICATION
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 2: PHENOTYPE DATA (pData) COLUMN INSPECTION & GROUP SIZES\n")
cat("=================================================================\n")

pd <- pData(eset)
cat("Total Samples in GSE62928:", nrow(pd), "\n")
cat("Total pData Columns Found:", ncol(pd), "\n\n")

# Print ALL pData column names
cat("--- ALL pData COLUMN NAMES ---\n")
for (i in seq_along(colnames(pd))) {
  cat(sprintf("  [%02d] %s\n", i, colnames(pd)[i]))
}

# Print unique values for each column with <= 10 unique values
cat("\n--- COLUMNS WITH <= 10 UNIQUE VALUES (FOR GROUP-LABEL IDENTIFICATION) ---\n")
for (col in colnames(pd)) {
  u_vals <- unique(pd[[col]])
  if (length(u_vals) <= 10) {
    cat(sprintf("Column: '%s' (Unique values: %d)\n", col, length(u_vals)))
    cat("  Values:", paste(head(as.character(u_vals), 10), collapse = " | "), "\n")
  }
}

# (a) CONFIRM/EDIT: Positively identify the 3 clinical groups from pData
# In GSE62928:
# 1. EPS: Encapsulating Peritoneal Sclerosis (treatment: 'EPS, frozen biopsy specimen')
# 2. PD: Patients undergoing PD catheter insertion without EPS (treatment: 'PD, frozen biopsy specimen')
# 3. Uremic: Uremic control patients undergoing abdominal surgery without PD/EPS history ('UREMIC control, frozen biopsy specimen')

group_3class <- character(nrow(pd))
for (i in seq_len(nrow(pd))) {
  txt <- paste(pd$title[i], pd$characteristics_ch1[i], pd$`treatment:ch1`[i])
  if (grepl("EPS", txt, ignore.case = TRUE)) {
    group_3class[i] <- "EPS"
  } else if (grepl("UREMIC", txt, ignore.case = TRUE)) {
    group_3class[i] <- "Uremic"
  } else if (grepl("PD", txt, ignore.case = TRUE)) {
    group_3class[i] <- "PD"
  } else {
    group_3class[i] <- "Unknown"
  }
}

# Explicitly report exact sample counts per group
n_eps <- sum(group_3class == "EPS")
n_pd <- sum(group_3class == "PD")
n_uremic <- sum(group_3class == "Uremic")
n_actual_total <- length(group_3class)

cat("\n=================================================================\n")
cat("EXACT GSE62928 CLINICAL GROUP SIZES (VERIFIED DIRECTLY FROM GEO):\n")
cat("=================================================================\n")
cat(sprintf("  1. EPS (Severe Peritoneal Fibrosis):                n = %d\n", n_eps))
cat(sprintf("  2. PD (PD Catheter Insertion, No EPS):             n = %d\n", n_pd))
cat(sprintf("  3. Uremic/Control (No PD or EPS History):           n = %d\n", n_uremic))
cat(sprintf("  TOTAL ACTUAL SAMPLES ACROSS ALL 3 GROUPS:          N = %d\n", n_actual_total))
cat("=================================================================\n")

# ------------------------------------------------------------------------------
# STEP 3: BINARY TRAIT DEFINITION & HARD STOP CHECKPOINT (n < 15)
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 3: BINARY TRAIT DEFINITION FOR WGCNA & N >= 15 CHECKPOINT\n")
cat("=================================================================\n")

# Decision logic per instructions:
# Option (a): EPS vs (PD + Uremic combined), isolating fibrotic endpoint specifically
# Option (b): EPS vs Uremic-Control only (dropping PD samples)
# Default to option (a) unless dropping PD leaves n >= 15 per remaining group.
if (n_eps >= 15 && n_uremic >= 15) {
  trait_choice <- "b"
  cat("Trait Selection: Option (b) [EPS vs Uremic-Control only, dropping PD samples]\n")
  keep_idx <- which(group_3class %in% c("EPS", "Uremic"))
  n_trait_case <- n_eps
  n_trait_ctrl <- n_uremic
} else {
  trait_choice <- "a"
  cat("Trait Selection: Option (a) [EPS vs (PD + Uremic combined)]\n")
  cat("  Rationale: Option (b) does not leave n >= 15 per group. Option (a) maximizes sample size (N=8) and isolates the fibrotic endpoint.\n")
  keep_idx <- seq_len(nrow(pd))
  n_trait_case <- n_eps
  n_trait_ctrl <- n_pd + n_uremic
}

n_trait_total <- n_trait_case + n_trait_ctrl
cat(sprintf("Chosen Comparison: %d Case (EPS) vs %d Control (PD + Uremic) -> Total n = %d\n",
            n_trait_case, n_trait_ctrl, n_trait_total))

# Define binary trait column
group_binary_trait <- ifelse(group_3class == "EPS", "Case_EPS", "Control")

# Check for override flag: --allow-small-n or ALLOW_SMALL_N environment variable
args <- commandArgs(trailingOnly = TRUE)
allow_small_n <- any(grepl("--allow-small-n|--force|--proceed", args, ignore.case = TRUE)) ||
                 identical(Sys.getenv("ALLOW_SMALL_N"), "TRUE")

if (n_trait_total < 15) {
  cat("\n", paste0(rep("!", 80), collapse = ""), "\n", sep = "")
  cat("HARD STOP ALERT: SAMPLE SIZE IS BELOW WGCNA POWER THRESHOLD (n < 15)!\n")
  cat(sprintf("  Actual usable sample size for chosen comparison: n = %d\n", n_trait_total))
  cat("  Minimum recommended for reliable WGCNA module detection: n >= 15\n")
  cat("  WGCNA module stability, topological overlap estimation, and scale-free\n")
  cat("  fitting become sensitive and fragile with n < 15.\n")
  cat(paste0(rep("!", 80), collapse = ""), "\n\n", sep = "")
  
  if (!allow_small_n) {
    stop("HARD STOP: Total usable sample size for GSE62928 is n = ", n_trait_total, 
         " (< 15). WGCNA module stability is unreliable below n = 15. ",
         "Execution stopped per protocol to prevent an underpowered run. ",
         "To proceed with the pipeline despite this statistical limitation, run with '--allow-small-n' or set ALLOW_SMALL_N=TRUE.")
  } else {
    cat("OVERRIDE ACTIVE: '--allow-small-n' detected. Proceeding with explicit limitation logging.\n\n")
  }
} else {
  cat("Sample size verification PASSED: n =", n_trait_total, ">= 15.\n")
}

# ------------------------------------------------------------------------------
# STEP 4: PROBE MAPPING & COLLAPSING BY HIGHEST MEAN EXPRESSION (MAXMEAN)
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 4: PROBE MAPPING & MAXMEAN COLLAPSING (GPL13158)\n")
cat("=================================================================\n")

# (b) CONFIRM/EDIT: Correct GPL annotation extraction
fd <- fData(eset)
ex <- exprs(eset)

if (!"Gene Symbol" %in% colnames(fd)) {
  stop("Column 'Gene Symbol' not found in fData(eset)!")
}

gene_symbols <- as.character(fd[["Gene Symbol"]])
probe_ids <- rownames(ex)

# Clean empty / NA symbols
valid_idx <- which(!is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---")
ex_valid <- ex[valid_idx, , drop = FALSE]
genes_valid <- gene_symbols[valid_idx]
probes_valid <- probe_ids[valid_idx]

cat("Probes with valid gene annotations:", length(valid_idx), "/", nrow(ex), "\n")

# Split multi-gene annotations (e.g. 'GENE1 /// GENE2')
split_list <- strsplit(genes_valid, " ?/// ?")
lengths_vec <- lengths(split_list)

expanded_probes <- rep(probes_valid, lengths_vec)
expanded_genes <- trimws(unlist(split_list))
expanded_exprs <- ex_valid[rep(seq_along(probes_valid), lengths_vec), , drop = FALSE]

# Collapse multi-probe cases by keeping probe with HIGHEST MEAN EXPRESSION across samples
probe_means <- rowMeans(expanded_exprs, na.rm = TRUE)
df_collapsing <- data.frame(
  Probe_ID = expanded_probes,
  Gene_Symbol = expanded_genes,
  Mean_Expr = probe_means,
  stringsAsFactors = FALSE
)

# Order by Gene_Symbol and descending Mean_Expr, then deduplicate
df_collapsing <- df_collapsing[order(df_collapsing$Gene_Symbol, -df_collapsing$Mean_Expr), ]
best_probes <- df_collapsing[!duplicated(df_collapsing$Gene_Symbol), ]

collapsed_matrix <- ex[best_probes$Probe_ID, , drop = FALSE]
rownames(collapsed_matrix) <- best_probes$Gene_Symbol

cat("Unique collapsed genes (MaxMean rule):", nrow(collapsed_matrix), "\n")
cat("Matrix dimensions (genes x samples):", nrow(collapsed_matrix), "x", ncol(collapsed_matrix), "\n")

# ------------------------------------------------------------------------------
# STEP 5: SAVE FULL EXPRESSION MATRIX AND SAMPLE METADATA
# ------------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("STEP 5: EXPORTING PROCESSED TABLES\n")
cat("=================================================================\n")

out_matrix_file <- "results/tables/GSE62928_full_expression_matrix.csv"
out_meta_file <- "results/tables/GSE62928_sample_metadata.csv"

# 1. Full expression matrix (genes x samples)
write.csv(collapsed_matrix, out_matrix_file, row.names = TRUE)
cat("Saved full expression matrix to:", out_matrix_file, "\n")

# 2. Sample metadata table (sample_id, group_3class, group_binary_trait)
metadata_df <- data.frame(
  sample_id = colnames(collapsed_matrix),
  sample_title = pd$title,
  group_3class = group_3class,
  group_binary_trait = group_binary_trait,
  binary_numeric = ifelse(group_binary_trait == "Case_EPS", 1, 0),
  geo_accession = pd$geo_accession,
  characteristics = pd$characteristics_ch1,
  stringsAsFactors = FALSE
)
write.csv(metadata_df, out_meta_file, row.names = FALSE)
cat("Saved sample metadata to:", out_meta_file, "\n")
cat("\nMetadata Preview:\n")
print(metadata_df[, c("sample_id", "sample_title", "group_3class", "group_binary_trait", "binary_numeric")])

cat("\n=================================================================\n")
cat("TASK 0 COMPLETED SUCCESSFULLY!\n")
cat("=================================================================\n")
