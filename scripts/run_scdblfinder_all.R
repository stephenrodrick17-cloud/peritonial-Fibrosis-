# Run scDblFinder across all 16 GSE248762 samples independently (hub-gene blind)
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

extract_dir <- "data/raw/GSE248762_extracted"
out_csv <- "results/tables/scDblFinder_per_barcode_calls.csv"

# Blocked hub genes (strict guard)
blocked_hubs <- c("ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", 
                  "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX")

files <- list.files(extract_dir, full.names = TRUE)
barcodes_files <- sort(grep("barcodes\\.tsv\\.gz$", files, value = TRUE))

results_list <- list()

cat("Starting scDblFinder across 16 samples...\n")
set.seed(42)

for (b_file in barcodes_files) {
  base <- sub("\\.barcodes\\.tsv\\.gz$", "", basename(b_file))
  parts <- strsplit(base, "_", fixed = TRUE)[[1]]
  gsm <- parts[1]
  sample_title <- paste(parts[-1], collapse = "_")
  
  feat_file <- file.path(extract_dir, paste0(base, ".features.tsv.gz"))
  mat_file <- file.path(extract_dir, paste0(base, ".matrix.mtx.gz"))
  
  # Read barcodes
  bcs <- readLines(gzfile(b_file))
  # Read features
  feats <- read.delim(gzfile(feat_file), header = FALSE, stringsAsFactors = FALSE)
  gene_symbols <- feats$V2
  
  # Verify hub gene exclusion
  hub_idx <- which(toupper(gene_symbols) %in% blocked_hubs)
  
  # Read matrix (genes x cells)
  mat <- readMM(gzfile(mat_file))
  rownames(mat) <- gene_symbols
  colnames(mat) <- paste0(gsm, "_", bcs)
  
  # Drop hub genes completely before passing to scDblFinder
  if (length(hub_idx) > 0) {
    mat <- mat[-hub_idx, ]
  }
  
  cat(sprintf("Running %s (%s): %d cells x %d genes (hub genes removed: %d)...\n", 
              gsm, sample_title, ncol(mat), nrow(mat), length(hub_idx)))
  
  sce <- SingleCellExperiment(assays = list(counts = mat))
  # Run scDblFinder with default settings
  sce <- scDblFinder(sce, verbose = FALSE)
  
  df <- data.frame(
    barcode = colnames(sce),
    gsm = gsm,
    sample = sample_title,
    scDblFinder_score = sce$scDblFinder.score,
    scDblFinder_class = as.character(sce$scDblFinder.class),
    stringsAsFactors = FALSE
  )
  results_list[[gsm]] <- df
  
  n_dbl <- sum(df$scDblFinder_class == "doublet")
  rate <- n_dbl / nrow(df) * 100
  cat(sprintf("  -> %s doublets: %d / %d (%.2f%%)\n", sample_title, n_dbl, nrow(df), rate))
}

all_res <- do.call(rbind, results_list)
write.csv(all_res, out_csv, row.names = FALSE)
cat(sprintf("scDblFinder complete! Total barcodes: %d, Total doublets: %d (%.2f%%). Saved to %s\n",
            nrow(all_res), sum(all_res$scDblFinder_class == "doublet"),
            sum(all_res$scDblFinder_class == "doublet") / nrow(all_res) * 100, out_csv))
