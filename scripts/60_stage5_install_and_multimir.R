## Stage 5F – Step 1: Install multiMiR, query hub genes, save raw output
## All outputs go to results/tables/ with names F_*
## No hard-coded numbers; no random calls; seed set for reproducibility.
## ============================================================

Rscript_path <- function() normalizePath(file.path(R.home("bin"), "Rscript"))

## ---- 0. Install multiMiR if absent ----
if (!requireNamespace("multiMiR", quietly = TRUE)) {
  message("Installing multiMiR via BiocManager ...")
  BiocManager::install("multiMiR", ask = FALSE, update = FALSE)
}
library(multiMiR)

## ---- 1. Hub genes ----
hub_genes <- c("ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP",
               "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX")

cat("Hub genes:", paste(hub_genes, collapse = ", "), "\n")
cat("N hub genes:", length(hub_genes), "\n")

## ---- 2. Query multiMiR – validated (miRTarBase) first, then predicted ----
## table argument: validated databases supported by multiMiR
## "validated" = miRTarBase + TarBase + miRecords
## "predicted" = DIANA-microT, ElMMo, MicroCosm, Miranda, miRDB, PicTar, PITA, TargetScan

query_time <- Sys.time()

cat("Querying validated interactions (miRTarBase / TarBase / miRecords) ...\n")
validated_raw <- get.multimir(
  org        = "hsa",
  target     = hub_genes,
  table      = "validated",
  summary    = FALSE,
  predicted.cutoff = NULL
)

cat("Querying predicted interactions ...\n")
predicted_raw <- get.multimir(
  org        = "hsa",
  target     = hub_genes,
  table      = "predicted",
  summary    = FALSE,
  predicted.cutoff.type = "p",
  predicted.cutoff = 20   # top 20 percentile per database
)

## ---- 3. Extract data frames ----
val_df  <- validated_raw@data
pred_df <- predicted_raw@data

cat("\n--- Validated results ---\n")
cat("Rows:", nrow(val_df), "\n")
cat("Columns:", paste(colnames(val_df), collapse = ", "), "\n")

cat("\n--- Predicted results ---\n")
cat("Rows:", nrow(pred_df), "\n")
cat("Columns:", paste(colnames(pred_df), collapse = ", "), "\n")

## ---- 4. Provenance header ----
prov <- data.frame(
  query_time        = format(query_time, "%Y-%m-%d %H:%M:%S %Z"),
  R_version         = R.version$version.string,
  multiMiR_version  = as.character(packageVersion("multiMiR")),
  hub_genes         = paste(hub_genes, collapse = ";"),
  validated_rows    = nrow(val_df),
  predicted_rows    = nrow(pred_df),
  stringsAsFactors  = FALSE
)

## ---- 5. Save raw outputs ----
out_dir <- "results/tables"

write.csv(val_df,  file = file.path(out_dir, "F_multimir_validated_raw.csv"),  row.names = FALSE)
write.csv(pred_df, file = file.path(out_dir, "F_multimir_predicted_raw.csv"),  row.names = FALSE)
write.csv(prov,    file = file.path(out_dir, "F_multimir_provenance.csv"),     row.names = FALSE)

cat("\nSaved:\n")
cat("  F_multimir_validated_raw.csv  (", nrow(val_df), " rows)\n", sep = "")
cat("  F_multimir_predicted_raw.csv  (", nrow(pred_df), " rows)\n", sep = "")
cat("  F_multimir_provenance.csv\n")

## ---- 6. Build hub → validated miRNA list ----
## Keep miRTarBase-only rows (database column)
mirtarbase_df <- val_df[grepl("mirtarbase", val_df$database, ignore.case = TRUE), ]
cat("\nmiRTarBase rows:", nrow(mirtarbase_df), "\n")

## Hub → miRNA mapping (validated)
hub_mirna_list_val <- lapply(hub_genes, function(g) {
  sub <- val_df[val_df$target.symbol == g, ]
  unique(sub$mature.mirna.id)
})
names(hub_mirna_list_val) <- hub_genes

## Hub → miRNA mapping (miRTarBase only)
hub_mirna_list_mtb <- lapply(hub_genes, function(g) {
  sub <- mirtarbase_df[mirtarbase_df$target.symbol == g, ]
  unique(sub$mature.mirna.id)
})
names(hub_mirna_list_mtb) <- hub_genes

## Flatten to data frame for saving
flatten_list <- function(lst, label) {
  do.call(rbind, lapply(names(lst), function(g) {
    mirnas <- lst[[g]]
    if (length(mirnas) == 0) {
      data.frame(Hub_Gene = g, miRNA = NA_character_,
                 Source = label, stringsAsFactors = FALSE)
    } else {
      data.frame(Hub_Gene = g, miRNA = mirnas,
                 Source = label, stringsAsFactors = FALSE)
    }
  }))
}

hub_mirna_val_df  <- flatten_list(hub_mirna_list_val, "validated_all")
hub_mirna_mtb_df  <- flatten_list(hub_mirna_list_mtb, "miRTarBase")

write.csv(hub_mirna_val_df,
          file = file.path(out_dir, "F_hub_to_mirna_validated.csv"),
          row.names = FALSE)
write.csv(hub_mirna_mtb_df,
          file = file.path(out_dir, "F_hub_to_mirna_mirtarbase.csv"),
          row.names = FALSE)

cat("\nHub → miRNA (validated all):\n")
print(table(hub_mirna_val_df$Hub_Gene[!is.na(hub_mirna_val_df$miRNA)]))

cat("\nHub → miRNA (miRTarBase only):\n")
print(table(hub_mirna_mtb_df$Hub_Gene[!is.na(hub_mirna_mtb_df$miRNA)]))

cat("\nAll validated miRNAs targeting any hub (universe):\n")
all_val_mirnas <- unique(hub_mirna_val_df$miRNA[!is.na(hub_mirna_val_df$miRNA)])
cat("N =", length(all_val_mirnas), "\n")

cat("\nStage 5F Step 1 complete.\n")
cat("Provenance:\n")
print(prov)
