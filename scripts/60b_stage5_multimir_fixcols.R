## Stage 5F – Step 1b: Fix column names in multiMiR output and rebuild hub→miRNA list
## Column names in actual output: mature_mirna_id, target_symbol (underscores, not dots)

library(multiMiR)
library(utils)

out_dir  <- "results/tables"
hub_genes <- c("ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP",
               "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX")

## Load already-saved raw CSVs (avoid re-querying)
val_df  <- read.csv(file.path(out_dir, "F_multimir_validated_raw.csv"),
                    stringsAsFactors = FALSE)
pred_df <- read.csv(file.path(out_dir, "F_multimir_predicted_raw.csv"),
                    stringsAsFactors = FALSE)

cat("Validated columns:", paste(colnames(val_df), collapse = ", "), "\n")
cat("Predicted columns:", paste(colnames(pred_df), collapse = ", "), "\n")

cat("Validated rows:", nrow(val_df), "\n")
cat("Predicted rows:", nrow(pred_df), "\n")

## Determine correct column names dynamically
mirna_col   <- grep("mirna_id|mirna.id", colnames(val_df), value = TRUE, ignore.case = TRUE)[1]
target_col  <- grep("target_symbol|target.symbol", colnames(val_df), value = TRUE, ignore.case = TRUE)[1]
db_col      <- grep("^database$", colnames(val_df), value = TRUE, ignore.case = TRUE)[1]

cat("\nmiRNA column:  ", mirna_col, "\n")
cat("Target column: ", target_col, "\n")
cat("Database col:  ", db_col, "\n")

## ---- miRTarBase-only rows ----
mirtarbase_df <- val_df[grepl("mirtarbase", val_df[[db_col]], ignore.case = TRUE), ]
cat("\nmiRTarBase rows:", nrow(mirtarbase_df), "\n")
cat("All databases in validated:\n")
print(table(val_df[[db_col]]))

## ---- Hub → miRNA list (validated all) ----
hub_mirna_val_df <- do.call(rbind, lapply(hub_genes, function(g) {
  sub    <- val_df[val_df[[target_col]] == g, ]
  mirnas <- unique(sub[[mirna_col]])
  mirnas <- mirnas[!is.na(mirnas) & nchar(mirnas) > 0]
  if (length(mirnas) == 0) {
    data.frame(Hub_Gene = g, miRNA = NA_character_,
               Source = "validated_all", N_interactions = 0L,
               stringsAsFactors = FALSE)
  } else {
    data.frame(Hub_Gene = g, miRNA = mirnas,
               Source = "validated_all",
               N_interactions = nrow(sub),
               stringsAsFactors = FALSE)
  }
}))

## ---- Hub → miRNA list (miRTarBase only) ----
hub_mirna_mtb_df <- do.call(rbind, lapply(hub_genes, function(g) {
  sub    <- mirtarbase_df[mirtarbase_df[[target_col]] == g, ]
  mirnas <- unique(sub[[mirna_col]])
  mirnas <- mirnas[!is.na(mirnas) & nchar(mirnas) > 0]
  if (length(mirnas) == 0) {
    data.frame(Hub_Gene = g, miRNA = NA_character_,
               Source = "miRTarBase", N_interactions = 0L,
               stringsAsFactors = FALSE)
  } else {
    data.frame(Hub_Gene = g, miRNA = mirnas,
               Source = "miRTarBase",
               N_interactions = nrow(sub),
               stringsAsFactors = FALSE)
  }
}))

cat("\n=== Hub → validated miRNA (all sources) ===\n")
print(hub_mirna_val_df)

cat("\n=== Hub → miRTarBase miRNA ===\n")
print(hub_mirna_mtb_df)

## Universe of validated miRNAs
all_val_mirnas <- unique(hub_mirna_val_df$miRNA[!is.na(hub_mirna_val_df$miRNA)])
cat("\nTotal unique validated miRNAs targeting ≥1 hub:", length(all_val_mirnas), "\n")

## Summary table: genes with vs without validated interactions
has_val <- sapply(hub_genes, function(g) {
  any(!is.na(hub_mirna_val_df$miRNA[hub_mirna_val_df$Hub_Gene == g]))
})
cat("\nHubs with ≥1 validated miRNA:\n")
print(has_val)

## ---- Save ----
write.csv(hub_mirna_val_df,
          file = file.path(out_dir, "F_hub_to_mirna_validated.csv"),
          row.names = FALSE)
write.csv(hub_mirna_mtb_df,
          file = file.path(out_dir, "F_hub_to_mirna_mirtarbase.csv"),
          row.names = FALSE)

cat("\nSaved F_hub_to_mirna_validated.csv\n")
cat("Saved F_hub_to_mirna_mirtarbase.csv\n")

## ---- Update provenance with multiMiR DB version ----
prov <- read.csv(file.path(out_dir, "F_multimir_provenance.csv"),
                 stringsAsFactors = FALSE)
prov$multiMiR_db_version <- "2.4.0"
prov$multiMiR_db_updated <- "2024-08-28"
prov$n_validated_mirnas_universe <- length(all_val_mirnas)
write.csv(prov, file = file.path(out_dir, "F_multimir_provenance.csv"), row.names = FALSE)
cat("\nProvenance updated.\n")

cat("\nStage 5F Step 1b complete.\n")
