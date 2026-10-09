# ==============================================================================
# MASTER PIPELINE ORCHESTRATOR
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Discovery Pipeline
# Dataset: GSE62928 (Peritoneal Fibrosis Transcriptomics)
# Intersection: GSE62928 DEGs ∩ ECM genes all.xlsx (Naba et al. Matrisome)
# ==============================================================================

cat("==============================================================================\n")
cat("STARTING GSE62928 PERITONEAL FIBROSIS & ECM-DEG DISCOVERY PIPELINE\n")
cat("==============================================================================\n\n")

start_time <- Sys.time()

cat("\n>>> EXECUTING SCRIPT 01: LOAD GSE62928 & PREPROCESS <<<\n")
source("01_load_qc_preprocess.R")

cat("\n>>> EXECUTING SCRIPT 02: DIFFERENTIAL EXPRESSION <<<\n")
source("02_differential_expression.R")

cat("\n>>> EXECUTING SCRIPT 03: GSE62928 ∩ ECM MASTERLIST INTERSECTION <<<\n")
source("03_matrisome_filtering.R")

cat("\n>>> EXECUTING SCRIPT 04: FUNCTIONAL ENRICHMENT (GO/KEGG) <<<\n")
source("04_functional_enrichment.R")

end_time <- Sys.time()
elapsed  <- round(difftime(end_time, start_time, units = "mins"), 2)

cat("\n==============================================================================\n")
cat("PIPELINE COMPLETED SUCCESSFULLY IN", elapsed, "MINUTES!\n")
cat("==============================================================================\n")
