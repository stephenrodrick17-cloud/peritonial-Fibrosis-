# ==============================================================================
# MASTER PIPELINE ORCHESTRATOR
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Discovery Pipeline
# Dataset: GSE125498 (Platform GPL10558 Illumina HumanHT-12 V4.0)
# ==============================================================================

cat("==============================================================================\n")
cat("STARTING GSE125498 PERITONEAL FIBROSIS & ECM-DEG DISCOVERY PIPELINE\n")
cat("==============================================================================\n\n")

start_time <- Sys.time()

cat("\n>>> EXECUTING SCRIPT 01: LOAD, QC & PREPROCESS <<<\n")
source("01_load_qc_preprocess.R")

cat("\n>>> EXECUTING SCRIPT 02: DIFFERENTIAL EXPRESSION (LIMMA) <<<\n")
source("02_differential_expression.R")

cat("\n>>> EXECUTING SCRIPT 03: ECM & MATRISOME FILTERING <<<\n")
source("03_matrisome_filtering.R")

cat("\n>>> EXECUTING SCRIPT 04: FUNCTIONAL ENRICHMENT (CLUSTERPROFILER) <<<\n")
source("04_functional_enrichment.R")

end_time <- Sys.time()
elapsed <- round(difftime(end_time, start_time, units = "mins"), 2)

cat("\n==============================================================================\n")
cat("PIPELINE COMPLETED SUCCESSFULLY IN", elapsed, "MINUTES!\n")
cat("==============================================================================\n")
