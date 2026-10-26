args <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[[1]]), winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "_common.R"))
set.seed(2026)

tasks <- sprintf("%02d", 1:11)
for (id in tasks) {
  files <- list.files(script_dir, pattern = paste0("^", id, "_.*\\.R$"), full.names = TRUE)
  if (length(files) != 1) {
    message("Skipping task ", id, ": expected exactly one task script, found ", length(files))
    next
  }
  message("Running ", basename(files))
  tryCatch(
    sys.source(files, envir = .GlobalEnv),
    error = function(e) {
      detail <- paste("Runner caught task error:", conditionMessage(e))
      task_log(paste0("T", id, "_runner_error"), detail, "failed")
      append_results_text(paste0("T", id, " task error"), paste(detail, "See the matching file in results/logs/."))
      message(detail)
    }
  )
}
writeLines(capture.output(sessionInfo()), file.path(LOGS, "sessionInfo.txt"))
message("Finished task runner. Review results/logs/ for completed, partial, skipped, or failed task statuses.")
