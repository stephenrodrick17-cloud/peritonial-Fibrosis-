args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("jsonlite"), "T10")

hpa_rows <- list()
status <- character()
previous_timeout <- getOption("timeout")
options(timeout = 15)
for (gene in HUBS) {
  url <- paste0("https://www.proteinatlas.org/api/search_download.php?search=",
    utils::URLencode(gene, reserved = TRUE), "&format=json")
  tmp <- tempfile(fileext = ".json")
  ok <- tryCatch({
    utils::download.file(url, tmp, quiet = TRUE, mode = "wb", method = "libcurl")
    file.exists(tmp) && file.info(tmp)$size > 0
  }, error = function(e) {
    status <<- c(status, paste(gene, "HPA request failed:", conditionMessage(e)))
    FALSE
  })
  if (ok) {
    dat <- tryCatch(jsonlite::fromJSON(tmp, flatten = TRUE), error = function(e) NULL)
    if (!is.null(dat) && nrow(dat)) {
      dat$query_gene <- gene
      hpa_rows[[length(hpa_rows) + 1L]] <- dat
      status <- c(status, paste(gene, "retrieved", nrow(dat), "rows"))
    } else status <- c(status, paste(gene, "response empty or could not be parsed"))
  }
  unlink(tmp)
}
options(timeout = previous_timeout)
if (!length(hpa_rows)) {
  write.csv(data.frame(status = "skipped", reason = if (length(status)) paste(status, collapse = " | ") else "HPA returned no rows."),
            file.path(TABLES, "T10_HPA_status.csv"), row.names = FALSE)
  task_log("T10_external_resources", status, "skipped")
  append_results_text("T10. External-resource support", paste("Human Protein Atlas data could not be retrieved for the consensus panel during this run.",
    "No tissue or single-cell expression values were imputed."))
} else {
  dat <- do.call(rbind, hpa_rows)
  utils::write.csv(dat, file.path(TABLES, "T10_HPA_nine_hubs_raw.csv"), row.names = FALSE)
  task_log("T10_external_resources", c(status, "The HPA download is exploratory; no related single-cell atlas was pre-specified in local inputs."), "partial")
  append_results_text("T10. External-resource support", paste("Public Human Protein Atlas search data were retrieved for the genes with responses from the endpoint.",
    "No relevant pre-specified single-cell atlas was available locally, so cell-type-specific expression was not inferred."))
}
