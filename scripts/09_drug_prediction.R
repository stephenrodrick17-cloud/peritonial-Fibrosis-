args <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", args[[1]])), "_common.R"))
require_pkgs(c("jsonlite", "ggplot2"), "T09")

out <- file.path(TABLES, "T09_Enrichr_drug_prediction.csv")
status_path <- file.path(TABLES, "T09_drug_prediction_status.csv")
curl <- Sys.which("curl.exe")
if (!nzchar(curl)) {
  write.csv(data.frame(status = "not_run", reason = "curl.exe unavailable for Enrichr POST requests."), status_path, row.names = FALSE)
  write_skip("T09_drug_prediction", "Enrichr's required list-upload POST could not be sent because curl.exe was unavailable.")
} else {
  hubs <- read_csv_safe(file.path(RESULTS, "revision", "ml_39_matrisome_wgcna", "ml_39_gene_results.csv"))
  selected <- hubs[hubs$Consensus_At_Least_2_of_3 == 1, ]
  up <- selected$Gene[selected$Discovery_log2FC > 0]
  down <- selected$Gene[selected$Discovery_log2FC < 0]
  libraries <- c("DSigDB", "LINCS_L1000_Chem_Pert_down", "LINCS_L1000_Chem_Pert_up", "DrugMatrix")
  rows <- list(); notes <- character()
  query_enrichr <- function(genes, direction, library) {
    if (!length(genes)) return(NULL)
    f <- tempfile(fileext = ".txt")
    writeLines(genes, f)
    response_file <- tempfile(fileext = ".json")
    uploaded <- tryCatch(system2(curl, c("--silent", "--show-error", "--fail", "--max-time", "15",
      "--form", paste0("list=@", f), "--form", paste0("description=", direction, "_consensus_hubs"),
      "--output", response_file,
      "https://maayanlab.cloud/Enrichr/addList"), stdout = TRUE, stderr = TRUE), warning = function(w) character())
    response_text <- if (file.exists(response_file)) readLines(response_file, warn = FALSE) else character()
    unlink(f)
    unlink(response_file)
    id_obj <- tryCatch(jsonlite::fromJSON(paste(response_text, collapse = "")), error = function(e) NULL)
    if (is.null(id_obj$userListId)) stop("Enrichr list upload failed or returned no userListId. ", paste(uploaded, collapse = " "))
    url <- paste0("https://maayanlab.cloud/Enrichr/enrich?userListId=", id_obj$userListId, "&backgroundType=", utils::URLencode(library, reserved = TRUE))
    txt <- system2(curl, c("--silent", "--show-error", "--fail", "--max-time", "15", url), stdout = TRUE, stderr = TRUE)
    response <- jsonlite::fromJSON(paste(txt, collapse = ""), simplifyVector = FALSE)
    result <- response[[library]]
    if (!length(result)) return(NULL)
    do.call(rbind, lapply(result, function(r) data.frame(Rank = r[[1]], Term = r[[2]], P = r[[3]],
      Z_score = r[[4]], Combined_score = r[[5]], Genes = paste(r[[6]], collapse = ";"),
      Adjusted_P = r[[7]], Library = library, Hub_direction = direction)))
  }
  for (lib in libraries) for (direction in c("Up", "Down")) {
    genes <- if (direction == "Up") up else down
    tryCatch({
      ans <- query_enrichr(genes, direction, lib)
      if (!is.null(ans)) rows[[length(rows) + 1L]] <<- ans
      notes <<- c(notes, paste(direction, lib, "queried for", length(genes), "genes"))
    }, error = function(e) notes <<- c(notes, paste("Query failed:", direction, lib, conditionMessage(e))))
  }
  results <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(results)) {
    results$Reversal_direction <- ifelse(results$Hub_direction == "Up",
      grepl("down", results$Library, ignore.case = TRUE), grepl("up", results$Library, ignore.case = TRUE))
    results$BH_FDR_across_queries <- p.adjust(results$P, method = "BH")
    write.csv(results, out, row.names = FALSE)
    reversal <- results[results$Reversal_direction, ]
    top <- head(reversal[order(reversal$BH_FDR_across_queries), ], 10)
    if (nrow(top)) {
      p <- ggplot2::ggplot(top, ggplot2::aes(reorder(Term, -log10(pmax(BH_FDR_across_queries, 1e-300))), -log10(pmax(BH_FDR_across_queries, 1e-300)))) +
        ggplot2::geom_col(fill = "#0072B2") + ggplot2::coord_flip() +
        ggplot2::labs(title = "Top direction-reversal drug signatures", x = NULL, y = "-log10 BH-adjusted P") + theme_project()
      save_plot(p, "T09_drug_reversal_signatures", 9, 6)
      edges <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
        genes <- strsplit(top$Genes[i], ";", fixed = TRUE)[[1]]
        data.frame(Gene = genes, Drug_signature = top$Term[i], Library = top$Library[i], Adjusted_P = top$BH_FDR_across_queries[i])
      }))
      write.csv(edges, file.path(TABLES, "T09_drug_gene_edges.csv"), row.names = FALSE)
    }
  } else {
    write.csv(data.frame(status = "no_results", reason = paste(notes, collapse = " | ")), status_path, row.names = FALSE)
  }
  task_log("T09_drug_prediction", c(notes, paste("Discovery-up consensus hubs:", paste(up, collapse = ", ")),
    paste("Discovery-down consensus hubs:", if (length(down)) paste(down, collapse = ", ") else "none; down-list queries skipped."),
    "Enrichr results are exploratory signatures, not evidence of drug efficacy; drug approval or anti-fibrotic status was not asserted without curation."), if (length(rows)) "partial" else "skipped")
  append_results_text("T9. Drug-signature exploration", if (length(rows)) paste("Enrichr perturbation libraries were queried separately using discovery-up and discovery-down consensus genes where those lists were non-empty.",
    "Candidate reversal signatures were prioritized by BH-adjusted P-values; these are database associations and do not establish therapeutic benefit.") else
    paste("Public Enrichr queries produced no usable results in this run; no drug candidates were inferred.",
      "Drug perturbation signatures would be exploratory and would not establish therapeutic efficacy."))
}
