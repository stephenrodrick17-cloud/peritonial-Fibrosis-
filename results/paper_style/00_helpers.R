options(stringsAsFactors = FALSE)

ROOT <- normalizePath(getwd(), winslash = "\\", mustWork = TRUE)
OUT <- file.path(ROOT, "results", "paper_style")
TABLES <- file.path(ROOT, "results", "tables")
REV <- file.path(ROOT, "results", "revision")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

needed <- c("ggplot2", "dplyr", "tidyr")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

theme_paper <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E6E9ED", linewidth = 0.3),
      plot.title = element_text(face = "bold", colour = "#16324F"),
      axis.title = element_text(colour = "#16324F"),
      strip.text = element_text(face = "bold", colour = "#16324F"),
      legend.position = "bottom"
    )
}

new_path <- function(path) {
  if (file.exists(path)) stop("Refusing to overwrite existing output: ", path)
  invisible(path)
}

save_csv <- function(x, name) {
  path <- file.path(OUT, name)
  new_path(path)
  write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

save_plot_pair <- function(plot, stem, width = 9, height = 6) {
  png_path <- file.path(OUT, paste0(stem, ".png"))
  pdf_path <- file.path(OUT, paste0(stem, ".pdf"))
  new_path(png_path)
  new_path(pdf_path)
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  invisible(c(png_path, pdf_path))
}

read_project_csv <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

add_letter <- function(plot, letter) {
  plot + labs(tag = letter) +
    theme(plot.tag = element_text(face = "bold", size = 14, colour = "#16324F"))
}

write_status <- function(step, status, detail) {
  path <- file.path(OUT, "execution_status.csv")
  row <- data.frame(step = step, status = status, detail = detail)
  existing <- if (file.exists(path)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else data.frame(step = character(), status = character(), detail = character())
  existing <- existing[existing$step != step, , drop = FALSE]
  updated <- rbind(existing, row)
  updated <- updated[order(updated$step), , drop = FALSE]
  write.csv(updated, path, row.names = FALSE, na = "")
}
