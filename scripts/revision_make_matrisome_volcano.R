#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(MatrisomeAnalyzeR)
})

results_path <- file.path("results", "tables", "GSE125498_all_results.csv")
out_dir <- file.path("results", "revision", "matrisome_volcano")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(results_path)) {
  stop("Missing gene-level results table: ", results_path)
}
if (!requireNamespace("MatrisomeAnalyzeR", quietly = TRUE)) {
  stop("MatrisomeAnalyzeR is required. Install MatrisomeAnalyzeR v1.0.1 and retry.")
}

results <- read.csv(results_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("Gene", "logFC", "P.Value", "adj.P.Val", "DEG_Status")
if (!all(required_columns %in% names(results))) {
  stop("The GSE125498 results table is missing required gene/statistic columns.")
}

annotated <- MatrisomeAnalyzeR::matriannotate(
  data = results,
  gene.column = "Gene",
  species = "human"
)
matrisome <- annotated[
  annotated[["Annotated Matrisome Division"]] != "Non-matrisome",
  c(
    "Annotated Gene",
    "Annotated Matrisome Division",
    "Annotated Matrisome Category",
    "logFC",
    "P.Value",
    "adj.P.Val"
  ),
  drop = FALSE
]
names(matrisome)[1:3] <- c("Gene", "Matrisome_Division", "Matrisome_Category")
matrisome <- matrisome[!duplicated(matrisome$Gene), , drop = FALSE]
matrisome <- matrisome[is.finite(matrisome$logFC) & is.finite(matrisome$P.Value), , drop = FALSE]
matrisome$DEG_Status <- ifelse(
  matrisome$logFC >= 0.585 & matrisome$P.Value < 0.05,
  "Upregulated",
  ifelse(
    matrisome$logFC <= -0.585 & matrisome$P.Value < 0.05,
    "Downregulated",
    "Not significant"
  )
)
matrisome$NegLog10P <- -log10(pmax(matrisome$P.Value, .Machine$double.xmin))

up_count <- sum(matrisome$DEG_Status == "Upregulated")
down_count <- sum(matrisome$DEG_Status == "Downregulated")
if (nrow(matrisome) != 455L || up_count != 33L || down_count != 13L) {
  stop(
    "Matrisome volcano input did not reproduce the expected counts: ",
    nrow(matrisome), " genes; ", up_count, " up; ", down_count, " down."
  )
}

label_data <- do.call(rbind, lapply(c("Upregulated", "Downregulated"), function(status) {
  subset <- matrisome[matrisome$DEG_Status == status, , drop = FALSE]
  subset[order(subset$P.Value), , drop = FALSE][seq_len(min(6L, nrow(subset))), , drop = FALSE]
}))

volcano <- ggplot2::ggplot(
  matrisome,
  ggplot2::aes(x = logFC, y = NegLog10P, colour = DEG_Status)
) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  ggplot2::geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  ggplot2::geom_point(size = 2, alpha = 0.82) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = Gene),
    size = 3.2,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Downregulated" = "#3478B8",
      "Not significant" = "#AAB3BF",
      "Upregulated" = "#C9413A"
    ),
    breaks = c("Upregulated", "Not significant", "Downregulated")
  ) +
  ggplot2::labs(
    title = "Matrisome differential expression in GSE125498",
    subtitle = paste0(
      "MatrisomeAnalyzeR v1.0.1 annotation: 455 unique genes; ",
      up_count, " upregulated and ", down_count, " downregulated"
    ),
    x = expression(log[2] * " fold-change (long-term PD vs short-term PD)"),
    y = expression(-log[10] * " nominal P-value"),
    colour = NULL,
    caption = "Dashed thresholds: |log2FC| = 0.585 and nominal P = 0.05. P-values are not multiple-testing adjusted for this classification."
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(colour = "grey25"),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "top",
    plot.caption = ggplot2::element_text(hjust = 0, size = 8.5, colour = "grey30")
  )

ggplot2::ggsave(
  file.path(out_dir, "GSE125498_MatrisomeAnalyzeR_volcano.png"),
  plot = volcano,
  width = 8.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggplot2::ggsave(
  file.path(out_dir, "GSE125498_MatrisomeAnalyzeR_volcano.pdf"),
  plot = volcano,
  width = 8.5,
  height = 6.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)
write.csv(
  matrisome,
  file.path(out_dir, "GSE125498_matrisome_volcano_source_data.csv"),
  row.names = FALSE
)

cat("Unique MatrisomeAnalyzeR genes: ", nrow(matrisome), "\n", sep = "")
cat("Upregulated matrisome DEGs: ", up_count, "\n", sep = "")
cat("Downregulated matrisome DEGs: ", down_count, "\n", sep = "")
cat("Figures saved to: ", normalizePath(out_dir, winslash = "/"), "\n", sep = "")
