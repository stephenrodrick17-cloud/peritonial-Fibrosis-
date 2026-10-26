source(file.path("results", "paper_style", "00_helpers.R"))

meta <- read_project_csv(file.path(TABLES, "GSE125498_sample_metadata.csv"))
stopifnot(nrow(meta) == 33L, sum(meta$group == "LPD") == 13L, sum(meta$group == "SPD") == 20L)

steps <- data.frame(
  step = factor(c(
    "Discovery cohort\nGSE125498\n13 LPD / 20 SPD",
    "Differential expression\nlimma; LPD vs SPD",
    "Matrisome / WGCNA\nfeature pool",
    "Three selectors\nLASSO / SVM-RFE / RF",
    "Seven-gene consensus\npanel",
    "Context analyses\npathways / immune / network",
    "LPD-only subtyping\nexploratory"
  ), levels = c(
    "Discovery cohort\nGSE125498\n13 LPD / 20 SPD",
    "Differential expression\nlimma; LPD vs SPD",
    "Matrisome / WGCNA\nfeature pool",
    "Three selectors\nLASSO / SVM-RFE / RF",
    "Seven-gene consensus\npanel",
    "Context analyses\npathways / immune / network",
    "LPD-only subtyping\nexploratory"
  )),
  x = seq_len(7), y = 1
)
edges <- data.frame(x = 1:6, xend = 2:7, y = 1, yend = 1)
plot <- ggplot(steps, aes(x, y)) +
  geom_segment(data = edges, aes(x = x + 0.25, xend = xend - 0.25, y = y, yend = yend),
               inherit.aes = FALSE, colour = "#4A718F",
               arrow = grid::arrow(length = grid::unit(0.12, "inches"))) +
  geom_label(aes(label = as.character(step)), fill = "#EAF2F8", colour = "#16324F",
             label.size = 0.45, label.padding = grid::unit(0.18, "lines"), size = 3) +
  scale_x_continuous(limits = c(0.5, 7.5), breaks = NULL) +
  scale_y_continuous(limits = c(0.75, 1.25), breaks = NULL) +
  labs(title = "Study workflow", x = NULL, y = NULL) +
  theme_paper() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

source_rows <- data.frame(
  stage = seq_len(7),
  description = c("GSE125498 discovery", "limma differential expression",
                  "matrisome/WGCNA candidate pool", "LASSO, SVM-RFE, random forest",
                  "seven-gene majority-vote consensus", "existing context analyses",
                  "LPD-only exploratory clustering"),
  n_LPD = c(13, 13, 13, 13, 13, 13, 13),
  n_SPD = c(20, 20, 20, 20, 20, 20, NA),
  status = c("reused", "reused", "reused", "reused", "reused", "reused", "new")
)
save_csv(source_rows, "Fig1_source_data.csv")
save_plot_pair(plot, "Fig1_workflow", 15, 4)
write_status("Fig1", "complete", "Workflow schematic; cohort counts verified from saved sample metadata.")
