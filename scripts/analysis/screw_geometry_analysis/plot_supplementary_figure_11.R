#!/usr/bin/env Rscript

# Rebuild the three-panel descriptive geometry figure from the corrected
# specimen-level shape-geometry table. Axial pitch is derived from the other
# two quantities and is labelled accordingly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_supplementary_figure_11.R <plot-data.csv> <specimen-key.csv> <output.png> <output.pdf>")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

read_robust <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (grepl("^sep\\s*=", lines[[1]], ignore.case = TRUE)) lines <- lines[-1]
  header <- lines[[1]]
  sep <- if (length(gregexpr(";", header, fixed = TRUE)[[1]]) >=
             length(gregexpr(",", header, fixed = TRUE)[[1]])) ";" else ","
  con <- textConnection(lines)
  on.exit(close(con), add = TRUE)
  read.table(con, header = TRUE, sep = sep, quote = "\"", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE, dec = ".")
}

to_num <- function(x) as.numeric(gsub(",", ".", as.character(x), fixed = TRUE))

plot_data <- read_robust(args[[1]])
specimens <- read_robust(args[[2]])

plot_data$angle_abs <- to_num(plot_data$angle_abs)
plot_data$axial_metric <- to_num(plot_data$axial_metric)
plot_data$axial_pitch <- to_num(plot_data$axial_pitch)

family_column <- if ("Family" %in% names(specimens)) "Family" else "family"
metadata <- specimens[, c("specimen_id", family_column)]
names(metadata) <- c("specimen_id", "family")
plot_data <- left_join(plot_data, metadata, by = "specimen_id")

family_order <- c("Anthribidae", "Attelabidae", "Belidae", "Brentidae",
                  "Caridae", "Curculionidae", "Nemonychidae")
family_colours <- c(
  Anthribidae = "#7F7F7F",
  Attelabidae = "#E76F51",
  Belidae = "#4C78A8",
  Brentidae = "#72B7B2",
  Caridae = "#E9C46A",
  Curculionidae = "#59A14F",
  Nemonychidae = "#B07AA1"
)
plot_data$family <- factor(plot_data$family, levels = family_order)

theme_joint <- theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 8),
    plot.tag = element_text(face = "bold", size = 8),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 6.5),
    legend.key.width = grid::unit(2.2, "mm"),
    legend.spacing.x = grid::unit(1.2, "mm"),
    legend.margin = margin(t = 1, unit = "mm")
  )

make_panel <- function(x, y, title, x_label, y_label) {
  ggplot(plot_data, aes(x = .data[[x]], y = .data[[y]], colour = family)) +
    geom_point(size = 1.8, alpha = 0.9) +
    geom_smooth(aes(group = 1), method = "lm", colour = "#4C9BB0",
                fill = "#9FD3DD", linewidth = 0.65, alpha = 0.28) +
    scale_colour_manual(
      values = family_colours,
      drop = TRUE,
      name = "Family",
      guide = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    labs(title = title, x = x_label, y = y_label) +
    theme_joint
}

p_a <- make_panel("angle_abs", "axial_metric", "Axial span vs winding angle",
                  "Absolute winding angle (deg)", "Axial span")
p_b <- make_panel("angle_abs", "axial_pitch", "Derived pitch vs winding angle",
                  "Absolute winding angle (deg)", "Derived axial pitch")
p_c <- make_panel("axial_metric", "axial_pitch", "Derived pitch vs axial span",
                  "Axial span", "Derived axial pitch")

figure <- (p_a + p_b + p_c) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

ggsave(args[[3]], figure, width = 180, height = 78, units = "mm", dpi = 400,
       bg = "white")
ggsave(args[[4]], figure, width = 180, height = 78, units = "mm", bg = "white")
