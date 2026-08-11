#!/usr/bin/env Rscript

# Publication-style supplementary PGLS scatterplots for the robust helix fit.
# The regression line is an ordinary visual guide; inferential values are read
# from the corresponding PGLS results table and printed explicitly.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: plot_robust_supplementary_pgls.R TIP_DATA PGLS_RESULTS OUTPUT_DIR")
}

tip_file <- args[[1L]]
pgls_file <- args[[2L]]
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tips <- read.csv2(tip_file, check.names = FALSE, stringsAsFactors = FALSE)
pgls <- read.csv2(pgls_file, check.names = FALSE, stringsAsFactors = FALSE)
# UTF-8 BOMs are present in the released semicolon-delimited tables.
names(tips)[1L] <- "tree_label"
names(pgls)[1L] <- "term"

specs <- list(
  list(
    predictor = "abs_winding_angle_deg",
    title = "PC1 and fitted winding angle",
    x = "Absolute fitted winding angle (degrees)",
    filename = "Supplementary_Fig_16_robust_PGLS_angle.jpg"
  ),
  list(
    predictor = "axial_pitch_360",
    title = "PC1 and fitted axial pitch",
    x = "Fitted axial pitch per 360-degree turn",
    filename = "Supplementary_Fig_17_robust_PGLS_pitch.jpg"
  )
)

for (spec in specs) {
  model_row <- pgls[
    pgls$response == "PC1" &
      pgls$predictor == spec$predictor &
      pgls$term == spec$predictor,
    ,
    drop = FALSE
  ]
  if (nrow(model_row) != 1L) {
    stop("Expected exactly one PGLS slope row for PC1 ~ ", spec$predictor)
  }

  plot_df <- tips[
    is.finite(tips$PC1) & is.finite(tips[[spec$predictor]]),
    c("tree_tip", "PC1", spec$predictor),
    drop = FALSE
  ]
  subtitle <- sprintf(
    "Exploratory PGLS: slope = %.3g, P = %.3f, FDR-adjusted P = %.3f; n = %d proxy tips",
    model_row$estimate,
    model_row$p_value,
    model_row$fdr_p_value,
    model_row$n_taxa
  )

  p <- ggplot(plot_df, aes(x = .data[[spec$predictor]], y = PC1)) +
    geom_smooth(method = "lm", se = TRUE, colour = "#356AE6", fill = "grey75") +
    geom_point(size = 2.8, colour = "black") +
    geom_text_repel(aes(label = tree_tip), size = 3.2, max.overlaps = Inf, seed = 1) +
    labs(title = spec$title, subtitle = subtitle, x = spec$x, y = "PC1 score") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    filename = file.path(output_dir, spec$filename),
    plot = p,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
}
