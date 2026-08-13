#!/usr/bin/env Rscript

# Re-render all robust screw-geometry figures in the established manuscript
# design language. The script only consumes checked-in source data and is
# intentionally separate from the inferential workflows: it changes visual
# presentation, not analyses or numerical results.

suppressPackageStartupMessages({
  library(ape)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(phytools)
  library(readr)
  library(scales)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- if (length(args) >= 2) args[[2]] else file.path(repo_root, "data", "screw_geometry", "figures")
manuscript_root <- if (length(args) >= 3 && nzchar(args[[3]])) normalizePath(args[[3]], winslash = "/", mustWork = TRUE) else NA_character_
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_dir <- file.path(repo_root, "data", "screw_geometry", "figure_source_data")
screw_dir <- file.path(repo_root, "data", "screw_geometry")
tree_file <- file.path(repo_root, "data", "phylogeny", "P01_Trees", "01_primary_tree_calibrated_grafen.tre")
style_file <- file.path(repo_root, "scripts", "analysis", "figure_rendering", "publication_style.R")
source(style_file)
set.seed(20260811)

required <- c(
  file.path(source_dir, "Figure_5_geometry_schematic.png"),
  file.path(source_dir, "allometry_merged_table.csv"),
  file.path(source_dir, "phylogenetic_signal_plot_data.csv"),
  file.path(source_dir, "asr_standardized_tip_data.csv"),
  file.path(source_dir, "evolutionary_models_univariate.csv"),
  file.path(source_dir, "evolutionary_models_multivariate.csv"),
  file.path(source_dir, "pcm_tip_level_data.csv"),
  file.path(source_dir, "pgls_tree_variant_detail.csv"),
  file.path(source_dir, "pgls_leave_one_out_detail.csv"),
  file.path(source_dir, "ecology_tip_level_data.csv"),
  file.path(screw_dir, "shape_geometry_analysis_dataset.csv"),
  file.path(screw_dir, "robust_helix_metrics.csv"),
  tree_file
)
missing_files <- required[!file.exists(required)]
if (length(missing_files) > 0) stop("Missing required source files:\n", paste(missing_files, collapse = "\n"))

read_mixed <- function(path) {
  probe <- readLines(path, n = 2, warn = FALSE, encoding = "UTF-8")
  delim <- if (length(probe) > 0 && grepl(";", probe[[1]], fixed = TRUE)) ";" else ","
  decimal_mark <- if (delim == ";" && length(probe) > 1 && grepl("[0-9],[0-9]", probe[[2]])) "," else "."
  suppressMessages(read_delim(
    path,
    delim = delim,
    locale = locale(decimal_mark = decimal_mark),
    trim_ws = TRUE,
    show_col_types = FALSE,
    progress = FALSE
  ))
}

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) return(format(x, scientific = TRUE, digits = 2))
  sprintf("%.3f", x)
}

stat_line <- function(row, adjusted = TRUE) {
  if (is.null(row) || nrow(row) == 0) return("")
  r2 <- if ("r_squared" %in% names(row)) row$r_squared[[1]] else NA_real_
  p <- if ("p_value" %in% names(row)) row$p_value[[1]] else NA_real_
  padj <- if (adjusted && "p_value_adjusted" %in% names(row)) row$p_value_adjusted[[1]] else NA_real_
  pieces <- c(if (!is.na(r2)) sprintf("R2 = %.3f", r2), if (!is.na(p)) paste0("P = ", fmt_p(p)))
  if (!is.na(padj)) pieces <- c(pieces, paste0("Holm P = ", fmt_p(padj)))
  paste(pieces, collapse = "; ")
}

lm_stat_line <- function(data, response, predictor) {
  fit <- lm(reformulate(predictor, response), data = data)
  sm <- summary(fit)
  sprintf("Point model: R2 = %.3f; P = %s", sm$r.squared, fmt_p(coef(sm)[2, 4]))
}

add_stat <- function(p, label) {
  if (!nzchar(label)) return(p)
  p + annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.05, vjust = 1.3, size = 2.7, colour = "#607487")
}

lm_panel <- function(data, x, y, xlab, ylab, title = NULL, colour_family = TRUE, stats = "", log_y = FALSE, stats_outside = TRUE) {
  aes_points <- if (colour_family && "Family" %in% names(data)) aes(x = .data[[x]], y = .data[[y]], fill = Family) else aes(x = .data[[x]], y = .data[[y]])
  p <- ggplot(data, aes_points) +
    geom_smooth(aes(x = .data[[x]], y = .data[[y]]), method = "lm", se = TRUE, inherit.aes = FALSE, colour = teal, fill = teal_light, linewidth = 1.0, alpha = 0.48) +
    {
      if (colour_family && "Family" %in% names(data)) family_point(2.35)
      else geom_point(size = 2.15, alpha = 1, colour = ink)
    } +
    labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = if (stats_outside && nzchar(stats)) stats else NULL
    ) +
    theme_pub()
  if (colour_family && "Family" %in% names(data)) {
    present_families <- family_order[family_order %in% unique(as.character(data$Family))]
    p <- p + scale_fill_manual(
      values = family_cols,
      limits = family_order,
      breaks = present_families,
      drop = TRUE,
      name = "Family"
    )
  }
  if (log_y) p <- p + scale_y_log10()
  if (stats_outside) {
    return(
      p + theme(
        plot.subtitle = element_text(
          colour = "#607487",
          size = rel(0.82),
          hjust = 0,
          margin = margin(l = 18, b = 8)
        ),
        plot.margin = margin(12, 8, 9, 8)
      )
    )
  }
  add_stat(p, stats)
}

save_repo <- function(plot, filename, width, height, dpi = 420) {
  path <- file.path(out_dir, filename)
  ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(path)
}

save_canonical <- function(plot, rel_base, width, height, tiff = FALSE) {
  if (is.na(manuscript_root)) return(invisible(NULL))
  base <- file.path(manuscript_root, rel_base)
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(base, ".png"), plot = plot, width = width, height = height, units = "in", dpi = 500, bg = "white", limitsize = FALSE)
  ggsave(paste0(base, ".pdf"), plot = plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
  if (tiff) {
    ggsave(paste0(base, ".tif"), plot = plot, width = width, height = height, units = "in", dpi = 600, bg = "white", compression = "lzw", limitsize = FALSE)
  }
  invisible(base)
}

allom <- read_mixed(file.path(source_dir, "allometry_merged_table.csv"))
shape <- read_mixed(file.path(screw_dir, "shape_geometry_analysis_dataset.csv"))
shape <- shape %>%
  left_join(allom %>% select(specimen_id, Family, logCS), by = "specimen_id")
metrics <- read_mixed(file.path(screw_dir, "robust_helix_metrics.csv")) %>%
  left_join(allom %>% select(specimen_id, Family), by = "specimen_id")
uni_stats <- read_mixed(file.path(source_dir, "allometry_univariate_PC1_to_PC5_results.csv"))
geom_stats <- read_mixed(file.path(source_dir, "allometry_continuous_traits_results.csv"))
rrpp_stats <- read_mixed(file.path(source_dir, "allometry_rrpp_multivariate_results.csv"))

# Main Figure 5 ----------------------------------------------------------------
schematic_file <- file.path(source_dir, "Figure_5_geometry_schematic.png")
schematic_composite <- png::readPNG(schematic_file)
schematic_crop <- schematic_composite[
  , seq_len(round(dim(schematic_composite)[[2]] * 0.505)), , drop = FALSE
]
# Remove the panel letters embedded in the original composite. Both labels are
# added to the final two-panel layout with one shared style below.
schematic_crop[seq_len(round(dim(schematic_crop)[[1]] * 0.12)), , ] <- 1
# The original composite's panel-b y-axis title begins just to the right of the
# complete trochanter tip. Mask that isolated carry-over without narrowing the
# crop, which would remove the rounded distal end again.
label_rows <- seq(
  round(dim(schematic_crop)[[1]] * 0.38),
  round(dim(schematic_crop)[[1]] * 0.47)
)
label_cols <- seq(
  round(dim(schematic_crop)[[2]] * 0.965),
  dim(schematic_crop)[[2]]
)
schematic_crop[label_rows, label_cols, ] <- 1
figure5_tag_theme <- theme(
  plot.tag = element_text(
    family = "Arial",
    face = "bold",
    size = 12,
    colour = ink,
    hjust = 0,
    vjust = 1
  ),
  plot.tag.position = c(0.015, 0.985)
)
p5a <- ggplot() +
  annotation_custom(
    grid::rasterGrob(schematic_crop, interpolate = TRUE),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  labs(tag = "a") +
  theme_void() +
  figure5_tag_theme +
  theme(plot.margin = margin(7, 5, 4, 5))
p5b <- ggplot(shape, aes(PC1, PC2, fill = angle_abs, size = axial_pitch)) +
  geom_point(shape = 21, colour = "white", stroke = 0.45, alpha = 1) +
  scale_fill_viridis_c(option = "plasma", name = "Winding angle (degrees)") +
  scale_size_continuous(range = c(1.7, 5.0), guide = "none") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = grid::unit(3.0, "cm"), barheight = grid::unit(0.22, "cm"))) +
  labs(x = "PC1", y = "PC2", tag = "b") +
  theme_pub() +
  figure5_tag_theme +
  theme(legend.position = "bottom", plot.margin = margin(7, 5, 4, 5))
fig5 <- (p5a | p5b) +
  plot_layout(widths = c(1, 1))
save_repo(fig5, "Figure_5_robust_shape_geometry.png", 12.4, 6.9)
save_canonical(fig5, "02_Main_Figures/Figure_5_screw_geometry_and_morphospace_180mm", 7.09, 3.94, tiff = TRUE)

# Extended Data Figure 3 -------------------------------------------------------
shape_score <- read_mixed(file.path(source_dir, "allometry_full_shape_regression_scores.csv"))
rrpp_row <- rrpp_stats %>% filter(Term == "logCS") %>% slice(1)
rrpp_label <- if (nrow(rrpp_row) == 1) sprintf("Full shape: R2 = %.3f; F = %.2f; Z = %.2f; P = %s", rrpp_row$Rsq, rrpp_row$F, rrpp_row$Z, fmt_p(rrpp_row$`Pr(>F)`)) else ""

ed3a <- lm_panel(shape_score, "logCS", "full_shape_regression_score", "log centroid size", "Multivariate allometry score", stats = rrpp_label, stats_outside = TRUE)
ed3b <- lm_panel(allom, "logCS", "PC1", "log centroid size", "PC1", stats = stat_line(uni_stats %>% filter(trait == "PC1")), stats_outside = TRUE)
ed3c <- lm_panel(allom, "logCS", "PC2", "log centroid size", "PC2", stats = stat_line(uni_stats %>% filter(trait == "PC2")), stats_outside = TRUE)
ed3d <- lm_panel(allom %>% filter(!is.na(abs_winding_angle_deg)), "logCS", "abs_winding_angle_deg", "log centroid size", "Absolute winding angle (degrees)", stats = stat_line(geom_stats %>% filter(trait == "abs_winding_angle_deg")), stats_outside = TRUE)
ed3e <- lm_panel(allom %>% filter(!is.na(axial_span)), "logCS", "axial_span", "log centroid size", "Fitted axial span", stats = stat_line(geom_stats %>% filter(trait == "axial_span")), stats_outside = TRUE)
ed3_legend <- get_legend(
  ed3b +
    guides(fill = guide_legend(title = "Family", ncol = 2, byrow = TRUE, override.aes = list(shape = 21, colour = "white", size = 3))) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box.just = "center"
    )
)
ed3_legend_panel <- wrap_elements(full = ed3_legend)
ed3 <- wrap_plots(
  ed3a + theme(legend.position = "none"),
  ed3b + theme(legend.position = "none"),
  ed3c + theme(legend.position = "none"),
  ed3d + theme(legend.position = "none"),
  ed3e + theme(legend.position = "none"),
  ed3_legend_panel,
  ncol = 2
) + plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e", "")))
save_repo(ed3, "Extended_Data_Fig_3_robust_allometry.png", 12.8, 14.4)
save_canonical(ed3, "03_Extended_Data_Figures/Extended_Data_Figure_3_allometry_180mm", 7.09, 7.95)

# Extended Data Figure 4 -------------------------------------------------------
signal_df <- read_mixed(file.path(source_dir, "phylogenetic_signal_plot_data.csv")) %>%
  mutate(
    trait_label = recode(trait_label, "Absolute winding angle" = "Absolute winding angle", "Axial span" = "Fitted axial span"),
    trait_label = factor(trait_label, levels = rev(c("PC1", "PC2", "PC3", "PC4", "PC5", "Absolute winding angle", "Fitted axial span"))),
    method = factor(method, levels = c("Blomberg's K", "Pagel's lambda"))
  )
ed4 <- ggplot(signal_df, aes(estimate, trait_label, fill = method, alpha = signal_status)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = ink, linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = bluegrey, linewidth = 0.65) +
  geom_text(data = signal_df %>% filter(!is.na(sig_label), sig_label != ""), aes(label = sig_label), position = position_dodge(width = 0.72), hjust = -0.45, colour = ink, size = 3.2) +
  scale_fill_manual(values = c("Blomberg's K" = bluegrey, "Pagel's lambda" = coral)) +
  scale_alpha_manual(values = c("No signal detected" = 0.72, "Signal detected" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.09))) +
  labs(x = "Signal estimate", y = NULL, fill = NULL) +
  theme_pub(10) + theme(legend.position = "top")
save_repo(ed4, "Extended_Data_Fig_4_robust_phylogenetic_signal.jpg", 12.2, 6.1)
save_canonical(ed4, "03_Extended_Data_Figures/Extended_Data_Figure_4_phylogenetic_signal_180mm", 7.09, 3.7)

# Extended Data Figure 5 (base phytools rendering) -----------------------------
asr_input <- read_mixed(file.path(source_dir, "asr_standardized_tip_data.csv"))
asr_tree <- read.tree(tree_file)
asr_input <- asr_input[complete.cases(asr_input[, c("PC1", "PC2", "PC3", "PC4", "PC5", "abs_winding_angle_deg", "axial_span")]), ]
asr_tree <- drop.tip(asr_tree, setdiff(asr_tree$tip.label, asr_input$tree_label))
asr_input <- asr_input[match(asr_tree$tip.label, asr_input$tree_label), ]
asr_traits <- c("PC1", "PC2", "PC3", "PC4", "PC5", "abs_winding_angle_deg", "axial_span")
asr_titles <- c(PC1 = "PC1", PC2 = "PC2", PC3 = "PC3", PC4 = "PC4", PC5 = "PC5", abs_winding_angle_deg = "Winding angle", axial_span = "Fitted axial span")
asr_maps <- lapply(asr_traits, function(v) {
  x <- setNames(asr_input[[v]], asr_input$tree_label)
  obj <- contMap(asr_tree, x, plot = FALSE, lims = c(-2.2, 2.2), res = 120)
  obj$cols <- setNames(colorRampPalette(continuous_cols)(length(obj$cols)), names(obj$cols))
  obj
})
names(asr_maps) <- asr_traits

short_tip <- function(x) {
  fam <- sub("___.*$", "", x)
  genus <- sub("^.*___", "", x)
  paste0(genus, " (", fam, ")")
}

draw_ed5 <- function(device_path, kind = c("png", "jpeg", "pdf"), width = 12.8, height = 11.2, dpi = 420) {
  kind <- match.arg(kind)
  if (kind == "pdf") pdf(device_path, width = width, height = height, bg = "white", useDingbats = FALSE)
  if (kind == "png") png(device_path, width = width, height = height, units = "in", res = dpi, bg = "white", type = "cairo")
  if (kind == "jpeg") jpeg(device_path, width = width, height = height, units = "in", res = dpi, bg = "white", quality = 96)
  on.exit(dev.off(), add = TRUE)
  layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow = TRUE), widths = c(1.08,1,1), heights = c(1,1,1))
  par(mar = c(1.0, 0.4, 2.4, 0.4), oma = c(0.5, 0.2, 0.5, 0.2), family = "sans", fg = ink, col.axis = ink, col.lab = ink)
  par(mar = c(0.5, 0.2, 1.7, 0.2))
  plot(asr_tree, show.tip.label = TRUE, tip.color = "#000000", cex = 0.58, label.offset = 2.5, edge.color = "#A7B8C2", no.margin = FALSE)
  tiplabels(pch = 16, col = "#000000", cex = 0.62)
  title(main = "a   Reference phylogeny", adj = 0, line = 0.25, cex.main = 0.86, col.main = ink)
  for (i in seq_along(asr_maps)) {
    plot(asr_maps[[i]], legend = FALSE, ftype = "off", lwd = 3.0, outline = FALSE, mar = c(0.5,0.2,1.7,0.2), direction = "rightwards")
    title(main = paste0(letters[[i + 1]], "   ", asr_titles[[names(asr_maps)[[i]]]]), adj = 0, line = 0.25, cex.main = 0.86, col.main = ink)
  }
  plot.new()
  par(xpd = NA)
  legend_cols <- colorRampPalette(continuous_cols)(120)
  legend_x <- seq(0.18, 0.82, length.out = length(legend_cols) + 1)
  for (j in seq_along(legend_cols)) {
    rect(legend_x[[j]], 0.46, legend_x[[j + 1]], 0.51, col = legend_cols[[j]], border = NA)
  }
  text(0.50, 0.64, "Standardized trait score", cex = 0.72, col = ink)
  text(0.18, 0.36, "-2.2\nlow", cex = 0.66, col = ink)
  text(0.82, 0.36, "2.2\nhigh", cex = 0.66, col = ink)
}

draw_ed5(file.path(out_dir, "Extended_Data_Fig_5_robust_ASR.jpg"), "jpeg")
if (!is.na(manuscript_root)) {
  ed5base <- file.path(manuscript_root, "03_Extended_Data_Figures", "Extended_Data_Figure_5_ancestral_state_reconstructions_180mm")
  draw_ed5(paste0(ed5base, ".png"), "png", width = 7.09, height = 6.2, dpi = 500)
  draw_ed5(paste0(ed5base, ".pdf"), "pdf", width = 7.09, height = 6.2)
}

# Extended Data Figure 6 -------------------------------------------------------
uni_models <- read_mixed(file.path(source_dir, "evolutionary_models_univariate.csv")) %>%
  mutate(set_label = factor(set_label, levels = rev(c("PC1", "PC2", "PC3", "PC4", "PC5"))))
mv_models <- read_mixed(file.path(source_dir, "evolutionary_models_multivariate.csv")) %>%
  mutate(set_label = factor(set_label, levels = rev(c("PC1-PC5", "PC1-PC4"))))
model_cols <- c(BM = bluegrey, OU = ink, EB = coral)
model_shapes <- c(BM = 16, OU = 17, EB = 15)
ed6a <- ggplot(uni_models, aes(delta_plot, set_label, colour = model, shape = model, group = set_label)) +
  geom_segment(data = uni_models %>% filter(fit_status == "converged") %>% group_by(set_label) %>% summarise(xmin = min(delta_plot), xmax = max(delta_plot), .groups = "drop"), aes(x = xmin, xend = xmax, y = set_label, yend = set_label), inherit.aes = FALSE, colour = grid_col, linewidth = 1.5) +
  geom_point(size = 3.2) + geom_vline(xintercept = 2, linetype = "dashed", colour = bluegrey) +
  scale_colour_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
  labs(x = expression(Delta*AIC[c]), y = NULL, colour = NULL, shape = NULL) + theme_pub() + theme(legend.position = "top")
ed6b <- ggplot(mv_models, aes(delta_plot, set_label, colour = model, shape = model)) +
  geom_segment(data = mv_models %>% filter(fit_status == "converged") %>% group_by(set_label) %>% summarise(xmin = min(delta_plot), xmax = max(delta_plot), .groups = "drop"), aes(x = xmin, xend = xmax, y = set_label, yend = set_label), inherit.aes = FALSE, colour = grid_col, linewidth = 1.5) +
  geom_point(data = mv_models %>% filter(fit_status == "converged"), size = 3.2) +
  geom_point(data = mv_models %>% filter(fit_status == "failed"), aes(x = 11.8), shape = 4, colour = coral, size = 4, stroke = 1.1) +
  geom_text(data = mv_models %>% filter(fit_status == "failed"), aes(x = 12.3, label = "not reliable"), colour = coral, hjust = 0, size = 3) +
  geom_vline(xintercept = c(2, 10), linetype = "dashed", colour = bluegrey) +
  scale_colour_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
  scale_x_continuous(limits = c(-0.2, 20), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(Delta*AIC), y = NULL, colour = NULL, shape = NULL) + theme_pub() + theme(legend.position = "none")
ed6 <- (ed6a / ed6b) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.25, 0.85))
save_repo(ed6, "Extended_Data_Fig_6_robust_models.jpg", 12.2, 8.2)
save_canonical(ed6, "03_Extended_Data_Figures/Extended_Data_Figure_6_evolutionary_model_support_180mm", 7.09, 4.85)

# Supplementary Figure 10 ------------------------------------------------------
s10a <- ggplot(shape, aes(angle_abs)) + geom_histogram(bins = 14, fill = bluegrey, colour = ink, linewidth = 0.3) + labs(title = "a", x = "Absolute winding angle (degrees)", y = "Count") + theme_pub()
s10b <- lm_panel(shape, "fit_rms", "angle_abs", "Helix RMS", "Winding angle (degrees)", title = "b", colour_family = FALSE)
s10c <- ggplot(shape, aes(PC1, PC2, fill = fit_rms)) +
  geom_point(shape = 21, size = 2.4, colour = "white", stroke = 0.45) +
  scale_fill_gradient(low = "#55BED0", high = teal, name = "Helix RMS") +
  guides(fill = guide_colourbar(title.position = "top", barwidth = grid::unit(28, "mm"), barheight = grid::unit(3, "mm"))) +
  labs(title = "c", x = "PC1", y = "PC2") + theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical")
s10d <- ggplot(shape, aes(PC1, PC2, fill = angle_abs, size = fit_rms)) +
  geom_point(shape = 21, colour = "white", stroke = 0.45, alpha = 1) +
  scale_fill_gradientn(colours = continuous_cols, name = "Winding angle (degrees)") +
  scale_size_continuous(range = c(1.4, 4.5), name = "Helix RMS") +
  guides(
    fill = guide_colourbar(title.position = "top", barwidth = grid::unit(28, "mm"), barheight = grid::unit(3, "mm")),
    size = guide_legend(title.position = "top", nrow = 1)
  ) +
  labs(title = "d", x = "PC1", y = "PC2") + theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical")
s10e <- lm_panel(shape, "fit_rms_rel", "angle_abs", "Helix RMS / radius", "Winding angle (degrees)", title = "e", colour_family = FALSE)
s10f <- lm_panel(shape %>% filter(shape_regime == "main_region"), "PC1", "angle_abs", "PC1", "Winding angle (degrees)", title = "f", colour_family = FALSE)
s10 <- wrap_plots(s10a,s10b,s10c,s10d,s10e,s10f,ncol=2,guides="keep") &
  theme(
    legend.margin = margin(t = 3, b = 5),
    plot.margin = margin(9, 9, 10, 10)
  )
save_repo(s10, "Supplementary_Fig_10_robust_geometry_QC.png", 11.8, 13.2)
save_canonical(s10, "04_Supplementary_Figures/Supplementary_Fig_10_PC1_PC5_vs_winding_angle_180mm", 7.09, 8.3)

# Supplementary Figure 11 ------------------------------------------------------
axial_df <- allom %>% filter(!is.na(abs_winding_angle_deg), !is.na(axial_span), !is.na(axial_pitch_360))
s11a <- lm_panel(axial_df, "abs_winding_angle_deg", "axial_span", "Absolute winding angle (degrees)", "Fitted axial span (log10 scale)", title = "a", stats = lm_stat_line(axial_df, "axial_span", "abs_winding_angle_deg"), log_y = TRUE)
s11b <- lm_panel(axial_df, "abs_winding_angle_deg", "axial_pitch_360", "Absolute winding angle (degrees)", "Fitted axial pitch per\n360-degree turn (log10 scale)", title = "b", stats = lm_stat_line(axial_df, "axial_pitch_360", "abs_winding_angle_deg"), log_y = TRUE)
s11 <- (s11a | s11b) + plot_layout(guides="collect") & theme(legend.position="bottom")
save_repo(s11, "Supplementary_Fig_11_robust_axial_relationships.png", 11.8, 5.7)
save_canonical(s11, "04_Supplementary_Figures/Supplementary_Fig_11_axial_screw_geometry_relationships_180mm", 7.09, 3.6)

# Supplementary Figure 12 ------------------------------------------------------
tip_df <- read_mixed(file.path(source_dir, "pcm_tip_level_data.csv"))
pgls_span <- read_mixed(file.path(source_dir, "pgls_core_axial_span.csv"))
pgls_label <- function(response, predictor) {
  r <- pgls_span %>% filter(.data$response == .env$response, .data$predictor == .env$predictor) %>% slice(1)
  if (nrow(r) == 0) return("")
  sprintf("PGLS beta = %.3f; P = %s\nFDR P = %s; lambda = %.2f; n = %d", r$estimate, fmt_p(r$p_value), fmt_p(r$fdr_p_value), r$lambda, r$n_taxa)
}
s12a <- lm_panel(tip_df, "axial_span", "PC1", "Fitted axial span", "PC1", title = "a", stats = pgls_label("PC1", "axial_span"))
s12b <- lm_panel(tip_df, "axial_span", "PC2", "Fitted axial span", "PC2", title = "b", stats = pgls_label("PC2", "axial_span"))
s12 <- (s12a | s12b) + plot_layout(guides="collect") & theme(legend.position="bottom")
save_repo(s12, "Supplementary_Fig_12_robust_PGLS.jpg", 11.8, 5.7)
save_canonical(s12, "04_Supplementary_Figures/Supplementary_Fig_12_PGLS_shape_vs_axial_span_180mm", 7.09, 3.6)

# Supplementary Figure 13 ------------------------------------------------------
jt <- read_mixed(file.path(source_dir, "joint_type_plot_data.csv")) %>% filter(joint_type_strict %in% c("True screw-nut joint", "Unopposed screw configuration"))
jt_cols <- c("True screw-nut joint" = "#527F9A", "Unopposed screw configuration" = "#E47B3F")
box_panel <- function(y, label, panel_label) {
  ggplot(jt, aes(joint_type_strict, .data[[y]], fill = joint_type_strict)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.88, colour = ink, linewidth = 0.55) +
    geom_point(position = position_jitter(width = 0.12, height = 0, seed = 20260811), size = 1.55, alpha = 1, colour = ink) +
    scale_fill_manual(values = jt_cols) + labs(title = panel_label, x = NULL, y = label) + theme_pub() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1, size = 7.5))
}
s13 <- box_panel("abs_winding_angle_deg", "Absolute winding angle\n(degrees)", "a") |
  box_panel("axial_pitch_360", "Fitted axial pitch\nper 360-degree turn", "b") |
  box_panel("axial_span", "Fitted axial span", "c")
save_repo(s13, "Supplementary_Fig_13_robust_joint_type.png", 13.2, 4.8)
save_canonical(s13, "04_Supplementary_Figures/Supplementary_Fig_13_screw_geometry_by_joint_type_180mm", 7.09, 2.9)

# Supplementary Figure 14 ------------------------------------------------------
supp14_specs <- tribble(
  ~var, ~label, ~log_y,
  "PC3", "PC3", FALSE,
  "PC4", "PC4", FALSE,
  "PC5", "PC5", FALSE,
  "n_turns_abs", "Absolute number\nof turns", FALSE,
  "axial_pitch_360", "Fitted axial pitch\nper 360-degree turn", FALSE,
  "start_end_dist", "Start-to-end\ndistance", FALSE,
  "fit_radius", "Fitted radius", FALSE,
  "fit_rms", "Helix RMS", TRUE,
  "axial_span", "Fitted axial span", FALSE,
  "ratio_axial_span_fit_radius", "Axial span /\nfitted radius", FALSE,
  "ratio_start_end_fit_radius", "Start-to-end distance\n/ radius", FALSE,
  "ratio_fit_rms_fit_radius", "Helix RMS /\nfitted radius", TRUE
)
s14_plots <- lapply(seq_len(nrow(supp14_specs)), function(i) {
  spec <- supp14_specs[i,]
  st <- if (startsWith(spec$var, "PC")) uni_stats %>% filter(trait == spec$var) else geom_stats %>% filter(trait == spec$var)
  lm_panel(
    allom %>% filter(!is.na(.data[[spec$var]])),
    "logCS", spec$var, "log centroid size", spec$label,
    title = NULL,
    stats = sub("; Holm", "\nHolm", stat_line(st), fixed = TRUE),
    log_y = spec$log_y
  ) +
    labs(tag = letters[[i]]) +
    theme(
      legend.position = "none",
      plot.tag = element_text(
        colour = ink,
        face = "bold",
        size = rel(1.05),
        hjust = 0,
        vjust = 1
      ),
      plot.tag.position = c(0, 1),
      plot.subtitle = element_text(
        colour = "#607487",
        size = rel(0.72),
        hjust = 0,
        lineheight = 0.92,
        margin = margin(l = 18, b = 0)
      ),
      plot.margin = margin(8, 8, 7, 11)
    )
})
s14 <- wrap_plots(s14_plots, ncol = 3)
save_repo(s14, "Supplementary_Fig_14_robust_allometry.png", 13.2, 15.5)
save_canonical(s14, "04_Supplementary_Figures/Supplementary_Fig_14_additional_specimen_allometry_180mm", 7.09, 8.3)

# Supplementary Figures 16 and 17 ---------------------------------------------
pgls_all <- read_mixed(file.path(screw_dir, "sensitivity", "pgls_primary_adequate.csv"))
pgls_stat <- function(response, predictor) {
  r <- pgls_all %>% filter(.data$response == .env$response, .data$predictor == .env$predictor, term == .env$predictor) %>% slice(1)
  if (nrow(r) == 0) return("")
  sprintf("PGLS beta = %.3g; P = %s\nFDR P = %s; lambda = %.2f; n = %d", r$estimate, fmt_p(r$p_value), fmt_p(r$fdr_p_value), r$lambda, r$n_taxa)
}
s16 <- lm_panel(tip_df, "abs_winding_angle_deg", "PC1", "Absolute fitted winding angle (degrees)", "PC1", stats = pgls_stat("PC1","abs_winding_angle_deg")) + theme(legend.position="bottom")
s17 <- lm_panel(tip_df, "axial_pitch_360", "PC1", "Fitted axial pitch per 360-degree turn", "PC1", stats = pgls_stat("PC1","axial_pitch_360")) + theme(legend.position="bottom")
save_repo(s16, "Supplementary_Fig_16_robust_PGLS_angle.jpg", 7.0, 5.4)
save_repo(s17, "Supplementary_Fig_17_robust_PGLS_pitch.jpg", 7.0, 5.4)
save_canonical(s16, "04_Supplementary_Figures/Supplementary_Fig_16_additional_PGLS_traits_I_180mm", 7.09, 5.1)
save_canonical(s17, "04_Supplementary_Figures/Supplementary_Fig_17_additional_PGLS_traits_II_180mm", 7.09, 5.1)

# Supplementary Figure 18 ------------------------------------------------------
tree_detail <- read_mixed(file.path(source_dir, "pgls_tree_variant_detail.csv")) %>%
  filter(model_type == "tree_variant", response == "PC1", predictor %in% c("abs_winding_angle_deg", "axial_span"), !is.na(estimate)) %>%
  mutate(
    predictor_label = recode(predictor, abs_winding_angle_deg = "PC1 ~ winding angle", axial_span = "PC1 ~ fitted axial span"),
    tree_label = gsub("_", " ", tree_id),
    tree_class = recode(variant_group, working = "Working tree", existing_file = "Existing tree", generated_branch_lengths = "Generated branch lengths")
  )
tree_class_cols <- c("Working tree" = teal, "Existing tree" = bluegrey, "Generated branch lengths" = coral)
s18 <- ggplot(tree_detail, aes(estimate, tree_label, colour = tree_class)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = bluegrey) +
  geom_point(size = 2.45, alpha = 1) +
  facet_wrap(~predictor_label, scales = "free", ncol = 2) +
  scale_colour_manual(values = tree_class_cols, na.value = "#8B9AA5") +
  labs(x = "PGLS slope estimate", y = NULL, colour = "Tree class") + theme_pub() +
  theme(legend.position="bottom", axis.text.y = element_text(size = 6.8))
save_repo(s18, "Supplementary_Fig_18_robust_tree_PGLS.jpg", 12.8, 7.2)
save_canonical(s18, "04_Supplementary_Figures/Supplementary_Fig_18_tree_sensitivity_PGLS_slopes_180mm", 7.09, 4.4)

# Supplementary Figures 19 and 20 ---------------------------------------------
loo <- read_mixed(file.path(source_dir, "pgls_leave_one_out_detail.csv")) %>%
  filter(predictor == "axial_span", response %in% c("PC1", "PC2"), !is.na(estimate)) %>%
  mutate(
    dropped_label = gsub("_", " ", sub("^.*___", "", dropped_taxon)),
    response = factor(response, levels = c("PC1","PC2"))
  )
s19 <- ggplot(loo, aes(estimate, dropped_label, colour = response)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = bluegrey) + geom_point(size = 2.2) +
  facet_wrap(~response, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(PC1 = teal, PC2 = coral), guide = "none") +
  labs(x = "Slope after dropping one proxy tip", y = NULL) + theme_pub() + theme(axis.text.y = element_text(size=7))
s20 <- ggplot(loo, aes(lambda, dropped_label, colour = response)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = bluegrey) + geom_point(size = 2.2) +
  facet_wrap(~response, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(PC1 = teal, PC2 = coral), guide = "none") +
  labs(x = "Estimated Pagel's lambda after dropping one proxy tip", y = NULL) + theme_pub() + theme(axis.text.y = element_text(size=7))
save_repo(s19, "Supplementary_Fig_19_robust_LOO_slopes.jpg", 9.0, 9.5)
save_repo(s20, "Supplementary_Fig_20_robust_LOO_lambda.jpg", 9.0, 9.5)
save_canonical(s19, "04_Supplementary_Figures/Supplementary_Fig_19_leave_one_out_PGLS_slopes_180mm", 7.09, 7.1)
save_canonical(s20, "04_Supplementary_Figures/Supplementary_Fig_20_leave_one_out_Pagels_lambda_180mm", 7.09, 7.1)

# Supplementary Figures 24 and 25 ---------------------------------------------
eco <- read_mixed(file.path(source_dir, "ecology_tip_level_data.csv"))
eco_panel <- function(response, predictor, title, ylab, fill_cols) {
  dat <- eco %>% filter(!is.na(.data[[response]]), !is.na(.data[[predictor]]))
  ggplot(dat, aes(.data[[predictor]], .data[[response]], fill = .data[[predictor]])) +
    geom_boxplot(width = 0.62, alpha = 0.88, outlier.shape = NA, colour = ink, linewidth = 0.6) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 20260811), size = 2.05, colour = ink, alpha = 1) +
    scale_fill_manual(values = fill_cols) + labs(x = NULL, y = ylab, title = title) + theme_pub() +
    theme(legend.position="none", axis.text.x = element_text(angle = 12, hjust = 1))
}
eco_figure <- function(response, ylab) {
  p1 <- eco_panel(response,"woody_association_broad","Woody association",ylab,ecology_palettes$woody_association_broad)
  p2 <- eco_panel(response,"larval_lifestyle_broad","Larval lifestyle",ylab,ecology_palettes$larval_lifestyle_broad)
  p3 <- eco_panel(response,"fungal_association_broad","Fungal association",ylab,ecology_palettes$fungal_association_broad)
  (p1 | p2) / (p3 | plot_spacer()) + plot_annotation(tag_levels="a")
}
s24 <- eco_figure("abs_winding_angle_deg", "Absolute winding angle (degrees)")
s25 <- eco_figure("axial_span", "Fitted axial span")
save_repo(s24, "Supplementary_Fig_24_robust_ecology_angle.png", 11.8, 9.0)
save_repo(s25, "Supplementary_Fig_25_robust_ecology_span.png", 11.8, 9.0)
save_canonical(s24, "04_Supplementary_Figures/Supplementary_Fig_24_ecology_absolute_winding_angle_180mm", 7.09, 5.4)
save_canonical(s25, "04_Supplementary_Figures/Supplementary_Fig_25_ecology_axial_span_180mm", 7.09, 5.4)

# Supplementary Figure 28 ------------------------------------------------------
short_specimen <- function(x) {
  z <- sub("^[0-9]+_", "", x)
  z <- sub("_trochanter_aligned$", "", z)
  gsub("_", " ", z)
}
metrics <- metrics %>% mutate(
  label = short_specimen(specimen_id),
  quality_class = factor(
    recode(quality_class,
      good = "Good",
      caution = "Caution",
      limited_identifiability = "Limited identifiability"
    ),
    levels = c("Good", "Caution", "Limited identifiability")
  )
)
label_angle <- metrics %>% slice_max(abs(robust_minus_released_angle_deg), n = 5, with_ties = FALSE)
label_pitch <- metrics %>% filter(is.finite(released_endpoint_equivalent_pitch_360), is.finite(fitted_pitch_360), released_endpoint_equivalent_pitch_360 > 0, fitted_pitch_360 > 0) %>% slice_max(abs(log10(fitted_pitch_360 / released_endpoint_equivalent_pitch_360)), n = 5, with_ties = FALSE)
label_quality <- metrics %>% slice_min(axial_angle_r_squared, n = 4, with_ties = FALSE)
label_rms <- metrics %>% slice_max(helix_rms_relative_to_radius, n = 4, with_ties = FALSE)
qcols <- c("Good" = teal, "Caution" = gold, "Limited identifiability" = coral)
s28a <- ggplot(metrics, aes(released_abs_winding_angle_deg, abs_winding_angle_deg, colour = quality_class)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour=bluegrey) + geom_point(size=2) +
  geom_text_repel(data=label_angle,aes(label=label),size=2.2,max.overlaps=Inf,box.padding=.25,seed=20260811) +
  scale_colour_manual(values=qcols, name="Fit quality", drop=FALSE) + labs(title="a", x="Released winding angle (degrees)",y="Robust fitted angle (degrees)") + theme_pub() + theme(legend.position="none")
s28b <- ggplot(metrics, aes(released_endpoint_equivalent_pitch_360, fitted_pitch_360, colour = quality_class)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour=bluegrey) + geom_point(size=2) +
  geom_text_repel(data=label_pitch,aes(label=label),size=2.2,max.overlaps=Inf,box.padding=.25,seed=20260811) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values=qcols, name="Fit quality", drop=FALSE) +
  labs(title="b", x="Released endpoint-equivalent pitch",y="Robust fitted pitch") + theme_pub() + theme(legend.position="none")
s28c <- ggplot(metrics, aes(abs_winding_angle_deg, axial_angle_r_squared, colour = quality_class)) + geom_point(size=2) +
  geom_text_repel(data=label_quality,aes(label=label),size=2.2,max.overlaps=Inf,box.padding=.25,seed=20260811) +
  scale_colour_manual(values=qcols, name="Fit quality", drop=FALSE) + labs(title="c", x="Robust fitted angle (degrees)",y="Axial-angular R2") + theme_pub() + theme(legend.position="none")
s28d <- ggplot(metrics, aes(abs_winding_angle_deg, helix_rms_relative_to_radius, colour = quality_class)) + geom_point(size=2) +
  geom_hline(yintercept=.10,linetype="dashed",colour=coral) +
  geom_text_repel(data=label_rms,aes(label=label),size=2.2,max.overlaps=Inf,box.padding=.25,seed=20260811) +
  scale_colour_manual(values=qcols, name="Fit quality", drop=FALSE) + labs(title="d", x="Robust fitted angle (degrees)",y="Helix RMS / fitted radius") + theme_pub() + theme(legend.position="none")
s28_legend_plot <- ggplot(
  data.frame(
    x = seq_along(qcols),
    y = 1,
    quality_class = factor(names(qcols), levels = names(qcols))
  ),
  aes(x, y, colour = quality_class)
) +
  geom_point(size = 2.7) +
  scale_colour_manual(values = qcols, name = "Fit quality", drop = FALSE) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", colour = ink, size = 9),
    legend.text = element_text(colour = ink, size = 8.5),
    legend.margin = margin(0, 0, 0, 0)
  )
s28_legend <- cowplot::get_legend(s28_legend_plot)
s28 <- ((s28a | s28b) / (s28c | s28d) / wrap_elements(full = s28_legend)) +
  plot_layout(heights = c(1, 1, 0.075))
save_repo(s28, "Supplementary_Fig_28_robust_fit_audit.png", 12.2, 10.4)
save_canonical(s28, "04_Supplementary_Figures/Supplementary_Fig_28_screw_geometry_quality_control_180mm", 7.09, 6.0)

message("Publication-style robust figures written to: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
