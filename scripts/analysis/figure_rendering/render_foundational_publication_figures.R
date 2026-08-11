suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- if (length(args) >= 2) args[[2]] else file.path(repo_root, "data", "foundational_figures")
manuscript_root <- if (length(args) >= 3 && nzchar(args[[3]])) normalizePath(args[[3]], winslash = "/", mustWork = TRUE) else NA_character_
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "scripts", "analysis", "figure_rendering", "publication_style.R"))
set.seed(20260811)

read_mixed <- function(path) {
  probe <- readLines(path, n = 2, warn = FALSE, encoding = "UTF-8")
  delim <- if (length(probe) && grepl(";", probe[[1]], fixed = TRUE)) ";" else ","
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
  if (x < 0.001) return("< 0.001")
  sprintf("%.3f", x)
}

save_plot_set <- function(plot, repo_name, canonical_base, width, height, dpi = 600) {
  repo_png <- file.path(out_dir, paste0(repo_name, ".png"))
  repo_pdf <- file.path(out_dir, paste0(repo_name, ".pdf"))
  ggsave(repo_png, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  ggsave(repo_pdf, plot = plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
  if (!is.na(manuscript_root)) {
    target <- file.path(manuscript_root, canonical_base)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    file.copy(repo_png, paste0(target, ".png"), overwrite = TRUE)
    file.copy(repo_pdf, paste0(target, ".pdf"), overwrite = TRUE)
  }
  invisible(c(repo_png, repo_pdf))
}

source_root <- file.path(repo_root, "data", "supplementary_source_data")

# Supplementary Figure 4: method-dependent clustering outcomes.
cluster_data <- read_mixed(file.path(source_root, "S01_PCA_and_Morphospace", "clusters_all_methods.csv")) %>%
  mutate(across(c(PC1, PC2), as.numeric)) %>%
  pivot_longer(
    cols = c(cluster_hclust, cluster_kmeans, cluster_mclust),
    names_to = "method",
    values_to = "cluster"
  ) %>%
  mutate(
    method = recode(method,
      cluster_hclust = "Hierarchical",
      cluster_kmeans = "k-means",
      cluster_mclust = "mclust"
    ),
    method = factor(method, levels = c("Hierarchical", "k-means", "mclust")),
    cluster = factor(paste("cluster", cluster), levels = names(cluster_cols))
  )

s4 <- ggplot(cluster_data, aes(PC1, PC2, fill = cluster)) +
  geom_point(shape = 21, size = 2.5, colour = "white", stroke = 0.48, alpha = 1) +
  facet_wrap(~method, nrow = 1) +
  scale_fill_manual(values = cluster_cols, drop = FALSE, name = NULL) +
  labs(x = "PC1", y = "PC2") +
  theme_pub(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(4.5, "mm"),
    legend.spacing.x = grid::unit(2.5, "mm"),
    panel.spacing.x = grid::unit(4.5, "mm"),
    strip.text = element_text(size = 9.5)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.7, colour = "white")))

save_plot_set(
  s4,
  "Supplementary_Fig_4_clustering_outcomes_alternative_methods",
  "04_Supplementary_Figures/Supplementary_Fig_4_clustering_outcomes_alternative_methods_180mm",
  7.09, 3.78
)

# Supplementary Figure 6: coxal wall thickness and opening state.
coxa <- read_mixed(file.path(source_root, "S04_Coxa_Thickness", "coxa_3d_wall_thickness_analysis_dataset.csv")) %>%
  mutate(
    opening_factor = factor(opening_factor, levels = c("Absent", "Present")),
    relative_thickness = as.numeric(relative_thickness)
  )
opening_cols <- c(Absent = "#168AA6", Present = "#EF7C35")

s6a <- ggplot(coxa, aes(bbox_diag_um, median_3d_thickness_um, fill = opening_factor)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, colour = teal, fill = teal_light, linewidth = 1, alpha = 0.48) +
  geom_point(shape = 21, size = 2.25, colour = "white", stroke = 0.45, alpha = 1) +
  scale_fill_manual(values = opening_cols, name = "Coxal wall opening") +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "Wall thickness vs coxa size",
    subtitle = "LM R2 = 0.534; P < 0.001",
    x = "Coxa size\n(bounding-box diagonal, micrometres)",
    y = "Median 3D wall thickness\n(micrometres)"
  ) + theme_pub(base_size = 6.5)

s6b <- ggplot(coxa, aes(opening_factor, relative_thickness, fill = opening_factor)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = muted_ink, linewidth = 0.7) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.88, colour = ink, linewidth = 0.6) +
  geom_point(position = position_jitter(width = 0.07, height = 0, seed = 20260811), shape = 21, size = 1.75, colour = "white", stroke = 0.38, alpha = 1) +
  scale_fill_manual(values = opening_cols, guide = "none") +
  labs(
    title = "Size-corrected wall thickness",
    subtitle = "Adjusted opening effect: P = 0.116",
    x = "Coxal wall opening",
    y = "Relative wall thickness\n(observed / size-predicted)"
  ) + theme_pub(base_size = 6.5)

s6c <- ggplot(coxa, aes(opening_factor, bbox_diag_um, fill = opening_factor)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.88, colour = ink, linewidth = 0.6) +
  geom_point(position = position_jitter(width = 0.07, height = 0, seed = 20260811), shape = 21, size = 1.75, colour = "white", stroke = 0.38, alpha = 1) +
  scale_fill_manual(values = opening_cols, guide = "none") +
  labs(
    title = "Coxa size",
    subtitle = "Opening vs coxa size: P = 0.803",
    x = "Coxal wall opening",
    y = "Coxa size\n(bounding-box diagonal, micrometres)"
  ) + theme_pub(base_size = 6.5)

s6 <- (s6a | s6b | s6c) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", widths = c(1.05, 1, 1)) &
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 6.8),
    axis.text = element_text(size = 6.2),
    plot.title = element_text(size = 7.2),
    plot.subtitle = element_text(size = 6.2),
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.3),
    plot.tag.position = c(0.035, 0.965),
    plot.margin = margin(11, 7, 8, 10)
  )

save_plot_set(
  s6,
  "Supplementary_Fig_6_coxal_wall_thickness_size_opening",
  "04_Supplementary_Figures/Supplementary_Fig_6_coxal_wall_thickness_size_opening_180mm",
  7.09, 3.15
)

# Supplementary Figure 8: trait contrasts across coxal-opening states.
opening_traits <- read_mixed(file.path(source_root, "S10_PCM_ASR_and_Disparity", "r28_coxal_wall_opening_plot_data.csv")) %>%
  mutate(
    coxal_wall_opening = factor(coxal_wall_opening, levels = c("absent", "present"), labels = c("Absent", "Present")),
    trait_label = factor(trait_label, levels = c(
      "PC1", "PC2", "Absolute winding angle (deg)", "Axial span", "Centroid size", "Absolute number of turns"
    ))
  )

s8 <- ggplot(opening_traits, aes(coxal_wall_opening, value, fill = coxal_wall_opening)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.88, colour = ink, linewidth = 0.6) +
  geom_point(position = position_jitter(width = 0.06, height = 0, seed = 20260811), shape = 21, size = 2, colour = "white", stroke = 0.42, alpha = 1) +
  facet_wrap(~trait_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = opening_cols, guide = "none") +
  labs(x = "Coxal wall opening", y = "Trait value") +
  theme_pub(base_size = 9) +
  theme(panel.spacing = grid::unit(4.5, "mm"), strip.text = element_text(size = 9))

save_plot_set(
  s8,
  "Supplementary_Fig_8_coxal_opening_trait_contrasts",
  "04_Supplementary_Figures/Supplementary_Fig_8_coxal_opening_trait_contrasts_180mm",
  7.09, 5.05
)

# Supplementary Figures 21-23: the same saturated predictor palettes are used
# for all response variables so colour encodes predictor consistently.
ecology <- read_mixed(file.path(source_root, "S06_Ecology_Matrix", "ecology_analysis_input_merged.csv"))
pgls <- read_mixed(file.path(source_root, "S07_Ecology_Tests", "ecology_pgls_factor_results.csv"))
phyanova <- read_mixed(file.path(source_root, "S07_Ecology_Tests", "ecology_phylogenetic_anova_results.csv"))

predictor_specs <- list(
  host_lineage_broad = list(title = "Host lineage", levels = c("angiosperm", "gymnosperm"), labels = c("Angiosperm", "Gymnosperm")),
  woody_association_broad = list(title = "Woody association", levels = c("nonwoody", "woody"), labels = c("Non-woody", "Woody")),
  larval_lifestyle_broad = list(title = "Larval lifestyle", levels = c("internal", "other"), labels = c("Internal", "Other")),
  fungal_association_broad = list(title = "Fungal association", levels = c("no", "yes"), labels = c("No", "Yes"))
)

ecology_panel <- function(response, predictor, ylab) {
  spec <- predictor_specs[[predictor]]
  dat <- ecology %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[predictor]])) %>%
    mutate(group = factor(.data[[predictor]], levels = spec$levels, labels = spec$labels)) %>%
    filter(!is.na(group))
  counts <- dat %>% count(group) %>% mutate(axis_label = paste0(group, "\nn = ", n))
  axis_labels <- setNames(counts$axis_label, counts$group)
  pgls_p <- pgls %>% filter(.data$response == .env$response, .data$predictor == .env$predictor) %>% summarise(value = min(p_value, na.rm = TRUE)) %>% pull(value)
  anova_p <- phyanova %>% filter(.data$response == .env$response, .data$predictor == .env$predictor) %>% slice(1) %>% pull(p_value)
  pal <- ecology_palettes[[predictor]]
  names(pal) <- spec$labels

  ggplot(dat, aes(group, .data[[response]], fill = group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.88, colour = ink, linewidth = 0.62) +
    geom_point(position = position_jitter(width = 0.065, height = 0, seed = 20260811), shape = 21, size = 2.15, colour = "white", stroke = 0.43, alpha = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = axis_labels) +
    labs(
      title = spec$title,
      subtitle = paste0("PGLS P ", fmt_p(pgls_p), "; phyANOVA P ", fmt_p(anova_p)),
      x = NULL,
      y = ylab
    ) +
    theme_pub(base_size = 9) +
    theme(plot.subtitle = element_text(size = rel(0.82), colour = muted_ink))
}

ecology_figure <- function(response, ylab) {
  p1 <- ecology_panel(response, "host_lineage_broad", ylab)
  p2 <- ecology_panel(response, "woody_association_broad", ylab)
  p3 <- ecology_panel(response, "larval_lifestyle_broad", ylab)
  p4 <- ecology_panel(response, "fungal_association_broad", ylab)
  ((p1 | p2) / (p3 | p4)) + plot_annotation(tag_levels = "a")
}

s21 <- ecology_figure("PC2", "PC2")
s22 <- ecology_figure("centroid_size", "Centroid size")
s23 <- ecology_figure("PC1", "PC1")

save_plot_set(s21, "Supplementary_Fig_21_ecology_PC2", "04_Supplementary_Figures/Supplementary_Fig_21_ecology_PC2_180mm", 7.09, 7.09)
save_plot_set(s22, "Supplementary_Fig_22_ecology_centroid_size", "04_Supplementary_Figures/Supplementary_Fig_22_ecology_centroid_size_180mm", 7.09, 7.09)
save_plot_set(s23, "Supplementary_Fig_23_ecology_PC1", "04_Supplementary_Figures/Supplementary_Fig_23_ecology_PC1_180mm", 7.09, 7.09)

message("Foundational publication-style figures written to: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
