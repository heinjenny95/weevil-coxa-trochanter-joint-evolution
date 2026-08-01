options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("Missing value for --", key)
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("metrics", "key", "out-dir", "figure-dir")
missing_args <- setdiff(required, names(args))
if (length(missing_args) > 0) {
  stop("Missing argument(s): ", paste(paste0("--", missing_args), collapse = ", "))
}

metrics_file <- normalizePath(args[["metrics"]], winslash = "/", mustWork = TRUE)
key_file <- normalizePath(args[["key"]], winslash = "/", mustWork = TRUE)
out_dir <- args[["out-dir"]]
figure_dir <- args[["figure-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_table_robust <- function(path) {
  comma <- tryCatch(read.csv(path, check.names = FALSE), error = function(e) NULL)
  if (!is.null(comma) && ncol(comma) > 1) return(comma)
  semicolon <- tryCatch(read.csv2(path, check.names = FALSE), error = function(e) NULL)
  if (!is.null(semicolon) && ncol(semicolon) > 1) return(semicolon)
  stop("Could not parse table: ", path)
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  value <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("TRUE", "T", "1", "YES", "Y", "JA", "X")] <- TRUE
  out[value %in% c("FALSE", "F", "0", "NO", "N", "NEIN")] <- FALSE
  out
}

model_p <- function(model) {
  f <- summary(model)$fstatistic
  unname(pf(f[[1]], f[[2]], f[[3]], lower.tail = FALSE))
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

extract_coefficients <- function(model, model_name, metric, subset_name) {
  tab <- as.data.frame(summary(model)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  tab$model <- model_name
  tab$metric <- metric
  tab$subset <- subset_name
  tab[, c("subset", "metric", "model", "term", "estimate", "std_error", "statistic", "p_value")]
}

extract_lm_summary <- function(model, model_name, metric, subset_name) {
  s <- summary(model)
  f <- s$fstatistic
  data.frame(
    subset = subset_name,
    metric = metric,
    model = model_name,
    n = nobs(model),
    r_squared = s$r.squared,
    adjusted_r_squared = s$adj.r.squared,
    f_statistic = unname(f[[1]]),
    df1 = unname(f[[2]]),
    df2 = unname(f[[3]]),
    p_value = model_p(model),
    aic = AIC(model),
    stringsAsFactors = FALSE
  )
}

metrics <- read_table_robust(metrics_file)
key <- read_table_robust(key_file)

required_metrics <- c(
  "specimen", "bbox_diag_um", "median_3d_thickness_um", "p10_3d_thickness_um",
  "touches_original_stack_boundary", "analysis_scale"
)
missing_metrics <- setdiff(required_metrics, names(metrics))
if (length(missing_metrics) > 0) {
  stop("Metrics table is missing: ", paste(missing_metrics, collapse = ", "))
}

required_key <- c("specimen_key", "Family", "Coxal wall hole")
missing_key <- setdiff(required_key, names(key))
if (length(missing_key) > 0) {
  stop("Specimen key is missing: ", paste(missing_key, collapse = ", "))
}

key_small <- key[, required_key]
dat <- merge(metrics, key_small, by.x = "specimen", by.y = "specimen_key", all.x = TRUE)
dat$opening <- as_bool(dat$`Coxal wall hole`)
dat$opening_factor <- factor(dat$opening, levels = c(FALSE, TRUE), labels = c("Absent", "Present"))
dat$touches_original_stack_boundary <- as_bool(dat$touches_original_stack_boundary)
dat$log_size <- log10(dat$bbox_diag_um)
dat$log_median_thickness <- log10(dat$median_3d_thickness_um)
dat$log_p10_thickness <- log10(dat$p10_3d_thickness_um)

complete <- dat[
  is.finite(dat$log_size) & is.finite(dat$log_median_thickness) & !is.na(dat$opening),
]
if (nrow(complete) < 10) stop("Too few complete observations after merging.")

fit_suite <- function(df, response, metric_name, subset_name) {
  size_formula <- as.formula(paste(response, "~ log_size"))
  additive_formula <- as.formula(paste(response, "~ log_size + opening_factor"))
  interaction_formula <- as.formula(paste(response, "~ log_size * opening_factor"))

  size_model <- lm(size_formula, data = df)
  additive_model <- lm(additive_formula, data = df)
  interaction_model <- lm(interaction_formula, data = df)
  interaction_lrt <- anova(additive_model, interaction_model)

  size_only_opening <- glm(opening ~ log_size, data = df, family = binomial())
  size_thickness_opening <- glm(
    as.formula(paste("opening ~ log_size +", response)),
    data = df,
    family = binomial()
  )
  logistic_lrt <- anova(size_only_opening, size_thickness_opening, test = "Chisq")

  summaries <- rbind(
    extract_lm_summary(size_model, "thickness_vs_coxa_size", metric_name, subset_name),
    extract_lm_summary(additive_model, "thickness_vs_size_plus_opening", metric_name, subset_name),
    extract_lm_summary(interaction_model, "thickness_vs_size_by_opening", metric_name, subset_name)
  )
  coefficients <- rbind(
    extract_coefficients(size_model, "thickness_vs_coxa_size", metric_name, subset_name),
    extract_coefficients(additive_model, "thickness_vs_size_plus_opening", metric_name, subset_name),
    extract_coefficients(interaction_model, "thickness_vs_size_by_opening", metric_name, subset_name)
  )

  extra <- data.frame(
    subset = subset_name,
    metric = metric_name,
    analysis = c(
      "opening_by_size_wilcoxon",
      "opening_by_size_logistic",
      "size_by_opening_interaction",
      "opening_by_size_plus_thickness"
    ),
    statistic = c(
      unname(wilcox.test(bbox_diag_um ~ opening_factor, data = df, exact = FALSE)$statistic),
      unname(coef(summary(size_only_opening))["log_size", "z value"]),
      interaction_lrt$F[[2]],
      logistic_lrt$Deviance[[2]]
    ),
    df = c(
      NA_real_,
      NA_real_,
      interaction_lrt$Df[[2]],
      logistic_lrt$Df[[2]]
    ),
    p_value = c(
      wilcox.test(bbox_diag_um ~ opening_factor, data = df, exact = FALSE)$p.value,
      coef(summary(size_only_opening))["log_size", "Pr(>|z|)"],
      interaction_lrt$`Pr(>F)`[[2]],
      logistic_lrt$`Pr(>Chi)`[[2]]
    ),
    stringsAsFactors = FALSE
  )

  list(
    size_model = size_model,
    additive_model = additive_model,
    summaries = summaries,
    coefficients = coefficients,
    extra = extra
  )
}

primary <- fit_suite(complete, "log_median_thickness", "median_3d_thickness_um", "all_masks")
no_boundary <- complete[!complete$touches_original_stack_boundary, ]
boundary_sensitivity <- fit_suite(
  no_boundary,
  "log_median_thickness",
  "median_3d_thickness_um",
  "exclude_boundary_touching_masks"
)
p10_sensitivity <- fit_suite(complete, "log_p10_thickness", "p10_3d_thickness_um", "all_masks")

all_summaries <- rbind(primary$summaries, boundary_sensitivity$summaries, p10_sensitivity$summaries)
all_coefficients <- rbind(primary$coefficients, boundary_sensitivity$coefficients, p10_sensitivity$coefficients)
all_extra <- rbind(primary$extra, boundary_sensitivity$extra, p10_sensitivity$extra)

complete$relative_thickness <- 10^(
  complete$log_median_thickness - predict(primary$size_model, newdata = complete)
)

group_desc <- do.call(rbind, lapply(split(complete, complete$opening_factor), function(x) {
  data.frame(
    coxal_wall_opening = as.character(unique(x$opening_factor)),
    n = nrow(x),
    median_coxa_size_um = median(x$bbox_diag_um),
    median_whole_volume_thickness_um = median(x$median_3d_thickness_um),
    median_relative_thickness = median(x$relative_thickness),
    stringsAsFactors = FALSE
  )
}))

write.csv(complete, file.path(out_dir, "coxa_3d_wall_thickness_analysis_dataset.csv"), row.names = FALSE)
write.csv(all_summaries, file.path(out_dir, "coxa_3d_wall_thickness_model_stats.csv"), row.names = FALSE)
write.csv(all_coefficients, file.path(out_dir, "coxa_3d_wall_thickness_coefficients.csv"), row.names = FALSE)
write.csv(all_extra, file.path(out_dir, "coxa_3d_wall_thickness_sensitivity_tests.csv"), row.names = FALSE)
write.csv(group_desc, file.path(out_dir, "coxa_3d_wall_thickness_group_descriptives.csv"), row.names = FALSE)

primary_opening_p <- primary$coefficients$p_value[
  primary$coefficients$model == "thickness_vs_size_plus_opening" &
    primary$coefficients$term == "opening_factorPresent"
]
primary_size_p <- primary$extra$p_value[primary$extra$analysis == "opening_by_size_logistic"]
primary_size_r2 <- primary$summaries$r_squared[primary$summaries$model == "thickness_vs_coxa_size"]
primary_size_lm_p <- primary$summaries$p_value[primary$summaries$model == "thickness_vs_coxa_size"]

compact <- data.frame(
  analysis = c(
    "Whole-volume coxal wall thickness versus coxa size",
    "Size-corrected whole-volume coxal wall thickness versus coxal wall opening",
    "Coxal wall opening versus coxa size",
    "Opening-by-size interaction for whole-volume coxal wall thickness"
  ),
  test = c("linear model", "linear model", "binomial GLM", "nested-model F test"),
  n = nrow(complete),
  effect = c(
    primary$coefficients$estimate[
      primary$coefficients$model == "thickness_vs_coxa_size" & primary$coefficients$term == "log_size"
    ],
    primary$coefficients$estimate[
      primary$coefficients$model == "thickness_vs_size_plus_opening" &
        primary$coefficients$term == "opening_factorPresent"
    ],
    coef(glm(opening ~ log_size, data = complete, family = binomial()))[["log_size"]],
    NA_real_
  ),
  statistic = c(
    primary$summaries$f_statistic[primary$summaries$model == "thickness_vs_coxa_size"],
    primary$coefficients$statistic[
      primary$coefficients$model == "thickness_vs_size_plus_opening" &
        primary$coefficients$term == "opening_factorPresent"
    ],
    primary$extra$statistic[primary$extra$analysis == "opening_by_size_logistic"],
    primary$extra$statistic[primary$extra$analysis == "size_by_opening_interaction"]
  ),
  p_value = c(
    primary_size_lm_p,
    primary_opening_p,
    primary_size_p,
    primary$extra$p_value[primary$extra$analysis == "size_by_opening_interaction"]
  ),
  r_squared = c(primary_size_r2, primary$summaries$r_squared[
    primary$summaries$model == "thickness_vs_size_plus_opening"
  ], NA_real_, NA_real_),
  stringsAsFactors = FALSE
)
write.csv(compact, file.path(out_dir, "coxa_size_association_summary_stats.csv"), row.names = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
  stop("Packages ggplot2 and patchwork are required to create the figure.")
}

library(ggplot2)
library(patchwork)

opening_colors <- c("Absent" = "#8DB9C9", "Present" = "#E09A6B")
theme_nature <- theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 7, face = "bold", margin = margin(b = 1.5)),
    plot.subtitle = element_text(size = 6, margin = margin(b = 3)),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(4, "mm"),
    plot.margin = margin(2, 3, 2, 2)
  )

p1 <- ggplot(complete, aes(bbox_diag_um, median_3d_thickness_um, color = opening_factor)) +
  geom_smooth(
    data = complete,
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#222222",
    fill = "#BDBDBD",
    linewidth = 0.55,
    inherit.aes = TRUE
  ) +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = opening_colors, name = "Coxal wall opening") +
  labs(
    title = "Wall thickness vs coxa size",
    subtitle = paste0("LM R2 = ", sprintf("%.3f", primary_size_r2), ", P ", fmt_p(primary_size_lm_p)),
    x = "Coxa size\n(bounding-box diagonal, micrometres)",
    y = "Median 3D wall thickness (micrometres)"
  ) +
  theme_nature +
  theme(legend.position = "bottom", legend.direction = "horizontal")

p2 <- ggplot(complete, aes(opening_factor, relative_thickness, fill = opening_factor)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#777777", linewidth = 0.45) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, color = "#333333", linewidth = 0.5) +
  geom_jitter(width = 0.10, size = 1.25, alpha = 0.75, color = "#333333") +
  scale_fill_manual(values = opening_colors, guide = "none") +
  labs(
    title = "Size-corrected wall thickness",
    subtitle = paste0("Adjusted opening effect: P ", fmt_p(primary_opening_p)),
    x = "Coxal wall opening",
    y = "Relative wall thickness\n(observed / size-predicted)"
  ) +
  theme_nature

p3 <- ggplot(complete, aes(opening_factor, bbox_diag_um, fill = opening_factor)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, color = "#333333", linewidth = 0.5) +
  geom_jitter(width = 0.10, size = 1.25, alpha = 0.75, color = "#333333") +
  scale_fill_manual(values = opening_colors, guide = "none") +
  labs(
    title = "Coxa size",
    subtitle = paste0("Opening vs coxa size: P ", fmt_p(primary_size_p)),
    x = "Coxal wall opening",
    y = "Coxa size\n(bounding-box diagonal, micrometres)"
  ) +
  theme_nature

combined <- (p1 + p2 + p3) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 7))

figure_pdf <- file.path(figure_dir, "Supplementary_Fig_5_coxal_wall_thickness_size_opening_180mm.pdf")
figure_png <- file.path(figure_dir, "Supplementary_Fig_5_coxal_wall_thickness_size_opening_180mm.png")
ggsave(figure_pdf, combined, width = 180, height = 61, units = "mm", device = cairo_pdf)
ggsave(figure_png, combined, width = 180, height = 61, units = "mm", dpi = 600, bg = "white")

cat("Analysed specimens:", nrow(complete), "\n")
cat("Boundary-touching masks:", sum(complete$touches_original_stack_boundary), "\n")
cat("Whole-volume thickness vs size: R2 =", primary_size_r2, ", P =", primary_size_lm_p, "\n")
cat("Adjusted opening effect: P =", primary_opening_p, "\n")
cat("Opening vs size: P =", primary_size_p, "\n")
cat("Figure:", figure_pdf, "\n")

