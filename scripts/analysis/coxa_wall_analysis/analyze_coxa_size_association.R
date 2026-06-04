options(stringsAsFactors = FALSE)

root <- "<MANUSCRIPT_PROJECT_ROOT>"
input_file <- file.path(root, "05_Results", "05_Coxa_Thickness", "coxa_thickness_hole_merged_debug.csv")
key_file <- file.path(root, "03_Data_and_Inputs", "Analysis_Input", "specimen_key.csv")
out_dir <- file.path(root, "05_Results", "05_Coxa_Thickness", "Coxa_Size_Association")
supp_dir <- file.path(root, "06_Supplementary_Materials", "Supplementary_Tables_CSV", "S04_Coxa_Thickness")
code_dir <- file.path(root, "04_Analysis_Code", "Coxa_Thickness")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_dir, recursive = TRUE, showWarnings = FALSE)

read_csv2_safe <- function(path) {
  read.csv2(path, check.names = FALSE, na.strings = c("", "NA", "NaN"))
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- trimws(as.character(x))
  out <- rep(NA, length(x))
  out[toupper(x) == "TRUE"] <- TRUE
  out[toupper(x) == "FALSE"] <- FALSE
  out
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, formatC(p, format = "e", digits = 3), sprintf("%.4f", p)))
}

dat <- read_csv2_safe(input_file)
needed <- c("specimen", "bbox_diag_um", "median_thickness_um", "log_size", "log_thickness", "hole")
missing_cols <- setdiff(needed, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

dat$hole <- as_bool(dat$hole)
dat$hole_factor <- factor(dat$hole, levels = c(FALSE, TRUE), labels = c("Absent", "Present"))
dat$log_size <- log10(dat$bbox_diag_um)
dat$log_thickness <- log10(dat$median_thickness_um)

if (file.exists(key_file)) {
  key <- read_csv2_safe(key_file)
  if (all(c("specimen_key", "Family") %in% names(key))) {
    dat <- merge(dat, key[, c("specimen_key", "Family")], by.x = "specimen", by.y = "specimen_key", all.x = TRUE)
  }
}
if (!"Family" %in% names(dat)) dat$Family <- NA_character_

size_complete <- dat[is.finite(dat$bbox_diag_um) & is.finite(dat$median_thickness_um), ]
hole_complete <- dat[is.finite(dat$bbox_diag_um) & !is.na(dat$hole), ]

desc <- do.call(rbind, lapply(split(hole_complete, hole_complete$hole_factor), function(x) {
  data.frame(
    group = as.character(unique(x$hole_factor)),
    n = nrow(x),
    coxa_size_mean_um = mean(x$bbox_diag_um),
    coxa_size_sd_um = sd(x$bbox_diag_um),
    coxa_size_median_um = median(x$bbox_diag_um),
    coxa_size_min_um = min(x$bbox_diag_um),
    coxa_size_max_um = max(x$bbox_diag_um),
    wall_thickness_median_um = median(x$median_thickness_um),
    stringsAsFactors = FALSE
  )
}))

glm_size <- glm(hole ~ log_size, data = hole_complete, family = binomial())
glm_null <- glm(hole ~ 1, data = hole_complete, family = binomial())
glm_sum <- summary(glm_size)
coef_tab <- as.data.frame(coef(glm_sum))
coef_tab$term <- rownames(coef_tab)
rownames(coef_tab) <- NULL
names(coef_tab)[1:4] <- c("estimate_log_odds", "std_error", "z_value", "p_value")
coef_tab$odds_ratio <- exp(coef_tab$estimate_log_odds)
ci <- suppressMessages(confint.default(glm_size))
coef_tab$odds_ratio_low95 <- exp(ci[, 1])
coef_tab$odds_ratio_high95 <- exp(ci[, 2])
coef_tab <- coef_tab[, c("term", "estimate_log_odds", "std_error", "z_value", "p_value",
                         "odds_ratio", "odds_ratio_low95", "odds_ratio_high95")]

lrt <- anova(glm_null, glm_size, test = "Chisq")
model_tab <- data.frame(
  model = "hole_presence_vs_coxa_size",
  response = "Coxal wall hole presence",
  predictor = "log10(bbox_diag_um)",
  n = nrow(hole_complete),
  n_hole_absent = sum(!hole_complete$hole),
  n_hole_present = sum(hole_complete$hole),
  null_deviance = glm_null$deviance,
  residual_deviance = glm_size$deviance,
  lr_chisq = lrt$Deviance[2],
  lr_df = lrt$Df[2],
  lr_p_value = lrt$`Pr(>Chi)`[2],
  aic_null = AIC(glm_null),
  aic_model = AIC(glm_size),
  mcfadden_pseudo_r2 = 1 - as.numeric(logLik(glm_size) / logLik(glm_null)),
  stringsAsFactors = FALSE
)

wilcox_size <- wilcox.test(bbox_diag_um ~ hole_factor, data = hole_complete, exact = FALSE)
ttest_size <- t.test(log_size ~ hole_factor, data = hole_complete)
size_group_tests <- data.frame(
  test = c("Wilcoxon rank-sum", "Welch t-test"),
  response = c("bbox_diag_um", "log10(bbox_diag_um)"),
  group = "Coxal wall hole presence",
  statistic = c(unname(wilcox_size$statistic), unname(ttest_size$statistic)),
  df = c(NA_real_, unname(ttest_size$parameter)),
  p_value = c(wilcox_size$p.value, ttest_size$p.value),
  median_absent = c(median(hole_complete$bbox_diag_um[!hole_complete$hole]), median(hole_complete$log_size[!hole_complete$hole])),
  median_present = c(median(hole_complete$bbox_diag_um[hole_complete$hole]), median(hole_complete$log_size[hole_complete$hole])),
  mean_absent = c(mean(hole_complete$bbox_diag_um[!hole_complete$hole]), mean(hole_complete$log_size[!hole_complete$hole])),
  mean_present = c(mean(hole_complete$bbox_diag_um[hole_complete$hole]), mean(hole_complete$log_size[hole_complete$hole])),
  stringsAsFactors = FALSE
)

lm_thickness <- lm(log_thickness ~ log_size, data = size_complete)
lm_sum <- summary(lm_thickness)
lm_coef <- as.data.frame(coef(lm_sum))
lm_coef$term <- rownames(lm_coef)
rownames(lm_coef) <- NULL
names(lm_coef)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
lm_coef$model <- "wall_thickness_vs_coxa_size"
lm_coef <- lm_coef[, c("model", "term", "estimate", "std_error", "t_value", "p_value")]

lm_model <- data.frame(
  model = "wall_thickness_vs_coxa_size",
  response = "log10(median_thickness_um)",
  predictor = "log10(bbox_diag_um)",
  n = nrow(size_complete),
  r_squared = lm_sum$r.squared,
  adj_r_squared = lm_sum$adj.r.squared,
  f_statistic = unname(lm_sum$fstatistic[1]),
  df1 = unname(lm_sum$fstatistic[2]),
  df2 = unname(lm_sum$fstatistic[3]),
  model_p_value = pf(lm_sum$fstatistic[1], lm_sum$fstatistic[2], lm_sum$fstatistic[3], lower.tail = FALSE),
  residual_se = lm_sum$sigma,
  stringsAsFactors = FALSE
)

pearson <- cor.test(size_complete$bbox_diag_um, size_complete$median_thickness_um, method = "pearson")
spearman <- suppressWarnings(cor.test(size_complete$bbox_diag_um, size_complete$median_thickness_um, method = "spearman", exact = FALSE))
cor_tab <- data.frame(
  test = c("Pearson", "Spearman"),
  x = "bbox_diag_um",
  y = "median_thickness_um",
  n = nrow(size_complete),
  estimate = c(unname(pearson$estimate), unname(spearman$estimate)),
  statistic = c(unname(pearson$statistic), unname(spearman$statistic)),
  p_value = c(pearson$p.value, spearman$p.value),
  stringsAsFactors = FALSE
)

write.csv2(dat, file.path(out_dir, "coxa_size_association_dataset.csv"), row.names = FALSE)
write.csv2(desc, file.path(out_dir, "coxa_size_by_wall_hole_descriptives.csv"), row.names = FALSE)
write.csv2(coef_tab, file.path(out_dir, "coxa_size_wall_hole_glm_coefficients.csv"), row.names = FALSE)
write.csv2(model_tab, file.path(out_dir, "coxa_size_wall_hole_model_stats.csv"), row.names = FALSE)
write.csv2(size_group_tests, file.path(out_dir, "coxa_size_wall_hole_group_tests.csv"), row.names = FALSE)
write.csv2(lm_coef, file.path(out_dir, "coxa_thickness_size_lm_coefficients.csv"), row.names = FALSE)
write.csv2(lm_model, file.path(out_dir, "coxa_thickness_size_model_stats.csv"), row.names = FALSE)
write.csv2(cor_tab, file.path(out_dir, "coxa_thickness_size_correlations.csv"), row.names = FALSE)

# Supplement-friendly compact table.
supp_table <- rbind(
  data.frame(
    analysis = "Coxal wall hole presence vs coxa size",
    test = "binomial GLM likelihood-ratio test",
    response = "hole",
    predictor = "log10(bbox_diag_um)",
    n = model_tab$n,
    estimate = coef_tab$estimate_log_odds[coef_tab$term == "log_size"],
    statistic = model_tab$lr_chisq,
    df = model_tab$lr_df,
    p_value = model_tab$lr_p_value,
    effect_summary = paste0("odds ratio = ", signif(coef_tab$odds_ratio[coef_tab$term == "log_size"], 4)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    analysis = "Coxal wall thickness vs coxa size",
    test = "linear model",
    response = "log10(median_thickness_um)",
    predictor = "log10(bbox_diag_um)",
    n = lm_model$n,
    estimate = lm_coef$estimate[lm_coef$term == "log_size"],
    statistic = lm_model$f_statistic,
    df = lm_model$df2,
    p_value = lm_model$model_p_value,
    effect_summary = paste0("R2 = ", signif(lm_model$r_squared, 4)),
    stringsAsFactors = FALSE
  )
)
write.csv2(supp_table, file.path(out_dir, "coxa_size_association_summary_stats.csv"), row.names = FALSE)
write.csv2(supp_table, file.path(supp_dir, "coxa_size_association_summary_stats.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  theme_set(theme_classic(base_size = 11))

  p_size <- ggplot(hole_complete, aes(x = hole_factor, y = bbox_diag_um, fill = hole_factor)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65, color = "grey25") +
    geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8, color = "grey20") +
    scale_fill_manual(values = c("Absent" = "#8FB6C9", "Present" = "#D98C5F"), guide = "none") +
    labs(x = "Coxal wall hole", y = "Coxa size (bounding-box diagonal, um)") +
    annotate("text", x = 1.5, y = max(hole_complete$bbox_diag_um, na.rm = TRUE),
             label = paste0("Wilcoxon p = ", fmt_p(wilcox_size$p.value)),
             vjust = 1.2, size = 3.4)

  pred <- data.frame(log_size = seq(min(hole_complete$log_size), max(hole_complete$log_size), length.out = 200))
  pred$bbox_diag_um <- 10^pred$log_size
  pred$fit <- predict(glm_size, newdata = pred, type = "response")
  pred_link <- predict(glm_size, newdata = pred, type = "link", se.fit = TRUE)
  pred$low <- plogis(pred_link$fit - 1.96 * pred_link$se.fit)
  pred$high <- plogis(pred_link$fit + 1.96 * pred_link$se.fit)

  p_glm <- ggplot(hole_complete, aes(x = bbox_diag_um, y = as.numeric(hole))) +
    geom_ribbon(data = pred, aes(x = bbox_diag_um, ymin = low, ymax = high), inherit.aes = FALSE,
                fill = "#D98C5F", alpha = 0.2) +
    geom_line(data = pred, aes(x = bbox_diag_um, y = fit), inherit.aes = FALSE,
              color = "#A94F2B", linewidth = 0.9) +
    geom_point(aes(color = hole_factor), size = 2.1, alpha = 0.85,
               position = position_jitter(height = 0.035, width = 0)) +
    scale_color_manual(values = c("Absent" = "#4E7F99", "Present" = "#A94F2B"), guide = "none") +
    scale_y_continuous(breaks = c(0, 1), labels = c("Absent", "Present"), limits = c(-0.08, 1.08)) +
    labs(x = "Coxa size (bounding-box diagonal, um)", y = "Coxal wall hole") +
    annotate("text", x = min(hole_complete$bbox_diag_um), y = 1.06,
             hjust = 0, vjust = 1,
             label = paste0("GLM LR p = ", fmt_p(model_tab$lr_p_value)), size = 3.4)

  p_thick <- ggplot(size_complete, aes(x = bbox_diag_um, y = median_thickness_um)) +
    geom_point(aes(color = Family), size = 2, alpha = 0.8, show.legend = FALSE) +
    geom_smooth(method = "lm", se = TRUE, color = "grey15", fill = "grey70", linewidth = 0.8) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Coxa size (bounding-box diagonal, um)",
         y = "Median coxal wall thickness (um)") +
    annotate("text", x = min(size_complete$bbox_diag_um), y = max(size_complete$median_thickness_um),
             hjust = 0, vjust = 1,
             label = paste0("LM R2 = ", signif(lm_model$r_squared, 3),
                            ", p = ", fmt_p(lm_model$model_p_value)), size = 3.4)

  ggsave(file.path(out_dir, "coxa_size_by_wall_hole_boxplot.pdf"), p_size, width = 90, height = 80, units = "mm")
  ggsave(file.path(out_dir, "coxa_size_by_wall_hole_boxplot.png"), p_size, width = 90, height = 80, units = "mm", dpi = 450)
  ggsave(file.path(out_dir, "coxa_size_wall_hole_logistic.pdf"), p_glm, width = 95, height = 80, units = "mm")
  ggsave(file.path(out_dir, "coxa_size_wall_hole_logistic.png"), p_glm, width = 95, height = 80, units = "mm", dpi = 450)
  ggsave(file.path(out_dir, "coxa_thickness_vs_coxa_size.pdf"), p_thick, width = 95, height = 80, units = "mm")
  ggsave(file.path(out_dir, "coxa_thickness_vs_coxa_size.png"), p_thick, width = 95, height = 80, units = "mm", dpi = 450)
}

script_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (!is.na(script_file)) {
  file.copy(script_file, file.path(code_dir, "coxa_size_association_analysis.R"), overwrite = TRUE)
}

cat("Input rows:", nrow(dat), "\n")
cat("Hole model rows:", nrow(hole_complete), "\n")
cat("Hole present/absent:", sum(hole_complete$hole), "/", sum(!hole_complete$hole), "\n")
cat("GLM LR p:", model_tab$lr_p_value, "\n")
cat("Thickness-size LM R2:", lm_model$r_squared, "\n")
cat("Thickness-size LM p:", lm_model$model_p_value, "\n")
cat("Outputs:", out_dir, "\n")
