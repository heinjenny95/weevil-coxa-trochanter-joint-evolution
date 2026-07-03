# ============================================================
# Coxa thickness (median) - exploratory plots (no hole variable)
# Combined 2x2 figure + CSV export of model statistics
# ============================================================

library(tidyverse)
library(ggrepel)
library(patchwork)

# ---- paths ----
in_csv    <- "<BEETLE_JOINTS_ROOT>/Processed/Coxa/coxa_combined_metrics.csv"
out_dir   <- "<BEETLE_JOINTS_ROOT>/Processed/Coxa"
out_png   <- file.path(out_dir, "coxa_thickness_2x2.png")
out_pdf   <- file.path(out_dir, "coxa_thickness_2x2.pdf")
stats_csv <- file.path(out_dir, "coxa_thickness_lm_stats.csv")

# ---- read ----
df <- read.csv2(in_csv, stringsAsFactors = FALSE)  # sep=";" dec=","

# ---- basic sanity ----
needed <- c("specimen","bbox_diag_um","median_thickness_um","first_slice","last_slice","mid_slice")
stopifnot(all(needed %in% names(df)))

df <- df %>%
  mutate(
    specimen = as.character(specimen),
    bbox_diag_um = as.numeric(bbox_diag_um),
    median_thickness_um = as.numeric(median_thickness_um),
    first_slice = as.numeric(first_slice),
    last_slice = as.numeric(last_slice),
    mid_slice = as.numeric(mid_slice),
    span_slices = last_slice - first_slice,
    mid_rel = (mid_slice - first_slice) / pmax(span_slices, 1)  # ideally around 0.5
  ) %>%
  filter(is.finite(bbox_diag_um), is.finite(median_thickness_um))

# ---- quick summary in console ----
cat("\nSummary: bbox_diag_um\n")
print(summary(df$bbox_diag_um))

cat("\nSummary: median_thickness_um\n")
print(summary(df$median_thickness_um))

# ---- helper for outlier labels (top/bottom N by thickness) ----
N_lab <- 6
label_df <- bind_rows(
  df %>% arrange(desc(median_thickness_um)) %>% slice_head(n = N_lab),
  df %>% arrange(median_thickness_um) %>% slice_head(n = N_lab)
) %>%
  distinct(specimen, .keep_all = TRUE)

# ============================================================
# PLOTS
# ============================================================

# 1) log-log thickness vs size + linear fit
p1 <- ggplot(df, aes(x = bbox_diag_um, y = median_thickness_um)) +
  geom_point(alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Median cuticle thickness vs Coxa size",
    x = "Coxa size (bbox diagonal, m) [log10]",
    y = "Median cuticle thickness (m) [log10]"
  ) +
  theme_classic()

# 2) thickness distribution
p2 <- ggplot(df, aes(x = median_thickness_um)) +
  geom_histogram(bins = 25) +
  labs(
    title = "Distribution of median cuticle thickness",
    x = "Median cuticle thickness (m)",
    y = "Count"
  ) +
  theme_classic()

# 3) same as (1) but label extremes
p3 <- ggplot(df, aes(x = bbox_diag_um, y = median_thickness_um)) +
  geom_point(alpha = 0.7) +
  geom_point(data = label_df, size = 2) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = specimen),
    max.overlaps = Inf,
    size = 3
  ) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Thickness vs size (labels = extreme thickness values)",
    x = "Coxa size (bbox diagonal, m) [log10]",
    y = "Median cuticle thickness (m) [log10]"
  ) +
  theme_classic()

# 4) slice QC: is mid_slice really mid between first/last?
p4 <- ggplot(df, aes(x = span_slices, y = mid_rel)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.8) +
  labs(
    title = "Mid-slice QC",
    x = "Coxa extent in Z (last_slice - first_slice)",
    y = "Relative mid position (0..1; ideally ~0.5)"
  ) +
  theme_classic()

# ============================================================
# COMBINED PLOT (2x2 GRID)
# ============================================================

combined_plot <- (p1 + p2) / (p3 + p4)

# ============================================================
# EXPORT FIGURE
# ============================================================

ggsave(
  filename = out_png,
  plot = combined_plot,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  filename = out_pdf,
  plot = combined_plot,
  width = 12,
  height = 10
)

cat("\nWrote figure PNG:", out_png, "\n")
cat("Wrote figure PDF:", out_pdf, "\n")

# ============================================================
# LINEAR MODEL
# ============================================================

model <- lm(log10(median_thickness_um) ~ log10(bbox_diag_um), data = df)
model_sum <- summary(model)

cat("\nLinear model summary:\n")
print(model_sum)

# coefficient table
coef_tab <- as.data.frame(model_sum$coefficients)
coef_tab$term <- rownames(coef_tab)
rownames(coef_tab) <- NULL

coef_tab <- coef_tab %>%
  select(term, everything()) %>%
  rename(
    estimate = Estimate,
    std_error = `Std. Error`,
    t_value = `t value`,
    p_value = `Pr(>|t|)`
  ) %>%
  mutate(type = "coefficient")

# model-level statistics
fstat <- unname(model_sum$fstatistic)

stats_tab <- data.frame(
  type = "model",
  term = "overall_model",
  estimate = NA_real_,
  std_error = NA_real_,
  t_value = NA_real_,
  p_value = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE),
  r_squared = model_sum$r.squared,
  adj_r_squared = model_sum$adj.r.squared,
  f_statistic = fstat[1],
  df1 = fstat[2],
  df2 = fstat[3],
  residual_se = model_sum$sigma,
  n = nrow(df)
)

# add empty model-stat columns to coef table so bind_rows works cleanly
coef_tab <- coef_tab %>%
  mutate(
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    f_statistic = NA_real_,
    df1 = NA_real_,
    df2 = NA_real_,
    residual_se = NA_real_,
    n = NA_real_
  )

# combine and export
out_stats <- bind_rows(coef_tab, stats_tab)

write.csv2(out_stats, stats_csv, row.names = FALSE)

cat("Wrote stats CSV:", stats_csv, "\n")


################################################################################

# ============================================================
# Coxa thickness vs Coxal wall hole
# sinnvoll reduzierte Analyse:
# - kein Wilcoxon
# - keine inferenzstatistische Analyse fuer Einbuchtung
# - Fokus auf size-corrected effect of hole
# - zusaetzlich: Interaktionstest und Residualplots
# ============================================================

library(tidyverse)
library(patchwork)

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
metrics_csv <- "<BEETLE_JOINTS_ROOT>/Processed/Coxa/coxa_combined_metrics.csv"
key_csv     <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_key.csv"
out_dir     <- "<BEETLE_JOINTS_ROOT>/Processed/Coxa"

out_stats_csv   <- file.path(out_dir, "coxa_thickness_hole_model_stats.csv")
out_counts_csv  <- file.path(out_dir, "coxa_trait_group_counts.csv")
out_merge_csv   <- file.path(out_dir, "coxa_thickness_hole_merged_debug.csv")
out_plot_png    <- file.path(out_dir, "coxa_thickness_hole_2x2.png")
out_plot_pdf    <- file.path(out_dir, "coxa_thickness_hole_2x2.pdf")

# ------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------
clean_names_basic <- function(x) {
  x %>%
    trimws() %>%
    gsub("\\s+", " ", .)
}

to_logical_trait <- function(x) {
  if (is.logical(x)) return(x)
  
  x <- as.character(x)
  x <- trimws(x)
  x_low <- tolower(x)
  
  case_when(
    x_low %in% c("true", "t", "1", "yes", "y", "ja", "x") ~ TRUE,
    x_low %in% c("false", "f", "0", "no", "n", "nein", "") ~ FALSE,
    is.na(x) ~ NA,
    TRUE ~ NA
  )
}

read_table_robust <- function(path) {
  x <- tryCatch(
    read.csv2(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(x) && ncol(x) > 1) return(x)
  
  x <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(x) && ncol(x) > 1) return(x)
  
  stop(paste("Konnte Datei nicht sinnvoll einlesen:", path))
}

extract_coefs <- function(model, model_name) {
  s <- summary(model)
  
  coef_tab <- as.data.frame(s$coefficients)
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL
  
  coef_tab %>%
    rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      statistic = `t value`,
      p_value = `Pr(>|t|)`
    ) %>%
    mutate(
      model = model_name,
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      aic = NA_real_,
      residual_se = NA_real_,
      df_model = NA_real_,
      df_resid = NA_real_,
      comparison = NA_character_,
      comparison_p = NA_real_
    ) %>%
    select(
      model, term, estimate, std_error, statistic, p_value,
      r_squared, adj_r_squared, aic, residual_se,
      df_model, df_resid, comparison, comparison_p
    )
}

extract_model_summary <- function(model, model_name) {
  s <- summary(model)
  
  data.frame(
    model = model_name,
    term = "overall_model",
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = pf(
      s$fstatistic[1],
      s$fstatistic[2],
      s$fstatistic[3],
      lower.tail = FALSE
    ),
    r_squared = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    aic = AIC(model),
    residual_se = s$sigma,
    df_model = unname(s$fstatistic[2]),
    df_resid = unname(s$fstatistic[3]),
    comparison = NA_character_,
    comparison_p = NA_real_
  )
}

extract_model_comparison <- function(model_small, model_large, comparison_name) {
  a <- anova(model_small, model_large)
  
  data.frame(
    model = "model_comparison",
    term = comparison_name,
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = a$F[2],
    p_value = a$`Pr(>F)`[2],
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    aic = NA_real_,
    residual_se = NA_real_,
    df_model = a$Df[2],
    df_resid = a$Res.Df[2],
    comparison = comparison_name,
    comparison_p = a$`Pr(>F)`[2]
  )
}

# ------------------------------------------------------------
# READ DATA
# ------------------------------------------------------------
df_metrics <- read.csv2(metrics_csv, stringsAsFactors = FALSE, check.names = FALSE)
df_key     <- read_table_robust(key_csv)

names(df_metrics) <- clean_names_basic(names(df_metrics))
names(df_key)     <- clean_names_basic(names(df_key))

cat("\n--- Column names in metrics ---\n")
print(names(df_metrics))

cat("\n--- Column names in specimen_key ---\n")
print(names(df_key))

# ------------------------------------------------------------
# CHECK REQUIRED COLUMNS
# ------------------------------------------------------------
needed_metrics <- c("specimen", "bbox_diag_um", "median_thickness_um")
stopifnot(all(needed_metrics %in% names(df_metrics)))

needed_key <- c("specimen_key", "Coxal wall hole", "Coxal Socket")
stopifnot(all(needed_key %in% names(df_key)))

# ------------------------------------------------------------
# CLEAN DATA
# ------------------------------------------------------------
df_metrics <- df_metrics %>%
  mutate(
    specimen = as.character(specimen),
    bbox_diag_um = as.numeric(bbox_diag_um),
    median_thickness_um = as.numeric(median_thickness_um)
  ) %>%
  filter(
    is.finite(bbox_diag_um),
    is.finite(median_thickness_um),
    bbox_diag_um > 0,
    median_thickness_um > 0
  ) %>%
  mutate(
    log_size = log10(bbox_diag_um),
    log_thickness = log10(median_thickness_um)
  )

df_traits <- df_key %>%
  mutate(
    specimen_key = as.character(specimen_key),
    hole = to_logical_trait(`Coxal wall hole`),
    indentation = to_logical_trait(`Coxal Socket`)
  ) %>%
  select(specimen_key, hole, indentation)

# ------------------------------------------------------------
# MATCHING CHECK
# ------------------------------------------------------------
metrics_ids <- unique(df_metrics$specimen)
key_ids     <- unique(df_traits$specimen_key)
n_match     <- sum(metrics_ids %in% key_ids)

cat("\n--- Matching summary ---\n")
cat("Unique metric specimen IDs: ", length(metrics_ids), "\n", sep = "")
cat("Unique specimen_key IDs:    ", length(key_ids), "\n", sep = "")
cat("Matched IDs:                ", n_match, "\n", sep = "")

if (n_match == 0) {
  stop("0 matching IDs between metrics table and specimen_key.csv")
}

# ------------------------------------------------------------
# MERGE
# ------------------------------------------------------------
df <- df_metrics %>%
  left_join(df_traits, by = c("specimen" = "specimen_key"))

write.csv2(df, out_merge_csv, row.names = FALSE)
cat("Wrote debug merged table: ", out_merge_csv, "\n", sep = "")

cat("\n--- Trait availability ---\n")
cat("Rows total:                ", nrow(df), "\n", sep = "")
cat("Rows with hole info:       ", sum(!is.na(df$hole)), "\n", sep = "")
cat("Rows with indentation info:", sum(!is.na(df$indentation)), "\n", sep = "")

# ------------------------------------------------------------
# GROUP COUNTS
# ------------------------------------------------------------
counts_tab <- bind_rows(
  df %>%
    mutate(group = case_when(
      is.na(hole) ~ "NA",
      hole ~ "TRUE",
      !hole ~ "FALSE"
    )) %>%
    group_by(group) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(trait = "hole"),
  
  df %>%
    mutate(group = case_when(
      is.na(indentation) ~ "NA",
      indentation ~ "TRUE",
      !indentation ~ "FALSE"
    )) %>%
    group_by(group) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(trait = "indentation")
) %>%
  select(trait, group, n)

write.csv2(counts_tab, out_counts_csv, row.names = FALSE)
cat("Wrote counts CSV: ", out_counts_csv, "\n", sep = "")

# ------------------------------------------------------------
# ANALYSIS DATA
# ------------------------------------------------------------
df_hole <- df %>%
  filter(!is.na(hole))

cat("\n--- Hole dataset summary ---\n")
cat("n total:  ", nrow(df_hole), "\n", sep = "")
cat("n TRUE:   ", sum(df_hole$hole), "\n", sep = "")
cat("n FALSE:  ", sum(!df_hole$hole), "\n", sep = "")

if (nrow(df_hole) < 10) {
  stop("Too few rows with non-missing hole information.")
}

# ------------------------------------------------------------
# MODELS
# ------------------------------------------------------------
# baseline: only size
m0 <- lm(log_thickness ~ log_size, data = df_hole)

# additive hole effect
m1 <- lm(log_thickness ~ log_size + hole, data = df_hole)

# interaction: different scaling by hole
m2 <- lm(log_thickness ~ log_size * hole, data = df_hole)

# residuals from size-only model
df_hole <- df_hole %>%
  mutate(
    resid_size_only = resid(m0),
    fitted_size_only = fitted(m0),
    fitted_additive = fitted(m1),
    fitted_interaction = fitted(m2)
  )

# ------------------------------------------------------------
# EXPORT STATS
# ------------------------------------------------------------
stats_tab <- bind_rows(
  extract_coefs(m0, "m0_size_only"),
  extract_model_summary(m0, "m0_size_only"),
  
  extract_coefs(m1, "m1_size_plus_hole"),
  extract_model_summary(m1, "m1_size_plus_hole"),
  
  extract_coefs(m2, "m2_size_hole_interaction"),
  extract_model_summary(m2, "m2_size_hole_interaction"),
  
  extract_model_comparison(m0, m1, "m0_vs_m1_add_hole"),
  extract_model_comparison(m1, m2, "m1_vs_m2_add_interaction")
)

write.csv2(stats_tab, out_stats_csv, row.names = FALSE)
cat("Wrote stats CSV: ", out_stats_csv, "\n", sep = "")

cat("\n--- Model summaries ---\n")
print(summary(m0))
print(summary(m1))
print(summary(m2))

cat("\n--- Model comparisons ---\n")
print(anova(m0, m1))
print(anova(m1, m2))

# ------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------

# A) scatter with separate fits by hole
p1 <- ggplot(df_hole, aes(x = bbox_diag_um, y = median_thickness_um, shape = hole)) +
  geom_point(alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Thickness vs size by Coxal wall hole",
    x = "Coxa size (bbox diagonal, m) [log10]",
    y = "Median thickness (m) [log10]",
    shape = "Coxal wall hole"
  ) +
  theme_classic()

# B) raw thickness boxplot by hole
p2 <- ggplot(df_hole, aes(x = factor(hole, levels = c(FALSE, TRUE)), y = median_thickness_um)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7) +
  labs(
    title = "Raw median thickness by Coxal wall hole",
    x = "Coxal wall hole",
    y = "Median thickness (m)"
  ) +
  theme_classic()

# C) size-corrected residuals by hole
p3 <- ggplot(df_hole, aes(x = factor(hole, levels = c(FALSE, TRUE)), y = resid_size_only)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7) +
  labs(
    title = "Size-corrected thickness residuals by Coxal wall hole",
    x = "Coxal wall hole",
    y = "Residual log10(thickness) from size-only model"
  ) +
  theme_classic()

# D) model diagnostic: residuals vs fitted
p4 <- ggplot(
  data.frame(
    fitted = fitted(m1),
    resid = resid(m1)
  ),
  aes(x = fitted, y = resid)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.8) +
  labs(
    title = "Diagnostic plot for additive model",
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_classic()

combined_plot <- (p1 + p2) / (p3 + p4)

ggsave(
  filename = out_plot_png,
  plot = combined_plot,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  filename = out_plot_pdf,
  plot = combined_plot,
  width = 12,
  height = 10
)

cat("Wrote plot PNG: ", out_plot_png, "\n", sep = "")
cat("Wrote plot PDF: ", out_plot_pdf, "\n", sep = "")

# ------------------------------------------------------------
# OPTIONAL: concise console interpretation helper
# ------------------------------------------------------------
hole_p <- summary(m1)$coefficients["holeTRUE", "Pr(>|t|)"]
interaction_p <- summary(m2)$coefficients["log_size:holeTRUE", "Pr(>|t|)"]

cat("\n--- Quick interpretation helper ---\n")
cat("p(hole additive effect)     = ", signif(hole_p, 4), "\n", sep = "")
cat("p(hole:size interaction)    = ", signif(interaction_p, 4), "\n", sep = "")
