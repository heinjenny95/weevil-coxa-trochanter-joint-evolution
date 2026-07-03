# ============================================================
# SCREW JOINT GEOMETRY vs SHAPE
# Unified analysis script
#
# This script does ALL of the following:
# 1) reads PCA and screw-geometry CSV files robustly
# 2) cleans and merges the data
# 3) applies the biologically motivated cutoff angle_abs >= 30
# 4) computes axial pitch
# 5) creates:
#    - MAIN TEXT FIGURE: 1x2 grid
#      (A) Morphospace colored by winding angle, sized by axial pitch
#      (B) Global relationship between PC1 and winding angle
#    - SUPPLEMENT FIGURE: multi-panel QC/interpretation grid
# 6) fits the key regression models
# 7) exports a clean regression summary table as CSV
# 8) exports a specimen-level data table used for plotting
#
# Core biological question:
# Is screw joint geometry (winding angle, axial pitch) coupled to
# overall trochanter shape, or does it vary partly independently?
# ============================================================

rm(list = ls())

# ------------------- INPUT -------------------
pca_file  <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"
geom_file <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/winding_metrics_excel.csv"

# ------------------- OUTPUT FOLDER -------------------
out_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/Shape_vs_ScrewGeometry"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------- OUTPUT FILES -------------------
main_fig_png <- file.path(out_dir, "Figure_main_shape_vs_screw_geometry.png")
main_fig_pdf <- file.path(out_dir, "Figure_main_shape_vs_screw_geometry.pdf")

supp_fig_png <- file.path(out_dir, "Figure_supplement_shape_vs_screw_geometry_QC.png")
supp_fig_pdf <- file.path(out_dir, "Figure_supplement_shape_vs_screw_geometry_QC.pdf")

regression_csv <- file.path(out_dir, "Table_regression_summary.csv")
plot_data_csv   <- file.path(out_dir, "Table_plotting_dataset.csv")
subset_ids_csv  <- file.path(out_dir, "Table_subset_membership.csv")

# ------------------- PACKAGES -------------------
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# ============================================================
# 1) HELPER FUNCTIONS
# ============================================================

# -------------------
# clean_str:
# Cleans string values and headers.
# Why:
# Hidden whitespace / BOM characters often break merges and column matching.
# -------------------
clean_str <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("^\ufeff", "", x)
  trimws(x)
}

# -------------------
# to_num:
# Converts numeric-looking strings safely.
# Why:
# Handles German decimal comma and prevents silent coercion issues.
# -------------------
to_num <- function(x) {
  x <- clean_str(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

# -------------------
# find_col:
# Finds a column from a list of acceptable names.
# Why:
# Allows robust handling of slightly different CSV schemas.
# -------------------
find_col <- function(df, candidates, label = NULL, required = TRUE) {
  nms <- names(df)
  nms_clean <- tolower(trimws(nms))
  cand_clean <- tolower(trimws(candidates))
  hit <- which(nms_clean %in% cand_clean)
  
  if (length(hit) == 0) {
    if (required) {
      if (is.null(label)) label <- paste(candidates, collapse = ", ")
      stop(
        "Could not find column for: ", label,
        "\nAvailable columns are:\n",
        paste(names(df), collapse = " | ")
      )
    } else {
      return(NULL)
    }
  }
  nms[hit[1]]
}

# -------------------
# read_csv_robust:
# Reads CSV files robustly, including Excel-DE style files with "sep=;".
# Why:
# Standard read.csv/read.csv2 often fails silently or inconsistently.
# -------------------
read_csv_robust <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
  
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) stop("File is empty: ", path)
  
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) stop("File contains only empty lines: ", path)
  
  lines[1] <- sub("^\ufeff", "", lines[1])
  
  # Case 1: Excel-style separator declaration
  if (grepl("^sep\\s*=\\s*.", lines[1], ignore.case = TRUE)) {
    sep <- sub("^sep\\s*=\\s*", "", lines[1], ignore.case = TRUE)
    sep <- substr(sep, 1, 1)
    
    tmp <- tempfile(fileext = ".csv")
    writeLines(lines[-1], tmp, useBytes = TRUE)
    
    x <- tryCatch(
      read.table(
        tmp,
        sep = sep,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        dec = ",",
        quote = "\"",
        comment.char = ""
      ),
      error = function(e) NULL
    )
    if (!is.null(x) && ncol(x) > 1) {
      names(x) <- clean_str(names(x))
      return(x)
    }
    
    x <- tryCatch(
      read.table(
        tmp,
        sep = sep,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        dec = ".",
        quote = "\"",
        comment.char = ""
      ),
      error = function(e) NULL
    )
    if (!is.null(x) && ncol(x) > 1) {
      names(x) <- clean_str(names(x))
      return(x)
    }
  }
  
  # Case 2: semicolon
  x <- tryCatch(
    read.table(
      path,
      sep = ";",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      dec = ",",
      quote = "\"",
      comment.char = "",
      fileEncoding = "UTF-8"
    ),
    error = function(e) NULL
  )
  if (!is.null(x) && ncol(x) > 1) {
    names(x) <- clean_str(names(x))
    return(x)
  }
  
  # Case 3: comma
  x <- tryCatch(
    read.table(
      path,
      sep = ",",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      dec = ".",
      quote = "\"",
      comment.char = "",
      fileEncoding = "UTF-8"
    ),
    error = function(e) NULL
  )
  if (!is.null(x) && ncol(x) > 1) {
    names(x) <- clean_str(names(x))
    return(x)
  }
  
  preview <- paste(utils::head(lines, 5), collapse = "\n")
  stop("Could not read CSV: ", path, "\nFirst lines were:\n", preview)
}

# -------------------
# write_csv_de:
# Exports CSV with semicolon separator.
# Why:
# Convenient for Excel in German locale.
# -------------------
write_csv_de <- function(df, path) {
  write.table(
    df,
    file = path,
    sep = ";",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    fileEncoding = "UTF-8"
  )
}

# -------------------
# make_quantile_bins:
# Creates 4 quantile-based size classes.
# Why:
# Makes point-size legend interpretable without being dominated by outliers.
# -------------------
make_quantile_bins <- function(x, labels = c("very low", "low", "high", "very high")) {
  x_ok <- x[is.finite(x) & !is.na(x)]
  if (length(x_ok) < 4) {
    stop("Not enough valid values for binning.")
  }
  
  qs <- quantile(x_ok, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 7)
  qs <- as.numeric(qs)
  
  if (length(unique(qs)) < 5) {
    rng <- range(x_ok, na.rm = TRUE)
    qs <- seq(rng[1], rng[2], length.out = 5)
  }
  
  cut(
    x,
    breaks = qs,
    include.lowest = TRUE,
    labels = labels
  )
}

# -------------------
# lm_summary_row:
# Extracts clean regression summary information into one row.
# Why:
# Produces a publication-friendly summary table.
# -------------------
lm_summary_row <- function(model, model_name, subset_name, n_obs) {
  sm <- summary(model)
  data.frame(
    model = model_name,
    subset = subset_name,
    n = n_obs,
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    f_statistic = unname(sm$fstatistic[1]),
    df1 = unname(sm$fstatistic[2]),
    df2 = unname(sm$fstatistic[3]),
    p_model = pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3], lower.tail = FALSE),
    intercept_estimate = coef(sm)[1, "Estimate"],
    intercept_p = coef(sm)[1, "Pr(>|t|)"],
    PC1_estimate = if ("PC1" %in% rownames(coef(sm))) coef(sm)["PC1", "Estimate"] else NA,
    PC1_p = if ("PC1" %in% rownames(coef(sm))) coef(sm)["PC1", "Pr(>|t|)"] else NA,
    PC2_estimate = if ("PC2" %in% rownames(coef(sm))) coef(sm)["PC2", "Estimate"] else NA,
    PC2_p = if ("PC2" %in% rownames(coef(sm))) coef(sm)["PC2", "Pr(>|t|)"] else NA,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 2) READ AND PREPARE DATA
# ============================================================

# -------------------
# Read the two key datasets:
# - PCA scores = overall shape
# - geometry metrics = winding angle, pitch-related quantities, fit quality
# -------------------
pca  <- read_csv_robust(pca_file)
geom <- read_csv_robust(geom_file)

names(pca)  <- clean_str(names(pca))
names(geom) <- clean_str(names(geom))

# -------------------
# Identify relevant columns robustly.
# Why:
# Input files may have slightly different names depending on previous export steps.
# -------------------
id_pca  <- find_col(pca,  c("specimen_id"), "specimen_id in PCA")
id_geom <- find_col(geom, c("specimen_id"), "specimen_id in geometry")

pc1 <- find_col(pca, c("PC1"), "PC1")
pc2 <- find_col(pca, c("PC2"), "PC2")

# Important:
# We prefer ABSOLUTE winding angle for biological interpretation of "how much screw".
angle_col <- find_col(
  geom,
  c("abs_winding_angle_deg", "winding_angle_deg"),
  "absolute winding angle"
)

signed_angle_col <- find_col(
  geom,
  c("signed_winding_angle_deg"),
  "signed winding angle",
  required = FALSE
)

axial_span_col <- find_col(
  geom,
  c("axial_span"),
  "axial_span",
  required = FALSE
)

start_end_col <- find_col(
  geom,
  c("start_end_dist"),
  "start_end_dist",
  required = FALSE
)

fit_rms_col <- find_col(
  geom,
  c("fit_rms"),
  "fit_rms",
  required = FALSE
)

fit_radius_col <- find_col(
  geom,
  c("fit_radius"),
  "fit_radius",
  required = FALSE
)

if (is.null(axial_span_col) && is.null(start_end_col)) {
  stop("Neither axial_span nor start_end_dist found in geometry file.")
}

# -------------------
# Clean specimen IDs before merging.
# Why:
# Invisible whitespace mismatches can kill the merge.
# -------------------
pca[[id_pca]]   <- clean_str(pca[[id_pca]])
geom[[id_geom]] <- clean_str(geom[[id_geom]])

# -------------------
# Merge PCA and geometry data.
# Why:
# This brings shape and screw-geometry into one analysis table.
# -------------------
geom_keep <- c(id_geom, angle_col)
if (!is.null(signed_angle_col)) geom_keep <- c(geom_keep, signed_angle_col)
if (!is.null(axial_span_col))   geom_keep <- c(geom_keep, axial_span_col)
if (!is.null(start_end_col))    geom_keep <- c(geom_keep, start_end_col)
if (!is.null(fit_rms_col))      geom_keep <- c(geom_keep, fit_rms_col)
if (!is.null(fit_radius_col))   geom_keep <- c(geom_keep, fit_radius_col)

geom_keep <- unique(geom_keep)

df <- merge(
  pca[, c(id_pca, pc1, pc2)],
  geom[, geom_keep],
  by.x = id_pca,
  by.y = id_geom,
  all = FALSE
)

cat("Merged rows:", nrow(df), "\n")

if (nrow(df) == 0) {
  stop("Merge resulted in 0 rows. Check specimen_id formatting.")
}

# -------------------
# Convert key columns to numeric.
# -------------------
num_cols <- c(pc1, pc2, angle_col, signed_angle_col, axial_span_col, start_end_col, fit_rms_col, fit_radius_col)
num_cols <- unique(num_cols[!is.null(num_cols)])

for (cc in num_cols) {
  df[[cc]] <- to_num(df[[cc]])
}

# -------------------
# Define core geometry variables.
# Why:
# angle_abs = amount of winding
# axial_pitch = axial travel normalized per 360
# -------------------
df$angle_abs <- abs(df[[angle_col]])

if (!is.null(signed_angle_col)) {
  df$angle_signed <- df[[signed_angle_col]]
} else {
  df$angle_signed <- NA_real_
}

df$axial_metric <- if (!is.null(axial_span_col)) {
  df[[axial_span_col]]
} else {
  df[[start_end_col]]
}

df$axial_pitch <- ifelse(
  is.na(df$angle_abs) | df$angle_abs == 0,
  NA,
  df$axial_metric * 360 / df$angle_abs
)

if (!is.null(fit_rms_col)) {
  df$fit_rms <- df[[fit_rms_col]]
} else {
  df$fit_rms <- NA_real_
}

if (!is.null(fit_radius_col)) {
  df$fit_radius <- df[[fit_radius_col]]
  df$fit_rms_rel <- ifelse(
    is.na(df$fit_radius) | df$fit_radius == 0,
    NA,
    df$fit_rms / df$fit_radius
  )
} else {
  df$fit_radius <- NA_real_
  df$fit_rms_rel <- NA_real_
}

# ============================================================
# 3) BIOLOGICALLY MOTIVATED FILTERING
# ============================================================

# -------------------
# Why angle_abs >= 30?
# Very small winding angles do not represent a meaningful screw-like geometry
# and produce unstable pitch estimates.
# -------------------
df <- df %>%
  filter(
    !is.na(.data[[pc1]]),
    !is.na(.data[[pc2]]),
    !is.na(angle_abs),
    !is.na(axial_pitch),
    is.finite(axial_pitch),
    axial_pitch > 0,
    angle_abs >= 30
  )

cat("Rows after final filtering:", nrow(df), "\n")

if (nrow(df) < 6) {
  stop("Too few rows left after filtering.")
}

# -------------------
# Define a coarse "main regime" subset for interpretation.
# Why:
# Global correlation may be driven by specimens occupying a distinct PC1 region.
# This split is NOT used to claim clusters; it is used to test whether
# the full-dataset trend is regime-driven.
# -------------------
df$shape_regime <- ifelse(df[[pc1]] < 0.1, "main_region", "right_region")

df_left  <- df %>% filter(shape_regime == "main_region")
df_right <- df %>% filter(shape_regime == "right_region")

# -------------------
# Create size classes for the main morphospace plot.
# Why:
# Easier to interpret than a raw continuous size legend for pitch.
# -------------------
df$axial_class <- make_quantile_bins(df$axial_pitch)
df$axial_class <- factor(
  df$axial_class,
  levels = c("very low", "low", "high", "very high")
)

size_values <- c(
  "very low"  = 2.5,
  "low"       = 4.2,
  "high"      = 5.8,
  "very high" = 7.2
)

# ============================================================
# 4) STATISTICAL MODELS
# ============================================================

# -------------------
# Full dataset:
# Tests whether geometry is associated with shape globally.
# -------------------
model_angle_full <- lm(angle_abs ~ PC1 + PC2, data = df)
model_pitch_full <- lm(axial_pitch ~ PC1 + PC2, data = df)

# -------------------
# Main-region subset:
# Tests whether the global trend persists within the dominant shape regime.
# This is the key decoupling test.
# -------------------
model_angle_left <- lm(angle_abs ~ PC1 + PC2, data = df_left)

# Optional right-region model if enough data
model_angle_right <- NULL
if (nrow(df_right) >= 5) {
  model_angle_right <- lm(angle_abs ~ PC1 + PC2, data = df_right)
}

# -------------------
# Build regression summary table.
# Why:
# This is the clean exportable stats table for the manuscript.
# -------------------
reg_table <- bind_rows(
  lm_summary_row(model_angle_full, "angle_abs ~ PC1 + PC2", "full_dataset", nrow(df)),
  lm_summary_row(model_pitch_full, "axial_pitch ~ PC1 + PC2", "full_dataset", nrow(df)),
  lm_summary_row(model_angle_left, "angle_abs ~ PC1 + PC2", "main_region_PC1_lt_0.1", nrow(df_left)),
  if (!is.null(model_angle_right)) lm_summary_row(model_angle_right, "angle_abs ~ PC1 + PC2", "right_region_PC1_ge_0.1", nrow(df_right))
)

# -------------------
# Correlation summary added for convenience.
# Why:
# Helps quick interpretation without opening model objects.
# -------------------
corr_table <- data.frame(
  metric = c("angle_abs ~ PC1", "angle_abs ~ PC2", "axial_pitch ~ PC1", "axial_pitch ~ PC2"),
  correlation = c(
    cor(df$angle_abs, df$PC1, use = "complete.obs"),
    cor(df$angle_abs, df$PC2, use = "complete.obs"),
    cor(df$axial_pitch, df$PC1, use = "complete.obs"),
    cor(df$axial_pitch, df$PC2, use = "complete.obs")
  ),
  stringsAsFactors = FALSE
)

# We combine both into one export by stacking a block label.
# Why:
# Regression models and simple correlations are both useful, but they have
# different columns. We therefore harmonize them before binding rows.

corr_export <- corr_table %>%
  mutate(
    table_block = "correlations",
    subset = NA_character_,
    n = NA_real_,
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    f_statistic = NA_real_,
    df1 = NA_real_,
    df2 = NA_real_,
    p_model = NA_real_,
    intercept_estimate = NA_real_,
    intercept_p = NA_real_,
    PC1_estimate = NA_real_,
    PC1_p = NA_real_,
    PC2_estimate = NA_real_,
    PC2_p = NA_real_
  ) %>%
  rename(model = metric) %>%
  mutate(correlation_value = correlation) %>%
  select(
    table_block, model, subset, n,
    r_squared, adj_r_squared, f_statistic, df1, df2, p_model,
    intercept_estimate, intercept_p,
    PC1_estimate, PC1_p,
    PC2_estimate, PC2_p,
    correlation_value
  )

reg_export_models <- reg_table %>%
  mutate(
    table_block = "regression_models",
    correlation_value = NA_real_
  ) %>%
  select(
    table_block, model, subset, n,
    r_squared, adj_r_squared, f_statistic, df1, df2, p_model,
    intercept_estimate, intercept_p,
    PC1_estimate, PC1_p,
    PC2_estimate, PC2_p,
    correlation_value
  )

reg_export <- bind_rows(
  reg_export_models,
  corr_export
)
# ============================================================
# 5) MAIN-TEXT FIGURES
# ============================================================

# -------------------
# MAIN FIGURE PANEL A
# Morphospace:
# color = winding angle
# size  = axial pitch class
#
# Why:
# This is the main visual overview of how screw geometry is distributed
# across overall shape space.
# -------------------
p_main_A <- ggplot(df, aes(x = .data[[pc1]], y = .data[[pc2]])) +
  geom_point(
    aes(size = axial_class, color = angle_abs),
    alpha = 0.9
  ) +
  scale_size_manual(
    values = size_values,
    drop = FALSE,
    name = "Axial pitch per 360"
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Winding angle ()"
  ) +
  labs(
    title = "A. Morphospace of screw joint geometry",
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  )

# -------------------
# MAIN FIGURE PANEL B
# Global relationship between shape and winding angle.
#
# Why:
# This shows the apparent global shape-geometry association and motivates
# the regime-specific follow-up interpretation.
# -------------------
p_main_B <- ggplot(df, aes(x = .data[[pc1]], y = angle_abs)) +
  geom_point(alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    title = "B. Global relationship between PC1 and winding angle",
    x = "PC1",
    y = "Absolute winding angle ()"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

main_grob <- arrangeGrob(
  p_main_A, p_main_B,
  ncol = 2
)

# ============================================================
# 6) SUPPLEMENT FIGURES
# ============================================================

# -------------------
# SUPP PANEL 1
# Distribution of winding angle
#
# Why:
# Shows whether geometry is spread smoothly or concentrated in broad regimes.
# -------------------
p_supp_1 <- ggplot(df, aes(x = angle_abs)) +
  geom_histogram(bins = 20) +
  labs(
    title = "S1. Distribution of absolute winding angle",
    x = "Absolute winding angle ()",
    y = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# -------------------
# SUPP PANEL 2
# Winding angle vs fit RMS
#
# Why:
# Checks whether extreme angles are merely fit artefacts.
# -------------------
p_supp_2 <- ggplot(df, aes(x = angle_abs, y = fit_rms)) +
  geom_point(alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  labs(
    title = "S2. Winding angle vs fit RMS",
    x = "Absolute winding angle ()",
    y = "fit_rms"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# -------------------
# SUPP PANEL 3
# Morphospace colored by fit RMS
#
# Why:
# Tests whether particular regions of morphospace are dominated by poor fits.
# -------------------
p_supp_3 <- ggplot(df, aes(x = .data[[pc1]], y = .data[[pc2]])) +
  geom_point(aes(color = fit_rms), size = 4.0, alpha = 0.9) +
  scale_color_viridis_c(
    option = "magma",
    name = "fit_rms"
  ) +
  labs(
    title = "S3. Morphospace colored by fit RMS",
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# -------------------
# SUPP PANEL 4
# Morphospace with angle + fit quality
#
# Why:
# Joint overview of biological signal (angle) and fit quality (size).
# -------------------
p_supp_4 <- ggplot(df, aes(x = .data[[pc1]], y = .data[[pc2]])) +
  geom_point(aes(color = angle_abs, size = fit_rms), alpha = 0.9) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Winding angle ()"
  ) +
  scale_size_continuous(
    name = "fit_rms"
  ) +
  labs(
    title = "S4. Morphospace: winding angle and fit RMS",
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  )

# -------------------
# SUPP PANEL 5
# Relative fit RMS
#
# Why:
# Raw RMS can depend on scale. Relative RMS is often more informative.
# -------------------
p_supp_5 <- ggplot(df, aes(x = angle_abs, y = fit_rms_rel)) +
  geom_point(alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  labs(
    title = "S5. Winding angle vs relative fit RMS",
    x = "Absolute winding angle ()",
    y = "fit_rms / fit_radius"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# -------------------
# SUPP PANEL 6
# Main-region-only relationship
#
# Why:
# This is the key decoupling test:
# does winding angle still track shape within the main regime?
# -------------------
p_supp_6 <- ggplot(df_left, aes(x = .data[[pc1]], y = angle_abs)) +
  geom_point(alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    title = "S6. Main-region relationship: PC1 vs winding angle",
    x = "PC1",
    y = "Absolute winding angle ()"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

supp_grob <- arrangeGrob(
  p_supp_1, p_supp_2,
  p_supp_3, p_supp_4,
  p_supp_5, p_supp_6,
  ncol = 2
)

# ============================================================
# 7) EXPORT TABLES (FINAL - Excel-safe & clean)
# ============================================================

# -------------------
# Helper: safe CSV export for Excel (DE)
# Why:
# - uses semicolon separator (Excel expects this)
# - keeps decimal point (scientific correctness)
# -------------------
write_csv_clean <- function(df, path) {
  write.table(
    df,
    file = path,
    sep = ";",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    dec = ".",
    fileEncoding = "UTF-8"
  )
}

# -------------------
# OPTIONAL: round values for readability (paper-ready)
# Why:
# Raw numbers are ugly and hard to read
# -------------------
reg_export_clean <- reg_export %>%
  mutate(
    r_squared = round(r_squared, 3),
    adj_r_squared = round(adj_r_squared, 3),
    f_statistic = round(f_statistic, 3),
    p_model = signif(p_model, 3),
    
    intercept_estimate = round(intercept_estimate, 3),
    intercept_p = signif(intercept_p, 3),
    
    PC1_estimate = round(PC1_estimate, 3),
    PC1_p = signif(PC1_p, 3),
    
    PC2_estimate = round(PC2_estimate, 3),
    PC2_p = signif(PC2_p, 3),
    
    correlation_value = signif(correlation_value, 3)
  )

# -------------------
# 1) Regression summary table
# -------------------
write_csv_clean(reg_export_clean, regression_csv)

# -------------------
# 2) Specimen-level plotting dataset
# -------------------
plot_export <- df %>%
  select(
    specimen_id = all_of(id_pca),
    PC1 = all_of(pc1),
    PC2 = all_of(pc2),
    angle_abs,
    angle_signed,
    axial_metric,
    axial_pitch,
    fit_rms,
    fit_radius,
    fit_rms_rel,
    axial_class,
    shape_regime
  )

write_csv_clean(plot_export, plot_data_csv)

# -------------------
# 3) Subset membership (debug / interpretation)
# -------------------
subset_export <- df %>%
  select(
    specimen_id = all_of(id_pca),
    PC1 = all_of(pc1),
    PC2 = all_of(pc2),
    shape_regime
  ) %>%
  arrange(shape_regime, PC1)

write_csv_clean(subset_export, subset_ids_csv)

# -------------------
# Console confirmation
# -------------------
cat("\nTables successfully exported:\n")
cat(regression_csv, "\n")
cat(plot_data_csv, "\n")
cat(subset_ids_csv, "\n")

# ============================================================
# 8) SAVE FIGURES
# ============================================================

ggsave(main_fig_png, main_grob, width = 14, height = 6.5, dpi = 400)
ggsave(main_fig_pdf, main_grob, width = 14, height = 6.5)

ggsave(supp_fig_png, supp_grob, width = 14, height = 18, dpi = 400)
ggsave(supp_fig_pdf, supp_grob, width = 14, height = 18)

# ============================================================
# 9) CONSOLE OUTPUT FOR INTERPRETATION
# ============================================================

cat("\n====================\n")
cat("KEY MODEL SUMMARIES\n")
cat("====================\n\n")

cat("--- Full dataset: angle_abs ~ PC1 + PC2 ---\n")
print(summary(model_angle_full))

cat("\n--- Full dataset: axial_pitch ~ PC1 + PC2 ---\n")
print(summary(model_pitch_full))

cat("\n--- Main region only: angle_abs ~ PC1 + PC2 ---\n")
print(summary(model_angle_left))

if (!is.null(model_angle_right)) {
  cat("\n--- Right region only: angle_abs ~ PC1 + PC2 ---\n")
  print(summary(model_angle_right))
}

cat("\n====================\n")
cat("FILES WRITTEN\n")
cat("====================\n")
cat(main_fig_png, "\n")
cat(main_fig_pdf, "\n")
cat(supp_fig_png, "\n")
cat(supp_fig_pdf, "\n")
cat(regression_csv, "\n")
cat(plot_data_csv, "\n")
cat(subset_ids_csv, "\n")



##########################################################################

# ============================================================
# 7.2) EXPORT MAIN TABLE (CLEAN & PAPER-READY)
# ============================================================

# -------------------
# Build minimal main-table
# Why:
# Only the core results needed for the manuscript
# -------------------

main_table <- bind_rows(
  data.frame(
    model = "Winding angle ~ shape",
    subset = "Full dataset",
    n = nrow(df),
    r_squared = summary(model_angle_full)$r.squared,
    p_value = summary(model_angle_full)$coefficients[2,4],
    PC1_effect = summary(model_angle_full)$coefficients[2,1],
    PC1_p = summary(model_angle_full)$coefficients[2,4]
  ),
  
  data.frame(
    model = "Axial pitch ~ shape",
    subset = "Full dataset",
    n = nrow(df),
    r_squared = summary(model_pitch_full)$r.squared,
    p_value = summary(model_pitch_full)$coefficients[2,4],
    PC1_effect = summary(model_pitch_full)$coefficients[2,1],
    PC1_p = summary(model_pitch_full)$coefficients[2,4]
  ),
  
  data.frame(
    model = "Winding angle ~ shape",
    subset = "Main region (PC1 < 0.1)",
    n = nrow(df_left),
    r_squared = summary(model_angle_left)$r.squared,
    p_value = summary(model_angle_left)$coefficients[2,4],
    PC1_effect = summary(model_angle_left)$coefficients[2,1],
    PC1_p = summary(model_angle_left)$coefficients[2,4]
  )
)

# -------------------
# Round values for readability
# -------------------
main_table <- main_table %>%
  mutate(
    r_squared = round(r_squared, 3),
    p_value = signif(p_value, 3),
    PC1_effect = round(PC1_effect, 2),
    PC1_p = signif(PC1_p, 3)
  )

# -------------------
# Export clean table (Excel-safe)
# -------------------
write.table(
  main_table,
  file = file.path(out_dir, "Table_main_results.csv"),
  sep = ";",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  dec = ".",
  fileEncoding = "UTF-8"
)

cat("\nMain results table exported.\n")


######################################################################################



# ============================================================
# EXTRA FIGURE: axial span & pitch vs winding angle
# ============================================================

# -------------------
# Why:
# - explores relationship between raw axial displacement (span)
#   and normalized displacement (pitch)
# - helps interpret how different geometry components interact
# -------------------

library(scales)

# optional: log-transform (empfohlen, weil Werte oft stark streuen)
df$log_axial_span  <- log10(df$axial_metric)
df$log_axial_pitch <- log10(df$axial_pitch)

# -------------------
# PANEL A: axial span vs winding angle
# -------------------
p_span <- ggplot(df, aes(x = angle_abs, y = log_axial_span)) +
  geom_point(alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    title = "A. Axial span vs winding angle",
    x = "Winding angle ()",
    y = "Axial span (log10)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# -------------------
# PANEL B: axial pitch vs winding angle
# -------------------
p_pitch <- ggplot(df, aes(x = angle_abs, y = log_axial_pitch)) +
  geom_point(alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    title = "B. Axial pitch vs winding angle",
    x = "Winding angle ()",
    y = "Axial pitch per 360 (log10)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# -------------------
# Combine into 1x2 figure
# -------------------
fig_axial <- gridExtra::arrangeGrob(
  p_span, p_pitch,
  ncol = 2
)

# -------------------
# Save
# -------------------
ggsave(
  file.path(out_dir, "Figure_supplement_axial_relationships.png"),
  fig_axial,
  width = 12,
  height = 5.5,
  dpi = 400
)

ggsave(
  file.path(out_dir, "Figure_supplement_axial_relationships.pdf"),
  fig_axial,
  width = 12,
  height = 5.5
)

cat("\nAxial span/pitch figure saved.\n")
