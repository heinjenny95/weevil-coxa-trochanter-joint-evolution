# ============================================================
# Joint type Ã— screw geometry (NON-phylogenetic) â€” UPDATED VERSION
#
# Compatible with new geometry file containing:
#   specimen_id
#   signed_winding_angle_deg
#   abs_winding_angle_deg
#   n_turns_signed
#   n_turns_abs
#   start_end_dist
#   axial_span
#   fit_radius
#   fit_rms
#
# And joint type file containing at least:
#   specimen_id
#   screw_state
#   joint_type
#   joint_type_strict   (optional)
#
# Computes:
# - axial_pitch_360 = axial_span * 360 / abs_winding_angle_deg
#
# Analyses (ONLY for screw joints):
# - Combined 1x3 figure:
#     abs_winding_angle_deg, axial_pitch_360, axial_span by joint_type
# - Kruskal-Wallis + post hoc
# - PERMANOVA (Euclidean) on scaled geometry vector ~ joint_type
# - Optional: screw_state tests
#
# Important:
# - Uses ABSOLUTE winding angle for all biological analyses
# - Filters out angle < 30Â° by default, because those cases do not
#   represent robust screw-like geometry and destabilize pitch estimates
# ============================================================

rm(list = ls())

# ------------------- PATHS -------------------
geom_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/winding_metrics_excel.csv"
type_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_joint_types.csv"

out_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/JointType_Geometry"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_fig_png <- file.path(out_dir, "joint_type_screw_geometry_1x3.png")
out_fig_pdf <- file.path(out_dir, "joint_type_screw_geometry_1x3.pdf")

out_stats_csv <- file.path(out_dir, "joint_type_screw_geometry_stats.csv")
out_plotdata_csv <- file.path(out_dir, "joint_type_screw_geometry_plotdata.csv")

# ------------------- SETTINGS -------------------
ANGLE_CUTOFF_DEG <- 30

# ------------------- PACKAGES -------------------
pkgs <- c("ggplot2", "FSA", "vegan", "patchwork", "dplyr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------- HELPERS -------------------
clean_str <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("^\ufeff", "", x)
  trimws(x)
}

to_num <- function(x) {
  x <- clean_str(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

read_csv_robust <- function(path) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) stop("File is empty: ", path)
  
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) stop("File contains only empty lines: ", path)
  
  lines[1] <- sub("^\ufeff", "", lines[1])
  
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
  
  stop("Could not read CSV: ", path)
}

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

is_screw_joint <- function(jt) {
  jt2 <- tolower(clean_str(jt))
  grepl("screw", jt2)
}

posthoc_nonparam <- function(df, response, group) {
  g <- factor(df[[group]])
  ng <- nlevels(g)
  
  if (ng < 2) {
    cat("Post-hoc skipped: <2 groups.\n")
    return(invisible(NULL))
  }
  
  if (ng == 2) {
    cat("\nTwo groups -> Wilcoxon rank-sum:", response, "\n")
    form <- as.formula(paste(response, "~", group))
    print(wilcox.test(form, data = df, exact = FALSE))
    return(invisible(NULL))
  }
  
  cat("\n>2 groups -> Dunn test (BH):", response, "\n")
  form <- as.formula(paste(response, "~", group))
  print(FSA::dunnTest(form, data = df, method = "bh"))
  invisible(NULL)
}

make_boxplot <- function(data, yvar, ylab, title_text, palette_vals) {
  ggplot(data, aes(x = joint_type, y = .data[[yvar]], fill = joint_type, color = joint_type)) +
    geom_boxplot(
      width = 0.6,
      alpha = 0.45,
      outlier.shape = NA,
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.12,
      size = 2,
      alpha = 0.8,
      stroke = 0
    ) +
    scale_fill_manual(values = palette_vals) +
    scale_color_manual(values = palette_vals) +
    labs(
      title = title_text,
      x = NULL,
      y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.position = "none"
    )
}

# small helper for stats export
extract_kw <- function(test_obj, response, grouping) {
  data.frame(
    analysis = "kruskal_wallis",
    response = response,
    grouping = grouping,
    statistic = unname(test_obj$statistic),
    df = unname(test_obj$parameter),
    p_value = unname(test_obj$p.value),
    stringsAsFactors = FALSE
  )
}

extract_wilcox <- function(test_obj, response, grouping) {
  data.frame(
    analysis = "wilcoxon",
    response = response,
    grouping = grouping,
    statistic = unname(test_obj$statistic),
    df = NA_real_,
    p_value = unname(test_obj$p.value),
    stringsAsFactors = FALSE
  )
}

# ------------------- READ -------------------
geom <- read_csv_robust(geom_path)
typ  <- read_csv_robust(type_path)

names(geom) <- clean_str(names(geom))
names(typ)  <- clean_str(names(typ))

# ------------------- CHECK REQUIRED COLUMNS -------------------
req_geom <- c("specimen_id", "abs_winding_angle_deg", "start_end_dist", "axial_span")
req_typ  <- c("specimen_id", "joint_type", "screw_state")

missing_geom <- setdiff(req_geom, names(geom))
missing_typ  <- setdiff(req_typ, names(typ))

if (length(missing_geom) > 0) stop("Missing in geometry file: ", paste(missing_geom, collapse = ", "))
if (length(missing_typ)  > 0) stop("Missing in specimen_joint_types.csv: ", paste(missing_typ, collapse = ", "))

# ------------------- CLEAN + NUMERIC -------------------
geom$specimen_id <- clean_str(geom$specimen_id)
typ$specimen_id  <- clean_str(typ$specimen_id)

geom$abs_winding_angle_deg    <- to_num(geom$abs_winding_angle_deg)
geom$signed_winding_angle_deg <- if ("signed_winding_angle_deg" %in% names(geom)) to_num(geom$signed_winding_angle_deg) else NA_real_
geom$start_end_dist           <- to_num(geom$start_end_dist)
geom$axial_span               <- to_num(geom$axial_span)
geom$fit_radius               <- if ("fit_radius" %in% names(geom)) to_num(geom$fit_radius) else NA_real_
geom$fit_rms                  <- if ("fit_rms" %in% names(geom)) to_num(geom$fit_rms) else NA_real_

typ$joint_type <- clean_str(typ$joint_type)
typ$screw_state <- clean_str(typ$screw_state)
typ$joint_type_strict <- if ("joint_type_strict" %in% names(typ)) clean_str(typ$joint_type_strict) else NA_character_

# ------------------- MERGE -------------------
df_all <- merge(
  geom,
  typ[, c("specimen_id", "joint_type", "joint_type_strict", "screw_state")],
  by = "specimen_id",
  all = FALSE
)

if (nrow(df_all) == 0) stop("Merge produced 0 rows. specimen_id mismatch?")

# ------------------- COMPUTE AXIAL PITCH -------------------
# Uses ABSOLUTE winding angle to avoid sign artefacts
df_all$axial_pitch_360 <- NA_real_
ok <- !is.na(df_all$axial_span) & !is.na(df_all$abs_winding_angle_deg) & df_all$abs_winding_angle_deg != 0
df_all$axial_pitch_360[ok] <- df_all$axial_span[ok] * 360 / df_all$abs_winding_angle_deg[ok]

# ------------------- BASIC GEOMETRY FILTER -------------------
df_all <- df_all %>%
  filter(
    !is.na(abs_winding_angle_deg),
    !is.na(axial_pitch_360),
    !is.na(joint_type)
  )

cat("\n==============================\n")
cat("MERGED DATA SUMMARY (geometry + joint types)\n")
cat("==============================\n")
cat("Rows after merge + NA filtering:", nrow(df_all), "\n\n")

cat("Joint types present in merged data:\n")
print(sort(unique(df_all$joint_type)))

cat("\nCounts by joint_type (merged geometry subset):\n")
print(table(df_all$joint_type))

# ------------------- FILTER: ONLY SCREW JOINTS -------------------
df <- df_all %>%
  filter(is_screw_joint(joint_type))

cat("\n==============================\n")
cat("SCREW-GEOMETRY ANALYSIS SUBSET (ONLY screw joints)\n")
cat("==============================\n")
cat("Rows kept before angle cutoff:", nrow(df), "\n")

excluded <- setdiff(unique(df_all$joint_type), unique(df$joint_type))
cat("\nExcluded joint types (non-screw) that had geometry entries:\n")
if (length(excluded) == 0) cat("None.\n") else print(excluded)

# ------------------- ANGLE CUTOFF -------------------
# Important:
# Angle < 30Â° does not represent robust screw-like geometry and makes
# pitch unstable. These cases are removed here.
df <- df %>%
  filter(abs_winding_angle_deg >= ANGLE_CUTOFF_DEG)

cat("\nRows kept after angle cutoff (>= ", ANGLE_CUTOFF_DEG, "Â°): ", nrow(df), "\n", sep = "")

if (nrow(df) < 3) stop("Too few screw-joint rows after filtering.")

# ------------------- FACTORS -------------------
df$joint_type  <- factor(df$joint_type)
df$screw_state <- factor(df$screw_state)

cat("\nCounts by joint_type (screw subset):\n")
print(table(df$joint_type))

cat("\nCounts by screw_state (screw subset):\n")
print(table(df$screw_state))

preferred_order <- c("True screwâ€“nut joint", "True screw-nut joint", "Unopposed screw joint")
present_levels <- levels(df$joint_type)
matched <- preferred_order[preferred_order %in% present_levels]
remaining <- setdiff(present_levels, matched)
new_order <- c(matched, remaining)
df$joint_type <- factor(df$joint_type, levels = new_order)

n_groups <- nlevels(df$joint_type)
palette_vals <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860")[seq_len(n_groups)]
names(palette_vals) <- levels(df$joint_type)

# ------------------- EXPORT PLOTTING DATA -------------------
plot_export <- df %>%
  select(
    specimen_id,
    joint_type,
    joint_type_strict,
    screw_state,
    abs_winding_angle_deg,
    signed_winding_angle_deg,
    axial_span,
    axial_pitch_360,
    start_end_dist,
    fit_radius,
    fit_rms
  )

write_csv_clean(plot_export, out_plotdata_csv)

# ============================================================
# COMBINED FIGURE (1 x 3)
# ============================================================

p1 <- make_boxplot(
  data = df,
  yvar = "abs_winding_angle_deg",
  ylab = "Absolute winding angle (Â°)",
  title_text = "Winding angle",
  palette_vals = palette_vals
)

p2 <- make_boxplot(
  data = df,
  yvar = "axial_pitch_360",
  ylab = "Axial pitch per 360Â°",
  title_text = "Axial pitch",
  palette_vals = palette_vals
)

if (sum(!is.na(df$axial_span)) > 0) {
  p3 <- make_boxplot(
    data = df,
    yvar = "axial_span",
    ylab = "Axial span",
    title_text = "Axial span",
    palette_vals = palette_vals
  )
} else {
  p3 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No axial span data", size = 5) +
    theme_void() +
    labs(title = "Axial span")
}

p_combined <- p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

print(p_combined)

ggsave(out_fig_png, p_combined, width = 15, height = 5.5, dpi = 400)
ggsave(out_fig_pdf, p_combined, width = 15, height = 5.5)

cat("\nCombined figure written to:\n", out_fig_png, "\n", out_fig_pdf, "\n")

# ============================================================
# 1) UNIVARIATE TESTS
# ============================================================

cat("\n==============================\n")
cat("1) Non-parametric tests (screw joints only)\n")
cat("==============================\n")

stats_list <- list()

kw_angle <- kruskal.test(abs_winding_angle_deg ~ joint_type, data = df)
kw_pitch <- kruskal.test(axial_pitch_360 ~ joint_type, data = df)

cat("\nKruskal-Wallis: abs_winding_angle_deg ~ joint_type\n")
print(kw_angle)
stats_list[[length(stats_list) + 1]] <- extract_kw(kw_angle, "abs_winding_angle_deg", "joint_type")

cat("\nKruskal-Wallis: axial_pitch_360 ~ joint_type\n")
print(kw_pitch)
stats_list[[length(stats_list) + 1]] <- extract_kw(kw_pitch, "axial_pitch_360", "joint_type")

if (kw_angle$p.value < 0.05) posthoc_nonparam(df, "abs_winding_angle_deg", "joint_type")
if (kw_pitch$p.value < 0.05) posthoc_nonparam(df, "axial_pitch_360", "joint_type")

if (sum(!is.na(df$axial_span)) > 0) {
  cat("\nKruskal-Wallis: axial_span ~ joint_type\n")
  kw_span <- kruskal.test(axial_span ~ joint_type, data = df)
  print(kw_span)
  stats_list[[length(stats_list) + 1]] <- extract_kw(kw_span, "axial_span", "joint_type")
  
  if (kw_span$p.value < 0.05) posthoc_nonparam(df, "axial_span", "joint_type")
}

# ============================================================
# 2) MULTIVARIATE: PERMANOVA
# ============================================================

cat("\n==============================\n")
cat("2) PERMANOVA (Euclidean; scaled geometry) ~ joint_type\n")
cat("==============================\n")

geom_cols <- c("abs_winding_angle_deg", "axial_pitch_360")
if (sum(!is.na(df$axial_span)) > 0) geom_cols <- c(geom_cols, "axial_span")

geom_mat <- scale(df[, geom_cols, drop = FALSE])
dist_euc <- dist(geom_mat, method = "euclidean")

disp <- vegan::betadisper(dist_euc, df$joint_type)
cat("\nHomogeneity of dispersion (betadisper ANOVA):\n")
disp_anova <- anova(disp)
print(disp_anova)

cat("\nPERMANOVA (adonis2, Euclidean):\n")
perm <- vegan::adonis2(dist_euc ~ joint_type, data = df, permutations = 999)
print(perm)

# add simple export rows
stats_list[[length(stats_list) + 1]] <- data.frame(
  analysis = "betadisper_anova",
  response = paste(geom_cols, collapse = " + "),
  grouping = "joint_type",
  statistic = disp_anova$`F value`[1],
  df = disp_anova$Df[1],
  p_value = disp_anova$`Pr(>F)`[1],
  stringsAsFactors = FALSE
)

stats_list[[length(stats_list) + 1]] <- data.frame(
  analysis = "permanova_adonis2",
  response = paste(geom_cols, collapse = " + "),
  grouping = "joint_type",
  statistic = perm$F[1],
  df = perm$Df[1],
  p_value = perm$`Pr(>F)`[1],
  stringsAsFactors = FALSE
)

# ============================================================
# 3) OPTIONAL: screw_state within screw joints
# ============================================================

cat("\n==============================\n")
cat("3) screw_state (clear vs ambiguous) within screw joints\n")
cat("==============================\n")

if (length(unique(df$screw_state[!is.na(df$screw_state)])) >= 2) {
  cat("\nWilcoxon: abs_winding_angle_deg ~ screw_state\n")
  w1 <- wilcox.test(abs_winding_angle_deg ~ screw_state, data = df, exact = FALSE)
  print(w1)
  stats_list[[length(stats_list) + 1]] <- extract_wilcox(w1, "abs_winding_angle_deg", "screw_state")
  
  cat("\nWilcoxon: axial_pitch_360 ~ screw_state\n")
  w2 <- wilcox.test(axial_pitch_360 ~ screw_state, data = df, exact = FALSE)
  print(w2)
  stats_list[[length(stats_list) + 1]] <- extract_wilcox(w2, "axial_pitch_360", "screw_state")
  
  if (sum(!is.na(df$axial_span)) > 0) {
    cat("\nWilcoxon: axial_span ~ screw_state\n")
    w3 <- wilcox.test(axial_span ~ screw_state, data = df, exact = FALSE)
    print(w3)
    stats_list[[length(stats_list) + 1]] <- extract_wilcox(w3, "axial_span", "screw_state")
  }
} else {
  cat("Only one screw_state present in screw subset; skipping.\n")
}

# ============================================================
# 4) EXPORT STATS
# ============================================================

stats_export <- bind_rows(stats_list) %>%
  mutate(
    statistic = round(statistic, 3),
    p_value = signif(p_value, 3)
  )

write_csv_clean(stats_export, out_stats_csv)

cat("\nStats written to:\n", out_stats_csv, "\n")
cat("Plotting data written to:\n", out_plotdata_csv, "\n")
cat("\n=== DONE ===\n")