# ============================================================
# Centroid Size (Landmark .txt) -> merge into key + PCA
# + Allometry tests on Deformetrica morphospace
#
# Robust ID matching:
#   1) full specimen_id match (including leading numeric ID)
#   2) fallback match by taxon-name only (strip leading numeric ID)
#
# Inputs (German CSV; sep=';' dec=','):
#   - specimen_key.csv (must contain column: specimen_id)
#   - PCA_scores_with_specimen_id.csv (must contain column: specimen_id)
#
# Landmark files:
#   - folder with .txt files like "15_Mononychus_punctumalbum_trochanter.txt"
#
# Outputs (all in folder "Allometry"):
#   - specimen_key_with_centroid_size.csv
#   - PCA_scores_with_specimen_id_with_centroid_size.csv
#   - allometry_univariate_PC1_to_PC5_results.csv
#   - allometry_rrpp_multivariate_results.csv
#   - allometry_procD_lm_results.csv
#   - allometry_univariate_PC1_to_PC5.png/.pdf
#   - morphospace_PC1_PC2_logCS.png/.pdf
#   - allometry_combined.png/.pdf
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(RRPP)
  library(geomorph)
  library(viridis)
  library(patchwork)
})

# ------------------- PATHS -------------------
lm_dir   <- "<BEETLE_JOINTS_ROOT>/Processed/Curculionoidea/Landmarks"
key_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_key.csv"
pca_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"

# base output directory
base_results_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results"
out_dir <- file.path(base_results_dir, "Allometry")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# outputs in Allometry folder
out_key_csv            <- file.path(out_dir, "specimen_key_with_centroid_size.csv")
out_pca_csv            <- file.path(out_dir, "PCA_scores_with_specimen_id_with_centroid_size.csv")
out_univar_csv         <- file.path(out_dir, "allometry_univariate_PC1_to_PC5_results.csv")
out_rrpp_csv           <- file.path(out_dir, "allometry_rrpp_multivariate_results.csv")
out_procD_csv          <- file.path(out_dir, "allometry_procD_lm_results.csv")

out_plot_uni_png       <- file.path(out_dir, "allometry_univariate_PC1_to_PC5.png")
out_plot_uni_pdf       <- file.path(out_dir, "allometry_univariate_PC1_to_PC5.pdf")
out_plot_morph_png     <- file.path(out_dir, "morphospace_PC1_PC2_logCS.png")
out_plot_morph_pdf     <- file.path(out_dir, "morphospace_PC1_PC2_logCS.pdf")
out_plot_combined_png  <- file.path(out_dir, "allometry_combined.png")
out_plot_combined_pdf  <- file.path(out_dir, "allometry_combined.pdf")

# ------------------- SETTINGS -------------------
pcs_use          <- paste0("PC", 1:5)
n_permutations   <- 9999
p_adjust_method  <- "holm"

# ------------------- HELPERS -------------------
stop_with_hint <- function(msg) {
  stop(paste0("\n ", msg, "\n"), call. = FALSE)
}

norm_id <- function(x) {
  x <- as.character(x)
  x <- tolower(x)
  x <- sub("\\.txt$", "", x, ignore.case = TRUE)
  x <- sub("\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("\\.vtk$", "", x, ignore.case = TRUE)
  x <- sub("_trochanter.*$", "", x, ignore.case = TRUE)
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

strip_leading_number <- function(x) {
  gsub("^[0-9]+_", "", x, perl = TRUE)
}

fmt_p <- function(x) {
  ifelse(is.na(x), NA_character_, format.pval(x, digits = 3, eps = 0.001))
}

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

read_cs_from_txt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  i_raw <- which(trimws(lines) == "[rawpoints]")
  if (length(i_raw) == 0) {
    stop_with_hint(paste0("No [rawpoints] in: ", basename(path)))
  }
  
  after <- lines[(i_raw[1] + 1):length(lines)]
  after <- after[nchar(trimws(after)) > 0]
  
  if (length(after) > 0 && grepl("^'?\\#?1\\s*$", trimws(after[1]))) {
    after <- after[-1]
  }
  
  is_xyz <- grepl(
    "^[[:space:]]*[-+0-9\\.eE]+[[:space:]]+[-+0-9\\.eE]+[[:space:]]+[-+0-9\\.eE]+[[:space:]]*$",
    after
  )
  xyz_lines <- after[is_xyz]
  
  if (length(xyz_lines) < 3) {
    stop_with_hint(paste0("Not enough coordinate lines in: ", basename(path)))
  }
  
  xyz_chr <- do.call(rbind, strsplit(trimws(xyz_lines), "\\s+"))
  xyz <- matrix(as.numeric(xyz_chr), ncol = 3, byrow = TRUE)
  
  centroid <- colMeans(xyz)
  sqrt(sum(rowSums((xyz - matrix(centroid, nrow(xyz), 3, byrow = TRUE))^2)))
}

# ------------------- 1) CENTROID SIZE TABLE -------------------
lm_files <- list.files(lm_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(lm_files) == 0) {
  stop_with_hint(paste0("No .txt landmark files found in: ", lm_dir))
}

lm_names  <- basename(lm_files)
spec_norm <- vapply(lm_names, norm_id, FUN.VALUE = character(1))
spec_norm <- as.character(unlist(spec_norm, use.names = FALSE))

cs_tbl <- tibble(
  lm_file = lm_files,
  lm_name = lm_names,
  specimen_id_norm = spec_norm,
  centroid_size = purrr::map_dbl(lm_files, read_cs_from_txt)
) %>%
  dplyr::mutate(
    specimen_id_norm   = as.character(specimen_id_norm),
    specimen_name_norm = strip_leading_number(specimen_id_norm)
  )

if (is.list(cs_tbl$specimen_id_norm) || is.list(cs_tbl$specimen_name_norm)) {
  stop_with_hint("specimen_id_norm/specimen_name_norm ist als LIST gespeichert. Bitte poste: str(cs_tbl)")
}

dup_full <- cs_tbl$specimen_id_norm[duplicated(cs_tbl$specimen_id_norm)]
dup_name <- cs_tbl$specimen_name_norm[duplicated(cs_tbl$specimen_name_norm)]

if (length(dup_full) > 0) {
  message(" Duplicate specimen_id_norm in landmark table (check filenames):")
  print(sort(table(cs_tbl$specimen_id_norm), decreasing = TRUE)[1:min(10, length(table(cs_tbl$specimen_id_norm)))])
}
if (length(dup_name) > 0) {
  message(" Duplicate specimen_name_norm in landmark table (name-only key not unique):")
  print(sort(table(cs_tbl$specimen_name_norm), decreasing = TRUE)[1:min(10, length(table(cs_tbl$specimen_name_norm)))])
  message("If name-only duplicates exist, fallback matching may be ambiguous for those taxa.")
}

# ------------------- 2) READ KEY + PCA -------------------
key <- read.csv2(key_path, stringsAsFactors = FALSE, check.names = FALSE)
pca <- read.csv2(pca_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!("specimen_id" %in% names(key))) {
  stop_with_hint("specimen_key.csv has no column named 'specimen_id'.")
}
if (!("specimen_id" %in% names(pca))) {
  stop_with_hint("PCA CSV has no column named 'specimen_id'.")
}

key2 <- key %>%
  dplyr::mutate(
    specimen_id_norm   = as.character(vapply(specimen_id, norm_id, FUN.VALUE = character(1))),
    specimen_name_norm = strip_leading_number(specimen_id_norm)
  )

pca2 <- pca %>%
  dplyr::mutate(
    specimen_id_norm   = as.character(vapply(specimen_id, norm_id, FUN.VALUE = character(1))),
    specimen_name_norm = strip_leading_number(specimen_id_norm)
  )

# ------------------- 3) MERGE centroid_size -------------------
key2 <- key2 %>%
  dplyr::left_join(
    cs_tbl %>% dplyr::select(specimen_id_norm, centroid_size),
    by = "specimen_id_norm"
  ) %>%
  dplyr::left_join(
    cs_tbl %>% dplyr::select(specimen_name_norm, centroid_size_name = centroid_size),
    by = "specimen_name_norm"
  ) %>%
  dplyr::mutate(centroid_size = dplyr::coalesce(centroid_size, centroid_size_name)) %>%
  dplyr::select(-centroid_size_name)

pca2 <- pca2 %>%
  dplyr::left_join(
    cs_tbl %>% dplyr::select(specimen_id_norm, centroid_size),
    by = "specimen_id_norm"
  ) %>%
  dplyr::left_join(
    cs_tbl %>% dplyr::select(specimen_name_norm, centroid_size_name = centroid_size),
    by = "specimen_name_norm"
  ) %>%
  dplyr::mutate(centroid_size = dplyr::coalesce(centroid_size, centroid_size_name)) %>%
  dplyr::select(-centroid_size_name)

# missing report
missing_key <- key2 %>% dplyr::filter(is.na(centroid_size)) %>% dplyr::pull(specimen_id)
missing_pca <- pca2 %>% dplyr::filter(is.na(centroid_size)) %>% dplyr::pull(specimen_id)

cat("\n--- Missing centroid_size after merge ---\n")
cat("key missing:", length(missing_key), "\n")
if (length(missing_key) > 0) print(missing_key)
cat("pca missing:", length(missing_pca), "\n")
if (length(missing_pca) > 0) print(missing_pca)

# ------------------- 4) WRITE merged outputs -------------------
write.csv2(
  key2 %>% dplyr::select(-specimen_id_norm, -specimen_name_norm),
  out_key_csv,
  row.names = FALSE
)

write.csv2(
  pca2 %>% dplyr::select(-specimen_id_norm, -specimen_name_norm),
  out_pca_csv,
  row.names = FALSE
)

cat("\n Wrote merged files:\n", out_key_csv, "\n", out_pca_csv, "\n")
cat("\n--- Centroid size summary (PCA table) ---\n")
print(summary(pca2$centroid_size))

# ============================================================
# 5) ALLOMETRY TESTS
# ============================================================

df <- pca2

df$centroid_size <- suppressWarnings(as.numeric(df$centroid_size))
if (any(is.na(df$centroid_size))) {
  stop_with_hint("centroid_size contains NA after numeric conversion.")
}
if (any(df$centroid_size <= 0)) {
  stop_with_hint("centroid_size must be > 0 for log().")
}

pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
if (length(pc_cols) < 5) {
  stop_with_hint(paste0("Found fewer than 5 PC columns. Detected: ", paste(pc_cols, collapse = ", ")))
}

missing_required_pcs <- setdiff(pcs_use, names(df))
if (length(missing_required_pcs) > 0) {
  stop_with_hint(paste0("Missing required PCs: ", paste(missing_required_pcs, collapse = ", ")))
}

for (cc in pcs_use) {
  df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
}

df0 <- df %>%
  dplyr::mutate(logCS = log(centroid_size)) %>%
  dplyr::filter(!is.na(logCS)) %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(pcs_use), ~ !is.na(.x)))

if (nrow(df0) < 10) {
  stop_with_hint(paste0("Too few rows after filtering NA. Remaining: ", nrow(df0)))
}

cat("\n=== Allometry data ===\n")
cat("Rows:", nrow(df0), "\n")
cat("Using PCs:", paste(pcs_use, collapse = ", "), "\n")
cat("Permutations:", n_permutations, "\n")
cat("P-value adjustment for univariate tests:", p_adjust_method, "\n\n")

# ============================================================
# 5.1 UNIVARIATE TESTS: PC1-PC5 ~ logCS
# ============================================================

cat("=== Univariat: PC1-PC5 ~ logCS ===\n")

univar_results <- lapply(pcs_use, function(pc) {
  form <- stats::as.formula(paste(pc, "~ logCS"))
  fit  <- stats::lm(form, data = df0)
  sm   <- summary(fit)
  
  coef_tab <- coef(sm)
  slope_row <- rownames(coef_tab) == "logCS"
  
  tibble(
    trait = pc,
    estimate = coef_tab[slope_row, "Estimate"],
    std_error = coef_tab[slope_row, "Std. Error"],
    t_value = coef_tab[slope_row, "t value"],
    p_value = coef_tab[slope_row, "Pr(>|t|)"],
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    f_statistic = unname(sm$fstatistic[1]),
    df_model = unname(sm$fstatistic[2]),
    df_residual = unname(sm$fstatistic[3]),
    n = nrow(df0)
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    p_value_adjusted = stats::p.adjust(p_value, method = p_adjust_method),
    significance_raw = sig_label(p_value),
    significance_adjusted = sig_label(p_value_adjusted),
    p_adjust_method = p_adjust_method
  )

print(univar_results)

write.table(
  univar_results,
  file = out_univar_csv,
  sep = ";",
  dec = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

# plot: univariate panels
plot_list <- lapply(pcs_use, function(pc) {
  row_i <- univar_results %>% dplyr::filter(trait == pc)
  
  ggplot(df0, aes_string(x = "logCS", y = pc)) +
    geom_point(
      size = 2.4,
      shape = 21,
      stroke = 0.5,
      fill = "white",
      color = "black",
      alpha = 0.9
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 1.0,
      color = "#1f1f1f"
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      hjust = 1.05, vjust = 1.35,
      size = 3.4,
      label = paste0(
        "R = ", format(round(row_i$r_squared, 3), nsmall = 3), "\n",
        "p = ", fmt_p(row_i$p_value), "\n",
        "Holm p = ", fmt_p(row_i$p_value_adjusted)
      )
    ) +
    labs(
      title = paste0(pc, " ~ log(Centroid size)"),
      x = "log(Centroid size)",
      y = paste0(pc, " score")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
})

p_uni_grid <- patchwork::wrap_plots(plotlist = plot_list, ncol = 2) +
  patchwork::plot_annotation(
    title = "Univariate allometry tests for PC1-PC5",
    subtitle = paste0("Multiple-testing correction: ", tools::toTitleCase(p_adjust_method))
  )

print(p_uni_grid)

ggsave(
  filename = out_plot_uni_png,
  plot = p_uni_grid,
  width = 12, height = 14, dpi = 400
)

ggsave(
  filename = out_plot_uni_pdf,
  plot = p_uni_grid,
  width = 12, height = 14
)

# ============================================================
# 5.2 MULTIVARIATE (RRPP): (PC1..PC5) ~ logCS
# ============================================================

cat("=== Multivariat (RRPP): (PC1..PC5) ~ logCS ===\n")
Y <- as.matrix(df0[, pcs_use, drop = FALSE])

fit_rrpp <- RRPP::lm.rrpp(
  Y ~ logCS,
  data = df0,
  iter = n_permutations,
  print.progress = FALSE
)

rrpp_anova <- anova(fit_rrpp)
rrpp_summary <- summary(fit_rrpp)

cat("\n--- ANOVA (RRPP) ---\n")
print(rrpp_anova)

cat("\n--- Summary (RRPP) ---\n")
print(rrpp_summary)

# robust RRPP table extraction
rrpp_tab <- NULL

if (!is.null(rrpp_anova$table)) {
  rrpp_tab <- rrpp_anova$table
} else if (!is.null(rrpp_anova$ANOVA)) {
  rrpp_tab <- rrpp_anova$ANOVA
} else if (inherits(rrpp_anova, "data.frame")) {
  rrpp_tab <- rrpp_anova
} else if (inherits(rrpp_anova, "matrix")) {
  rrpp_tab <- as.data.frame(rrpp_anova)
}

if (is.null(rrpp_tab)) {
  stop_with_hint("Could not extract RRPP ANOVA table. Run: str(rrpp_anova) and paste the output.")
}

rrpp_tab <- as.data.frame(rrpp_tab)
rrpp_tab$Term <- rownames(rrpp_tab)
rownames(rrpp_tab) <- NULL
rrpp_tab <- rrpp_tab[, c("Term", setdiff(names(rrpp_tab), "Term")), drop = FALSE]

write.table(
  rrpp_tab,
  file = out_rrpp_csv,
  sep = ";",
  dec = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n Wrote RRPP multivariate table:\n", out_rrpp_csv, "\n")
# ============================================================
# 5.2b MORPHOSPACE PLOT: PC1-PC2 colored by logCS
# ============================================================

p_morph <- ggplot(df0, aes(x = PC1, y = PC2, color = logCS)) +
  geom_point(size = 3, alpha = 0.95) +
  scale_color_viridis_c(
    option = "plasma",
    end = 0.95,
    name = "log(CS)"
  ) +
  labs(
    title = "Morphospace PC1-PC2 colored by log(Centroid size)",
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black")
  )

print(p_morph)

ggsave(
  filename = out_plot_morph_png,
  plot = p_morph,
  width = 7.2, height = 5.4, dpi = 400
)

ggsave(
  filename = out_plot_morph_pdf,
  plot = p_morph,
  width = 7.2, height = 5.4
)

# ============================================================
# 5.2c COMBINED PLOT
# ============================================================

p_combined <- p_morph / p_uni_grid +
  patchwork::plot_annotation(tag_levels = "A")

print(p_combined)

ggsave(
  filename = out_plot_combined_png,
  plot = p_combined,
  width = 12, height = 18, dpi = 400
)

ggsave(
  filename = out_plot_combined_pdf,
  plot = p_combined,
  width = 12, height = 18
)

# ============================================================
# 5.3 geomorph procD.lm
# ============================================================

cat("\n=== geomorph: procD.lm on PC1-PC5 ===\n")
fit_allo <- geomorph::procD.lm(
  Y ~ logCS,
  data = df0,
  iter = n_permutations
)

fit_allo_summary <- summary(fit_allo)
print(fit_allo_summary)

aov_obj <- anova(fit_allo)

tab <- NULL
if (!is.null(aov_obj$ANOVA)) tab <- aov_obj$ANOVA
if (is.null(tab) && !is.null(aov_obj$aov.table)) tab <- aov_obj$aov.table
if (is.null(tab) && !is.null(aov_obj$table)) tab <- aov_obj$table

if (is.null(tab)) {
  stop_with_hint("Could not find ANOVA table inside aov_obj. Run: str(aov_obj) and paste the output.")
}

procD_tab <- as.data.frame(tab)
procD_tab$Term <- rownames(procD_tab)
rownames(procD_tab) <- NULL

wanted <- c("Term", "Df", "SS", "MS", "Rsq", "F", "Z", "Pr(>F)")
have <- intersect(wanted, names(procD_tab))
procD_tab <- procD_tab[, c(have, setdiff(names(procD_tab), have)), drop = FALSE]

write.table(
  procD_tab,
  file = out_procD_csv,
  sep = ";",
  dec = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n Wrote procD.lm multivariate table:\n", out_procD_csv, "\n")

# ============================================================
# FINAL MESSAGE
# ============================================================

cat("\n========================================\n")
cat("DONE. All outputs were written to:\n")
cat(out_dir, "\n")
cat("========================================\n")
