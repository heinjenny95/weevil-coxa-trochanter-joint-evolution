# ============================================================
# EXTENDED ALLOMETRY PIPELINE
# FINAL DEBUG VERSION WITH HEADER-RECOVERY FOR GEOMETRY
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(RRPP)
  library(geomorph)
  library(viridis)
  library(patchwork)
  library(nnet)
  library(data.table)
})

# ------------------- PATHS -------------------
lm_dir   <- "<BEETLE_JOINTS_ROOT>/Processed/Curculionoidea/Landmarks"
key_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_key.csv"
pca_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"
geom_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/winding_metrics_excel.csv"

base_results_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results"
out_dir <- file.path(base_results_dir, "Allometry")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------- OUTPUT FILES -------------------
out_key_csv                 <- file.path(out_dir, "specimen_key_with_centroid_size.csv")
out_pca_csv                 <- file.path(out_dir, "PCA_scores_with_specimen_id_with_centroid_size.csv")
out_merged_csv              <- file.path(out_dir, "allometry_merged_table.csv")

out_univar_pc_csv           <- file.path(out_dir, "allometry_univariate_PC1_to_PC5_results.csv")
out_rrpp_csv                <- file.path(out_dir, "allometry_rrpp_multivariate_results.csv")
out_procD_csv               <- file.path(out_dir, "allometry_procD_lm_results.csv")
out_cont_csv                <- file.path(out_dir, "allometry_continuous_traits_results.csv")
out_group_csv               <- file.path(out_dir, "allometry_group_size_tests_results.csv")
out_binom_csv               <- file.path(out_dir, "allometry_binary_trait_glm_results.csv")
out_multinom_csv            <- file.path(out_dir, "allometry_multinomial_joint_type_results.csv")

out_manifest_txt            <- file.path(out_dir, "allometry_manifest_detected_columns.txt")
out_log_txt                 <- file.path(out_dir, "allometry_log.txt")
out_geom_match_csv          <- file.path(out_dir, "geometry_matching_diagnostics.csv")
out_geom_unmatched_pca_csv  <- file.path(out_dir, "geometry_unmatched_pca_ids.csv")
out_geom_unmatched_geom_csv <- file.path(out_dir, "geometry_unmatched_geometry_ids.csv")

out_plot_uni_png            <- file.path(out_dir, "allometry_univariate_PC1_to_PC5.png")
out_plot_uni_pdf            <- file.path(out_dir, "allometry_univariate_PC1_to_PC5.pdf")
out_plot_morph_png          <- file.path(out_dir, "morphospace_PC1_PC2_logCS.png")
out_plot_morph_pdf          <- file.path(out_dir, "morphospace_PC1_PC2_logCS.pdf")
out_plot_summary_png        <- file.path(out_dir, "allometry_summary_effects.png")
out_plot_summary_pdf        <- file.path(out_dir, "allometry_summary_effects.pdf")
out_plot_geom_png           <- file.path(out_dir, "allometry_geometry_selected_panels.png")
out_plot_geom_pdf           <- file.path(out_dir, "allometry_geometry_selected_panels.pdf")

# ------------------- SETTINGS -------------------
pcs_use                <- paste0("PC", 1:5)
n_permutations         <- 9999
p_adjust_method        <- "holm"
min_n_for_group_test   <- 3
max_panels_geometry    <- 6
epsilon_ratio          <- 1e-12

# ------------------- LOGGING -------------------
log_lines <- character()

log_msg <- function(...) {
  txt <- paste0(...)
  message(txt)
  log_lines <<- c(log_lines, txt)
}

# ------------------- HELPERS -------------------
stop_with_hint <- function(msg) {
  stop(paste0("\n ", msg, "\n"), call. = FALSE)
}

norm_id <- function(x) {
  x <- as.character(x)
  x <- enc2utf8(x)
  x <- trimws(x)
  x <- sub("^\ufeff", "", x)
  x <- tolower(x)
  x <- sub("\\.txt$", "", x, ignore.case = TRUE)
  x <- sub("\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("\\.vtk$", "", x, ignore.case = TRUE)
  x <- sub("_trochanter.*$", "", x, ignore.case = TRUE)
  x <- gsub("_aligned$", "", x, ignore.case = TRUE)
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

strip_leading_number <- function(x) {
  gsub("^[0-9]+_", "", x, perl = TRUE)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

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

safe_bind_rows <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (length(lst) == 0) return(tibble::tibble())
  dplyr::bind_rows(lst)
}

looks_like_default_vnames <- function(nm) {
  all(grepl("^V[0-9]+$", nm))
}

recover_header_from_first_row <- function(df, filename = "") {
  if (nrow(df) < 2) return(df)
  if (!looks_like_default_vnames(names(df))) return(df)
  
  first_row <- as.character(unlist(df[1, , drop = TRUE]))
  first_row <- enc2utf8(first_row)
  first_row <- trimws(first_row)
  
  header_keywords <- c(
    "specimen_id",
    "signed_winding_angle_deg",
    "abs_winding_angle_deg",
    "n_turns_signed",
    "n_turns_abs",
    "start_end_dist",
    "axial_span",
    "fit_radius",
    "fit_rms"
  )
  
  keyword_hits <- sum(first_row %in% header_keywords)
  
  if (keyword_hits >= 3) {
    names(df) <- first_row
    df <- df[-1, , drop = FALSE]
    rownames(df) <- NULL
    log_msg(" Recovered header from first data row for ", filename, ".")
  }
  
  df
}

read_fread_robust <- function(path) {
  if (!file.exists(path)) stop_with_hint(paste0("File does not exist: ", path))
  
  # plain fread first, like your successful manual test
  df <- try(
    data.table::fread(path, data.table = FALSE),
    silent = TRUE
  )
  
  if (inherits(df, "try-error") || !is.data.frame(df) || ncol(df) < 2) {
    # fallback with a few extra args
    df <- try(
      data.table::fread(
        path,
        data.table = FALSE,
        encoding = "UTF-8",
        fill = TRUE,
        quote = "\""
      ),
      silent = TRUE
    )
  }
  
  if (inherits(df, "try-error") || !is.data.frame(df) || ncol(df) < 2) {
    stop_with_hint(paste0("Could not read file with fread: ", path))
  }
  
  names(df) <- enc2utf8(names(df))
  names(df) <- sub("^\ufeff", "", names(df))
  names(df) <- trimws(names(df))
  
  df <- recover_header_from_first_row(df, basename(path))
  
  names(df) <- enc2utf8(names(df))
  names(df) <- sub("^\ufeff", "", names(df))
  names(df) <- trimws(names(df))
  
  log_msg(" fread read: ", basename(path), " | ncol=", ncol(df), " nrow=", nrow(df))
  log_msg("   columns: ", paste(names(df), collapse = " | "))
  
  df
}

force_geometry_id_col <- function(df) {
  nm <- names(df)
  nm_clean <- tolower(trimws(nm))
  nm_clean <- sub("^\ufeff", "", nm_clean)
  
  if ("specimen_id" %in% nm_clean) {
    names(df)[match("specimen_id", nm_clean)] <- "specimen_id"
    return(df)
  }
  
  alias_idx <- which(nm_clean %in% c(
    "specimen id", "specimen", "id", "sample_id", "sample id",
    "filename", "file_name", "name"
  ))
  if (length(alias_idx) > 0) {
    old_name <- nm[alias_idx[1]]
    names(df)[alias_idx[1]] <- "specimen_id"
    log_msg(" Geometry ID column was not exact. Renamed '", old_name, "' to 'specimen_id'.")
    return(df)
  }
  
  old_name <- nm[1]
  names(df)[1] <- "specimen_id"
  log_msg(" No explicit geometry ID column found. Forced first column ('", old_name, "') to 'specimen_id'.")
  df
}

rename_first_match <- function(df, candidates, new_name) {
  hit <- intersect(candidates, names(df))
  if (length(hit) > 0) {
    names(df)[names(df) == hit[1]] <- new_name
  }
  df
}

canonicalize_geometry_cols <- function(df) {
  df <- rename_first_match(df, c("signed_winding_angle_deg", "winding_angle_deg", "signed_winding_angle"), "signed_winding_angle_deg")
  df <- rename_first_match(df, c("abs_winding_angle_deg", "absolute_winding_angle", "winding_angle_abs"), "abs_winding_angle_deg")
  df <- rename_first_match(df, c("n_turns_signed", "turns_signed"), "n_turns_signed")
  df <- rename_first_match(df, c("n_turns_abs", "turns_abs", "absolute_turn_number"), "n_turns_abs")
  df <- rename_first_match(df, c("start_end_dist", "start_end_distance"), "start_end_dist")
  df <- rename_first_match(df, c("axial_span"), "axial_span")
  df <- rename_first_match(df, c("fit_radius", "fitted_radius", "radius"), "fit_radius")
  df <- rename_first_match(df, c("fit_rms", "rms", "fit_error"), "fit_rms")
  df
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

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

run_lm_table <- function(df, response, predictor = "logCS", log_response = FALSE) {
  tmp <- df %>%
    dplyr::select(dplyr::all_of(c(response, predictor))) %>%
    dplyr::mutate(
      response_raw = .data[[response]],
      predictor_raw = .data[[predictor]]
    )
  
  if (log_response) {
    tmp <- tmp %>% dplyr::mutate(response_used = ifelse(response_raw > 0, log(response_raw), NA_real_))
  } else {
    tmp <- tmp %>% dplyr::mutate(response_used = response_raw)
  }
  
  tmp <- tmp %>%
    dplyr::mutate(predictor_used = predictor_raw) %>%
    dplyr::filter(!is.na(response_used), !is.na(predictor_used))
  
  if (nrow(tmp) < 10) {
    return(tibble::tibble(
      trait = response,
      model_type = ifelse(log_response, "lm_log_response", "lm"),
      n = nrow(tmp),
      estimate = NA_real_,
      std_error = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      f_statistic = NA_real_,
      df_model = NA_real_,
      df_residual = NA_real_,
      response_transform = ifelse(log_response, "log", "none")
    ))
  }
  
  fit <- stats::lm(response_used ~ predictor_used, data = tmp)
  sm  <- summary(fit)
  cf  <- coef(sm)
  
  tibble::tibble(
    trait = response,
    model_type = ifelse(log_response, "lm_log_response", "lm"),
    n = nrow(tmp),
    estimate = cf["predictor_used", "Estimate"],
    std_error = cf["predictor_used", "Std. Error"],
    t_value = cf["predictor_used", "t value"],
    p_value = cf["predictor_used", "Pr(>|t|)"],
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    f_statistic = unname(sm$fstatistic[1]),
    df_model = unname(sm$fstatistic[2]),
    df_residual = unname(sm$fstatistic[3]),
    response_transform = ifelse(log_response, "log", "none")
  )
}

run_group_size_test <- function(df, group_col, size_col = "logCS", min_n = 3) {
  tmp <- df %>%
    dplyr::select(dplyr::all_of(c(group_col, size_col))) %>%
    dplyr::rename(group = dplyr::all_of(group_col), size = dplyr::all_of(size_col)) %>%
    dplyr::filter(!is.na(group), !is.na(size)) %>%
    dplyr::mutate(group = as.factor(group))
  
  group_sizes <- table(tmp$group)
  keep_levels <- names(group_sizes[group_sizes >= min_n])
  tmp <- tmp %>% dplyr::filter(group %in% keep_levels) %>% droplevels()
  
  if (nlevels(tmp$group) < 2 || nrow(tmp) < 10) {
    return(tibble::tibble(
      grouping = group_col,
      test = NA_character_,
      n = nrow(tmp),
      n_levels = nlevels(tmp$group),
      statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_
    ))
  }
  
  if (nlevels(tmp$group) == 2) {
    wt <- stats::wilcox.test(size ~ group, data = tmp, exact = FALSE)
    tibble::tibble(
      grouping = group_col,
      test = "wilcox_logCS_by_group",
      n = nrow(tmp),
      n_levels = nlevels(tmp$group),
      statistic = unname(wt$statistic),
      df = NA_real_,
      p_value = wt$p.value
    )
  } else {
    kt <- stats::kruskal.test(size ~ group, data = tmp)
    tibble::tibble(
      grouping = group_col,
      test = "kruskal_logCS_by_group",
      n = nrow(tmp),
      n_levels = nlevels(tmp$group),
      statistic = unname(kt$statistic),
      df = unname(kt$parameter),
      p_value = kt$p.value
    )
  }
}

run_binomial_glm <- function(df, trait_col, predictor = "logCS") {
  tmp <- df %>%
    dplyr::select(dplyr::all_of(c(trait_col, predictor))) %>%
    dplyr::rename(y = dplyr::all_of(trait_col), x = dplyr::all_of(predictor)) %>%
    dplyr::filter(!is.na(y), !is.na(x))
  
  if (is.logical(tmp$y)) {
    tmp$y <- as.integer(tmp$y)
  } else if (is.numeric(tmp$y)) {
    uniq <- sort(unique(tmp$y))
    if (!all(uniq %in% c(0, 1))) return(NULL)
  } else {
    y_chr <- tolower(trimws(as.character(tmp$y)))
    map <- c("true"=1,"false"=0,"yes"=1,"no"=0,"present"=1,"absent"=0)
    if (!all(y_chr %in% names(map))) return(NULL)
    tmp$y <- unname(map[y_chr])
  }
  
  if (length(unique(tmp$y)) < 2 || nrow(tmp) < 10) return(NULL)
  
  fit <- try(stats::glm(y ~ x, data = tmp, family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  sm <- summary(fit)
  cf <- coef(sm)
  
  tibble::tibble(
    trait = trait_col,
    n = nrow(tmp),
    estimate = cf["x", "Estimate"],
    std_error = cf["x", "Std. Error"],
    z_value = cf["x", "z value"],
    p_value = cf["x", "Pr(>|z|)"],
    odds_ratio = exp(cf["x", "Estimate"])
  )
}

run_multinom_joint_type <- function(df, trait_col, predictor = "logCS", min_n = 3) {
  tmp <- df %>%
    dplyr::select(dplyr::all_of(c(trait_col, predictor))) %>%
    dplyr::rename(y = dplyr::all_of(trait_col), x = dplyr::all_of(predictor)) %>%
    dplyr::filter(!is.na(y), !is.na(x)) %>%
    dplyr::mutate(y = as.factor(y))
  
  ytab <- table(tmp$y)
  keep_levels <- names(ytab[ytab >= min_n])
  tmp <- tmp %>% dplyr::filter(y %in% keep_levels) %>% droplevels()
  
  if (nlevels(tmp$y) < 3 || nrow(tmp) < 15) return(NULL)
  
  fit <- try(nnet::multinom(y ~ x, data = tmp, trace = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  sm <- summary(fit)
  est <- sm$coefficients
  se  <- sm$standard.errors
  
  if (is.null(dim(est))) {
    z <- est["x"] / se["x"]
    p <- 2 * (1 - pnorm(abs(z)))
    out <- tibble::tibble(
      trait = trait_col,
      comparison = rownames(as.matrix(est))[1],
      n = nrow(tmp),
      estimate = unname(est["x"]),
      std_error = unname(se["x"]),
      z_value = unname(z),
      p_value = unname(p)
    )
  } else {
    out <- tibble::tibble(
      trait = trait_col,
      comparison = rownames(est),
      n = nrow(tmp),
      estimate = est[, "x"],
      std_error = se[, "x"],
      z_value = est[, "x"] / se[, "x"],
      p_value = 2 * (1 - pnorm(abs(est[, "x"] / se[, "x"])))
    )
  }
  
  out
}

# ============================================================
# 1) CENTROID SIZE TABLE
# ============================================================
lm_files <- list.files(lm_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(lm_files) == 0) stop_with_hint(paste0("No .txt landmark files found in: ", lm_dir))

lm_names  <- basename(lm_files)
spec_norm <- vapply(lm_names, norm_id, FUN.VALUE = character(1))

cs_tbl <- tibble::tibble(
  lm_file = lm_files,
  lm_name = lm_names,
  specimen_id_norm = as.character(spec_norm),
  centroid_size = purrr::map_dbl(lm_files, read_cs_from_txt)
) %>%
  dplyr::mutate(
    specimen_name_norm = strip_leading_number(specimen_id_norm)
  )

# ============================================================
# 2) READ INPUT TABLES
# ============================================================
key <- read_fread_robust(key_path)
pca <- read_fread_robust(pca_path)
geom_df <- read_fread_robust(geom_path)
geom_df <- force_geometry_id_col(geom_df)
geom_df <- canonicalize_geometry_cols(geom_df)

if (!("specimen_id" %in% names(key)))     stop_with_hint("specimen_key.csv has no column 'specimen_id'.")
if (!("specimen_id" %in% names(pca)))     stop_with_hint("PCA file has no column 'specimen_id'.")
if (!("specimen_id" %in% names(geom_df))) stop_with_hint("Geometry file still has no usable ID column after forced handling.")

# ============================================================
# 3) NORMALIZE IDS
# ============================================================
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

geom2 <- geom_df %>%
  dplyr::mutate(
    specimen_id_norm   = as.character(vapply(specimen_id, norm_id, FUN.VALUE = character(1))),
    specimen_name_norm = strip_leading_number(specimen_id_norm)
  )

# ============================================================
# 4) MERGE CENTROID SIZE
# ============================================================
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

# ============================================================
# 5) MERGE KEY TRAITS INTO PCA
# ============================================================
key_trait_candidates <- c(
  "Family", "family",
  "Coxal wall hole", "coxal_wall_hole",
  "Coxal Socket", "coxal_socket",
  "Windung Coxa", "windung_coxa",
  "joint_type", "screw_state",
  "Schraube", "screw"
)

key_traits_present <- intersect(key_trait_candidates, names(key2))

key_traits_tbl <- key2 %>%
  dplyr::select(specimen_id_norm, specimen_name_norm, dplyr::all_of(key_traits_present)) %>%
  dplyr::distinct()

pca2 <- pca2 %>%
  dplyr::left_join(key_traits_tbl, by = c("specimen_id_norm", "specimen_name_norm"))

# ============================================================
# 6) GEOMETRY MATCHING DIAGNOSTICS
# ============================================================
geom_trait_candidates <- c(
  "signed_winding_angle_deg", "abs_winding_angle_deg",
  "n_turns_signed", "n_turns_abs",
  "start_end_dist", "axial_span", "fit_radius", "fit_rms"
)
geom_traits_present <- intersect(geom_trait_candidates, names(geom2))

full_id_matches <- sum(pca2$specimen_id_norm %in% geom2$specimen_id_norm)
name_matches    <- sum(pca2$specimen_name_norm %in% geom2$specimen_name_norm)

log_msg("Geometry traits present in geometry table: ", paste(geom_traits_present, collapse = ", "))
log_msg("PCA rows matching geometry by full normalized ID: ", full_id_matches, " / ", nrow(pca2))
log_msg("PCA rows matching geometry by stripped name ID: ", name_matches, " / ", nrow(pca2))

geom_match_diag <- pca2 %>%
  dplyr::select(specimen_id, specimen_id_norm, specimen_name_norm) %>%
  dplyr::mutate(
    match_full_id = specimen_id_norm %in% geom2$specimen_id_norm,
    match_name    = specimen_name_norm %in% geom2$specimen_name_norm
  )

write.csv2(geom_match_diag, out_geom_match_csv, row.names = FALSE)

unmatched_pca <- geom_match_diag %>%
  dplyr::filter(!match_full_id & !match_name)
write.csv2(unmatched_pca, out_geom_unmatched_pca_csv, row.names = FALSE)

unmatched_geom <- geom2 %>%
  dplyr::select(specimen_id, specimen_id_norm, specimen_name_norm) %>%
  dplyr::mutate(
    used_full = specimen_id_norm %in% pca2$specimen_id_norm,
    used_name = specimen_name_norm %in% pca2$specimen_name_norm
  ) %>%
  dplyr::filter(!used_full & !used_name)
write.csv2(unmatched_geom, out_geom_unmatched_geom_csv, row.names = FALSE)

# ============================================================
# 7) MERGE GEOMETRY INTO PCA WITH FULL-ID + NAME FALLBACK
# ============================================================
geom_full_tbl <- geom2 %>%
  dplyr::select(specimen_id_norm, dplyr::all_of(geom_traits_present)) %>%
  dplyr::distinct()

geom_name_tbl <- geom2 %>%
  dplyr::select(specimen_name_norm, dplyr::all_of(geom_traits_present)) %>%
  dplyr::distinct()

pca2 <- pca2 %>%
  dplyr::left_join(
    geom_full_tbl,
    by = "specimen_id_norm"
  )

geom_name_tbl_renamed <- geom_name_tbl
names(geom_name_tbl_renamed)[names(geom_name_tbl_renamed) %in% geom_traits_present] <- paste0(geom_traits_present, "_name")

pca2 <- pca2 %>%
  dplyr::left_join(
    geom_name_tbl_renamed,
    by = "specimen_name_norm"
  )

for (gc in geom_traits_present) {
  gc_name <- paste0(gc, "_name")
  if (gc_name %in% names(pca2)) {
    pca2[[gc]] <- dplyr::coalesce(pca2[[gc]], pca2[[gc_name]])
  }
}

pca2 <- pca2 %>%
  dplyr::select(-dplyr::any_of(paste0(geom_traits_present, "_name")))

geom_non_na_counts <- sapply(geom_traits_present, function(v) sum(!is.na(pca2[[v]])))
for (v in names(geom_non_na_counts)) {
  log_msg("Non-NA values after geometry merge for ", v, ": ", geom_non_na_counts[[v]])
}

# ============================================================
# 8) WRITE MERGED TABLES
# ============================================================
write.csv2(key2 %>% dplyr::select(-specimen_id_norm, -specimen_name_norm), out_key_csv, row.names = FALSE)
write.csv2(pca2 %>% dplyr::select(-specimen_id_norm, -specimen_name_norm), out_pca_csv, row.names = FALSE)

# ============================================================
# 9) PREP MERGED DATA
# ============================================================
df <- pca2

df$centroid_size <- safe_numeric(df$centroid_size)
if (any(is.na(df$centroid_size))) stop_with_hint("centroid_size contains NA after numeric conversion.")
if (any(df$centroid_size <= 0)) stop_with_hint("centroid_size must be > 0 for log().")

missing_required_pcs <- setdiff(pcs_use, names(df))
if (length(missing_required_pcs) > 0) stop_with_hint(paste0("Missing required PCs: ", paste(missing_required_pcs, collapse = ", ")))

for (cc in pcs_use) df[[cc]] <- safe_numeric(df[[cc]])

geom_aliases <- list(
  abs_winding_angle_deg    = c("abs_winding_angle_deg"),
  signed_winding_angle_deg = c("signed_winding_angle_deg"),
  axial_span               = c("axial_span"),
  start_end_dist           = c("start_end_dist"),
  fit_radius               = c("fit_radius"),
  n_turns_abs              = c("n_turns_abs"),
  fit_rms                  = c("fit_rms")
)

detected_geom_cols <- purrr::map_chr(geom_aliases, ~ find_first_col(df, .x))
detected_geom_cols <- detected_geom_cols[!is.na(detected_geom_cols)]

trait_aliases <- list(
  joint_type         = c("joint_type"),
  screw_state        = c("screw_state", "Schraube", "screw"),
  coxal_wall_hole          = c("Coxal wall hole", "coxal_wall_hole"),
  coxal_socket   = c("Coxal Socket", "coxal_socket"),
  windung_coxa       = c("Windung Coxa", "windung_coxa"),
  family             = c("Family", "family")
)

detected_trait_cols <- purrr::map_chr(trait_aliases, ~ find_first_col(df, .x))
detected_trait_cols <- detected_trait_cols[!is.na(detected_trait_cols)]

for (nm in unique(unname(detected_geom_cols))) {
  df[[nm]] <- safe_numeric(df[[nm]])
}

df <- df %>%
  dplyr::mutate(
    logCS = log(centroid_size),
    ratio_axial_span_fit_radius = if (all(c("axial_span","fit_radius") %in% names(df))) axial_span / (fit_radius + epsilon_ratio) else NA_real_,
    ratio_start_end_fit_radius  = if (all(c("start_end_dist","fit_radius") %in% names(df))) start_end_dist / (fit_radius + epsilon_ratio) else NA_real_,
    ratio_fit_rms_fit_radius    = if (all(c("fit_rms","fit_radius") %in% names(df))) fit_rms / (fit_radius + epsilon_ratio) else NA_real_
  )

df_shape <- df %>%
  dplyr::filter(!is.na(logCS)) %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(pcs_use), ~ !is.na(.x)))

if (nrow(df_shape) < 10) stop_with_hint(paste0("Too few rows after filtering NA for shape analyses. Remaining: ", nrow(df_shape)))

write.csv2(df %>% dplyr::select(-specimen_id_norm, -specimen_name_norm), out_merged_csv, row.names = FALSE)

# ============================================================
# 10) MANIFEST
# ============================================================
manifest_lines <- c(
  "=== DETECTED GEOMETRY COLUMNS ===",
  if (length(detected_geom_cols) > 0) paste(names(detected_geom_cols), "->", unname(detected_geom_cols)) else "none",
  "",
  "=== DETECTED DISCRETE TRAIT COLUMNS ===",
  if (length(detected_trait_cols) > 0) paste(names(detected_trait_cols), "->", unname(detected_trait_cols)) else "none",
  "",
  "=== USED PCs ===",
  paste(pcs_use, collapse = ", ")
)
writeLines(manifest_lines, out_manifest_txt)

# ============================================================
# 11) UNIVARIATE PC TESTS
# ============================================================
univar_pc_results <- lapply(pcs_use, function(pc) {
  run_lm_table(df_shape, response = pc, predictor = "logCS", log_response = FALSE)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    p_value_adjusted = stats::p.adjust(p_value, method = p_adjust_method),
    significance_raw = sig_label(p_value),
    significance_adjusted = sig_label(p_value_adjusted),
    p_adjust_method = p_adjust_method,
    trait_group = "shape_PC"
  )

write.table(
  univar_pc_results,
  file = out_univar_pc_csv,
  sep = ";", dec = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 12) MULTIVARIATE SHAPE ALLOMETRY
# ============================================================
Y <- as.matrix(df_shape[, pcs_use, drop = FALSE])

fit_rrpp <- RRPP::lm.rrpp(Y ~ logCS, data = df_shape, iter = n_permutations, print.progress = FALSE)
rrpp_anova <- anova(fit_rrpp)
rrpp_tab <- if (!is.null(rrpp_anova$table)) rrpp_anova$table else as.data.frame(rrpp_anova)
rrpp_tab <- as.data.frame(rrpp_tab)
rrpp_tab$Term <- rownames(rrpp_tab)
rownames(rrpp_tab) <- NULL
rrpp_tab <- rrpp_tab[, c("Term", setdiff(names(rrpp_tab), "Term")), drop = FALSE]

write.table(
  rrpp_tab,
  file = out_rrpp_csv,
  sep = ";", dec = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE,
  fileEncoding = "UTF-8"
)

fit_allo <- geomorph::procD.lm(Y ~ logCS, data = df_shape, iter = n_permutations)
aov_obj <- anova(fit_allo)
tab <- NULL
if (!is.null(aov_obj$ANOVA)) tab <- aov_obj$ANOVA
if (is.null(tab) && !is.null(aov_obj$aov.table)) tab <- aov_obj$aov.table
if (is.null(tab) && !is.null(aov_obj$table)) tab <- aov_obj$table
if (is.null(tab)) stop_with_hint("Could not extract procD.lm ANOVA table.")

procD_tab <- as.data.frame(tab)
procD_tab$Term <- rownames(procD_tab)
rownames(procD_tab) <- NULL
wanted <- c("Term", "Df", "SS", "MS", "Rsq", "F", "Z", "Pr(>F)")
have <- intersect(wanted, names(procD_tab))
procD_tab <- procD_tab[, c(have, setdiff(names(procD_tab), have)), drop = FALSE]

write.table(
  procD_tab,
  file = out_procD_csv,
  sep = ";", dec = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 13) CONTINUOUS GEOMETRY TRAITS ~ logCS
# ============================================================
continuous_traits_raw <- unique(c(
  unname(detected_geom_cols),
  "ratio_axial_span_fit_radius",
  "ratio_start_end_fit_radius",
  "ratio_fit_rms_fit_radius"
))
continuous_traits_raw <- continuous_traits_raw[continuous_traits_raw %in% names(df)]

log_candidate_traits <- unique(c(
  unname(detected_geom_cols[c("axial_span", "start_end_dist", "fit_radius", "fit_rms")]),
  "ratio_axial_span_fit_radius",
  "ratio_start_end_fit_radius",
  "ratio_fit_rms_fit_radius"
))
log_candidate_traits <- log_candidate_traits[!is.na(log_candidate_traits)]
log_candidate_traits <- log_candidate_traits[log_candidate_traits %in% names(df)]

cont_results_none <- safe_bind_rows(lapply(continuous_traits_raw, function(tr) {
  run_lm_table(df, response = tr, predictor = "logCS", log_response = FALSE)
}))

cont_results_log <- safe_bind_rows(lapply(log_candidate_traits, function(tr) {
  run_lm_table(df, response = tr, predictor = "logCS", log_response = TRUE)
}))

cont_results <- dplyr::bind_rows(cont_results_none, cont_results_log)

if (nrow(cont_results) > 0) {
  cont_results <- cont_results %>%
    dplyr::mutate(
      p_value_adjusted = stats::p.adjust(p_value, method = p_adjust_method),
      significance_raw = sig_label(p_value),
      significance_adjusted = sig_label(p_value_adjusted),
      p_adjust_method = p_adjust_method,
      trait_group = dplyr::case_when(
        grepl("^ratio_", trait) ~ "geometry_ratio",
        TRUE ~ "geometry_continuous"
      )
    ) %>%
    dplyr::arrange(trait_group, p_value_adjusted, desc(abs(r_squared)))
} else {
  cont_results <- tibble::tibble()
}

write.table(
  cont_results,
  file = out_cont_csv,
  sep = ";", dec = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 14) GROUP TESTS
# ============================================================
group_test_cols <- unique(unname(detected_trait_cols))
group_test_cols <- group_test_cols[group_test_cols %in% names(df)]

group_results <- safe_bind_rows(lapply(group_test_cols, function(gc) {
  run_group_size_test(df, group_col = gc, size_col = "logCS", min_n = min_n_for_group_test)
}))

if (nrow(group_results) > 0) {
  group_results <- group_results %>%
    dplyr::mutate(
      p_value_adjusted = stats::p.adjust(p_value, method = p_adjust_method),
      significance_raw = sig_label(p_value),
      significance_adjusted = sig_label(p_value_adjusted),
      p_adjust_method = p_adjust_method
    )
}

write.table(
  group_results,
  file = out_group_csv,
  sep = ";", dec = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 15) OPTIONAL BINOMIAL GLMS
# ============================================================
binom_results <- safe_bind_rows(lapply(group_test_cols, function(gc) {
  run_binomial_glm(df, trait_col = gc, predictor = "logCS")
}))

if (nrow(binom_results) > 0) {
  binom_results <- binom_results %>%
    dplyr::mutate(
      p_value_adjusted = stats::p.adjust(p_value, method = p_adjust_method),
      significance_raw = sig_label(p_value),
      significance_adjusted = sig_label(p_value_adjusted),
      p_adjust_method = p_adjust_method
    )
}

write.table(
  binom_results,
  file = out_binom_csv,
  sep = ";", dec = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 16) OPTIONAL MULTINOMIAL
# ============================================================
multinom_results <- tibble::tibble()

if ("joint_type" %in% names(detected_trait_cols)) {
  jt_col <- unname(detected_trait_cols["joint_type"])
  mm <- run_multinom_joint_type(df, trait_col = jt_col, predictor = "logCS", min_n = min_n_for_group_test)
  if (!is.null(mm)) {
    multinom_results <- mm %>%
      dplyr::mutate(
        p_value_adjusted = stats::p.adjust(p_value, method = p_adjust_method),
        significance_raw = sig_label(p_value),
        significance_adjusted = sig_label(p_value_adjusted),
        p_adjust_method = p_adjust_method
      )
  }
}

write.table(
  multinom_results,
  file = out_multinom_csv,
  sep = ";", dec = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 17) PLOTS
# ============================================================
plot_list <- lapply(pcs_use, function(pc) {
  row_i <- univar_pc_results %>% dplyr::filter(trait == pc)
  
  ggplot2::ggplot(df_shape, ggplot2::aes_string(x = "logCS", y = pc)) +
    ggplot2::geom_point(size = 2.4, shape = 21, stroke = 0.5, fill = "white", color = "black", alpha = 0.9) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.0, color = "#1f1f1f") +
    ggplot2::annotate(
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
    ggplot2::labs(
      title = paste0(pc, " ~ log(Centroid size)"),
      x = "log(Centroid size)",
      y = paste0(pc, " score")
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black")
    )
})

p_uni_grid <- patchwork::wrap_plots(plotlist = plot_list, ncol = 2) +
  patchwork::plot_annotation(
    title = "Univariate allometry tests for PC1-PC5",
    subtitle = paste0("Multiple-testing correction: ", tools::toTitleCase(p_adjust_method))
  )

ggplot2::ggsave(out_plot_uni_png, p_uni_grid, width = 12, height = 14, dpi = 400)
ggplot2::ggsave(out_plot_uni_pdf, p_uni_grid, width = 12, height = 14)

p_morph <- ggplot2::ggplot(df_shape, ggplot2::aes(x = PC1, y = PC2, color = logCS)) +
  ggplot2::geom_point(size = 3, alpha = 0.95) +
  viridis::scale_color_viridis(option = "plasma", end = 0.95, name = "log(CS)", discrete = FALSE) +
  ggplot2::labs(
    title = "Morphospace PC1-PC2 colored by log(Centroid size)",
    x = "PC1",
    y = "PC2"
  ) +
  ggplot2::theme_classic(base_size = 13)

ggplot2::ggsave(out_plot_morph_png, p_morph, width = 7.2, height = 5.4, dpi = 400)
ggplot2::ggsave(out_plot_morph_pdf, p_morph, width = 7.2, height = 5.4)

summary_plot_df <- dplyr::bind_rows(
  univar_pc_results %>% dplyr::select(trait, estimate, p_value_adjusted, r_squared, trait_group, response_transform),
  cont_results %>% dplyr::select(trait, estimate, p_value_adjusted, r_squared, trait_group, response_transform)
) %>%
  dplyr::filter(!is.na(estimate), !is.na(p_value_adjusted), !is.na(r_squared))

if (nrow(summary_plot_df) > 0) {
  summary_plot_df <- summary_plot_df %>%
    dplyr::mutate(
      trait = forcats::fct_reorder(trait, r_squared),
      sig = p_value_adjusted < 0.05
    )
  
  p_summary <- ggplot2::ggplot(summary_plot_df, ggplot2::aes(x = r_squared, y = trait)) +
    ggplot2::geom_point(ggplot2::aes(shape = sig, size = abs(estimate), fill = trait_group), color = "black") +
    ggplot2::scale_shape_manual(values = c(21, 24), name = "adj. p < 0.05") +
    ggplot2::labs(
      title = "Allometry summary across shape and screw-geometry traits",
      subtitle = "x-axis = R, point size = absolute slope estimate",
      x = expression(R^2),
      y = NULL,
      fill = "Trait group",
      size = "|Slope|"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  ggplot2::ggsave(out_plot_summary_png, p_summary, width = 10, height = 8, dpi = 400)
  ggplot2::ggsave(out_plot_summary_pdf, p_summary, width = 10, height = 8)
}

geom_candidates_for_panels <- cont_results %>%
  dplyr::filter(trait_group %in% c("geometry_continuous", "geometry_ratio")) %>%
  dplyr::group_by(trait) %>%
  dplyr::slice_min(order_by = p_value_adjusted, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p_value_adjusted, desc(r_squared)) %>%
  dplyr::slice_head(n = max_panels_geometry)

if (nrow(geom_candidates_for_panels) > 0) {
  panel_traits <- geom_candidates_for_panels$trait
  
  geom_panel_list <- lapply(panel_traits, function(tr) {
    tmp <- df %>%
      dplyr::select(logCS, dplyr::all_of(tr)) %>%
      dplyr::rename(y = dplyr::all_of(tr)) %>%
      dplyr::filter(!is.na(logCS), !is.na(y))
    
    if (nrow(tmp) < 5) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(tr))
    
    ggplot2::ggplot(tmp, ggplot2::aes(x = logCS, y = y)) +
      ggplot2::geom_point(size = 2.3, shape = 21, fill = "white", color = "black") +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "#1f1f1f") +
      ggplot2::labs(title = tr, x = "log(Centroid size)", y = tr) +
      ggplot2::theme_classic(base_size = 11)
  })
  
  p_geom_grid <- patchwork::wrap_plots(geom_panel_list, ncol = 2) +
    patchwork::plot_annotation(
      title = "Selected geometry-vs-size panels",
      subtitle = "Most interesting exploratory candidates based on adjusted p and R"
    )
  
  ggplot2::ggsave(out_plot_geom_png, p_geom_grid, width = 12, height = 14, dpi = 400)
  ggplot2::ggsave(out_plot_geom_pdf, p_geom_grid, width = 12, height = 14)
}

# ============================================================
# 18) LOG FILE
# ============================================================
writeLines(log_lines, out_log_txt)

cat("\n========================================\n")
cat("DONE. All outputs were written to:\n")
cat(out_dir, "\n")
cat("========================================\n")
cat("\nMost important debug files:\n")
cat("- ", out_log_txt, "\n")
cat("- ", out_geom_match_csv, "\n")
cat("- ", out_geom_unmatched_pca_csv, "\n")
cat("- ", out_geom_unmatched_geom_csv, "\n")
cat("- ", out_manifest_txt, "\n")




###############################################################################

# ============================================================
# MAIN TEXT FIGURE: ALLOMETRY (2 x 2 GRID)
# FINAL CLEAN VERSION
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(viridis)
  library(data.table)
})

# ------------------- PATHS -------------------
allom_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/Allometry"

merged_table_path <- file.path(allom_dir, "allometry_merged_table.csv")
univar_path <- file.path(allom_dir, "allometry_univariate_PC1_to_PC5_results.csv")
cont_path   <- file.path(allom_dir, "allometry_continuous_traits_results.csv")

out_png <- file.path(allom_dir, "Figure_Main_Allometry_2x2_final_clean.png")
out_pdf <- file.path(allom_dir, "Figure_Main_Allometry_2x2_final_clean.pdf")

# ------------------- CONFIRMED MULTIVARIATE VALUES -------------------
multivar_R2 <- 0.104
multivar_F  <- 7.668
multivar_p  <- 0.0001

# ------------------- HELPERS -------------------
stop_with_hint <- function(msg) {
  stop(paste0("\n ", msg, "\n"), call. = FALSE)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

fmt_p_line <- function(label, x, digits = 3, threshold = 0.001) {
  if (is.na(x)) return(paste0(label, " = NA"))
  if (x < threshold) return(paste0(label, " < ", format(threshold, nsmall = 3)))
  paste0(label, " = ", format(round(x, digits), nsmall = digits))
}

read_csv_auto <- function(path) {
  if (!file.exists(path)) stop_with_hint(paste0("File not found: ", path))
  df <- data.table::fread(path, data.table = FALSE, encoding = "UTF-8")
  names(df) <- trimws(sub("^\ufeff", "", enc2utf8(names(df))))
  df
}

as_num_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

panel_theme_clean <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(color = "black"),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

make_reg_panel <- function(df, xvar, yvar, title, ylab, stats_label,
                           xlim_vals,
                           ylim_vals = NULL,
                           point_size = 2.7,
                           x_text = Inf, y_text = NULL) {
  if (is.null(ylim_vals)) {
    ylim_vals <- range(df[[yvar]], na.rm = TRUE)
    ypad <- diff(ylim_vals) * 0.05
    ylim_vals <- c(ylim_vals[1] - ypad, ylim_vals[2] + ypad)
  }
  
  if (is.null(y_text)) {
    y_text <- ylim_vals[2] - 0.08 * diff(ylim_vals)
  }
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]])) +
    ggplot2::geom_point(
      size = point_size,
      shape = 21,
      stroke = 0.45,
      fill = "white",
      color = "black",
      alpha = 0.95
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 1.0,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = x_text, y = y_text,
      hjust = 1.02, vjust = 1,
      label = stats_label,
      size = 3.7,
      lineheight = 1.05
    ) +
    ggplot2::labs(
      title = title,
      x = "log(Centroid size)",
      y = ylab
    ) +
    panel_theme_clean(base_size = 12) +
    ggplot2::coord_cartesian(
      xlim = xlim_vals,
      ylim = ylim_vals,
      clip = "off"
    )
}

# ------------------- READ INPUT -------------------
df <- read_csv_auto(merged_table_path)
univar_tab <- read_csv_auto(univar_path)
cont_tab <- read_csv_auto(cont_path)

required_merged_columns <- c("PC1", "PC2", "PC3", "PC4", "PC5", "centroid_size", "abs_winding_angle_deg")
missing_merged_columns <- setdiff(required_merged_columns, names(df))
if (length(missing_merged_columns) > 0) {
  stop_with_hint(paste0(
    "Missing required columns in allometry_merged_table.csv: ",
    paste(missing_merged_columns, collapse = ", ")
  ))
}

for (cc in c("PC1", "PC2", "PC3", "PC4", "PC5", "centroid_size", "abs_winding_angle_deg")) {
  df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
}

df <- df %>%
  dplyr::mutate(logCS = log(centroid_size))

df_shape <- df %>%
  dplyr::filter(!is.na(logCS), !is.na(PC1), !is.na(PC2), !is.na(PC3), !is.na(PC4), !is.na(PC5))

df_wind <- df %>%
  dplyr::filter(!is.na(logCS), !is.na(abs_winding_angle_deg))

if (nrow(df_shape) < 10) stop_with_hint("Too few complete rows for shape plotting.")
if (nrow(df_wind) < 10) stop_with_hint("Too few complete rows for winding-angle plotting.")

# ------------------- COMMON AXES -------------------
x_range <- range(c(df_shape$logCS, df_wind$logCS), na.rm = TRUE)
x_pad <- diff(x_range) * 0.05
xlim_common <- c(x_range[1] - x_pad, x_range[2] + x_pad)

pc1_ylim <- range(df_shape$PC1, na.rm = TRUE)
pc1_ylim <- pc1_ylim + c(-1, 1) * diff(pc1_ylim) * 0.05

pc2_ylim <- range(df_shape$PC2, na.rm = TRUE)
pc2_ylim <- pc2_ylim + c(-1, 1) * diff(pc2_ylim) * 0.05

wind_ylim <- range(df_wind$abs_winding_angle_deg, na.rm = TRUE)
wind_ylim <- wind_ylim + c(-1, 1) * diff(wind_ylim) * 0.05

# ------------------- EXTRACT STATS -------------------
get_univar_row <- function(tab, trait_name) {
  out <- tab %>% dplyr::filter(trait == trait_name)
  if (nrow(out) == 0) stop_with_hint(paste0("Trait '", trait_name, "' not found in univariate table."))
  out[1, , drop = FALSE]
}

pc1_row <- get_univar_row(univar_tab, "PC1")
pc2_row <- get_univar_row(univar_tab, "PC2")

pc1_p <- as_num_safe(pc1_row$p_value)
pc1_holm <- as_num_safe(pc1_row$p_value_adjusted)
pc2_p <- as_num_safe(pc2_row$p_value)
pc2_holm <- as_num_safe(pc2_row$p_value_adjusted)

pc1_label <- paste(
  paste0("R = ", fmt_num(as_num_safe(pc1_row$r_squared))),
  fmt_p_line("p", pc1_p),
  fmt_p_line("Holm p", pc1_holm),
  sep = "\n"
)

pc2_label <- paste(
  paste0("R = ", fmt_num(as_num_safe(pc2_row$r_squared))),
  fmt_p_line("p", pc2_p),
  fmt_p_line("Holm p", pc2_holm),
  sep = "\n"
)

wind_row <- cont_tab %>%
  dplyr::filter(trait == "abs_winding_angle_deg", model_type == "lm")

if (nrow(wind_row) == 0) {
  wind_row <- cont_tab %>% dplyr::filter(trait == "abs_winding_angle_deg")
}
if (nrow(wind_row) == 0) {
  stop_with_hint("Trait 'abs_winding_angle_deg' not found in continuous trait table.")
}
wind_row <- wind_row[1, , drop = FALSE]

wind_p <- as_num_safe(wind_row$p_value)
wind_holm <- as_num_safe(wind_row$p_value_adjusted)

wind_label <- paste(
  paste0("R = ", fmt_num(as_num_safe(wind_row$r_squared))),
  fmt_p_line("p", wind_p),
  fmt_p_line("Holm p", wind_holm),
  sep = "\n"
)

multivar_label <- paste(
  "PC1-PC5 ~ log(CS)",
  paste0("R = ", fmt_num(multivar_R2)),
  paste0("F = ", fmt_num(multivar_F)),
  "p < 0.001",
  sep = "\n"
)

# ------------------- PANEL A AXIS RANGE -------------------
pcA_xrange <- range(df_shape$PC1, na.rm = TRUE)
pcA_yrange <- range(df_shape$PC2, na.rm = TRUE)
pcA_xpad <- diff(pcA_xrange) * 0.10
pcA_ypad <- diff(pcA_yrange) * 0.08
pcA_xlim <- c(pcA_xrange[1] - pcA_xpad, pcA_xrange[2] + pcA_xpad)
pcA_ylim <- c(pcA_yrange[1] - pcA_ypad, pcA_yrange[2] + pcA_ypad)

panel_stats_y_frac <- 0.92
panelA_stats_y <- pcA_ylim[1] + panel_stats_y_frac * diff(pcA_ylim)
panelA_stats_x <- pcA_xlim[2] - 0.01 * diff(pcA_xlim)

# smaller, more subtle legend inside panel A
legend_pos_x <- 0.88
legend_pos_y <- 0.47

# ------------------- BUILD PANELS -------------------
pA <- ggplot2::ggplot(df_shape, ggplot2::aes(x = PC1, y = PC2, fill = logCS)) +
  ggplot2::geom_point(
    size = 3.0,
    alpha = 0.98,
    shape = 21,
    color = "black",
    stroke = 0.20
  ) +
  ggplot2::annotate(
    "text",
    x = panelA_stats_x,
    y = panelA_stats_y,
    hjust = 1.02,
    vjust = 1,
    label = multivar_label,
    size = 3.7,
    lineheight = 1.05
  ) +
  viridis::scale_fill_viridis(
    option = "plasma",
    end = 0.95,
    discrete = FALSE,
    name = "log(CS)"
  ) +
  ggplot2::labs(
    title = NULL,
    subtitle = NULL,
    x = "PC1",
    y = "PC2"
  ) +
  panel_theme_clean(base_size = 12) +
  ggplot2::coord_cartesian(
    xlim = pcA_xlim,
    ylim = pcA_ylim,
    clip = "off"
  ) +
  ggplot2::theme(
    legend.position = c(legend_pos_x, legend_pos_y),
    legend.direction = "vertical",
    legend.background = ggplot2::element_blank(),
    legend.key.height = grid::unit(0.65, "cm"),
    legend.key.width = grid::unit(0.22, "cm"),
    legend.title = ggplot2::element_text(size = 10, face = "bold"),
    legend.text = ggplot2::element_text(size = 8.5)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      barheight = grid::unit(2.4, "cm"),
      barwidth = grid::unit(0.28, "cm"),
      ticks = TRUE
    )
  )

pB <- make_reg_panel(
  df = df_shape,
  xvar = "logCS",
  yvar = "PC1",
  title = "PC1",
  ylab = "PC1 score",
  stats_label = pc1_label,
  xlim_vals = xlim_common,
  ylim_vals = pc1_ylim,
  y_text = pc1_ylim[1] + panel_stats_y_frac * diff(pc1_ylim)
)

pC <- make_reg_panel(
  df = df_shape,
  xvar = "logCS",
  yvar = "PC2",
  title = "PC2",
  ylab = "PC2 score",
  stats_label = pc2_label,
  xlim_vals = xlim_common,
  ylim_vals = pc2_ylim,
  y_text = pc2_ylim[1] + panel_stats_y_frac * diff(pc2_ylim)
)

pD <- make_reg_panel(
  df = df_wind,
  xvar = "logCS",
  yvar = "abs_winding_angle_deg",
  title = NULL,
  ylab = "Absolute winding angle ()",
  stats_label = wind_label,
  xlim_vals = xlim_common,
  ylim_vals = wind_ylim,
  y_text = wind_ylim[1] + 0.88 * diff(wind_ylim)
)

# ------------------- COMBINE -------------------
final_plot <- (pA + pB) / (pC + pD) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.99)
  )

# ------------------- SAVE -------------------
ggplot2::ggsave(
  filename = out_png,
  plot = final_plot,
  width = 13.5,
  height = 10.5,
  dpi = 500,
  bg = "white"
)

ggplot2::ggsave(
  filename = out_pdf,
  plot = final_plot,
  width = 13.5,
  height = 10.5,
  bg = "white"
)

cat("\n========================================\n")
cat("DONE. Figure written to:\n")
cat(out_png, "\n")
cat(out_pdf, "\n")
cat("========================================\n")




################################################################################


# ============================================================
# SUPPLEMENTARY ALLOMETRY PLOTS
# Creates all remaining allometry plots not included in the main figure
#
# Output folder:
#   .../Allometry/Supplementary_Plots
#
# Included:
#   - PC3, PC4, PC5
#   - all remaining continuous geometry traits except abs_winding_angle_deg
# ============================================================

# ------------------- PATHS -------------------
supp_dir <- file.path(allom_dir, "Supplementary_Plots")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

out_pc_supp_png   <- file.path(supp_dir, "Supp_Allometry_PC3_PC5.png")
out_pc_supp_pdf   <- file.path(supp_dir, "Supp_Allometry_PC3_PC5.pdf")
out_geom_supp_png <- file.path(supp_dir, "Supp_Allometry_Geometry_Remaining.png")
out_geom_supp_pdf <- file.path(supp_dir, "Supp_Allometry_Geometry_Remaining.pdf")

# ------------------- HELPERS -------------------
format_trait_label <- function(x) {
  dplyr::case_when(
    x == "signed_winding_angle_deg"   ~ "Signed winding angle ()",
    x == "abs_winding_angle_deg"      ~ "Absolute winding angle ()",
    x == "n_turns_signed"             ~ "Signed number of turns",
    x == "n_turns_abs"                ~ "Absolute number of turns",
    x == "start_end_dist"             ~ "Start-end distance",
    x == "axial_span"                 ~ "Axial span",
    x == "fit_radius"                 ~ "Fitted radius",
    x == "fit_rms"                    ~ "Fit RMS",
    x == "ratio_axial_span_fit_radius" ~ "Axial span / fitted radius",
    x == "ratio_start_end_fit_radius"  ~ "Start-end distance / fitted radius",
    x == "ratio_fit_rms_fit_radius"    ~ "Fit RMS / fitted radius",
    TRUE ~ x
  )
}

make_stats_label <- function(r2, p, padj = NA_real_) {
  lines <- c(
    paste0("R = ", fmt_num(r2)),
    fmt_p_line("p", p)
  )
  if (!is.na(padj)) {
    lines <- c(lines, fmt_p_line("Holm p", padj))
  }
  paste(lines, collapse = "\n")
}

make_simple_reg_panel <- function(df_plot, yvar, title, ylab, stats_label, xlim_vals) {
  yvals <- df_plot[[yvar]]
  ylim_vals <- range(yvals, na.rm = TRUE)
  ypad <- diff(ylim_vals) * 0.05
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.05
  ylim_vals <- c(ylim_vals[1] - ypad, ylim_vals[2] + ypad)
  
  y_text <- ylim_vals[1] + 0.92 * diff(ylim_vals)
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = logCS, y = .data[[yvar]])) +
    ggplot2::geom_point(
      size = 2.5,
      shape = 21,
      stroke = 0.4,
      fill = "white",
      color = "black",
      alpha = 0.95
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.95,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = Inf, y = y_text,
      hjust = 1.02, vjust = 1,
      label = stats_label,
      size = 3.5,
      lineheight = 1.05
    ) +
    ggplot2::labs(
      title = title,
      x = "log(Centroid size)",
      y = ylab
    ) +
    panel_theme_clean(base_size = 11) +
    ggplot2::coord_cartesian(
      xlim = xlim_vals,
      ylim = ylim_vals,
      clip = "off"
    )
}

choose_cont_row <- function(tab, trait_name) {
  # Prefer raw lm row for plotting raw data
  out <- tab %>%
    dplyr::filter(trait == trait_name)
  
  if (nrow(out) == 0) return(NULL)
  
  out_lm <- out %>% dplyr::filter(model_type == "lm")
  if (nrow(out_lm) > 0) return(out_lm[1, , drop = FALSE])
  
  out[1, , drop = FALSE]
}

# ------------------- COMMON X LIMITS -------------------
# Use same x-limits as main regression panels
x_range_supp <- range(df$logCS, na.rm = TRUE)
x_pad_supp <- diff(x_range_supp) * 0.05
xlim_supp <- c(x_range_supp[1] - x_pad_supp, x_range_supp[2] + x_pad_supp)

# ============================================================
# 1) SUPPLEMENT: PC3-PC5
# ============================================================
pcs_supp <- c("PC3", "PC4", "PC5")

pc_supp_plots <- lapply(pcs_supp, function(pc) {
  row_i <- univar_tab %>% dplyr::filter(trait == pc)
  
  if (nrow(row_i) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle(pc)
    )
  }
  
  df_plot <- df %>%
    dplyr::filter(!is.na(logCS), !is.na(.data[[pc]]))
  
  stats_label <- make_stats_label(
    r2   = as_num_safe(row_i$r_squared[1]),
    p    = as_num_safe(row_i$p_value[1]),
    padj = as_num_safe(row_i$p_value_adjusted[1])
  )
  
  make_simple_reg_panel(
    df_plot = df_plot,
    yvar = pc,
    title = pc,
    ylab = paste0(pc, " score"),
    stats_label = stats_label,
    xlim_vals = xlim_supp
  )
})

p_pc_supp <- patchwork::wrap_plots(pc_supp_plots, ncol = 2) +
  patchwork::plot_annotation(
    title = "Supplementary allometry plots: PC3-PC5",
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 13),
    plot.tag.position = c(0.01, 0.99)
  )

ggplot2::ggsave(
  filename = out_pc_supp_png,
  plot = p_pc_supp,
  width = 11,
  height = 8.5,
  dpi = 500,
  bg = "white"
)

ggplot2::ggsave(
  filename = out_pc_supp_pdf,
  plot = p_pc_supp,
  width = 11,
  height = 8.5,
  bg = "white"
)

# ============================================================
# 2) SUPPLEMENT: REMAINING GEOMETRY TRAITS
# ============================================================
# Pull remaining continuous traits directly from the results table
remaining_geom_traits <- cont_tab %>%
  dplyr::distinct(trait) %>%
  dplyr::pull(trait) %>%
  setdiff("abs_winding_angle_deg")

# Keep only traits that actually exist in merged table
remaining_geom_traits <- remaining_geom_traits[remaining_geom_traits %in% names(df)]

geom_supp_plots <- lapply(remaining_geom_traits, function(tr) {
  row_i <- choose_cont_row(cont_tab, tr)
  
  if (is.null(row_i)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle(tr)
    )
  }
  
  df_plot <- df %>%
    dplyr::filter(!is.na(logCS), !is.na(.data[[tr]]))
  
  if (nrow(df_plot) < 5) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle(format_trait_label(tr))
    )
  }
  
  stats_label <- make_stats_label(
    r2   = as_num_safe(row_i$r_squared[1]),
    p    = as_num_safe(row_i$p_value[1]),
    padj = as_num_safe(row_i$p_value_adjusted[1])
  )
  
  make_simple_reg_panel(
    df_plot = df_plot,
    yvar = tr,
    title = format_trait_label(tr),
    ylab = format_trait_label(tr),
    stats_label = stats_label,
    xlim_vals = xlim_supp
  )
})

# choose sensible ncol depending on number of panels
geom_n <- length(geom_supp_plots)
geom_ncol <- ifelse(geom_n <= 4, 2, 3)

p_geom_supp <- patchwork::wrap_plots(geom_supp_plots, ncol = geom_ncol) +
  patchwork::plot_annotation(
    title = "Supplementary allometry plots: remaining screw-geometry traits",
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 13),
    plot.tag.position = c(0.01, 0.99)
  )

# dynamic height
geom_height <- ifelse(geom_n <= 4, 8.5,
                      ifelse(geom_n <= 6, 11,
                             ifelse(geom_n <= 9, 14, 17)))

ggplot2::ggsave(
  filename = out_geom_supp_png,
  plot = p_geom_supp,
  width = 14,
  height = geom_height,
  dpi = 500,
  bg = "white"
)

ggplot2::ggsave(
  filename = out_geom_supp_pdf,
  plot = p_geom_supp,
  width = 14,
  height = geom_height,
  bg = "white"
)

cat("\n========================================\n")
cat("SUPPLEMENTARY ALLOMETRY PLOTS DONE\n")
cat("Written to:\n")
cat(out_pc_supp_png, "\n")
cat(out_pc_supp_pdf, "\n")
cat(out_geom_supp_png, "\n")
cat(out_geom_supp_pdf, "\n")
cat("========================================\n")



################################################################################
################################################################################
################################################################################
# ============================================================


# ============================================================
# PHYLOGENETICALLY INFORMED ALLOMETRY (PGLS)
# AGGREGATED-TO-TREE-TIP VERSION
#
# Fits univariate PGLS for:
#   - PC1 ~ logCS
#   - PC2 ~ logCS
#   - abs_winding_angle_deg ~ logCS
#
# Key fixes:
# - Aggregates specimen-level data to one value per tree tip
# - Robustly matches tree labels against tree_tip and Family___tree_tip
# - Removes internal node labels from tree to avoid caper duplication error
# - Moves stats text lower for abs_winding_angle_deg panel
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ape)
  library(caper)
  library(ggplot2)
  library(patchwork)
})

# ------------------- PATHS -------------------
merged_table_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/Allometry/allometry_merged_table.csv"
key_path    <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_key.csv"
tree_path   <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/curc_fig1_withCaridae_calibrated_Grafen.tre"

base_results_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results"
out_dir <- file.path(base_results_dir, "Allometry_Phylogenetic")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_input_specimen_csv <- file.path(out_dir, "pgls_input_specimen_level.csv")
out_input_agg_csv      <- file.path(out_dir, "pgls_input_aggregated_to_tree_tip.csv")
out_tree_used          <- file.path(out_dir, "pgls_tree_used.tre")
out_results_csv        <- file.path(out_dir, "pgls_results_main_traits.csv")
out_summary_txt        <- file.path(out_dir, "pgls_model_summaries.txt")
out_log_txt            <- file.path(out_dir, "pgls_log.txt")
out_match_diag_csv     <- file.path(out_dir, "pgls_tree_label_matching_diagnostics.csv")

out_pc1_png            <- file.path(out_dir, "pgls_PC1.png")
out_pc1_pdf            <- file.path(out_dir, "pgls_PC1.pdf")
out_pc2_png            <- file.path(out_dir, "pgls_PC2.png")
out_pc2_pdf            <- file.path(out_dir, "pgls_PC2.pdf")
out_wind_png           <- file.path(out_dir, "pgls_abs_winding_angle.png")
out_wind_pdf           <- file.path(out_dir, "pgls_abs_winding_angle.pdf")
out_combined_png       <- file.path(out_dir, "pgls_main_traits_combined.png")
out_combined_pdf       <- file.path(out_dir, "pgls_main_traits_combined.pdf")

# ------------------- SETTINGS -------------------
traits_to_test <- c("PC1", "PC2", "abs_winding_angle_deg")
lambda_bounds  <- c(1e-6, 1.5)
min_tips_needed <- 10

# ------------------- LOGGING -------------------
log_lines <- character()

log_msg <- function(...) {
  txt <- paste0(...)
  message(txt)
  log_lines <<- c(log_lines, txt)
}

# ------------------- HELPERS -------------------
stop_with_hint <- function(msg) {
  stop(paste0("\n ", msg, "\n"), call. = FALSE)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", format.pval(x, digits = 3, eps = 0.001))
}

read_csv_auto <- function(path) {
  if (!file.exists(path)) stop_with_hint(paste0("File not found: ", path))
  df <- data.table::fread(path, data.table = FALSE, encoding = "UTF-8")
  names(df) <- trimws(sub("^\ufeff", "", enc2utf8(names(df))))
  df
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

norm_id <- function(x) {
  x <- as.character(x)
  x <- enc2utf8(x)
  x <- trimws(x)
  x <- sub("^\ufeff", "", x)
  x <- tolower(x)
  x <- sub("\\.txt$", "", x, ignore.case = TRUE)
  x <- sub("\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("\\.vtk$", "", x, ignore.case = TRUE)
  x <- sub("_trochanter.*$", "", x, ignore.case = TRUE)
  x <- gsub("_aligned$", "", x, ignore.case = TRUE)
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

norm_tree_label <- function(x) {
  x <- as.character(x)
  x <- enc2utf8(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

fix_known_tree_labels <- function(x) {
  x <- gsub("^caridae___neydus$", "caridae___nedyus", x, ignore.case = TRUE)
  x <- gsub("^caridae_neydus$", "caridae_nedyus", x, ignore.case = TRUE)
  x <- gsub("^curculionidae___neydus$", "curculionidae___nedyus", x, ignore.case = TRUE)
  x <- gsub("^curculionidae_neydus$", "curculionidae_nedyus", x, ignore.case = TRUE)
  x <- gsub("^neydus$", "nedyus", x, ignore.case = TRUE)
  x <- gsub("^belidae___belidae$", "belidae___agnesiotis", x, ignore.case = TRUE)
  x <- gsub("^belidae_belidae$", "belidae_agnesiotis", x, ignore.case = TRUE)
  x <- gsub("^caridae___caridae$", "caridae___car", x, ignore.case = TRUE)
  x <- gsub("^caridae_caridae$", "caridae_car", x, ignore.case = TRUE)
  x <- gsub("^belidae$", "agnesiotis", x, ignore.case = TRUE)
  x <- gsub("^caridae$", "car", x, ignore.case = TRUE)
  x
}

build_tree_label_family_tip <- function(family, tree_tip) {
  out <- paste(family, tree_tip, sep = "___")
  out <- norm_tree_label(out)
  out <- fix_known_tree_labels(out)
  out
}

sanitize_tree_for_caper <- function(tr) {
  tr$tip.label <- norm_tree_label(tr$tip.label)
  tr$tip.label <- fix_known_tree_labels(tr$tip.label)
  
  if (!is.null(tr$node.label)) {
    n_nonempty <- sum(!is.na(tr$node.label) & trimws(tr$node.label) != "")
    if (n_nonempty > 0) {
      log_msg(" Removing ", n_nonempty, " internal node labels from tree for caper compatibility.")
    }
    tr$node.label <- NULL
  }
  
  tr
}

panel_theme_clean <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

get_coef_row <- function(coef_tab, predictor = "logCS") {
  rn <- rownames(coef_tab)
  idx <- which(rn == predictor)
  if (length(idx) == 0) return(NULL)
  coef_tab[idx[1], , drop = FALSE]
}

extract_pgls_stats <- function(model, trait_name) {
  sm <- summary(model)
  cf <- sm$coefficients
  row_i <- get_coef_row(cf, "logCS")
  
  if (is.null(row_i)) {
    return(tibble::tibble(
      trait = trait_name,
      n = nobs(model),
      intercept = NA_real_,
      slope = NA_real_,
      std_error = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      lambda = if (!is.null(model$param["lambda"])) unname(model$param["lambda"]) else NA_real_,
      logLik = as.numeric(logLik(model)),
      AIC = AIC(model)
    ))
  }
  
  tibble::tibble(
    trait = trait_name,
    n = nobs(model),
    intercept = cf["(Intercept)", "Estimate"],
    slope = row_i[1, "Estimate"],
    std_error = row_i[1, "Std. Error"],
    t_value = row_i[1, "t value"],
    p_value = row_i[1, "Pr(>|t|)"],
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    lambda = if (!is.null(model$param["lambda"])) unname(model$param["lambda"]) else NA_real_,
    logLik = as.numeric(logLik(model)),
    AIC = AIC(model)
  )
}

fit_pgls_with_fallback <- function(formula_obj, comp_obj, lambda_bounds = c(1e-06, 1.5)) {
  attempts <- list(
    list(label = "lambda_ML", lambda = "ML", bounds = list(lambda = lambda_bounds)),
    list(label = "lambda_fixed_1", lambda = 1, bounds = NULL),
    list(label = "lambda_fixed_0", lambda = 0, bounds = NULL)
  )
  
  last_error <- NULL
  
  for (att in attempts) {
    fit_try <- tryCatch(
      {
        if (identical(att$lambda, "ML")) {
          caper::pgls(
            formula = formula_obj,
            data = comp_obj,
            lambda = att$lambda,
            bounds = att$bounds
          )
        } else {
          caper::pgls(
            formula = formula_obj,
            data = comp_obj,
            lambda = att$lambda
          )
        }
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    
    if (!is.null(fit_try)) {
      return(list(model = fit_try, method = att$label, error = NULL))
    }
  }
  
  list(model = NULL, method = NA_character_, error = last_error)
}

plot_pgls_panel <- function(df_plot, yvar, title, ylab, stats_label, outfile_png, outfile_pdf,
                            stats_y_frac = 0.92) {
  x_range <- range(df_plot$logCS, na.rm = TRUE)
  x_pad <- diff(x_range) * 0.05
  xlim_vals <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  
  y_range <- range(df_plot[[yvar]], na.rm = TRUE)
  y_pad <- diff(y_range) * 0.05
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.05
  ylim_vals <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  
  y_text <- ylim_vals[1] + stats_y_frac * diff(ylim_vals)
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = logCS, y = .data[[yvar]])) +
    ggplot2::geom_point(
      size = 2.7,
      shape = 21,
      stroke = 0.45,
      fill = "white",
      color = "black",
      alpha = 0.95
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 1.0,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = Inf, y = y_text,
      hjust = 1.02, vjust = 1,
      label = stats_label,
      size = 3.6,
      lineheight = 1.05
    ) +
    ggplot2::labs(
      title = title,
      x = "log(Centroid size)",
      y = ylab
    ) +
    panel_theme_clean(base_size = 12) +
    ggplot2::coord_cartesian(
      xlim = xlim_vals,
      ylim = ylim_vals,
      clip = "off"
    )
  
  ggplot2::ggsave(outfile_png, p, width = 6.8, height = 5.4, dpi = 500, bg = "white")
  ggplot2::ggsave(outfile_pdf, p, width = 6.8, height = 5.4, bg = "white")
  
  p
}

# ------------------- READ INPUTS -------------------
merged_data <- read_csv_auto(merged_table_path)
key <- read_csv_auto(key_path)

if (!file.exists(tree_path)) stop_with_hint(paste0("Tree file not found: ", tree_path))
tree <- ape::read.tree(tree_path)
tree <- sanitize_tree_for_caper(tree)

log_msg("Merged-table rows: ", nrow(merged_data))
log_msg("Key rows: ", nrow(key))
log_msg("Tree tips: ", length(tree$tip.label))

# ------------------- BASIC CHECKS -------------------
required_key <- c("specimen_id", "Family", "tree_tip")
missing_key <- setdiff(required_key, names(key))
if (length(missing_key) > 0) {
  stop_with_hint(paste0("specimen_key.csv is missing required columns: ",
                        paste(missing_key, collapse = ", ")))
}

required_merged_columns <- c("specimen_id", "centroid_size", traits_to_test)
missing_merged_columns <- setdiff(required_merged_columns, names(merged_data))
if (length(missing_merged_columns) > 0) {
  stop_with_hint(paste0("allometry_merged_table.csv is missing required columns: ",
                        paste(missing_merged_columns, collapse = ", ")))
}

# ------------------- PREP TREE LABELS -------------------
merged_data2 <- merged_data %>%
  dplyr::mutate(
    specimen_id_norm = vapply(specimen_id, norm_id, FUN.VALUE = character(1)),
    centroid_size = safe_num(centroid_size),
    logCS = log(centroid_size)
  )

for (tr in traits_to_test) {
  merged_data2[[tr]] <- safe_num(merged_data2[[tr]])
}

key2 <- key %>%
  dplyr::mutate(
    specimen_id_norm = vapply(specimen_id, norm_id, FUN.VALUE = character(1)),
    tree_tip_norm = fix_known_tree_labels(norm_tree_label(tree_tip)),
    tree_label_family_tip = build_tree_label_family_tip(Family, tree_tip),
    tree_label_tip_only   = tree_tip_norm
  ) %>%
  dplyr::select(
    specimen_id_norm, specimen_id, Family, tree_tip,
    tree_label_family_tip, tree_label_tip_only
  ) %>%
  dplyr::distinct()

dat <- merged_data2 %>%
  dplyr::left_join(key2, by = "specimen_id_norm") %>%
  dplyr::filter(!is.na(logCS))

# choose the matching strategy with higher overlap
n_match_family_tip <- sum(unique(dat$tree_label_family_tip) %in% tree$tip.label, na.rm = TRUE)
n_match_tip_only   <- sum(unique(dat$tree_label_tip_only)   %in% tree$tip.label, na.rm = TRUE)

log_msg("Unique matches using Family___tree_tip: ", n_match_family_tip)
log_msg("Unique matches using tree_tip only: ", n_match_tip_only)

if (n_match_family_tip >= n_match_tip_only) {
  dat <- dat %>% dplyr::mutate(tree_label = tree_label_family_tip)
  chosen_strategy <- "Family___tree_tip"
} else {
  dat <- dat %>% dplyr::mutate(tree_label = tree_label_tip_only)
  chosen_strategy <- "tree_tip_only"
}

log_msg("Chosen tree-label strategy: ", chosen_strategy)

dat <- dat %>%
  dplyr::filter(!is.na(tree_label))

write.csv2(dat, out_input_specimen_csv, row.names = FALSE)

# diagnostics
match_diag <- dat %>%
  dplyr::mutate(matches_tree = tree_label %in% tree$tip.label) %>%
  dplyr::count(tree_label, matches_tree, name = "n_specimens") %>%
  dplyr::arrange(dplyr::desc(matches_tree), dplyr::desc(n_specimens), tree_label)

write.csv2(match_diag, out_match_diag_csv, row.names = FALSE)

# retain only matched labels
dat_matched <- dat %>%
  dplyr::filter(tree_label %in% tree$tip.label)

if (nrow(dat_matched) == 0) {
  stop_with_hint("No matched rows between data and tree after trying both tree-label strategies.")
}

log_msg("Matched specimen rows: ", nrow(dat_matched))
log_msg("Unique matched tree labels before aggregation: ", dplyr::n_distinct(dat_matched$tree_label))

# ------------------- AGGREGATE TO TREE TIP -------------------
agg_dat <- dat_matched %>%
  dplyr::group_by(tree_label) %>%
  dplyr::summarise(
    n_specimens = dplyr::n(),
    logCS = mean(logCS, na.rm = TRUE),
    centroid_size_mean = mean(centroid_size, na.rm = TRUE),
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    abs_winding_angle_deg = mean(abs_winding_angle_deg, na.rm = TRUE),
    .groups = "drop"
  )

write.csv2(agg_dat, out_input_agg_csv, row.names = FALSE)

log_msg("Rows after aggregation to one row per tree tip: ", nrow(agg_dat))

keep_tips <- intersect(tree$tip.label, agg_dat$tree_label)
if (length(keep_tips) < min_tips_needed) {
  stop_with_hint(paste0("Too few overlapping tree tips after aggregation. Found: ", length(keep_tips)))
}

tree_used <- ape::drop.tip(tree, setdiff(tree$tip.label, keep_tips))
tree_used <- sanitize_tree_for_caper(tree_used)

agg_used <- agg_dat %>%
  dplyr::filter(tree_label %in% tree_used$tip.label)

agg_used <- agg_used[match(tree_used$tip.label, agg_used$tree_label), , drop = FALSE]

if (!all(agg_used$tree_label == tree_used$tip.label)) {
  stop_with_hint("Aggregated data order could not be matched to tree order.")
}

ape::write.tree(tree_used, file = out_tree_used)

log_msg("Final tree tips used: ", length(tree_used$tip.label))
log_msg("Final aggregated rows used: ", nrow(agg_used))

# ------------------- FIT PGLS MODELS -------------------
results_list <- list()
summary_text <- character()
plot_list <- list()
failed_traits <- character()

for (tr in traits_to_test) {
  df_tr <- agg_used %>%
    dplyr::filter(!is.na(.data[[tr]]), !is.na(logCS))
  
  log_msg("Trait ", tr, ": aggregated rows with complete data = ", nrow(df_tr))
  
  if (nrow(df_tr) < min_tips_needed) {
    log_msg(" Trait ", tr, " skipped: too few rows.")
    next
  }
  
  tree_tr <- ape::drop.tip(tree_used, setdiff(tree_used$tip.label, df_tr$tree_label))
  tree_tr <- sanitize_tree_for_caper(tree_tr)
  
  df_tr <- df_tr[match(tree_tr$tip.label, df_tr$tree_label), , drop = FALSE]
  
  comp <- caper::comparative.data(
    phy = tree_tr,
    data = as.data.frame(df_tr),
    names.col = "tree_label",
    vcv = TRUE,
    warn.dropped = TRUE
  )
  
  form <- stats::as.formula(paste(tr, "~ logCS"))
  
  fit_info <- fit_pgls_with_fallback(
    formula_obj = form,
    comp_obj = comp,
    lambda_bounds = lambda_bounds
  )
  fit <- fit_info$model
  
  if (is.null(fit)) {
    log_msg(" Trait ", tr, " skipped: PGLS failed across fallback strategies. Last error: ", fit_info$error)
    failed_traits <- c(failed_traits, tr)
    next
  }
  
  log_msg("Trait ", tr, ": PGLS fitted with strategy ", fit_info$method)
  
  res_row <- extract_pgls_stats(fit, tr)
  res_row$model_fit_strategy <- fit_info$method
  results_list[[tr]] <- res_row
  
  sm <- capture.output(summary(fit))
  summary_text <- c(
    summary_text,
    paste0("============================================================"),
    paste0("Trait: ", tr),
    paste0("Fit strategy: ", fit_info$method),
    paste0("============================================================"),
    sm,
    ""
  )
  
  stats_label <- paste(
    "PGLS (lambda ML)",
    paste0("R = ", fmt_num(res_row$r_squared)),
    paste0("p = ", fmt_p(res_row$p_value)),
    paste0("lambda = ", fmt_num(res_row$lambda)),
    sep = "\n"
  )
  
  title_i <- dplyr::case_when(
    tr == "PC1" ~ "Phylogenetic allometry: PC1",
    tr == "PC2" ~ "Phylogenetic allometry: PC2",
    tr == "abs_winding_angle_deg" ~ "Phylogenetic allometry: Absolute winding angle",
    TRUE ~ paste("Phylogenetic allometry:", tr)
  )
  
  ylab_i <- dplyr::case_when(
    tr == "PC1" ~ "PC1 score",
    tr == "PC2" ~ "PC2 score",
    tr == "abs_winding_angle_deg" ~ "Absolute winding angle ()",
    TRUE ~ tr
  )
  
  out_png_i <- switch(
    tr,
    "PC1" = out_pc1_png,
    "PC2" = out_pc2_png,
    "abs_winding_angle_deg" = out_wind_png,
    file.path(out_dir, paste0("pgls_", tr, ".png"))
  )
  
  out_pdf_i <- switch(
    tr,
    "PC1" = out_pc1_pdf,
    "PC2" = out_pc2_pdf,
    "abs_winding_angle_deg" = out_wind_pdf,
    file.path(out_dir, paste0("pgls_", tr, ".pdf"))
  )
  
  stats_y_frac_i <- if (tr == "abs_winding_angle_deg") 0.74 else 0.92
  
  plot_list[[tr]] <- plot_pgls_panel(
    df_plot = df_tr,
    yvar = tr,
    title = title_i,
    ylab = ylab_i,
    stats_label = stats_label,
    outfile_png = out_png_i,
    outfile_pdf = out_pdf_i,
    stats_y_frac = stats_y_frac_i
  )
}

# ------------------- WRITE RESULTS -------------------
results_tab <- dplyr::bind_rows(results_list)

if (nrow(results_tab) == 0) {
  stop_with_hint("No PGLS models could be fitted.")
}

results_tab <- results_tab %>%
  dplyr::mutate(
    p_value_adjusted = stats::p.adjust(p_value, method = "holm")
  )

write.csv2(results_tab, out_results_csv, row.names = FALSE)
writeLines(summary_text, out_summary_txt)
if (length(failed_traits) > 0) {
  log_msg("Traits skipped after all PGLS fallback attempts: ", paste(unique(failed_traits), collapse = ", "))
}
writeLines(log_lines, out_log_txt)

# ------------------- COMBINED FIGURE -------------------
plots_to_combine <- plot_list[c("PC1", "PC2", "abs_winding_angle_deg")]
plots_to_combine <- plots_to_combine[!vapply(plots_to_combine, is.null, logical(1))]

if (length(plots_to_combine) > 0) {
  p_combined <- patchwork::wrap_plots(plots_to_combine, ncol = 2) +
    patchwork::plot_annotation(
      title = "Phylogenetically informed allometry (PGLS)",
      tag_levels = "A"
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99)
    )
  
  ggplot2::ggsave(out_combined_png, p_combined, width = 12.5, height = 8.5, dpi = 500, bg = "white")
  ggplot2::ggsave(out_combined_pdf, p_combined, width = 12.5, height = 8.5, bg = "white")
}

cat("\n============================================================\n")
cat("DONE: PHYLOGENETICALLY INFORMED ALLOMETRY (PGLS)\n")
cat("Output folder:\n")
cat(out_dir, "\n")
cat("Main result table:\n")
cat(out_results_csv, "\n")
cat("============================================================\n")



