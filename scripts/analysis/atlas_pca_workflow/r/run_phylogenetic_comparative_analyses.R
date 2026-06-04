# ============================================================
# Phylogenetic comparative analysis workflow
# Project: Curculionoidea coxa-trochanter joint paper
# Purpose:
#   Exploratory phylogenetic comparative analyses (PCM) on
#   genus-level aggregated shape and joint-geometry data.
#
# IMPORTANT:
#   This script treats PCM as exploratory / trend-seeking.
#   Results should be interpreted cautiously.
#
# KEY DATA LOGIC:
#   - specimen-level data are merged by specimen_id
#   - multiple specimens per phylogenetic tip are averaged
#   - the phylogeny uses labels in the form Family___Genus
#   - therefore this script builds a matching column:
#         tree_label = paste(Family, tree_tip, sep = "___")
#   - all tree-aware analyses use tree_label, NOT tree_tip alone
# ============================================================

# ============================================================
# SECTION 0: PACKAGE SETUP
# ============================================================

required_packages <- c(
  "ape", "phytools", "geiger", "caper", "nlme", "mvMORPH",
  "tidyverse", "data.table", "ggrepel", "patchwork",
  "RColorBrewer", "scales"
)

install_missing_packages <- FALSE

missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
  message("Missing packages detected: ", paste(missing_packages, collapse = ", "))
  if (isTRUE(install_missing_packages)) {
    install.packages(missing_packages, dependencies = TRUE)
  } else {
    stop(
      "Please install the missing packages before running the script:\n",
      paste(missing_packages, collapse = ", ")
    )
  }
}

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(geiger)
  library(caper)
  library(nlme)
  library(mvMORPH)
  library(tidyverse)
  library(data.table)
  library(ggrepel)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
})

options(stringsAsFactors = FALSE)
theme_set(ggplot2::theme_bw(base_size = 12))

# ============================================================
# SECTION 1: USER SETTINGS
# ============================================================

tree_file <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/curc_fig1_withCaridae_calibrated_Grafen.tre"
tree_variant_dir <- "<MANUSCRIPT_PROJECT_ROOT>/Phylogeny/Dataset_S4_Supermatrices_partitions/Supermatrix for Fig. 1/AA"
pc_file <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"
specimen_key_file <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_key.csv"
geometry_file <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/winding_metrics_excel.csv"
centroid_file <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id_with_centroid_size.csv"
output_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/PCM"

tree_variant_files <- c(
  ml_unrooted = "curc_fig1.treefile",
  consensus_unrooted = "curc_fig1.contree",
  rooted_ml = "curc_fig1_rooted.tre",
  ultrametric = "curc_fig1_ultrametric.tre",
  grafen = "curc_fig1_grafen.tre",
  with_caridae_ml = "curc_fig1_withCaridae_ML.tre",
  with_caridae_calibrated = "curc_fig1_withCaridae_calibrated.tre",
  with_caridae_calibrated_grafen = "curc_fig1_withCaridae_calibrated_Grafen.tre",
  ultrametric_with_caridae = "curc_fig1_ultrametric_withCaridae_correct.tre",
  grafen_with_caridae = "curc_fig1_grafen_withCaridae_correct.tre"
)

tree_outgroup_label <- "Nemonychidae___Rhynchitomacerinus"
tree_topology_stress_n <- 0

pc_vars_main <- paste0("PC", 1:5)

geometry_vars_main <- c(
  "abs_winding_angle_deg",
  "n_turns_abs",
  "start_end_dist",
  "axial_span",
  "fit_radius"
)

size_var <- "centroid_size"

discrete_vars_main <- c(
  "Coxal wall hole",
  "Coxal Socket",
  "Family"
)

min_n_group_univariate <- 2
min_n_group_disparity <- 3
min_n_family_multivariate <- 3

jpeg_width <- 2400
jpeg_height <- 1800
jpeg_res <- 300

set.seed(123)

# Manual fixes for known label mismatches between specimen table and tree.
# Add more if needed.
tree_tip_corrections <- c(
  "Neydus" = "Nedyus",
  "Belidae" = "Agnesiotis",
  "Caridae" = "Car"
)

# ============================================================
# SECTION 2: OUTPUT FOLDER STRUCTURE
# ============================================================

subdirs <- c(
  "01_QC",
  "02_Phylogenetic_signal",
  "03_Evolutionary_models",
  "04_PGLS",
  "05_Allometry",
  "06_Group_tests",
  "07_Disparity",
  "08_ASR",
  "09_Morphospace_and_tree_plots",
  "10_Logs",
  "11_Tree_robustness"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
purrr::walk(file.path(output_dir, subdirs), ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# ============================================================
# SECTION 3: HELPER FUNCTIONS
# ============================================================

read_delim_guess <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  delim <- if (grepl(";", first_line)) ";" else ","
  df <- data.table::fread(
    path,
    sep = delim,
    dec = ",",
    encoding = "UTF-8",
    data.table = FALSE
  )
  tibble::as_tibble(df)
}

clean_names_basic <- function(x) {
  x <- gsub("^\\ufeff", "", x)
  x <- trimws(x)
  x
}

save_plot_dual <- function(plot_obj, filename_base, outdir, width = 10, height = 8, units = "in", dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(outdir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    units = units
  )
  ggplot2::ggsave(
    filename = file.path(outdir, paste0(filename_base, ".jpg")),
    plot = plot_obj,
    width = width,
    height = height,
    units = units,
    dpi = dpi
  )
}

write_clean_csv <- function(df, path) {
  readr::write_excel_csv2(df, path)
}

add_skip_log <- function(log_df, analysis_block, response, predictor, reason) {
  dplyr::bind_rows(
    log_df,
    tibble::tibble(
      analysis_block = analysis_block,
      response = response,
      predictor = predictor,
      reason = reason
    )
  )
}

mode_single_value <- function(x) {
  ux <- unique(stats::na.omit(x))
  if (length(ux) == 1) ux else NA
}

p_adjust_if_possible <- function(p) {
  if (length(p) == 0) return(numeric(0))
  if (all(is.na(p))) return(rep(NA_real_, length(p)))
  stats::p.adjust(p, method = "fdr")
}

make_tree_label <- function(family, genus_tip) {
  family <- trimws(as.character(family))
  genus_tip <- trimws(as.character(genus_tip))
  family <- gsub("^\\ufeff", "", family)
  genus_tip <- gsub("^\\ufeff", "", genus_tip)
  dplyr::if_else(
    !is.na(family) & !is.na(genus_tip) & family != "" & genus_tip != "",
    paste(family, genus_tip, sep = "___"),
    NA_character_
  )
}

safe_phylo_signal <- function(tree, trait_vec, trait_name) {
  dat <- tibble::tibble(tree_label = names(trait_vec), value = as.numeric(trait_vec)) %>%
    dplyr::filter(!is.na(value), !is.na(tree_label))

  if (nrow(dat) < 5) return(NULL)

  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
  vec <- dat$value
  names(vec) <- dat$tree_label
  vec <- vec[tr$tip.label]

  k_res <- tryCatch({
    tmp <- phytools::phylosig(tr, vec, method = "K", test = TRUE)
    tibble::tibble(
      trait = trait_name,
      method = "Blomberg_K",
      estimate = unname(tmp$K),
      statistic = NA_real_,
      p_value = unname(tmp$P)
    )
  }, error = function(e) NULL)

  lambda_res <- tryCatch({
    tmp <- phytools::phylosig(tr, vec, method = "lambda", test = TRUE)
    tibble::tibble(
      trait = trait_name,
      method = "Pagel_lambda",
      estimate = unname(tmp$lambda),
      statistic = unname(tmp$logL),
      p_value = unname(tmp$P)
    )
  }, error = function(e) NULL)

  dplyr::bind_rows(k_res, lambda_res)
}

safe_fit_continuous_models <- function(tree, trait_vec, trait_name) {
  dat <- tibble::tibble(tree_label = names(trait_vec), value = as.numeric(trait_vec)) %>%
    dplyr::filter(!is.na(value), !is.na(tree_label))

  if (nrow(dat) < 5) return(NULL)

  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
  vec <- dat$value
  names(vec) <- dat$tree_label
  vec <- vec[tr$tip.label]

  models <- c("BM", "OU", "EB")
  res <- vector("list", length(models))
  names(res) <- models

  for (m in models) {
    fit <- tryCatch(
      suppressWarnings(geiger::fitContinuous(tr, vec, model = m)),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      res[[m]] <- tibble::tibble(
        trait = trait_name,
        model = m,
        logLik = as.numeric(fit$opt$lnL),
        AIC = as.numeric(fit$opt$aic),
        AICc = as.numeric(fit$opt$aicc),
        converged = TRUE
      )
    } else {
      res[[m]] <- tibble::tibble(
        trait = trait_name,
        model = m,
        logLik = NA_real_,
        AIC = NA_real_,
        AICc = NA_real_,
        converged = FALSE
      )
    }
  }

  out <- dplyr::bind_rows(res) %>%
    dplyr::group_by(trait) %>%
    dplyr::mutate(
      delta_AICc = AICc - min(AICc, na.rm = TRUE),
      best_model = model[which.min(AICc)]
    ) %>%
    dplyr::ungroup()

  out
}

safe_pgls <- function(tree, df, response, predictor, predictor_is_factor = FALSE) {
  needed <- c("tree_label", response, predictor)
  dat <- df %>%
    dplyr::select(dplyr::all_of(needed)) %>%
    dplyr::filter(!dplyr::if_any(dplyr::all_of(needed), is.na)) %>%
    dplyr::distinct()

  if (nrow(dat) < 5) {
    return(list(result = NULL, reason = "Too few complete taxa (<5)."))
  }

  if (predictor_is_factor) {
    dat[[predictor]] <- as.factor(dat[[predictor]])
    if (nlevels(dat[[predictor]]) < 2) {
      return(list(result = NULL, reason = "Predictor has <2 levels after filtering."))
    }
    level_counts <- table(dat[[predictor]])
    if (any(level_counts < min_n_group_univariate)) {
      return(list(result = NULL, reason = "At least one group has n below threshold."))
    }
  }

  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
  dat <- dat %>% dplyr::filter(tree_label %in% tr$tip.label)
  dat <- dat[match(tr$tip.label, dat$tree_label), , drop = FALSE]
  dat <- as.data.frame(dat)

  if (anyDuplicated(dat$tree_label) > 0) {
    return(list(result = NULL, reason = "Duplicated tree_label values in model data."))
  }

  form <- stats::as.formula(paste0("`", response, "` ~ `", predictor, "`"))

  cor_struct <- tryCatch(
    ape::corPagel(value = 0.5, phy = tr, fixed = FALSE, form = ~tree_label),
    error = function(e) NULL
  )
  if (is.null(cor_struct)) {
    return(list(result = NULL, reason = "Failed to build corPagel structure."))
  }

  fit <- tryCatch(
    suppressWarnings(
      nlme::gls(
        model = form,
        data = dat,
        correlation = cor_struct,
        method = "ML",
        na.action = na.omit
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(result = NULL, reason = "gls PGLS fit failed."))
  }

  sm <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(sm) || is.null(sm$tTable)) {
    return(list(result = NULL, reason = "Could not extract GLS summary table."))
  }

  coef_tab <- as.data.frame(sm$tTable) %>%
    tibble::rownames_to_column("term") %>%
    tibble::as_tibble()

  lambda_est <- tryCatch(
    {
      coef_val <- stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE)
      as.numeric(coef_val[[1]])
    },
    error = function(e) NA_real_
  )

  anova_tab <- tryCatch(
    as.data.frame(stats::anova(fit)) %>%
      tibble::rownames_to_column("term") %>%
      tibble::as_tibble(),
    error = function(e) NULL
  )

  result <- coef_tab %>%
    dplyr::mutate(
      response = response,
      predictor = predictor,
      n_taxa = nrow(dat),
      lambda = lambda_est,
      predictor_type = ifelse(predictor_is_factor, "factor", "continuous")
    )

  list(result = result, anova = anova_tab, model = fit, reason = NA_character_)
}

safe_phyl_anova <- function(tree, df, response, group_var) {
  needed <- c("tree_label", response, group_var)
  dat <- df %>%
    dplyr::select(dplyr::all_of(needed)) %>%
    dplyr::filter(!dplyr::if_any(dplyr::all_of(needed), is.na)) %>%
    dplyr::distinct()

  if (nrow(dat) < 5) return(list(result = NULL, reason = "Too few complete taxa (<5)."))

  dat[[group_var]] <- as.factor(dat[[group_var]])
  if (nlevels(dat[[group_var]]) < 2) {
    return(list(result = NULL, reason = "Grouping variable has <2 levels."))
  }

  group_counts <- table(dat[[group_var]])
  if (any(group_counts < min_n_group_univariate)) {
    return(list(result = NULL, reason = "At least one group has n below threshold."))
  }

  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
  dat <- dat %>% dplyr::filter(tree_label %in% tr$tip.label)
  dat <- dat[match(tr$tip.label, dat$tree_label), , drop = FALSE]

  dat <- as.data.frame(dat)

  x <- dat[[response]]
  names(x) <- dat$tree_label
  g <- dat[[group_var]]
  names(g) <- dat$tree_label

  res <- tryCatch(
    suppressWarnings(phytools::phylANOVA(tr, g = g, x = x, nsim = 1000)),
    error = function(e) NULL
  )

  if (is.null(res)) return(list(result = NULL, reason = "phylANOVA failed."))

  out <- tibble::tibble(
    response = response,
    predictor = group_var,
    method = "phylANOVA",
    statistic = unname(res$F),
    p_value = unname(res$Pf),
    n_taxa = nrow(dat),
    n_groups = nlevels(g)
  )

  list(result = out, posthoc = res$Ps, reason = NA_character_)
}

safe_fastAnc <- function(tree, trait_vec, trait_name) {
  dat <- tibble::tibble(tree_label = names(trait_vec), value = as.numeric(trait_vec)) %>%
    dplyr::filter(!is.na(value), !is.na(tree_label))

  if (nrow(dat) < 5) return(NULL)

  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
  vec <- dat$value
  names(vec) <- dat$tree_label
  vec <- vec[tr$tip.label]

  anc <- tryCatch(suppressWarnings(phytools::fastAnc(tr, vec, vars = TRUE, CI = TRUE)), error = function(e) NULL)
  if (is.null(anc)) return(NULL)

  tibble::tibble(
    trait = trait_name,
    node = names(anc$ace),
    estimate = as.numeric(anc$ace),
    variance = as.numeric(anc$var),
    CI95_lower = as.numeric(anc$CI95[, 1]),
    CI95_upper = as.numeric(anc$CI95[, 2])
  )
}

standardize_tree_labels <- function(tree) {
  tree$tip.label <- trimws(tree$tip.label)
  tree$tip.label <- gsub("^\\ufeff", "", tree$tip.label)
  tree
}

safe_root_tree <- function(tree, outgroup_label) {
  tree <- standardize_tree_labels(tree)
  if (!outgroup_label %in% tree$tip.label) return(tree)
  tryCatch(
    ape::root(tree, outgroup = outgroup_label, resolve.root = TRUE),
    error = function(e) tree
  )
}

prepare_tree_for_pcm <- function(tree, target_tips, outgroup_label = NULL) {
  tree <- standardize_tree_labels(tree)
  if (!is.null(outgroup_label)) tree <- safe_root_tree(tree, outgroup_label)

  shared <- intersect(tree$tip.label, target_tips)
  if (length(shared) < 5) return(NULL)

  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, shared))
  tr <- ape::keep.tip(tr, shared)
  tr$tip.label <- trimws(tr$tip.label)
  if (is.null(tr$edge.length)) {
    tr <- tryCatch(ape::compute.brlen(tr, method = "Grafen"), error = function(e) tr)
  }
  if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0)) return(NULL)
  tr
}

build_tree_variant_set <- function(
    tree_variant_dir,
    tree_variant_files,
    reference_tree,
    target_tips,
    outgroup_label = NULL,
    topology_stress_n = 0
) {
  variants <- list()
  meta <- list()
  idx <- 1
  tree_signatures <- character()

  add_variant <- function(tree_obj, variant_id, source_path, variant_group, note = NA_character_) {
    tr <- prepare_tree_for_pcm(tree_obj, target_tips = target_tips, outgroup_label = outgroup_label)
    if (is.null(tr)) return(NULL)

    signature <- tryCatch(ape::write.tree(tr), error = function(e) NA_character_)
    if (!is.na(signature) && signature %in% tree_signatures) return(NULL)
    tree_signatures <<- c(tree_signatures, signature)

    tree_height <- tryCatch(max(ape::node.depth.edgelength(tr)), error = function(e) NA_real_)
    total_branch_length <- tryCatch(sum(tr$edge.length), error = function(e) NA_real_)
    rf_distance <- tryCatch(
      as.numeric(ape::dist.topo(ape::unroot(reference_tree), ape::unroot(tr)))[1],
      error = function(e) NA_real_
    )

    variants[[variant_id]] <<- tr
    meta[[idx]] <<- tibble::tibble(
      tree_id = variant_id,
      source_path = source_path,
      variant_group = variant_group,
      note = note,
      n_tips = ape::Ntip(tr),
      is_rooted = ape::is.rooted(tr),
      is_ultrametric = ape::is.ultrametric(tr),
      tree_height = tree_height,
      total_branch_length = total_branch_length,
      rf_distance_vs_reference = rf_distance
    )
    idx <<- idx + 1
    invisible(NULL)
  }

  add_variant(reference_tree, "working_tree", tree_file, "working", "Primary tree used in main PCM workflow.")

  for (variant_id in names(tree_variant_files)) {
    path_i <- file.path(tree_variant_dir, tree_variant_files[[variant_id]])
    if (!file.exists(path_i)) next
    tree_i <- tryCatch(ape::read.tree(path_i), error = function(e) NULL)
    if (is.null(tree_i)) next
    add_variant(tree_i, variant_id, path_i, "existing_file", "Imported from phylogeny workflow outputs.")
  }

  ref_topology <- reference_tree
  ref_topology$edge.length <- rep(1, nrow(ref_topology$edge))

  for (pw in c(0.5, 1, 2)) {
    tr_tmp <- tryCatch(
      ape::compute.brlen(ref_topology, method = "Grafen", power = pw),
      error = function(e) NULL
    )
    if (!is.null(tr_tmp)) {
      add_variant(
        tr_tmp,
        paste0("grafen_power_", gsub("\\.", "_", as.character(pw))),
        "generated_in_script",
        "generated_branch_lengths",
        paste0("Grafen branch lengths with power = ", pw)
      )
    }
  }

  if (topology_stress_n > 0) {
    for (i in seq_len(topology_stress_n)) {
      tr_tmp <- tryCatch(ape::rtree(n = ape::Ntip(reference_tree), tip.label = reference_tree$tip.label), error = function(e) NULL)
      if (is.null(tr_tmp)) next
      tr_tmp <- tryCatch(ape::compute.brlen(tr_tmp, method = "Grafen", power = 1), error = function(e) NULL)
      if (is.null(tr_tmp)) next
      add_variant(
        tr_tmp,
        paste0("topology_stress_random_", i),
        "generated_in_script",
        "topology_stress",
        "Random topology stress test using the same observed tips; interpret only as a sensitivity stress test, not as empirical tree support."
      )
    }
  }

  meta_df <- dplyr::bind_rows(meta) %>%
    dplyr::distinct(tree_id, .keep_all = TRUE)

  list(
    trees = variants,
    metadata = meta_df
  )
}

# ============================================================
# SECTION 4: DATA IMPORT
# ============================================================

tree <- ape::read.tree(tree_file)
tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- gsub("^\\ufeff", "", tree$tip.label)

pc_df <- read_delim_guess(pc_file)
specimen_key_df <- read_delim_guess(specimen_key_file)
geometry_df <- read_delim_guess(geometry_file)
centroid_df <- read_delim_guess(centroid_file)

names(pc_df) <- clean_names_basic(names(pc_df))
names(specimen_key_df) <- clean_names_basic(names(specimen_key_df))
names(geometry_df) <- clean_names_basic(names(geometry_df))
names(centroid_df) <- clean_names_basic(names(centroid_df))

pc_keep <- c("specimen_id", pc_vars_main)
pc_df <- pc_df %>% dplyr::select(dplyr::any_of(pc_keep))

centroid_df <- centroid_df %>%
  dplyr::select(dplyr::any_of(c("specimen_id", size_var))) %>%
  dplyr::distinct()

specimen_key_df <- specimen_key_df %>%
  dplyr::select(dplyr::any_of(c("specimen_id", "tree_tip", "Family", "Coxal wall hole", "Coxal Socket", "taxon_binomial"))) %>%
  dplyr::distinct()

geometry_df <- geometry_df %>%
  dplyr::select(dplyr::any_of(c("specimen_id", geometry_vars_main))) %>%
  dplyr::distinct()

for (v in c("Coxal wall hole", "Coxal Socket")) {
  if (v %in% names(specimen_key_df)) {
    specimen_key_df[[v]] <- dplyr::case_when(
      specimen_key_df[[v]] %in% c(TRUE, "TRUE", "True", "true", 1, "1") ~ TRUE,
      specimen_key_df[[v]] %in% c(FALSE, "FALSE", "False", "false", 0, "0") ~ FALSE,
      TRUE ~ NA
    )
  }
}

for (v in intersect(pc_vars_main, names(pc_df))) pc_df[[v]] <- as.numeric(pc_df[[v]])
for (v in intersect(geometry_vars_main, names(geometry_df))) geometry_df[[v]] <- as.numeric(geometry_df[[v]])
if (size_var %in% names(centroid_df)) centroid_df[[size_var]] <- as.numeric(centroid_df[[size_var]])

# ============================================================
# SECTION 5: SPECIMEN-LEVEL MERGE
# ============================================================

specimen_merged <- specimen_key_df %>%
  dplyr::left_join(pc_df, by = "specimen_id") %>%
  dplyr::left_join(geometry_df, by = "specimen_id") %>%
  dplyr::left_join(centroid_df, by = "specimen_id")

specimen_merged <- specimen_merged %>%
  dplyr::mutate(
    tree_tip = trimws(tree_tip),
    Family = trimws(Family),
    tree_tip = dplyr::recode(tree_tip, !!!tree_tip_corrections),
    tree_label = make_tree_label(Family, tree_tip)
  )

merge_summary <- tibble::tibble(
  metric = c(
    "tree_tip_count_in_tree",
    "rows_specimen_key",
    "rows_pc",
    "rows_geometry",
    "rows_centroid",
    "rows_specimen_merged",
    "unique_tree_tip_values_in_specimen_data",
    "unique_tree_label_values_in_specimen_data"
  ),
  value = c(
    length(tree$tip.label),
    nrow(specimen_key_df),
    nrow(pc_df),
    nrow(geometry_df),
    nrow(centroid_df),
    nrow(specimen_merged),
    dplyr::n_distinct(specimen_merged$tree_tip),
    dplyr::n_distinct(specimen_merged$tree_label)
  )
)

write_clean_csv(merge_summary, file.path(output_dir, "01_QC", "qc_merge_summary.csv"))

# ============================================================
# SECTION 6: TIP-LEVEL AGGREGATION
# ============================================================

continuous_vars <- c(pc_vars_main, geometry_vars_main, size_var)
continuous_vars <- continuous_vars[continuous_vars %in% names(specimen_merged)]

tip_counts <- specimen_merged %>%
  dplyr::filter(!is.na(tree_label)) %>%
  dplyr::count(tree_label, name = "n_specimens")

discrete_consistency <- specimen_merged %>%
  dplyr::filter(!is.na(tree_label)) %>%
  dplyr::group_by(tree_label) %>%
  dplyr::summarise(
    tree_tip_values = paste(sort(unique(stats::na.omit(tree_tip))), collapse = " | "),
    Family_values = paste(sort(unique(stats::na.omit(Family))), collapse = " | "),
    Family_consistent = dplyr::n_distinct(stats::na.omit(Family)) <= 1,
    Coxal_wall_hole_values = paste(sort(unique(stats::na.omit(`Coxal wall hole`))), collapse = " | "),
    Coxal_wall_hole_consistent = dplyr::n_distinct(stats::na.omit(`Coxal wall hole`)) <= 1,
    Coxal_Socket_values = paste(sort(unique(stats::na.omit(`Coxal Socket`))), collapse = " | "),
    Coxal_Socket_consistent = dplyr::n_distinct(stats::na.omit(`Coxal Socket`)) <= 1,
    .groups = "drop"
  )

tip_level_df <- specimen_merged %>%
  dplyr::filter(!is.na(tree_label)) %>%
  dplyr::group_by(tree_label) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(continuous_vars), ~ mean(.x, na.rm = TRUE)),
    tree_tip = mode_single_value(tree_tip),
    Family = mode_single_value(Family),
    `Coxal wall hole` = mode_single_value(`Coxal wall hole`),
    `Coxal Socket` = mode_single_value(`Coxal Socket`),
    n_specimens = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(continuous_vars), ~ ifelse(is.nan(.x), NA, .x))
  )

write_clean_csv(tip_counts, file.path(output_dir, "01_QC", "qc_specimens_per_tree_label.csv"))
write_clean_csv(discrete_consistency, file.path(output_dir, "01_QC", "qc_tree_label_aggregation_summary.csv"))

# ============================================================
# SECTION 7: TREE MATCHING
# ============================================================

tip_level_df$tree_label <- trimws(tip_level_df$tree_label)
tip_level_df$tree_label <- gsub("^\\ufeff", "", tip_level_df$tree_label)

shared_tips <- intersect(tree$tip.label, tip_level_df$tree_label)

if (length(shared_tips) == 0) {
  stop(
    "No shared tree tips found between tree$tip.label and tip_level_df$tree_label.\n",
    "Check Family/tree_tip formatting and spelling."
  )
}

tips_in_tree_not_data <- setdiff(tree$tip.label, tip_level_df$tree_label)
tips_in_data_not_tree <- setdiff(tip_level_df$tree_label, tree$tip.label)

tree_pruned <- ape::drop.tip(tree, tips_in_tree_not_data)

tip_level_df <- tip_level_df %>%
  dplyr::filter(tree_label %in% shared_tips) %>%
  dplyr::arrange(match(tree_label, tree_pruned$tip.label))

tree_match_summary <- tibble::tibble(
  metric = c(
    "tips_in_tree",
    "tips_in_tip_level_data_before_match",
    "shared_tips",
    "tips_in_tree_not_data",
    "tips_in_data_not_tree"
  ),
  value = c(
    length(tree$tip.label),
    length(unique(specimen_merged$tree_label)),
    length(shared_tips),
    length(tips_in_tree_not_data),
    length(tips_in_data_not_tree)
  )
)

write_clean_csv(tree_match_summary, file.path(output_dir, "01_QC", "qc_tree_matching_summary.csv"))
write_clean_csv(tibble::tibble(tree_label = tips_in_tree_not_data), file.path(output_dir, "01_QC", "qc_tree_tips_without_data.csv"))
write_clean_csv(tibble::tibble(tree_label = tips_in_data_not_tree), file.path(output_dir, "01_QC", "qc_data_tips_without_tree.csv"))

tree_variant_bundle <- build_tree_variant_set(
  tree_variant_dir = tree_variant_dir,
  tree_variant_files = tree_variant_files,
  reference_tree = tree_pruned,
  target_tips = tip_level_df$tree_label,
  outgroup_label = tree_outgroup_label,
  topology_stress_n = tree_topology_stress_n
)

tree_variant_list <- tree_variant_bundle$trees
tree_variant_metadata <- tree_variant_bundle$metadata

if (is.data.frame(tree_variant_metadata) && nrow(tree_variant_metadata) > 0) {
  write_clean_csv(tree_variant_metadata, file.path(output_dir, "01_QC", "qc_tree_variant_metadata.csv"))
}

# ============================================================
# SECTION 8: QC TABLES
# ============================================================

missingness_tbl <- tibble::tibble(
  variable = names(tip_level_df),
  n_missing = sapply(tip_level_df, function(x) sum(is.na(x))),
  prop_missing = sapply(tip_level_df, function(x) mean(is.na(x)))
)

group_counts_tbl <- dplyr::bind_rows(
  tip_level_df %>%
    dplyr::filter(!is.na(Family)) %>%
    dplyr::count(Family, name = "n") %>%
    dplyr::transmute(variable = "Family", level = as.character(Family), n = n),
  tip_level_df %>%
    dplyr::filter(!is.na(`Coxal wall hole`)) %>%
    dplyr::count(`Coxal wall hole`, name = "n") %>%
    dplyr::transmute(variable = "Coxal wall hole", level = as.character(`Coxal wall hole`), n = n),
  tip_level_df %>%
    dplyr::filter(!is.na(`Coxal Socket`)) %>%
    dplyr::count(`Coxal Socket`, name = "n") %>%
    dplyr::transmute(variable = "Coxal Socket", level = as.character(`Coxal Socket`), n = n)
)

write_clean_csv(missingness_tbl, file.path(output_dir, "01_QC", "qc_missingness.csv"))
write_clean_csv(group_counts_tbl, file.path(output_dir, "01_QC", "qc_group_counts.csv"))

# ============================================================
# SECTION 9: QC PLOTS
# ============================================================

qc_cont_vars <- c(pc_vars_main, geometry_vars_main, size_var)
qc_cont_vars <- qc_cont_vars[qc_cont_vars %in% names(tip_level_df)]

qc_long <- tip_level_df %>%
  dplyr::select(tree_label, dplyr::all_of(qc_cont_vars)) %>%
  tidyr::pivot_longer(-tree_label, names_to = "variable", values_to = "value")

p_hist <- ggplot2::ggplot(qc_long, ggplot2::aes(x = value)) +
  ggplot2::geom_histogram(bins = 20) +
  ggplot2::facet_wrap(~ variable, scales = "free") +
  ggplot2::labs(
    title = "Distributions of continuous variables (tip level)",
    x = "Value",
    y = "Count"
  )

save_plot_dual(p_hist, "qc_histograms_continuous", file.path(output_dir, "01_QC"), width = 12, height = 10)

if (length(qc_cont_vars) >= 2) {
  cor_tbl <- tip_level_df %>%
    dplyr::select(dplyr::all_of(qc_cont_vars)) %>%
    stats::cor(use = "pairwise.complete.obs", method = "pearson")

  cor_long <- as.data.frame(as.table(cor_tbl)) %>%
    tibble::as_tibble() %>%
    dplyr::rename(var1 = Var1, var2 = Var2, correlation = Freq)

  p_cor <- ggplot2::ggplot(cor_long, ggplot2::aes(var1, var2, fill = correlation)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = round(correlation, 2)), size = 3) +
    ggplot2::scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
    ggplot2::coord_equal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Correlation matrix of continuous variables", x = NULL, y = NULL)

  save_plot_dual(p_cor, "qc_correlation_matrix", file.path(output_dir, "01_QC"), width = 12, height = 10)
  write_clean_csv(cor_long, file.path(output_dir, "01_QC", "qc_correlation_matrix.csv"))
}

# ============================================================
# SECTION 10: PHYLOGENETIC SIGNAL
# ============================================================

phylo_signal_results <- purrr::map_dfr(qc_cont_vars, function(v) {
  vec <- tip_level_df[[v]]
  names(vec) <- tip_level_df$tree_label
  safe_phylo_signal(tree_pruned, vec, v)
})

if (nrow(phylo_signal_results) > 0) {
  phylo_signal_results <- phylo_signal_results %>%
    dplyr::mutate(fdr_p_value = p_adjust_if_possible(p_value))
}
write_clean_csv(phylo_signal_results, file.path(output_dir, "02_Phylogenetic_signal", "phylogenetic_signal_continuous.csv"))

key_trait_maps <- intersect(c("PC1", "abs_winding_angle_deg", "axial_span", "centroid_size"), qc_cont_vars)
for (v in key_trait_maps) {
  dat <- tip_level_df %>%
    dplyr::select(tree_label, dplyr::all_of(v)) %>%
    dplyr::filter(!is.na(.data[[v]]))

  if (nrow(dat) >= 5) {
    tr <- ape::drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, dat$tree_label))
    vec <- dat[[v]]
    names(vec) <- dat$tree_label
    vec <- vec[tr$tip.label]

    grDevices::pdf(file.path(output_dir, "02_Phylogenetic_signal", paste0("trait_map_", v, ".pdf")), width = 10, height = 10)
    try({
      phytools::contMap(tr, vec, plot = TRUE, legend = 0.7, fsize = 0.8)
      title(main = paste("Trait map:", v))
    }, silent = TRUE)
    grDevices::dev.off()

    grDevices::jpeg(file.path(output_dir, "02_Phylogenetic_signal", paste0("trait_map_", v, ".jpg")),
                    width = jpeg_width, height = jpeg_height, res = jpeg_res)
    try({
      phytools::contMap(tr, vec, plot = TRUE, legend = 0.7, fsize = 0.8)
      title(main = paste("Trait map:", v))
    }, silent = TRUE)
    grDevices::dev.off()
  }
}

# ============================================================
# SECTION 11: UNIVARIATE EVOLUTIONARY MODEL FITS
# ============================================================

evol_model_results <- purrr::map_dfr(qc_cont_vars, function(v) {
  vec <- tip_level_df[[v]]
  names(vec) <- tip_level_df$tree_label
  safe_fit_continuous_models(tree_pruned, vec, v)
})

write_clean_csv(evol_model_results, file.path(output_dir, "03_Evolutionary_models", "evolutionary_model_fits_univariate.csv"))

# ============================================================
# SECTION 12: PGLS - CONTINUOUS VS CONTINUOUS
# ============================================================

skip_log <- tibble::tibble(
  analysis_block = character(),
  response = character(),
  predictor = character(),
  reason = character()
)

geometry_pairs <- combn(geometry_vars_main, 2, simplify = FALSE)
geometry_models <- purrr::map(geometry_pairs, ~ tibble::tibble(response = .x[1], predictor = .x[2])) %>%
  dplyr::bind_rows()

shape_geometry_models <- tidyr::expand_grid(response = pc_vars_main, predictor = geometry_vars_main) %>%
  dplyr::filter(response %in% names(tip_level_df), predictor %in% names(tip_level_df))

size_models <- dplyr::bind_rows(
  tibble::tibble(response = pc_vars_main, predictor = size_var),
  tibble::tibble(response = geometry_vars_main, predictor = size_var)
) %>%
  dplyr::filter(response %in% names(tip_level_df), predictor %in% names(tip_level_df))

all_cont_models <- dplyr::bind_rows(
  geometry_models,
  shape_geometry_models,
  size_models
)

pgls_results <- list()
pgls_anova_results <- list()
pgls_model_index <- 1

for (i in seq_len(nrow(all_cont_models))) {
  resp <- all_cont_models$response[i]
  pred <- all_cont_models$predictor[i]

  tmp <- safe_pgls(tree_pruned, tip_level_df, resp, pred, predictor_is_factor = FALSE)

  if (is.null(tmp$result)) {
    skip_log <- add_skip_log(skip_log, "PGLS_continuous", resp, pred, tmp$reason)
  } else {
    pgls_results[[pgls_model_index]] <- tmp$result
    if (!is.null(tmp$anova)) {
      pgls_anova_results[[pgls_model_index]] <- tmp$anova %>%
        dplyr::mutate(response = resp, predictor = pred)
    }

    plot_df <- tip_level_df %>%
      dplyr::select(tree_label, dplyr::all_of(resp), dplyr::all_of(pred)) %>%
      dplyr::filter(!dplyr::if_any(c(dplyr::all_of(resp), dplyr::all_of(pred)), is.na))

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[pred]], y = .data[[resp]], label = tree_label)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_smooth(method = "lm", se = TRUE) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
      ggplot2::labs(
        title = paste(resp, "vs", pred, "(tip level)"),
        subtitle = "Visual trend only; inferential result is the PGLS model",
        x = pred,
        y = resp
      )

    save_plot_dual(
      p,
      paste0("pgls_scatter_", resp, "_vs_", pred),
      file.path(output_dir, "04_PGLS"),
      width = 10,
      height = 8
    )

    pgls_model_index <- pgls_model_index + 1
  }
}

pgls_results_df <- dplyr::bind_rows(pgls_results)

if (nrow(pgls_results_df) > 0) {
  pgls_results_df <- pgls_results_df %>%
    dplyr::rename_with(~ "estimate", dplyr::any_of(c("Estimate", "Value"))) %>%
    dplyr::rename_with(~ "std_error", dplyr::any_of(c("Std. Error", "Std.Error"))) %>%
    dplyr::rename_with(~ "statistic", dplyr::any_of(c("t value", "t-value"))) %>%
    dplyr::rename_with(~ "p_value", dplyr::any_of(c("Pr(>|t|)", "p-value"))) %>%
    dplyr::mutate(
      fdr_p_value = if ("p_value" %in% names(.)) p_adjust_if_possible(p_value) else NA_real_
    )
}
write_clean_csv(pgls_results_df, file.path(output_dir, "04_PGLS", "pgls_continuous_vs_continuous.csv"))

if (length(pgls_anova_results) > 0) {
  pgls_anova_df <- dplyr::bind_rows(pgls_anova_results)
  write_clean_csv(pgls_anova_df, file.path(output_dir, "04_PGLS", "pgls_continuous_vs_continuous_anova.csv"))
}

# ============================================================
# SECTION 13: ALLOMETRY BLOCK
# ============================================================

allometry_df <- if (nrow(pgls_results_df) > 0 && "predictor" %in% names(pgls_results_df)) {
  pgls_results_df %>% dplyr::filter(predictor == size_var)
} else {
  tibble::tibble()
}

write_clean_csv(allometry_df, file.path(output_dir, "05_Allometry", "allometry_results.csv"))

# ============================================================
# SECTION 14: GROUP TESTS (DISCRETE VARIABLES)
# ============================================================

group_test_vars <- c("Coxal wall hole", "Coxal Socket", "Family")
group_test_vars <- group_test_vars[group_test_vars %in% names(tip_level_df)]

format_group_label <- function(x) {
  dplyr::case_when(
    x == "Coxal wall hole" ~ "coxal wall opening",
    x == "Coxal Socket" ~ "coxal socket",
    TRUE ~ x
  )
}

group_res_list <- list()
group_model_index <- 1

for (gvar in group_test_vars) {
  gvar_label <- format_group_label(gvar)

  for (resp in qc_cont_vars) {
    dat_tmp <- tip_level_df %>%
      dplyr::select(tree_label, dplyr::all_of(resp), dplyr::all_of(gvar)) %>%
      dplyr::filter(!dplyr::if_any(dplyr::all_of(c(resp, gvar)), is.na))

    if (nrow(dat_tmp) < 5) {
      skip_log <- add_skip_log(skip_log, "Group_test", resp, gvar, "Too few complete taxa (<5).")
      next
    }

    counts <- table(dat_tmp[[gvar]])
    if (length(counts) < 2) {
      skip_log <- add_skip_log(skip_log, "Group_test", resp, gvar, "Grouping variable has <2 levels.")
      next
    }
    if (any(counts < min_n_group_univariate)) {
      skip_log <- add_skip_log(skip_log, "Group_test", resp, gvar, "At least one group has n below threshold.")
      next
    }

    pgls_tmp <- safe_pgls(tree_pruned, tip_level_df, resp, gvar, predictor_is_factor = TRUE)
    if (is.null(pgls_tmp$result)) {
      skip_log <- add_skip_log(skip_log, "Group_test_PGLS", resp, gvar, pgls_tmp$reason)
    } else {
      group_res_list[[group_model_index]] <- pgls_tmp$result %>%
        dplyr::mutate(method = "factor_PGLS")
      group_model_index <- group_model_index + 1
    }

    phyl_tmp <- safe_phyl_anova(tree_pruned, tip_level_df, resp, gvar)
    if (is.null(phyl_tmp$result)) {
      skip_log <- add_skip_log(skip_log, "Group_test_phylANOVA", resp, gvar, phyl_tmp$reason)
    } else {
      group_res_list[[group_model_index]] <- phyl_tmp$result
      group_model_index <- group_model_index + 1
    }

    dat_plot <- dat_tmp
    if (gvar == "Coxal wall hole") {
      dat_plot[[gvar]] <- factor(dat_plot[[gvar]], levels = c(FALSE, TRUE), labels = c("absent", "present"))
    }

    p <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = .data[[gvar]], y = .data[[resp]], fill = .data[[gvar]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.9) +
      ggplot2::labs(
        title = paste(resp, "by", gvar_label),
        subtitle = "Visual group plot; inferential output is exported separately",
        x = gvar_label,
        y = resp
      ) +
      ggplot2::theme(legend.position = "none")

    save_plot_dual(
      p,
      paste0("groupplot_", resp, "_by_", make.names(gvar)),
      file.path(output_dir, "06_Group_tests"),
      width = 10,
      height = 8
    )
  }
}

group_test_results_df <- dplyr::bind_rows(group_res_list)
if (nrow(group_test_results_df) > 0 && "p_value" %in% names(group_test_results_df)) {
  group_test_results_df <- group_test_results_df %>%
    dplyr::mutate(fdr_p_value = p_adjust_if_possible(p_value))
}

write_clean_csv(group_test_results_df, file.path(output_dir, "06_Group_tests", "group_tests_results.csv"))

# ============================================================
# SECTION 15: DISPARITY ANALYSES
# ============================================================

pc_for_disparity <- pc_vars_main[pc_vars_main %in% names(tip_level_df)]
disparity_group_vars <- c("Family", "Coxal wall hole", "Coxal Socket")
disparity_group_vars <- disparity_group_vars[disparity_group_vars %in% names(tip_level_df)]

disparity_results <- list()
disp_i <- 1

for (gvar in disparity_group_vars) {
  dat <- tip_level_df %>%
    dplyr::select(tree_label, dplyr::all_of(gvar), dplyr::all_of(pc_for_disparity)) %>%
    dplyr::filter(!dplyr::if_any(dplyr::all_of(c(gvar, pc_for_disparity)), is.na))

  if (nrow(dat) < 5) {
    skip_log <- add_skip_log(skip_log, "Disparity", "PC_space", gvar, "Too few complete taxa.")
    next
  }

  counts <- table(dat[[gvar]])
  keep_levels <- names(counts[counts >= min_n_group_disparity])
  dat <- dat %>% dplyr::filter(.data[[gvar]] %in% keep_levels)

  if (length(unique(dat[[gvar]])) < 2) {
    skip_log <- add_skip_log(skip_log, "Disparity", "PC_space", gvar, "Fewer than 2 groups with n >= disparity threshold.")
    next
  }

  group_split <- split(dat, dat[[gvar]])
  disp_tbl <- purrr::map_dfr(names(group_split), function(grp) {
    sub <- group_split[[grp]]
    mat <- as.matrix(sub[, pc_for_disparity, drop = FALSE])
    center <- colMeans(mat, na.rm = TRUE)
    d <- sqrt(rowSums((mat - matrix(center, nrow(mat), ncol(mat), byrow = TRUE))^2))
    tibble::tibble(
      group_var = gvar,
      group = grp,
      n_taxa = nrow(sub),
      disparity_mean_distance = mean(d, na.rm = TRUE),
      disparity_sd_distance = stats::sd(d, na.rm = TRUE)
    )
  })

  disparity_results[[disp_i]] <- disp_tbl
  disp_i <- disp_i + 1

  p <- ggplot2::ggplot(disp_tbl, ggplot2::aes(x = reorder(group, disparity_mean_distance), y = disparity_mean_distance)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = disparity_mean_distance - disparity_sd_distance,
        ymax = disparity_mean_distance + disparity_sd_distance
      ),
      width = 0.2
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Morphospace disparity by", gvar),
      x = gvar,
      y = "Mean distance to group centroid"
    )

  save_plot_dual(
    p,
    paste0("disparity_", make.names(gvar)),
    file.path(output_dir, "07_Disparity"),
    width = 10,
    height = 8
  )
}

disparity_results_df <- dplyr::bind_rows(disparity_results)
write_clean_csv(disparity_results_df, file.path(output_dir, "07_Disparity", "disparity_results.csv"))

# ============================================================
# SECTION 16: MULTIVARIATE EVOLUTIONARY MODEL FITS
# ============================================================

mv_model_results <- tibble::tibble()

mv_pc_sets <- list(
  PC1_5 = pc_vars_main[pc_vars_main %in% names(tip_level_df)],
  PC1_4 = intersect(paste0("PC", 1:4), names(tip_level_df))
)

for (set_name in names(mv_pc_sets)) {
  pcs <- mv_pc_sets[[set_name]]
  if (length(pcs) < 2) {
    skip_log <- add_skip_log(skip_log, "Multivariate_Evol", set_name, "BM/OU/EB", "Too few PCs available.")
    next
  }

  dat <- tip_level_df %>%
    dplyr::select(tree_label, dplyr::all_of(pcs)) %>%
    dplyr::filter(!dplyr::if_any(dplyr::all_of(pcs), is.na))

  if (nrow(dat) < 6) {
    skip_log <- add_skip_log(skip_log, "Multivariate_Evol", set_name, "BM/OU/EB", "Too few complete taxa.")
    next
  }

  tr <- ape::drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, dat$tree_label))
  dat <- dat[match(tr$tip.label, dat$tree_label), , drop = FALSE]
  dat_df <- as.data.frame(dat)
  Y <- as.matrix(dat_df[, pcs, drop = FALSE])

  fit_mv_model <- function(model_name) {
    tryCatch({
      fit <- switch(
        model_name,
        "BM" = suppressWarnings(mvMORPH::mvBM(tree = tr, data = Y, model = "BM1")),
        "OU" = suppressWarnings(mvMORPH::mvOU(tree = tr, data = Y, model = "OU1")),
        "EB" = suppressWarnings(mvMORPH::mvEB(tree = tr, data = Y, model = "EB"))
      )
      tibble::tibble(
        pc_set = set_name,
        model = model_name,
        logLik = as.numeric(stats::logLik(fit)),
        AIC = stats::AIC(fit),
        converged = TRUE
      )
    }, error = function(e) {
      skip_log <<- add_skip_log(skip_log, "Multivariate_Evol", set_name, model_name, "Model failed or unstable.")
      tibble::tibble(
        pc_set = set_name,
        model = model_name,
        logLik = NA_real_,
        AIC = NA_real_,
        converged = FALSE
      )
    })
  }

  mv_tmp <- dplyr::bind_rows(
    fit_mv_model("BM"),
    fit_mv_model("OU"),
    fit_mv_model("EB")
  ) %>%
    dplyr::group_by(pc_set) %>%
    dplyr::mutate(
      delta_AIC = AIC - min(AIC, na.rm = TRUE),
      best_model = model[which.min(AIC)]
    ) %>%
    dplyr::ungroup()

  mv_model_results <- dplyr::bind_rows(mv_model_results, mv_tmp)
}

write_clean_csv(mv_model_results, file.path(output_dir, "03_Evolutionary_models", "evolutionary_model_fits_multivariate.csv"))

# ============================================================
# SECTION 17: ANCESTRAL STATE RECONSTRUCTION (EXPLORATORY)
# ============================================================

asr_cont_vars <- intersect(c("abs_winding_angle_deg", "axial_span", "PC1", "centroid_size"), names(tip_level_df))
asr_cont_results <- purrr::map_dfr(asr_cont_vars, function(v) {
  vec <- tip_level_df[[v]]
  names(vec) <- tip_level_df$tree_label
  safe_fastAnc(tree_pruned, vec, v)
})
write_clean_csv(asr_cont_results, file.path(output_dir, "08_ASR", "asr_continuous_fastAnc.csv"))

discrete_asr_vars <- intersect(c("Coxal wall hole", "Coxal Socket"), names(tip_level_df))
discrete_asr_summary <- list()
dai <- 1

for (v in discrete_asr_vars) {
  dat <- tip_level_df %>%
    dplyr::select(tree_label, dplyr::all_of(v)) %>%
    dplyr::filter(!is.na(.data[[v]]))

  if (nrow(dat) < 5) {
    skip_log <- add_skip_log(skip_log, "ASR_discrete", v, "simmap", "Too few taxa with non-missing states.")
    next
  }
  if (dplyr::n_distinct(dat[[v]]) < 2) {
    skip_log <- add_skip_log(skip_log, "ASR_discrete", v, "simmap", "Trait has fewer than 2 states.")
    next
  }

  tr <- ape::drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, dat$tree_label))
  dat_df <- as.data.frame(dat)
  states <- dat_df[[v]]
  names(states) <- dat_df$tree_label
  states <- states[tr$tip.label]

  smap <- tryCatch(
    suppressWarnings(phytools::make.simmap(tr, states, model = "ER", nsim = 50)),
    error = function(e) NULL
  )

  if (is.null(smap)) {
    skip_log <- add_skip_log(skip_log, "ASR_discrete", v, "simmap", "make.simmap failed.")
    next
  }

  summary_obj <- summary(smap)
  if (!is.null(summary_obj$ace)) {
    tmp <- as.data.frame(summary_obj$ace) %>%
      tibble::rownames_to_column("node") %>%
      tidyr::pivot_longer(-node, names_to = "state", values_to = "posterior_probability") %>%
      dplyr::mutate(trait = v)

    discrete_asr_summary[[dai]] <- tmp
    dai <- dai + 1
  }

  grDevices::pdf(file.path(output_dir, "08_ASR", paste0("simmap_", make.names(v), ".pdf")), width = 10, height = 10)
  try({
    plot(summary_obj, fsize = 0.8)
    title(main = paste("Stochastic character mapping:", v))
  }, silent = TRUE)
  grDevices::dev.off()

  grDevices::jpeg(file.path(output_dir, "08_ASR", paste0("simmap_", make.names(v), ".jpg")),
                  width = jpeg_width, height = jpeg_height, res = jpeg_res)
  try({
    plot(summary_obj, fsize = 0.8)
    title(main = paste("Stochastic character mapping:", v))
  }, silent = TRUE)
  grDevices::dev.off()
}

discrete_asr_summary_df <- dplyr::bind_rows(discrete_asr_summary)
write_clean_csv(discrete_asr_summary_df, file.path(output_dir, "08_ASR", "asr_discrete_simmap_summary.csv"))

# ============================================================
# SECTION 18: MORPHOSPACE AND TREE PLOTS
# ============================================================

pc_plot_pairs <- list(
  c("PC1", "PC2"),
  c("PC1", "PC3"),
  c("PC2", "PC3")
)

family_n <- tip_level_df %>%
  dplyr::filter(!is.na(Family)) %>%
  dplyr::count(Family, sort = TRUE)

family_levels <- family_n$Family
n_fams <- length(family_levels)
palette_use <- if (n_fams <= 12) {
  RColorBrewer::brewer.pal(max(3, min(12, n_fams)), "Set3")[seq_len(n_fams)]
} else {
  scales::hue_pal()(n_fams)
}
family_colors <- setNames(palette_use, family_levels)

for (pair in pc_plot_pairs) {
  if (!all(pair %in% names(tip_level_df))) next
  dat <- tip_level_df %>%
    dplyr::select(tree_label, Family, `Coxal wall hole`, dplyr::all_of(pair)) %>%
    dplyr::filter(!dplyr::if_any(dplyr::all_of(pair), is.na))

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data[[pair[1]]],
      y = .data[[pair[2]]],
      color = Family,
      shape = `Coxal wall hole`,
      label = tree_label
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    ggplot2::scale_color_manual(values = family_colors, na.translate = TRUE) +
    ggplot2::labs(
      title = paste(pair[1], "vs", pair[2], "(tip-level morphospace)"),
      subtitle = "Color = Family, shape = Coxal wall hole",
      x = pair[1],
      y = pair[2]
    )

  save_plot_dual(
    p,
    paste0("morphospace_", pair[1], "_vs_", pair[2]),
    file.path(output_dir, "09_Morphospace_and_tree_plots"),
    width = 10,
    height = 8
  )
}

if (all(c("PC1", "PC2") %in% names(tip_level_df))) {
  dat <- tip_level_df %>%
    dplyr::select(tree_label, PC1, PC2) %>%
    dplyr::filter(!dplyr::if_any(dplyr::all_of(c("PC1", "PC2")), is.na))

  if (nrow(dat) >= 5) {
    dat_df <- as.data.frame(dat)
    tr <- ape::drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, dat_df$tree_label))
    xy_df <- dat_df[match(tr$tip.label, dat_df$tree_label), , drop = FALSE]
    xy <- as.matrix(xy_df[, c("PC1", "PC2"), drop = FALSE])
    rownames(xy) <- xy_df$tree_label

    grDevices::pdf(file.path(output_dir, "09_Morphospace_and_tree_plots", "phylomorphospace_PC1_PC2.pdf"), width = 10, height = 8)
    try({
      phytools::phylomorphospace(tr, xy, label = "off", xlab = "PC1", ylab = "PC2")
      graphics::points(xy[, 1], xy[, 2], pch = 19, cex = 1.2)
      graphics::text(xy[, 1], xy[, 2], labels = rownames(xy), pos = 3, cex = 0.7)
      graphics::title(main = "Phylomorphospace: PC1 vs PC2")
    }, silent = TRUE)
    grDevices::dev.off()

    grDevices::jpeg(file.path(output_dir, "09_Morphospace_and_tree_plots", "phylomorphospace_PC1_PC2.jpg"),
                    width = jpeg_width, height = jpeg_height, res = jpeg_res)
    try({
      phytools::phylomorphospace(tr, xy, label = "off", xlab = "PC1", ylab = "PC2")
      graphics::points(xy[, 1], xy[, 2], pch = 19, cex = 1.2)
      graphics::text(xy[, 1], xy[, 2], labels = rownames(xy), pos = 3, cex = 0.7)
      graphics::title(main = "Phylomorphospace: PC1 vs PC2")
    }, silent = TRUE)
    grDevices::dev.off()
  }
}

# ============================================================
# SECTION 19: FINAL LOG EXPORT
# ============================================================

write_clean_csv(tip_level_df, file.path(output_dir, "10_Logs", "tip_level_dataset_used_for_PCM.csv"))
write_clean_csv(skip_log, file.path(output_dir, "10_Logs", "skipped_analyses_log.csv"))

# ============================================================
# SECTION 20: CONSOLE SUMMARY
# ============================================================

message("============================================================")
message("PHYLOGENETIC COMPARATIVE ANALYSIS FINISHED")
message("Results written to: ", output_dir)
message("Shared tree tips analysed: ", nrow(tip_level_df))
message("Continuous variables analysed: ", paste(qc_cont_vars, collapse = ", "))
message("Grouping variables considered: ", paste(group_test_vars, collapse = ", "))
message("IMPORTANT: Treat all PCM results as exploratory trend analyses.")
message("============================================================")



# ============================================================
# SECTION 21: FIGURE - PHYLOGENETIC SIGNAL (CLEAN GROUPED BARPLOT)
# ============================================================

message("SECTION 21: Building clean grouped phylogenetic signal barplot ...")

phylo_signal_plot_dir <- file.path(output_dir, "09_Morphospace_and_tree_plots")
dir.create(phylo_signal_plot_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1) Daten holen
# ------------------------------------------------------------
phylo_signal_plot_df <- NULL

if (exists("phylo_signal_results") && is.data.frame(phylo_signal_results) && nrow(phylo_signal_results) > 0) {
  phylo_signal_plot_df <- phylo_signal_results
} else {
  phylo_signal_csv <- file.path(output_dir, "02_Phylogenetic_signal", "phylogenetic_signal_continuous.csv")
  if (file.exists(phylo_signal_csv)) {
    phylo_signal_plot_df <- read_delim_guess(phylo_signal_csv)
  } else {
    stop("No phylogenetic signal data found.")
  }
}

required_cols <- c("trait", "method", "estimate", "p_value")
missing_cols <- setdiff(required_cols, names(phylo_signal_plot_df))
if (length(missing_cols) > 0) {
  stop(paste("Phylogenetic signal table is missing required columns:", paste(missing_cols, collapse = ", ")))
}

if (!"fdr_p_value" %in% names(phylo_signal_plot_df)) {
  phylo_signal_plot_df <- phylo_signal_plot_df %>%
    dplyr::mutate(fdr_p_value = p_adjust_if_possible(p_value))
}

# ------------------------------------------------------------
# 2) Nur Main-Text-Traits
# ------------------------------------------------------------
traits_keep <- c(
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "abs_winding_angle_deg",
  "axial_span"
)

trait_label_lookup <- c(
  "PC1" = "PC1",
  "PC2" = "PC2",
  "PC3" = "PC3",
  "PC4" = "PC4",
  "PC5" = "PC5",
  "abs_winding_angle_deg" = "Absolute winding angle",
  "axial_span" = "Axial span"
)

trait_order <- c(
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "abs_winding_angle_deg", "axial_span"
)

phylo_signal_plot_df <- phylo_signal_plot_df %>%
  dplyr::filter(
    !is.na(trait),
    !is.na(method),
    !is.na(estimate),
    trait %in% traits_keep
  ) %>%
  dplyr::mutate(
    trait = as.character(trait),
    trait_label = unname(trait_label_lookup[trait]),
    method = dplyr::case_when(
      method == "Blomberg_K" ~ "Blomberg's K",
      method == "Pagel_lambda" ~ "Pagel's lambda",
      TRUE ~ as.character(method)
    ),
    signal_status = dplyr::case_when(
      !is.na(fdr_p_value) & fdr_p_value < 0.05 ~ "Signal detected",
      TRUE ~ "No signal detected"
    ),
    sig_label = dplyr::case_when(
      !is.na(fdr_p_value) & fdr_p_value < 0.001 ~ "***",
      !is.na(fdr_p_value) & fdr_p_value < 0.01  ~ "**",
      !is.na(fdr_p_value) & fdr_p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::mutate(
    trait = factor(trait, levels = rev(trait_order)),
    trait_label = factor(trait_label, levels = rev(unname(trait_label_lookup[trait_order]))),
    method = factor(method, levels = c("Blomberg's K", "Pagel's lambda")),
    signal_status = factor(signal_status, levels = c("No signal detected", "Signal detected"))
  ) %>%
  dplyr::mutate(
    y_num = as.numeric(trait_label),
    y_text = dplyr::if_else(method == "Blomberg's K", y_num - 0.18, y_num + 0.18)
  )

if (nrow(phylo_signal_plot_df) == 0) {
  stop("Phylogenetic signal table is empty after filtering.")
}

# ------------------------------------------------------------
# 3) Farben + Positionen
# ------------------------------------------------------------
method_colors <- c(
  "Blomberg's K" = "#3B82F6",
  "Pagel's lambda" = "#E68613"
)

bar_dodge <- ggplot2::position_dodge(width = 0.72)

x_max <- max(phylo_signal_plot_df$estimate, na.rm = TRUE)
x_upper <- max(1.2, x_max * 1.12)

# ------------------------------------------------------------
# 4) Plot
# ------------------------------------------------------------
p_phylo_signal_bar <- ggplot2::ggplot(
  phylo_signal_plot_df,
  ggplot2::aes(x = estimate, y = trait_label, fill = method, alpha = signal_status)
) +
  ggplot2::geom_vline(
    xintercept = 1.0,
    linewidth = 0.55,
    linetype = "dashed",
    color = "grey45"
  ) +
  ggplot2::geom_hline(
    yintercept = 2.5,
    linewidth = 0.45,
    color = "grey60"
  ) +
  ggplot2::geom_col(
    position = bar_dodge,
    width = 0.64,
    color = "black",
    linewidth = 0.28
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(phylo_signal_plot_df, sig_label != ""),
    ggplot2::aes(
      x = estimate + 0.03,
      y = y_text,
      label = sig_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 3.4,
    color = "black"
  ) +
  ggplot2::scale_fill_manual(values = method_colors) +
  ggplot2::scale_alpha_manual(
    values = c(
      "No signal detected" = 0.28,
      "Signal detected" = 0.95
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, x_upper),
    expand = ggplot2::expansion(mult = c(0, 0.03))
  ) +
  ggplot2::labs(
    title = "Phylogenetic signal in key shape and joint-geometry traits",
    subtitle = "Color indicates metric",
    x = "Signal estimate",
    y = NULL,
    fill = NULL,
    alpha = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(color = "grey88", linewidth = 0.35),
    axis.text.y = ggplot2::element_text(size = 10, color = "black"),
    axis.text.x = ggplot2::element_text(size = 10, color = "black"),
    legend.position = "top",
    legend.box = "horizontal",
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10, color = "grey20")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1)),
    alpha = "none"
  )

save_plot_dual(
  p_phylo_signal_bar,
  "figure_phylogenetic_signal_grouped_clean_final",
  phylo_signal_plot_dir,
  width = 10.4,
  height = 6.8
)

write_clean_csv(
  phylo_signal_plot_df %>%
    dplyr::select(trait, trait_label, method, estimate, statistic, p_value, fdr_p_value, signal_status, sig_label),
  file.path(phylo_signal_plot_dir, "figure_phylogenetic_signal_grouped_clean_final_table.csv")
)

message("Clean grouped phylogenetic signal barplot saved to: ", phylo_signal_plot_dir)
print(p_phylo_signal_bar)


#####################################################################################

# ============================================================
# SECTION 22: MAIN-TEXT FIGURE - PGLS CORE ASSOCIATIONS
# PC1 ~ axial_span
# PC2 ~ axial_span
# ============================================================

message("SECTION 22: Building main-text PGLS figure for PC1/PC2 vs axial span ...")

pgls_plot_dir <- file.path(output_dir, "04_PGLS")
dir.create(pgls_plot_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1) Settings
# ------------------------------------------------------------
pgls_main_pairs <- tibble::tribble(
  ~response, ~predictor,    ~panel_title,
  "PC1",     "axial_span",  "PC1 vs axial span",
  "PC2",     "axial_span",  "PC2 vs axial span"
)

# optional: use Family colors if available, otherwise simple black points
use_family_colors <- TRUE

# ------------------------------------------------------------
# 2) Minimal checks
# ------------------------------------------------------------
needed_vars <- unique(c(pgls_main_pairs$response, pgls_main_pairs$predictor, "tree_label"))
missing_needed <- setdiff(needed_vars, names(tip_level_df))
if (length(missing_needed) > 0) {
  stop("Missing required variables in tip_level_df: ", paste(missing_needed, collapse = ", "))
}

if (!exists("pgls_results_df") || !is.data.frame(pgls_results_df) || nrow(pgls_results_df) == 0) {
  message("pgls_results_df not found or empty. Figure will be made without model annotation text.")
  pgls_results_df <- tibble::tibble()
}

# ------------------------------------------------------------
# 3) Family colors (only if available and desired)
# ------------------------------------------------------------
family_colors_local <- NULL
if (use_family_colors && "Family" %in% names(tip_level_df)) {
  fams_present <- tip_level_df %>%
    dplyr::filter(!is.na(Family)) %>%
    dplyr::pull(Family) %>%
    unique() %>%
    sort()
  
  if (length(fams_present) > 0) {
    pal <- if (length(fams_present) <= 12) {
      RColorBrewer::brewer.pal(max(3, min(12, length(fams_present))), "Set3")[seq_along(fams_present)]
    } else {
      scales::hue_pal()(length(fams_present))
    }
    family_colors_local <- setNames(pal, fams_present)
  }
}

# ------------------------------------------------------------
# 4) Helper to extract annotation from PGLS results
# ------------------------------------------------------------
extract_pgls_label <- function(res_df, response_name, predictor_name) {
  if (!is.data.frame(res_df) || nrow(res_df) == 0) return("PGLS summary unavailable")
  
  tmp <- res_df %>%
    dplyr::filter(
      response == response_name,
      predictor == predictor_name
    )
  
  if (nrow(tmp) == 0) return("PGLS summary unavailable")
  
  # keep non-intercept term if present
  if ("term" %in% names(tmp)) {
    tmp_nonint <- tmp %>% dplyr::filter(term != "(Intercept)")
    if (nrow(tmp_nonint) > 0) tmp <- tmp_nonint
  }
  
  tmp <- tmp[1, , drop = FALSE]
  
  est_txt <- if ("estimate" %in% names(tmp) && !is.na(tmp$estimate)) {
    paste0("slope = ", formatC(tmp$estimate, format = "f", digits = 3))
  } else {
    "slope = n/a"
  }
  
  p_txt <- if ("p_value" %in% names(tmp) && !is.na(tmp$p_value)) {
    if (tmp$p_value < 0.001) "p < 0.001" else paste0("p = ", formatC(tmp$p_value, format = "f", digits = 3))
  } else {
    "p = n/a"
  }
  
  fdr_txt <- if ("fdr_p_value" %in% names(tmp) && !is.na(tmp$fdr_p_value)) {
    if (tmp$fdr_p_value < 0.001) "FDR < 0.001" else paste0("FDR = ", formatC(tmp$fdr_p_value, format = "f", digits = 3))
  } else {
    "FDR = n/a"
  }
  
  lambda_txt <- if ("lambda" %in% names(tmp) && !is.na(tmp$lambda)) {
    paste0("lambda = ", formatC(tmp$lambda, format = "f", digits = 2))
  } else {
    "lambda = n/a"
  }
  
  paste(est_txt, p_txt, fdr_txt, lambda_txt, sep = "\n")
}

# ------------------------------------------------------------
# 5) Build individual panels
# ------------------------------------------------------------
pgls_main_plots <- vector("list", nrow(pgls_main_pairs))

for (i in seq_len(nrow(pgls_main_pairs))) {
  resp <- pgls_main_pairs$response[i]
  pred <- pgls_main_pairs$predictor[i]
  panel_title <- pgls_main_pairs$panel_title[i]
  
  plot_df <- tip_level_df %>%
    dplyr::select(dplyr::any_of(c("tree_label", "Family", resp, pred))) %>%
    dplyr::filter(!is.na(.data[[resp]]), !is.na(.data[[pred]]))
  
  ann_label <- extract_pgls_label(pgls_results_df, resp, pred)
  
  x_rng <- range(plot_df[[pred]], na.rm = TRUE)
  y_rng <- range(plot_df[[resp]], na.rm = TRUE)
  
  ann_x <- x_rng[1] + 0.03 * diff(x_rng)
  ann_y <- y_rng[2] - 0.04 * diff(y_rng)
  
  if (isTRUE(use_family_colors) && !is.null(family_colors_local) && "Family" %in% names(plot_df)) {
    p_tmp <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[pred]], y = .data[[resp]], color = Family)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::geom_smooth(
        method = "lm",
        se = TRUE,
        linewidth = 0.8,
        color = "black",
        fill = "grey75"
      ) +
      ggplot2::scale_color_manual(values = family_colors_local, na.translate = FALSE)
  } else {
    p_tmp <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[pred]], y = .data[[resp]])
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9, color = "black") +
      ggplot2::geom_smooth(
        method = "lm",
        se = TRUE,
        linewidth = 0.8,
        color = "black",
        fill = "grey75"
      )
  }
  
  p_tmp <- p_tmp +
    ggplot2::annotate(
      "label",
      x = ann_x,
      y = ann_y,
      label = ann_label,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      label.size = 0.25,
      fill = "white"
    ) +
    ggplot2::labs(
      title = panel_title,
      subtitle = "Trend line shown for visualisation; inference based on PGLS",
      x = "Axial span",
      y = resp,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      legend.position = if (i == 1 && isTRUE(use_family_colors) && !is.null(family_colors_local)) "right" else "none",
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 10)
    )
  
  pgls_main_plots[[i]] <- p_tmp
}

# ------------------------------------------------------------
# 6) Combine and save
# ------------------------------------------------------------
p_pgls_main <- pgls_main_plots[[1]] + pgls_main_plots[[2]] +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Phylogenetically informed associations between shape and axial span",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )
  )

save_plot_dual(
  p_pgls_main,
  "figure_pgls_core_pc1_pc2_vs_axial_span",
  pgls_plot_dir,
  width = 12,
  height = 6
)

# ------------------------------------------------------------
# 7) Export compact table used for panel annotation
# ------------------------------------------------------------
pgls_main_table <- pgls_results_df %>%
  dplyr::filter(
    response %in% c("PC1", "PC2"),
    predictor == "axial_span"
  ) %>%
  dplyr::filter(if ("term" %in% names(.)) term != "(Intercept)" else TRUE)

write_clean_csv(
  pgls_main_table,
  file.path(pgls_plot_dir, "figure_pgls_core_pc1_pc2_vs_axial_span_table.csv")
)

message("Main-text PGLS figure saved to: ", pgls_plot_dir)
print(p_pgls_main)


###############################################################################

# ============================================================
# SECTION 23: ROBUSTNESS CHECKS FOR KEY PGLS MODELS
# ============================================================

message("SECTION 23: Running robustness checks for key PGLS models ...")

robust_dir <- file.path(output_dir, "04_PGLS", "robustness_checks")
dir.create(robust_dir, showWarnings = FALSE, recursive = TRUE)

robust_pairs <- tibble::tribble(
  ~response, ~predictor,
  "PC1", "axial_span",
  "PC2", "axial_span"
)

# ------------------------------------------------------------
# Helper: extract compact result from safe_pgls output
# ------------------------------------------------------------
extract_compact_pgls <- function(pgls_obj, response, predictor, dropped_taxon = NA_character_, model_type = "original") {
  if (is.null(pgls_obj$result)) {
    return(tibble::tibble(
      response = response,
      predictor = predictor,
      model_type = model_type,
      dropped_taxon = dropped_taxon,
      term = NA_character_,
      estimate = NA_real_,
      p_value = NA_real_,
      lambda = NA_real_,
      n_taxa = NA_real_,
      reason = pgls_obj$reason
    ))
  }
  
  res <- pgls_obj$result
  res <- res %>%
    dplyr::rename_with(~ "estimate", dplyr::any_of(c("Estimate", "Value"))) %>%
    dplyr::rename_with(~ "std_error", dplyr::any_of(c("Std. Error", "Std.Error"))) %>%
    dplyr::rename_with(~ "statistic", dplyr::any_of(c("t value", "t-value"))) %>%
    dplyr::rename_with(~ "p_value", dplyr::any_of(c("Pr(>|t|)", "p-value")))
  
  if ("term" %in% names(res)) {
    res_nonint <- res %>% dplyr::filter(term != "(Intercept)")
    if (nrow(res_nonint) > 0) res <- res_nonint
  }
  
  res <- res[1, , drop = FALSE]
  
  tibble::tibble(
    response = response,
    predictor = predictor,
    model_type = model_type,
    dropped_taxon = dropped_taxon,
    term = if ("term" %in% names(res)) res$term else NA_character_,
    estimate = if ("estimate" %in% names(res)) res$estimate else NA_real_,
    p_value = if ("p_value" %in% names(res)) res$p_value else NA_real_,
    lambda = if ("lambda" %in% names(res)) res$lambda else NA_real_,
    n_taxa = if ("n_taxa" %in% names(res)) res$n_taxa else NA_real_,
    reason = NA_character_
  )
}

# ------------------------------------------------------------
# 1) Original models
# ------------------------------------------------------------
robust_original <- purrr::map_dfr(seq_len(nrow(robust_pairs)), function(i) {
  resp <- robust_pairs$response[i]
  pred <- robust_pairs$predictor[i]
  
  fit <- safe_pgls(tree_pruned, tip_level_df, resp, pred, predictor_is_factor = FALSE)
  extract_compact_pgls(fit, resp, pred, dropped_taxon = NA_character_, model_type = "original")
})

write_clean_csv(
  robust_original,
  file.path(robust_dir, "robustness_original_models.csv")
)

# ------------------------------------------------------------
# 2) Leave-one-out
# ------------------------------------------------------------
robust_loo <- purrr::map_dfr(seq_len(nrow(robust_pairs)), function(i) {
  resp <- robust_pairs$response[i]
  pred <- robust_pairs$predictor[i]
  
  dat_base <- tip_level_df %>%
    dplyr::select(tree_label, dplyr::all_of(resp), dplyr::all_of(pred)) %>%
    dplyr::filter(!is.na(.data[[resp]]), !is.na(.data[[pred]]))
  
  taxa <- dat_base$tree_label
  
  purrr::map_dfr(taxa, function(tx) {
    dat_minus_one <- tip_level_df %>%
      dplyr::filter(tree_label != tx)
    
    fit <- safe_pgls(tree_pruned, dat_minus_one, resp, pred, predictor_is_factor = FALSE)
    extract_compact_pgls(fit, resp, pred, dropped_taxon = tx, model_type = "leave_one_out")
  })
})

robust_loo <- robust_loo %>%
  dplyr::group_by(response, predictor) %>%
  dplyr::mutate(fdr_p_value = p_adjust_if_possible(p_value)) %>%
  dplyr::ungroup()

write_clean_csv(
  robust_loo,
  file.path(robust_dir, "robustness_leave_one_out.csv")
)

# ------------------------------------------------------------
# 3) Summary of leave-one-out stability
# ------------------------------------------------------------
robust_summary <- robust_loo %>%
  dplyr::group_by(response, predictor) %>%
  dplyr::summarise(
    n_models = dplyr::n(),
    n_failed = sum(is.na(estimate)),
    sign_consistent_positive = {
      est_ok <- estimate[!is.na(estimate)]
      length(est_ok) > 0 && all(est_ok > 0)
    },
    sign_consistent_negative = {
      est_ok <- estimate[!is.na(estimate)]
      length(est_ok) > 0 && all(est_ok < 0)
    },
    min_estimate = if (all(is.na(estimate))) NA_real_ else min(estimate, na.rm = TRUE),
    max_estimate = if (all(is.na(estimate))) NA_real_ else max(estimate, na.rm = TRUE),
    median_estimate = if (all(is.na(estimate))) NA_real_ else stats::median(estimate, na.rm = TRUE),
    min_p = if (all(is.na(p_value))) NA_real_ else min(p_value, na.rm = TRUE),
    max_p = if (all(is.na(p_value))) NA_real_ else max(p_value, na.rm = TRUE),
    n_nominal_sig = sum(p_value < 0.05, na.rm = TRUE),
    n_fdr_sig = sum(fdr_p_value < 0.05, na.rm = TRUE),
    min_lambda = if (all(is.na(lambda))) NA_real_ else min(lambda, na.rm = TRUE),
    max_lambda = if (all(is.na(lambda))) NA_real_ else max(lambda, na.rm = TRUE),
    median_lambda = if (all(is.na(lambda))) NA_real_ else stats::median(lambda, na.rm = TRUE),
    .groups = "drop"
  )

write_clean_csv(
  robust_summary,
  file.path(robust_dir, "robustness_leave_one_out_summary.csv")
)

# ------------------------------------------------------------
# 4) Simple diagnostic plots
# ------------------------------------------------------------
if (nrow(robust_loo) > 0) {
  robust_loo_est <- robust_loo %>% dplyr::filter(!is.na(estimate))
  robust_loo_lambda <- robust_loo %>% dplyr::filter(!is.na(lambda))
  
  if (nrow(robust_loo_est) > 0) {
    p_slope <- ggplot2::ggplot(
      robust_loo_est,
      ggplot2::aes(x = dropped_taxon, y = estimate)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_grid(response ~ predictor, scales = "free_y") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Leave-one-out robustness of PGLS slope estimates",
        x = "Dropped taxon",
        y = "Slope estimate"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 7)
      )

    save_plot_dual(
      p_slope,
      "robustness_leave_one_out_slopes",
      robust_dir,
      width = 10,
      height = 8
    )
  }
  
  if (nrow(robust_loo_lambda) > 0) {
    p_lambda <- ggplot2::ggplot(
      robust_loo_lambda,
      ggplot2::aes(x = dropped_taxon, y = lambda)
    ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_grid(response ~ predictor, scales = "free_y") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Leave-one-out robustness of estimated lambda",
        x = "Dropped taxon",
        y = "Estimated lambda"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 7)
      )

    save_plot_dual(
      p_lambda,
      "robustness_leave_one_out_lambda",
      robust_dir,
      width = 10,
      height = 8
    )
  }
}

message("SECTION 23 finished. Robustness outputs written to: ", robust_dir)


###############################################################################
###############################################################################
###############################################################################

# ============================================================
# SECTION 23B: TREE / TIME-CALIBRATION ROBUSTNESS ACROSS PCM
# ============================================================

message("SECTION 23B: Running tree-variant robustness checks across PCM blocks ...")

tree_robust_dir <- file.path(output_dir, "11_Tree_robustness")
dir.create(tree_robust_dir, showWarnings = FALSE, recursive = TRUE)

if (!exists("tree_variant_list") || length(tree_variant_list) == 0) {
  tree_variant_list <- list(working_tree = tree_pruned)
}

if (!exists("tree_variant_metadata") || !is.data.frame(tree_variant_metadata)) {
  tree_variant_metadata <- tibble::tibble(tree_id = names(tree_variant_list))
}

extract_root_fastAnc <- function(tree, trait_vec, trait_name) {
  tmp <- safe_fastAnc(tree, trait_vec, trait_name)
  if (is.null(tmp) || nrow(tmp) == 0) {
    return(tibble::tibble(
      trait = trait_name,
      node = NA_character_,
      estimate = NA_real_,
      variance = NA_real_,
      CI95_lower = NA_real_,
      CI95_upper = NA_real_,
      root_node_expected = as.character(ape::Ntip(tree) + 1),
      matched_root = FALSE
    ))
  }

  root_node_expected <- as.character(ape::Ntip(tree) + 1)
  root_row <- tmp %>%
    dplyr::filter(node == root_node_expected)

  if (nrow(root_row) == 0) root_row <- tmp[1, , drop = FALSE]

  root_row %>%
    dplyr::mutate(
      root_node_expected = root_node_expected,
      matched_root = node == root_node_expected
    )
}

robust_phylo_signal <- purrr::map_dfr(names(tree_variant_list), function(tree_id) {
  tr <- tree_variant_list[[tree_id]]
  purrr::map_dfr(qc_cont_vars, function(v) {
    vec <- tip_level_df[[v]]
    names(vec) <- tip_level_df$tree_label
    tmp <- safe_phylo_signal(tr, vec, v)
    if (is.null(tmp) || nrow(tmp) == 0) {
      return(tibble::tibble(
        tree_id = tree_id,
        trait = v,
        method = NA_character_,
        estimate = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        reason = "Signal estimation failed or too few taxa."
      ))
    }
    tmp %>%
      dplyr::mutate(tree_id = tree_id, reason = NA_character_) %>%
      dplyr::select(tree_id, dplyr::everything())
  })
}) %>%
  dplyr::left_join(tree_variant_metadata, by = "tree_id")

write_clean_csv(
  robust_phylo_signal,
  file.path(tree_robust_dir, "robustness_phylogenetic_signal_across_trees.csv")
)

robust_phylo_signal_summary <- robust_phylo_signal %>%
  dplyr::filter(!is.na(method)) %>%
  dplyr::group_by(trait, method) %>%
  dplyr::summarise(
    n_trees = dplyr::n(),
    estimate_min = min(estimate, na.rm = TRUE),
    estimate_max = max(estimate, na.rm = TRUE),
    estimate_median = stats::median(estimate, na.rm = TRUE),
    p_min = min(p_value, na.rm = TRUE),
    p_max = max(p_value, na.rm = TRUE),
    n_nominal_sig = sum(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write_clean_csv(
  robust_phylo_signal_summary,
  file.path(tree_robust_dir, "robustness_phylogenetic_signal_summary.csv")
)

robust_evol_models <- purrr::map_dfr(names(tree_variant_list), function(tree_id) {
  tr <- tree_variant_list[[tree_id]]
  purrr::map_dfr(qc_cont_vars, function(v) {
    vec <- tip_level_df[[v]]
    names(vec) <- tip_level_df$tree_label
    tmp <- safe_fit_continuous_models(tr, vec, v)
    if (is.null(tmp) || nrow(tmp) == 0) {
      return(tibble::tibble(
        tree_id = tree_id,
        trait = v,
        model = NA_character_,
        logLik = NA_real_,
        AIC = NA_real_,
        delta_AIC = NA_real_,
        best_model = NA_character_,
        reason = "Evolutionary model fit failed or too few taxa."
      ))
    }
    tmp %>%
      dplyr::mutate(tree_id = tree_id, reason = NA_character_) %>%
      dplyr::select(tree_id, dplyr::everything())
  })
}) %>%
  dplyr::left_join(tree_variant_metadata, by = "tree_id")

write_clean_csv(
  robust_evol_models,
  file.path(tree_robust_dir, "robustness_evolutionary_models_across_trees.csv")
)

robust_evol_summary <- robust_evol_models %>%
  dplyr::filter(!is.na(best_model)) %>%
  dplyr::count(trait, best_model, name = "n_trees") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(prop_trees = n_trees / sum(n_trees)) %>%
  dplyr::ungroup()

write_clean_csv(
  robust_evol_summary,
  file.path(tree_robust_dir, "robustness_evolutionary_models_best_model_frequency.csv")
)

robust_pgls_trees <- purrr::map_dfr(names(tree_variant_list), function(tree_id) {
  tr <- tree_variant_list[[tree_id]]
  purrr::map_dfr(seq_len(nrow(all_cont_models)), function(i) {
    resp <- all_cont_models$response[i]
    pred <- all_cont_models$predictor[i]
    fit <- safe_pgls(tr, tip_level_df, resp, pred, predictor_is_factor = FALSE)
    extract_compact_pgls(fit, resp, pred, dropped_taxon = NA_character_, model_type = "tree_variant") %>%
      dplyr::mutate(tree_id = tree_id, .before = 1)
  })
}) %>%
  dplyr::left_join(tree_variant_metadata, by = "tree_id")

robust_pgls_trees <- robust_pgls_trees %>%
  dplyr::group_by(response, predictor) %>%
  dplyr::mutate(fdr_p_value = p_adjust_if_possible(p_value)) %>%
  dplyr::ungroup()

write_clean_csv(
  robust_pgls_trees,
  file.path(tree_robust_dir, "robustness_pgls_across_trees.csv")
)

robust_pgls_summary <- robust_pgls_trees %>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::group_by(response, predictor) %>%
  dplyr::summarise(
    n_trees = dplyr::n(),
    estimate_min = min(estimate, na.rm = TRUE),
    estimate_max = max(estimate, na.rm = TRUE),
    estimate_median = stats::median(estimate, na.rm = TRUE),
    sign_consistent_positive = all(estimate > 0),
    sign_consistent_negative = all(estimate < 0),
    p_min = min(p_value, na.rm = TRUE),
    p_max = max(p_value, na.rm = TRUE),
    n_nominal_sig = sum(p_value < 0.05, na.rm = TRUE),
    n_fdr_sig = sum(fdr_p_value < 0.05, na.rm = TRUE),
    lambda_min = min(lambda, na.rm = TRUE),
    lambda_max = max(lambda, na.rm = TRUE),
    lambda_median = stats::median(lambda, na.rm = TRUE),
    .groups = "drop"
  )

write_clean_csv(
  robust_pgls_summary,
  file.path(tree_robust_dir, "robustness_pgls_summary.csv")
)

robust_asr_root <- purrr::map_dfr(names(tree_variant_list), function(tree_id) {
  tr <- tree_variant_list[[tree_id]]
  purrr::map_dfr(asr_cont_vars, function(v) {
    vec <- tip_level_df[[v]]
    names(vec) <- tip_level_df$tree_label
    extract_root_fastAnc(tr, vec, v) %>%
      dplyr::mutate(tree_id = tree_id, .before = 1)
  })
}) %>%
  dplyr::left_join(tree_variant_metadata, by = "tree_id")

write_clean_csv(
  robust_asr_root,
  file.path(tree_robust_dir, "robustness_asr_root_estimates_across_trees.csv")
)

robust_asr_root_summary <- robust_asr_root %>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::group_by(trait) %>%
  dplyr::summarise(
    n_trees = dplyr::n(),
    estimate_min = min(estimate, na.rm = TRUE),
    estimate_max = max(estimate, na.rm = TRUE),
    estimate_median = stats::median(estimate, na.rm = TRUE),
    ci_lower_min = min(CI95_lower, na.rm = TRUE),
    ci_upper_max = max(CI95_upper, na.rm = TRUE),
    .groups = "drop"
  )

write_clean_csv(
  robust_asr_root_summary,
  file.path(tree_robust_dir, "robustness_asr_root_estimates_summary.csv")
)

if (nrow(robust_pgls_trees) > 0) {
  robust_pgls_plot_df <- robust_pgls_trees %>% dplyr::filter(!is.na(estimate), !is.na(tree_id))

  if (nrow(robust_pgls_plot_df) > 0) {
  p_tree_robust <- ggplot2::ggplot(
    robust_pgls_plot_df,
    ggplot2::aes(x = reorder(tree_id, estimate), y = estimate, color = variant_group)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::facet_grid(response ~ predictor, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "PGLS sensitivity across available tree variants",
      x = "Tree variant",
      y = "Slope estimate",
      color = "Tree class"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7)
    )

  save_plot_dual(
    p_tree_robust,
    "robustness_pgls_across_trees_slopes",
    tree_robust_dir,
    width = 12,
    height = 10
  )
  }
}

message("SECTION 23B finished. Tree robustness outputs written to: ", tree_robust_dir)


# ============================================================
# SECTION 24: CONTINUOUS ASR / CONTMAP FOR SHAPE + SCREW-JOINT GEOMETRY
# Combined figure with shared legend
# BLACK BACKGROUND VERSION - REFINED
# ============================================================

message("SECTION 24: Building combined contMap figure for PC1-PC5 plus screw-joint geometry ...")

asr_plot_dir <- file.path(output_dir, "08_ASR")
dir.create(asr_plot_dir, showWarnings = FALSE, recursive = TRUE)

asr_main_traits <- paste0("PC", 1:5)
asr_geom_traits <- c("abs_winding_angle_deg", "axial_span")
asr_figure_traits <- c(asr_main_traits, asr_geom_traits)
asr_figure_traits <- asr_figure_traits[asr_figure_traits %in% names(tip_level_df)]

if (length(asr_figure_traits) == 0) {
  stop("None of the requested ASR figure traits are available in tip_level_df.")
}

build_contmap_obj <- function(tree, df, trait_name) {
  dat <- df %>%
    dplyr::select(tree_label, dplyr::all_of(trait_name)) %>%
    dplyr::filter(!is.na(.data[[trait_name]]), !is.na(tree_label))

  if (nrow(dat) < 5) return(NULL)

  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
  vec <- dat[[trait_name]]
  names(vec) <- dat$tree_label
  vec <- vec[tr$tip.label]

  tryCatch(
    suppressWarnings(
      phytools::contMap(
        tree = tr,
        x = vec,
        plot = FALSE
      )
    ),
    error = function(e) {
      message("contMap failed for ", trait_name, ": ", conditionMessage(e))
      NULL
    }
  )
}

draw_shared_contmap_legend_black <- function(cols, zlim, title = "Standardized trait score") {
  plot.new()
  par(usr = c(0, 1, 0, 1))

  rect(0, 0, 1, 1, col = "black", border = NA)
  text(0.5, 0.80, labels = title, cex = 1.55, font = 2, col = "white")

  x_left <- 0.12
  x_right <- 0.88
  y_bottom <- 0.40
  y_top <- 0.58

  n_cols <- length(cols)
  xs <- seq(x_left, x_right, length.out = n_cols + 1)

  for (i in seq_len(n_cols)) {
    rect(xs[i], y_bottom, xs[i + 1], y_top, col = cols[i], border = cols[i])
  }

  rect(x_left, y_bottom, x_right, y_top, border = "grey85", lwd = 1.0)

  tick_vals <- pretty(zlim, n = 5)
  tick_vals <- tick_vals[tick_vals >= zlim[1] & tick_vals <= zlim[2]]

  if (length(tick_vals) >= 2 && diff(zlim) > 0) {
    tick_pos <- x_left + (tick_vals - zlim[1]) / diff(zlim) * (x_right - x_left)
    segments(tick_pos, y_bottom - 0.04, tick_pos, y_bottom, col = "grey85", lwd = 1.0)
    text(
      tick_pos,
      y_bottom - 0.10,
      labels = formatC(tick_vals, format = "f", digits = 1),
      cex = 1.25,
      col = "grey95"
    )
  }

  text(x_left, 0.20, labels = "Low", adj = c(0, 0.5), cex = 1.30, col = "white", font = 2)
  text(x_right, 0.20, labels = "High", adj = c(1, 0.5), cex = 1.30, col = "white", font = 2)
}

draw_empty_black_panel <- function() {
  plot.new()
  par(usr = c(0, 1, 0, 1))
  rect(0, 0, 1, 1, col = "black", border = NA)
}

trait_display_labels <- c(
  setNames(asr_main_traits, asr_main_traits),
  "abs_winding_angle_deg" = "Winding angle",
  "axial_span" = "Axial span"
)

asr_df_scaled <- tip_level_df %>%
  dplyr::select(tree_label, dplyr::all_of(asr_figure_traits))

for (trait_name in asr_figure_traits) {
  vals <- asr_df_scaled[[trait_name]]
  if (sum(is.finite(vals)) >= 2 && stats::sd(vals, na.rm = TRUE) > 0) {
    asr_df_scaled[[trait_name]] <- as.numeric(scale(vals))
  } else {
    asr_df_scaled[[trait_name]] <- NA_real_
  }
}

contmap_list <- purrr::map(asr_figure_traits, ~ build_contmap_obj(tree_pruned, asr_df_scaled, .x))
names(contmap_list) <- asr_figure_traits

failed_traits <- names(contmap_list)[vapply(contmap_list, is.null, logical(1))]
if (length(failed_traits) > 0) {
  warning("Could not build contMap for: ", paste(failed_traits, collapse = ", "))
}

contmap_list <- contmap_list[!vapply(contmap_list, is.null, logical(1))]

if (length(contmap_list) == 0) {
  stop("No contMap objects could be built.")
}

all_z_vals <- unlist(
  lapply(names(contmap_list), function(trait_name) asr_df_scaled[[trait_name]]),
  use.names = FALSE
)
all_z_vals <- all_z_vals[is.finite(all_z_vals)]

if (length(all_z_vals) == 0) {
  stop("No finite standardized values available for shared legend.")
}

global_zlim <- range(all_z_vals, na.rm = TRUE)

shared_palette_fun <- colorRampPalette(
  c("#46B1E8", "#7FD3E6", "#D9F0E6", "#F2E85C", "#F7B267", "#F4845F", "#D7263D")
)

contmap_list <- lapply(contmap_list, function(cm) {
  new_cols <- shared_palette_fun(length(cm$cols))
  names(new_cols) <- names(cm$cols)
  cm$cols <- new_cols
  cm
})

shared_cols <- shared_palette_fun(260)
plot_order <- names(contmap_list)

n_cols_layout <- 3
n_panels_total <- length(plot_order) + 1
n_rows_layout <- ceiling(n_panels_total / n_cols_layout)
blank_panel_id <- length(plot_order) + 2
panel_ids <- c(seq_along(plot_order), length(plot_order) + 1)
if (length(panel_ids) < n_rows_layout * n_cols_layout) {
  panel_ids <- c(panel_ids, rep(blank_panel_id, n_rows_layout * n_cols_layout - length(panel_ids)))
}
layout_mat <- matrix(panel_ids, nrow = n_rows_layout, byrow = TRUE)

combined_pdf <- file.path(asr_plot_dir, "figure_contASR_shape_geometry_combined_black_refined.pdf")
combined_jpg <- file.path(asr_plot_dir, "figure_contASR_shape_geometry_combined_black_refined.jpg")

draw_combined_contmap_panel_set <- function() {
  layout(layout_mat)
  par(
    oma = c(0, 0, 2.6, 0),
    xpd = NA,
    bg = "black",
    fg = "white",
    col = "white",
    col.main = "white",
    col.lab = "white",
    col.axis = "white"
  )

  for (trait_name in plot_order) {
    par(
      mar = c(1.5, 1.3, 1.3, 1.3),
      bg = "black",
      fg = "white",
      col = "white",
      col.main = "white",
      col.lab = "white",
      col.axis = "white"
    )
    plot(
      contmap_list[[trait_name]],
      type = "phylogram",
      ftype = "off",
      lwd = 4.2,
      legend = 0,
      fsize = 0.72,
      outline = FALSE
    )
    graphics::mtext(
      text = dplyr::coalesce(unname(trait_display_labels[[trait_name]]), trait_name),
      side = 1,
      line = -1.2,
      adj = 0.02,
      cex = 1.2,
      font = 2,
      col = "white"
    )
  }

  par(
    mar = c(1.8, 1.6, 1.8, 1.6),
    bg = "black",
    fg = "white",
    col = "white",
    col.main = "white",
    col.lab = "white",
    col.axis = "white"
  )
  draw_shared_contmap_legend_black(
    cols = shared_cols,
    zlim = global_zlim,
    title = "Standardized trait score"
  )

  if (any(layout_mat == blank_panel_id)) {
    n_blank <- sum(layout_mat[, ] == blank_panel_id)
    for (i in seq_len(n_blank)) {
      par(
        mar = c(0, 0, 0, 0),
        bg = "black",
        fg = "white",
        col = "white"
      )
      draw_empty_black_panel()
    }
  }

  mtext(
    "Continuous ancestral reconstructions of shape and screw-joint geometry",
    outer = TRUE,
    cex = 1.12,
    font = 2,
    line = 0.9,
    col = "white"
  )
}

grDevices::pdf(combined_pdf, width = 16, height = 13, bg = "black")
op <- par(no.readonly = TRUE)
draw_combined_contmap_panel_set()
par(op)
grDevices::dev.off()

grDevices::jpeg(combined_jpg, width = 4800, height = 3900, res = 300, bg = "black")
op <- par(no.readonly = TRUE)
draw_combined_contmap_panel_set()
par(op)
grDevices::dev.off()

for (trait_name in names(contmap_list)) {
  pdf_ind <- file.path(asr_plot_dir, paste0("figure_contASR_", trait_name, "_single_black_refined.pdf"))
  jpg_ind <- file.path(asr_plot_dir, paste0("figure_contASR_", trait_name, "_single_black_refined.jpg"))

  grDevices::pdf(pdf_ind, width = 8.2, height = 8.2, bg = "black")
  try({
    par(
      mar = c(2.0, 2.0, 3.0, 2.0),
      bg = "black",
      fg = "white",
      col = "white",
      col.main = "white",
      col.lab = "white",
      col.axis = "white"
    )
    plot(
      contmap_list[[trait_name]],
      type = "phylogram",
      ftype = "off",
      lwd = 4.2,
      legend = 0.9,
      fsize = 0.75,
      outline = FALSE
    )
    title(
      main = dplyr::coalesce(unname(trait_display_labels[[trait_name]]), trait_name),
      line = 1,
      cex.main = 1.12,
      font.main = 2,
      col.main = "white"
    )
  }, silent = TRUE)
  grDevices::dev.off()

  grDevices::jpeg(jpg_ind, width = 2400, height = 2400, res = 300, bg = "black")
  try({
    par(
      mar = c(2.0, 2.0, 3.0, 2.0),
      bg = "black",
      fg = "white",
      col = "white",
      col.main = "white",
      col.lab = "white",
      col.axis = "white"
    )
    plot(
      contmap_list[[trait_name]],
      type = "phylogram",
      ftype = "off",
      lwd = 4.2,
      legend = 0.9,
      fsize = 0.75,
      outline = FALSE
    )
    title(
      main = dplyr::coalesce(unname(trait_display_labels[[trait_name]]), trait_name),
      line = 1,
      cex.main = 1.12,
      font.main = 2,
      col.main = "white"
    )
  }, silent = TRUE)
  grDevices::dev.off()
}

write_clean_csv(
  asr_df_scaled,
  file.path(asr_plot_dir, "figure_contASR_shape_geometry_standardized_input.csv")
)

if (exists("asr_cont_results") && is.data.frame(asr_cont_results) && nrow(asr_cont_results) > 0) {
  asr_main_table <- asr_cont_results %>%
    dplyr::filter(trait %in% names(contmap_list))

  write_clean_csv(
    asr_main_table,
    file.path(asr_plot_dir, "figure_contASR_shape_geometry_table.csv")
  )
}

message("Combined contASR figure for shape + screw-joint geometry saved to: ", asr_plot_dir)


################################################################################
################################################################################
################################################################################

# ============================================================
# SECTION 26: MAIN-TEXT FIGURE - UNIVARIATE + MULTIVARIATE
# EVOLUTIONARY MODEL SUPPORT (CLEAN DUMBBELL VERSION)
# ============================================================

message("SECTION 26: Building combined univariate + multivariate evolutionary model figure (clean dumbbell version) ...")

evol_plot_dir <- file.path(output_dir, "03_Evolutionary_models")
dir.create(evol_plot_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1) Checks
# ------------------------------------------------------------
if (!exists("evol_model_results") || !is.data.frame(evol_model_results) || nrow(evol_model_results) == 0) {
  stop("evol_model_results not found or empty.")
}

if (!exists("mv_model_results") || !is.data.frame(mv_model_results) || nrow(mv_model_results) == 0) {
  stop("mv_model_results not found or empty.")
}

# ------------------------------------------------------------
# 2) Prepare univariate data
# ------------------------------------------------------------
uni_df <- evol_model_results %>%
  dplyr::filter(trait %in% paste0("PC", 1:5)) %>%
  dplyr::mutate(
    analysis_level = "Univariate",
    set_label = as.character(trait),
    model = factor(model, levels = c("BM", "OU", "EB")),
    converged = as.logical(converged)
  )

uni_conv_df <- uni_df %>%
  dplyr::filter(converged, !is.na(AICc)) %>%
  dplyr::group_by(set_label) %>%
  dplyr::mutate(
    delta_plot = AICc - min(AICc, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

uni_plot_df <- uni_df %>%
  dplyr::left_join(
    uni_conv_df %>% dplyr::select(set_label, model, delta_plot),
    by = c("set_label", "model")
  ) %>%
  dplyr::mutate(
    fit_status = dplyr::case_when(
      converged & !is.na(AICc) ~ "converged",
      TRUE ~ "failed"
    ),
    set_label = factor(set_label, levels = rev(paste0("PC", 1:5)))
  )

# ------------------------------------------------------------
# 3) Prepare multivariate data
# ------------------------------------------------------------
pc_set_labels <- c(
  "PC1_5" = "PC1-PC5",
  "PC1_4" = "PC1-PC4"
)

mv_df <- mv_model_results %>%
  dplyr::mutate(
    analysis_level = "Multivariate",
    set_label = dplyr::if_else(
      as.character(pc_set) %in% names(pc_set_labels),
      pc_set_labels[as.character(pc_set)],
      as.character(pc_set)
    ),
    model = factor(model, levels = c("BM", "OU", "EB")),
    converged = as.logical(converged)
  )

mv_conv_df <- mv_df %>%
  dplyr::filter(converged, !is.na(AIC)) %>%
  dplyr::group_by(set_label) %>%
  dplyr::mutate(
    delta_plot = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

mv_plot_df <- mv_df %>%
  dplyr::left_join(
    mv_conv_df %>% dplyr::select(set_label, model, delta_plot),
    by = c("set_label", "model")
  ) %>%
  dplyr::mutate(
    fit_status = dplyr::case_when(
      converged & !is.na(AIC) ~ "converged",
      TRUE ~ "failed"
    ),
    set_label = factor(set_label, levels = rev(c("PC1-PC5", "PC1-PC4")))
  )

# ------------------------------------------------------------
# 4) Common styling data
# ------------------------------------------------------------
model_cols <- c(
  "BM" = "#1b9e77",
  "OU" = "#7570b3",
  "EB" = "#d95f02"
)

model_shapes <- c(
  "BM" = 16,
  "OU" = 17,
  "EB" = 15
)

failed_col <- "#d7301f"

uni_segment_df <- uni_plot_df %>%
  dplyr::filter(fit_status == "converged", !is.na(delta_plot)) %>%
  dplyr::group_by(set_label) %>%
  dplyr::summarise(
    x_min = min(delta_plot, na.rm = TRUE),
    x_max = max(delta_plot, na.rm = TRUE),
    .groups = "drop"
  )

mv_segment_df <- mv_plot_df %>%
  dplyr::filter(fit_status == "converged", !is.na(delta_plot)) %>%
  dplyr::group_by(set_label) %>%
  dplyr::summarise(
    x_min = min(delta_plot, na.rm = TRUE),
    x_max = max(delta_plot, na.rm = TRUE),
    .groups = "drop"
  )

failed_x <- 11.6

# ------------------------------------------------------------
# 5) Panel A: Univariate
# ------------------------------------------------------------
p_uni <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, color = "grey35", linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = 2, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  ggplot2::geom_segment(
    data = uni_segment_df,
    ggplot2::aes(x = x_min, xend = x_max, y = set_label, yend = set_label),
    color = "grey82",
    linewidth = 1.0
  ) +
  ggplot2::geom_point(
    data = uni_plot_df %>% dplyr::filter(fit_status == "converged", !is.na(delta_plot)),
    ggplot2::aes(x = delta_plot, y = set_label, color = model, shape = model),
    size = 4.2,
    stroke = 0.8
  ) +
  ggplot2::scale_color_manual(values = model_cols) +
  ggplot2::scale_shape_manual(values = model_shapes) +
  ggplot2::scale_x_continuous(
    limits = c(-0.2, 7.2),
    breaks = c(0, 2, 5),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Univariate model support",
    x = expression(Delta*AIC[c]),
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold", size = 13)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(ncol = 1, byrow = TRUE),
    shape = ggplot2::guide_legend(ncol = 1, byrow = TRUE)
  )

# ------------------------------------------------------------
# 6) Panel B: Multivariate
# ------------------------------------------------------------
p_mv <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, color = "grey35", linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = 2, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  ggplot2::geom_vline(xintercept = 10, color = "grey75", linewidth = 0.4, linetype = "dashed") +
  ggplot2::geom_segment(
    data = mv_segment_df,
    ggplot2::aes(x = x_min, xend = x_max, y = set_label, yend = set_label),
    color = "grey82",
    linewidth = 1.0
  ) +
  ggplot2::geom_point(
    data = mv_plot_df %>% dplyr::filter(fit_status == "converged", !is.na(delta_plot)),
    ggplot2::aes(x = delta_plot, y = set_label, color = model, shape = model),
    size = 4.2,
    stroke = 0.8
  ) +
  ggplot2::geom_point(
    data = mv_plot_df %>% dplyr::filter(fit_status == "failed"),
    ggplot2::aes(x = failed_x, y = set_label),
    shape = 4,
    size = 4.8,
    stroke = 1.2,
    color = failed_col
  ) +
  ggplot2::geom_text(
    data = mv_plot_df %>% dplyr::filter(fit_status == "failed"),
    ggplot2::aes(x = failed_x + 0.35, y = set_label, label = "failed"),
    hjust = 0,
    vjust = 0.5,
    size = 3.4,
    color = failed_col
  ) +
  ggplot2::scale_color_manual(values = model_cols) +
  ggplot2::scale_shape_manual(values = model_shapes) +
  ggplot2::scale_x_continuous(
    limits = c(-0.2, 14),
    breaks = c(0, 2, 5, 10),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Multivariate model support",
    x = expression(Delta*AIC),
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", size = 13)
  )

# ------------------------------------------------------------
# 7) Combined figure
# ------------------------------------------------------------
p_evol_combined <- p_uni / p_mv +
  patchwork::plot_layout(guides = "collect", heights = c(1.35, 0.9)) +
  patchwork::plot_annotation(
    title = "Evolutionary model support across univariate and multivariate shape representations",
    subtitle = "Lower values indicate stronger support; dashed lines mark AIC = 2 and 10 where shown",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
  )

save_plot_dual(
  p_evol_combined,
  "figure_univariate_multivariate_evolutionary_models_dumbbell",
  evol_plot_dir,
  width = 10.6,
  height = 8.2
)

# ------------------------------------------------------------
# 8) Export figure tables
# ------------------------------------------------------------
uni_fig_table <- uni_plot_df %>%
  dplyr::select(analysis_level, set_label, model, fit_status, delta_plot, converged, AIC, AICc, logLik)

mv_fig_table <- mv_plot_df %>%
  dplyr::select(analysis_level, set_label, model, fit_status, delta_plot, converged, AIC, logLik)

write_clean_csv(
  uni_fig_table,
  file.path(evol_plot_dir, "figure_univariate_evolutionary_models_dumbbell_table.csv")
)

write_clean_csv(
  mv_fig_table,
  file.path(evol_plot_dir, "figure_multivariate_evolutionary_models_dumbbell_table.csv")
)

message("Combined univariate + multivariate evolutionary model dumbbell figure saved to: ", evol_plot_dir)
print(p_evol_combined)
