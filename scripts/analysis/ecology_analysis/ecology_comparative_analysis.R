# Broad ecological association tests used in the manuscript.
#
# Set project_root to the local directory containing the analysis_data input
# and results folders before running this script.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ape)
  library(caper)
  library(phytools)
})

# ------------------- PATHS -------------------
project_root <- "<MANUSCRIPT_PROJECT_ROOT>"
ecology_path <- file.path(project_root, "analysis_data", "Input", "ecology_tree_tip_matrix.csv")
pcm_tip_path <- file.path(project_root, "analysis_data", "Results", "PCM", "tip_level_dataset_used_for_PCM.csv")
tree_path <- file.path(project_root, "analysis_data", "Input", "curculionoidea_primary_tree.tre")
out_dir <- file.path(project_root, "analysis_data", "Results", "Ecology")
plot_dir <- file.path(out_dir, "plots")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ------------------- OUTPUTS -------------------
out_input_merged <- file.path(out_dir, "ecology_analysis_input_merged.csv")
out_factor_counts <- file.path(out_dir, "ecology_broad_factor_counts.csv")
out_group_summary <- file.path(out_dir, "ecology_trait_group_summary.csv")
out_nonphylo <- file.path(out_dir, "ecology_nonphylo_group_tests.csv")
out_phylanova <- file.path(out_dir, "ecology_phylogenetic_anova_results.csv")
out_pgls <- file.path(out_dir, "ecology_pgls_factor_results.csv")
out_manifest <- file.path(out_dir, "ecology_analysis_manifest.txt")

# ------------------- SETTINGS -------------------
response_traits <- c("PC1", "PC2", "abs_winding_angle_deg", "axial_span", "centroid_size")
min_total_tips <- 5
min_group_n <- 2
phylanova_nsim <- 1000
set.seed(20260627)

# ------------------- HELPERS -------------------
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_factor <- function(x, levels_keep = NULL) {
  x <- as.character(x)
  x[trimws(x) == ""] <- NA_character_
  if (!is.null(levels_keep)) {
    x[!(x %in% levels_keep)] <- NA_character_
    return(factor(x, levels = levels_keep))
  }
  factor(x)
}

fmt_msg <- function(...) cat(paste0(..., "\n"))

write_output <- function(x, path) {
  write.csv2(x, path, row.names = FALSE)
  write.csv(x, sub("\\.csv$", "_gsheets.csv", path), row.names = FALSE)
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
          caper::pgls(formula = formula_obj, data = comp_obj, lambda = att$lambda, bounds = att$bounds)
        } else {
          caper::pgls(formula = formula_obj, data = comp_obj, lambda = att$lambda)
        }
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(fit_try)) {
      return(list(model = fit_try, method = att$label, error = NA_character_))
    }
  }
  list(model = NULL, method = NA_character_, error = last_error)
}

extract_pgls_factor_rows <- function(model, response, predictor, fit_method) {
  sm <- summary(model)
  cf <- as.data.frame(sm$coefficients)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  cf <- cf[cf$term != "(Intercept)", , drop = FALSE]
  if (nrow(cf) == 0) {
    return(data.frame())
  }
  data.frame(
    response = response,
    predictor = predictor,
    term = cf$term,
    estimate = cf$Estimate,
    std_error = cf$`Std. Error`,
    t_value = cf$`t value`,
    p_value = cf$`Pr(>|t|)`,
    lambda = if (!is.null(model$param["lambda"])) unname(model$param["lambda"]) else NA_real_,
    logLik = as.numeric(logLik(model)),
    AIC = AIC(model),
    n = nobs(model),
    model_fit_strategy = fit_method,
    stringsAsFactors = FALSE
  )
}

sanitize_tree_labels <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}

# ------------------- READ INPUT -------------------
ecology <- read.csv2(ecology_path, stringsAsFactors = FALSE, check.names = FALSE)
pcm <- read.csv2(pcm_tip_path, stringsAsFactors = FALSE, check.names = FALSE)
tree <- ape::read.tree(tree_path)
tree$tip.label <- sanitize_tree_labels(tree$tip.label)
tree$node.label <- NULL

for (v in response_traits) {
  pcm[[v]] <- safe_num(pcm[[v]])
}

# ------------------- BROAD ECOLOGY FACTORS -------------------
ecology2 <- ecology %>%
  dplyr::mutate(
    # Mixed composite tips are deliberately unresolved in broad tests.
    host_lineage_broad = dplyr::case_when(
      host_lineage_simple == "angiosperm" ~ "angiosperm",
      host_lineage_simple == "gymnosperm" ~ "gymnosperm",
      TRUE ~ NA_character_
    ),
    woody_association_broad = dplyr::case_when(
      woody_association_simple %in% c("woody", "woody_or_dead_plant_associated") ~ "woody",
      woody_association_simple %in% c("mostly_nonwoody", "nonwoody_or_shrub_associated") ~ "nonwoody",
      TRUE ~ NA_character_
    ),
    larval_lifestyle_broad = dplyr::case_when(
      larval_lifestyle_simple == "internal" ~ "internal",
      larval_lifestyle_simple %in% c("mixed", "external_root_feeding") ~ "other",
      TRUE ~ NA_character_
    ),
    fungal_association_broad = dplyr::case_when(
      fungal_symbiosis_simple %in% c("yes", "yes_or_common") ~ "yes",
      fungal_symbiosis_simple == "no" ~ "no",
      TRUE ~ NA_character_
    )
  )

broad_predictors <- c(
  "host_lineage_broad",
  "woody_association_broad",
  "larval_lifestyle_broad",
  "fungal_association_broad"
)

# ------------------- MERGE -------------------
analysis_df <- pcm %>%
  dplyr::select(tree_label, dplyr::all_of(response_traits), Family, tree_tip, n_specimens) %>%
  dplyr::left_join(
    ecology2 %>%
      dplyr::select(
        tree_label, representative_specimens, coding_confidence, ecology_notes,
        trophic_mode_simple, host_lineage_simple, substrate_simple, larval_lifestyle_simple,
        woody_association_simple, fungal_symbiosis_simple, aquatic_association_simple,
        dplyr::all_of(broad_predictors)
      ),
    by = "tree_label"
  )

write_output(analysis_df, out_input_merged)

# ------------------- FACTOR COUNTS -------------------
factor_counts <- dplyr::bind_rows(lapply(broad_predictors, function(pred) {
  analysis_df %>%
    dplyr::count(.data[[pred]], name = "n_tips") %>%
    dplyr::mutate(predictor = pred, level = .data[[pred]]) %>%
    dplyr::select(predictor, level, n_tips)
}))
write_output(factor_counts, out_factor_counts)

# ------------------- DESCRIPTIVE GROUP SUMMARIES -------------------
group_summary <- dplyr::bind_rows(lapply(broad_predictors, function(pred) {
  dplyr::bind_rows(lapply(response_traits, function(resp) {
    dat <- analysis_df %>%
      dplyr::select(tree_label, dplyr::all_of(pred), dplyr::all_of(resp)) %>%
      dplyr::filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
    if (nrow(dat) == 0) return(data.frame())
    dat %>%
      dplyr::group_by(.data[[pred]]) %>%
      dplyr::summarise(
        n_tips = n(),
        mean = mean(.data[[resp]], na.rm = TRUE),
        sd = sd(.data[[resp]], na.rm = TRUE),
        median = median(.data[[resp]], na.rm = TRUE),
        min = min(.data[[resp]], na.rm = TRUE),
        max = max(.data[[resp]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        predictor = pred,
        response = resp,
        level = .data[[pred]]
      ) %>%
      dplyr::select(predictor, response, level, n_tips, mean, sd, median, min, max)
  }))
}))
write_output(group_summary, out_group_summary)

# ------------------- NON-PHYLO GROUP TESTS -------------------
nonphylo_results <- dplyr::bind_rows(lapply(broad_predictors, function(pred) {
  dplyr::bind_rows(lapply(response_traits, function(resp) {
    dat <- analysis_df %>%
      dplyr::select(tree_label, dplyr::all_of(pred), dplyr::all_of(resp)) %>%
      dplyr::filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
    if (nrow(dat) < min_total_tips) return(data.frame())
    counts <- table(dat[[pred]])
    if (length(counts) != 2 || any(counts < min_group_n)) return(data.frame())
    wt <- tryCatch(wilcox.test(dat[[resp]] ~ dat[[pred]], exact = FALSE), error = function(e) NULL)
    if (is.null(wt)) return(data.frame())
    data.frame(
      predictor = pred,
      response = resp,
      n = nrow(dat),
      group_levels = paste(names(counts), collapse = " vs "),
      group_n = paste(as.integer(counts), collapse = " vs "),
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }))
}))
if (nrow(nonphylo_results) > 0) {
  nonphylo_results$p_adj_fdr <- p.adjust(nonphylo_results$p_value, method = "fdr")
}
write_output(nonphylo_results, out_nonphylo)

# ------------------- PHYLOGENETIC ANOVA -------------------
phylanova_results <- dplyr::bind_rows(lapply(broad_predictors, function(pred) {
  dplyr::bind_rows(lapply(response_traits, function(resp) {
    dat <- analysis_df %>%
      dplyr::select(tree_label, dplyr::all_of(pred), dplyr::all_of(resp)) %>%
      dplyr::filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
    if (nrow(dat) < min_total_tips) return(data.frame())
    counts <- table(dat[[pred]])
    if (length(counts) != 2 || any(counts < min_group_n)) return(data.frame())
    tr <- drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
    dat <- dat[match(tr$tip.label, dat$tree_label), , drop = FALSE]
    grp <- as.character(dat[[pred]])
    names(grp) <- dat$tree_label
    y <- dat[[resp]]
    names(y) <- dat$tree_label
    pa <- tryCatch(
      phytools::phylANOVA(tree = tr, x = grp, y = y, nsim = phylanova_nsim),
      error = function(e) NULL
    )
    if (is.null(pa)) return(data.frame())
    data.frame(
      predictor = pred,
      response = resp,
      n = nrow(dat),
      group_levels = paste(names(counts), collapse = " vs "),
      group_n = paste(as.integer(counts), collapse = " vs "),
      f_statistic = if (!is.null(pa$F)) unname(pa$F[1]) else NA_real_,
      p_value = if (!is.null(pa$Pf)) unname(pa$Pf[1]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))
if (nrow(phylanova_results) > 0) {
  phylanova_results$p_adj_fdr <- p.adjust(phylanova_results$p_value, method = "fdr")
}
write_output(phylanova_results, out_phylanova)

# ------------------- PGLS WITH ECOLOGY FACTORS -------------------
pgls_results <- dplyr::bind_rows(lapply(broad_predictors, function(pred) {
  dplyr::bind_rows(lapply(response_traits, function(resp) {
    dat <- analysis_df %>%
      dplyr::select(tree_label, dplyr::all_of(pred), dplyr::all_of(resp)) %>%
      dplyr::filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
    if (nrow(dat) < min_total_tips) return(data.frame())
    counts <- table(dat[[pred]])
    if (length(counts) != 2 || any(counts < min_group_n)) return(data.frame())
    tr <- drop.tip(tree, setdiff(tree$tip.label, dat$tree_label))
    dat <- dat[match(tr$tip.label, dat$tree_label), , drop = FALSE]
    dat[[pred]] <- factor(dat[[pred]])
    comp <- tryCatch(
      caper::comparative.data(
        phy = tr,
        data = as.data.frame(dat),
        names.col = "tree_label",
        vcv = TRUE,
        warn.dropped = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(comp)) return(data.frame())
    form <- stats::as.formula(paste(resp, "~", pred))
    fit_info <- fit_pgls_with_fallback(formula_obj = form, comp_obj = comp)
    if (is.null(fit_info$model)) {
      return(data.frame(
        response = resp,
        predictor = pred,
        term = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        t_value = NA_real_,
        p_value = NA_real_,
        lambda = NA_real_,
        logLik = NA_real_,
        AIC = NA_real_,
        n = nrow(dat),
        model_fit_strategy = NA_character_,
        error = fit_info$error,
        stringsAsFactors = FALSE
      ))
    }
    out <- extract_pgls_factor_rows(fit_info$model, response = resp, predictor = pred, fit_method = fit_info$method)
    out$error <- NA_character_
    out
  }))
}))
if (nrow(pgls_results) > 0) {
  pgls_results$p_adj_fdr <- p.adjust(pgls_results$p_value, method = "fdr")
}
write_output(pgls_results, out_pgls)

# ------------------- PLOTS -------------------
for (pred in broad_predictors) {
  for (resp in response_traits) {
    dat <- analysis_df %>%
      dplyr::select(tree_label, dplyr::all_of(pred), dplyr::all_of(resp)) %>%
      dplyr::filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
    counts <- table(dat[[pred]])
    if (nrow(dat) < min_total_tips || length(counts) < 2 || any(counts < min_group_n)) next
    p <- ggplot(dat, aes(x = .data[[pred]], y = .data[[resp]], fill = .data[[pred]])) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.40, color = "black") +
      geom_jitter(width = 0.08, height = 0, size = 2.4, alpha = 0.9, color = "black") +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      labs(
        title = paste(resp, "by", pred),
        x = pred,
        y = resp
      )
    out_png <- file.path(plot_dir, paste0("ecology_boxplot_", resp, "_by_", pred, ".png"))
    out_pdf <- file.path(plot_dir, paste0("ecology_boxplot_", resp, "_by_", pred, ".pdf"))
    ggsave(out_png, p, width = 6.2, height = 4.8, dpi = 400, bg = "white")
    ggsave(out_pdf, p, width = 6.2, height = 4.8, bg = "white")
  }
}

# ------------------- COMBINED PLOTS BY RESPONSE TRAIT -------------------
predictor_labels <- c(
  host_lineage_broad = "Host lineage",
  woody_association_broad = "Woody association",
  larval_lifestyle_broad = "Larval lifestyle",
  fungal_association_broad = "Fungal association"
)

response_labels <- c(
  PC1 = "PC1",
  PC2 = "PC2",
  abs_winding_angle_deg = "Absolute winding angle (deg)",
  axial_span = "Axial span",
  centroid_size = "Centroid size"
)

for (resp in response_traits) {
  combined_dat <- dplyr::bind_rows(lapply(broad_predictors, function(pred) {
    dat <- analysis_df %>%
      dplyr::select(tree_label, dplyr::all_of(pred), dplyr::all_of(resp)) %>%
      dplyr::filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
    counts <- table(dat[[pred]])
    if (nrow(dat) < min_total_tips || length(counts) < 2 || any(counts < min_group_n)) {
      return(data.frame())
    }
    data.frame(
      tree_label = dat$tree_label,
      predictor = pred,
      predictor_label = unname(predictor_labels[[pred]]),
      group = dat[[pred]],
      response = dat[[resp]],
      stringsAsFactors = FALSE
    )
  }))

  if (nrow(combined_dat) == 0) next

  combined_dat$predictor_label <- factor(
    combined_dat$predictor_label,
    levels = unname(predictor_labels[broad_predictors])
  )

  p_combined <- ggplot(combined_dat, aes(x = group, y = response, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.40, color = "black") +
    geom_jitter(width = 0.08, height = 0, size = 2.4, alpha = 0.9, color = "black") +
    facet_wrap(~ predictor_label, scales = "free_x", ncol = 2) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    ) +
    labs(
      title = paste("Ecology contrasts for", response_labels[[resp]]),
      x = NULL,
      y = response_labels[[resp]]
    )

  out_png_combined <- file.path(plot_dir, paste0("ecology_combined_boxplot_", resp, ".png"))
  out_pdf_combined <- file.path(plot_dir, paste0("ecology_combined_boxplot_", resp, ".pdf"))
  ggsave(out_png_combined, p_combined, width = 9.2, height = 7.0, dpi = 400, bg = "white")
  ggsave(out_pdf_combined, p_combined, width = 9.2, height = 7.0, bg = "white")
}

# ------------------- MANIFEST -------------------
manifest_lines <- c(
  "ECOLOGY ANALYSIS SCRIPT OUTPUTS",
  paste("Ecology table:", ecology_path),
  paste("PCM tip-level input:", pcm_tip_path),
  paste("Tree:", tree_path),
  "",
  "Broad ecology predictors tested:",
  paste("-", broad_predictors),
  "",
  "Response traits tested:",
  paste("-", response_traits),
  "",
  "Files written:",
  paste("-", out_input_merged),
  paste("-", out_factor_counts),
  paste("-", out_group_summary),
  paste("-", out_nonphylo),
  paste("-", out_phylanova),
  paste("-", out_pgls),
  paste("-", plot_dir)
)
writeLines(manifest_lines, out_manifest)

fmt_msg("Done. Ecology outputs written to: ", out_dir)
