#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(caper)
  library(dplyr)
  library(readr)
})

root_dir <- "<MANUSCRIPT_PROJECT_ROOT>"
pgls_dir <- file.path(root_dir, "05_Results", "03_Allometry", "Allometry_Phylogenetic")

specimen_file <- file.path(pgls_dir, "pgls_input_specimen_level.csv")
aggregate_file <- file.path(pgls_dir, "pgls_input_aggregated_to_tree_tip.csv")
results_file <- file.path(pgls_dir, "pgls_results_main_traits.csv")
summary_file <- file.path(pgls_dir, "pgls_model_summaries.txt")
log_file <- file.path(pgls_dir, "pgls_log.txt")
tree_file <- file.path(pgls_dir, "pgls_tree_used.tre")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
backup_dir <- file.path(pgls_dir, "Backup", paste0("before_axial_span_phylo_allometry_", timestamp))
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(aggregate_file, results_file, summary_file, log_file)) {
  if (file.exists(path)) file.copy(path, backup_dir, overwrite = TRUE)
}

read_csv2_local <- function(path) {
  readr::read_delim(
    path,
    delim = ";",
    locale = readr::locale(decimal_mark = ",", grouping_mark = ""),
    show_col_types = FALSE,
    trim_ws = TRUE
  ) |>
    as.data.frame(check.names = FALSE)
}

write_csv2_local <- function(x, path) {
  utils::write.csv2(x, path, row.names = FALSE)
}

fit_pgls_row <- function(agg, tree, trait) {
  df <- agg |>
    dplyr::filter(is.finite(logCS), is.finite(.data[[trait]])) |>
    as.data.frame()

  tree_trait <- ape::drop.tip(tree, setdiff(tree$tip.label, df$tree_label))
  df <- df[match(tree_trait$tip.label, df$tree_label), , drop = FALSE]

  comp <- caper::comparative.data(
    phy = tree_trait,
    data = df,
    names.col = "tree_label",
    vcv = TRUE,
    warn.dropped = TRUE
  )

  fit <- caper::pgls(stats::as.formula(paste(trait, "~ logCS")), data = comp, lambda = 1)
  sm <- summary(fit)
  coef_tab <- sm$coefficients

  data.frame(
    trait = trait,
    n = nrow(df),
    intercept = unname(coef_tab["(Intercept)", "Estimate"]),
    slope = unname(coef_tab["logCS", "Estimate"]),
    std_error = unname(coef_tab["logCS", "Std. Error"]),
    t_value = unname(coef_tab["logCS", "t value"]),
    p_value = unname(coef_tab["logCS", "Pr(>|t|)"]),
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    lambda = 1,
    logLik = as.numeric(stats::logLik(fit)),
    AIC = stats::AIC(fit),
    model_fit_strategy = "lambda_fixed_1",
    stringsAsFactors = FALSE
  )
}

specimen <- read_csv2_local(specimen_file)
for (nm in c("logCS", "centroid_size", "PC1", "PC2", "abs_winding_angle_deg", "axial_span")) {
  if (nm %in% names(specimen)) specimen[[nm]] <- as.numeric(specimen[[nm]])
}

aggregate <- specimen |>
  dplyr::filter(tree_label %in% ape::read.tree(tree_file)$tip.label) |>
  dplyr::group_by(tree_label) |>
  dplyr::summarise(
    n_specimens = dplyr::n(),
    logCS = mean(logCS, na.rm = TRUE),
    centroid_size_mean = mean(centroid_size, na.rm = TRUE),
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    abs_winding_angle_deg = mean(abs_winding_angle_deg, na.rm = TRUE),
    axial_span = mean(axial_span, na.rm = TRUE),
    .groups = "drop"
  ) |>
  as.data.frame()

write_csv2_local(aggregate, aggregate_file)

tree <- ape::read.tree(tree_file)
axial_row <- fit_pgls_row(aggregate, tree, "axial_span")

results <- read_csv2_local(results_file)
for (nm in c("n", "intercept", "slope", "std_error", "t_value", "p_value", "r_squared", "adj_r_squared", "lambda", "logLik", "AIC", "p_value_adjusted")) {
  if (nm %in% names(results)) results[[nm]] <- as.numeric(results[[nm]])
}

results <- results |>
  dplyr::filter(trait != "axial_span") |>
  dplyr::bind_rows(axial_row) |>
  dplyr::mutate(p_value_adjusted = stats::p.adjust(p_value, method = "holm"))

write_csv2_local(results, results_file)

summary_lines <- c(
  if (file.exists(summary_file)) readLines(summary_file, warn = FALSE, encoding = "UTF-8") else character(),
  "",
  "============================================================",
  "Trait: axial_span",
  "Fit strategy: lambda_fixed_1",
  "============================================================",
  capture.output({
    df <- aggregate |>
      dplyr::filter(is.finite(logCS), is.finite(axial_span)) |>
      as.data.frame()
    tree_axial <- ape::drop.tip(tree, setdiff(tree$tip.label, df$tree_label))
    df <- df[match(tree_axial$tip.label, df$tree_label), , drop = FALSE]
    comp <- caper::comparative.data(tree_axial, df, names.col = "tree_label", vcv = TRUE, warn.dropped = TRUE)
    print(summary(caper::pgls(axial_span ~ logCS, data = comp, lambda = 1)))
  }),
  ""
)
writeLines(summary_lines, summary_file, useBytes = TRUE)

log_lines <- c(
  if (file.exists(log_file)) readLines(log_file, warn = FALSE, encoding = "UTF-8") else character(),
  paste0("Trait axial_span: aggregated rows with complete data = ", axial_row$n),
  "Trait axial_span: PGLS fitted with strategy lambda_fixed_1"
)
writeLines(log_lines, log_file, useBytes = TRUE)

message("Added axial_span to phylogenetic allometry.")
message("Backup: ", backup_dir)
print(axial_row)
