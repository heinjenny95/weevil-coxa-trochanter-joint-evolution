#!/usr/bin/env Rscript

# Family-level disparity in PC1-PC5.
#
# Outputs observed mean squared distances to family centroids, bootstrap
# confidence intervals, a global PERMDISP test, FDR-adjusted pairwise PERMDISP
# tests, and a rarefaction sensitivity analysis at the smallest retained family
# sample size. Families represented by fewer than three specimens are excluded.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(vegan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    paste(
      "Usage: Rscript calculate_family_disparity.R",
      "<pca_scores.csv> <specimen_key.csv> <output_dir>"
    )
  )
}

pca_path <- normalizePath(args[[1]], mustWork = TRUE)
key_path <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

k_pcs <- 5L
min_n <- 3L
n_boot <- 5000L
n_perm <- 9999L
n_rarefy <- 5000L
seed <- 20260729L

read_any <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  delim <- if (grepl(";", first, fixed = TRUE)) ";" else ","
  decimal <- if (delim == ";") "," else "."
  read_delim(
    path,
    delim = delim,
    locale = locale(decimal_mark = decimal, grouping_mark = ""),
    show_col_types = FALSE,
    trim_ws = TRUE
  )
}

write_table <- function(x, filename) {
  write_csv(x, file.path(out_dir, filename), na = "")
}

disparity_msd <- function(x) {
  centre <- colMeans(x)
  mean(rowSums(sweep(x, 2, centre, "-")^2))
}

bootstrap_disparity <- function(x, n_boot) {
  n <- nrow(x)
  replicate(
    n_boot,
    disparity_msd(x[sample.int(n, n, replace = TRUE), , drop = FALSE])
  )
}

rarefied_disparity <- function(x, target_n, n_iter) {
  if (nrow(x) == target_n) {
    return(rep(disparity_msd(x), n_iter))
  }
  replicate(
    n_iter,
    disparity_msd(x[sample.int(nrow(x), target_n, replace = FALSE), , drop = FALSE])
  )
}

extract_permdisp <- function(x, group, permutations) {
  fit <- betadisper(
    dist(x),
    group = droplevels(factor(group)),
    type = "centroid",
    bias.adjust = TRUE
  )
  test <- permutest(fit, permutations = permutations)
  tab <- as.data.frame(test$tab)
  tibble(
    F = unname(tab[1, "F"]),
    p_value = unname(tab[1, "Pr(>F)"]),
    permutations = permutations
  )
}

pca <- read_any(pca_path)
key <- read_any(key_path)

if (!"specimen_id" %in% names(pca)) stop("PCA table lacks specimen_id.")
if (!"specimen_id" %in% names(key)) stop("Specimen key lacks specimen_id.")

family_column <- intersect(c("Family", "family"), names(key))
if (!"Family" %in% names(pca)) {
  if (!length(family_column)) stop("No family column found in PCA table or key.")
  pca <- pca |>
    left_join(
      key |> transmute(specimen_id, Family = .data[[family_column[[1]]]]),
      by = "specimen_id"
    )
}

pc_columns <- grep("^PC[0-9]+$", names(pca), value = TRUE)
pc_columns <- pc_columns[order(as.integer(sub("^PC", "", pc_columns)))]
if (length(pc_columns) < k_pcs) stop("Fewer than five PC columns found.")
pc_columns <- pc_columns[seq_len(k_pcs)]

data <- pca |>
  mutate(
    Family = trimws(as.character(Family)),
    across(all_of(pc_columns), as.numeric)
  ) |>
  filter(
    !is.na(Family),
    Family != "",
    if_all(all_of(pc_columns), is.finite)
  )

family_counts <- data |>
  count(Family, name = "n") |>
  arrange(desc(n), Family)

retained_families <- family_counts |>
  filter(n >= min_n) |>
  pull(Family)

analysis_data <- data |>
  filter(Family %in% retained_families) |>
  mutate(Family = droplevels(factor(Family)))

if (nlevels(analysis_data$Family) < 2) {
  stop("Fewer than two families meet the minimum sample-size criterion.")
}

x <- as.matrix(analysis_data[, pc_columns])

set.seed(seed)
summary_table <- analysis_data |>
  group_by(Family) |>
  group_modify(\(group_data, family_key) {
    family_x <- as.matrix(group_data[, pc_columns])
    boot <- bootstrap_disparity(family_x, n_boot)
    tibble(
      n = nrow(family_x),
      disparity = disparity_msd(family_x),
      ci_low = unname(quantile(boot, 0.025)),
      ci_high = unname(quantile(boot, 0.975)),
      bootstrap_replicates = n_boot
    )
  }) |>
  ungroup() |>
  arrange(desc(disparity))

set.seed(seed)
global_test <- extract_permdisp(
  x,
  analysis_data$Family,
  permutations = n_perm
) |>
  mutate(
    test = "PERMDISP on Euclidean distances to family centroids",
    dimensions = paste0("PC1-PC", k_pcs),
    families = nlevels(analysis_data$Family),
    specimens = nrow(analysis_data),
    minimum_family_n = min_n,
    .before = 1
  )

family_pairs <- combn(levels(analysis_data$Family), 2, simplify = FALSE)
set.seed(seed)
pairwise_tests <- bind_rows(lapply(family_pairs, function(pair) {
  subset_data <- analysis_data |>
    filter(Family %in% pair) |>
    mutate(Family = droplevels(Family))
  pair_result <- extract_permdisp(
    as.matrix(subset_data[, pc_columns]),
    subset_data$Family,
    permutations = n_perm
  )
  pair_result |>
    mutate(
      family_1 = pair[[1]],
      family_2 = pair[[2]],
      n_1 = sum(subset_data$Family == pair[[1]]),
      n_2 = sum(subset_data$Family == pair[[2]]),
      .before = 1
    )
})) |>
  mutate(
    p_adjusted_fdr = p.adjust(p_value, method = "BH"),
    significant_fdr_0_05 = p_adjusted_fdr < 0.05
  ) |>
  arrange(p_adjusted_fdr, p_value)

target_n <- min(summary_table$n)
set.seed(seed)
rarefied_table <- analysis_data |>
  group_by(Family) |>
  group_modify(\(group_data, family_key) {
    family_x <- as.matrix(group_data[, pc_columns])
    draws <- rarefied_disparity(family_x, target_n, n_rarefy)
    tibble(
      original_n = nrow(family_x),
      rarefied_n = target_n,
      rarefied_disparity_mean = mean(draws),
      rarefied_disparity_median = median(draws),
      rarefied_ci_low = unname(quantile(draws, 0.025)),
      rarefied_ci_high = unname(quantile(draws, 0.975)),
      rarefaction_replicates = n_rarefy
    )
  }) |>
  ungroup() |>
  arrange(desc(rarefied_disparity_mean))

scope <- tibble(
  analysis = "Family-level disparity",
  response = "PC1-PC5",
  disparity_metric = "Mean squared Euclidean distance to the family centroid",
  global_test = "PERMDISP with centroid distances and bias adjustment",
  pairwise_test = "Pairwise PERMDISP with Benjamini-Hochberg FDR correction",
  sensitivity_analysis = paste0(
    "Rarefaction to ", target_n,
    " specimens per family without replacement"
  ),
  minimum_family_n = min_n,
  bootstrap_replicates = n_boot,
  permutation_replicates = n_perm,
  rarefaction_replicates = n_rarefy,
  random_seed = seed
)

write_table(family_counts, "family_disparity_sample_counts.csv")
write_table(summary_table, "family_disparity_summary_pc1_pc5.csv")
write_table(global_test, "family_disparity_global_test_pc1_pc5.csv")
write_table(pairwise_tests, "family_disparity_pairwise_tests_pc1_pc5.csv")
write_table(rarefied_table, "family_disparity_rarefied_n4_pc1_pc5.csv")
write_table(scope, "family_disparity_analysis_scope.csv")

# Compatibility filename used by the publication figure generator.
write_csv(
  summary_table,
  file.path(out_dir, "disparity_by_family_k5_nGT2.csv"),
  na = ""
)

cat("Retained families:\n")
print(summary_table)
cat("\nGlobal PERMDISP:\n")
print(global_test)
cat("\nPairwise PERMDISP:\n")
print(pairwise_tests)
cat("\nRarefied disparity:\n")
print(rarefied_table)
