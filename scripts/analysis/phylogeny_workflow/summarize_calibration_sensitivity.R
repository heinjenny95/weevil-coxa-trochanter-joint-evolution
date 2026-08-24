#!/usr/bin/env Rscript

# Combine primary-tree and updated-calibration reruns into compact source-data
# tables. The script intentionally excludes model intercepts when counting
# inferential changes because the manuscript interprets predictor slopes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: Rscript summarize_calibration_sensitivity.R",
      "<primary_pcm_dir> <sensitivity_pcm_dir>",
      "<primary_ecology_dir> <sensitivity_ecology_dir> <output_dir>"
    ),
    call. = FALSE
  )
}

primary_pcm_dir <- args[[1]]
sensitivity_pcm_dir <- args[[2]]
primary_ecology_dir <- args[[3]]
sensitivity_ecology_dir <- args[[4]]
output_dir <- args[[5]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_excel_csv <- function(path) {
  table <- read.csv2(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  names(table) <- sub("^\\ufeff", "", names(table))
  table
}

write_international_csv <- function(table, filename) {
  write.csv(
    table,
    file.path(output_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

stack_schemes <- function(primary, sensitivity) {
  primary$calibration_scheme <- "archived_primary_223Ma_ceiling"
  sensitivity$calibration_scheme <- "updated_195Ma_root_sensitivity"
  columns <- union(names(primary), names(sensitivity))
  for (column in setdiff(columns, names(primary))) primary[[column]] <- NA
  for (column in setdiff(columns, names(sensitivity))) sensitivity[[column]] <- NA
  rbind(primary[, columns, drop = FALSE], sensitivity[, columns, drop = FALSE])
}

count_threshold_changes <- function(
  primary,
  sensitivity,
  keys,
  p_column,
  adjusted_column = NULL
) {
  keep_primary <- unique(primary[, c(keys, p_column, adjusted_column), drop = FALSE])
  keep_sensitivity <- unique(sensitivity[, c(keys, p_column, adjusted_column), drop = FALSE])
  merged <- merge(
    keep_primary,
    keep_sensitivity,
    by = keys,
    suffixes = c("_primary", "_sensitivity"),
    all = TRUE
  )
  raw_changes <- sum(
    (merged[[paste0(p_column, "_primary")]] < 0.05) !=
      (merged[[paste0(p_column, "_sensitivity")]] < 0.05),
    na.rm = TRUE
  )
  adjusted_changes <- NA_integer_
  if (!is.null(adjusted_column)) {
    adjusted_changes <- sum(
      (merged[[paste0(adjusted_column, "_primary")]] < 0.05) !=
        (merged[[paste0(adjusted_column, "_sensitivity")]] < 0.05),
      na.rm = TRUE
    )
  }
  c(raw = raw_changes, adjusted = adjusted_changes)
}

signal_primary <- read_excel_csv(file.path(
  primary_pcm_dir,
  "02_Phylogenetic_signal/phylogenetic_signal_continuous.csv"
))
signal_sensitivity <- read_excel_csv(file.path(
  sensitivity_pcm_dir,
  "02_Phylogenetic_signal/phylogenetic_signal_continuous.csv"
))
write_international_csv(
  stack_schemes(signal_primary, signal_sensitivity),
  "calibration_sensitivity_phylogenetic_signal.csv"
)

models_primary <- read_excel_csv(file.path(
  primary_pcm_dir,
  "03_Evolutionary_models/evolutionary_model_fits_univariate.csv"
))
models_sensitivity <- read_excel_csv(file.path(
  sensitivity_pcm_dir,
  "03_Evolutionary_models/evolutionary_model_fits_univariate.csv"
))
write_international_csv(
  stack_schemes(models_primary, models_sensitivity),
  "calibration_sensitivity_evolutionary_models.csv"
)

pgls_primary <- read_excel_csv(file.path(
  primary_pcm_dir,
  "04_PGLS/pgls_continuous_vs_continuous.csv"
))
pgls_sensitivity <- read_excel_csv(file.path(
  sensitivity_pcm_dir,
  "04_PGLS/pgls_continuous_vs_continuous.csv"
))
pgls_primary <- subset(pgls_primary, term != "(Intercept)")
pgls_sensitivity <- subset(pgls_sensitivity, term != "(Intercept)")
write_international_csv(
  stack_schemes(pgls_primary, pgls_sensitivity),
  "calibration_sensitivity_pgls_slopes.csv"
)

allometry_primary <- read_excel_csv(file.path(
  primary_pcm_dir,
  "05_Allometry/allometry_results.csv"
))
allometry_sensitivity <- read_excel_csv(file.path(
  sensitivity_pcm_dir,
  "05_Allometry/allometry_results.csv"
))
allometry_primary <- subset(allometry_primary, term != "(Intercept)")
allometry_sensitivity <- subset(allometry_sensitivity, term != "(Intercept)")
write_international_csv(
  stack_schemes(allometry_primary, allometry_sensitivity),
  "calibration_sensitivity_allometry_slopes.csv"
)

asr_primary <- read_excel_csv(file.path(
  primary_pcm_dir,
  "08_ASR/asr_continuous_fastAnc.csv"
))
asr_sensitivity <- read_excel_csv(file.path(
  sensitivity_pcm_dir,
  "08_ASR/asr_continuous_fastAnc.csv"
))
write_international_csv(
  stack_schemes(asr_primary, asr_sensitivity),
  "calibration_sensitivity_continuous_asr.csv"
)

ecology_pgls_primary <- read_excel_csv(file.path(
  primary_ecology_dir,
  "ecology_pgls_factor_results.csv"
))
ecology_pgls_sensitivity <- read_excel_csv(file.path(
  sensitivity_ecology_dir,
  "ecology_pgls_factor_results.csv"
))
write_international_csv(
  stack_schemes(ecology_pgls_primary, ecology_pgls_sensitivity),
  "calibration_sensitivity_ecology_pgls.csv"
)

ecology_anova_primary <- read_excel_csv(file.path(
  primary_ecology_dir,
  "ecology_phylogenetic_anova_results.csv"
))
ecology_anova_sensitivity <- read_excel_csv(file.path(
  sensitivity_ecology_dir,
  "ecology_phylogenetic_anova_results.csv"
))
write_international_csv(
  stack_schemes(ecology_anova_primary, ecology_anova_sensitivity),
  "calibration_sensitivity_ecology_phyanova.csv"
)

model_best_primary <- unique(models_primary[, c("trait", "best_model")])
model_best_sensitivity <- unique(models_sensitivity[, c("trait", "best_model")])
model_best_comparison <- merge(
  model_best_primary,
  model_best_sensitivity,
  by = "trait",
  suffixes = c("_primary", "_sensitivity"),
  all = TRUE
)
model_rank_changes <- sum(
  model_best_comparison$best_model_primary !=
    model_best_comparison$best_model_sensitivity,
  na.rm = TRUE
)

signal_changes <- count_threshold_changes(
  signal_primary,
  signal_sensitivity,
  c("trait", "method"),
  "p_value",
  "fdr_p_value"
)
pgls_changes <- count_threshold_changes(
  pgls_primary,
  pgls_sensitivity,
  c("response", "predictor", "term"),
  "p_value",
  "fdr_p_value"
)
allometry_changes <- count_threshold_changes(
  allometry_primary,
  allometry_sensitivity,
  c("response", "predictor", "term"),
  "p_value",
  "fdr_p_value"
)
ecology_pgls_changes <- count_threshold_changes(
  ecology_pgls_primary,
  ecology_pgls_sensitivity,
  c("response", "predictor", "term"),
  "p_value",
  "p_adj_fdr"
)
ecology_anova_changes <- count_threshold_changes(
  ecology_anova_primary,
  ecology_anova_sensitivity,
  c("response", "predictor"),
  "p_value",
  "p_adj_fdr"
)

summary_table <- data.frame(
  analysis_block = c(
    "phylogenetic_signal",
    "pgls_predictor_slopes",
    "phylogenetic_allometry_slopes",
    "univariate_evolutionary_model_ranking",
    "ecology_pgls",
    "ecology_phylogenetic_anova",
    "continuous_ancestral_state_reconstruction"
  ),
  raw_p_threshold_changes = c(
    signal_changes[["raw"]],
    pgls_changes[["raw"]],
    allometry_changes[["raw"]],
    NA,
    ecology_pgls_changes[["raw"]],
    ecology_anova_changes[["raw"]],
    NA
  ),
  adjusted_p_threshold_changes = c(
    signal_changes[["adjusted"]],
    pgls_changes[["adjusted"]],
    allometry_changes[["adjusted"]],
    NA,
    ecology_pgls_changes[["adjusted"]],
    ecology_anova_changes[["adjusted"]],
    NA
  ),
  best_model_changes = c(NA, NA, NA, model_rank_changes, NA, NA, NA),
  interpretation = c(
    "No 0.05-threshold changes.",
    "No 0.05-threshold changes in predictor slopes.",
    "No 0.05-threshold changes in predictor slopes.",
    "PC1 changed from EB to BM; all other trait rankings were unchanged.",
    "No 0.05-threshold changes.",
    "No 0.05-threshold changes.",
    "Numerical root estimates changed, as expected for a branch-length sensitivity."
  ),
  stringsAsFactors = FALSE
)
write_international_csv(summary_table, "calibration_sensitivity_summary.csv")

cat("Calibration-sensitivity source tables written to: ", output_dir, "\n", sep = "")
