#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript export_final_allometry_tables.R ",
    "<bounded_pgls_allometry.csv> <output_directory>"
  )
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = TRUE)

results <- read.csv2(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c(
  "term", "estimate", "std_error", "statistic", "p_value",
  "response", "predictor", "n_taxa", "lambda", "fdr_p_value"
)
missing_columns <- setdiff(required, names(results))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

traits <- c("PC1", "PC2", "PC3", "PC4", "PC5", "abs_winding_angle_deg", "axial_span")
summary <- results[
  results$term == results$predictor &
    results$predictor == "centroid_size" &
    results$response %in% traits,
  required
]

summary <- summary[match(traits, summary$response), ]
if (anyNA(summary$response)) {
  stop("One or more expected allometry responses are absent from the input table.")
}

summary$trait <- c(
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "Absolute winding angle", "Axial span"
)
summary$significant_raw <- summary$p_value < 0.05
summary$significant_fdr <- summary$fdr_p_value < 0.05

summary <- summary[, c(
  "trait", "response", "predictor", "n_taxa", "estimate", "std_error",
  "statistic", "p_value", "fdr_p_value", "lambda",
  "significant_raw", "significant_fdr"
)]
names(summary)[names(summary) == "statistic"] <- "t_value"

write.csv2(
  summary,
  file.path(output_dir, "pgls_results_main_traits.csv"),
  row.names = FALSE,
  quote = TRUE
)

file.copy(
  input_file,
  file.path(output_dir, "pgls_allometry_all_traits.csv"),
  overwrite = TRUE
)

message("Wrote final bounded-lambda allometry tables to: ", output_dir)
