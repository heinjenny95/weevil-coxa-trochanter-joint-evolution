#!/usr/bin/env Rscript

# Rebuild the concise robust shape--geometry result tables from the canonical
# primary (n = 63) and strict (n = 53) regression summaries. Model-level
# P values are kept distinct from the PC1 coefficient P values.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: refresh_shape_geometry_tables.R <repository-root>")
}

repo_root <- normalizePath(args[[1]], mustWork = TRUE)
data_dir <- file.path(repo_root, "data", "screw_geometry")

read_sc <- function(path) {
  read.table(
    path, header = TRUE, sep = ";", quote = "\"", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8"
  )
}

build_main <- function(regression_path, expected_full_n, expected_main_n) {
  reg <- read_sc(regression_path)
  keep <- reg$table_block == "regression_models" &
    reg$model %in% c("angle_abs ~ PC1 + PC2", "axial_pitch ~ PC1 + PC2") &
    reg$subset %in% c("full_dataset", "main_region_PC1_lt_0.1")
  out <- reg[keep, c(
    "model", "subset", "n", "r_squared", "p_model", "PC1_estimate", "PC1_p"
  )]
  names(out) <- c(
    "model", "subset", "n", "r_squared", "p_value", "PC1_effect", "PC1_p"
  )
  out$model <- sub("angle_abs", "Winding angle", out$model, fixed = TRUE)
  out$model <- sub("axial_pitch", "Axial pitch", out$model, fixed = TRUE)
  out$model <- sub("PC1 + PC2", "shape", out$model, fixed = TRUE)
  out$subset <- sub("full_dataset", "Full dataset", out$subset, fixed = TRUE)
  out$subset <- sub(
    "main_region_PC1_lt_0.1", "Main region (PC1 < 0.1)", out$subset,
    fixed = TRUE
  )
  expected_n <- ifelse(
    out$subset == "Main region (PC1 < 0.1)", expected_main_n, expected_full_n
  )
  if (!all(out$n == expected_n)) {
    stop("Unexpected sample size in ", regression_path, ": ",
         paste(unique(out$n), collapse = ", "))
  }
  out
}

primary <- build_main(file.path(data_dir, "regression_summary.csv"), 63, 55)
strict <- build_main(file.path(data_dir, "regression_summary_strict_good.csv"), 53, 53)

write.csv(
  primary, file.path(data_dir, "main_results.csv"),
  row.names = FALSE, na = "NA", fileEncoding = "UTF-8"
)
write.csv(
  strict, file.path(data_dir, "main_results_strict_good.csv"),
  row.names = FALSE, na = "NA", fileEncoding = "UTF-8"
)

cat("Primary rows:", nrow(primary), "(n = 63)\n")
cat("Strict rows:", nrow(strict), "(n = 53)\n")
