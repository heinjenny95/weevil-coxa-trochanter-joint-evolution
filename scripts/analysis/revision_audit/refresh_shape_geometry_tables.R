#!/usr/bin/env Rscript

# Build the two concise supplementary tables from the corrected 60-specimen
# geometry analysis. This prevents stale 59-specimen summaries from surviving
# in the submission package.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop(
  "Usage: refresh_shape_geometry_tables.R <shape-data.csv> <winding.csv> <joint-types.csv> <regression-summary.csv> <specimen-output.csv> <main-results-output.csv>"
)

read_sc <- function(path, dec = ".") {
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  read.table(
    path, header = TRUE, sep = ";", quote = "\"", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE, dec = dec,
    fileEncoding = "UTF-8", skip = ifelse(grepl("^sep=", first), 1, 0)
  )
}

shape <- read_sc(args[[1]])
winding <- read_sc(args[[2]], dec = ",")
joint <- read_sc(args[[3]])
reg <- read_sc(args[[4]])

alias <- c("308_lisshorhoptrus_oryzophilus_aligned" = "308_lisshorhoptrus_oryzophilus_trochanter_aligned")
hit <- winding$specimen_id %in% names(alias)
winding$specimen_id[hit] <- unname(alias[winding$specimen_id[hit]])

specimen <- merge(shape, winding[, c("specimen_id", "start_end_dist")], by = "specimen_id", all.x = TRUE, sort = FALSE)
specimen <- merge(specimen, joint[, c("specimen_id", "joint_type", "joint_type_strict", "screw_state")], by = "specimen_id", all.x = TRUE, sort = FALSE)
specimen <- specimen[, c(
  "specimen_id", "joint_type", "joint_type_strict", "screw_state",
  "angle_abs", "angle_signed", "axial_metric", "axial_pitch",
  "start_end_dist", "fit_radius", "fit_rms"
)]
names(specimen) <- c(
  "specimen_id", "joint_type", "joint_type_strict", "screw_state",
  "abs_winding_angle_deg", "signed_winding_angle_deg", "axial_span",
  "endpoint_equivalent_axial_pitch_360", "start_end_dist", "fit_radius", "radial_circle_fit_rms"
)

keep <- reg$table_block == "regression_models" & reg$model %in% c(
  "angle_abs ~ PC1 + PC2", "axial_pitch ~ PC1 + PC2"
)
main <- reg[keep, c(
  "model", "subset", "n", "r_squared", "adj_r_squared", "f_statistic",
  "df1", "df2", "p_model", "PC1_estimate", "PC1_p", "PC2_estimate", "PC2_p"
)]
main$quantity_note <- ifelse(
  grepl("axial_pitch", main$model),
  "endpoint-equivalent pitch is algebraically derived from angle and span",
  "measured winding angle"
)

dir.create(dirname(args[[5]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args[[6]]), recursive = TRUE, showWarnings = FALSE)
write.csv(specimen, args[[5]], row.names = FALSE, na = "")
write.csv(main, args[[6]], row.names = FALSE, na = "")

cat("Specimen rows:", nrow(specimen), "\n")
cat("Main-result rows:", nrow(main), "\n")
