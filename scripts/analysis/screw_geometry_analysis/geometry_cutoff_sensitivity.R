#!/usr/bin/env Rscript

# Sensitivity of the principal shape--geometry models to the minimum absolute
# winding-angle threshold.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: geometry_cutoff_sensitivity.R <PCA-and-geometry.csv> <output.csv>")

d <- read.table(
  args[[1]], header = TRUE, sep = ";", quote = "\"", comment.char = "",
  stringsAsFactors = FALSE, check.names = FALSE, dec = ".", fileEncoding = "UTF-8"
)
if ("abs_winding_angle_deg" %in% names(d)) names(d)[names(d) == "abs_winding_angle_deg"] <- "angle_abs"
if ("axial_span" %in% names(d)) names(d)[names(d) == "axial_span"] <- "axial_metric"
numeric_cols <- c("PC1", "PC2", "angle_abs", "axial_metric")
d[numeric_cols] <- lapply(d[numeric_cols], function(x) as.numeric(gsub(",", ".", x, fixed = TRUE)))
if (!"axial_pitch" %in% names(d)) d$axial_pitch <- d$axial_metric * 360 / d$angle_abs

summarise_model <- function(model, threshold, response, formula_text) {
  s <- summary(model)
  data.frame(
    threshold_deg = threshold,
    response = response,
    formula = formula_text,
    n = nobs(model),
    r_squared = unname(s$r.squared),
    adjusted_r_squared = unname(s$adj.r.squared),
    F = unname(s$fstatistic[[1]]),
    df1 = unname(s$fstatistic[[2]]),
    df2 = unname(s$fstatistic[[3]]),
    model_p = pf(s$fstatistic[[1]], s$fstatistic[[2]], s$fstatistic[[3]], lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (threshold in c(0, 10, 20, 30, 45)) {
  x <- d[is.finite(d$angle_abs) & d$angle_abs >= threshold, ]
  rows[[length(rows) + 1]] <- summarise_model(lm(angle_abs ~ PC1 + PC2, data = x), threshold, "absolute_winding_angle", "angle_abs ~ PC1 + PC2")
  rows[[length(rows) + 1]] <- summarise_model(lm(axial_metric ~ angle_abs, data = x), threshold, "axial_span", "axial_metric ~ angle_abs")
  rows[[length(rows) + 1]] <- summarise_model(lm(axial_pitch ~ PC1 + PC2, data = x), threshold, "endpoint_equivalent_axial_pitch", "axial_pitch ~ PC1 + PC2")
}

out <- do.call(rbind, rows)
dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
write.table(out, args[[2]], sep = ";", row.names = FALSE, quote = FALSE, na = "", fileEncoding = "UTF-8")
