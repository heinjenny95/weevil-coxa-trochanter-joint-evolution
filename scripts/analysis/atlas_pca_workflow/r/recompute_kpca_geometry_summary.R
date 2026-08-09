#!/usr/bin/env Rscript

# Recompute the kPCA shape--geometry sensitivity table from released matched
# kernel scores while applying the same winding-angle eligibility criterion as
# the main 60-specimen analysis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop(
  "Usage: recompute_kpca_geometry_summary.R <shape-data.csv> <linear-scores.csv> <legacy-scores.csv> <adaptive-scores.csv> <output.csv> <min-angle>"
)

shape <- read.table(args[[1]], header = TRUE, sep = ";", quote = "\"", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
linear <- read.table(args[[2]], header = TRUE, sep = ";", quote = "\"", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, dec = ",")
legacy <- read.csv(args[[3]], stringsAsFactors = FALSE, check.names = FALSE)
adaptive <- read.csv(args[[4]], stringsAsFactors = FALSE, check.names = FALSE)
cutoff <- as.numeric(args[[6]])

names(shape)[names(shape) == "angle_abs"] <- "abs_winding_angle_deg"
names(shape)[names(shape) == "axial_metric"] <- "axial_span"
numeric_shape <- c("PC1", "PC2", "abs_winding_angle_deg", "axial_span", "axial_pitch")
shape[numeric_shape] <- lapply(shape[numeric_shape], as.numeric)

score_cols <- function(x) {
  preferred <- grep("^(matched_kPC|kPC|PC)[0-9]+$", names(x), value = TRUE)
  if (length(preferred) < 5) stop("Could not identify five score columns")
  preferred[seq_len(5)]
}

merge_scores <- function(base, scores, prefix) {
  cols <- score_cols(scores)
  z <- scores[, c("specimen_id", cols), drop = FALSE]
  names(z)[-1] <- paste0(prefix, seq_along(cols))
  merge(base, z, by = "specimen_id", all.x = TRUE, sort = FALSE)
}

d <- merge(shape[, c("specimen_id", "abs_winding_angle_deg", "axial_span", "axial_pitch")], linear[, c("specimen_id", paste0("PC", 1:5))], by = "specimen_id", all.x = TRUE, sort = FALSE)
d <- merge_scores(d, legacy, "L")
d <- merge_scores(d, adaptive, "A")
d <- d[is.finite(d$abs_winding_angle_deg) & d$abs_winding_angle_deg >= cutoff, ]

methods <- list(
  "Linear PCA" = c("PC1", "PC2", "PC3", "PC4", "PC5"),
  "Legacy kPCA" = paste0("L", 1:5),
  "Adaptive kPCA" = paste0("A", 1:5)
)

one <- function(response, method, axes, n_axes) {
  f <- reformulate(axes[seq_len(n_axes)], response = response)
  fit <- lm(f, data = d)
  s <- summary(fit)
  data.frame(
    response = response,
    method = method,
    predictor_axes = n_axes,
    minimum_abs_winding_angle_deg = cutoff,
    n = nobs(fit),
    R2 = unname(s$r.squared),
    F = unname(s$fstatistic[[1]]),
    df1 = unname(s$fstatistic[[2]]),
    df2 = unname(s$fstatistic[[3]]),
    P = pf(s$fstatistic[[1]], s$fstatistic[[2]], s$fstatistic[[3]], lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (response in c("abs_winding_angle_deg", "axial_span", "axial_pitch")) {
  for (method in names(methods)) {
    for (n_axes in c(2, 5)) rows[[length(rows) + 1]] <- one(response, method, methods[[method]], n_axes)
  }
}
out <- do.call(rbind, rows)
dir.create(dirname(args[[5]]), recursive = TRUE, showWarnings = FALSE)
write.csv(out, args[[5]], row.names = FALSE)
