#!/usr/bin/env Rscript

# Descriptive Hopkins clusterability analysis for the final PC1-PC5
# morphospace representation. The repeated randomized evaluations sample
# observed points and uniform reference points within each call; they do not
# bootstrap specimens or create additional independent observations.

suppressPackageStartupMessages(library(hopkins))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript run_hopkins_clusterability.R <PCA_scores.csv> <output_dir>"
  )
}

pca_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

get_H <- function(x) {
  if (is.list(x)) {
    if (!is.null(x$H)) return(as.numeric(x$H))
    return(as.numeric(unlist(x)[1]))
  }
  as.numeric(x)
}

df <- read.csv(
  pca_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
pc_names <- paste0("PC", 1:5)
if (!all(pc_names %in% names(df))) {
  stop("Input must contain PC1-PC5 columns")
}

X <- scale(as.matrix(df[, pc_names, drop = FALSE]))
n_specimens <- nrow(X)
m_primary <- floor(0.1 * n_specimens)
n_primary <- 100L
primary_seed <- 1L

# Preserve the historical RNG sequence: the archived workflow calculated one
# single evaluation before the 100 repeated evaluations.
set.seed(primary_seed)
single_evaluation <- get_H(hopkins::hopkins(X, m = m_primary))
primary_values <- replicate(
  n_primary,
  get_H(hopkins::hopkins(X, m = m_primary))
)

randomized_values <- data.frame(
  evaluation = seq_along(primary_values),
  hopkins = primary_values
)

primary_summary <- data.frame(
  dimensions = "PC1-PC5",
  n_specimens = n_specimens,
  m = m_primary,
  randomized_evaluations = n_primary,
  seed = primary_seed,
  single_evaluation = single_evaluation,
  mean = mean(primary_values),
  median = median(primary_values),
  sd = sd(primary_values),
  min = min(primary_values),
  max = max(primary_values),
  q025 = unname(quantile(primary_values, 0.025)),
  q975 = unname(quantile(primary_values, 0.975)),
  fraction_above_0_75 = mean(primary_values > 0.75)
)

m_sensitivity <- do.call(
  rbind,
  lapply(4:10, function(m_candidate) {
    sensitivity_seed <- 1000L + m_candidate
    set.seed(sensitivity_seed)
    values <- replicate(
      1000L,
      get_H(hopkins::hopkins(X, m = m_candidate))
    )
    data.frame(
      m = m_candidate,
      sample_fraction = m_candidate / n_specimens,
      randomized_evaluations = length(values),
      seed = sensitivity_seed,
      mean = mean(values),
      median = median(values),
      sd = sd(values),
      q025 = unname(quantile(values, 0.025)),
      q25 = unname(quantile(values, 0.25)),
      q75 = unname(quantile(values, 0.75)),
      q975 = unname(quantile(values, 0.975)),
      fraction_above_0_75 = mean(values > 0.75)
    )
  })
)

write.csv(
  randomized_values,
  file.path(out_dir, "hopkins_randomized_values.csv"),
  row.names = FALSE
)
write.csv(
  primary_summary,
  file.path(out_dir, "hopkins_summary.csv"),
  row.names = FALSE
)
write.csv(
  m_sensitivity,
  file.path(out_dir, "hopkins_m_sensitivity.csv"),
  row.names = FALSE
)

cat(
  sprintf(
    paste0(
      "PC1-PC5: n=%d, m=%d, evaluations=%d, seed=%d, ",
      "mean=%.15f, median=%.15f, range=%.15f-%.15f\n"
    ),
    n_specimens,
    m_primary,
    n_primary,
    primary_seed,
    mean(primary_values),
    median(primary_values),
    min(primary_values),
    max(primary_values)
  )
)
cat(
  sprintf(
    "m-sensitivity means %.3f-%.3f; fraction H > 0.75 %.3f-%.3f\n",
    min(m_sensitivity$mean),
    max(m_sensitivity$mean),
    min(m_sensitivity$fraction_above_0_75),
    max(m_sensitivity$fraction_above_0_75)
  )
)
