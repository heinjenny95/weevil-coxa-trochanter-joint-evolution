# Phylogenetically informed allometry of the full atlas-shape score.
#
# The score projects all non-zero atlas PC axes onto the normalized
# specimen-level regression vector for log centroid size. The specimen-level
# RRPP analysis remains the primary full-shape significance test.

suppressPackageStartupMessages({
  library(ape)
  library(caper)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript fit_multivariate_allometry_score_pgls.R",
      "<specimen_table.csv> <tip_table.csv> <tree.tre> <output.csv>"
    ),
    call. = FALSE
  )
}

specimen_path <- args[[1]]
tip_path <- args[[2]]
tree_path <- args[[3]]
output_path <- args[[4]]

read_table_auto <- function(path) {
  header <- readLines(path, n = 1, warn = FALSE)
  semicolon <- grepl(";", header, fixed = TRUE)
  data <- utils::read.table(
    path,
    header = TRUE,
    sep = if (semicolon) ";" else ",",
    dec = if (semicolon) "," else ".",
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  names(data) <- sub("^\\ufeff", "", names(data))
  names(data)[grepl("tree_label$", names(data))] <- "tree_label"
  names(data)[grepl("specimen_id$", names(data))] <- "specimen_id"
  data
}

as_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}

normalize_tree_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- gsub("neydus", "nedyus", x, fixed = TRUE)
  x
}

prepare_table <- function(data, require_label = FALSE) {
  if (!"logCS" %in% names(data)) {
    if (!"centroid_size" %in% names(data)) {
      stop("Input requires logCS or centroid_size.", call. = FALSE)
    }
    data$logCS <- log(as_numeric(data$centroid_size))
  } else {
    data$logCS <- as_numeric(data$logCS)
  }

  if (require_label && !"tree_label" %in% names(data)) {
    if ("tree_tip" %in% names(data)) {
      data$tree_label <- data$tree_tip
    } else {
      stop("Tip-level input requires tree_label or tree_tip.", call. = FALSE)
    }
  }
  data
}

specimen <- prepare_table(read_table_auto(specimen_path), require_label = TRUE)
tip_data <- prepare_table(read_table_auto(tip_path), require_label = TRUE)
tree <- ape::read.tree(tree_path)
tree$node.label <- NULL
tree$tip.label <- normalize_tree_label(tree$tip.label)
tip_data$tree_label <- normalize_tree_label(tip_data$tree_label)

pc_names <- grep("^PC[0-9]+$", names(specimen), value = TRUE)
pc_names <- pc_names[order(as.integer(sub("^PC", "", pc_names)))]
if (length(pc_names) < 5) {
  stop("Fewer than five PC columns were detected.", call. = FALSE)
}
for (name in pc_names) specimen[[name]] <- as_numeric(specimen[[name]])
pc_variance <- vapply(specimen[pc_names], stats::var, numeric(1), na.rm = TRUE)
pc_names <- pc_names[
  is.finite(pc_variance) & pc_variance > sqrt(.Machine$double.eps)
]

specimen <- specimen[
  stats::complete.cases(specimen[, c(pc_names, "logCS", "tree_label")]),
]
tip_data <- tip_data[
  stats::complete.cases(tip_data[, c("logCS", "tree_label")]),
]

pc_matrix <- as.matrix(specimen[, pc_names])
pc_fit <- stats::lm(pc_matrix ~ specimen$logCS)
regression_vector <- stats::coef(pc_fit)[2, ]
regression_norm <- sqrt(sum(regression_vector^2))
if (!is.finite(regression_norm) || regression_norm == 0) {
  stop("The specimen-level regression vector has zero or invalid length.", call. = FALSE)
}
regression_vector <- regression_vector / regression_norm
specimen_center <- colMeans(pc_matrix)

project_score <- function(data) {
  centered <- sweep(as.matrix(data[, pc_names]), 2, specimen_center, "-")
  as.vector(centered %*% regression_vector)
}

specimen$multivariate_score <- project_score(specimen)
tip_scores <- stats::aggregate(
  specimen[, c("logCS", "multivariate_score")],
  by = list(tree_label = specimen$tree_label),
  FUN = mean,
  na.rm = TRUE
)
tip_data <- merge(
  tip_data[, c("tree_label", "logCS"), drop = FALSE],
  tip_scores[, c("tree_label", "multivariate_score"), drop = FALSE],
  by = "tree_label",
  all = FALSE
)

common_tips <- intersect(tree$tip.label, tip_data$tree_label)
if (length(common_tips) < 4) {
  stop("Fewer than four tree tips could be matched to the tip-level table.", call. = FALSE)
}
tree <- ape::keep.tip(tree, common_tips)
tip_data <- tip_data[match(tree$tip.label, tip_data$tree_label), ]
pgls_data <- tip_data[, c("tree_label", "logCS", "multivariate_score")]

comparative <- caper::comparative.data(
  phy = tree,
  data = pgls_data,
  names.col = "tree_label",
  vcv = TRUE,
  warn.dropped = TRUE
)
fit <- caper::pgls(multivariate_score ~ logCS, data = comparative, lambda = 1)
fit_summary <- summary(fit)
coefficient <- as.data.frame(fit_summary$coefficients)["logCS", , drop = FALSE]

results <- data.frame(
  response = "full_atlas_shape_regression_score",
  analysis_scope = "all_nonzero_atlas_PCs",
  n_pc_axes = length(pc_names),
  n = stats::nobs(fit),
  estimate = coefficient$Estimate,
  std_error = coefficient$`Std. Error`,
  t_value = coefficient$`t value`,
  p_value = coefficient$`Pr(>|t|)`,
  r_squared = fit_summary$r.squared,
  adjusted_r_squared = fit_summary$adj.r.squared,
  lambda = unname(fit$param["lambda"]),
  model = "PGLS_lambda_fixed_1",
  stringsAsFactors = FALSE
)

scores <- rbind(
  data.frame(level = "specimen", label = seq_len(nrow(specimen)),
             logCS = specimen$logCS, multivariate_score = specimen$multivariate_score),
  data.frame(level = "phylogenetic_tip", label = tip_data$tree_label,
             logCS = tip_data$logCS, multivariate_score = tip_data$multivariate_score)
)
vector_table <- data.frame(
  principal_component = pc_names,
  normalized_regression_loading = as.numeric(regression_vector),
  specimen_mean = as.numeric(specimen_center)
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv2(results, output_path, row.names = FALSE)
utils::write.csv2(scores, sub("\\.csv$", "_scores.csv", output_path), row.names = FALSE)
utils::write.csv2(vector_table, sub("\\.csv$", "_projection_vector.csv", output_path), row.names = FALSE)

message("Wrote PGLS results and projection diagnostics to: ", dirname(output_path))
