# Full atlas-shape allometry against an isometric null.
#
# Atlas shapes were normalized for scale before PCA. Under isometry, the
# multivariate regression of shape on log centroid size therefore has a zero
# slope. The primary test uses every non-zero atlas PC axis, which is
# equivalent to the complete atlas-derived shape representation. PC1-PC5 are
# retained only as an anatomical sensitivity analysis.

suppressPackageStartupMessages({
  library(RRPP)
  library(geomorph)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    paste(
      "Usage: Rscript test_full_shape_allometry.R",
      "<pca_with_centroid_size.csv> <output_directory>"
    ),
    call. = FALSE
  )
}

input_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_permutations <- 9999
p_adjust_method <- "holm"
interpretation_pcs <- paste0("PC", 1:5)

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
  data
}

as_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}

sort_pc_names <- function(x) {
  x[order(as.integer(sub("^PC", "", x)))]
}

extract_rrpp_table <- function(fit, scope, pc_names, variance_percent) {
  fit_anova <- anova(fit)
  tab <- if (!is.null(fit_anova$table)) fit_anova$table else as.data.frame(fit_anova)
  tab <- as.data.frame(tab)
  tab$Term <- rownames(tab)
  rownames(tab) <- NULL
  tab <- tab[, c("Term", setdiff(names(tab), "Term")), drop = FALSE]
  tab$analysis_scope <- scope
  tab$n_pc_axes <- length(pc_names)
  tab$pc_axes <- paste(pc_names, collapse = ",")
  tab$variance_represented_percent <- variance_percent
  tab$null_hypothesis <- "zero multivariate shape slope under isometry"
  tab$n_permutations <- n_permutations
  tab
}

extract_procd_table <- function(fit, scope, pc_names, variance_percent) {
  fit_anova <- anova(fit)
  tab <- NULL
  if (!is.null(fit_anova$ANOVA)) tab <- fit_anova$ANOVA
  if (is.null(tab) && !is.null(fit_anova$aov.table)) tab <- fit_anova$aov.table
  if (is.null(tab) && !is.null(fit_anova$table)) tab <- fit_anova$table
  if (is.null(tab)) stop("Could not extract procD.lm ANOVA table.", call. = FALSE)
  tab <- as.data.frame(tab)
  tab$Term <- rownames(tab)
  rownames(tab) <- NULL
  wanted <- c("Term", "Df", "SS", "MS", "Rsq", "F", "Z", "Pr(>F)")
  tab <- tab[, c(intersect(wanted, names(tab)), setdiff(names(tab), wanted)), drop = FALSE]
  tab$analysis_scope <- scope
  tab$n_pc_axes <- length(pc_names)
  tab$pc_axes <- paste(pc_names, collapse = ",")
  tab$variance_represented_percent <- variance_percent
  tab$null_hypothesis <- "zero multivariate shape slope under isometry"
  tab$n_permutations <- n_permutations
  tab
}

data <- read_table_auto(input_path)
if (!"centroid_size" %in% names(data) && !"logCS" %in% names(data)) {
  stop("Input requires centroid_size or logCS.", call. = FALSE)
}

pc_names <- sort_pc_names(grep("^PC[0-9]+$", names(data), value = TRUE))
if (length(pc_names) < 5) {
  stop("Fewer than five PC columns were detected.", call. = FALSE)
}

for (name in pc_names) data[[name]] <- as_numeric(data[[name]])
if ("logCS" %in% names(data)) {
  data$logCS <- as_numeric(data$logCS)
} else {
  data$centroid_size <- as_numeric(data$centroid_size)
  if (any(data$centroid_size <= 0, na.rm = TRUE)) {
    stop("centroid_size must be positive.", call. = FALSE)
  }
  data$logCS <- log(data$centroid_size)
}

pc_variance <- vapply(data[pc_names], stats::var, numeric(1), na.rm = TRUE)
nonzero_pcs <- pc_names[is.finite(pc_variance) & pc_variance > sqrt(.Machine$double.eps)]
if (length(nonzero_pcs) == 0) {
  stop("No non-zero PC axes were detected.", call. = FALSE)
}
if (!all(interpretation_pcs %in% nonzero_pcs)) {
  stop("PC1-PC5 are not all present among the non-zero axes.", call. = FALSE)
}

complete_rows <- stats::complete.cases(data[, c("logCS", nonzero_pcs), drop = FALSE])
analysis_data <- data[complete_rows, , drop = FALSE]
if (nrow(analysis_data) < 10) {
  stop("Fewer than ten complete specimens remain.", call. = FALSE)
}

pc_variance_complete <- vapply(
  analysis_data[nonzero_pcs],
  stats::var,
  numeric(1),
  na.rm = TRUE
)
total_variance <- sum(pc_variance_complete)
pc5_variance_percent <- 100 * sum(pc_variance_complete[interpretation_pcs]) / total_variance

full_shape <- as.matrix(analysis_data[, nonzero_pcs, drop = FALSE])
pc5_shape <- as.matrix(analysis_data[, interpretation_pcs, drop = FALSE])

set.seed(20260728)
fit_full <- RRPP::lm.rrpp(
  full_shape ~ logCS,
  data = analysis_data,
  iter = n_permutations,
  print.progress = FALSE
)
fit_pc5 <- RRPP::lm.rrpp(
  pc5_shape ~ logCS,
  data = analysis_data,
  iter = n_permutations,
  print.progress = FALSE
)

full_rrpp <- extract_rrpp_table(
  fit_full,
  "full_atlas_shape",
  nonzero_pcs,
  100
)
pc5_rrpp <- extract_rrpp_table(
  fit_pc5,
  "PC1-PC5_sensitivity",
  interpretation_pcs,
  pc5_variance_percent
)

set.seed(20260728)
fit_procd <- geomorph::procD.lm(
  full_shape ~ logCS,
  data = analysis_data,
  iter = n_permutations,
  print.progress = FALSE
)
full_procd <- extract_procd_table(
  fit_procd,
  "full_atlas_shape",
  nonzero_pcs,
  100
)

univariate <- lapply(nonzero_pcs, function(pc) {
  fit <- stats::lm(analysis_data[[pc]] ~ analysis_data$logCS)
  fit_summary <- summary(fit)
  coefficient <- fit_summary$coefficients[2, ]
  data.frame(
    trait = pc,
    estimate = coefficient[["Estimate"]],
    std_error = coefficient[["Std. Error"]],
    t_value = coefficient[["t value"]],
    p_value = coefficient[["Pr(>|t|)"]],
    r_squared = fit_summary$r.squared,
    adjusted_r_squared = fit_summary$adj.r.squared,
    n = stats::nobs(fit),
    stringsAsFactors = FALSE
  )
})
univariate <- do.call(rbind, univariate)
univariate$p_value_adjusted <- stats::p.adjust(
  univariate$p_value,
  method = p_adjust_method
)
univariate$p_adjust_method <- p_adjust_method
univariate$analysis_scope <- "all_nonzero_atlas_PCs"

pc5_univariate <- univariate[
  match(interpretation_pcs, univariate$trait),
  ,
  drop = FALSE
]

regression_fit <- stats::lm(full_shape ~ analysis_data$logCS)
regression_vector_raw <- stats::coef(regression_fit)[2, ]
regression_norm <- sqrt(sum(regression_vector_raw^2))
regression_vector <- regression_vector_raw / regression_norm
shape_center <- colMeans(full_shape)
centered_shape <- sweep(full_shape, 2, shape_center, "-")
regression_scores <- as.vector(centered_shape %*% regression_vector)

fitted_ss_by_pc <- regression_vector_raw^2
fitted_ss_total <- sum(fitted_ss_by_pc)
contributions <- data.frame(
  principal_component = nonzero_pcs,
  regression_coefficient = as.numeric(regression_vector_raw),
  normalized_regression_loading = as.numeric(regression_vector),
  fitted_effect_percent = 100 * fitted_ss_by_pc / fitted_ss_total,
  total_shape_variance_percent = 100 * pc_variance_complete[nonzero_pcs] / total_variance,
  stringsAsFactors = FALSE
)
contributions <- contributions[
  order(contributions$fitted_effect_percent, decreasing = TRUE),
  ,
  drop = FALSE
]

id_name <- names(analysis_data)[grepl("specimen_id", names(analysis_data))][1]
score_output <- data.frame(
  specimen_id = if (!is.na(id_name)) analysis_data[[id_name]] else seq_len(nrow(analysis_data)),
  logCS = analysis_data$logCS,
  full_shape_regression_score = regression_scores,
  stringsAsFactors = FALSE
)
if ("Family" %in% names(analysis_data)) score_output$Family <- analysis_data$Family
if ("Family.x" %in% names(analysis_data)) score_output$Family <- analysis_data$Family.x
if ("tree_label" %in% names(analysis_data)) score_output$tree_label <- analysis_data$tree_label

analysis_scope <- data.frame(
  primary_analysis = "RRPP regression of all non-zero atlas PCs on log centroid size",
  isometric_null = "zero multivariate shape slope because atlas shape is scale-normalized",
  n_specimens = nrow(analysis_data),
  n_pc_axes_primary = length(nonzero_pcs),
  primary_variance_percent = 100,
  sensitivity_analysis = "RRPP regression of PC1-PC5 on log centroid size",
  sensitivity_variance_percent = pc5_variance_percent,
  pc1_pc5_fitted_effect_percent = sum(
    contributions$fitted_effect_percent[
      contributions$principal_component %in% interpretation_pcs
    ]
  ),
  n_permutations = n_permutations,
  stringsAsFactors = FALSE
)

write_output <- function(data, filename) {
  utils::write.csv2(
    data,
    file.path(output_dir, filename),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

write_output(full_rrpp, "allometry_rrpp_multivariate_results.csv")
write_output(full_procd, "allometry_procD_lm_results.csv")
write_output(pc5_rrpp, "allometry_rrpp_PC1_to_PC5_sensitivity_results.csv")
write_output(univariate, "allometry_univariate_all_PC_results.csv")
write_output(pc5_univariate, "allometry_univariate_PC1_to_PC5_results.csv")
write_output(contributions, "allometry_PC_allometric_contributions.csv")
write_output(score_output, "allometry_full_shape_regression_scores.csv")
write_output(analysis_scope, "allometry_analysis_scope.csv")

message("Full-shape allometry outputs written to: ", output_dir)
