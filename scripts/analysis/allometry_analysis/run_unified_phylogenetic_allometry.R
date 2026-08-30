#!/usr/bin/env Rscript

# Unified univariate phylogenetic allometry for the final manuscript.
#
# The workflow uses one prespecified predictor and one fitting strategy:
#   - predictor: log of the arithmetic mean centroid size at each proxy tip;
#   - shape responses: all 68 atlas specimens aggregated to 15 proxy tips;
#   - geometry responses: the 63 main-dataset specimens matched before
#     aggregation to 14 proxy tips;
#   - Pagel's lambda: maximum likelihood constrained to [0, 1].

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: Rscript run_unified_phylogenetic_allometry.R",
      "<allometry_merged_table.csv> <specimen_key.csv>",
      "<primary_tree.tre> <output_directory>"
    ),
    call. = FALSE
  )
}

allometry_path <- normalizePath(args[[1]], mustWork = TRUE)
key_path <- normalizePath(args[[2]], mustWork = TRUE)
tree_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- normalizePath(args[[4]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

shape_responses <- paste0("PC", 1:5)
geometry_responses <- c(
  "abs_winding_angle_deg",
  "n_turns_abs",
  "start_end_dist",
  "axial_span",
  "axial_pitch_360",
  "fit_radius"
)
predictor_name <- "log_centroid_size"

read_table_auto <- function(path) {
  header <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
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

normalize_genus_tip <- function(x) {
  x <- trimws(as.character(x))
  x[x == "Neydus"] <- "Nedyus"
  x[x == "Belidae"] <- "Agnesiotis"
  x[x == "Caridae"] <- "Car"
  x
}

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

aggregate_tip_data <- function(data, responses, analysis_set) {
  needed <- c("tree_label", "centroid_size", responses)
  complete <- stats::complete.cases(data[, needed, drop = FALSE]) &
    data$centroid_size > 0
  data <- data[complete, needed, drop = FALSE]

  if (nrow(data) == 0L) {
    stop("No complete specimen rows for analysis set: ", analysis_set)
  }

  values <- stats::aggregate(
    data[, c("centroid_size", responses), drop = FALSE],
    by = list(tree_label = data$tree_label),
    FUN = mean_or_na
  )
  names(values)[names(values) == "centroid_size"] <- "centroid_size_mean"

  counts <- stats::aggregate(
    rep(1L, nrow(data)),
    by = list(tree_label = data$tree_label),
    FUN = sum
  )
  names(counts)[2] <- "n_specimens"

  values <- merge(values, counts, by = "tree_label", all.x = TRUE, sort = FALSE)
  values[[predictor_name]] <- log(values$centroid_size_mean)
  values$analysis_set <- analysis_set
  values
}

fit_bounded_pgls <- function(tree, tip_data, response, analysis_set) {
  needed <- c("tree_label", response, predictor_name)
  data <- tip_data[stats::complete.cases(tip_data[, needed, drop = FALSE]), needed]
  tree_used <- ape::keep.tip(tree, intersect(tree$tip.label, data$tree_label))
  data <- data[match(tree_used$tip.label, data$tree_label), , drop = FALSE]

  if (nrow(data) < 5L || anyNA(data$tree_label)) {
    stop("Too few matched tips for ", response, ".")
  }

  formula <- stats::as.formula(paste(response, "~", predictor_name))

  fit_at_lambda <- function(lambda) {
    correlation <- tryCatch(
      ape::corPagel(
        value = lambda,
        phy = tree_used,
        fixed = TRUE,
        form = ~tree_label
      ),
      error = function(e) NULL
    )
    if (is.null(correlation)) return(NULL)

    tryCatch(
      suppressWarnings(
        nlme::gls(
          model = formula,
          data = data,
          correlation = correlation,
          method = "ML",
          na.action = na.omit
        )
      ),
      error = function(e) NULL
    )
  }

  negative_log_likelihood <- function(lambda) {
    fit <- fit_at_lambda(lambda)
    if (is.null(fit)) return(Inf)
    value <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
    if (is.finite(value)) -value else Inf
  }

  optimum <- tryCatch(
    stats::optimize(negative_log_likelihood, interval = c(0, 1), tol = 1e-9),
    error = function(e) NULL
  )
  candidates <- c(0, 1)
  if (!is.null(optimum) && is.finite(optimum$minimum)) {
    candidates <- c(candidates, optimum$minimum)
  }
  objective <- vapply(candidates, negative_log_likelihood, numeric(1))
  if (all(!is.finite(objective))) {
    stop("No valid bounded-lambda PGLS fit for ", response, ".")
  }

  lambda <- candidates[which.min(objective)]
  fit <- fit_at_lambda(lambda)
  coefficient_table <- as.data.frame(summary(fit)$tTable)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL

  data.frame(
    analysis_set = analysis_set,
    term = coefficient_table$term,
    estimate = coefficient_table$Value,
    std_error = coefficient_table$`Std.Error`,
    statistic = coefficient_table$`t-value`,
    p_value = coefficient_table$`p-value`,
    response = response,
    predictor = predictor_name,
    n_taxa = nrow(data),
    lambda = lambda,
    lambda_fit_strategy = "maximum_likelihood_bounded_0_1",
    log_likelihood = as.numeric(stats::logLik(fit)),
    AIC = stats::AIC(fit),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

allometry <- read_table_auto(allometry_path)
key <- read_table_auto(key_path)
required_allometry <- c("specimen_id", "centroid_size", shape_responses, geometry_responses)
required_key <- c("specimen_id", "Family", "tree_tip")
missing_allometry <- setdiff(required_allometry, names(allometry))
missing_key <- setdiff(required_key, names(key))
if (length(missing_allometry) > 0L) {
  stop("Allometry table is missing: ", paste(missing_allometry, collapse = ", "))
}
if (length(missing_key) > 0L) {
  stop("Specimen key is missing: ", paste(missing_key, collapse = ", "))
}

for (name in c("centroid_size", shape_responses, geometry_responses)) {
  allometry[[name]] <- as_numeric(allometry[[name]])
}

key <- key[, required_key, drop = FALSE]
names(key)[names(key) == "Family"] <- "Family_key"
key$tree_tip <- normalize_genus_tip(key$tree_tip)
data <- merge(allometry, key, by = "specimen_id", all.x = TRUE, sort = FALSE)
family <- if ("Family" %in% names(data)) data$Family else data$Family_key
family[is.na(family) | !nzchar(family)] <- data$Family_key[is.na(family) | !nzchar(family)]
data$tree_label <- paste(trimws(family), data$tree_tip, sep = "___")

shape_tip_data <- aggregate_tip_data(data, shape_responses, "shape_only")
geometry_tip_data <- aggregate_tip_data(data, geometry_responses, "main_dataset")

tree <- ape::read.tree(tree_path)
tree$tip.label <- trimws(sub("^\\ufeff", "", tree$tip.label))
tree$node.label <- NULL

shape_tip_data <- shape_tip_data[shape_tip_data$tree_label %in% tree$tip.label, ]
geometry_tip_data <- geometry_tip_data[geometry_tip_data$tree_label %in% tree$tip.label, ]
shape_tip_data <- shape_tip_data[match(intersect(tree$tip.label, shape_tip_data$tree_label), shape_tip_data$tree_label), ]
geometry_tip_data <- geometry_tip_data[match(intersect(tree$tip.label, geometry_tip_data$tree_label), geometry_tip_data$tree_label), ]

if (nrow(shape_tip_data) != 15L) {
  stop("Expected 15 shape proxy tips, found ", nrow(shape_tip_data), ".")
}
if (nrow(geometry_tip_data) != 14L) {
  stop("Expected 14 geometry proxy tips, found ", nrow(geometry_tip_data), ".")
}

shape_results <- do.call(
  rbind,
  lapply(shape_responses, function(response) {
    fit_bounded_pgls(tree, shape_tip_data, response, "shape_only")
  })
)
geometry_results <- do.call(
  rbind,
  lapply(geometry_responses, function(response) {
    fit_bounded_pgls(tree, geometry_tip_data, response, "main_dataset")
  })
)
results <- rbind(shape_results, geometry_results)

results$fdr_family <- ifelse(
  results$analysis_set == "shape_only",
  "shape_PC1_to_PC5",
  "main_geometry_traits"
)
results$fdr_p_value <- NA_real_
is_slope <- results$term == results$predictor
for (family_name in unique(results$fdr_family)) {
  index <- which(is_slope & results$fdr_family == family_name)
  results$fdr_p_value[index] <- stats::p.adjust(results$p_value[index], method = "BH")
}

summary_responses <- c(shape_responses, "abs_winding_angle_deg", "axial_span")
summary_labels <- c(
  shape_responses,
  "Absolute winding angle",
  "Axial span"
)
summary_table <- results[
  is_slope & results$response %in% summary_responses,
  c(
    "analysis_set", "response", "predictor", "n_taxa", "estimate",
    "std_error", "statistic", "p_value", "fdr_p_value", "lambda",
    "lambda_fit_strategy", "fdr_family"
  )
]
summary_table <- summary_table[match(summary_responses, summary_table$response), ]
summary_table$trait <- summary_labels
summary_table$significant_raw <- summary_table$p_value < 0.05
summary_table$significant_fdr <- summary_table$fdr_p_value < 0.05
summary_table <- summary_table[, c(
  "trait", "analysis_set", "response", "predictor", "n_taxa",
  "estimate", "std_error", "statistic", "p_value", "fdr_p_value",
  "lambda", "lambda_fit_strategy", "fdr_family",
  "significant_raw", "significant_fdr"
)]
names(summary_table)[names(summary_table) == "estimate"] <- "slope"
names(summary_table)[names(summary_table) == "statistic"] <- "t_value"

geometry_table <- results[results$analysis_set == "main_dataset", ]

write_csv_lf <- function(data, path) {
  connection <- file(path, open = "wb")
  on.exit(close(connection))
  utils::write.csv(data, connection, row.names = FALSE, na = "", eol = "\n")
}

write_csv_lf(shape_tip_data, file.path(output_dir, "pgls_allometry_tip_data_shape.csv"))
write_csv_lf(geometry_tip_data, file.path(output_dir, "pgls_allometry_tip_data_geometry.csv"))
write_csv_lf(results, file.path(output_dir, "pgls_allometry_all_traits.csv"))
write_csv_lf(summary_table, file.path(output_dir, "pgls_results_main_traits.csv"))
write_csv_lf(geometry_table, file.path(output_dir, "pgls_geometry_main_traits_main_dataset.csv"))

message(
  "Unified phylogenetic allometry completed: ",
  nrow(shape_tip_data), " shape tips, ",
  nrow(geometry_tip_data), " geometry tips; ",
  sum(summary_table$significant_fdr), " FDR-supported main-trait slopes."
)
