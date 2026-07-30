# Direct phylogenetic multivariate allometry across the complete atlas-PC space.
#
# The inferential test uses all non-constant atlas PC axes as a multivariate
# response and the phylogenetic covariance matrix among matched tree tips.
# A one-dimensional projection is produced only for visualization and is not
# used as the response in the significance test.

suppressPackageStartupMessages({
  library(ape)
  library(geomorph)
  library(RRPP)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript run_phylogenetic_multivariate_allometry.R",
      "<specimen_table.csv> <specimen_key.csv> <tree.tre> <output_directory>"
    ),
    call. = FALSE
  )
}

specimen_path <- args[[1]]
key_path <- args[[2]]
tree_path <- args[[3]]
output_dir <- args[[4]]

set.seed(20260729)
permutations <- 9999

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

normalize_genus_tip <- function(x) {
  x <- trimws(as.character(x))
  x[x == "Neydus"] <- "Nedyus"
  x[x == "Belidae"] <- "Agnesiotis"
  x[x == "Caridae"] <- "Car"
  x
}

extract_anova_row <- function(fit, method_name) {
  table <- anova(fit)$table
  effect <- table["logCS", , drop = FALSE]
  data.frame(
    method = method_name,
    response = "all_nonconstant_atlas_PCs",
    n_tips = nrow(fit$LM$data),
    n_pc_axes = fit$LM$p,
    permutations = permutations,
    df = effect$Df,
    sum_squares = effect$SS,
    mean_square = effect$MS,
    r_squared = effect$Rsq,
    f_value = effect$F,
    z_score = effect$Z,
    p_value = effect$`Pr(>F)`,
    stringsAsFactors = FALSE
  )
}

specimen <- read_table_auto(specimen_path)
key <- read_table_auto(key_path)

required_specimen <- c("specimen_id", "Family", "centroid_size")
missing_specimen <- setdiff(required_specimen, names(specimen))
if (length(missing_specimen) > 0) {
  stop(
    "Specimen table is missing: ",
    paste(missing_specimen, collapse = ", "),
    call. = FALSE
  )
}
if (!all(c("specimen_id", "tree_tip") %in% names(key))) {
  stop("Specimen key requires specimen_id and tree_tip.", call. = FALSE)
}

pc_names <- grep("^PC[0-9]+$", names(specimen), value = TRUE)
pc_names <- pc_names[order(as.integer(sub("^PC", "", pc_names)))]
if (length(pc_names) < 5) {
  stop("Fewer than five atlas PC columns were detected.", call. = FALSE)
}
for (name in pc_names) specimen[[name]] <- as_numeric(specimen[[name]])
specimen$centroid_size <- as_numeric(specimen$centroid_size)

key_map <- key[, c("specimen_id", "tree_tip")]
names(key_map)[2] <- "tree_tip_key"
data <- merge(specimen, key_map, by = "specimen_id", all.x = TRUE)
data$tree_tip_key <- normalize_genus_tip(data$tree_tip_key)
data$tree_label <- paste(data$Family, data$tree_tip_key, sep = "___")

complete <- stats::complete.cases(
  data[, c(pc_names, "centroid_size", "tree_label")]
)
data <- data[complete & data$centroid_size > 0, , drop = FALSE]

pc_variance <- vapply(data[pc_names], stats::var, numeric(1), na.rm = TRUE)
informative_pc_names <- pc_names[
  is.finite(pc_variance) & pc_variance > sqrt(.Machine$double.eps)
]
if (length(informative_pc_names) < 2) {
  stop("Fewer than two non-constant PC axes remain.", call. = FALSE)
}

tip_data <- stats::aggregate(
  data[, c(informative_pc_names, "centroid_size")],
  by = list(tree_label = data$tree_label),
  FUN = mean,
  na.rm = TRUE
)
tip_counts <- stats::aggregate(
  data$specimen_id,
  by = list(tree_label = data$tree_label),
  FUN = length
)
names(tip_counts)[2] <- "n_specimens"
tip_data <- merge(tip_data, tip_counts, by = "tree_label", all.x = TRUE)

tree <- ape::read.tree(tree_path)
tree$node.label <- NULL
matched_tips <- intersect(tree$tip.label, tip_data$tree_label)
if (length(matched_tips) < 4) {
  stop("Fewer than four phylogenetic tips could be matched.", call. = FALSE)
}
tree <- ape::keep.tip(tree, matched_tips)
tip_data <- tip_data[match(tree$tip.label, tip_data$tree_label), , drop = FALSE]

Y <- as.matrix(tip_data[, informative_pc_names, drop = FALSE])
rownames(Y) <- tip_data$tree_label
logCS <- setNames(log(tip_data$centroid_size), tip_data$tree_label)
phylogenetic_covariance <- ape::vcv.phylo(tree, corr = FALSE)

rrpp_fit <- RRPP::lm.rrpp(
  Y ~ logCS,
  Cov = phylogenetic_covariance,
  iter = permutations,
  RRPP = TRUE,
  print.progress = FALSE
)
geomorph_fit <- geomorph::procD.pgls(
  Y ~ logCS,
  phy = tree,
  iter = permutations,
  print.progress = FALSE
)

results <- rbind(
  extract_anova_row(rrpp_fit, "RRPP_lm.rrpp_with_phylogenetic_covariance"),
  extract_anova_row(geomorph_fit, "geomorph_procD.pgls")
)
results$pc_axes_available <- length(pc_names)
results$pc_variance_represented_percent <- 100
results$predictor <- "log_centroid_size"
results$interpretation <- paste(
  "Direct multivariate phylogenetic test across the complete non-constant",
  "atlas-PC representation; visualization scores are descriptive only."
)

# Obtain the phylogenetic generalized least-squares direction in PC space.
X <- cbind(intercept = 1, logCS = as.numeric(logCS))
C_inverse <- solve(phylogenetic_covariance)
beta <- solve(t(X) %*% C_inverse %*% X) %*%
  t(X) %*% C_inverse %*% Y
allometric_direction <- as.numeric(beta["logCS", ])
direction_norm <- sqrt(sum(allometric_direction^2))
if (!is.finite(direction_norm) || direction_norm == 0) {
  stop("The phylogenetic allometric direction has zero length.", call. = FALSE)
}
allometric_direction <- allometric_direction / direction_norm
tip_center <- colMeans(Y)

project_score <- function(matrix) {
  centered <- sweep(matrix, 2, tip_center, "-")
  as.vector(centered %*% allometric_direction)
}

tip_scores <- project_score(Y)
specimen_matrix <- as.matrix(data[, informative_pc_names, drop = FALSE])
specimen_scores <- project_score(specimen_matrix)

family_by_tip <- stats::aggregate(
  data$Family,
  by = list(tree_label = data$tree_label),
  FUN = function(x) names(sort(table(x), decreasing = TRUE))[1]
)
names(family_by_tip)[2] <- "Family"
tip_score_table <- merge(
  data.frame(
    level = "phylogenetic_tip",
    label = tip_data$tree_label,
    tree_label = tip_data$tree_label,
    log_centroid_size = logCS,
    multivariate_allometry_score = tip_scores,
    n_specimens = tip_data$n_specimens,
    stringsAsFactors = FALSE
  ),
  family_by_tip,
  by = "tree_label",
  all.x = TRUE
)
specimen_score_table <- data.frame(
  level = "specimen",
  label = data$specimen_id,
  tree_label = data$tree_label,
  log_centroid_size = log(data$centroid_size),
  multivariate_allometry_score = specimen_scores,
  n_specimens = 1,
  Family = data$Family,
  stringsAsFactors = FALSE
)
score_table <- rbind(
  specimen_score_table[, names(tip_score_table)],
  tip_score_table
)

projection_table <- data.frame(
  principal_component = informative_pc_names,
  phylogenetic_gls_direction = allometric_direction,
  tip_mean = tip_center,
  stringsAsFactors = FALSE
)

matching_table <- data.frame(
  tree_tip = tree$tip.label,
  matched = tree$tip.label %in% tip_data$tree_label,
  n_specimens = tip_data$n_specimens,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  results,
  file.path(output_dir, "phylogenetic_multivariate_allometry_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  tip_data,
  file.path(output_dir, "phylogenetic_multivariate_allometry_tip_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  score_table,
  file.path(output_dir, "phylogenetic_multivariate_allometry_visualization_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  projection_table,
  file.path(output_dir, "phylogenetic_multivariate_allometry_projection.csv"),
  row.names = FALSE
)
utils::write.csv(
  matching_table,
  file.path(output_dir, "phylogenetic_multivariate_allometry_matching.csv"),
  row.names = FALSE
)

message(
  "Direct phylogenetic multivariate allometry completed for ",
  nrow(tip_data),
  " tips and ",
  length(informative_pc_names),
  " PC axes."
)
