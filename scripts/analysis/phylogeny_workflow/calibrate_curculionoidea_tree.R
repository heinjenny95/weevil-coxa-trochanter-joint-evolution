#!/usr/bin/env Rscript

# Calibrate the Curculionoidea phylogeny used for comparative analyses.
#
# This script starts from an inferred maximum-likelihood tree, roots it on the
# specified nemonychid outgroup, inserts Caridae at the backbone position used
# in the manuscript, assigns Grafen branch lengths, and time-calibrates the tree
# with explicitly source-separated fossil minima and upper bounds.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript calibrate_curculionoidea_tree.R <input_ml_tree> <output_tree>",
      "[outgroup_tip] [caridae_tip] [curculionoidea_max_ma]",
      "[curculionidae_max_ma] [calibration_table_csv]",
      "[curculionoidea_min_ma] [curculionidae_min_ma]",
      "[curculionoidea_start_ma] [curculionidae_start_ma]",
      "[lambda] [model] [dual_iter_max]"
    ),
    call. = FALSE
  )
}

input_tree <- args[[1]]
output_tree <- args[[2]]
outgroup_tip <- ifelse(
  length(args) >= 3,
  args[[3]],
  "Nemonychidae___Rhynchitomacerinus"
)
caridae_tip <- ifelse(length(args) >= 4, args[[4]], "Caridae___Car")
curculionoidea_max_ma <- ifelse(length(args) >= 5, as.numeric(args[[5]]), 195.0)
curculionidae_max_ma <- ifelse(
  length(args) >= 6,
  as.numeric(args[[6]]),
  151.0
)
calibration_table_csv <- ifelse(
  length(args) >= 7,
  args[[7]],
  paste0(output_tree, ".calibration.csv")
)

curculionoidea_min_ma <- ifelse(length(args) >= 8, as.numeric(args[[8]]), 195.0)
curculionidae_min_ma <- ifelse(length(args) >= 9, as.numeric(args[[9]]), 113.0)
curculionoidea_start_ma <- ifelse(
  length(args) >= 10,
  as.numeric(args[[10]]),
  mean(c(curculionoidea_min_ma, curculionoidea_max_ma))
)
curculionidae_start_ma <- ifelse(
  length(args) >= 11,
  as.numeric(args[[11]]),
  mean(c(curculionidae_min_ma, curculionidae_max_ma))
)
chronos_lambda <- ifelse(length(args) >= 12, as.numeric(args[[12]]), 1)
chronos_model <- ifelse(length(args) >= 13, args[[13]], "correlated")
dual_iter_max <- ifelse(length(args) >= 14, as.integer(args[[14]]), 1000L)

if (
  !is.finite(curculionoidea_max_ma) ||
    !is.finite(curculionidae_max_ma) ||
    !is.finite(curculionoidea_min_ma) ||
    !is.finite(curculionidae_min_ma) ||
    !is.finite(curculionoidea_start_ma) ||
    !is.finite(curculionidae_start_ma) ||
    !is.finite(chronos_lambda) ||
    !is.finite(dual_iter_max) ||
    curculionoidea_max_ma < curculionoidea_min_ma ||
    curculionidae_max_ma < curculionidae_min_ma ||
    curculionidae_max_ma > curculionoidea_max_ma ||
    curculionoidea_start_ma < curculionoidea_min_ma ||
    curculionoidea_start_ma > curculionoidea_max_ma ||
    curculionidae_start_ma < curculionidae_min_ma ||
    curculionidae_start_ma > curculionidae_max_ma ||
    curculionidae_start_ma >= curculionoidea_start_ma ||
    chronos_lambda <= 0 ||
    dual_iter_max < 1
) {
  stop("Invalid or phylogenetically inconsistent calibration bounds.", call. = FALSE)
}

chronos_model <- match.arg(
  tolower(chronos_model),
  c("correlated", "relaxed", "discrete", "clock")
)

tree_ml <- read.tree(input_tree)

if (!outgroup_tip %in% tree_ml$tip.label) {
  stop("Outgroup tip not found in input tree: ", outgroup_tip, call. = FALSE)
}

tree_rooted <- root(tree_ml, outgroup = outgroup_tip, resolve.root = TRUE)
tips <- tree_rooted$tip.label

brentidae_tips <- grep("Brentidae", tips, value = TRUE)
curculionidae_tips <- grep("Curculionidae", tips, value = TRUE)

if (length(brentidae_tips) == 0 || length(curculionidae_tips) == 0) {
  stop("Could not identify Brentidae and/or Curculionidae tips.", call. = FALSE)
}

node_brentidae_curculionidae <- getMRCA(
  tree_rooted,
  c(brentidae_tips, curculionidae_tips)
)

# Caridae is absent from the source alignment and is therefore inserted at the
# backbone position supported by the phylogenetic literature used in the study.
tree_with_caridae <- bind.tip(
  tree_rooted,
  tip.label = caridae_tip,
  where = node_brentidae_curculionidae,
  position = 0.0001
)

tree_grafen <- compute.brlen(tree_with_caridae, method = "Grafen")
tips <- tree_grafen$tip.label

node_curculionoidea <- getMRCA(tree_grafen, tips)
node_curculionidae <- getMRCA(
  tree_grafen,
  grep("Curculionidae", tips, value = TRUE)
)

# Minimum ages follow the fossil-bearing-stratum ages tabulated by McKenna et
# al. (2019; PNAS, doi:10.1073/pnas.1909655116, Supplementary Table S5).
# The 223 Ma Curculionoidea ceiling was previously used as a conservative upper
# bound in weevil divergence dating by Letsch et al. (2020; Systematic
# Entomology, doi:10.1111/syen.12396), based on the then-interpreted oldest
# polyphagan fossil Leehermania prorova. Leehermania was subsequently placed in
# Myxophaga (Fikacek et al. 2020, doi:10.1111/syen.12386), so 223 Ma is used
# only when explicitly supplied as a historical sensitivity ceiling. The 195 Ma
# Curculionoidea ceiling follows current Belidae dating practice (Li et al.
# 2024, doi:10.7554/eLife.97552.3), and 151 Ma reproduces the tighter
# Curculionidae maximum described by Letsch et al.
calibration <- makeChronosCalib(
  tree_grafen,
  node = c(node_curculionoidea, node_curculionidae),
  age.min = c(curculionoidea_min_ma, curculionidae_min_ma),
  age.max = c(curculionoidea_max_ma, curculionidae_max_ma)
)

calibration$clade <- c("Curculionoidea", "Curculionidae")
calibration$age.start <- c(curculionoidea_start_ma, curculionidae_start_ma)
calibration$bound_role <- c(
  if (identical(curculionoidea_min_ma, curculionoidea_max_ma)) {
    "fixed root age for calibration sensitivity"
  } else {
    "root minimum and maximum"
  },
  if (identical(curculionidae_max_ma, curculionoidea_max_ma)) {
    "minimum with non-restrictive root ceiling"
  } else {
    "minimum and sensitivity maximum"
  }
)

control <- chronos.control(
  tol = 1e-8,
  iter.max = 10000,
  eval.max = 10000,
  dual.iter.max = dual_iter_max
)

tree_calibrated <- chronos(
  tree_grafen,
  calibration = calibration,
  lambda = chronos_lambda,
  model = chronos_model,
  control = control
)

write.tree(tree_calibrated, output_tree)
write.csv(calibration, calibration_table_csv, row.names = FALSE)

diagnostic_table_csv <- paste0(output_tree, ".fit_diagnostics.csv")
write.csv(
  data.frame(
    input_tree = basename(input_tree),
    output_tree = basename(output_tree),
    model = chronos_model,
    lambda = chronos_lambda,
    dual_iter_max = dual_iter_max,
    converged = isTRUE(attr(tree_calibrated, "convergence")),
    optimizer_message = as.character(attr(tree_calibrated, "message")),
    dual_iterations = as.integer(attr(tree_calibrated, "niter")),
    stringsAsFactors = FALSE
  ),
  diagnostic_table_csv,
  row.names = FALSE
)

node_depth <- node.depth.edgelength(tree_calibrated)
root_age <- max(node_depth[seq_len(Ntip(tree_calibrated))])
curculionidae_age <- root_age - node_depth[node_curculionidae]

message("Wrote calibrated tree: ", output_tree)
message("Wrote calibration table: ", calibration_table_csv)
message("Wrote fit diagnostics: ", diagnostic_table_csv)
message("Ultrametric: ", is.ultrametric(tree_calibrated))
message("Curculionoidea root age (Ma): ", signif(root_age, 8))
message("Curculionidae crown age (Ma): ", signif(curculionidae_age, 8))
