#!/usr/bin/env Rscript

# Calibrate the Curculionoidea phylogeny used for comparative analyses.
#
# This script starts from an inferred maximum-likelihood tree, roots it on the
# specified nemonychid outgroup, inserts Caridae at the backbone position used
# in the manuscript, assigns Grafen branch lengths, and time-calibrates the tree
# with fossil constraints from the McKenna et al. source dataset.

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
      "[outgroup_tip] [caridae_tip]"
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

calibration <- makeChronosCalib(
  tree_grafen,
  node = c(node_curculionoidea, node_curculionidae),
  age.min = c(157.3, 113.0),
  age.max = c(223.0, 223.0)
)

control <- chronos.control(tol = 1e-8, iter.max = 1000)

tree_calibrated <- chronos(
  tree_grafen,
  calibration = calibration,
  lambda = 1,
  model = "correlated",
  control = control
)

write.tree(tree_calibrated, output_tree)

message("Wrote calibrated tree: ", output_tree)
message("Ultrametric: ", is.ultrametric(tree_calibrated))
