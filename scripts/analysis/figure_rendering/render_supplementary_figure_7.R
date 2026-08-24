suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
manuscript_root <- if (length(args) >= 2 && nzchar(args[[2]])) normalizePath(args[[2]], winslash = "/", mustWork = TRUE) else NA_character_

specimen_file <- file.path(repo_root, "data", "supplementary_source_data", "S00_Specimen_Metadata", "specimen_key.csv")
tree_file <- file.path(repo_root, "data", "phylogeny", "P01_Trees", "01_primary_tree_grafen.tre")
repo_out_dir <- file.path(repo_root, "data", "foundational_figures")
dir.create(repo_out_dir, recursive = TRUE, showWarnings = FALSE)

specimens <- read.csv(specimen_file, check.names = FALSE, stringsAsFactors = FALSE)
specimens$tree_tip <- trimws(specimens$tree_tip)
tip_corrections <- c("Neydus" = "Nedyus", "Belidae" = "Agnesiotis", "Caridae" = "Car")
replace_idx <- specimens$tree_tip %in% names(tip_corrections)
specimens$tree_tip[replace_idx] <- unname(tip_corrections[specimens$tree_tip[replace_idx]])
specimens$tree_label <- ifelse(
  is.na(specimens$Family) | is.na(specimens$tree_tip) | specimens$Family == "" | specimens$tree_tip == "",
  NA_character_, paste(specimens$Family, specimens$tree_tip, sep = "___")
)

opening <- specimens[["Coxal wall hole"]]
opening <- ifelse(is.na(opening) | trimws(as.character(opening)) == "", NA, toupper(trimws(as.character(opening))) == "TRUE")
specimens$opening <- opening
tip_groups <- split(specimens$opening, specimens$tree_label)
tip_states <- vapply(tip_groups, function(x) {
  values <- unique(stats::na.omit(x))
  if (length(values) == 1) values[[1]] else NA
}, logical(1))
tip_states <- tip_states[!is.na(tip_states)]

tree <- ape::read.tree(tree_file)
shared <- intersect(tree$tip.label, names(tip_states))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, shared))
tip_states <- tip_states[tree$tip.label]
if (length(tip_states) != 10) stop("Expected 10 unambiguous coxal-wall-opening proxy tips; found ", length(tip_states))

display_labels <- sub("^.*___", "", tree$tip.label)
state_factor <- factor(ifelse(tip_states, "Present", "Absent"), levels = c("Absent", "Present"))
names(state_factor) <- display_labels
tree$tip.label <- display_labels

set.seed(123)
smap <- suppressWarnings(phytools::make.simmap(tree, state_factor, model = "ER", nsim = 50))
summary_map <- summary(smap)
state_colours <- c("Absent" = "#2A6F83", "Present" = "#D96857")

draw_map <- function() {
  layout(matrix(c(1, 2), nrow = 2), heights = c(9.2, 0.8))
  par(mar = c(0.6, 0.8, 0.8, 7.8), family = "sans", fg = "#203447", col.axis = "#203447", xpd = NA, cex = 0.95)
  plot(summary_map, colors = state_colours, fsize = 0.82, ftype = "i", lwd = 2.2, split.vertical = TRUE)
  par(mar = c(0, 0, 0, 0), family = "sans", fg = "#203447")
  plot.new()
  legend("center", legend = names(state_colours), pch = 21, pt.bg = unname(state_colours), col = "#203447", pt.cex = 1.35, cex = 0.9, bty = "n", horiz = TRUE, x.intersp = 0.7)
}

repo_pdf <- file.path(repo_out_dir, "Supplementary_Fig_7_coxal_wall_opening_stochastic_map.pdf")
repo_png <- file.path(repo_out_dir, "Supplementary_Fig_7_coxal_wall_opening_stochastic_map.png")
grDevices::cairo_pdf(repo_pdf, width = 7.09, height = 4.35, family = "Arial"); draw_map(); grDevices::dev.off()
grDevices::png(repo_png, width = 7.09, height = 4.35, units = "in", res = 600, bg = "white"); draw_map(); grDevices::dev.off()

state_table <- data.frame(tree_label = names(tip_states), display_label = display_labels, coxal_wall_opening = ifelse(tip_states, "Present", "Absent"), stringsAsFactors = FALSE)
write.csv(state_table, file.path(repo_out_dir, "Supplementary_Fig_7_tip_states.csv"), row.names = FALSE)

if (!is.na(manuscript_root)) {
  canonical_base <- file.path(manuscript_root, "04_Supplementary_Figures", "Supplementary_Fig_7_coxal_wall_opening_stochastic_map_180mm")
  file.copy(repo_pdf, paste0(canonical_base, ".pdf"), overwrite = TRUE)
  file.copy(repo_png, paste0(canonical_base, ".png"), overwrite = TRUE)
}

message("Rendered Supplementary Figure 7 on the Grafen primary tree from ", length(tip_states), " unambiguous proxy tips.")
