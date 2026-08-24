suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
pcm_root <- if (length(args) >= 2) normalizePath(args[[2]], winslash = "/", mustWork = TRUE) else stop("Pass the main-dataset PCM output directory as argument 2.")
manuscript_root <- if (length(args) >= 3 && nzchar(args[[3]])) normalizePath(args[[3]], winslash = "/", mustWork = TRUE) else NA_character_

source(file.path(repo_root, "scripts", "analysis", "figure_rendering", "publication_style.R"))
tree <- ape::read.tree(file.path(repo_root, "data", "phylogeny", "P01_Trees", "01_primary_tree_grafen.tre"))
tip_data <- read.csv2(file.path(pcm_root, "10_Logs", "tip_level_dataset_used_for_PCM.csv"), check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
label_map <- read.csv(file.path(repo_root, "data", "supplementary_source_data", "S05_Joint_Typology", "joint_type_phylogeny_tip_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
names(label_map)[names(label_map) == "tree_tip"] <- "tree_label"

tip_data <- tip_data[stats::complete.cases(tip_data[, c("tree_label", "PC1", "PC2", "Family")]), ]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tip_data$tree_label))
tip_data <- tip_data[match(tree$tip.label, tip_data$tree_label), ]
if (!identical(tree$tip.label, tip_data$tree_label)) stop("Tree and tip-level PCA data could not be aligned.")

tip_xy <- data.frame(node = seq_along(tree$tip.label), x = tip_data$PC1, y = tip_data$PC2, tree_label = tree$tip.label, Family = tip_data$Family, stringsAsFactors = FALSE)
tip_xy$display_label <- label_map$display_label[match(tip_xy$tree_label, label_map$tree_label)]
tip_xy$display_label[is.na(tip_xy$display_label)] <- sub("^.*___", "", tip_xy$tree_label[is.na(tip_xy$display_label)])
tip_xy$nudge_x <- ifelse(tip_xy$display_label == "Curculio + 12", -0.025, 0)
tip_xy$nudge_y <- ifelse(tip_xy$display_label == "Curculio + 12", -0.014, 0)

anc_x <- phytools::fastAnc(tree, setNames(tip_xy$x, tree$tip.label))
anc_y <- phytools::fastAnc(tree, setNames(tip_xy$y, tree$tip.label))
n_tip <- length(tree$tip.label)
all_xy <- data.frame(
  node = c(seq_len(n_tip), as.integer(names(anc_x))),
  x = c(tip_xy$x, as.numeric(anc_x)),
  y = c(tip_xy$y, as.numeric(anc_y))
)
edges <- data.frame(parent = tree$edge[, 1], child = tree$edge[, 2])
edges$x <- all_xy$x[match(edges$parent, all_xy$node)]
edges$y <- all_xy$y[match(edges$parent, all_xy$node)]
edges$xend <- all_xy$x[match(edges$child, all_xy$node)]
edges$yend <- all_xy$y[match(edges$child, all_xy$node)]
internal_xy <- all_xy[all_xy$node > n_tip, ]

p <- ggplot() +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend), colour = "#AAB6BE", linewidth = 0.7, lineend = "round") +
  geom_point(data = internal_xy, aes(x = x, y = y), shape = 21, size = 2.2, fill = "white", colour = "#97A7B0", stroke = 0.6) +
  geom_point(data = tip_xy, aes(x = x, y = y, fill = Family), shape = 21, size = 3.2, colour = "white", stroke = 0.55) +
  ggrepel::geom_text_repel(data = tip_xy, aes(x = x, y = y, label = display_label, colour = Family), family = "Arial", fontface = "italic", size = 3.0, seed = 123, nudge_x = tip_xy$nudge_x, nudge_y = tip_xy$nudge_y, box.padding = 0.4, point.padding = 0.3, force = 3, force_pull = 0.15, max.time = 5, max.iter = 20000, segment.colour = "#AAB6BE", segment.size = 0.25, max.overlaps = Inf, min.segment.length = 0) +
  scale_fill_manual(values = family_cols, limits = family_order, drop = TRUE, name = "Family") +
  scale_colour_manual(values = family_cols, limits = family_order, drop = TRUE, guide = "none") +
  labs(x = "PC1", y = "PC2") +
  theme_pub(base_size = 9) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_text(face = "bold"), legend.key.width = grid::unit(3.5, "mm"), legend.spacing.x = grid::unit(1.1, "mm"), plot.margin = margin(8, 12, 7, 8)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3.2)))

repo_out_dir <- file.path(repo_root, "data", "foundational_figures")
dir.create(repo_out_dir, recursive = TRUE, showWarnings = FALSE)
repo_pdf <- file.path(repo_out_dir, "Supplementary_Fig_15_phylomorphospace_PC1_PC2.pdf")
repo_png <- file.path(repo_out_dir, "Supplementary_Fig_15_phylomorphospace_PC1_PC2.png")
ggsave(repo_pdf, p, width = 7.09, height = 4.15, units = "in", device = cairo_pdf)
ggsave(repo_png, p, width = 7.09, height = 4.15, units = "in", dpi = 600, bg = "white")

if (!is.na(manuscript_root)) {
  canonical_base <- file.path(manuscript_root, "04_Supplementary_Figures", "Supplementary_Fig_15_phylomorphospace_PC1_PC2_180mm")
  file.copy(repo_pdf, paste0(canonical_base, ".pdf"), overwrite = TRUE)
  file.copy(repo_png, paste0(canonical_base, ".png"), overwrite = TRUE)
}

message("Rendered Supplementary Figure 15 on the Grafen primary tree with ", n_tip, " proxy tips.")
