# Shared visual design for all manuscript-facing quantitative figures.
#
# The design follows the sharp, saturated appearance of the kPCA sensitivity
# figure: opaque marks, white point contours, dark typography, and a precise
# light grid.  These tokens change presentation only; they do not alter data,
# model fits, statistics, scales, or panel content.

family_order <- c(
  "Anthribidae", "Attelabidae", "Belidae", "Brentidae", "Caridae",
  "Curculionidae", "Nemonychidae"
)

family_cols <- c(
  Anthribidae = "#7E5AA6",
  Attelabidae = "#D96B5F",
  Belidae = "#7FA2D6",
  Brentidae = "#C2A532",
  Caridae = "#E58A2E",
  Curculionidae = "#55A88F",
  Nemonychidae = "#78B84B"
)

ink <- "#263746"
muted_ink <- "#52687A"
teal <- "#137F9B"
teal_light <- "#B9DDE6"
coral <- "#DF644E"
bluegrey <- "#66869A"
gold <- "#C69420"
grid_col <- "#D9E2E8"
border_col <- "#B8C5CE"

continuous_cols <- c("#164F7A", "#168AA6", "#F0CF3F", "#EF7C35", "#C43838")

cluster_cols <- c(
  "cluster 1" = "#55A88F",
  "cluster 2" = "#E7784F",
  "cluster 3" = "#718FC7",
  "cluster 4" = "#8B5AA7",
  "cluster 5" = "#C79E26",
  "cluster 6" = "#338DA0",
  "cluster 7" = "#C0608A",
  "cluster 8" = "#6E7880"
)

ecology_palettes <- list(
  host_lineage_broad = c(
    angiosperm = "#167A5C",
    gymnosperm = "#62BE9E"
  ),
  wood_association_broad = c(
    "non-wood" = "#B25C00",
    wood = "#F09A3F"
  ),
  woody_association_broad = c(
    nonwoody = "#B25C00",
    woody = "#F09A3F"
  ),
  larval_lifestyle_broad = c(
    internal = "#2D6FAF",
    other = "#79A9DB"
  ),
  fungal_association_broad = c(
    no = "#693E98",
    yes = "#A876CC"
  )
)

theme_pub <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = grid_col, linewidth = 0.45),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = border_col, linewidth = 0.5),
      axis.ticks = ggplot2::element_line(colour = border_col, linewidth = 0.45),
      axis.ticks.length = grid::unit(1.8, "mm"),
      axis.title = ggplot2::element_text(colour = ink, face = "plain"),
      axis.text = ggplot2::element_text(colour = ink),
      plot.title = ggplot2::element_text(
        colour = ink, face = "bold", size = ggplot2::rel(1.05), hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        colour = muted_ink, size = ggplot2::rel(0.86), hjust = 0
      ),
      legend.title = ggplot2::element_text(colour = ink, face = "bold"),
      legend.text = ggplot2::element_text(colour = ink),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      strip.text = ggplot2::element_text(colour = ink, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "white", colour = border_col, linewidth = 0.5
      ),
      plot.tag = ggplot2::element_text(face = "bold", colour = ink, size = ggplot2::rel(1.25)),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(8, 9, 8, 9)
    )
}

family_point <- function(size = 2.25) {
  ggplot2::geom_point(
    shape = 21,
    size = size,
    colour = "white",
    stroke = 0.45,
    alpha = 1
  )
}

scale_family_fill <- function(drop = FALSE, name = "Family") {
  ggplot2::scale_fill_manual(values = family_cols, limits = family_order, drop = drop, name = name)
}

scale_family_colour <- function(drop = FALSE, name = "Family") {
  ggplot2::scale_colour_manual(values = family_cols, limits = family_order, drop = drop, name = name)
}
