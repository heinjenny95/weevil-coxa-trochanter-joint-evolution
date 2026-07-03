# Phylogeny workflow scripts

This folder contains study-specific scripts used to prepare phylogenetic inputs
for the comparative analyses.

- `subset_phylip_alignment.py` subsets the source sequential PHYLIP alignment to
  the taxa retained for tree inference.
- `calibrate_curculionoidea_tree.R` roots the inferred tree, inserts Caridae at
  the manuscript backbone position, applies Grafen branch lengths, and performs
  fossil-calibrated time scaling with `ape::chronos`.

The original source alignment, partition files and source-study inference
outputs are not redistributed in this code repository. They should be obtained
from the McKenna et al. source dataset cited in the manuscript.
