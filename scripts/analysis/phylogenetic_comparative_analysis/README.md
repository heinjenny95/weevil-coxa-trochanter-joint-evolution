# Phylogenetic comparative analyses

This folder contains the phylogenetic comparative analyses for the manuscript.

The main workflow is implemented in:

`run_phylogenetic_comparative_analyses.R`

That script contains the primary analyses of phylogenetic signal, evolutionary
model fitting, PGLS, phylogenetic ANOVA, ancestral-state reconstruction and
tree-sensitivity checks.

Additional focused scripts in this folder are:

- `run_phylogenetic_multivariate_allometry.R`
  performs the direct phylogenetically informed multivariate RRPP analysis
  across all non-zero atlas PC axes. The fitted one-dimensional projection is
  exported only for visualization.
- `test_phylogenetic_allometry_axial_span.R`
  adds the phylogenetically informed allometry test for axial span after the
  main allometry workflow had already been established.
