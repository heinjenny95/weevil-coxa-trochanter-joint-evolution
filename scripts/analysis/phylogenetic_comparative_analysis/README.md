# Phylogenetic comparative analyses

This folder contains the phylogenetic comparative analyses for the manuscript.

The main workflow is implemented in:

`run_phylogenetic_comparative_analyses.R`

That script contains the primary analyses of phylogenetic signal, evolutionary
model fitting, PGLS, phylogenetic ANOVA, ancestral-state reconstruction and
tree-sensitivity checks.

Additional focused scripts in this folder are:

- `fit_multivariate_allometry_score_pgls.R`
  fits the phylogenetically informed visualization model for the multivariate
  PC1-PC5 allometry score.
- `test_phylogenetic_allometry_axial_span.R`
  adds the phylogenetically informed allometry test for axial span after the
  main allometry workflow had already been established.
