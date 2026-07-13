# Targeted phylogenetic comparative tests

This folder contains targeted follow-up analyses that were added separately
from the main phylogenetic comparative workflow.

The main phylogenetic comparative analyses for the manuscript are implemented
in:

`../atlas_pca_workflow/r/run_phylogenetic_comparative_analyses.R`

That script contains the primary analyses of phylogenetic signal, evolutionary
model fitting, PGLS, phylogenetic ANOVA, ancestral-state reconstruction and
tree-sensitivity checks.

The script in this folder is a focused add-on:

- `test_phylogenetic_allometry_axial_span.R`
  adds the phylogenetically informed allometry test for axial span after the
  main allometry workflow had already been established.

It should therefore be read as an additional targeted test, not as the complete
PCM workflow.
