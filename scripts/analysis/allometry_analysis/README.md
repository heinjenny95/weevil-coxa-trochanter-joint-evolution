# Allometry analysis

This folder contains the specimen-level allometry workflow used to test
whether trochanter shape and screw-geometry variables scale with centroid
size.

- `allometry_full_workflow.R`
  merges analysis identifiers, PCA scores, centroid size and screw-geometry
  measurements; runs univariate analyses and the primary multivariate RRPP
  test across all non-zero atlas PC axes; and exports diagnostic tables and
  plots. The PC1-PC5 model is retained as an anatomical sensitivity analysis,
  not as the full-shape test. Its former embedded PGLS block has been removed
  from the active workflow and is preserved in version history and the dated
  project backup.

- `run_unified_phylogenetic_allometry.R`
  is the sole current univariate phylogenetic-allometry workflow. It uses the
  log of proxy-tip mean centroid size, all 15 shape tips, the 14 geometry tips
  obtained by matching the 63 main-dataset specimens before aggregation, and
  maximum-likelihood Pagel lambda constrained to the interval [0, 1]. It
  exports the complete coefficient table, the concise main-trait table and
  both analysis-ready tip datasets.

The direct phylogenetically informed multivariate RRPP analysis across all
non-zero PC axes is stored with the other phylogenetic comparative analyses in
`../phylogenetic_comparative_analysis/`.
