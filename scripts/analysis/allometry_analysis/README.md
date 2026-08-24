# Allometry analysis

This folder contains the specimen-level allometry workflow used to test
whether trochanter shape and screw-geometry variables scale with centroid
size.

- `allometry_full_workflow.R`
  merges analysis identifiers, PCA scores, centroid size and screw-geometry
  measurements; runs univariate analyses and the primary multivariate RRPP
  test across all non-zero atlas PC axes; and exports diagnostic tables and
  plots. The PC1-PC5 model is retained as an anatomical sensitivity analysis,
  not as the full-shape test.

The direct phylogenetically informed multivariate RRPP analysis across all
non-zero PC axes is stored with the other phylogenetic comparative analyses in
`../phylogenetic_comparative_analysis/`.
