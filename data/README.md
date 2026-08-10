# Revision source data

This directory contains the compact tables needed to audit the manuscript
corrections made during the final consistency review.

## `metadata`

- `specimen_key.csv`: corrected specimen metadata. Human-readable taxonomy uses
  *Rhynchites cupreus* and the phylogeny display label uses *Nedyus*.
- `taxonomic_proxy_mapping.csv`: explicit mapping of 68 specimens to the 15
  taxonomic proxy tips used by the comparative analyses. Five specimens have an
  exact sampled-genus match; the other 63 use a broader taxonomic proxy. The
  table retains the source proxy label, exact analysis-tree tip, composite
  family--tip label and any applied tree-label correction.

## `screw_geometry`

- `robust_helix_metrics.csv`: all 64 robust circular-helix fits, bootstrap
  intervals, fit diagnostics, quality flags and analysis-set membership.
- `robust_helix_bootstrap_draws.csv`: 200 successful conditional-residual
  bootstrap draws for each of the 64 trajectories.
- `robust_helix_point_residuals.csv`: point-level radial, axial and combined
  residuals for the fitted trajectory.
- `robust_geometry_all.csv`, `robust_geometry_primary_adequate.csv` and
  `robust_geometry_strict_good.csv`: downstream-compatible tables for all 64,
  the primary 63 and the strict 53 trajectories.
- `geometry_sample_flow.csv`, `shape_geometry_analysis_dataset.csv`,
  `regression_summary.csv`, `main_results.csv`,
  `joint_type_geometry_stats.csv` and
  `specimen_level_screw_joint_geometry.csv`: primary-set downstream outputs.
  Files with the `_strict_good` suffix provide the corresponding quality
  sensitivity where applicable.
- `shape_model_uncertainty_summary.csv`,
  `shape_coefficients_rubin_summary.csv`,
  `joint_type_uncertainty_summary.csv`,
  `pgls_shape_geometry_rubin_summary.csv` and
  `ecology_pgls_rubin_summary.csv`: downstream propagation of conditional
  measurement uncertainty.
- `ROBUST_REANALYSIS_SUMMARY.md`: interpretation audit across quality sets,
  bootstrap draws, phylogenetic trees, leave-one-out checks, allometry and
  ecology.
- `figure_source_data/`: compact, checked-in plotting inputs for the final
  publication-style main, Extended Data and Supplementary figures. The
  phylogenetic comparative tables use the 14-tip `pcm_matched` primary set
  reported in the manuscript.
- `figures/`: publication-style robust-analysis figures plus primary and
  strict shape-geometry sensitivity outputs.
- `sensitivity/`: matched-tip PGLS, alternative-tree and leave-one-out
  summaries plus primary/strict allometry and ecology tables.
- `winding_metrics.csv`, `winding_metrics_legacy_endpoint_axis.csv` and
  `geometry_cutoff_sensitivity.csv`: legacy endpoint-axis outputs retained for
  numerical audit, not the current primary geometry analysis.

## `morphospace`

- `kPCA_sensitivity_summary.csv` and
  `kPCA_shape_geometry_regressions.csv`: legacy threshold-aligned ordinary-PCA
  and kernel-PCA sensitivity results using the >=30-degree, n = 60 subset. The
  code now accepts robust fitted pitch and an upstream-filtered zero-degree
  cutoff, but numerical regeneration requires the source `Atlas_Momentas.txt`,
  which is not included in this compact repository.

## `phylogeny`

- `P01_Trees/`: primary and source tree variants supplied with the study. The
  comparative results remain conditional on the documented 15-tip taxonomic
  proxy mapping.

The identifiers for *Lissorhoptrus oryzophilus*, *Ormiscus saltator* and
*Tropiphorus pseudonasutus* are normalized to the atlas/PCA spellings. Current
pitch is the fitted axial rise per complete turn using all trajectory points;
the algebraic endpoint quotient is retained only in the legacy tables.
