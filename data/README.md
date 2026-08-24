# Revision source data

This directory contains the complete supplementary source-data payload and
the detailed tables needed to audit the manuscript analyses and final
consistency corrections.

## `supplementary_source_data`

- Complete source-data payload for Supplementary Tables 1-40. The
  `_manifest.csv` file maps each numbered table to its CSV file(s), caption,
  byte size and SHA-256 checksum. These tables complement the detailed robust
  helix-fit outputs below.

## `metadata`

- `specimen_key.csv`: analysis identifiers and coded analytical variables.
  Human-readable taxonomy uses
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
  `robust_geometry_strict_good.csv`: downstream-compatible internal tables for
  all 64 fits, the 63-specimen main dataset and the 53-specimen
  high-confidence subset. The historical internal suffixes are mapped to the
  manuscript terminology in the released supplementary-source tables.
- `geometry_sample_flow.csv`, `shape_geometry_analysis_dataset.csv`,
  `regression_summary.csv`, `main_results.csv`,
  `joint_type_geometry_stats.csv` and
  `specimen_level_screw_joint_geometry.csv`: main-dataset downstream outputs.
  Internal files with the `_strict_good` suffix provide the corresponding
  high-confidence sensitivity where applicable.
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
  phylogenetic comparative tables use the 14-tip `pcm_matched` main dataset
  reported in the manuscript. `Figure_5_geometry_schematic.png` preserves the
  original pre-reanalysis Figure 5 composite at full resolution; the renderer
  extracts the complete descriptive trochanter/helix panel, removes embedded
  original labeling and places it beside the robust regenerated morphospace.
- `figures/`: publication-style robust-analysis figures plus primary and
  strict shape-geometry sensitivity outputs.
- `sensitivity/`: matched-tip PGLS, alternative-tree and leave-one-out
  summaries plus primary/strict allometry and ecology tables.
- `winding_metrics.csv` and `winding_metrics_legacy_endpoint_axis.csv`: legacy
  endpoint-axis outputs retained for numerical audit, not the current robust
  geometry analysis.

## `morphospace`

- `kPCA_sensitivity_summary.csv` and
  `kPCA_shape_geometry_regressions.csv`: ordination-sensitivity outputs for
  ordinary PCA and kernel PCA. Geometry-related historical outputs are kept
  separate from the current robust helix-fit evidence. The repository includes
  the kPCA workflow, ordered atlas dataset and exact Deformetrica configuration;
  regenerated momenta can be produced by rerunning that atlas workflow on the
  deposited meshes.

## `phylogeny`

- `P01_Trees/`: primary and source tree variants supplied with the study. The
  comparative results remain conditional on the documented 15-tip taxonomic
  proxy mapping. `15_calibration_sensitivity_root195_curculionidae151.tre` is a
  deterministic updated-calibration sensitivity tree; its constraint and fit-
  diagnostic CSVs are stored beside it. It is not substituted silently for the
  archived primary tree.

The identifiers for *Lissorhoptrus oryzophilus*, *Ormiscus saltator* and
*Tropiphorus pseudonasutus* are normalized to the atlas/PCA spellings. Current
pitch is the fitted axial rise per complete turn using all trajectory points;
the algebraic endpoint quotient is retained only in the legacy tables.
