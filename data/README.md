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

- `winding_metrics.csv`: 64 exported traced surface trajectories before the 30-degree
  winding-angle eligibility cutoff.
- `geometry_sample_flow.csv`: explicit accounting of the four excluded and 60
  eligible trajectories.
- `geometry_cutoff_sensitivity.csv`: model summaries at 0, 10, 20, 30 and 45
  degree minimum winding-angle thresholds (n = 64 down to n = 58).
- `shape_geometry_analysis_dataset.csv`: the corrected 60-specimen dataset used
  for shape--geometry regressions.
- `regression_summary.csv`: revised OLS results, including the three simple
  geometry-on-geometry regressions.
- `joint_type_geometry_stats.csv`: joint-type comparison results for the 58
  specimens in recurrent screw-bearing joint categories.
- `main_results.csv`: corrected model-level and PC1/PC2 effect statistics for
  the 60-specimen analysis and 56-specimen main-region sensitivity subset.
- `specimen_level_screw_joint_geometry.csv`: all 60 eligible trajectories,
  including the two specimens outside the recurrent joint-type comparison.

## `morphospace`

- `kPCA_sensitivity_summary.csv` and
  `kPCA_shape_geometry_regressions.csv`: threshold-aligned ordinary-PCA and
  kernel-PCA sensitivity results; all geometry models use the >=30-degree,
  n = 60 subset.

## `phylogeny`

- `P01_Trees/`: primary and source tree variants supplied with the study. The
  comparative results remain conditional on the documented 15-tip taxonomic
  proxy mapping.

The identifier for *Lissorhoptrus oryzophilus* is normalized to
`308_lisshorhoptrus_oryzophilus_trochanter_aligned`, matching the atlas/PCA
dataset. Pitch is algebraically derived from winding angle and axial span and
must not be interpreted as an independent measured trait.
