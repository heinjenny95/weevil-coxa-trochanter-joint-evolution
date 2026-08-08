# Revision source data

This directory contains the compact tables needed to audit the manuscript
corrections made during the final consistency review.

## `metadata`

- `specimen_key.csv`: corrected specimen metadata. Human-readable taxonomy uses
  *Rhynchites cupreus* and the phylogeny display label uses *Nedyus*.
- `taxonomic_proxy_mapping.csv`: explicit mapping of 68 specimens to the 15
  taxonomic proxy tips used by the comparative analyses. Five specimens have an
  exact sampled-genus match; the other 63 use a broader taxonomic proxy.

## `screw_geometry`

- `winding_metrics.csv`: 64 exported helical trajectories before the 30-degree
  winding-angle eligibility cutoff.
- `geometry_sample_flow.csv`: explicit accounting of the four excluded and 60
  eligible trajectories.
- `shape_geometry_analysis_dataset.csv`: the corrected 60-specimen dataset used
  for shape--geometry regressions.
- `regression_summary.csv`: revised OLS results, including the three simple
  geometry-on-geometry regressions.
- `joint_type_geometry_stats.csv`: joint-type comparison results for the 58
  specimens in recurrent screw-bearing joint categories.

The identifier for *Lissorhoptrus oryzophilus* is normalized to
`308_lisshorhoptrus_oryzophilus_trochanter_aligned`, matching the atlas/PCA
dataset. Pitch is algebraically derived from winding angle and axial span and
must not be interpreted as an independent measured trait.
