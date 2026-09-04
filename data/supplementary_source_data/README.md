# Supplementary source data

This directory contains the complete CSV payload for Supplementary Tables 1-40
of the associated manuscript. Files are grouped by analysis domain. The
`_manifest.csv` records each unique payload once, including its primary table
group, caption, relative path, role, file size and SHA-256 checksum.
`_table_file_map.csv` is the normalized authoritative table-to-file map, with
one row per table-file relationship. It therefore records the three shared
source files for Tables 10 and 11 as separate rows instead of a combined
`10;11` value. A human-readable workbook containing Tables 1-40 is attached to
the corresponding GitHub release.

Table 1 includes `radar_filename`, the exact deposited basename used to bridge
the 12 genus-only `_sp` files and the legacy `Rynchites_cupreus` filename to
the normalized analytical identifiers. The captions for Tables 24 and 33 state
explicitly that the continuous ancestral-state summaries are conditional on
Brownian motion.

Robust-geometry inference uses the 63-specimen main dataset without an angular
cutoff. Other fit warnings remain descriptive diagnostics and do not define an
exclusion subset because point count and angular coverage partly reflect the
biological extent of the helix. Phylogenetic geometry tables use the 14 matched
proxy tips from this main dataset. Atlas-only shape analyses retain their
documented 15-tip scope.

`S01_PCA_and_Morphospace/Kernel_Width_Selection/` contains the quantitative
audit supporting Supplementary Fig. 25. The five archived kernel-width runs
contain 67 specimens and were performed during workflow development. Kernel
width 0.1 was then retained without further tuning for the final 68-specimen
atlas, which additionally included *Rhynchites cupreus*. The audit is therefore
evidence for parameter selection, not a sensitivity rerun of the final atlas.

`S11_Robustness/Calibration_Sensitivity/` documents why the historical 223 Ma
tree was replaced by a topology-based Grafen primary working tree and contains
comparisons across calibrated and alternative branch-length representations
for phylogenetic signal, evolutionary models, predictor slopes, continuous
ancestral states and ecological tests. These files supplement the numbered
robustness tables without creating an additional Supplementary Table number.
`calibration_sensitivity_manifest.csv` records the repository-relative path,
size and SHA-256 checksum of every added calibration tree, diagnostic and rerun
table.

Processed tomograms and coxa and trochanter joint-surface meshes are publicly
accessible through RADAR4KIT under
https://doi.org/10.35097/9p77hjk7wa656d6k and are not duplicated here as bulk
imaging or mesh files.
