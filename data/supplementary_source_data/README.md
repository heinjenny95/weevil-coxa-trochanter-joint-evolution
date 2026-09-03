# Supplementary source data

This directory contains the complete CSV payload for Supplementary Tables 1-40
of the associated manuscript. Files are grouped by analysis domain. The
`_manifest.csv` file records the table number, caption, relative path, role,
file size and SHA-256 checksum for every payload.

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
