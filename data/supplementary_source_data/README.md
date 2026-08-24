# Supplementary source data

This directory contains the complete CSV payload for Supplementary Tables 1-40
of the associated manuscript. Files are grouped by analysis domain. The
`_manifest.csv` file records the table number, caption, relative path, role,
file size and SHA-256 checksum for every payload.

Robust-geometry tables distinguish the 63-specimen main dataset
from the 53-specimen high-confidence subset. Files that combine both definitions include
an explicit `quality_set` column; phylogenetic geometry tables likewise
distinguish the 14-tip main and 12-tip high-confidence matched datasets. Atlas-only
shape analyses retain their documented 15-tip scope.

`S11_Robustness/Calibration_Sensitivity/` documents the provenance of the
historical 223 Ma calibration ceiling and contains primary-versus-195 Ma reruns
for phylogenetic signal, evolutionary models, predictor slopes, continuous
ancestral states and ecological tests. These files supplement the numbered
robustness tables without creating an additional Supplementary Table number.
`calibration_sensitivity_manifest.csv` records the repository-relative path,
size and SHA-256 checksum of every added calibration tree, diagnostic and rerun
table.

Processed tomograms, coxa and trochanter joint-surface meshes and corresponding
specimen metadata will be made publicly accessible through RADAR4KIT upon
publication under https://doi.org/10.35097/9p77hjk7wa656d6k and are not
duplicated here as bulk imaging or mesh files.
