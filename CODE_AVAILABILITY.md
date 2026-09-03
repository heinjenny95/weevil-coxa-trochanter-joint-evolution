# Code availability

This repository contains the analysis code used for the manuscript
**Evolutionary diversification of the biological screw joint in weevils**.
The current workflows are located in `scripts/analysis`.

The release includes mesh alignment, the exact Deformetrica atlas
configuration and ordered study dataset, morphometric and allometric analyses,
robust 3D helix fitting and downstream joint-geometry analysis, coxal-wall
measurements, ecological tests, and phylogenetic comparative analyses. It also
includes the publication-figure renderers, shared visual style, panel assembly
code and compact figure-source tables used for the final main, Extended Data
and Supplementary figures. Image-only descriptive panels are retained as
source assets alongside the corresponding scripted assembly workflows.

The current standard-Python circular 3D helix fitter, conditional bootstrap and
quality-set export are included in
`scripts/analysis/screw_geometry_analysis/fit_helical_paths.py` and the
associated downstream scripts. The original Cinema 4D/Python implementation
used for helical-axis and circle fitting is retained as
`scripts/analysis/screw_geometry_analysis/helical_path_metrics_cinema4d.py`.
The historical script reads ordered three-dimensional semilandmark CSV files, performs the
multi-stage axis search and projected circle fit, and exports the winding-angle,
axial-span, radius and fit-error measurements used for numerical legacy
validation. Running that historical script requires a compatible Cinema 4D
Python environment; the current robust fitter does not.

All analysis data other than the processed tomograms and coxa/trochanter meshes
are included under `data/`. This includes the
complete source-data payload for Supplementary Tables 1-40, derived analysis
tables, analysis-ready tabular inputs, quality-control outputs and the detailed
robust geometry audit files. The imaging data and joint-surface meshes are
publicly accessible through RADAR4KIT
under https://doi.org/10.35097/9p77hjk7wa656d6k. Local absolute paths have been
replaced by explicit placeholder roots or command-line arguments so that no
workstation-specific paths are published.
