# Code availability

This repository contains the analysis code used for the manuscript
**Evolutionary diversification of the biological screw joint in weevils**.
The current workflows are located in `scripts/analysis`.

The release includes the three-landmark mesh-alignment script, the Deformetrica
atlas configuration and ordered study dataset, morphometric and allometric analyses,
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

The analysis-ready tabular inputs and derived outputs needed for the downstream
statistical analyses are included under `data/`. This includes the
complete source-data payload for Supplementary Tables 1-40, derived analysis
tables, analysis-ready tabular inputs, quality-control outputs and the detailed
robust geometry audit files. The imaging data and joint-surface meshes are
publicly accessible through RADAR4KIT
under https://doi.org/10.35097/9p77hjk7wa656d6k. Local absolute paths have been
replaced by explicit placeholder roots or command-line arguments so that no
workstation-specific paths are published.

The release does not contain every raw or intermediate input used in the
raw-to-derived workflow. Not redistributed here are the study-specific Cinema
4D `Obj_Processing.py` batch-remeshing script, the original Checkpoint landmark
files, aligned VTK meshes, `initial_template.vtk`, binary coxa segmentation
masks, the original ordered per-specimen semilandmark CSV files, or the source
alignment, partition files and IQ-TREE command log/report. The point-residual
table retains the observed semilandmark XYZ coordinates, but not the original
input-file packaging. Consequently, the public release supports exact
regeneration of downstream statistics, tables and scripted quantitative
figures from the released derived inputs; exact reruns of mesh preprocessing,
atlas estimation, wall-thickness extraction and source-tree inference require
the corresponding non-repository inputs.
