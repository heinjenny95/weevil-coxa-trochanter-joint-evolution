# Code availability

This repository contains the analysis code used for the manuscript
**Evolutionary diversification of the biological screw joint in weevils**.
The current workflows are located in `scripts/analysis`.

The release includes mesh alignment, the exact Deformetrica atlas
configuration and ordered study dataset, morphometric and allometric analyses,
downstream analysis of exported joint-geometry measurements, coxal-wall
measurements, ecological tests, and phylogenetic comparative analyses. It also
includes the script used to generate the revised Supplementary Figure 11.
Final journal-figure layout, colour adjustment, panel assembly and visual
polishing scripts are otherwise excluded; some analysis scripts still create
diagnostic plots as workflow byproducts.

The original Cinema 4D/Python implementation used for helical-axis and circle
fitting is included as
`scripts/analysis/screw_geometry_analysis/helical_path_metrics_cinema4d.py`.
It reads ordered three-dimensional semilandmark CSV files, performs the
multi-stage axis search and projected circle fit, and exports the winding-angle,
axial-span, radius and fit-error measurements used by the downstream workflows.
Running it requires a compatible Cinema 4D Python environment.

Selected compact data used to audit the revised mapping and geometry analyses
are included under `data/`; the remaining study data are distributed separately
as described in the manuscript. Local absolute paths have been replaced by
explicit placeholder roots or command-line arguments so that no
workstation-specific paths are published.
