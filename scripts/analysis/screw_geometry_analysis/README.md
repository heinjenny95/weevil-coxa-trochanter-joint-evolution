# Screw-geometry analysis

This folder contains the scripts used to process and analyse exported
screw-related geometry from the coxa-trochanteral joint.

The helical-axis and circle fits were originally produced from manually placed
semilandmarks with a custom Cinema 4D/Python workflow. That source file and the
raw semilandmark coordinates were not retained in a portable form. The scripts
here therefore begin with the exported measurement table; they do not recreate
the measurements from semilandmarks.

Recommended order:

1. `screw_geometry_extraction_workflow.R`
   reads PCA and screw-geometry measurements, applies geometry filters,
   computes derived variables such as axial pitch and exports shape-geometry
   regression summaries. Despite its historical filename, this is a downstream
   measurement-processing script rather than the original helical-fit extractor.
2. `analyze_joint_type_screw_geometry.R`
   compares screw-geometry variables among joint-type categories using
   non-phylogenetic group tests and multivariate geometry summaries.
3. `plot_supplementary_figure_11.R`
   rebuilds the three-panel specimen-level relationships figure from the
   corrected plotting table and labels pitch explicitly as a derived quantity.

Both scripts accept input and output paths as their first command-line
arguments. The shape-geometry workflow uses `<PCA CSV> <geometry CSV>
<output directory>`; the joint-type workflow uses `<geometry CSV> <joint-type
CSV> <output directory>`.

Axial pitch is calculated as `axial span * 360 / absolute winding angle`.
Regressions of pitch against either component are exported for descriptive
completeness only and are not independent tests of biological trait coupling.

The final sample flow is written to `geometry_sample_flow.csv`. In the current
dataset, 64 trajectories were exported, four fell below the documented 30-degree
cutoff and 60 were matched to the atlas after correcting the legacy
Lissorhoptrus identifier.

The joint-type PERMANOVA uses 999 permutations with the fixed seed `20260808`
so its tabulated P value is exactly reproducible.
