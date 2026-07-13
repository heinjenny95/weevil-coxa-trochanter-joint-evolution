# Screw-geometry analysis

This folder contains the scripts used to extract and analyse screw-related
geometry from the coxa-trochanteral joint.

Recommended order:

1. `screw_geometry_extraction_workflow.R`
   reads PCA and screw-geometry measurements, applies geometry filters,
   computes derived variables such as axial pitch and exports shape-geometry
   regression summaries.
2. `analyze_joint_type_screw_geometry.R`
   compares screw-geometry variables among joint-type categories using
   non-phylogenetic group tests and multivariate geometry summaries.
