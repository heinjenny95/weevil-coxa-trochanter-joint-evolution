# Coxa-wall analyses

This folder contains the reproducible workflow used to quantify coxal wall
thickness across the complete three-dimensional coxa masks and to analyse its
association with coxa size and the coxal wall opening. The primary thickness
measure is the median local thickness across all foreground voxels in the
largest connected component of each mask. Local thickness is defined as the
diameter of the largest sphere that contains a foreground voxel and remains
inside the coxa mask.

For computational efficiency, the complete foreground-cropped mask is analysed
at an isotropic scale of 0.75. Comparison with full-resolution estimates for a
validation subset changed the median thickness by no more than 0.5%. Distances
are converted back to original-voxel units and then to micrometres.

Recommended order:

1. `update_coxal_wall_opening_coding.R`
   Standardizes the binary coxal wall opening character used in the downstream
   analyses.
2. `update_coxa_wall_character_coding.R`
   Updates the broader coxa-wall character table used for joint-character
   summaries.
3. `extract_3d_coxa_wall_thickness.py`
   Extracts whole-volume local-thickness summaries, coxa size and mask-quality
   diagnostics from paired binary TIFF masks and OBJ meshes.
4. `analyse_3d_coxa_wall_thickness.R`
   Tests allometric scaling of whole-volume wall thickness and whether coxal
   wall opening is associated with size-corrected thickness. It also runs
   boundary-exclusion and lower-tail-thickness sensitivity analyses and creates
   the supplementary figure and result tables.

The older central-section scripts are retained only to document the exploratory
analysis that preceded the whole-volume workflow. They are not used for the
reported results.

The scripts assume the input tables distributed with the study and use
placeholder project roots rather than workstation-specific paths.
