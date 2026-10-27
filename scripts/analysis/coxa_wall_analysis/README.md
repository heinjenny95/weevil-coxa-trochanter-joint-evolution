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
five-specimen validation subset changed the median thickness by 0.46-5.21%.
Distances
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

The extraction script requires the original paired binary coxa masks and OBJ
meshes. The [68 original masks](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/download/v1.4.0/Coxa_original_masks_68_specimens.zip) are attached to
release v1.4.0; derived thickness tables remain in the source tree.
Re-extraction uses those masks and the corresponding meshes on RADAR4KIT. Scripts use placeholder
project roots rather than workstation-specific paths.

## Recovered boundary QC in the reviewer supplement

The v1.3.2 metrics table omitted `touches_original_stack_boundary`, although
the downstream R script required it for boundary-exclusion sensitivity. The
supplement supplies `data/reviewer_support/coxa_mask_boundary_qc.csv`, recovered
from the six faces of each of the 68 original TIFF masks. Stack dimensions
were checked against the released metrics and SHA-256 checksums identify the
source masks. Three masks touch a stack boundary.

Pass `--boundary-qc data/reviewer_support/coxa_mask_boundary_qc.csv` to
`analyse_3d_coxa_wall_thickness.R` with its other required arguments. The script
matches unique specimen IDs and rejects missing or invalid flags. It retains
an existing boundary column when one is already supplied. Main and lower-tail
statistics on all masks reproduce the released values; additional
boundary-exclusion output uses 65 masks and is not present in the v1.3.2
reference tables. No published numerical result was overwritten.

To independently recover the flags from the original masks, use:

```text
python recover_mask_boundary_qc.py metrics.csv ORIGINAL_MASK_DIRECTORY recovered_qc.csv
```

This optional recovery requires numpy, tifffile and an appropriate TIFF codec
(imagecodecs was used in the validation). It is not needed to run the bundled
R reviewer example.
