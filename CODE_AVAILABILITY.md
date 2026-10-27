# Code and associated input availability

Study-specific R and Python analysis scripts for **Evolutionary diversification
of the biological screw joint in weevils** are released under the MIT License:

- [v1.4.0: installation, reviewer example and expanded input access](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/tag/v1.4.0)
- [v1.3.2: unchanged original submission release](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/tag/v1.3.2)

The source tree includes preprocessing/alignment code, the current robust
helix fitter and legacy implementation, statistical and morphometric scripts,
phylogenetic comparative workflows, figure renderers and derived source tables.
An archived copy of `Obj_Processing.py` has been restored; its input path was
made portable and its syntax checked. Cinema 4D execution was not repeated.

## Reviewer entry point

The root README provides installation instructions, tested operating system
and package versions, hardware context, measured installation/demo times,
expected results and commands for other input files. The R example verifies
25 reference tables from seven workflows. Test settings and results are in
`docs/tested_environment`; the Nature documentation map is
`docs/NATURE_CODE_CHECKLIST.md`. CI runs the R example and Python tests.

## Associated inputs

Release v1.4.0 supplies three checksum-indexed data assets:

1. [Atlas inputs/outputs for all 68 specimens](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/download/v1.4.0/Atlas_inputs_and_outputs_68_specimens.zip): every referenced aligned VTK mesh, initial/estimated templates, momenta, control points, configuration and logs. PCA reconstructed from the momenta matches all released numerical scores and variance values after an explicitly recorded historical identifier correction.
2. [Original coxa masks](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/download/v1.4.0/Coxa_original_masks_68_specimens.zip): 68 original uncropped binary masks and per-specimen boundary QC.
3. [Phylogenetic inference inputs/logs](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/download/v1.4.0/Phylogeny_inference_inputs_and_logs.zip): the study alignment, partition files, inferred tree and original IQ-TREE run records. The upstream phylogenomic source is attributed in the package.

The repository also provides 64 helix-input CSV files reconstructed from the
observed coordinates already retained in the released point-residual table.
They preserve those numbers and point order; they do not recreate the original
export packaging or normal vectors. The original Checkpoint landmark files
used for the earlier alignment step are not bundled. The aligned atlas inputs
make those files unnecessary for starting at the atlas stage.

Tomographic volumes and joint-surface meshes are deposited on RADAR4KIT:
https://doi.org/10.35097/9p77hjk7wa656d6k. The data repository's terms apply;
the code MIT licence does not replace third-party data, article or media terms.

## What has and has not been tested

The public validation covers the R example, atlas-output-to-PCA conversion,
all 64 helix point refits without bootstrap, the Python tests and selected
additional analysis/rendering entry points. It does not claim a new atlas
optimization, IQ-TREE inference, full-volume thickness extraction, independent
manual tracing or re-execution of all measurement-uncertainty analyses.
Historical hardware/environment uncertainties are recorded in the README.

The versioned release and commit identify the software snapshot. CITATION.cff
preserves all manuscript authors and identifies the associated preprint DOI;
that DOI is not a software DOI. A DOI-minting archive of the final code release
is recommended before publication and has not been created by this release.
