# Atlas ordination workflow

This directory contains the scripts used to summarize the subject-specific Deformetrica momenta and assess the robustness of the resulting morphospace.

## Primary analysis

The `r/` scripts reproduce the ordinary-PCA workflow used for the main analyses, including PCA variance, specimen scores, morphospace clustering, disparity and allometry.

## Nonlinear sensitivity analysis

`python/run_kpca_sensitivity.py` compares ordinary PCA with two radial-basis-function kernel PCA ordinations:

1. `gamma = 0.25`, retained from the adapted exploratory notebook.
2. A data-adaptive kernel with `gamma = 1 / (2 * median squared pairwise distance)`.

The script matches and sign-aligns the first five kernel axes to PC1-PC5 and compares global morphospace geometry, family and joint-type structure, allometry, and associations with screw joint geometry. When a robust fitted-pitch column is present, it is used in preference to the legacy endpoint quotient. Current shape-geometry analyses use the upstream quality-filtered robust-helix tables without an angular cutoff; historical endpoint-derived outputs are not used for current inference. The script writes tables, publication-ready figures and an input-file checksum manifest to the selected output directory.

Run from the repository root with study-specific input paths:

```powershell
python scripts/analysis/atlas_pca_workflow/python/run_kpca_sensitivity.py `
  --momenta path/to/Atlas_Momentas.txt `
  --linear-scores path/to/PCA_scores_with_specimen_id.csv `
  --analysis-data path/to/PCA_scores_with_specimen_id_with_centroid_size.csv `
  --joint-data path/to/specimen_joint_types.csv `
  --output-dir path/to/kpca_sensitivity_output
```

Kernel eigenvalue fractions describe inertia in kernel feature space and should not be interpreted as percentages of the original anatomical shape variance.

## Atlas kernel-width selection audit

`python/analyze_kernel_width_stability.py` quantifies agreement among the five
archived Deformetrica runs used during parameter selection (kernel widths
0.025, 0.075, 0.1, 0.15 and 0.2). These development-stage runs contain the 67
specimens available at that time. They preceded the final 68-specimen atlas,
which added *Rhynchites cupreus* while retaining kernel width 0.1 without
further tuning. The archived runs therefore document parameter selection; they
are not final-dataset sensitivity reruns and do not imply that 0.1 is uniquely
optimal.

The audit reports rank concordance of all pairwise specimen distances,
PC1-PC5 distance concordance, Procrustes similarity and disparity, nearest-
neighbour overlap, repeated 80% subsampling intervals and 999-permutation
Mantel P values. Effect-size concordance, rather than Mantel significance, is
used to assess stability.

```powershell
python scripts/analysis/atlas_pca_workflow/python/analyze_kernel_width_stability.py `
  --kernel-root path/to/archived_kernel_width_runs `
  --final-scores data/supplementary_source_data/S01_PCA_and_Morphospace/PCA_scores_with_specimen_id.csv `
  --output-dir data/supplementary_source_data/S01_PCA_and_Morphospace/Kernel_Width_Selection
```
