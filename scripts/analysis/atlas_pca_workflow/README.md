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
