# Screw-geometry analysis

This folder contains the scripts used to extract, process and analyse
screw-related geometry from the coxa-trochanteral joint.

The current workflow fits manually placed, ordered three-dimensional
semilandmarks with `fit_helical_paths.py`, a robust circular 3D helix model.
The original Cinema 4D/Python implementation,
`helical_path_metrics_cinema4d.py`, is retained as an auditable legacy method
and is reproduced internally by the new script. The legacy method estimates an
axis by minimizing radial circle-fit RMS but does not regress axial position on
angular position.

Recommended order:

1. `fit_helical_paths.py`
   fits the robust 3D helix, reproduces the legacy Cinema 4D grid search for
   validation, writes specimen- and point-level diagnostics, and performs the
   conditional moving-block residual bootstrap.
2. `prepare_robust_geometry_inputs.py`
   writes downstream-compatible tables for all fits, the predeclared primary
   adequacy set and the strict quality-sensitivity set.
3. `screw_geometry_extraction_workflow.R`
   reads PCA and screw-geometry measurements, applies geometry filters,
   uses fitted pitch when present and exports shape-geometry regression
   summaries. Despite its historical filename, this is a downstream
   measurement-processing script rather than the helical-fit extractor.
4. `analyze_joint_type_screw_geometry.R`
   compares screw-geometry variables among joint-type categories using
   non-phylogenetic group tests and multivariate geometry summaries.
5. `propagate_robust_geometry_uncertainty.R`
   propagates all 200 specimen-level bootstrap draws through shape, joint-type,
   PGLS and broad ecology models.
6. `plot_supplementary_figure_11.R`
   rebuilds the three-panel specimen-level relationships figure from the
   corrected plotting table.
7. `plot_robust_supplementary_pgls.R`
   rebuilds the publication-labelled primary-set PC1 plots for fitted winding
   angle and fitted axial pitch from the released proxy-tip and PGLS tables.
8. `render_publication_style_figures.R`
   rebuilds the manuscript-facing main, Extended Data and Supplementary
   figures from the checked-in robust-analysis tables while retaining the
   established manuscript palette, panel lettering and white-background
   layout. Main Figure 5 extracts the complete original descriptive
   trochanter/helix schematic as panel a, gives panels a and b equal widths and
   identical programmatic labels, and regenerates the robust morphospace as
   panel b; the separate PC1--winding-angle regression remains in
   Supplementary Figure 10.
   The renderer replaces raw diagnostic plotting templates and does not alter
   any numerical analysis results. In Extended Data Figure 5, the reference-
   phylogeny tips in panel a are intentionally black; the continuous trait
   palette is reserved for the ancestral-state maps in panels b--h.

The downstream scripts accept input and output paths as command-line
arguments. The shape-geometry workflow uses `<PCA CSV> <geometry CSV>
<output directory>`; the joint-type workflow uses `<geometry CSV> <joint-type
CSV> <output directory>`. The supplementary PGLS plotter uses
`<tip-level CSV> <PGLS-results CSV> <output directory>`.

For robust input, axial pitch is the fitted axial rise per full turn and uses
all ordered points. For legacy input only, endpoint-equivalent pitch is
calculated as `axial span * 360 / absolute winding angle`; regressions of that
legacy quotient against either component are descriptive rather than
independent tests.

The final sample flow is written to `geometry_sample_flow.csv`. The robust
analysis fitted all 64 trajectories. The primary adequacy set contains 63
trajectories, and the strict quality-sensitivity set contains 53.

The joint-type PERMANOVA uses 999 permutations with the fixed seed `20260808`
so its tabulated P value is exactly reproducible.

Rebuild the publication-style figure set from the repository root with:

```bash
Rscript scripts/analysis/screw_geometry_analysis/render_publication_style_figures.R \
  . \
  data/screw_geometry/figures
```

An optional third argument may point to the manuscript asset directory. When
provided, matching PNG, PDF and (for main Figure 5) TIFF deliverables are
written there in addition to the repository rasters.

## Robust three-dimensional helix analysis

`fit_helical_paths.py` is the current extraction workflow for the original,
ordered Cinema 4D landmark exports. The legacy measurements remain distributed
separately so that old-versus-new differences can be audited.

The script performs two fits for every specimen:

1. A direct Python port of `helical_path_metrics_cinema4d.py`. This retains the
   three-stage angular grid search and exists to verify that the raw landmark
   files reproduce the released measurements.
2. A continuous circular-helix fit. For a candidate axis, each point is
   represented by radius, unwrapped angular position and axial position. Axis
   direction, axis position, radius, axial intercept and axial rise per radian
   are optimized together. Radial residuals and deviations from the linear
   axial-position-versus-angle relationship are minimized with SciPy's
   `soft_l1` robust loss.

The fitted pitch is `2 * pi * abs(axial rise per radian)`. Unlike the released
endpoint-equivalent pitch, it uses all ordered points. The output retains both
quantities so that their difference remains explicit.

Install the Python dependencies with:

```bash
python -m pip install -r requirements-python.txt
```

Example batch run:

```bash
python fit_helical_paths.py \
  "/path/to/Landmarks Windungen" \
  "/path/to/new_validation_output" \
  --released-metrics "/path/to/winding_metrics_excelDE.csv" \
  --bootstrap 200 \
  --seed 20260810
```

The output directory contains:

- `robust_helix_metrics.csv`: specimen-level parameters, direct comparison
  with the released values, provisional quality flags and bootstrap intervals;
- `robust_helix_bootstrap_draws.csv`: every successful specimen-level
  bootstrap draw for uncertainty propagation through downstream analyses;
- `robust_helix_point_residuals.csv`: point-level fitted coordinates and radial,
  axial and combined residuals;
- `specimen_qc/*.png`: observed points, fitted helix and diagnostic panels;
- `robust_helix_summary_qc.png`: old-versus-new comparisons and fit-quality
  overview;
- `analysis_manifest.json`: software versions, settings, input checksums and
  exact uncertainty interpretation;
- `fit_failures.csv`: preserved batch errors, if any.

The moving-block residual bootstrap is conditional on the traced landmark
path. Its intervals measure fit sensitivity to the observed deviations from a
circular helix; they do **not** measure manual placement repeatability. That
would require repeated independent landmark placement. Quality flags are
diagnostics, not automatic exclusion rules. In particular, short arcs can
admit several geometrically plausible axes even when their residuals are low.

For downstream validation, the script records three analysis sets without
modifying the raw fits: `all` contains every successful fit;
`primary_adequate` requires helix RMS divided by fitted radius to be no greater
than 0.10; and `strict_good` retains only traces without provisional geometry
warnings. The primary threshold is based only on geometric model adequacy and
is fixed before inspecting downstream biological associations.

Run the synthetic recovery and parser tests with:

```bash
python -m unittest -v test_fit_helical_paths.py
```
