# Weevil coxa-trochanteral joint evolution

Analysis code accompanying the manuscript **Evolutionary diversification of
the biological screw joint in weevils**.

The repository contains the preprocessing, Deformetrica atlas, statistical
and morphometric analysis workflows used to study trochanter shape,
screw-joint geometry, joint typology, phylogenetic history, allometry and
broad ecological associations across Curculionoidea. Final journal-layout,
colour-adjustment and panel-assembly scripts are intentionally excluded.

## Repository contents

- `scripts/analysis/atlas_pca_workflow`: PCA extraction, morphospace
  clustering, family disparity, allometry and the main phylogenetic
  comparative workflow.
- `scripts/analysis/deformetrica_workflow`: exact deterministic-atlas
  configuration, ordered 68-subject dataset XML, portable runner and an
  output-free Jupyter workflow.
- `scripts/analysis/mesh_preprocessing`: three-landmark GPA alignment used to
  prepare trochanter meshes for atlas construction.
- `scripts/analysis/coxa_wall_analysis`: coxal size, wall-thickness and opening
  analyses.
- `scripts/analysis/ecology_analysis`: exploratory ecological association
  tests.
- `scripts/analysis/phylogeny_workflow`: alignment subsetting and preparation
  of the study phylogeny.
- `scripts/analysis/phylogenetic_comparative_analysis`: targeted add-on
  phylogenetically informed tests not covered by the main workflow.
- `scripts/analysis/standalone_workflows`: joint typology, screw geometry and
  combined analysis workflows.
- `scripts/code_manifest.csv`: file sizes and SHA-256 checksums for the
  released scripts.
- `r_package_inventory.txt` and `python_requirements.txt`: detected software
  dependencies.

## Running the analyses

The scripts are a transparent analysis record rather than a containerized
one-command pipeline. Most workflows expect the input tables, trees, meshes or
segmentation-derived measurements distributed with the study.

1. Clone this repository.
2. Install the R and Python dependencies listed at the repository root.
3. Replace placeholder roots such as `<MANUSCRIPT_PROJECT_ROOT>`,
   `<KERNEL_WIDTH_PROJECT_ROOT>` and `<BEETLE_JOINTS_ROOT>` with local paths, or
   pass the command-line inputs documented by the individual script.
4. Run the workflow-specific scripts from `scripts/analysis`. The
   Deformetrica workflow has its own detailed README and notebook entry point.

The intended generic project layout is:

```text
<MANUSCRIPT_PROJECT_ROOT>/
  analysis_data/
    Input/
    Results/
```

Some scripts create diagnostic graphics as analysis byproducts. These are not
the final publication layouts.

## Data availability

Tomographic data, processed meshes, morphometric scores, geometry
measurements, phylogenetic trees and source-data tables are distributed
separately as described in the manuscript's Data availability statement. Raw
and processed research data are not duplicated in this code repository.

## Reproducibility scope

The code was developed in R 4.x and Python 3.x; atlas construction used
Deformetrica 4.3.0. Syntax and repository hygiene were checked for the v1.0.0
release. Exact numerical reproduction additionally depends on the study data,
external software described in the Methods and the package versions used for
the original analyses.

## Citation

Please cite the associated manuscript and this repository. Citation metadata
are provided in `CITATION.cff`; the manuscript DOI can be added after
publication without changing the analysis history.

## License

The code is released under the MIT License. See `LICENSE`.
