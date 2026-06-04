# Weevil coxa-trochanteral joint evolution

Analysis code for the manuscript:

**Evolutionary diversification of the coxa-trochanteral joint morphology in weevils (Coleoptera: Curculionoidea)**

This repository is prepared as a public code companion for the manuscript. It contains analysis scripts and minimal repository metadata. Final journal-figure layout and polishing scripts are intentionally omitted.

## Repository Status

This code-only export was created from the current local manuscript project on 2026-06-04. The original working supplement folder was left untouched as a backup.

Important caveat: scripts that originally depended on local absolute paths have been sanitized with placeholder roots such as `<MANUSCRIPT_PROJECT_ROOT>`, `<KERNEL_WIDTH_PROJECT_ROOT>` and `<BEETLE_JOINTS_ROOT>`. Before rerunning the full workflow on another computer, replace these placeholders with local project paths or parameterize the scripts.

## Contents

- `scripts/analysis`: current analysis scripts, grouped by task.

## Manifests

- `scripts/code_manifest.csv`: code-only manifest.
- `r_package_inventory.txt`: automatically detected R package calls in exported R scripts.

## Quick Start

1. Browse the analysis scripts in `scripts/analysis`.
2. Start with the workflow-specific subfolders for atlas PCA, coxa-wall analysis, phylogenetic comparative analyses and standalone analysis workflows.
3. Replace placeholder roots such as `<MANUSCRIPT_PROJECT_ROOT>`, `<KERNEL_WIDTH_PROJECT_ROOT>` and `<BEETLE_JOINTS_ROOT>` with local paths before rerunning scripts.

## Reproducibility Notes

The repository is designed as a transparent code record for review and publication. It is not a fully containerized rerun environment. Several scripts were written during exploratory analysis and therefore assume a project-specific folder layout. Some analysis pipelines create diagnostic plots as byproducts, but final figure layout, colour-adjustment, panel-assembly and journal-formatting scripts are not included.

## Citation

If you use this code, please cite the associated manuscript and the archived GitHub or Zenodo release. Before public release, update `CITATION.cff` with the final author list, DOI and repository URL.

## Before Public Release

Before making the repository public, check `RELEASE_CHECKLIST.md`. The main remaining manual edits are author names, final manuscript citation, Zenodo DOI and any desired path-parameterization of scripts.

## License

Code is released under the MIT License, unless otherwise noted. See `LICENSE`.
