# Deformetrica atlas workflow

This directory records the study-specific workflow used to construct the
trochanter atlas in Deformetrica 4.3.0 and to derive the momenta used for
principal-component analysis.

## Provenance

The workflow was developed from the public *Landmark-Free Morphometry*
tutorial by Toussaint and colleagues
([GitLab repository](https://gitlab.com/ntoussaint/landmark-free-morphometry),
commit `cc10f96`) and adapted to the weevil dataset. The checked upstream
revision did not contain an explicit licence file, so its notebooks and mouse
demonstration code are not copied here. Instead, this directory contains an
original, study-specific implementation, the exact atlas configuration and an
attribution to the upstream workflow.

## Final study configuration

- Deformetrica version: 4.3.0
- model: deterministic atlas in three dimensions
- object type: `SurfaceMesh`
- attachment: current
- template noise standard deviation: 0.05
- template and deformation kernel: KeOps, width 0.1
- kernel device used for the archived run: CPU
- deformation trajectory: 30 time points
- optimizer: gradient ascent
- initial step size: 0.01
- maximum iterations: 1,000
- final atlas sample: 68 aligned trochanter meshes

The exact ordered subject list is preserved in `study_data_set.xml`. Meshes
are distributed separately with the study data and are not duplicated in
this code repository.

## Files

- `deformetrica_atlas_workflow.ipynb`: output-free notebook entry point.
- `build_dataset_xml.py`: builds a deterministic-atlas dataset file from a
  directory of aligned VTK meshes.
- `run_deformetrica_atlas.py`: runs the atlas estimation, records the log and
  gives the principal outputs stable filenames.
- `config/model.xml`: exact final model parameters.
- `config/optimization_parameters.xml`: exact final optimizer parameters.
- `study_data_set.xml`: exact subject order used in the 68-specimen run.

The resulting `Atlas_Momentas.txt`, `Atlas_ControlPoints.txt` and
`data_set.xml` are passed to
`../atlas_pca_workflow/r/extract_pca_scores_with_specimen_ids.R`, which
performs the unscaled, centred PCA of subject momenta used in the manuscript.

## Command-line use

Place the aligned VTK meshes and `initial_template.vtk` in one working
directory, then run:

```bash
python build_dataset_xml.py /path/to/atlasing --output data_set.xml
cp config/model.xml config/optimization_parameters.xml /path/to/atlasing/
python run_deformetrica_atlas.py /path/to/atlasing
Rscript ../atlas_pca_workflow/r/extract_pca_scores_with_specimen_ids.R \
  /path/to/atlasing/Atlas_Momentas.txt \
  /path/to/atlasing/data_set.xml \
  /path/to/atlasing/pca_outputs
```

Use `study_data_set.xml` instead of rebuilding it when reproducing the exact
archived subject order.
