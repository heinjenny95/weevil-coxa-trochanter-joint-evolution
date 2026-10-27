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
- saved XML kernel device: CPU (runtime log reports `GpuMode.KERNEL`; historical
  hardware use is not established by the XML alone)
- deformation trajectory: 30 time points
- optimizer: gradient ascent
- initial step size: 0.01
- maximum iterations: 1,000
- final atlas sample: 68 aligned trochanter meshes

The v1.3.2 optimization configuration records the optimizer, initial step
size and maximum iteration limit, without a convergence tolerance or a
distributed execution log. The reviewer supplement subsequently located the
matching archive under `Testcode/data/atlasing`. Its log records convergence
tolerance 0.0001, final reported iteration 51, stopping after the tolerance
threshold was met, and estimation time 32 minutes 44 seconds. These recovered logs and inputs accompany v1.4.0 as a release asset;
the v1.3.2 tag remains unchanged.

The exact ordered subject list is preserved in `study_data_set.xml`. The
aligned VTK meshes referenced by that file, the Checkpoint landmarks used for
alignment and `initial_template.vtk` are not distributed in this code
source tree. The aligned meshes and template are provided in the
versioned atlas release asset below; surface meshes are also on RADAR4KIT.

The [v1.4.0 atlas asset](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/download/v1.4.0/Atlas_inputs_and_outputs_68_specimens.zip) packages
the 68 referenced aligned meshes, initial/estimated templates, momenta,
control points, XML files and recovered logs. PCA reconstructed from those
momenta matches all released numerical scores and variance values exactly.
The extractor explicitly records the historical `t_pseudonastus` to
`t_pseudonasutus` identifier correction. No atlas optimization was rerun in
this audit. The runtime log reports `GpuMode.KERNEL` despite `kernel-device`
being `cpu` in the XML; actual historical hardware use remains to be resolved.

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

After obtaining or reconstructing the aligned VTK meshes and
`initial_template.vtk`, place them in one working directory, then run:

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
