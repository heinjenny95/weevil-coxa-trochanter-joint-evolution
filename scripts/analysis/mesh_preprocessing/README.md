# Three-landmark GPA alignment

`gpa_alignment_cinema4d.py` is the study-specific Cinema 4D script used to
align trochanter meshes from three homologous landmarks digitized in
Checkpoint. It pairs OBJ meshes and Morphologika landmark files by basename,
performs generalized Procrustes alignment with centroid-size normalization,
and exports aligned OBJ meshes.

Run the script from the Cinema 4D Script Manager. It prompts for separate
mesh, landmark and output directories. The landmark files must contain a
`[rawpoints]` section, and the first three points are used in a fixed
homologous order.

Cinema 4D's Python environment provides the `c4d` module; this dependency is
therefore not installable through the repository's Python requirements file.

Release v1.4.0 restores `Obj_Processing.py` from an
archived pre-submission copy. It matches the documented target of 10,000
polygons, three smoothing iterations at alpha 0.40 and disabled hole closing.
Set `WEEVIL_RAW_MESH_DIR` before starting Cinema 4D, or edit `BASE_PATH` in
the script, then use the Script Manager. It writes cleaned OBJ/PLY files to
named subdirectories and a processing log. This audit parsed its syntax but
did not execute or validate the Cinema 4D API calls. The historical execution
log and exact original software environment were not recovered. Original Checkpoint/Morphologika landmark files are not bundled.
Aligned VTK inputs are in the [atlas release asset](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/download/v1.4.0/Atlas_inputs_and_outputs_68_specimens.zip);
processed surface meshes are on RADAR4KIT. Reproducing the raw mesh-preprocessing sequence
therefore requires those non-repository inputs in addition to this alignment
script.
