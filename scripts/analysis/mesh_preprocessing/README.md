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
