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

This folder documents and releases the GPA alignment step only. The
study-specific `Obj_Processing.py` script used earlier for batch remeshing,
smoothing and mesh cleanup is not included in the public repository. The
original Checkpoint/Morphologika landmark files and the processed or aligned
mesh inputs are likewise not included here; surface meshes are deposited
separately on RADAR4KIT. Reproducing the raw mesh-preprocessing sequence
therefore requires those non-repository inputs in addition to this alignment
script.
