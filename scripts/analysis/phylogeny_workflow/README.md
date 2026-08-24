# Phylogeny workflow scripts

This folder contains study-specific scripts used to prepare phylogenetic inputs
for the comparative analyses.

- `subset_phylip_alignment.py` subsets the source sequential PHYLIP alignment to
  the taxa retained for tree inference.
- `calibrate_curculionoidea_tree.R` roots the inferred tree, inserts Caridae at
  the manuscript backbone position, applies Grafen branch lengths, and performs
  explicitly parameterized time scaling with `ape::chronos`. Calibration bounds,
  start ages, model, smoothing parameter and maximum dual iterations are command
  line arguments, and the script writes both the node constraints and optimizer
  diagnostics beside the output tree. With only the required input and output
  paths, the defaults reproduce the documented 195 Ma calibration-sensitivity
  specification (root fixed to 195 Ma; Curculionidae constrained to 113-151 Ma).
  The historical 223 Ma ceiling is used only when supplied explicitly.
- `summarize_calibration_sensitivity.R` combines primary-tree and updated-
  calibration reruns into compact international CSV source-data tables.

## Calibration provenance and sensitivity

The exact downstream primary tree is archived as
`data/phylogeny/P01_Trees/01_primary_tree_calibrated_grafen.tre`. Its fossil
minima (157.3 Ma for Curculionoidea and 113 Ma for Curculionidae) follow
McKenna et al. (2019; doi:10.1073/pnas.1909655116, Supplementary Table S5).
The 223 Ma ceiling follows Letsch et al. (2020;
doi:10.1111/syen.12396), who used the then-interpreted oldest polyphagan fossil
*Leehermania prorova*. *Leehermania* was subsequently reassigned to Myxophaga
(Fikáček et al. 2020; doi:10.1111/syen.12386), so 223 Ma is documented as a
historical conservative ceiling rather than direct fossil evidence for
Curculionoidea.

The parameterized script does not claim to recreate the earlier primary tree
bit-for-bit. That archived tree remains the exact input for the reported
primary analyses. Instead, the script regenerates an explicit deterministic
calibration sensitivity tree using a 195 Ma root age, the current hard maximum
used in Belidae dating (Li et al. 2024; doi:10.7554/eLife.97552.3), and the
tighter 151 Ma Curculionidae ceiling used by Letsch et al. The sensitivity tree,
constraints, optimizer diagnostics and rerun summaries are distributed in
`data/phylogeny/P01_Trees` and
`data/supplementary_source_data/S11_Robustness/Calibration_Sensitivity`.

The 195 Ma sensitivity reruns produced no raw or adjusted 0.05-threshold changes
for phylogenetic signal, predictor slopes, ecological PGLS or phylogenetic ANOVA.
The best univariate model for PC1 changed from early burst to Brownian motion;
all other univariate rankings were unchanged. Numerical ancestral-state root
estimates shifted, as expected for a branch-length sensitivity.

The original source alignment, partition files and source-study inference
outputs are not redistributed in this code repository. They should be obtained
from the McKenna et al. source dataset cited in the manuscript.
