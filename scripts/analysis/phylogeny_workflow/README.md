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
  diagnostics beside the output tree. Calibrated outputs are sensitivity trees;
  none is used as the primary working phylogeny. With only the required input
  and output paths, the defaults reproduce the documented 195 Ma calibration-
  sensitivity specification (root fixed to 195 Ma; Curculionidae constrained
  to 113-151 Ma). The historical 223 Ma ceiling is used only when supplied
  explicitly.
- `summarize_calibration_sensitivity.R` combines primary-tree and updated-
  calibration reruns into compact international CSV source-data tables.

## Calibration provenance and sensitivity

The primary working tree is
`data/phylogeny/P01_Trees/01_primary_tree_grafen.tre`. It preserves the inferred
topology, the literature-supported Caridae insertion and transparent Grafen
branch lengths. These branch lengths are a topology-based working representation
and are not interpreted as divergence times. The study does not estimate node
ages.

The exact historical 223 Ma tree is retained as
`data/phylogeny/P01_Trees/16_historical_primary_tree_223ma.tre`. Its fossil
minima (157.3 Ma for Curculionoidea and 113 Ma for Curculionidae) follow McKenna
et al. (2019; doi:10.1073/pnas.1909655116, Supplementary Table S5).
The 223 Ma ceiling follows Letsch et al. (2020;
doi:10.1111/syen.12396), who used the then-interpreted oldest polyphagan fossil
*Leehermania prorova*. *Leehermania* was subsequently reassigned to Myxophaga
(Fikáček et al. 2020; doi:10.1111/syen.12386), so 223 Ma is documented as a
historical conservative ceiling rather than direct fossil evidence for
Curculionoidea. It is therefore a labelled sensitivity representation, not the
primary tree.

The parameterized script does not claim to recreate the historical tree bit for
bit or to identify a unique timescale for the reduced proxy-tip topology. It
generates explicit calibration sensitivity trees, including a fixed 195 Ma root
scenario following the hard maximum used in Belidae dating (Li et al. 2024;
doi:10.7554/eLife.97552.3) and interval scenarios retained for robustness
testing. The sensitivity trees, constraints, optimizer diagnostics and rerun
summaries are distributed in
`data/phylogeny/P01_Trees` and
`data/supplementary_source_data/S11_Robustness/Calibration_Sensitivity`.

Results are evaluated across the topology-based primary tree and alternative
branch-length representations. Findings that change with calibration or tree
choice are reported as sensitivity-dependent.

The original source alignment, partition files and source-study inference
outputs are not redistributed in this code repository. They should be obtained
from the McKenna et al. source dataset cited in the manuscript.
