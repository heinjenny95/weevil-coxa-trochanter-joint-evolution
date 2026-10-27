# Nature code and software documentation map

This map identifies the material supplied for the study-specific analysis
scripts in v1.4.0. It supports completion of the journal's checklist; it is
not an editorial certification or a claim that every upstream stage was rerun.

Reference: [Nature Code and Software Submission Checklist](https://www.nature.com/documents/nr-software-policy.pdf).

| Checklist item | Material supplied |
| --- | --- |
| Source code | `scripts/analysis/`, including the restored archived mesh-preprocessing source. |
| Small real or simulated dataset | Checked-in study tables used by `scripts/run_reviewer_demo.R`, three synthetic/parser tests and `data/reviewer_support/helix_traces/`. |
| Dependencies and OS with versions | Root README, `r_package_inventory.txt`, workflow-specific Python requirements and `docs/tested_environment/R_packages.csv` plus `requirements-helix-tested.txt`. Windows 11/R 4.5.2/Python 3.12.14 are the local validation environment. Legacy upstream applications use separate environments with the limits stated below. |
| Versions tested | Environment files, `validation_summary.json` and GitHub Actions logs. RRPP 2.1.2/geomorph 4.0.10 reproduce the released multivariate values. |
| Non-standard hardware | The tabular example and current helix fitter require no GPU or acquisition hardware. The README states the measured workstation and distinguishes upstream resources. Minimum RAM and every upstream resource requirement were not benchmarked. |
| Installation instructions | `Rscript --vanilla scripts/install_R_dependencies.R`; separate Python environment instructions in the README. The root Python inventory is not one combined environment. |
| Installation time | About two minutes for R dependencies and 20 seconds for helix dependencies on the stated workstation, excluding installation of R/Python. Exact measurements are in `validation_summary.json`; machine/network speed affect timing. |
| Demo instructions | `Rscript --vanilla scripts/run_reviewer_demo.R`; README also gives Python test and point-refit commands. |
| Expected output | 25 successful rows in `verification.csv`, timings, session information, logs, regenerated tables/figures and exit code zero. Python tests end with three tests and `OK`. |
| Demo runtime | About half a minute for the R example, seven seconds for cold Python tests and 74 seconds for 64 point refits on the documented workstation. |
| Other input data | Root/workflow READMEs specify argument order, columns and specimen/tip matching. New datasets need their own mappings and inclusion decisions. |
| Optional reproduction instructions | Root/workflow READMEs and checksum-indexed assets for the atlas, coxa masks and phylogenetic inference. The example compares published numbers. |
| Licence | MIT for source code; data and third-party material retain their stated terms. |
| Open repository/version | [v1.4.0](https://github.com/heinjenny95/weevil-coxa-trochanter-joint-evolution/releases/tag/v1.4.0), citation metadata and Git commit; original v1.3.2 is unchanged. |
| Functionality in manuscript | Submitted Methods describe processing and analyses; workflow READMEs/comments document executable entry points. Updating manuscript/PDF files is separate from this repository update. |

## Validation scope

Local validation verifies seven R workflows against 25 reference tables, the
68-subject atlas-output-to-PCA conversion, 64 helix point refits without
bootstrap, three Python tests and selected additional workflows/renderers.
The GitHub R check repeats installation and the example on a separate runner;
its status and logs are in the Actions tab.

No new Deformetrica optimization, IQ-TREE inference, full-volume thickness
extraction or full repeat of the 12,800 bootstrap fits is claimed. Original
Checkpoint landmarks for the earlier alignment are not bundled; aligned atlas
inputs are supplied. Cinema 4D preprocessing source was syntax checked but not
executed. Historical atlas hardware is uncertain because the XML says CPU and
the log reports `GpuMode.KERNEL`. A fully frozen historical upstream environment
is not claimed.

## Recommendations distinct from supplied materials

The checklist strongly recommends testing by a colleague unfamiliar with the
software. Automated validation is documented and is not represented as that
human test. The final code snapshot should also be preserved in a DOI-minting
archive before publication. No new software DOI is claimed; the bioRxiv DOI
in CITATION.cff identifies the associated manuscript.
