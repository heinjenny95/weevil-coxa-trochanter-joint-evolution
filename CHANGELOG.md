# Changes

## v1.4.0 - 27 October 2026

This release makes the study-specific analysis scripts easier to install and
check for peer review. The original v1.3.2 submission release and its numerical
reference tables remain available unchanged.

- Add a separate R installer, tested package-version records, a seven-workflow
  example and automatic comparison of 25 released tables. Preserve RRPP 2.1.2
  and geomorph 4.0.10 for the released multivariate results.
- Add a Windows R example check in GitHub Actions alongside the Python tests.
- Restore the archived Cinema 4D batch-preprocessing source with a portable
  input-path override. Its syntax was checked; Cinema 4D execution was not.
- Supply source-mask boundary flags for all 68 coxae. The downstream R script
  can read this QC table to run its existing boundary-exclusion sensitivity.
  The published all-mask and lower-tail results are unchanged.
- Make the standalone PC1-PC5 export use the five-test Holm family already
  used by the primary workflow and the released table. The separate 67-PC
  table retains its own adjustment family.
- Record the historical atlas identifier correction from `t_pseudonastus` to
  `t_pseudonasutus`. No numerical PCA values change.
- Provide all 68 aligned atlas meshes, initial/estimated templates, momenta,
  control points, XML settings and logs as a release asset. PCA recomputed
  from those momenta matches all 4,624 numerical scores and variance values.
- Provide the 68 original coxa masks and the phylogenetic inference archive
  as a separate release asset, with per-file SHA-256 manifests.
- Reconstruct all 64 ordered helix inputs from the already released observed
  coordinates. Original normal vectors and export packaging are not claimed.
- Update data/code availability, citation metadata and the Nature-checklist
  documentation map. A software DOI has not been minted by this release;
  the associated preprint DOI in CITATION.cff is not a software DOI.
