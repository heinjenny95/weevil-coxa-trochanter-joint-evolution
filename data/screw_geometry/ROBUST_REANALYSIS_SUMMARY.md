# Robust screw-geometry reanalysis

## Analysis sets and uncertainty scope

The ordered semilandmark trajectories were refitted with a robust circular
three-dimensional helix model. The fit estimates a continuous axis, centre,
radius, winding angle and axial slope from all trajectory points. Measurement
uncertainty was evaluated with 200 moving-block conditional-residual bootstrap
draws per specimen (12,800 successful draws in total).

- `all`: all 64 traced trajectories.
- `primary_adequate`: 63 trajectories with helix RMS / fitted radius <= 0.10.
  Only *Dryophthorus corticalis* was excluded (relative RMS = 0.180).
- `strict_good`: 53 trajectories without provisional fit-quality warnings.

The bootstrap is conditional on the traced semilandmarks. It does not estimate
repeatability of manual semilandmark placement.

## Findings that survive the audit

1. **Shape allometry is robust.** Full-atlas multivariate shape varies with
   centroid size (RRPP, 67 non-zero axes, R2 = 0.0702, p = 0.0001), with the
   strongest univariate contribution on PC2 (Holm-adjusted p = 0.000019).
2. **Axial pitch is not robustly associated with shape or size.** In the
   primary specimen-level model, pitch ~ PC1 + PC2 was only marginal
   (R2 = 0.0973, p = 0.0463) and was significant in 55.5% of measurement
   bootstrap draws. It disappeared in the strict set (R2 = 0.0021,
   p = 0.9487). Pitch was unrelated to centroid size in both sets
   (primary p = 0.611; strict p = 0.678).
3. **The winding-angle/shape association is quality-sensitive.** The primary
   specimen-level angle ~ PC1 + PC2 model was strong (R2 = 0.5119,
   p = 4.52e-10) and stable to conditional measurement uncertainty, but the
   strict-set model was not supported (R2 = 0.0591, p = 0.218). The analogous
   angle-size association survived Holm correction in the primary set
   (R2 = 0.1503, adjusted p = 0.0305) but not the strict set (R2 = 0.0522,
   adjusted p = 1.0).
4. **Phylogenetic shape-geometry results are exploratory trends, not robust
   positive tests.** With matched specimen sets aggregated to 14 primary or
   12 strict proxy tips, PC1 ~ winding angle gave p = 0.0543 in the primary
   tree and p = 0.150 in the strict set. The slope was negative across all
   primary tree variants, but none of the 13 tree-specific tests survived FDR
   correction; only 3/13 trees were nominally significant in the strict set.
   PC1 ~ axial span was supported on the primary tree (p = 0.0394) but not in
   the strict set (p = 0.137), and its significance changed under both tree
   choice and leave-one-tip-out sensitivity. PC1 ~ fitted pitch was unsupported
   on the primary and strict trees (p = 0.165 and 0.519).
5. **Joint-type comparisons are under-replicated.** The primary screw-joint
   subset contains 57 true screw-nut joints and only 3 unopposed screw
   configurations. Winding angle differs nominally in this unbalanced sample
   (Kruskal-Wallis p = 0.0202), whereas pitch and axial span do not. The strict
   set contains 51 versus 1 specimen, so inferential group tests are invalid.
   The joint-type result must therefore be described as exploratory and
   unstable rather than as a general group difference.
6. **No ecological association survives multiplicity correction.** Across
   broad host lineage, woody association, larval lifestyle and fungal
   association factors, all ecology PGLS and phylogenetic-ANOVA results were
   FDR-nonsignificant. A few nominal strict-set results were dependent on the
   response, model and uncertainty treatment and are not robust discoveries.

## Rejected or limited analyses

- Multivariate OU fits repeatedly reported non-convergence or unreliable
  Hessian solutions with only 12--14 proxy tips. They are not used for
  biological inference.
- The published kPCA sensitivity could not be rerun in this audit because the
  source `Atlas_Momentas.txt` file was not present in the accessible project
  directories. The script now prefers robust fitted pitch and accepts a zero
  angular cutoff for upstream quality-filtered tables, but the missing source
  momenta remain necessary to regenerate its numerical outputs.

## Interpretation for the manuscript

The defensible conclusion is that screw-joint geometry is heterogeneous and
that some primary-set associations suggest covariation with major trochanter
shape and size axes. Those associations are sensitive to fit-quality
definition, phylogenetic proxy sampling and tree choice. Fitted axial pitch,
joint-type contrasts and broad ecological predictors do not provide robust
confirmatory signals in the present data.
