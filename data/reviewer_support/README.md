# Inputs supporting the reviewer example

- `coxa_mask_boundary_qc.csv` records whether each of 68 original uncropped
  binary TIFF masks touches a stack boundary. Three do. It includes dimensions,
  filenames and SHA-256 checksums. Original masks accompany v1.4.0;
  `recover_mask_boundary_qc.py` can independently recover these flags.
- `helix_traces/` contains 64 ordered CSV files reconstructed from observed
  XYZ coordinates and point names in the released
  `data/screw_geometry/robust_helix_point_residuals.csv`. Units come from
  `robust_helix_metrics.csv`. Coordinate text and point order are preserved.
  These are not original Cinema 4D exports; normal vectors are not reconstructed
  or required by the current robust fitter.

Rebuild the helix inputs into an empty output directory:

```text
python scripts/analysis/screw_geometry_analysis/reconstruct_helix_traces.py data/screw_geometry/robust_helix_point_residuals.csv data/screw_geometry/robust_helix_metrics.csv new-trace-directory
```

Original v1.3.2 numerical reference tables have not been replaced. Boundary
metadata make an existing sensitivity branch runnable. Historical specimen
spelling aliases are explicit in the preparation scripts; they add no specimens.
