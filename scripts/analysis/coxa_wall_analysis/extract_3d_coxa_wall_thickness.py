"""Extract whole-volume 3D coxal wall-thickness metrics.

Local thickness is defined for each foreground voxel as the diameter of the
largest sphere that is fully contained in the segmented coxal cuticle and
contains that voxel. The specimen-level primary metric is the median of this
whole-volume local-thickness distribution.

Inputs
------
One directory containing paired binary TIFF masks and OBJ meshes with matching
basenames. TIFF arrays must be three-dimensional and use zero for background.

Outputs
-------
``coxa_3d_wall_thickness_metrics.csv`` and an optional XLSX copy. The table
contains whole-volume thickness summaries, coxa size, and mask-QC fields.
"""

from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path

import localthickness as lt
import numpy as np
import pandas as pd
import tifffile
import trimesh
from scipy import ndimage


DEFAULT_VOXEL_SIZE_UM = 1.22
DEFAULT_PADDING_VOXELS = 2
DEFAULT_ANALYSIS_SCALE = 0.75


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Measure whole-volume 3D local thickness in binary coxa masks."
    )
    parser.add_argument(
        "root_dir",
        nargs="?",
        default=".",
        help="Directory containing paired binary TIFF masks and OBJ meshes.",
    )
    parser.add_argument(
        "--voxel-size-um",
        type=float,
        default=DEFAULT_VOXEL_SIZE_UM,
        help="Isotropic voxel edge length in micrometres.",
    )
    parser.add_argument(
        "--padding-voxels",
        type=int,
        default=DEFAULT_PADDING_VOXELS,
        help="Background padding retained around the foreground bounding box.",
    )
    parser.add_argument(
        "--analysis-scale",
        type=float,
        default=DEFAULT_ANALYSIS_SCALE,
        help=(
            "Isotropic analysis scale applied to the complete cropped 3D mask. "
            "The default 0.75 was validated against full-resolution estimates "
            "and preserves whole-volume median thickness within 0.5% in the "
            "validation subset."
        ),
    )
    parser.add_argument(
        "--prefix",
        default="coxa_3d_wall_thickness_metrics",
        help="Output filename prefix.",
    )
    parser.add_argument(
        "--specimen",
        action="append",
        default=None,
        help="Optional specimen basename to process; may be supplied repeatedly.",
    )
    parser.add_argument(
        "--no-xlsx",
        action="store_true",
        help="Write only CSV output.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from an existing output CSV and skip completed specimens.",
    )
    return parser.parse_args()


def paired_files(root_dir: Path) -> list[tuple[str, Path, Path]]:
    obj_by_stem = {path.stem: path for path in root_dir.glob("*.obj")}
    tif_by_stem = {
        path.stem: path
        for pattern in ("*.tif", "*.tiff")
        for path in root_dir.glob(pattern)
    }
    return [
        (stem, obj_by_stem[stem], tif_by_stem[stem])
        for stem in sorted(obj_by_stem.keys() & tif_by_stem.keys())
    ]


def bbox_diagonal_from_obj(obj_path: Path) -> float:
    mesh = trimesh.load(obj_path, force="mesh", process=False)
    vertices = np.asarray(mesh.vertices, dtype=float)
    if vertices.size == 0:
        raise ValueError("OBJ mesh has no vertices.")
    return float(np.linalg.norm(vertices.max(axis=0) - vertices.min(axis=0)))


def crop_to_foreground(mask: np.ndarray, padding: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    occupied = np.where(mask)
    if not occupied[0].size:
        raise ValueError("TIFF mask contains no foreground voxels.")

    start = np.array([axis.min() for axis in occupied], dtype=int)
    stop = np.array([axis.max() + 1 for axis in occupied], dtype=int)
    start = np.maximum(start - padding, 0)
    stop = np.minimum(stop + padding, np.asarray(mask.shape))
    slices = tuple(slice(int(lo), int(hi)) for lo, hi in zip(start, stop))
    return mask[slices], start, stop


def retain_largest_component(mask: np.ndarray) -> tuple[np.ndarray, int, float]:
    structure = ndimage.generate_binary_structure(rank=3, connectivity=1)
    labels, component_count = ndimage.label(mask, structure=structure)
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    largest_label = int(counts.argmax())
    foreground_count = int(mask.sum())
    retained = labels == largest_label
    removed_fraction = 1.0 - (float(retained.sum()) / foreground_count)
    return retained, int(component_count), removed_fraction


def format_triplet(values: np.ndarray | tuple[int, ...]) -> str:
    return "x".join(str(int(value)) for value in values)


def resample_binary_mask(mask: np.ndarray, scale: float) -> np.ndarray:
    if not 0 < scale <= 1:
        raise ValueError("analysis_scale must be > 0 and <= 1.")
    if scale == 1:
        return mask

    reciprocal = 1.0 / scale
    integer_step = int(round(reciprocal))
    if np.isclose(reciprocal, integer_step):
        return mask[tuple(slice(None, None, integer_step) for _ in range(mask.ndim))]

    return ndimage.zoom(mask, zoom=scale, order=0, prefilter=False).astype(bool)


def measure_specimen(
    specimen: str,
    obj_path: Path,
    tif_path: Path,
    voxel_size_um: float,
    padding_voxels: int,
    analysis_scale: float,
) -> dict[str, object]:
    raw = tifffile.imread(tif_path)
    if raw.ndim != 3:
        raise ValueError(f"Expected a 3D TIFF mask, got shape {raw.shape}.")
    original_shape = np.asarray(raw.shape, dtype=int)

    mask = raw > 0
    occupied = np.where(mask)
    if not occupied[0].size:
        raise ValueError("TIFF mask contains no foreground voxels.")
    raw_start = np.array([axis.min() for axis in occupied], dtype=int)
    raw_stop = np.array([axis.max() + 1 for axis in occupied], dtype=int)
    touches_stack_boundary = bool(
        np.any(raw_start == 0) or np.any(raw_stop == np.asarray(mask.shape))
    )

    cropped, crop_start, crop_stop = crop_to_foreground(mask, padding_voxels)
    del raw, mask, occupied
    analysis_mask = resample_binary_mask(cropped, analysis_scale)
    del cropped
    coxa_mask, component_count, removed_fraction = retain_largest_component(analysis_mask)
    del analysis_mask

    # localthickness returns local sphere radii; multiply by two for diameter.
    local_radius_vox = lt.local_thickness(coxa_mask, scale=1)
    thickness_analysis_vox = 2.0 * local_radius_vox[coxa_mask]
    thickness_vox = thickness_analysis_vox / analysis_scale
    del local_radius_vox, thickness_analysis_vox
    if not thickness_vox.size:
        raise ValueError("No local-thickness values were produced.")

    bbox_diag_units = bbox_diagonal_from_obj(obj_path)
    quantiles = np.percentile(thickness_vox, [10, 25, 50, 75, 90])

    return {
        "specimen": specimen,
        "bbox_diag_units": bbox_diag_units,
        "bbox_diag_um": bbox_diag_units * voxel_size_um,
        "mean_3d_thickness_vox": float(np.mean(thickness_vox)),
        "median_3d_thickness_vox": float(quantiles[2]),
        "p10_3d_thickness_vox": float(quantiles[0]),
        "p25_3d_thickness_vox": float(quantiles[1]),
        "p75_3d_thickness_vox": float(quantiles[3]),
        "p90_3d_thickness_vox": float(quantiles[4]),
        "max_3d_thickness_vox": float(np.max(thickness_vox)),
        "mean_3d_thickness_um": float(np.mean(thickness_vox) * voxel_size_um),
        "median_3d_thickness_um": float(quantiles[2] * voxel_size_um),
        "p10_3d_thickness_um": float(quantiles[0] * voxel_size_um),
        "p25_3d_thickness_um": float(quantiles[1] * voxel_size_um),
        "p75_3d_thickness_um": float(quantiles[3] * voxel_size_um),
        "p90_3d_thickness_um": float(quantiles[4] * voxel_size_um),
        "max_3d_thickness_um": float(np.max(thickness_vox) * voxel_size_um),
        "voxel_size_um": voxel_size_um,
        "analysis_scale": analysis_scale,
        "effective_analysis_voxel_size_um": voxel_size_um / analysis_scale,
        "foreground_voxels_retained": int(coxa_mask.sum()),
        "connected_components_6n": component_count,
        "fraction_removed_outside_largest_component": removed_fraction,
        "touches_original_stack_boundary": touches_stack_boundary,
        "original_stack_shape_zyx": format_triplet(original_shape),
        "analysis_crop_start_zyx": format_triplet(crop_start),
        "analysis_crop_stop_zyx": format_triplet(crop_stop),
        "analysis_crop_shape_zyx": format_triplet(coxa_mask.shape),
        "obj_file": obj_path.name,
        "tif_file": tif_path.name,
        "thickness_definition": "diameter_of_largest_inscribed_sphere_containing_voxel",
        "analysis_domain": "complete_3d_coxa_mask_after_foreground_crop",
    }


def main() -> None:
    args = parse_args()
    root_dir = Path(args.root_dir).expanduser().resolve()
    if not root_dir.is_dir():
        raise NotADirectoryError(root_dir)
    if not 0 < args.analysis_scale <= 1:
        raise ValueError("--analysis-scale must be > 0 and <= 1.")

    pairs = paired_files(root_dir)
    if args.specimen:
        requested = set(args.specimen)
        pairs = [pair for pair in pairs if pair[0] in requested]
        missing = requested - {pair[0] for pair in pairs}
        if missing:
            raise ValueError(f"Requested specimens not found as paired files: {sorted(missing)}")
    if not pairs:
        raise RuntimeError("No paired OBJ and TIFF files with matching basenames found.")

    csv_path = root_dir / f"{args.prefix}.csv"
    rows: list[dict[str, object]] = []
    completed: set[str] = set()
    if args.resume and csv_path.exists():
        previous = pd.read_csv(csv_path)
        rows = previous.to_dict(orient="records")
        completed = set(previous["specimen"].astype(str))
        pairs = [pair for pair in pairs if pair[0] not in completed]
        print(f"Resuming after {len(completed)} completed specimens.", flush=True)
    errors: list[dict[str, str]] = []
    started = time.time()
    print(f"Processing {len(pairs)} paired specimens from {root_dir}")
    for index, (specimen, obj_path, tif_path) in enumerate(pairs, start=1):
        specimen_started = time.time()
        try:
            row = measure_specimen(
                specimen,
                obj_path,
                tif_path,
                args.voxel_size_um,
                args.padding_voxels,
                args.analysis_scale,
            )
            row["runtime_seconds"] = time.time() - specimen_started
            rows.append(row)
            print(
                f"[{index:02d}/{len(pairs):02d}] {specimen}: "
                f"median={row['median_3d_thickness_um']:.3f} um "
                f"({row['runtime_seconds']:.1f}s)",
                flush=True,
            )
            pd.DataFrame(rows).sort_values("specimen").to_csv(csv_path, index=False)
        except Exception as exc:  # continue to produce an explicit error manifest
            errors.append({"specimen": specimen, "error": str(exc)})
            print(f"[{index:02d}/{len(pairs):02d}] ERROR {specimen}: {exc}", flush=True)
        gc.collect()

    if not rows:
        raise RuntimeError("Every specimen failed; no output was written.")

    output = pd.DataFrame(rows).sort_values("specimen")
    output.to_csv(csv_path, index=False)
    if not args.no_xlsx:
        output.to_excel(root_dir / f"{args.prefix}.xlsx", index=False)
    if errors:
        pd.DataFrame(errors).to_csv(root_dir / f"{args.prefix}_errors.csv", index=False)

    print(f"Wrote {csv_path}")
    print(f"Completed {len(rows)}/{len(pairs)} specimens in {time.time() - started:.1f}s")


if __name__ == "__main__":
    main()

