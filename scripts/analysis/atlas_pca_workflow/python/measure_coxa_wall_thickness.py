# ============================================================
# Coxa metrics (combined):
# 1) Coxa size from OBJ (bounding-box diagonal)
# 2) Cuticle thickness from binary TIF label stacks:
#    - find FIRST and LAST slice that contain cuticle
#    - take mid-slice between them
#    - compute thickness on that 2D slice via EDT:
#         thickness_vox ≈ 2 * distance_transform_edt(cuticle_mask)
#    - summarize by MEDIAN (+ P10/P90)
#
# Outputs:
#  - Excel file (.xlsx): safest for Excel
#  - German CSV (.csv): sep=";" and decimal=","
#
# Requirements:
#   pip install numpy pandas trimesh tifffile scipy openpyxl
# ============================================================

import os
import time
import numpy as np
import pandas as pd
import trimesh
from tifffile import imread
from scipy.ndimage import distance_transform_edt


# ---------------- SETTINGS ----------------
ROOT_DIR = r"<BEETLE_JOINTS_ROOT>\Processed\Coxa"

OUT_XLSX = os.path.join(ROOT_DIR, "coxa_combined_metrics.xlsx")
OUT_CSV  = os.path.join(ROOT_DIR, "coxa_combined_metrics.csv")

VOXEL_SIZE_UM = 1.22

# Ignore tiny speckles/noise: slice must have at least this many foreground pixels
MIN_PIXELS_PER_SLICE = 200

# Round outputs to avoid Excel scientific notation and locale weirdness
ROUND_UM = 3      # µm columns
ROUND_VOX = 3     # voxel columns

# ------------------------------------------


def find_mid_slice(binary_stack: np.ndarray):
    """
    binary_stack: (Z, Y, X) boolean
    returns first_z, last_z, mid_z based on slices that contain enough cuticle pixels
    """
    counts = binary_stack.reshape(binary_stack.shape[0], -1).sum(axis=1)
    valid = np.where(counts >= MIN_PIXELS_PER_SLICE)[0]
    if len(valid) == 0:
        return None, None, None
    first_z = int(valid[0])
    last_z  = int(valid[-1])
    mid_z   = int(round((first_z + last_z) / 2))
    return first_z, last_z, mid_z


def bbox_diag_from_obj(obj_path: str) -> float:
    """
    Returns bounding-box diagonal in OBJ coordinate units.
    Robust even if mesh is not watertight.
    """
    mesh = trimesh.load(obj_path, force="mesh", process=False)
    verts = np.asarray(mesh.vertices, dtype=float)
    if verts.size == 0:
        raise ValueError("Empty mesh (no vertices)")

    vmin = verts.min(axis=0)
    vmax = verts.max(axis=0)
    extents = vmax - vmin
    diag = float(np.linalg.norm(extents))
    return diag


def thickness_stats_from_tif(tif_path: str):
    """
    Loads 3D label stack (Z,Y,X), finds mid-slice of structure,
    computes thickness per cuticle pixel via EDT on cropped slice,
    returns median/p10/p90 in voxels and slice indices.
    """
    stack = imread(tif_path)
    if stack.ndim != 3:
        raise ValueError(f"Expected 3D stack (Z,Y,X), got shape {stack.shape}")

    cuticle_3d = stack > 0
    first_z, last_z, mid_z = find_mid_slice(cuticle_3d)
    if mid_z is None:
        raise ValueError("No valid slices found (check MIN_PIXELS_PER_SLICE or labels)")

    sl = cuticle_3d[mid_z, :, :]
    n_pix = int(sl.sum())
    if n_pix == 0:
        raise ValueError("Mid-slice has no cuticle pixels")

    # Crop to object bounding box -> avoids 1-minute outliers
    ys, xs = np.where(sl)
    y0, y1 = ys.min(), ys.max()
    x0, x1 = xs.min(), xs.max()
    sl_crop = sl[y0:y1+1, x0:x1+1]

    # EDT (distance to background), thickness estimate = 2*distance
    dt = distance_transform_edt(sl_crop)
    thickness_vox = 2.0 * dt[sl_crop]

    med_vox = float(np.median(thickness_vox))
    p10_vox = float(np.percentile(thickness_vox, 10))
    p90_vox = float(np.percentile(thickness_vox, 90))

    return {
        "first_slice": first_z,
        "last_slice": last_z,
        "mid_slice": mid_z,
        "n_pixels_mid_slice": n_pix,
        "median_thickness_vox": med_vox,
        "p10_thickness_vox": p10_vox,
        "p90_thickness_vox": p90_vox,
    }


def safe_round(x, nd):
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return np.nan
        return float(np.round(x, nd))
    except Exception:
        return np.nan


def main():
    t_start = time.time()

    files = os.listdir(ROOT_DIR)
    obj_bases = {os.path.splitext(f)[0] for f in files if f.lower().endswith(".obj")}
    tif_bases = {os.path.splitext(f)[0] for f in files if f.lower().endswith(".tif") or f.lower().endswith(".tiff")}

    paired = sorted(list(obj_bases & tif_bases))
    if not paired:
        raise RuntimeError("No paired .obj + .tif found with matching base names in ROOT_DIR.")

    print(f"Found pairs: {len(paired)}")

    rows = []
    errors = []

    for base in paired:
        obj_path = os.path.join(ROOT_DIR, base + ".obj")
        tif_path = os.path.join(ROOT_DIR, base + ".tif")
        if not os.path.exists(tif_path):
            # fallback .tiff
            tif_path = os.path.join(ROOT_DIR, base + ".tiff")

        t0 = time.time()
        try:
            bbox_diag_units = bbox_diag_from_obj(obj_path)
            bbox_diag_um = bbox_diag_units * VOXEL_SIZE_UM  # your meshes appear voxel-based; keep this conversion

            thick = thickness_stats_from_tif(tif_path)

            med_um = thick["median_thickness_vox"] * VOXEL_SIZE_UM
            p10_um = thick["p10_thickness_vox"] * VOXEL_SIZE_UM
            p90_um = thick["p90_thickness_vox"] * VOXEL_SIZE_UM

            row = {
                "specimen": base,

                # size
                "bbox_diag_units": safe_round(bbox_diag_units, ROUND_VOX),
                "bbox_diag_um": safe_round(bbox_diag_um, ROUND_UM),

                # thickness (vox + µm)
                "median_thickness_vox": safe_round(thick["median_thickness_vox"], ROUND_VOX),
                "p10_thickness_vox": safe_round(thick["p10_thickness_vox"], ROUND_VOX),
                "p90_thickness_vox": safe_round(thick["p90_thickness_vox"], ROUND_VOX),

                "median_thickness_um": safe_round(med_um, ROUND_UM),
                "p10_thickness_um": safe_round(p10_um, ROUND_UM),
                "p90_thickness_um": safe_round(p90_um, ROUND_UM),

                # slice info
                "first_slice": thick["first_slice"],
                "last_slice": thick["last_slice"],
                "mid_slice": thick["mid_slice"],
                "n_pixels_mid_slice": thick["n_pixels_mid_slice"],

                # file paths (optional but helpful)
                "obj_file": os.path.basename(obj_path),
                "tif_file": os.path.basename(tif_path),
            }

            rows.append(row)
            print(f"✓ {base}  ({time.time()-t0:.2f}s)")

        except Exception as e:
            msg = str(e)
            print(f"✗ {base}: {msg}")
            errors.append({"specimen": base, "error": msg})
            continue

    df = pd.DataFrame(rows).sort_values("specimen")

    # Write XLSX (best for Excel)
    df.to_excel(OUT_XLSX, index=False)
    print(f"\nWrote XLSX: {OUT_XLSX}")

    # Write German CSV (semicolon + decimal comma)
    # rounding already prevents scientific notation in most cases
    df.to_csv(OUT_CSV, index=False, sep=";", decimal=",")
    print(f"Wrote CSV:  {OUT_CSV}")

    # Optional error log
    if errors:
        err_df = pd.DataFrame(errors).sort_values("specimen")
        err_path = os.path.join(ROOT_DIR, "coxa_combined_metrics_ERRORS.xlsx")
        err_df.to_excel(err_path, index=False)
        print(f"Wrote ERRORS: {err_path} ({len(errors)} specimens failed)")

    print(f"\nDone in {time.time()-t_start:.1f}s. Processed: {len(df)} / {len(paired)}")


if __name__ == "__main__":
    main()
