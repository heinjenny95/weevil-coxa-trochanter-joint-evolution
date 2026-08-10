#!/usr/bin/env python3
"""Create downstream-compatible geometry tables from robust helix fits."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Sequence


SPECIMEN_ID_ALIASES = {
    "308_lisshorhoptrus_oryzophilus_aligned": (
        "308_lisshorhoptrus_oryzophilus_trochanter_aligned"
    ),
    "310_ormiscus_saltator_trochanter_mirrored_aligned": (
        "310_ormiscus_saltator_trochanter_aligned"
    ),
    "t_pseudonastus_trochanter_aligned": (
        "t_pseudonasutus_trochanter_aligned"
    ),
}


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "").strip()
    if not value:
        return math.nan
    return float(value)


def as_bool(row: dict[str, str], key: str) -> bool:
    return row.get(key, "").strip().lower() in {"true", "1", "yes"}


def compatibility_record(row: dict[str, str]) -> dict[str, object]:
    signed_angle = as_float(row, "signed_winding_angle_deg")
    abs_angle = as_float(row, "abs_winding_angle_deg")
    endpoint_span = as_float(row, "endpoint_axial_span_new_axis")
    fitted_pitch = as_float(row, "fitted_pitch_360")
    helix_rms = as_float(row, "helix_rms")
    return {
        "specimen_id": SPECIMEN_ID_ALIASES.get(
            row["specimen_id"], row["specimen_id"]
        ),
        "raw_geometry_id": row["specimen_id"],
        "signed_winding_angle_deg": signed_angle,
        "abs_winding_angle_deg": abs_angle,
        "n_turns_signed": signed_angle / 360.0,
        "n_turns_abs": abs_angle / 360.0,
        "start_end_dist": as_float(row, "start_end_dist"),
        # Retain the historical meaning: observed endpoint separation projected
        # on the fitted axis. Do not replace it with fitted pitch * turns.
        "axial_span": endpoint_span,
        "fitted_axial_span": as_float(row, "fitted_axial_span"),
        # This is the robust slope through all ordered points, not an endpoint
        # quotient. Both aliases are supplied for existing downstream readers.
        "axial_pitch_360": fitted_pitch,
        "fitted_pitch_360": fitted_pitch,
        "fit_radius": as_float(row, "fit_radius"),
        # Compatibility alias now represents total 3D helix RMS. Radial and
        # axial components remain explicit in adjacent columns.
        "fit_rms": helix_rms,
        "radial_rms": as_float(row, "radial_rms"),
        "axial_rms": as_float(row, "axial_rms"),
        "helix_rms": helix_rms,
        "helix_rms_relative_to_radius": as_float(
            row, "helix_rms_relative_to_radius"
        ),
        "axial_angle_r_squared": as_float(row, "axial_angle_r_squared"),
        "quality_class": row.get("quality_class", ""),
        "quality_warnings": row.get("quality_warnings", ""),
        "analysis_set_all": as_bool(row, "analysis_set_all"),
        "analysis_set_primary_adequate": as_bool(
            row, "analysis_set_primary_adequate"
        ),
        "analysis_set_strict_good": as_bool(row, "analysis_set_strict_good"),
        "measurement_method": "robust_circular_3d_helix_soft_l1",
        "uncertainty_scope": row.get("uncertainty_interpretation", ""),
    }


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    with arguments.metrics_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        input_rows = list(csv.DictReader(handle))
    if not input_rows:
        raise ValueError("Robust metrics table is empty")
    records = [compatibility_record(row) for row in input_rows]
    sets = {
        "all": records,
        "primary_adequate": [
            record for record in records if record["analysis_set_primary_adequate"]
        ],
        "strict_good": [
            record for record in records if record["analysis_set_strict_good"]
        ],
    }
    expected = {"all": 64, "primary_adequate": 63, "strict_good": 53}
    for set_name, set_records in sets.items():
        if len(set_records) != expected[set_name]:
            raise ValueError(
                f"Unexpected {set_name} size: {len(set_records)}; "
                f"expected {expected[set_name]}"
            )
        output_path = arguments.output_dir / f"robust_geometry_{set_name}.csv"
        write_csv(output_path, set_records)
        print(f"{set_name}: {len(set_records)} -> {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
