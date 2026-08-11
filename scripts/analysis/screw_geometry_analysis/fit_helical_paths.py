#!/usr/bin/env python3
"""Robustly fit circular helices to ordered 3D semilandmark traces.

This workflow deliberately keeps the historical Cinema 4D measurement model
and the new model separate:

* ``legacy`` ports the released Cinema 4D grid-search algorithm so that its
  output can be reproduced outside Cinema 4D.
* ``robust_helix`` continuously optimizes a three-dimensional circular helix.
  It minimizes radial and axial residuals together with a robust loss.

The input CSV files are the Cinema 4D landmark exports used by this project.
The first line declares the number of landmarks, the second line contains
``X``, ``Y`` and ``Z`` columns, and subsequent rows are ordered along the
trajectory. Coordinates are assumed to be in millimetres unless the export
states otherwise.

The bootstrap intervals quantify conditional path-fit uncertainty. They do
not estimate manual landmark-placement repeatability, which requires repeated
independent tracings.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import re
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy import __version__ as scipy_version
from scipy.optimize import least_squares


ALGORITHM_VERSION = "1.0.0"
EPS = 1e-12


@dataclass
class LandmarkTrace:
    specimen_id: str
    source_path: Path
    point_names: list[str]
    points: np.ndarray
    normals: np.ndarray | None
    declared_count: int | None
    units: str


@dataclass
class LegacyFit:
    axis_point: np.ndarray
    axis: np.ndarray
    basis_u: np.ndarray
    basis_v: np.ndarray
    theta: np.ndarray
    signed_winding_angle_deg: float
    abs_winding_angle_deg: float
    axial_span: float
    start_end_dist: float
    radius: float
    radial_rms: float


@dataclass
class RobustHelixFit:
    success: bool
    message: str
    axis_point: np.ndarray
    axis: np.ndarray
    basis_e1: np.ndarray
    basis_e2: np.ndarray
    radius: float
    axial_intercept: float
    axial_rise_per_radian: float
    theta: np.ndarray
    radial_coordinates: np.ndarray
    axial_coordinates: np.ndarray
    predicted_points: np.ndarray
    radial_residuals: np.ndarray
    axial_residuals: np.ndarray
    point_distances: np.ndarray
    objective_cost: float
    nfev: int
    data_scale: float
    candidate_count: int
    second_best_cost_ratio: float
    near_optimal_solution_count: int
    near_optimal_angle_range_deg: float
    near_optimal_pitch_range: float
    near_optimal_axis_max_separation_deg: float
    warnings: list[str] = field(default_factory=list)

    @property
    def signed_winding_angle_deg(self) -> float:
        return float(np.degrees(self.theta[-1] - self.theta[0]))

    @property
    def abs_winding_angle_deg(self) -> float:
        return abs(self.signed_winding_angle_deg)

    @property
    def n_turns_abs(self) -> float:
        return self.abs_winding_angle_deg / 360.0

    @property
    def fitted_pitch_360(self) -> float:
        return 2.0 * math.pi * abs(self.axial_rise_per_radian)

    @property
    def fitted_axial_span(self) -> float:
        return abs(
            self.axial_rise_per_radian * (self.theta[-1] - self.theta[0])
        )

    @property
    def radial_rms(self) -> float:
        return float(np.sqrt(np.mean(np.square(self.radial_residuals))))

    @property
    def axial_rms(self) -> float:
        return float(np.sqrt(np.mean(np.square(self.axial_residuals))))

    @property
    def helix_rms(self) -> float:
        return float(np.sqrt(np.mean(np.square(self.point_distances))))


def normalize(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm <= EPS:
        raise ValueError("Cannot normalize a zero or non-finite vector")
    return vector / norm


def canonicalize_axis(
    axis: np.ndarray, e1: np.ndarray, e2: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Choose a deterministic sign while preserving a right-handed basis."""
    axis = normalize(axis)
    dominant = int(np.argmax(np.abs(axis)))
    if axis[dominant] < 0.0:
        return -axis, e1, -e2
    return axis, e1, e2


def orthonormal_basis(
    axis: np.ndarray, reference: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    axis = normalize(axis)
    if reference is None:
        candidates = np.eye(3)
        reference = candidates[int(np.argmin(np.abs(candidates @ axis)))]
    reference = np.asarray(reference, dtype=float)
    e1 = reference - axis * float(np.dot(reference, axis))
    if np.linalg.norm(e1) <= 1e-8:
        candidates = np.eye(3)
        reference = candidates[int(np.argmin(np.abs(candidates @ axis)))]
        e1 = reference - axis * float(np.dot(reference, axis))
    e1 = normalize(e1)
    e2 = normalize(np.cross(axis, e1))
    return e1, e2


def rotate_vector(vector: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    axis = normalize(axis)
    return (
        vector * math.cos(angle)
        + np.cross(axis, vector) * math.sin(angle)
        + axis * float(np.dot(axis, vector)) * (1.0 - math.cos(angle))
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_column(name: str) -> str:
    return re.sub(r"[\s_'\"]", "", name.strip().lower())


def read_landmark_csv(path: Path) -> LandmarkTrace:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))

    rows = [row for row in rows if any(cell.strip() for cell in row)]
    if len(rows) < 3:
        raise ValueError(f"Too few rows in {path}")

    declared_count: int | None = None
    units = "unknown"
    if rows[0] and rows[0][0].strip().lower() == "landmarks":
        if len(rows[0]) > 1:
            try:
                declared_count = int(rows[0][1])
            except ValueError:
                declared_count = None
        if len(rows[0]) > 3:
            units = rows[0][3].strip().strip("()") or "unknown"

    header_index = None
    indices: dict[str, int] = {}
    aliases = {
        "x": {"x", "pos.x", "posx", "px", "coordx"},
        "y": {"y", "pos.y", "posy", "py", "coordy"},
        "z": {"z", "pos.z", "posz", "pz", "coordz"},
    }
    for row_index, row in enumerate(rows[:12]):
        canonical = [canonical_column(cell) for cell in row]
        found: dict[str, int] = {}
        for coordinate, accepted in aliases.items():
            for column_index, value in enumerate(canonical):
                if value in accepted:
                    found[coordinate] = column_index
                    break
        if len(found) == 3:
            header_index = row_index
            indices.update(found)
            for optional in ("name", "nx", "ny", "nz"):
                if optional in canonical:
                    indices[optional] = canonical.index(optional)
            break

    if header_index is None:
        raise ValueError(f"Could not locate X/Y/Z columns in {path}")

    points: list[list[float]] = []
    normals: list[list[float]] = []
    names: list[str] = []
    normal_columns_present = all(key in indices for key in ("nx", "ny", "nz"))
    required_max = max(indices["x"], indices["y"], indices["z"])

    for row_number, row in enumerate(rows[header_index + 1 :], start=1):
        if len(row) <= required_max:
            continue
        try:
            point = [
                float(row[indices[coordinate]].strip().replace(",", "."))
                for coordinate in ("x", "y", "z")
            ]
        except ValueError:
            continue
        points.append(point)
        if "name" in indices and len(row) > indices["name"]:
            names.append(row[indices["name"]].strip() or f"point_{row_number:03d}")
        else:
            names.append(f"point_{row_number:03d}")
        if normal_columns_present:
            try:
                normals.append(
                    [
                        float(row[indices[coordinate]].strip().replace(",", "."))
                        for coordinate in ("nx", "ny", "nz")
                    ]
                )
            except ValueError:
                normals.append([math.nan, math.nan, math.nan])

    point_array = np.asarray(points, dtype=float)
    if point_array.ndim != 2 or point_array.shape[0] < 5 or point_array.shape[1] != 3:
        raise ValueError(f"Need at least five valid XYZ points in {path}")
    if not np.all(np.isfinite(point_array)):
        raise ValueError(f"Non-finite XYZ coordinates in {path}")
    if declared_count is not None and declared_count != point_array.shape[0]:
        raise ValueError(
            f"Declared {declared_count} landmarks but read {point_array.shape[0]} in {path}"
        )

    normal_array = np.asarray(normals, dtype=float) if normal_columns_present else None
    return LandmarkTrace(
        specimen_id=path.stem,
        source_path=path,
        point_names=names,
        points=point_array,
        normals=normal_array,
        declared_count=declared_count,
        units=units,
    )


def fit_circle_algebraic(coordinates: np.ndarray) -> tuple[float, float, float]:
    coordinates = np.asarray(coordinates, dtype=float)
    x = coordinates[:, 0]
    y = coordinates[:, 1]
    design = np.column_stack((x, y, np.ones_like(x)))
    response = -(np.square(x) + np.square(y))
    coefficients, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    if rank < 3:
        raise ValueError("Projected points do not define a circle")
    a_coefficient, b_coefficient, c_coefficient = coefficients
    center_x = -a_coefficient / 2.0
    center_y = -b_coefficient / 2.0
    radius_squared = (
        center_x * center_x + center_y * center_y - c_coefficient
    )
    if not np.isfinite(radius_squared) or radius_squared <= EPS:
        raise ValueError("Invalid fitted circle radius")
    return float(center_x), float(center_y), math.sqrt(float(radius_squared))


def legacy_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    axis = normalize(axis)
    temporary = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, temporary))) > 0.9:
        temporary = np.array([0.0, 1.0, 0.0])
    basis_u = normalize(temporary - axis * float(np.dot(axis, temporary)))
    basis_v = normalize(np.cross(axis, basis_u))
    return basis_u, basis_v


def evaluate_legacy_axis(points: np.ndarray, axis: np.ndarray) -> dict[str, object]:
    axis = normalize(axis)
    origin = np.mean(points, axis=0)
    basis_u, basis_v = legacy_basis(axis)
    centered = points - origin
    coordinates = np.column_stack((centered @ basis_u, centered @ basis_v))
    center_x, center_y, _ = fit_circle_algebraic(coordinates)
    axis_point = origin + center_x * basis_u + center_y * basis_v
    centered_on_axis = points - axis_point
    coordinates = np.column_stack(
        (centered_on_axis @ basis_u, centered_on_axis @ basis_v)
    )
    radii = np.linalg.norm(coordinates, axis=1)
    radius = float(np.mean(radii))
    rms = float(np.sqrt(np.mean(np.square(radii - radius))))
    return {
        "axis": axis,
        "axis_point": axis_point,
        "basis_u": basis_u,
        "basis_v": basis_v,
        "coordinates": coordinates,
        "radius": radius,
        "rms": rms,
    }


def refine_legacy_axis_grid(
    points: np.ndarray, base_axis: np.ndarray, range_degrees: float, step_degrees: float
) -> dict[str, object]:
    base_axis = normalize(base_axis)
    basis_1, basis_2 = legacy_basis(base_axis)
    values = np.arange(
        -range_degrees, range_degrees + step_degrees * 0.25, step_degrees
    )
    best: dict[str, object] | None = None
    for angle_1 in values:
        for angle_2 in values:
            candidate = rotate_vector(base_axis, basis_1, math.radians(float(angle_1)))
            candidate = rotate_vector(candidate, basis_2, math.radians(float(angle_2)))
            try:
                evaluated = evaluate_legacy_axis(points, candidate)
            except (ValueError, np.linalg.LinAlgError):
                continue
            if best is None or float(evaluated["rms"]) < float(best["rms"]):
                best = evaluated
    if best is None:
        raise RuntimeError("Legacy axis grid search failed")
    return best


def fit_legacy_cinema_method(points: np.ndarray) -> LegacyFit:
    centered = points - np.mean(points, axis=0)
    covariance = centered.T @ centered / float(points.shape[0])
    _, eigenvectors = np.linalg.eigh(covariance)
    eigenvectors = eigenvectors[:, ::-1]

    best: dict[str, object] | None = None
    for index in range(3):
        for sign in (1.0, -1.0):
            seed = sign * eigenvectors[:, index]
            try:
                candidate = refine_legacy_axis_grid(points, seed, 40.0, 8.0)
            except RuntimeError:
                continue
            if best is None or float(candidate["rms"]) < float(best["rms"]):
                best = candidate
    if best is None:
        raise RuntimeError("No usable legacy initial axis")

    best = refine_legacy_axis_grid(points, np.asarray(best["axis"]), 10.0, 2.0)
    best = refine_legacy_axis_grid(points, np.asarray(best["axis"]), 2.0, 0.5)
    coordinates = np.asarray(best["coordinates"])
    theta = np.unwrap(np.arctan2(coordinates[:, 1], coordinates[:, 0]))
    winding = float(np.degrees(theta[-1] - theta[0]))
    axis = normalize(np.asarray(best["axis"]))
    return LegacyFit(
        axis_point=np.asarray(best["axis_point"]),
        axis=axis,
        basis_u=np.asarray(best["basis_u"]),
        basis_v=np.asarray(best["basis_v"]),
        theta=theta,
        signed_winding_angle_deg=winding,
        abs_winding_angle_deg=abs(winding),
        axial_span=abs(float(np.dot(points[-1] - points[0], axis))),
        start_end_dist=float(np.linalg.norm(points[-1] - points[0])),
        radius=float(best["radius"]),
        radial_rms=float(best["rms"]),
    )


def data_normalization(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    center = np.mean(points, axis=0)
    centered = points - center
    scale = float(np.sqrt(np.mean(np.sum(np.square(centered), axis=1))))
    if not np.isfinite(scale) or scale <= EPS:
        raise ValueError("Landmark trace has zero spatial extent")
    return centered / scale, center, scale


def axis_from_tangent_parameters(
    seed_axis: np.ndarray, tangent_1: np.ndarray, tangent_2: np.ndarray, a: float, b: float
) -> np.ndarray:
    return normalize(seed_axis + a * tangent_1 + b * tangent_2)


def initial_helix_parameters(
    normalized_points: np.ndarray, seed_axis: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    seed_axis = normalize(seed_axis)
    tangent_1, tangent_2 = orthonormal_basis(seed_axis)
    projected = np.column_stack(
        (normalized_points @ tangent_1, normalized_points @ tangent_2)
    )
    center_x, center_y, _ = fit_circle_algebraic(projected)
    shifted = projected - np.array([center_x, center_y])
    radii = np.linalg.norm(shifted, axis=1)
    radius = max(float(np.median(radii)), 1e-6)
    theta = np.unwrap(np.arctan2(shifted[:, 1], shifted[:, 0]))
    theta_centered = theta - float(np.mean(theta))
    axial = normalized_points @ seed_axis
    if float(np.ptp(theta_centered)) > 1e-8:
        slope, intercept = np.polyfit(theta_centered, axial, 1)
    else:
        slope, intercept = 0.0, float(np.mean(axial))
    parameters = np.array(
        [0.0, 0.0, center_x, center_y, math.log(radius), intercept, slope],
        dtype=float,
    )
    return parameters, tangent_1


def helix_components(
    parameters: np.ndarray,
    normalized_points: np.ndarray,
    seed_axis: np.ndarray,
    tangent_1: np.ndarray,
) -> dict[str, np.ndarray | float]:
    tangent_2 = normalize(np.cross(seed_axis, tangent_1))
    axis = axis_from_tangent_parameters(
        seed_axis, tangent_1, tangent_2, parameters[0], parameters[1]
    )
    basis_e1, basis_e2 = orthonormal_basis(axis, tangent_1)
    axis_point = parameters[2] * basis_e1 + parameters[3] * basis_e2
    relative = normalized_points - axis_point
    x = relative @ basis_e1
    y = relative @ basis_e2
    axial = relative @ axis
    radius_coordinates = np.sqrt(np.square(x) + np.square(y))
    theta = np.unwrap(np.arctan2(y, x))
    theta_centered = theta - float(np.mean(theta))
    radius = math.exp(float(parameters[4]))
    fitted_axial = parameters[5] + parameters[6] * theta_centered
    radial_residuals = radius_coordinates - radius
    axial_residuals = axial - fitted_axial
    return {
        "axis": axis,
        "basis_e1": basis_e1,
        "basis_e2": basis_e2,
        "axis_point": axis_point,
        "x": x,
        "y": y,
        "axial": axial,
        "radius_coordinates": radius_coordinates,
        "theta": theta,
        "theta_centered": theta_centered,
        "radius": radius,
        "fitted_axial": fitted_axial,
        "radial_residuals": radial_residuals,
        "axial_residuals": axial_residuals,
    }


def robust_residual_vector(
    parameters: np.ndarray,
    normalized_points: np.ndarray,
    seed_axis: np.ndarray,
    tangent_1: np.ndarray,
) -> np.ndarray:
    components = helix_components(parameters, normalized_points, seed_axis, tangent_1)
    return np.concatenate(
        (
            np.asarray(components["radial_residuals"]),
            np.asarray(components["axial_residuals"]),
        )
    )


def provisional_quality_warnings(
    fit: RobustHelixFit, points: np.ndarray
) -> list[str]:
    warnings: list[str] = []
    n_points = points.shape[0]
    trajectory_span = float(np.linalg.norm(points[-1] - points[0]))
    relative_rms = fit.helix_rms / fit.radius if fit.radius > EPS else math.inf
    angle = fit.abs_winding_angle_deg
    theta_differences = np.diff(fit.theta)
    direction = np.sign(float(np.sum(theta_differences))) or 1.0
    backtracking = float(np.mean(theta_differences * direction < -1e-5))
    axial_total = float(np.sum(np.square(fit.axial_coordinates - np.mean(fit.axial_coordinates))))
    axial_r2 = (
        1.0 - float(np.sum(np.square(fit.axial_residuals))) / axial_total
        if axial_total > EPS
        else math.nan
    )

    if n_points < 10:
        warnings.append("very_few_points")
    elif n_points < 20:
        warnings.append("few_points")
    if angle < 30.0:
        warnings.append("angular_coverage_below_30_deg")
    elif angle < 180.0:
        warnings.append("angular_coverage_below_180_deg")
    if trajectory_span > EPS and fit.radius / trajectory_span > 5.0:
        warnings.append("large_radius_relative_to_trace")
    if relative_rms > 0.10:
        warnings.append("high_relative_helix_rms")
    if np.isfinite(axial_r2) and axial_r2 < 0.80:
        warnings.append("weak_axial_angle_linearity")
    if backtracking > 0.10:
        warnings.append("angular_backtracking")
    if fit.near_optimal_angle_range_deg > 20.0:
        warnings.append("near_optimal_angle_instability")
    if fit.near_optimal_axis_max_separation_deg > 15.0:
        warnings.append("near_optimal_axis_instability")
    return warnings


def fit_robust_helix(
    points: np.ndarray,
    initial_axes: Sequence[np.ndarray] | None = None,
    robust_loss: str = "soft_l1",
    f_scale: float = 0.02,
    max_nfev: int = 4000,
    include_pca_seeds: bool = True,
) -> RobustHelixFit:
    normalized_points, original_center, data_scale = data_normalization(points)
    covariance = normalized_points.T @ normalized_points / float(points.shape[0])
    _, eigenvectors = np.linalg.eigh(covariance)
    eigenvectors = eigenvectors[:, ::-1]

    seeds: list[np.ndarray] = []
    if initial_axes:
        seeds.extend(normalize(np.asarray(axis)) for axis in initial_axes)
    if include_pca_seeds or not seeds:
        for index in range(3):
            seeds.extend((eigenvectors[:, index], -eigenvectors[:, index]))

    unique_seeds: list[np.ndarray] = []
    for seed in seeds:
        seed = normalize(seed)
        if not any(abs(float(np.dot(seed, prior))) > 0.999999 for prior in unique_seeds):
            unique_seeds.append(seed)

    lower = np.array([-5.0, -5.0, -100.0, -100.0, -12.0, -20.0, -50.0])
    upper = np.array([5.0, 5.0, 100.0, 100.0, math.log(100.0), 20.0, 50.0])
    candidates: list[tuple[float, object, dict[str, np.ndarray | float], np.ndarray]] = []

    for seed in unique_seeds:
        try:
            initial, tangent_1 = initial_helix_parameters(normalized_points, seed)
            initial = np.minimum(np.maximum(initial, lower + 1e-8), upper - 1e-8)
            result = least_squares(
                robust_residual_vector,
                initial,
                args=(normalized_points, seed, tangent_1),
                bounds=(lower, upper),
                loss=robust_loss,
                f_scale=f_scale,
                x_scale="jac",
                max_nfev=max_nfev,
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
            components = helix_components(result.x, normalized_points, seed, tangent_1)
            residuals = np.concatenate(
                (
                    np.asarray(components["radial_residuals"]),
                    np.asarray(components["axial_residuals"]),
                )
            )
            score = float(result.cost)
            if np.isfinite(score):
                candidates.append((score, result, components, result.x))
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError):
            continue

    if not candidates:
        raise RuntimeError("All robust helix starting configurations failed")
    candidates.sort(key=lambda item: item[0])
    best_score, result, components, parameters = candidates[0]

    second_best_cost_ratio = (
        float(candidates[1][0] / best_score)
        if len(candidates) > 1 and best_score > EPS
        else math.nan
    )
    near_optimal = [
        candidate
        for candidate in candidates
        if candidate[0] <= best_score * 1.05 + 1e-12
    ]
    near_angles = np.asarray(
        [
            abs(
                math.degrees(
                    float(
                        np.asarray(candidate[2]["theta"])[-1]
                        - np.asarray(candidate[2]["theta"])[0]
                    )
                )
            )
            for candidate in near_optimal
        ],
        dtype=float,
    )
    near_pitches = np.asarray(
        [2.0 * math.pi * abs(float(candidate[3][6])) * data_scale for candidate in near_optimal],
        dtype=float,
    )
    best_axis_for_diagnostics = normalize(np.asarray(components["axis"]))
    near_axis_separations = np.asarray(
        [
            math.degrees(
                math.acos(
                    float(
                        np.clip(
                            abs(
                                np.dot(
                                    best_axis_for_diagnostics,
                                    normalize(np.asarray(candidate[2]["axis"])),
                                )
                            ),
                            -1.0,
                            1.0,
                        )
                    )
                )
            )
            for candidate in near_optimal
        ],
        dtype=float,
    )

    axis = np.asarray(components["axis"])
    basis_e1 = np.asarray(components["basis_e1"])
    basis_e2 = np.asarray(components["basis_e2"])
    axis_before_canonicalization = axis.copy()
    axis, basis_e1, basis_e2 = canonicalize_axis(axis, basis_e1, basis_e2)
    axis_was_flipped = float(np.dot(axis, axis_before_canonicalization)) < 0.0

    normalized_axis_point = np.asarray(components["axis_point"])
    axis_point = original_center + data_scale * normalized_axis_point
    relative = (points - axis_point) / data_scale
    x = relative @ basis_e1
    y = relative @ basis_e2
    axial_normalized = relative @ axis
    theta = np.unwrap(np.arctan2(y, x))
    theta_centered = theta - float(np.mean(theta))

    radius = float(components["radius"]) * data_scale
    # Preserve the robustly optimized axial regression under the deterministic
    # axis-sign transformation. Flipping both axis and angular direction changes
    # the intercept sign but leaves rise per radian unchanged.
    axial_intercept_normalized = (
        -float(parameters[5]) if axis_was_flipped else float(parameters[5])
    )
    axial_rise_normalized = float(parameters[6])
    fitted_axial = axial_intercept_normalized + axial_rise_normalized * theta_centered
    axial_intercept = axial_intercept_normalized * data_scale
    axial_rise_per_radian = axial_rise_normalized * data_scale
    radius_coordinates = np.sqrt(np.square(x) + np.square(y)) * data_scale
    radial_residuals = radius_coordinates - radius
    axial_residuals = (axial_normalized - fitted_axial) * data_scale

    radial_unit = np.column_stack(
        (np.cos(theta), np.sin(theta))
    )
    predicted_points = (
        axis_point
        + np.outer(fitted_axial * data_scale, axis)
        + radius
        * (
            np.outer(radial_unit[:, 0], basis_e1)
            + np.outer(radial_unit[:, 1], basis_e2)
        )
    )
    point_distances = np.linalg.norm(points - predicted_points, axis=1)

    fit = RobustHelixFit(
        success=bool(result.success),
        message=str(result.message),
        axis_point=axis_point,
        axis=axis,
        basis_e1=basis_e1,
        basis_e2=basis_e2,
        radius=radius,
        axial_intercept=axial_intercept,
        axial_rise_per_radian=axial_rise_per_radian,
        theta=theta,
        radial_coordinates=radius_coordinates,
        axial_coordinates=axial_normalized * data_scale,
        predicted_points=predicted_points,
        radial_residuals=radial_residuals,
        axial_residuals=axial_residuals,
        point_distances=point_distances,
        objective_cost=float(result.cost),
        nfev=int(result.nfev),
        data_scale=data_scale,
        candidate_count=len(candidates),
        second_best_cost_ratio=second_best_cost_ratio,
        near_optimal_solution_count=len(near_optimal),
        near_optimal_angle_range_deg=float(np.ptp(near_angles)) if len(near_angles) else math.nan,
        near_optimal_pitch_range=float(np.ptp(near_pitches)) if len(near_pitches) else math.nan,
        near_optimal_axis_max_separation_deg=float(np.max(near_axis_separations)) if len(near_axis_separations) else math.nan,
    )
    fit.warnings = provisional_quality_warnings(fit, points)
    return fit


def quality_metrics(fit: RobustHelixFit, points: np.ndarray) -> dict[str, float | str]:
    theta_differences = np.diff(fit.theta)
    direction = np.sign(float(np.sum(theta_differences))) or 1.0
    backtracking_fraction = float(np.mean(theta_differences * direction < -1e-5))
    axial_total = float(
        np.sum(np.square(fit.axial_coordinates - np.mean(fit.axial_coordinates)))
    )
    axial_r2 = (
        1.0 - float(np.sum(np.square(fit.axial_residuals))) / axial_total
        if axial_total > EPS
        else math.nan
    )
    radius_cv = (
        float(np.std(fit.radial_coordinates, ddof=1) / np.mean(fit.radial_coordinates))
        if len(fit.radial_coordinates) > 1 and np.mean(fit.radial_coordinates) > EPS
        else math.nan
    )
    trajectory_length = float(np.sum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
    chord = float(np.linalg.norm(points[-1] - points[0]))
    relative_rms = fit.helix_rms / fit.radius if fit.radius > EPS else math.nan
    if any(
        warning in fit.warnings
        for warning in (
            "very_few_points",
            "angular_coverage_below_30_deg",
            "angular_coverage_below_180_deg",
            "large_radius_relative_to_trace",
            "near_optimal_angle_instability",
            "near_optimal_axis_instability",
        )
    ):
        quality_class = "limited_identifiability"
    elif fit.warnings:
        quality_class = "caution"
    else:
        quality_class = "good"
    return {
        "axial_angle_r_squared": axial_r2,
        "radius_cv": radius_cv,
        "angular_backtracking_fraction": backtracking_fraction,
        "trajectory_length": trajectory_length,
        "start_end_dist": chord,
        "helix_rms_relative_to_radius": relative_rms,
        "quality_class": quality_class,
        "quality_warnings": "|".join(fit.warnings),
    }


def moving_block_indices(n: int, block_length: int, rng: np.random.Generator) -> np.ndarray:
    output: list[int] = []
    while len(output) < n:
        start = int(rng.integers(0, n))
        output.extend((start + offset) % n for offset in range(block_length))
    return np.asarray(output[:n], dtype=int)


def bootstrap_fit_intervals(
    points: np.ndarray,
    fit: RobustHelixFit,
    replicates: int,
    seed: int,
) -> tuple[dict[str, float | int], list[dict[str, float | int]]]:
    if replicates <= 0:
        return {"bootstrap_requested": 0, "bootstrap_successful": 0}, []
    rng = np.random.default_rng(seed)
    residuals = points - fit.predicted_points
    residuals = residuals - np.mean(residuals, axis=0)
    block_length = max(2, int(round(math.sqrt(points.shape[0]))))
    collected: dict[str, list[float]] = {
        "abs_winding_angle_deg": [],
        "fitted_pitch_360": [],
        "endpoint_axial_span": [],
        "fitted_axial_span": [],
        "radius": [],
        "helix_rms": [],
    }
    draws: list[dict[str, float | int]] = []
    for replicate_index in range(1, replicates + 1):
        indices = moving_block_indices(points.shape[0], block_length, rng)
        synthetic_points = fit.predicted_points + residuals[indices]
        try:
            bootstrap_fit = fit_robust_helix(
                synthetic_points,
                initial_axes=[fit.axis],
                max_nfev=1200,
                include_pca_seeds=False,
            )
        except RuntimeError:
            continue
        collected["abs_winding_angle_deg"].append(
            bootstrap_fit.abs_winding_angle_deg
        )
        collected["fitted_pitch_360"].append(bootstrap_fit.fitted_pitch_360)
        endpoint_axial_span = abs(
            float(
                np.dot(
                    synthetic_points[-1] - synthetic_points[0], bootstrap_fit.axis
                )
            )
        )
        collected["endpoint_axial_span"].append(endpoint_axial_span)
        collected["fitted_axial_span"].append(bootstrap_fit.fitted_axial_span)
        collected["radius"].append(bootstrap_fit.radius)
        collected["helix_rms"].append(bootstrap_fit.helix_rms)
        bootstrap_quality = quality_metrics(bootstrap_fit, synthetic_points)
        draws.append(
            {
                "bootstrap_replicate": replicate_index,
                "abs_winding_angle_deg": bootstrap_fit.abs_winding_angle_deg,
                "fitted_pitch_360": bootstrap_fit.fitted_pitch_360,
                "endpoint_axial_span": endpoint_axial_span,
                "fitted_axial_span": bootstrap_fit.fitted_axial_span,
                "start_end_dist": float(
                    np.linalg.norm(synthetic_points[-1] - synthetic_points[0])
                ),
                "fit_radius": bootstrap_fit.radius,
                "radial_rms": bootstrap_fit.radial_rms,
                "axial_rms": bootstrap_fit.axial_rms,
                "helix_rms": bootstrap_fit.helix_rms,
                "axial_angle_r_squared": float(
                    bootstrap_quality["axial_angle_r_squared"]
                ),
                "helix_rms_relative_to_radius": float(
                    bootstrap_quality["helix_rms_relative_to_radius"]
                ),
            }
        )

    successful = len(collected["radius"])
    output: dict[str, float | int] = {
        "bootstrap_requested": replicates,
        "bootstrap_successful": successful,
        "bootstrap_block_length": block_length,
    }
    if successful:
        for metric, values in collected.items():
            array = np.asarray(values, dtype=float)
            output[f"{metric}_ci_low"] = float(np.percentile(array, 2.5))
            output[f"{metric}_ci_high"] = float(np.percentile(array, 97.5))
            output[f"{metric}_bootstrap_sd"] = float(np.std(array, ddof=1)) if successful > 1 else 0.0
    return output, draws


def read_legacy_metrics(path: Path | None) -> dict[str, dict[str, float]]:
    if path is None:
        return {}
    with path.open("r", encoding="utf-8-sig") as handle:
        lines = handle.readlines()
    if lines and lines[0].strip().lower().startswith("sep="):
        delimiter = lines[0].strip()[4:5]
        lines = lines[1:]
    else:
        delimiter = ";" if lines and lines[0].count(";") > lines[0].count(",") else ","
    output: dict[str, dict[str, float]] = {}
    for row in csv.DictReader(lines, delimiter=delimiter):
        specimen_id = (row.get("specimen_id") or "").strip()
        if not specimen_id:
            continue
        values: dict[str, float] = {}
        for key, value in row.items():
            if key == "specimen_id" or value is None or not value.strip():
                continue
            try:
                values[key] = float(value.strip().replace(",", "."))
            except ValueError:
                continue
        output[specimen_id] = values
    return output


def finite_or_blank(value: object) -> object:
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return ""
    return value


def write_records(path: Path, records: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for record in records:
        for key in record:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow({key: finite_or_blank(value) for key, value in record.items()})


def specimen_summary_record(
    trace: LandmarkTrace,
    legacy: LegacyFit,
    robust: RobustHelixFit,
    old_metrics: dict[str, float],
    bootstrap: dict[str, float | int],
) -> dict[str, object]:
    quality = quality_metrics(robust, trace.points)
    old_angle = old_metrics.get("abs_winding_angle_deg", math.nan)
    old_axial_span = old_metrics.get("axial_span", math.nan)
    old_radius = old_metrics.get("fit_radius", math.nan)
    old_rms = old_metrics.get("fit_rms", math.nan)
    old_pitch = (
        old_axial_span * 360.0 / old_angle
        if np.isfinite(old_axial_span) and np.isfinite(old_angle) and old_angle > EPS
        else math.nan
    )
    endpoint_axial_span = abs(float(np.dot(trace.points[-1] - trace.points[0], robust.axis)))
    endpoint_pitch_new_axis = (
        endpoint_axial_span * 360.0 / robust.abs_winding_angle_deg
        if robust.abs_winding_angle_deg > EPS
        else math.nan
    )
    record: dict[str, object] = {
        "specimen_id": trace.specimen_id,
        # Keep released tables relocatable and avoid embedding workstation paths.
        "source_file": trace.source_path.name,
        "n_points": trace.points.shape[0],
        "units": trace.units,
        "robust_fit_success": robust.success,
        "robust_fit_message": robust.message,
        "robust_candidate_count": robust.candidate_count,
        "robust_nfev": robust.nfev,
        "second_best_cost_ratio": robust.second_best_cost_ratio,
        "near_optimal_solution_count": robust.near_optimal_solution_count,
        "near_optimal_angle_range_deg": robust.near_optimal_angle_range_deg,
        "near_optimal_pitch_range": robust.near_optimal_pitch_range,
        "near_optimal_axis_max_separation_deg": robust.near_optimal_axis_max_separation_deg,
        "signed_winding_angle_deg": robust.signed_winding_angle_deg,
        "abs_winding_angle_deg": robust.abs_winding_angle_deg,
        "n_turns_abs": robust.n_turns_abs,
        "fitted_pitch_360": robust.fitted_pitch_360,
        "endpoint_equivalent_pitch_new_axis": endpoint_pitch_new_axis,
        "endpoint_axial_span_new_axis": endpoint_axial_span,
        "fitted_axial_span": robust.fitted_axial_span,
        "fit_radius": robust.radius,
        "radial_rms": robust.radial_rms,
        "axial_rms": robust.axial_rms,
        "helix_rms": robust.helix_rms,
        "axis_x": robust.axis[0],
        "axis_y": robust.axis[1],
        "axis_z": robust.axis[2],
        "axis_point_x": robust.axis_point[0],
        "axis_point_y": robust.axis_point[1],
        "axis_point_z": robust.axis_point[2],
        **quality,
        "analysis_set_all": True,
        "analysis_set_primary_adequate": bool(
            float(quality["helix_rms_relative_to_radius"]) <= 0.10
        ),
        "analysis_set_strict_good": quality["quality_class"] == "good",
        "analysis_set_definition": (
            "primary: helix_rms/fit_radius <= 0.10; "
            "strict: no provisional geometry-quality warnings"
        ),
        "ported_legacy_abs_winding_angle_deg": legacy.abs_winding_angle_deg,
        "ported_legacy_axial_span": legacy.axial_span,
        "ported_legacy_fit_radius": legacy.radius,
        "ported_legacy_radial_rms": legacy.radial_rms,
        "released_abs_winding_angle_deg": old_angle,
        "released_axial_span": old_axial_span,
        "released_endpoint_equivalent_pitch_360": old_pitch,
        "released_fit_radius": old_radius,
        "released_radial_rms": old_rms,
        "ported_minus_released_angle_deg": legacy.abs_winding_angle_deg - old_angle,
        "ported_minus_released_axial_span": legacy.axial_span - old_axial_span,
        "ported_minus_released_radius": legacy.radius - old_radius,
        "ported_minus_released_rms": legacy.radial_rms - old_rms,
        "robust_minus_released_angle_deg": robust.abs_winding_angle_deg - old_angle,
        "robust_minus_released_pitch_360": robust.fitted_pitch_360 - old_pitch,
        "robust_minus_released_radius": robust.radius - old_radius,
        "uncertainty_interpretation": "conditional moving-block residual bootstrap; excludes manual tracing repeatability",
        **bootstrap,
    }
    return record


def residual_records(
    trace: LandmarkTrace, fit: RobustHelixFit
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for index, (point, predicted) in enumerate(zip(trace.points, fit.predicted_points)):
        records.append(
            {
                "specimen_id": trace.specimen_id,
                "point_index": index + 1,
                "point_name": trace.point_names[index],
                "x": point[0],
                "y": point[1],
                "z": point[2],
                "fit_x": predicted[0],
                "fit_y": predicted[1],
                "fit_z": predicted[2],
                "theta_deg": math.degrees(float(fit.theta[index])),
                "radial_coordinate": fit.radial_coordinates[index],
                "axial_coordinate": fit.axial_coordinates[index],
                "radial_residual": fit.radial_residuals[index],
                "axial_residual": fit.axial_residuals[index],
                "helix_distance": fit.point_distances[index],
            }
        )
    return records


def make_specimen_qc_plot(trace: LandmarkTrace, fit: RobustHelixFit, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(15, 5))
    axis_3d = figure.add_subplot(1, 3, 1, projection="3d")
    axis_3d.plot(*trace.points.T, color="#333333", linewidth=1.0, alpha=0.65)
    axis_3d.scatter(*trace.points.T, c=np.arange(trace.points.shape[0]), cmap="viridis", s=24)

    dense_theta = np.linspace(float(np.min(fit.theta)), float(np.max(fit.theta)), 400)
    centered_dense = dense_theta - float(np.mean(fit.theta))
    dense_axial = fit.axial_intercept + fit.axial_rise_per_radian * centered_dense
    dense_points = (
        fit.axis_point
        + np.outer(dense_axial, fit.axis)
        + fit.radius
        * (
            np.outer(np.cos(dense_theta), fit.basis_e1)
            + np.outer(np.sin(dense_theta), fit.basis_e2)
        )
    )
    axis_3d.plot(*dense_points.T, color="#d62728", linewidth=2.0, label="robust helix")
    for point, predicted in zip(trace.points, fit.predicted_points):
        axis_3d.plot(
            [point[0], predicted[0]],
            [point[1], predicted[1]],
            [point[2], predicted[2]],
            color="#999999",
            linewidth=0.5,
        )
    half_span = max(float(np.ptp(fit.axial_coordinates)), fit.radius) * 0.7
    axis_line = np.vstack(
        (fit.axis_point - half_span * fit.axis, fit.axis_point + half_span * fit.axis)
    )
    axis_3d.plot(*axis_line.T, color="#1f77b4", linewidth=2.0, label="fitted axis")
    axis_3d.set_title("3D fit")
    axis_3d.set_xlabel("X (mm)")
    axis_3d.set_ylabel("Y (mm)")
    axis_3d.set_zlabel("Z (mm)")
    axis_3d.legend(loc="best", fontsize=8)

    axis_radial = figure.add_subplot(1, 3, 2)
    theta_degrees = np.degrees(fit.theta - fit.theta[0])
    axis_radial.scatter(theta_degrees, fit.radial_coordinates, color="#333333", s=22)
    axis_radial.axhline(fit.radius, color="#d62728", linewidth=1.5)
    axis_radial.set_xlabel("Unwrapped angle from first point (deg)")
    axis_radial.set_ylabel("Radius (mm)")
    axis_radial.set_title(f"Radial RMS = {fit.radial_rms:.4g} mm")

    axis_axial = figure.add_subplot(1, 3, 3)
    axis_axial.scatter(theta_degrees, fit.axial_coordinates, color="#333333", s=22)
    fitted_axial = fit.axial_intercept + fit.axial_rise_per_radian * (
        fit.theta - float(np.mean(fit.theta))
    )
    axis_axial.plot(theta_degrees, fitted_axial, color="#d62728", linewidth=1.5)
    axis_axial.set_xlabel("Unwrapped angle from first point (deg)")
    axis_axial.set_ylabel("Axial coordinate (mm)")
    axis_axial.set_title(f"Pitch/360 = {fit.fitted_pitch_360:.4g} mm")

    warning_text = ", ".join(fit.warnings) if fit.warnings else "none"
    figure.suptitle(
        f"{trace.specimen_id} | n={trace.points.shape[0]} | "
        f"angle={fit.abs_winding_angle_deg:.2f} deg | warnings: {warning_text}",
        fontsize=10,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def annotate_outliers(axis: object, x: np.ndarray, y: np.ndarray, names: Sequence[str], count: int = 6) -> None:
    finite = np.isfinite(x) & np.isfinite(y)
    if not np.any(finite):
        return
    x_finite = x[finite]
    y_finite = y[finite]
    names_finite = np.asarray(names)[finite]
    if len(x_finite) <= count:
        indices = np.arange(len(x_finite))
    else:
        x_z = (x_finite - np.mean(x_finite)) / (np.std(x_finite) + EPS)
        y_z = (y_finite - np.mean(y_finite)) / (np.std(y_finite) + EPS)
        indices = np.argsort(np.square(x_z) + np.square(y_z))[-count:]
    for index in indices:
        axis.annotate(
            names_finite[index],
            (x_finite[index], y_finite[index]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=6,
        )


def make_summary_qc_plot(records: Sequence[dict[str, object]], output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_path.parent.mkdir(parents=True, exist_ok=True)
    names = [str(record["specimen_id"]) for record in records]

    def values(key: str) -> np.ndarray:
        return np.asarray([float(record.get(key, math.nan)) for record in records])

    figure, axes = plt.subplots(2, 2, figsize=(13, 11))
    plots = [
        (
            "released_abs_winding_angle_deg",
            "abs_winding_angle_deg",
            "Released angle (deg)",
            "Robust fitted angle (deg)",
        ),
        (
            "released_endpoint_equivalent_pitch_360",
            "fitted_pitch_360",
            "Released endpoint-equivalent pitch (mm)",
            "Robust fitted pitch (mm)",
        ),
        (
            "abs_winding_angle_deg",
            "axial_angle_r_squared",
            "Robust fitted angle (deg)",
            "Axial-angle R-squared",
        ),
        (
            "abs_winding_angle_deg",
            "helix_rms_relative_to_radius",
            "Robust fitted angle (deg)",
            "Helix RMS / radius",
        ),
    ]
    for axis, (x_key, y_key, x_label, y_label) in zip(axes.flat, plots):
        x = values(x_key)
        y = values(y_key)
        finite = np.isfinite(x) & np.isfinite(y)
        axis.scatter(x[finite], y[finite], s=28, color="#2c7fb8", alpha=0.8)
        if x_key.startswith("released") and np.any(finite):
            limits = [min(float(np.min(x[finite])), float(np.min(y[finite]))), max(float(np.max(x[finite])), float(np.max(y[finite])))]
            axis.plot(limits, limits, linestyle="--", color="#555555", linewidth=1.0)
        if "pitch" in x_key and np.any(finite & (x > 0.0) & (y > 0.0)):
            axis.set_xscale("log")
            axis.set_yscale("log")
        annotate_outliers(axis, x, y, names)
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.grid(alpha=0.2)
    figure.suptitle("Robust 3D helix fit: comparison and quality control")
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing landmark CSV files")
    parser.add_argument("output_dir", type=Path, help="New directory for validation outputs")
    parser.add_argument(
        "--released-metrics",
        type=Path,
        default=None,
        help="Optional released winding-metrics CSV for reproduction checks",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help="Conditional moving-block residual-bootstrap replicates per specimen",
    )
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--no-specimen-plots", action="store_true")
    parser.add_argument("--pattern", default="*.csv")
    return parser.parse_args(argv)


def run_batch(arguments: argparse.Namespace) -> int:
    input_dir = arguments.input_dir.resolve()
    output_dir = arguments.output_dir.resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    if output_dir == input_dir or input_dir in output_dir.parents:
        raise ValueError("Output directory must not be the input directory or one of its children")
    output_dir.mkdir(parents=True, exist_ok=True)

    released_path = arguments.released_metrics
    if released_path is None:
        candidate = input_dir / "winding_metrics_excelDE.csv"
        released_path = candidate if candidate.exists() else None
    released_metrics = read_legacy_metrics(released_path)

    input_files = sorted(
        path
        for path in input_dir.glob(arguments.pattern)
        if path.is_file() and "winding_metrics" not in path.name.lower()
    )
    if not input_files:
        raise FileNotFoundError(f"No landmark CSV files found in {input_dir}")

    summaries: list[dict[str, object]] = []
    residuals: list[dict[str, object]] = []
    bootstrap_draws: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    file_manifest: list[dict[str, object]] = []
    started = time.time()

    for file_index, path in enumerate(input_files, start=1):
        print(f"[{file_index:02d}/{len(input_files):02d}] {path.stem}", flush=True)
        file_manifest.append(
            {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
        try:
            trace = read_landmark_csv(path)
            legacy = fit_legacy_cinema_method(trace.points)
            robust = fit_robust_helix(trace.points, initial_axes=[legacy.axis])
            bootstrap, specimen_bootstrap_draws = bootstrap_fit_intervals(
                trace.points,
                robust,
                arguments.bootstrap,
                arguments.seed + file_index * 1009,
            )
            for draw in specimen_bootstrap_draws:
                bootstrap_draws.append(
                    {
                        "specimen_id": trace.specimen_id,
                        "specimen_seed": arguments.seed + file_index * 1009,
                        **draw,
                    }
                )
            summary = specimen_summary_record(
                trace,
                legacy,
                robust,
                released_metrics.get(trace.specimen_id, {}),
                bootstrap,
            )
            summaries.append(summary)
            residuals.extend(residual_records(trace, robust))
            if not arguments.no_specimen_plots:
                make_specimen_qc_plot(
                    trace, robust, output_dir / "specimen_qc" / f"{trace.specimen_id}.png"
                )
            print(
                f"    robust angle={robust.abs_winding_angle_deg:.2f} deg, "
                f"pitch={robust.fitted_pitch_360:.5f} mm, "
                f"RMS={robust.helix_rms:.5f} mm, "
                f"quality={summary['quality_class']}",
                flush=True,
            )
        except Exception as error:  # Continue batch while preserving the exact failure.
            failures.append(
                {
                    "specimen_id": path.stem,
                    "source_file": path.name,
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
            print(f"    FAILED: {type(error).__name__}: {error}", flush=True)

    write_records(output_dir / "robust_helix_metrics.csv", summaries)
    write_records(output_dir / "robust_helix_point_residuals.csv", residuals)
    write_records(output_dir / "robust_helix_bootstrap_draws.csv", bootstrap_draws)
    write_records(output_dir / "fit_failures.csv", failures)
    if summaries:
        make_summary_qc_plot(summaries, output_dir / "robust_helix_summary_qc.png")

    reproduction_fields = (
        "ported_minus_released_angle_deg",
        "ported_minus_released_axial_span",
        "ported_minus_released_radius",
        "ported_minus_released_rms",
    )
    reproduction_maxima: dict[str, float] = {}
    for field_name in reproduction_fields:
        finite_values = [
            abs(float(record[field_name]))
            for record in summaries
            if np.isfinite(float(record.get(field_name, math.nan)))
        ]
        reproduction_maxima[f"max_abs_{field_name}"] = max(finite_values) if finite_values else math.nan

    quality_counts: dict[str, int] = {}
    for record in summaries:
        quality = str(record["quality_class"])
        quality_counts[quality] = quality_counts.get(quality, 0) + 1
    analysis_set_counts = {
        "all": len(summaries),
        "primary_adequate": sum(
            bool(record["analysis_set_primary_adequate"]) for record in summaries
        ),
        "strict_good": sum(
            bool(record["analysis_set_strict_good"]) for record in summaries
        ),
    }
    manifest = {
        "algorithm": "robust circular 3D helix fit",
        "algorithm_version": ALGORITHM_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_directory": str(input_dir),
        "output_directory": str(output_dir),
        "released_metrics": released_path.name if released_path else None,
        "coordinate_units": "read from each Cinema 4D export; expected mm",
        "n_input_files": len(input_files),
        "n_successful": len(summaries),
        "n_failed": len(failures),
        "bootstrap_replicates_per_specimen": arguments.bootstrap,
        "random_seed": arguments.seed,
        "robust_loss": "soft_l1",
        "robust_f_scale_normalized": 0.02,
        "uncertainty_scope": "conditional moving-block residual bootstrap; does not include manual tracing repeatability",
        "quality_flags_are": "provisional diagnostics, not preregistered exclusion criteria",
        "quality_class_counts": quality_counts,
        "analysis_set_counts": analysis_set_counts,
        "analysis_set_definition": {
            "all": "all successful geometry fits",
            "primary_adequate": "helix RMS divided by fitted radius <= 0.10",
            "strict_good": "no provisional geometry-quality warnings",
        },
        "n_bootstrap_draw_rows": len(bootstrap_draws),
        "legacy_reproduction_max_absolute_differences": reproduction_maxima,
        "runtime_seconds": time.time() - started,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy_version,
        "input_files": file_manifest,
    }
    with (output_dir / "analysis_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    print(f"\nSuccessful fits: {len(summaries)}/{len(input_files)}")
    print(f"Failures: {len(failures)}")
    print(f"Outputs: {output_dir}")
    return 0 if not failures else 2


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    return run_batch(arguments)


if __name__ == "__main__":
    raise SystemExit(main())
