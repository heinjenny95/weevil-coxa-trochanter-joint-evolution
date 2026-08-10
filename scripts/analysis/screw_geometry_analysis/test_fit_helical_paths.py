import csv
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from fit_helical_paths import (
    fit_robust_helix,
    normalize,
    orthonormal_basis,
    read_landmark_csv,
)


def synthetic_helix(noise_sd=0.0, outlier=False, seed=1234):
    rng = np.random.default_rng(seed)
    axis = normalize(np.array([0.31, -0.42, 0.852]))
    e1, e2 = orthonormal_basis(axis)
    center = np.array([1.2, -0.7, 0.4])
    theta = np.linspace(0.25, 2.65 * math.pi, 48)
    radius = 0.37
    rise_per_radian = 0.052
    points = (
        center
        + np.outer(rise_per_radian * (theta - np.mean(theta)), axis)
        + radius
        * (
            np.outer(np.cos(theta), e1)
            + np.outer(np.sin(theta), e2)
        )
    )
    if noise_sd:
        points = points + rng.normal(scale=noise_sd, size=points.shape)
    if outlier:
        points[19] += np.array([0.08, -0.05, 0.06])
    return points, axis, radius, rise_per_radian, theta


class RobustHelixFitTests(unittest.TestCase):
    def test_clean_helix_parameter_recovery(self):
        points, axis, radius, rise, theta = synthetic_helix(noise_sd=0.0005)
        fit = fit_robust_helix(points)
        self.assertGreater(abs(float(np.dot(fit.axis, axis))), 0.999)
        self.assertAlmostEqual(fit.radius, radius, delta=0.002)
        self.assertAlmostEqual(
            fit.fitted_pitch_360, 2.0 * math.pi * rise, delta=0.004
        )
        self.assertAlmostEqual(
            fit.abs_winding_angle_deg,
            math.degrees(theta[-1] - theta[0]),
            delta=1.0,
        )
        self.assertLess(fit.helix_rms, 0.002)

    def test_soft_l1_fit_resists_single_outlier(self):
        points, _, radius, rise, _ = synthetic_helix(
            noise_sd=0.0005, outlier=True
        )
        fit = fit_robust_helix(points)
        self.assertAlmostEqual(fit.radius, radius, delta=0.01)
        self.assertAlmostEqual(
            fit.fitted_pitch_360, 2.0 * math.pi * rise, delta=0.02
        )

    def test_cinema_landmark_csv_parser(self):
        points, _, _, _, _ = synthetic_helix()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "example.csv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["Landmarks", len(points), "v1.0", "(mm)"])
                writer.writerow(["#", "Name", "Type", "X", "Y", "Z", "NX", "NY", "NZ"])
                for index, point in enumerate(points, start=1):
                    writer.writerow(
                        [
                            index,
                            f"C.1.{index:02d}",
                            "normal",
                            *point,
                            0.0,
                            0.0,
                            1.0,
                        ]
                    )
            trace = read_landmark_csv(path)
        self.assertEqual(trace.units, "mm")
        self.assertEqual(trace.declared_count, len(points))
        self.assertEqual(trace.point_names[0], "C.1.01")
        np.testing.assert_allclose(trace.points, points)
        self.assertEqual(trace.normals.shape, points.shape)


if __name__ == "__main__":
    unittest.main()
