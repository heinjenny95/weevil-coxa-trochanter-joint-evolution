from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


WIDTHS = (0.025, 0.075, 0.1, 0.15, 0.2)
SEED = 20260825


def normalize_id(value: str) -> str:
    value = value.strip().lower()
    for suffix in ("_trochanter_mirrored_aligned", "_trochanter_aligned"):
        if value.endswith(suffix):
            value = value[: -len(suffix)]
    # Historical filename typo; these identify the same specimen.
    return (
        value.replace("t_pseudonastus", "t_pseudonasutus")
        .replace("rynchites_cupreus", "rhynchites_cupreus")
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Quantify agreement among archived Deformetrica kernel-width runs. "
            "These runs were made during workflow development with 67 specimens "
            "and are not final-dataset sensitivity reruns."
        )
    )
    parser.add_argument(
        "--kernel-root",
        type=Path,
        required=True,
        help="Directory containing the 0.025, 0.075, 0.1, 0.15 and 0.2 run folders.",
    )
    parser.add_argument(
        "--final-scores",
        type=Path,
        required=True,
        help="Final 68-specimen PCA score table used only to audit specimen coverage.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def read_scores(kernel_root: Path, width: float) -> tuple[list[str], np.ndarray]:
    path = kernel_root / str(width) / "pca_scores_with_specimen_ids_DE.csv"
    df = pd.read_csv(path, sep=";", decimal=",")
    ids = [normalize_id(v) for v in df["specimen_id"].astype(str)]
    pc_cols = sorted(
        (c for c in df.columns if c.startswith("PC")),
        key=lambda c: int(c[2:]),
    )
    x = df[pc_cols].to_numpy(dtype=float)
    nonzero = np.nanstd(x, axis=0) > 1e-14
    return ids, x[:, nonzero]


def condensed(a: np.ndarray) -> np.ndarray:
    iu = np.triu_indices(a.shape[0], 1)
    return a[iu]


def distances(x: np.ndarray) -> np.ndarray:
    delta = x[:, None, :] - x[None, :, :]
    return np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))


def ranks(v: np.ndarray) -> np.ndarray:
    return pd.Series(v).rank(method="average").to_numpy(dtype=float)


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    aa = a - a.mean()
    bb = b - b.mean()
    denom = np.linalg.norm(aa) * np.linalg.norm(bb)
    return float(np.dot(aa, bb) / denom)


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return pearson(ranks(a), ranks(b))


def procrustes_similarity(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Return correlation-like similarity and squared disparity after rotation."""
    a = a - a.mean(axis=0, keepdims=True)
    b = b - b.mean(axis=0, keepdims=True)
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    u, singular_values, vt = np.linalg.svd(b.T @ a, full_matrices=False)
    rotation = u @ vt
    aligned = b @ rotation
    disparity = float(np.square(a - aligned).sum())
    similarity = float(singular_values.sum())
    return similarity, disparity


def knn_overlap(da: np.ndarray, db: np.ndarray, k: int) -> float:
    aa = np.argsort(da, axis=1)[:, 1 : k + 1]
    bb = np.argsort(db, axis=1)[:, 1 : k + 1]
    values = [len(set(x).intersection(y)) / k for x, y in zip(aa, bb)]
    return float(np.mean(values))


def mantel_permutation_p(
    da: np.ndarray,
    db: np.ndarray,
    observed: float,
    rng: np.random.Generator,
    permutations: int = 999,
) -> float:
    n = da.shape[0]
    a = ranks(condensed(da))
    exceed = 0
    for _ in range(permutations):
        idx = rng.permutation(n)
        permuted = db[np.ix_(idx, idx)]
        stat = pearson(a, ranks(condensed(permuted)))
        exceed += stat >= observed
    return (exceed + 1) / (permutations + 1)


def subsample_interval(
    da: np.ndarray,
    db: np.ndarray,
    rng: np.random.Generator,
    fraction: float = 0.8,
    repetitions: int = 1000,
) -> tuple[float, float, float]:
    n = da.shape[0]
    m = round(n * fraction)
    values = np.empty(repetitions, dtype=float)
    for i in range(repetitions):
        idx = np.sort(rng.choice(n, size=m, replace=False))
        values[i] = spearman(
            condensed(da[np.ix_(idx, idx)]),
            condensed(db[np.ix_(idx, idx)]),
        )
    q = np.quantile(values, [0.025, 0.5, 0.975])
    return tuple(float(v) for v in q)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    data: dict[float, dict[str, object]] = {}
    reference_ids: list[str] | None = None
    for width in WIDTHS:
        ids, x = read_scores(args.kernel_root, width)
        if reference_ids is None:
            reference_ids = ids
        if set(ids) != set(reference_ids):
            raise RuntimeError(f"Specimen mismatch at kernel width {width}")
        order = [ids.index(v) for v in reference_ids]
        x = x[order]
        data[width] = {
            "ids": reference_ids,
            "x": x,
            "d_full": distances(x),
            "d_pc5": distances(x[:, :5]),
        }

    current_df = pd.read_csv(args.final_scores)
    current_ids = {normalize_id(v) for v in current_df["specimen_id"].astype(str)}
    historical_ids = set(reference_ids or [])

    rows: list[dict[str, object]] = []
    for i, wa in enumerate(WIDTHS):
        for wb in WIDTHS[i + 1 :]:
            a = data[wa]
            b = data[wb]
            dfa = a["d_full"]
            dfb = b["d_full"]
            d5a = a["d_pc5"]
            d5b = b["d_pc5"]
            assert isinstance(dfa, np.ndarray)
            assert isinstance(dfb, np.ndarray)
            assert isinstance(d5a, np.ndarray)
            assert isinstance(d5b, np.ndarray)
            full_rho = spearman(condensed(dfa), condensed(dfb))
            pc5_rho = spearman(condensed(d5a), condensed(d5b))
            x_a = a["x"]
            x_b = b["x"]
            assert isinstance(x_a, np.ndarray)
            assert isinstance(x_b, np.ndarray)
            proc_similarity, proc_disparity = procrustes_similarity(x_a[:, :5], x_b[:, :5])
            lo, med, hi = subsample_interval(dfa, dfb, rng)
            rows.append(
                {
                    "kernel_width_a": wa,
                    "kernel_width_b": wb,
                    "n_specimens": len(reference_ids or []),
                    "full_space_distance_spearman_rho": full_rho,
                    "full_space_mantel_p_999": mantel_permutation_p(
                        dfa, dfb, full_rho, rng
                    ),
                    "full_space_rho_80pct_subsample_q025": lo,
                    "full_space_rho_80pct_subsample_median": med,
                    "full_space_rho_80pct_subsample_q975": hi,
                    "pc1_pc5_distance_spearman_rho": pc5_rho,
                    "pc1_pc5_procrustes_similarity": proc_similarity,
                    "pc1_pc5_procrustes_disparity": proc_disparity,
                    "knn_overlap_k5": knn_overlap(dfa, dfb, 5),
                    "knn_overlap_k10": knn_overlap(dfa, dfb, 10),
                }
            )

    result = pd.DataFrame(rows)
    result.to_csv(
        args.output_dir / "kernel_width_pairwise_stability.csv", index=False
    )

    control_points: list[dict[str, object]] = []
    for width in WIDTHS:
        header = (args.kernel_root / str(width) / "Atlas_Momentas.txt").read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()[0]
        subjects, points, dimensions = [int(v) for v in header.split()[:3]]
        control_points.append(
            {
                "kernel_width": width,
                "n_subjects": subjects,
                "n_control_points": points,
                "dimensions": dimensions,
                "workflow_stage": "parameter selection before final atlas",
                "included_in_final_atlas": False,
            }
        )
    run_summary = pd.DataFrame(control_points)
    run_summary.to_csv(args.output_dir / "kernel_width_run_summary.csv", index=False)

    audit = {
        "seed": args.seed,
        "historical_n": len(historical_ids),
        "current_n": len(current_ids),
        "missing_from_historical_sensitivity": sorted(current_ids - historical_ids),
        "absent_from_current_final": sorted(historical_ids - current_ids),
        "workflow_stage": "parameter selection before the final 68-specimen atlas",
        "selected_width": 0.1,
        "final_atlas_policy": (
            "Kernel width 0.1 was retained without further tuning when Rhynchites "
            "cupreus was added to the final 68-specimen atlas."
        ),
        "interpretation_limit": (
            "The five archived runs quantify stability of the development-stage "
            "67-specimen ordinations. They document parameter selection but are not "
            "sensitivity reruns of the final 68-specimen atlas and do not establish "
            "that width 0.1 is uniquely optimal."
        ),
    }
    (args.output_dir / "kernel_width_dataset_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )

    focus = result[
        (result["kernel_width_a"] == 0.1) | (result["kernel_width_b"] == 0.1)
    ]
    print("DATASET AUDIT")
    print(json.dumps(audit, indent=2))
    print("\nCOMPARISONS INVOLVING WIDTH 0.1")
    print(focus.to_string(index=False))
    print("\nADJACENT-WIDTH COMPARISONS")
    adjacent = result[
        np.isclose(
            result["kernel_width_b"] - result["kernel_width_a"],
            [0.05 if a == 0.025 else 0.025 if a == 0.075 else 0.05 for a in result["kernel_width_a"]],
        )
    ]
    # Explicit pairs avoid treating the uneven width grid as a regular sequence.
    adjacent = result[
        result.apply(
            lambda r: (r["kernel_width_a"], r["kernel_width_b"])
            in {(0.025, 0.075), (0.075, 0.1), (0.1, 0.15), (0.15, 0.2)},
            axis=1,
        )
    ]
    print(adjacent.to_string(index=False))


if __name__ == "__main__":
    main()
