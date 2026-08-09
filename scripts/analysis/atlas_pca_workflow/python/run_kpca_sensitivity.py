from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f as f_distribution
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import KernelPCA


SEED = 20260803
N_PERM = 9999
N_BOOT = 5000

OUT = Path("kPCA_sensitivity_output")
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
LOGS = OUT / "logs"

MOMENTA = Path("Atlas_Momentas.txt")
LINEAR_SCORES = Path("PCA_scores_with_specimen_id.csv")
ANALYSIS_DATA = Path("PCA_scores_with_specimen_id_with_centroid_size.csv")
JOINT_DATA = Path("specimen_joint_types.csv")

FAMILY_ORDER = [
    "Anthribidae",
    "Attelabidae",
    "Belidae",
    "Brentidae",
    "Caridae",
    "Curculionidae",
    "Nemonychidae",
]
FAMILY_COLORS = {
    "Anthribidae": "#7E5AA6",
    "Attelabidae": "#D96B5F",
    "Belidae": "#7FA2D6",
    "Brentidae": "#C2A532",
    "Caridae": "#E58A2E",
    "Curculionidae": "#55A88F",
    "Nemonychidae": "#78B84B",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare ordinary PCA with two RBF-kPCA parameterizations."
    )
    parser.add_argument("--momenta", type=Path, required=True, help="Deformetrica Atlas_Momentas.txt file.")
    parser.add_argument("--linear-scores", type=Path, required=True, help="Released ordinary-PCA score table.")
    parser.add_argument("--analysis-data", type=Path, required=True, help="Specimen-level PCA and trait table.")
    parser.add_argument("--joint-data", type=Path, required=True, help="Specimen-level joint-type table.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for tables, figures and logs.")
    parser.add_argument(
        "--geometry-min-angle",
        type=float,
        default=30.0,
        help="Minimum absolute winding angle (degrees) for all shape--geometry sensitivity models.",
    )
    return parser.parse_args()


def configure_paths(args: argparse.Namespace) -> None:
    global OUT, TABLES, FIGURES, LOGS, MOMENTA, LINEAR_SCORES, ANALYSIS_DATA, JOINT_DATA, GEOMETRY_MIN_ANGLE
    OUT = args.output_dir.resolve()
    TABLES = OUT / "tables"
    FIGURES = OUT / "figures"
    LOGS = OUT / "logs"
    MOMENTA = args.momenta.resolve()
    LINEAR_SCORES = args.linear_scores.resolve()
    ANALYSIS_DATA = args.analysis_data.resolve()
    JOINT_DATA = args.joint_data.resolve()
    GEOMETRY_MIN_ANGLE = float(args.geometry_min_angle)


def ensure_dirs() -> None:
    for directory in (TABLES, FIGURES, LOGS):
        directory.mkdir(parents=True, exist_ok=True)


def read_semicolon(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=";", dtype=str)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.replace(",", ".", regex=False), errors="coerce")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_momenta(path: Path) -> np.ndarray:
    with path.open("r", encoding="utf-8") as handle:
        header = [int(value) for value in handle.readline().split()]
    if len(header) != 3:
        raise ValueError(f"Unexpected momenta header: {header}")
    n_subjects, n_control_points, n_dimensions = header
    values = np.loadtxt(path, skiprows=1)
    expected = n_subjects * n_control_points * n_dimensions
    if values.size != expected:
        raise ValueError(f"Expected {expected} momenta values, found {values.size}")
    return values.reshape(n_subjects, n_control_points, n_dimensions).reshape(n_subjects, -1)


def load_inputs() -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    x = load_momenta(MOMENTA)

    score_df = read_semicolon(LINEAR_SCORES)
    pc_columns = sorted(
        [column for column in score_df.columns if column.startswith("PC")],
        key=lambda value: int(value[2:]),
    )
    for column in pc_columns:
        score_df[column] = numeric(score_df[column])

    analysis_df = read_semicolon(ANALYSIS_DATA)
    for column in analysis_df.columns:
        if column.startswith("PC") or column in {
            "centroid_size",
            "abs_winding_angle_deg",
            "n_turns_abs",
            "axial_span",
            "fit_radius",
            "fit_rms",
        }:
            analysis_df[column] = numeric(analysis_df[column])

    joint_df = read_semicolon(JOINT_DATA)

    if len(score_df) != x.shape[0]:
        raise ValueError("The momenta and final linear-score table have different specimen counts")
    if analysis_df["specimen_id"].duplicated().any():
        raise ValueError("Duplicate specimen IDs in the main analysis table")
    if score_df["specimen_id"].tolist() != analysis_df["specimen_id"].tolist():
        raise ValueError("Specimen order differs between final score and analysis tables")

    return x, score_df, analysis_df, joint_df, pc_columns


def fit_kpca(x: np.ndarray, gamma: float) -> tuple[np.ndarray, np.ndarray, KernelPCA]:
    model = KernelPCA(
        n_components=x.shape[0] - 1,
        kernel="rbf",
        gamma=gamma,
        eigen_solver="dense",
        remove_zero_eig=True,
    )
    scores = model.fit_transform(x)
    eigenvalues = np.asarray(model.eigenvalues_, dtype=float)
    positive = eigenvalues > np.finfo(float).eps * max(1.0, eigenvalues.max())
    return scores[:, positive], eigenvalues[positive], model


def align_axes(
    linear: np.ndarray,
    kernel: np.ndarray,
    n_target: int = 5,
    n_candidates: int = 12,
) -> tuple[np.ndarray, pd.DataFrame]:
    n_candidates = min(n_candidates, kernel.shape[1])
    correlations = np.empty((n_target, n_candidates), dtype=float)
    for i in range(n_target):
        for j in range(n_candidates):
            correlations[i, j] = pearsonr(linear[:, i], kernel[:, j]).statistic

    rows, cols = linear_sum_assignment(-np.abs(correlations))
    assignment = {row: col for row, col in zip(rows, cols)}
    aligned = np.empty((kernel.shape[0], n_target), dtype=float)
    records = []
    for i in range(n_target):
        j = assignment[i]
        sign = 1.0 if correlations[i, j] >= 0 else -1.0
        aligned[:, i] = kernel[:, j] * sign
        records.append(
            {
                "linear_axis": f"PC{i + 1}",
                "matched_native_kernel_axis": f"kPC{j + 1}",
                "sign_multiplier": int(sign),
                "pearson_r": correlations[i, j],
                "abs_pearson_r": abs(correlations[i, j]),
            }
        )
    return aligned, pd.DataFrame(records)


def pseudo_f_group(y: np.ndarray, groups: np.ndarray) -> tuple[float, float, float, int, int]:
    groups = np.asarray(groups)
    levels = pd.unique(groups)
    grand = y.mean(axis=0)
    ss_total = float(np.sum((y - grand) ** 2))
    ss_within = 0.0
    for level in levels:
        block = y[groups == level]
        ss_within += float(np.sum((block - block.mean(axis=0)) ** 2))
    ss_between = ss_total - ss_within
    df_between = len(levels) - 1
    df_within = len(y) - len(levels)
    statistic = (ss_between / df_between) / (ss_within / df_within)
    r2 = ss_between / ss_total
    return statistic, r2, ss_within, df_between, df_within


def permutation_group_test(
    y: np.ndarray,
    groups: np.ndarray,
    rng: np.random.Generator,
    n_perm: int = N_PERM,
) -> dict[str, float]:
    observed, r2, _, df1, df2 = pseudo_f_group(y, groups)
    exceed = 0
    for _ in range(n_perm):
        permuted = rng.permutation(groups)
        test_statistic, _, _, _, _ = pseudo_f_group(y, permuted)
        exceed += test_statistic >= observed - 1e-12
    return {
        "pseudo_F": observed,
        "R2": r2,
        "df1": df1,
        "df2": df2,
        "permutation_P": (exceed + 1) / (n_perm + 1),
        "permutations": n_perm,
    }


def multivariate_regression(y: np.ndarray, x: np.ndarray) -> tuple[float, float, int, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    design = np.column_stack([np.ones(len(x)), x])
    fitted = design @ np.linalg.lstsq(design, y, rcond=None)[0]
    grand = y.mean(axis=0)
    ss_total = float(np.sum((y - grand) ** 2))
    ss_error = float(np.sum((y - fitted) ** 2))
    ss_model = ss_total - ss_error
    df_model = design.shape[1] - 1
    df_error = len(y) - design.shape[1]
    pseudo_f = (ss_model / df_model) / (ss_error / df_error)
    return ss_model / ss_total, pseudo_f, df_model, df_error


def permutation_regression(
    y: np.ndarray,
    x: np.ndarray,
    rng: np.random.Generator,
    n_perm: int = N_PERM,
) -> dict[str, float]:
    r2, observed, df1, df2 = multivariate_regression(y, x)
    exceed = 0
    for _ in range(n_perm):
        _, statistic, _, _ = multivariate_regression(y, rng.permutation(x))
        exceed += statistic >= observed - 1e-12
    return {
        "R2": r2,
        "pseudo_F": observed,
        "df1": df1,
        "df2": df2,
        "permutation_P": (exceed + 1) / (n_perm + 1),
        "permutations": n_perm,
    }


def ols_summary(y: np.ndarray, x: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
    y = y[mask]
    x = x[mask]
    design = np.column_stack([np.ones(len(x)), x])
    fitted = design @ np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - fitted
    ss_total = float(np.sum((y - y.mean()) ** 2))
    ss_error = float(np.sum(residual**2))
    r2 = 1.0 - ss_error / ss_total
    df1 = x.shape[1]
    df2 = len(y) - design.shape[1]
    f_value = ((ss_total - ss_error) / df1) / (ss_error / df2)
    return {
        "n": len(y),
        "R2": r2,
        "F": f_value,
        "df1": df1,
        "df2": df2,
        "P": float(f_distribution.sf(f_value, df1, df2)),
    }


def family_disparity(y: np.ndarray, groups: np.ndarray) -> pd.DataFrame:
    records = []
    for family in FAMILY_ORDER:
        block = y[groups == family]
        if len(block) < 3:
            continue
        squared = np.sum((block - block.mean(axis=0)) ** 2, axis=1)
        records.append(
            {
                "Family": family,
                "n": len(block),
                "mean_squared_distance": float(squared.mean()),
                "median_squared_distance": float(np.median(squared)),
            }
        )
    return pd.DataFrame(records)


def dispersion_test(
    y: np.ndarray,
    groups: np.ndarray,
    rng: np.random.Generator,
    n_perm: int = N_PERM,
) -> dict[str, float]:
    def distances(labels: np.ndarray) -> np.ndarray:
        result = np.empty(len(labels), dtype=float)
        for level in pd.unique(labels):
            positions = labels == level
            centre = y[positions].mean(axis=0)
            result[positions] = np.linalg.norm(y[positions] - centre, axis=1)
        return result

    def one_way(values: np.ndarray, labels: np.ndarray) -> float:
        levels = pd.unique(labels)
        grand = values.mean()
        ss_between = sum(np.sum(labels == level) * (values[labels == level].mean() - grand) ** 2 for level in levels)
        ss_within = sum(np.sum((values[labels == level] - values[labels == level].mean()) ** 2) for level in levels)
        return (ss_between / (len(levels) - 1)) / (ss_within / (len(values) - len(levels)))

    observed_distances = distances(groups)
    observed = one_way(observed_distances, groups)
    exceed = 0
    for _ in range(n_perm):
        permuted = rng.permutation(groups)
        statistic = one_way(distances(permuted), permuted)
        exceed += statistic >= observed - 1e-12
    return {
        "F": observed,
        "df1": len(pd.unique(groups)) - 1,
        "df2": len(groups) - len(pd.unique(groups)),
        "permutation_P": (exceed + 1) / (n_perm + 1),
        "permutations": n_perm,
    }


def rarefied_disparity(
    y: np.ndarray,
    groups: np.ndarray,
    rng: np.random.Generator,
    sample_size: int = 4,
    n_boot: int = N_BOOT,
) -> pd.DataFrame:
    records = []
    for family in FAMILY_ORDER:
        block = y[groups == family]
        if len(block) < sample_size:
            continue
        estimates = []
        for _ in range(n_boot):
            indices = rng.choice(len(block), size=sample_size, replace=False)
            sampled = block[indices]
            estimates.append(np.mean(np.sum((sampled - sampled.mean(axis=0)) ** 2, axis=1)))
        records.append(
            {
                "Family": family,
                "original_n": len(block),
                "rarefied_n": sample_size,
                "mean": float(np.mean(estimates)),
                "lower_95": float(np.quantile(estimates, 0.025)),
                "upper_95": float(np.quantile(estimates, 0.975)),
                "resamples": n_boot,
            }
        )
    return pd.DataFrame(records)


def bh_adjust(values: pd.Series) -> pd.Series:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return pd.Series(output, index=values.index)


def pairwise_disparity_tests(
    y: np.ndarray,
    groups: np.ndarray,
    rng: np.random.Generator,
    n_perm: int = N_PERM,
) -> pd.DataFrame:
    eligible = [family for family in FAMILY_ORDER if np.sum(groups == family) >= 3]
    records = []
    for index, first in enumerate(eligible):
        for second in eligible[index + 1 :]:
            mask = np.isin(groups, [first, second])
            local_y = y[mask]
            local_groups = groups[mask]

            def disparity_difference(labels: np.ndarray) -> float:
                estimates = []
                for level in (first, second):
                    block = local_y[labels == level]
                    estimates.append(np.mean(np.sum((block - block.mean(axis=0)) ** 2, axis=1)))
                return estimates[0] - estimates[1]

            observed = disparity_difference(local_groups)
            exceed = 0
            for _ in range(n_perm):
                permuted = rng.permutation(local_groups)
                exceed += abs(disparity_difference(permuted)) >= abs(observed) - 1e-12
            records.append(
                {
                    "family_1": first,
                    "family_2": second,
                    "difference_family1_minus_family2": observed,
                    "permutation_P": (exceed + 1) / (n_perm + 1),
                    "permutations": n_perm,
                }
            )
    result = pd.DataFrame(records)
    if not result.empty:
        result["FDR_P"] = bh_adjust(result["permutation_P"])
    return result


def write_scores(
    specimen_ids: pd.Series,
    family: pd.Series,
    scores: np.ndarray,
    prefix: str,
    path: Path,
) -> None:
    frame = pd.DataFrame(scores, columns=[f"{prefix}{i + 1}" for i in range(scores.shape[1])])
    frame.insert(0, "Family", family.to_numpy())
    frame.insert(0, "specimen_id", specimen_ids.to_numpy())
    frame.to_csv(path, index=False)


def configure_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 6.5,
            "axes.labelsize": 6.5,
            "axes.titlesize": 7,
            "axes.linewidth": 0.6,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 5.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    sns.set_style("whitegrid", {"grid.color": "#E3E9ED", "grid.linewidth": 0.5})


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def scatter_panel(
    ax: plt.Axes,
    values: np.ndarray,
    families: np.ndarray,
    x_index: int,
    y_index: int,
    title: str,
    x_label: str,
    y_label: str,
) -> None:
    for family in FAMILY_ORDER:
        mask = families == family
        ax.scatter(
            values[mask, x_index],
            values[mask, y_index],
            s=13,
            color=FAMILY_COLORS[family],
            edgecolor="white",
            linewidth=0.35,
            alpha=0.92,
            label=family,
        )
    ax.set_title(title, pad=3, fontweight="bold")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def make_figures(
    linear: np.ndarray,
    legacy_aligned: np.ndarray,
    adaptive_aligned: np.ndarray,
    families: np.ndarray,
    linear_fraction: np.ndarray,
    legacy_fraction: np.ndarray,
    adaptive_fraction: np.ndarray,
    corr_legacy: np.ndarray,
    corr_adaptive: np.ndarray,
    disparity_long: pd.DataFrame,
) -> None:
    configure_plotting()

    fig, axes = plt.subplots(3, 2, figsize=(7.0866, 7.10))
    fig.subplots_adjust(
        left=0.095,
        right=0.985,
        top=0.985,
        bottom=0.085,
        wspace=0.24,
        hspace=0.34,
    )
    datasets = [
        (linear, "Linear PCA", linear_fraction, "PC"),
        (legacy_aligned, "RBF-kPCA, legacy gamma = 0.25", legacy_fraction, "matched kPC"),
        (adaptive_aligned, "RBF-kPCA, median-distance gamma", adaptive_fraction, "matched kPC"),
    ]
    panel_letters = iter("abcdef")
    for row, (values, title, fractions, label_prefix) in enumerate(datasets):
        x1 = f"{label_prefix}1 ({100 * fractions[0]:.1f}%)"
        y1 = f"{label_prefix}2 ({100 * fractions[1]:.1f}%)"
        x2 = f"{label_prefix}2 ({100 * fractions[1]:.1f}%)"
        y2 = f"{label_prefix}3 ({100 * fractions[2]:.1f}%)"
        scatter_panel(axes[row, 0], values, families, 0, 1, f"{title}: axes 1-2", x1, y1)
        scatter_panel(axes[row, 1], values, families, 1, 2, f"{title}: axes 2-3", x2, y2)
        for column in range(2):
            axes[row, column].set_box_aspect(0.72)
            axes[row, column].text(
                -0.15,
                1.08,
                next(panel_letters),
                transform=axes[row, column].transAxes,
                fontweight="bold",
                va="top",
            )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=7,
        frameon=False,
        bbox_to_anchor=(0.5, 0.012),
        columnspacing=1.15,
        handletextpad=0.35,
    )
    save_figure(fig, "kPCA_morphospace_comparison_180mm")

    fig, axes = plt.subplots(1, 2, figsize=(7.0866, 3.15), constrained_layout=True)
    for letter, ax, matrix, title in [
        ("a", axes[0], corr_legacy, "Legacy gamma = 0.25"),
        ("b", axes[1], corr_adaptive, "Median-distance gamma"),
    ]:
        sns.heatmap(
            matrix,
            cmap="vlag",
            center=0,
            vmin=-1,
            vmax=1,
            square=True,
            annot=True,
            fmt=".2f",
            annot_kws={"size": 5},
            cbar_kws={"label": "Pearson r", "shrink": 0.75},
            ax=ax,
        )
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Native kPC")
        ax.set_ylabel("Linear PC")
        ax.set_xticklabels([f"kPC{i}" for i in range(1, matrix.shape[1] + 1)], rotation=45, ha="right")
        ax.set_yticklabels([f"PC{i}" for i in range(1, matrix.shape[0] + 1)], rotation=0)
        ax.text(-0.17, 1.04, letter, transform=ax.transAxes, fontweight="bold", va="top")
    save_figure(fig, "kPCA_axis_correlations_180mm")

    fig, ax = plt.subplots(figsize=(7.0866, 3.25), constrained_layout=True)
    methods = ["Linear PCA", "Legacy kPCA", "Adaptive kPCA"]
    offsets = {method: (index - 1) * 0.22 for index, method in enumerate(methods)}
    x_positions = np.arange(len(FAMILY_ORDER))
    for method in methods:
        block = disparity_long[disparity_long["method"] == method].set_index("Family")
        y = [block.loc[family, "mean_squared_distance"] if family in block.index else np.nan for family in FAMILY_ORDER]
        ax.scatter(x_positions + offsets[method], y, s=30, label=method, zorder=3)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(FAMILY_ORDER, rotation=30, ha="right")
    ax.set_ylabel("Mean squared distance to family centroid")
    ax.set_xlabel("")
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.text(-0.045, 1.03, "a", transform=ax.transAxes, fontweight="bold", va="top")
    save_figure(fig, "kPCA_family_disparity_comparison_180mm")


def main() -> None:
    configure_paths(parse_args())
    ensure_dirs()
    rng = np.random.default_rng(SEED)
    x, score_df, analysis_df, joint_df, pc_columns = load_inputs()

    linear = score_df[pc_columns].to_numpy(dtype=float)
    linear_nonzero = linear[:, np.nanstd(linear, axis=0) > 1e-14]
    if linear_nonzero.shape[1] != x.shape[0] - 1:
        raise ValueError(f"Expected 67 non-zero linear PCs, found {linear_nonzero.shape[1]}")

    # Validate that the flattened momenta reproduce the released ordinary PCA.
    x_centered = x - x.mean(axis=0)
    _, singular_values, right_vectors = np.linalg.svd(x_centered, full_matrices=False)
    recomputed = x_centered @ right_vectors.T
    validation = np.array(
        [[pearsonr(linear_nonzero[:, i], recomputed[:, j]).statistic for j in range(10)] for i in range(10)]
    )
    diagonal_validation = np.abs(np.diag(validation))
    if not np.allclose(diagonal_validation, 1.0, atol=1e-8):
        raise ValueError("The source momenta do not reproduce the released linear PCA")

    squared_distances = squareform(pdist(x, metric="sqeuclidean"))
    off_diagonal = squared_distances[np.triu_indices_from(squared_distances, k=1)]
    median_squared_distance = float(np.median(off_diagonal))
    gamma_legacy = 0.25
    gamma_adaptive = 1.0 / (2.0 * median_squared_distance)

    legacy_scores, legacy_eigenvalues, _ = fit_kpca(x, gamma_legacy)
    adaptive_scores, adaptive_eigenvalues, _ = fit_kpca(x, gamma_adaptive)

    legacy_aligned, legacy_matching = align_axes(linear_nonzero, legacy_scores)
    adaptive_aligned, adaptive_matching = align_axes(linear_nonzero, adaptive_scores)
    legacy_matching.insert(0, "method", "legacy_gamma_0.25")
    adaptive_matching.insert(0, "method", "median_distance_gamma")
    pd.concat([legacy_matching, adaptive_matching], ignore_index=True).to_csv(
        TABLES / "kPCA_axis_matching.csv", index=False
    )

    write_scores(
        score_df["specimen_id"],
        analysis_df["Family"],
        legacy_scores,
        "kPC",
        TABLES / "kPCA_scores_legacy_gamma_0.25.csv",
    )
    write_scores(
        score_df["specimen_id"],
        analysis_df["Family"],
        adaptive_scores,
        "kPC",
        TABLES / "kPCA_scores_median_distance_gamma.csv",
    )
    write_scores(
        score_df["specimen_id"],
        analysis_df["Family"],
        legacy_aligned,
        "matched_kPC",
        TABLES / "kPCA_scores_legacy_first5_matched_to_linear.csv",
    )
    write_scores(
        score_df["specimen_id"],
        analysis_df["Family"],
        adaptive_aligned,
        "matched_kPC",
        TABLES / "kPCA_scores_adaptive_first5_matched_to_linear.csv",
    )

    linear_variance = np.var(linear_nonzero, axis=0, ddof=1)
    linear_fraction = linear_variance / linear_variance.sum()
    legacy_fraction_native = legacy_eigenvalues / legacy_eigenvalues.sum()
    adaptive_fraction_native = adaptive_eigenvalues / adaptive_eigenvalues.sum()

    def matched_fraction(matching: pd.DataFrame, native_fraction: np.ndarray) -> np.ndarray:
        indices = matching["matched_native_kernel_axis"].str.replace("kPC", "", regex=False).astype(int) - 1
        return native_fraction[indices.to_numpy()]

    legacy_fraction = matched_fraction(legacy_matching, legacy_fraction_native)
    adaptive_fraction = matched_fraction(adaptive_matching, adaptive_fraction_native)

    eigen_table = []
    for method, values in [
        ("linear_PCA", linear_variance),
        ("legacy_kPCA", legacy_eigenvalues),
        ("adaptive_kPCA", adaptive_eigenvalues),
    ]:
        fractions = values / values.sum()
        for index, (value, fraction) in enumerate(zip(values, fractions), start=1):
            eigen_table.append(
                {
                    "method": method,
                    "axis": index,
                    "eigenvalue_or_variance": value,
                    "fraction_of_method_total": fraction,
                    "cumulative_fraction": fractions[:index].sum(),
                }
            )
    pd.DataFrame(eigen_table).to_csv(TABLES / "kPCA_eigenspectra.csv", index=False)

    corr_legacy = np.array(
        [[pearsonr(linear_nonzero[:, i], legacy_scores[:, j]).statistic for j in range(10)] for i in range(10)]
    )
    corr_adaptive = np.array(
        [[pearsonr(linear_nonzero[:, i], adaptive_scores[:, j]).statistic for j in range(10)] for i in range(10)]
    )
    pd.DataFrame(corr_legacy, index=[f"PC{i}" for i in range(1, 11)], columns=[f"kPC{i}" for i in range(1, 11)]).to_csv(
        TABLES / "axis_correlations_legacy.csv"
    )
    pd.DataFrame(corr_adaptive, index=[f"PC{i}" for i in range(1, 11)], columns=[f"kPC{i}" for i in range(1, 11)]).to_csv(
        TABLES / "axis_correlations_adaptive.csv"
    )

    comparison_records = []
    for method, values in [
        ("legacy_kPCA", legacy_aligned),
        ("adaptive_kPCA", adaptive_aligned),
    ]:
        _, _, disparity = procrustes(linear_nonzero[:, :5], values)
        distance_rho = spearmanr(pdist(linear_nonzero[:, :5]), pdist(values)).statistic
        comparison_records.append(
            {
                "method": method,
                "procrustes_disparity_first5": disparity,
                "pairwise_distance_spearman_rho_first5": distance_rho,
            }
        )
    pd.DataFrame(comparison_records).to_csv(TABLES / "global_space_similarity.csv", index=False)

    families = analysis_df["Family"].to_numpy()
    family_methods = {
        "Linear PCA": linear_nonzero[:, :5],
        "Legacy kPCA": legacy_aligned,
        "Adaptive kPCA": adaptive_aligned,
    }
    family_test_records = []
    disparity_frames = []
    rarefaction_frames = []
    pairwise_frames = []
    for method, values in family_methods.items():
        test = permutation_group_test(values, families, rng)
        dispersion = dispersion_test(values, families, rng)
        family_test_records.append(
            {
                "method": method,
                "test": "family_location_PERMANOVA",
                **test,
            }
        )
        family_test_records.append(
            {
                "method": method,
                "test": "family_dispersion_PERMDISP",
                "pseudo_F": dispersion["F"],
                "R2": np.nan,
                "df1": dispersion["df1"],
                "df2": dispersion["df2"],
                "permutation_P": dispersion["permutation_P"],
                "permutations": dispersion["permutations"],
            }
        )
        local_disparity = family_disparity(values, families)
        local_disparity.insert(0, "method", method)
        disparity_frames.append(local_disparity)
        local_rarefaction = rarefied_disparity(values, families, rng)
        local_rarefaction.insert(0, "method", method)
        rarefaction_frames.append(local_rarefaction)
        local_pairwise = pairwise_disparity_tests(values, families, rng)
        local_pairwise.insert(0, "method", method)
        pairwise_frames.append(local_pairwise)

    pd.DataFrame(family_test_records).to_csv(TABLES / "family_structure_tests.csv", index=False)
    disparity_long = pd.concat(disparity_frames, ignore_index=True)
    disparity_long.to_csv(TABLES / "family_disparity_comparison.csv", index=False)
    pd.concat(rarefaction_frames, ignore_index=True).to_csv(TABLES / "family_disparity_rarefied_n4.csv", index=False)
    pd.concat(pairwise_frames, ignore_index=True).to_csv(TABLES / "family_disparity_pairwise_tests.csv", index=False)

    joint_merged = analysis_df[["specimen_id"]].merge(
        joint_df[["specimen_id", "joint_type_strict"]], on="specimen_id", how="left", validate="one_to_one"
    )
    joint_mask = joint_merged["joint_type_strict"].notna().to_numpy()
    joint_groups = joint_merged.loc[joint_mask, "joint_type_strict"].to_numpy()
    joint_records = []
    for method, values in family_methods.items():
        result = permutation_group_test(values[joint_mask], joint_groups, rng)
        joint_records.append({"method": method, **result})
    pd.DataFrame(joint_records).to_csv(TABLES / "joint_type_structure_tests.csv", index=False)

    geometry = analysis_df.copy()
    geometry["axial_pitch"] = geometry["axial_span"] / geometry["n_turns_abs"].replace(0, np.nan)
    geometry_eligible = (
        np.isfinite(geometry["abs_winding_angle_deg"].to_numpy(dtype=float))
        & (geometry["abs_winding_angle_deg"].to_numpy(dtype=float) >= GEOMETRY_MIN_ANGLE)
    )
    geometry_records = []
    for trait in ["abs_winding_angle_deg", "axial_span", "axial_pitch"]:
        y = geometry[trait].to_numpy(dtype=float)
        valid = geometry_eligible & np.isfinite(y)
        for method, values in family_methods.items():
            for n_axes in (2, 5):
                result = ols_summary(y[valid], values[valid, :n_axes])
                geometry_records.append(
                    {
                        "response": trait,
                        "method": method,
                        "predictor_axes": n_axes,
                        "minimum_abs_winding_angle_deg": GEOMETRY_MIN_ANGLE,
                        **result,
                    }
                )
    pd.DataFrame(geometry_records).to_csv(TABLES / "shape_geometry_regressions.csv", index=False)

    log_size = np.log(analysis_df["centroid_size"].to_numpy(dtype=float))
    allometry_records = []
    for method, values, space_note in [
        ("Linear PCA", linear_nonzero, "100% of original momenta variance"),
        ("Legacy kPCA", legacy_scores, "100% of positive kernel feature-space inertia"),
        ("Adaptive kPCA", adaptive_scores, "100% of positive kernel feature-space inertia"),
    ]:
        full_result = permutation_regression(values, log_size, rng)
        allometry_records.append(
            {
                "method": method,
                "space": "all_nonzero_axes",
                "axes": values.shape[1],
                "interpretation": space_note,
                **full_result,
            }
        )
    for method, values in family_methods.items():
        result = permutation_regression(values, log_size, rng)
        allometry_records.append(
            {
                "method": method,
                "space": "first5_or_matched_first5",
                "axes": 5,
                "interpretation": "major-shape sensitivity subspace",
                **result,
            }
        )
    pd.DataFrame(allometry_records).to_csv(TABLES / "allometry_sensitivity_tests.csv", index=False)

    kernel_diagnostics = {
        "n_specimens": int(x.shape[0]),
        "n_flattened_momenta_features": int(x.shape[1]),
        "median_squared_pairwise_distance": median_squared_distance,
        "legacy_gamma": gamma_legacy,
        "adaptive_gamma_formula": "1 / (2 * median squared pairwise Euclidean distance)",
        "adaptive_gamma": gamma_adaptive,
        "legacy_off_diagonal_kernel_similarity_quantiles": np.quantile(
            np.exp(-gamma_legacy * off_diagonal), [0, 0.25, 0.5, 0.75, 1]
        ).tolist(),
        "adaptive_off_diagonal_kernel_similarity_quantiles": np.quantile(
            np.exp(-gamma_adaptive * off_diagonal), [0, 0.25, 0.5, 0.75, 1]
        ).tolist(),
        "legacy_positive_kernel_axes": int(legacy_scores.shape[1]),
        "adaptive_positive_kernel_axes": int(adaptive_scores.shape[1]),
    }
    (TABLES / "kernel_parameters.json").write_text(json.dumps(kernel_diagnostics, indent=2), encoding="utf-8")

    manifest = pd.DataFrame(
        [
            {"file": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
            for path in (MOMENTA, LINEAR_SCORES, ANALYSIS_DATA, JOINT_DATA)
        ]
    )
    manifest.to_csv(LOGS / "input_manifest_sha256.csv", index=False)

    make_figures(
        linear_nonzero,
        legacy_aligned,
        adaptive_aligned,
        families,
        linear_fraction,
        legacy_fraction,
        adaptive_fraction,
        corr_legacy,
        corr_adaptive,
        disparity_long,
    )

    readme = f"""# RBF-kPCA sensitivity analysis

This directory is isolated from the released linear-PCA analysis. It does not overwrite manuscript inputs.

## Inputs and validation

- 68 subject-specific Deformetrica momenta vectors, each flattened to 2,700 values.
- The source momenta reproduce the released ordinary PCA exactly for the first ten axes (absolute Pearson r = 1.0 on each matched axis).
- Input file hashes are recorded in `logs/input_manifest_sha256.csv`.

## Kernel choices

1. Legacy workflow: RBF kernel with gamma = {gamma_legacy:.8g}, matching the adapted developer notebook.
2. Data-adaptive sensitivity: RBF kernel with gamma = 1 / (2 * median squared pairwise distance) = {gamma_adaptive:.8g}.

The legacy setting produces a broad, nearly linear kernel. The adaptive setting deliberately tests a stronger nonlinear representation. Kernel eigenvalue fractions describe inertia in the kernel feature space and must not be called percentages of original anatomical shape variance.

## Comparisons

- Morphospace geometry and axis correspondence with the released linear PCA.
- Family location and within-family disparity.
- Joint-type separation.
- Shape associations with winding angle, axial span and endpoint-equivalent axial pitch, consistently restricted to trajectories with absolute winding angle >= {GEOMETRY_MIN_ANGLE:g} degrees.
- Multivariate association with log centroid size.

For direct visual and inferential comparison, the first five kPCA axes were matched to PC1-PC5 by maximum absolute Pearson correlation and sign-aligned. Native kPCA scores are also retained.

Run with:

`python scripts/run_kpca_sensitivity.py --momenta path/to/Atlas_Momentas.txt --linear-scores path/to/PCA_scores_with_specimen_id.csv --analysis-data path/to/PCA_scores_with_specimen_id_with_centroid_size.csv --joint-data path/to/specimen_joint_types.csv --output-dir path/to/kpca_sensitivity_output`
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    print(json.dumps(kernel_diagnostics, indent=2))
    print(f"Outputs written to {OUT}")


if __name__ == "__main__":
    main()
