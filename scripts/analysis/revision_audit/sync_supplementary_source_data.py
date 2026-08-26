#!/usr/bin/env python3
"""Synchronize numbered supplementary-table payloads with robust analyses.

The script keeps the atlas-only 15-tip shape analyses separate from the
specimen-matched robust-geometry analyses (14 main-dataset tips), labels
geometry-dependent rows explicitly as ``main_dataset``, rebuilds the SHA-256
manifest, and removes workstation paths from released tree-robustness tables.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from pathlib import Path, PureWindowsPath


TREE_PATHS = {
    "curc_fig1_grafen_withCaridae_correct.tre": "data/phylogeny/P01_Trees/01_primary_tree_grafen.tre",
    "curc_fig1_withCaridae_calibrated_Grafen.tre": "data/phylogeny/P01_Trees/16_historical_primary_tree_223ma.tre",
    "candidate_historical_223.tre": "data/phylogeny/P01_Trees/16_historical_primary_tree_223ma.tre",
    "candidate_fixed_195.tre": "data/phylogeny/P01_Trees/15_calibration_sensitivity_root195_curculionidae151.tre",
    "candidate_interval_157_3_170.tre": "data/phylogeny/P01_Trees/17_calibration_sensitivity_interval_157_3_170.tre",
    "candidate_interval_157_3_195.tre": "data/phylogeny/P01_Trees/18_calibration_sensitivity_interval_157_3_195.tre",
    "candidate_interval_157_3_223.tre": "data/phylogeny/P01_Trees/19_calibration_sensitivity_interval_157_3_223.tre",
    "curc_fig1_grafen.tre": "data/phylogeny/P01_Trees/07_grafen.tre",
    "curc_fig1_rooted.tre": "data/phylogeny/P01_Trees/06_rooted_ml.tre",
    "curc_fig1_ultrametric_withCaridae_correct.tre": "data/phylogeny/P01_Trees/12_ultrametric_with_caridae_corrected.tre",
    "curc_fig1_ultrametric.tre": "data/phylogeny/P01_Trees/08_ultrametric.tre",
    "curc_fig1_withCaridae_calibrated.tre": "data/phylogeny/P01_Trees/02_primary_tree_calibrated.tre",
    "curc_fig1.contree": "data/phylogeny/P01_Trees/05_consensus_unrooted.contree",
    "curc_fig1.treefile": "data/phylogeny/P01_Trees/04_ml_unrooted.treefile",
}

DECIMAL_COMMA = re.compile(r"^[+-]?(?:\d+,\d*|\d*,\d+)(?:[eE][+-]?\d+)?$")


def delimiter_for(path: Path) -> str:
    sample = path.read_text(encoding="utf-8-sig")[:8192]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;").delimiter
    except csv.Error:
        header = sample.splitlines()[0]
        return ";" if header.count(";") > header.count(",") else ","


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter_for(path))
        return list(reader.fieldnames or []), list(reader)


def sanitize_value(value: str) -> str:
    if not value:
        return value
    if ":\\" in value or value.startswith("/"):
        name = PureWindowsPath(value).name
        return TREE_PATHS.get(name, name)
    if DECIMAL_COMMA.fullmatch(value):
        return value.replace(",", ".")
    return value


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: sanitize_value(row.get(key, "")) for key in fieldnames})


def copy_main(primary: Path, output: Path) -> None:
    p_fields, p_rows = read_csv(primary)
    fields = ["quality_set"] + p_fields
    rows = [dict(quality_set="main_dataset", **row) for row in p_rows]
    write_csv(output, fields, rows)


def copy(source: Path, target: Path) -> None:
    if source.suffix.lower() == ".csv":
        fields, rows = read_csv(source)
        write_csv(target, fields, rows)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def filter_shape_only(path: Path) -> None:
    fields, rows = read_csv(path)
    rows = [row for row in rows if row.get("trait", "").startswith("PC")]
    write_csv(path, fields, rows)


def sanitize_file(path: Path) -> None:
    fields, rows = read_csv(path)
    write_csv(path, fields, rows)


def rebuild_main_results(
    regression: Path, output: Path, expected_full_n: str, expected_main_n: str
) -> None:
    _, rows = read_csv(regression)
    wanted = {
        ("angle_abs ~ PC1 + PC2", "full_dataset"),
        ("axial_pitch ~ PC1 + PC2", "full_dataset"),
        ("angle_abs ~ PC1 + PC2", "main_region_PC1_lt_0.1"),
    }
    selected = [
        row for row in rows
        if row.get("table_block") == "regression_models"
        and (row.get("model"), row.get("subset")) in wanted
    ]
    expected = {
        ("angle_abs ~ PC1 + PC2", "full_dataset"): expected_full_n,
        ("axial_pitch ~ PC1 + PC2", "full_dataset"): expected_full_n,
        ("angle_abs ~ PC1 + PC2", "main_region_PC1_lt_0.1"): expected_main_n,
    }
    if any(row.get("n") != expected[(row["model"], row["subset"])] for row in selected) or len(selected) != 3:
        raise ValueError(f"Unexpected rows or sample size in {regression}")
    model_names = {
        "angle_abs ~ PC1 + PC2": "Winding angle ~ shape",
        "axial_pitch ~ PC1 + PC2": "Axial pitch ~ shape",
    }
    subset_names = {
        "full_dataset": "Full dataset",
        "main_region_PC1_lt_0.1": "Main region (PC1 < 0.1)",
    }
    output_rows = [{
        "model": model_names[row["model"]],
        "subset": subset_names[row["subset"]],
        "n": row["n"],
        "r_squared": row["r_squared"],
        "p_value": row["p_model"],
        "PC1_effect": row["PC1_estimate"],
        "PC1_p": row["PC1_p"],
    } for row in selected]
    write_csv(
        output,
        ["model", "subset", "n", "r_squared", "p_value", "PC1_effect", "PC1_p"],
        output_rows,
    )


def rebuild_manifest(source_root: Path, additions: dict[str, tuple[str, str]]) -> None:
    manifest_path = source_root / "_manifest.csv"
    fields, old_rows = read_csv(manifest_path)
    old_by_path = {row["relative_path"]: row for row in old_rows}
    retained: list[dict[str, str]] = []
    for rel, row in old_by_path.items():
        if (source_root / Path(rel)).exists():
            retained.append(row)
    retained_by_path = {row["relative_path"]: row for row in retained}
    for rel, (number, title) in additions.items():
        if rel not in retained_by_path:
            row = {
                "supplementary_table": number,
                "caption": title,
                "relative_path": rel,
                "role": "numbered_table_payload",
                "size_bytes": "",
                "sha256": "",
            }
            retained.append(row)
            retained_by_path[rel] = row
        else:
            row = retained_by_path[rel]
            row["supplementary_table"] = number
            row["caption"] = title
            row["role"] = "numbered_table_payload"

    # Keep the release manifest exhaustive. Calibration-sensitivity outputs are
    # auxiliary source data rather than numbered table payloads, but they still
    # need checksums so a checkout can be verified without consulting a second
    # inventory.
    retained_by_path = {row["relative_path"]: row for row in retained}
    for path in source_root.rglob("*.csv"):
        if path.name == "_manifest.csv":
            continue
        rel = path.relative_to(source_root).as_posix()
        if rel in retained_by_path:
            continue
        is_calibration = rel.startswith("S11_Robustness/Calibration_Sensitivity/")
        row = {
            "supplementary_table": "",
            "caption": (
                "Calibration-sensitivity analysis provenance and outputs."
                if is_calibration
                else "Auxiliary analysis source data."
            ),
            "relative_path": rel,
            "role": (
                "calibration_sensitivity_source_data"
                if is_calibration
                else "auxiliary_source_data"
            ),
            "size_bytes": "",
            "sha256": "",
        }
        retained.append(row)
        retained_by_path[rel] = row
    for row in retained:
        path = source_root / Path(row["relative_path"])
        row["size_bytes"] = str(path.stat().st_size)
        row["sha256"] = sha256(path)

    def sort_key(row: dict[str, str]) -> tuple[int, str]:
        number = row.get("supplementary_table", "").split(";", 1)[0]
        return (int(number) if number.isdigit() else 999, row["relative_path"])

    write_csv(manifest_path, fields, sorted(retained, key=sort_key))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo_root", type=Path)
    parser.add_argument("robust_point_estimates", type=Path)
    args = parser.parse_args()

    repo = args.repo_root.resolve()
    robust = args.robust_point_estimates.resolve()
    source = repo / "data" / "supplementary_source_data"
    screw = repo / "data" / "screw_geometry"
    primary = robust / "primary_adequate"

    # Rebuild the concise model table before copying it into the table payload.
    rebuild_main_results(
        screw / "regression_summary.csv", screw / "main_results.csv", "63", "55"
    )
    # Table 13: all successful fits plus the 63-fit main analysis dataset.
    t13 = source / "S02_Shape_Geometry"
    copy(screw / "robust_geometry_primary_adequate.csv", t13 / "specimen_level_screw_joint_geometry.csv")
    copy(screw / "robust_geometry_all.csv", t13 / "robust_geometry_all.csv")
    copy(screw / "robust_helix_metrics.csv", t13 / "robust_helix_metrics.csv")

    # Tables 14--17: canonical robust summaries and uncertainty propagation.
    for stale in (t13 / "Table_main_results.csv", t13 / "Table_regression_summary.csv"):
        if stale.exists():
            stale.unlink()
    copy(screw / "main_results.csv", t13 / "main_results.csv")
    copy(screw / "shape_model_uncertainty_summary.csv", t13 / "shape_model_uncertainty_summary.csv")
    copy(screw / "regression_summary.csv", t13 / "regression_summary.csv")

    t16 = source / "S09_PCM_PGLS"
    copy_main(
        primary / "pcm_matched/04_PGLS/figure_pgls_core_pc1_pc2_vs_axial_span_table.csv",
        t16 / "figure_pgls_core_pc1_pc2_vs_axial_span_table.csv",
    )
    copy(screw / "pgls_shape_geometry_rubin_summary.csv", t16 / "pgls_shape_geometry_rubin_summary.csv")

    t17 = source / "S05_Joint_Typology"
    copy(screw / "joint_type_geometry_stats.csv", t17 / "joint_type_screw_geometry_stats.csv")
    copy(screw / "joint_type_uncertainty_summary.csv", t17 / "joint_type_uncertainty_summary.csv")

    # Tables 19 and 22: specimen-level and phylogenetic allometry.
    t19 = source / "S03_Allometry"
    copy_main(
        primary / "allometry/allometry_continuous_traits_results.csv",
        t19 / "allometry_continuous_traits_results.csv",
    )
    filter_shape_only(t19 / "pgls_results_main_traits.csv")
    copy(
        primary / "allometry/Allometry_Phylogenetic/pgls_results_main_traits.csv",
        t19 / "pgls_geometry_main_traits_main_dataset.csv",
    )

    # Tables 23--34: matched-tip PCM outputs for the main dataset.
    matched_outputs = {
        "S08_PCM_Signal_and_Models/phylogenetic_signal_continuous.csv": "pcm_matched/02_Phylogenetic_signal/phylogenetic_signal_continuous.csv",
        "S08_PCM_Signal_and_Models/evolutionary_model_fits_univariate.csv": "pcm_matched/03_Evolutionary_models/evolutionary_model_fits_univariate.csv",
        "S08_PCM_Signal_and_Models/evolutionary_model_fits_multivariate.csv": "pcm_matched/03_Evolutionary_models/evolutionary_model_fits_multivariate.csv",
        "S09_PCM_PGLS/pgls_continuous_vs_continuous.csv": "pcm_matched/04_PGLS/pgls_continuous_vs_continuous.csv",
        "S09_PCM_PGLS/pgls_continuous_vs_continuous_anova.csv": "pcm_matched/04_PGLS/pgls_continuous_vs_continuous_anova.csv",
        "S10_PCM_ASR_and_Disparity/asr_continuous_fastAnc.csv": "pcm_matched/08_ASR/asr_continuous_fastAnc.csv",
        "S11_Robustness/robustness_pgls_summary.csv": "pcm_matched/11_Tree_robustness/robustness_pgls_summary.csv",
        "S11_Robustness/robustness_pgls_across_trees.csv": "pcm_matched/11_Tree_robustness/robustness_pgls_across_trees.csv",
        "S11_Robustness/robustness_phylogenetic_signal_across_trees.csv": "pcm_matched/11_Tree_robustness/robustness_phylogenetic_signal_across_trees.csv",
        "S11_Robustness/robustness_phylogenetic_signal_summary.csv": "pcm_matched/11_Tree_robustness/robustness_phylogenetic_signal_summary.csv",
        "S11_Robustness/robustness_evolutionary_models_across_trees.csv": "pcm_matched/11_Tree_robustness/robustness_evolutionary_models_across_trees.csv",
        "S11_Robustness/robustness_evolutionary_models_best_model_frequency.csv": "pcm_matched/11_Tree_robustness/robustness_evolutionary_models_best_model_frequency.csv",
        "S11_Robustness/robustness_asr_root_estimates_across_trees.csv": "pcm_matched/11_Tree_robustness/robustness_asr_root_estimates_across_trees.csv",
        "S11_Robustness/robustness_asr_root_estimates_summary.csv": "pcm_matched/11_Tree_robustness/robustness_asr_root_estimates_summary.csv",
        "S11_Robustness/robustness_leave_one_out_summary.csv": "pcm_matched/04_PGLS/robustness_checks/robustness_leave_one_out_summary.csv",
    }
    for destination, relative in matched_outputs.items():
        copy_main(primary / relative, source / destination)

    # Tables 35--39: ecology matrices and tests for the main dataset.
    ecology_outputs = {
        "S07_Ecology_Tests/ecology_phylogenetic_anova_results.csv": "ecology/ecology_phylogenetic_anova_results.csv",
        "S07_Ecology_Tests/ecology_pgls_factor_results.csv": "ecology/ecology_pgls_factor_results.csv",
        "S07_Ecology_Tests/ecology_trait_group_summary.csv": "ecology/ecology_trait_group_summary.csv",
        "S07_Ecology_Tests/ecology_nonphylo_group_tests.csv": "ecology/ecology_nonphylo_group_tests.csv",
        "S06_Ecology_Matrix/ecology_analysis_input_merged.csv": "ecology/ecology_analysis_input_merged.csv",
    }
    for destination, relative in ecology_outputs.items():
        copy_main(primary / relative, source / destination)

    # Remove the one remaining private workstation path from a released table.
    sanitize_file(screw / "figure_source_data/pgls_tree_variant_detail.csv")

    # Publish one machine-readable CSV convention throughout the package,
    # independent of the locale used by the originating R session.
    for path in source.rglob("*.csv"):
        if path.name != "_manifest.csv":
            sanitize_file(path)

    titles = {
        "13": "Table 13: Specimen-level screw joint geometry measurements.",
        "14": "Table 14: Main shape-geometry analysis results.",
        "15": "Table 15: Regression summary for shape-geometry relationships.",
        "16": "Table 16: Core PGLS relationships between shape and axial span.",
        "17": "Table 17: Screw joint geometry statistics by joint type.",
        "22": "Table 22: Phylogenetically informed allometry results.",
    }
    additions: dict[str, tuple[str, str]] = {}
    for name in ("robust_geometry_all.csv", "robust_helix_metrics.csv"):
        additions[f"S02_Shape_Geometry/{name}"] = ("13", titles["13"])
    for name in ("main_results.csv", "shape_model_uncertainty_summary.csv"):
        additions[f"S02_Shape_Geometry/{name}"] = ("14", titles["14"])
    for name in ("regression_summary.csv",):
        additions[f"S02_Shape_Geometry/{name}"] = ("15", titles["15"])
    additions["S09_PCM_PGLS/pgls_shape_geometry_rubin_summary.csv"] = ("16", titles["16"])
    additions["S05_Joint_Typology/joint_type_uncertainty_summary.csv"] = ("17", titles["17"])
    additions["S03_Allometry/pgls_geometry_main_traits_main_dataset.csv"] = ("22", titles["22"])
    rebuild_manifest(source, additions)

    print("Supplementary source data synchronized and manifest rebuilt.")


if __name__ == "__main__":
    main()
