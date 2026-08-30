#!/usr/bin/env python3
"""Publish the canonical 2026-08-24 Grafen-primary rebuild.

The script copies only machine-readable outputs produced by the completed
Grafen-primary main-dataset workflows. It normalizes released CSVs to
UTF-8/comma/decimal point, labels geometry-dependent rows as main_dataset,
strips workstation paths from tree-robustness outputs, and refreshes the
supplementary-table manifest. It does not run inferential analyses.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path, PureWindowsPath


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


def clean_value(value: str) -> str:
    if value is None:
        return ""
    if DECIMAL_COMMA.fullmatch(value):
        return value.replace(",", ".")
    if ":\\" in value or value.startswith("/"):
        name = PureWindowsPath(value).name
        if name.endswith((".tre", ".tree", ".treefile", ".contree")):
            return f"data/phylogeny/P01_Trees/{name}"
    return value


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean_value(row.get(field, "")) for field in fields})


def copy_csv(source: Path, target: Path) -> None:
    fields, rows = read_csv(source)
    write_csv(target, fields, rows)


def copy_main(main: Path, target: Path) -> None:
    main_fields, main_rows = read_csv(main)
    fields = ["quality_set"] + main_fields
    rows = [dict(quality_set="main_dataset", **row) for row in main_rows]
    write_csv(target, fields, rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def refresh_manifest(root: Path) -> None:
    manifest = root / "_manifest.csv"
    fields, rows = read_csv(manifest)
    kept = []
    for row in rows:
        path = root / Path(row["relative_path"])
        if not path.exists():
            continue
        row["size_bytes"] = str(path.stat().st_size)
        row["sha256"] = sha256(path)
        kept.append(row)

    known = {row["relative_path"] for row in kept}
    for path in root.rglob("*.csv"):
        if path.name == "_manifest.csv":
            continue
        rel = path.relative_to(root).as_posix()
        if rel in known:
            continue
        is_calibration = rel.startswith("S11_Robustness/Calibration_Sensitivity/")
        is_allometry_tip_data = rel in {
            "S03_Allometry/pgls_allometry_all_traits_anova.csv",
            "S03_Allometry/pgls_allometry_tip_data_geometry.csv",
            "S03_Allometry/pgls_allometry_tip_data_shape.csv",
        }
        kept.append({
            "supplementary_table": "22" if is_allometry_tip_data else "",
            "caption": (
                "Table 22: Phylogenetically informed allometry results."
                if is_allometry_tip_data
                else (
                    "Calibration-sensitivity analysis provenance and outputs."
                    if is_calibration
                    else "Auxiliary analysis source data."
                )
            ),
            "relative_path": rel,
            "role": (
                "numbered_table_payload"
                if is_allometry_tip_data
                else (
                    "calibration_sensitivity_source_data"
                    if is_calibration
                    else "auxiliary_source_data"
                )
            ),
            "size_bytes": str(path.stat().st_size),
            "sha256": sha256(path),
        })
    write_csv(manifest, fields, kept)


def mirror_csv_tree(source: Path, target: Path) -> None:
    for path in source.rglob("*.csv"):
        copy_csv(path, target / path.relative_to(source))
    for name in ("README.md",):
        src = source / name
        if src.exists():
            (target / name).write_bytes(src.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo_root", type=Path)
    parser.add_argument("rebuild_root", type=Path)
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()

    repo = args.repo_root.resolve()
    rebuild = args.rebuild_root.resolve()
    project = args.project_root.resolve()
    pcm_main = rebuild / "PCM_Main"
    eco_main = rebuild / "Ecology_Main"
    allom = rebuild / "Allometry_Main"
    phy_allom = rebuild / "Phylogenetic_Multivariate_Allometry"
    unified_allom = rebuild / "Unified_Phylogenetic_Allometry"

    fig = repo / "data/screw_geometry/figure_source_data"
    figure_map = {
        pcm_main / "09_Morphospace_and_tree_plots/figure_phylogenetic_signal_grouped_clean_final_table.csv": fig / "phylogenetic_signal_plot_data.csv",
        pcm_main / "08_ASR/figure_contASR_shape_geometry_standardized_input.csv": fig / "asr_standardized_tip_data.csv",
        pcm_main / "03_Evolutionary_models/figure_univariate_evolutionary_models_dumbbell_table.csv": fig / "evolutionary_models_univariate.csv",
        pcm_main / "03_Evolutionary_models/figure_multivariate_evolutionary_models_dumbbell_table.csv": fig / "evolutionary_models_multivariate.csv",
        pcm_main / "10_Logs/tip_level_dataset_used_for_PCM.csv": fig / "pcm_tip_level_data.csv",
        pcm_main / "04_PGLS/figure_pgls_core_pc1_pc2_vs_axial_span_table.csv": fig / "pgls_core_axial_span.csv",
        pcm_main / "11_Tree_robustness/robustness_pgls_across_trees.csv": fig / "pgls_tree_variant_detail.csv",
        pcm_main / "04_PGLS/robustness_checks/robustness_leave_one_out.csv": fig / "pgls_leave_one_out_detail.csv",
        eco_main / "ecology_analysis_input_merged.csv": fig / "ecology_tip_level_data.csv",
    }
    for name in (
        "allometry_merged_table.csv",
        "allometry_univariate_PC1_to_PC5_results.csv",
        "allometry_continuous_traits_results.csv",
        "allometry_rrpp_multivariate_results.csv",
        "allometry_full_shape_regression_scores.csv",
    ):
        figure_map[allom / name] = fig / name
    for source, target in figure_map.items():
        copy_csv(source, target)

    copy_csv(
        pcm_main / "04_PGLS/pgls_continuous_vs_continuous.csv",
        repo / "data/screw_geometry/sensitivity/pgls_primary_adequate.csv",
    )

    supp = repo / "data/supplementary_source_data"
    s03 = supp / "S03_Allometry"
    for name in (
        "allometry_analysis_scope.csv",
        "allometry_continuous_traits_results.csv",
        "allometry_full_shape_regression_scores.csv",
        "allometry_PC_allometric_contributions.csv",
        "allometry_procD_lm_results.csv",
        "allometry_rrpp_multivariate_results.csv",
        "allometry_rrpp_PC1_to_PC5_sensitivity_results.csv",
        "allometry_univariate_all_PC_results.csv",
        "allometry_univariate_PC1_to_PC5_results.csv",
    ):
        copy_csv(allom / name, s03 / name)
    for name in (
        "pgls_allometry_all_traits.csv",
        "pgls_allometry_all_traits_anova.csv",
        "pgls_allometry_tip_data_geometry.csv",
        "pgls_allometry_tip_data_shape.csv",
        "pgls_geometry_main_traits_main_dataset.csv",
        "pgls_results_main_traits.csv",
    ):
        copy_csv(unified_allom / name, s03 / name)
    for name in (
        "phylogenetic_multivariate_allometry_matching.csv",
        "phylogenetic_multivariate_allometry_projection.csv",
        "phylogenetic_multivariate_allometry_results.csv",
        "phylogenetic_multivariate_allometry_tip_data.csv",
        "phylogenetic_multivariate_allometry_visualization_scores.csv",
    ):
        copy_csv(phy_allom / name, s03 / name)

    combined_outputs = {
        "S08_PCM_Signal_and_Models/phylogenetic_signal_continuous.csv": "02_Phylogenetic_signal/phylogenetic_signal_continuous.csv",
        "S08_PCM_Signal_and_Models/evolutionary_model_fits_univariate.csv": "03_Evolutionary_models/evolutionary_model_fits_univariate.csv",
        "S08_PCM_Signal_and_Models/evolutionary_model_fits_multivariate.csv": "03_Evolutionary_models/evolutionary_model_fits_multivariate.csv",
        "S09_PCM_PGLS/figure_pgls_core_pc1_pc2_vs_axial_span_table.csv": "04_PGLS/figure_pgls_core_pc1_pc2_vs_axial_span_table.csv",
        "S09_PCM_PGLS/pgls_continuous_vs_continuous.csv": "04_PGLS/pgls_continuous_vs_continuous.csv",
        "S09_PCM_PGLS/pgls_continuous_vs_continuous_anova.csv": "04_PGLS/pgls_continuous_vs_continuous_anova.csv",
        "S10_PCM_ASR_and_Disparity/asr_continuous_fastAnc.csv": "08_ASR/asr_continuous_fastAnc.csv",
        "S11_Robustness/robustness_pgls_summary.csv": "11_Tree_robustness/robustness_pgls_summary.csv",
        "S11_Robustness/robustness_pgls_across_trees.csv": "11_Tree_robustness/robustness_pgls_across_trees.csv",
        "S11_Robustness/robustness_phylogenetic_signal_across_trees.csv": "11_Tree_robustness/robustness_phylogenetic_signal_across_trees.csv",
        "S11_Robustness/robustness_phylogenetic_signal_summary.csv": "11_Tree_robustness/robustness_phylogenetic_signal_summary.csv",
        "S11_Robustness/robustness_evolutionary_models_across_trees.csv": "11_Tree_robustness/robustness_evolutionary_models_across_trees.csv",
        "S11_Robustness/robustness_evolutionary_models_best_model_frequency.csv": "11_Tree_robustness/robustness_evolutionary_models_best_model_frequency.csv",
        "S11_Robustness/robustness_asr_root_estimates_across_trees.csv": "11_Tree_robustness/robustness_asr_root_estimates_across_trees.csv",
        "S11_Robustness/robustness_asr_root_estimates_summary.csv": "11_Tree_robustness/robustness_asr_root_estimates_summary.csv",
        "S11_Robustness/robustness_leave_one_out_summary.csv": "04_PGLS/robustness_checks/robustness_leave_one_out_summary.csv",
    }
    for destination, relative in combined_outputs.items():
        copy_main(pcm_main / relative, supp / destination)

    ecology_outputs = {
        "S06_Ecology_Matrix/ecology_analysis_input_merged.csv": "ecology_analysis_input_merged.csv",
        "S07_Ecology_Tests/ecology_phylogenetic_anova_results.csv": "ecology_phylogenetic_anova_results.csv",
        "S07_Ecology_Tests/ecology_pgls_factor_results.csv": "ecology_pgls_factor_results.csv",
        "S07_Ecology_Tests/ecology_trait_group_summary.csv": "ecology_trait_group_summary.csv",
        "S07_Ecology_Tests/ecology_nonphylo_group_tests.csv": "ecology_nonphylo_group_tests.csv",
    }
    for destination, name in ecology_outputs.items():
        copy_main(eco_main / name, supp / destination)

    for path in supp.rglob("*.csv"):
        if path.name != "_manifest.csv":
            copy_csv(path, path)
    refresh_manifest(supp)

    project_tables = project / "05_Supplementary_Tables"
    mirror_csv_tree(supp, project_tables)
    refresh_manifest(project_tables)
    print("Grafen-primary source data and supplementary tables synchronized.")


if __name__ == "__main__":
    main()
