"""Run a deterministic Deformetrica atlas and normalize principal outputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
CONVERGENCE_RX = re.compile(
    rf"Log-likelihood\s*=\s*({FLOAT}).*?attachment\s*=\s*({FLOAT})"
    rf".*?regularity\s*=\s*({FLOAT})",
    re.IGNORECASE,
)


def require_file(workdir: Path, name: str) -> Path:
    path = workdir / name
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return path


def move_unique(workdir: Path, pattern: str, destination: str) -> Path:
    matches = [path for path in workdir.glob(pattern) if path.name != destination]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one output matching {pattern!r}, found {len(matches)}"
        )
    target = workdir / destination
    if target.exists():
        target.unlink()
    return Path(shutil.move(str(matches[0]), str(target)))


def write_convergence(log_path: Path, output_path: Path) -> int:
    rows = []
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = CONVERGENCE_RX.search(line)
        if match:
            rows.append(" ".join(match.groups()))
    output_path.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--deformetrica", default="deformetrica")
    parser.add_argument("--model", default="model.xml")
    parser.add_argument("--dataset", default="data_set.xml")
    parser.add_argument("--optimization", default="optimization_parameters.xml")
    args = parser.parse_args()

    workdir = args.workdir.expanduser().resolve()
    require_file(workdir, args.model)
    require_file(workdir, args.dataset)
    require_file(workdir, args.optimization)
    require_file(workdir, "initial_template.vtk")

    command = [
        args.deformetrica,
        "estimate",
        args.model,
        args.dataset,
        "-p",
        args.optimization,
        "--output=.",
        "-v",
        "DEBUG",
    ]
    log_path = workdir / "deformetrica.log"
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            log.write(line)
        return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)

    momenta = move_unique(workdir, "*EstimatedParameters__Momenta.txt", "Atlas_Momentas.txt")
    control_points = move_unique(
        workdir, "*EstimatedParameters__ControlPoints.txt", "Atlas_ControlPoints.txt"
    )
    template = move_unique(
        workdir, "*_EstimatedParameters__Template_*.vtk", "Atlas_initial_template.vtk"
    )
    n_rows = write_convergence(log_path, workdir / "convergence.txt")
    print(f"Atlas outputs: {momenta.name}, {control_points.name}, {template.name}")
    print(f"Extracted {n_rows} convergence records")


if __name__ == "__main__":
    main()
