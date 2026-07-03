"""Build a Deformetrica dataset XML file from aligned VTK meshes."""

from __future__ import annotations

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET


def build_dataset(
    mesh_dir: Path,
    output: Path,
    pattern: str,
    object_id: str,
    template_name: str,
) -> list[Path]:
    meshes = sorted(
        path
        for path in mesh_dir.glob(pattern)
        if path.is_file()
        and path.name not in {template_name, "Atlas_initial_template.vtk"}
        and not path.name.startswith("DeterministicAtlas__")
    )
    if not meshes:
        raise RuntimeError(f"No meshes matching {pattern!r} found in {mesh_dir}")

    root = ET.Element("data_set")
    for mesh in meshes:
        subject = ET.SubElement(root, "subject", {"id": mesh.name})
        visit = ET.SubElement(subject, "visit", {"id": "experiment"})
        filename = ET.SubElement(visit, "filename", {"object_id": object_id})
        filename.text = mesh.name

    tree = ET.ElementTree(root)
    ET.indent(tree, space="    ")
    output.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output, encoding="utf-8", xml_declaration=True)
    return meshes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh_dir", type=Path, help="Directory containing aligned VTK meshes")
    parser.add_argument("--output", type=Path, default=Path("data_set.xml"))
    parser.add_argument("--pattern", default="*_aligned.vtk")
    parser.add_argument("--object-id", default="joint")
    parser.add_argument("--template-name", default="initial_template.vtk")
    args = parser.parse_args()

    mesh_dir = args.mesh_dir.expanduser().resolve()
    output = args.output
    if not output.is_absolute():
        output = mesh_dir / output
    meshes = build_dataset(
        mesh_dir,
        output.resolve(),
        args.pattern,
        args.object_id,
        args.template_name,
    )
    print(f"Wrote {output.resolve()} with {len(meshes)} subjects")


if __name__ == "__main__":
    main()
