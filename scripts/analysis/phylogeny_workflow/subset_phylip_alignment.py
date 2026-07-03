#!/usr/bin/env python3

"""Subset a sequential PHYLIP alignment by taxon name.

The script reads a one-line-per-taxon sequential PHYLIP file, keeps only taxa
listed in a plain-text file, and writes a new PHYLIP alignment with an updated
taxon count. Taxon names must match the first token of each PHYLIP record.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_keep_taxa(path: Path) -> set[str]:
    return {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def read_phylip(path: Path) -> tuple[int, list[tuple[str, str]]]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip().split()
        if len(header) < 2:
            raise ValueError("Input PHYLIP header must contain ntax and nchar.")

        try:
            nchar = int(header[1])
        except ValueError as exc:
            raise ValueError("Could not parse nchar from PHYLIP header.") from exc

        records: list[tuple[str, str]] = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            taxon = fields[0]
            sequence = "".join(fields[1:])
            records.append((taxon, sequence))

    return nchar, records


def write_phylip(path: Path, nchar: int, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(records)} {nchar}\n")
        for taxon, sequence in records:
            handle.write(f"{taxon} {sequence}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_phylip", type=Path)
    parser.add_argument("keep_taxa", type=Path)
    parser.add_argument("output_phylip", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    keep = read_keep_taxa(args.keep_taxa)
    nchar, records = read_phylip(args.input_phylip)
    filtered = [record for record in records if record[0] in keep]

    if not filtered:
        raise SystemExit(
            "ERROR: no taxa matched. Check that keep_taxa names match PHYLIP labels."
        )

    missing = sorted(keep - {taxon for taxon, _ in records})
    if missing:
        print(
            "WARNING: requested taxa not found in alignment: "
            + ", ".join(missing[:20])
            + (" ..." if len(missing) > 20 else "")
        )

    write_phylip(args.output_phylip, nchar, filtered)
    print(f"Wrote {len(filtered)} taxa to {args.output_phylip}")


if __name__ == "__main__":
    main()
