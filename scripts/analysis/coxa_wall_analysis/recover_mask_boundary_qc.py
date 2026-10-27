"""Recover the extraction script's missing boundary flag from original masks.

This reads the six faces of each unmodified original mask. It does not repeat
local-thickness extraction or estimate a boundary flag from a cropped mask.
"""
import argparse
import csv
import hashlib
from pathlib import Path
import numpy as np
import tifffile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('metrics_csv', type=Path)
    parser.add_argument('mask_directory', type=Path)
    parser.add_argument('output_csv', type=Path)
    args = parser.parse_args()
    with args.metrics_csv.open(encoding='utf-8-sig', newline='') as handle:
        metrics = list(csv.DictReader(handle))
    rows = []
    for row in metrics:
        path = args.mask_directory / Path(row['tif_file']).name
        raw = tifffile.imread(path)
        if raw.ndim != 3:
            raise ValueError(f'Expected 3D original mask: {path.name}')
        shape = 'x'.join(map(str, raw.shape))
        if shape != row['original_stack_shape_zyx']:
            raise ValueError(f'Original-stack shape differs: {path.name}')
        boundary = any(np.any(face > 0) for face in (
            raw[0], raw[-1], raw[:, 0], raw[:, -1], raw[:, :, 0], raw[:, :, -1]))
        rows.append({
            'specimen': row['specimen'],
            'touches_original_stack_boundary': bool(boundary),
            'mask_filename': path.name,
            'mask_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'original_stack_shape_zyx': shape,
        })
        del raw
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'Verified {len(rows)} original masks; {sum(r["touches_original_stack_boundary"] for r in rows)} touch a stack boundary.')


if __name__ == '__main__':
    main()
