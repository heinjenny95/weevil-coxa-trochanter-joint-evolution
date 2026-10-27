"""Repackage released observed XYZ points as ordered helix-input CSV files.

These are reconstructed inputs, not the original Cinema 4D exports. Original
normal vectors are not reconstructed or needed by the robust helix fitter.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('point_residuals_csv', type=Path)
    parser.add_argument('metrics_csv', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    groups = defaultdict(list)
    with args.point_residuals_csv.open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            groups[row['specimen_id']].append(row)
    with args.metrics_csv.open(encoding='utf-8-sig', newline='') as handle:
        metadata = {row['specimen_id']: row for row in csv.DictReader(handle)}
    if not groups or not set(groups).issubset(metadata):
        raise ValueError('Every trace needs a corresponding specimen metrics record.')
    args.output_directory.mkdir(parents=True, exist_ok=True)
    if any(args.output_directory.iterdir()):
        raise ValueError('Choose an empty output directory; existing inputs are preserved.')
    for specimen, rows in sorted(groups.items()):
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', specimen):
            raise ValueError(f'Unsafe specimen filename: {specimen!r}')
        rows.sort(key=lambda row: int(row['point_index']))
        if len(set(row['point_index'] for row in rows)) != len(rows):
            raise ValueError(f'Duplicate point indices for {specimen}')
        if int(metadata[specimen]['n_points']) != len(rows):
            raise ValueError(f'Point count differs from the released metrics for {specimen}')
        units = metadata[specimen].get('units') or 'unknown'
        with (args.output_directory / f'{specimen}.csv').open('w', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, lineterminator='\n')
            writer.writerow(['Landmarks', len(rows), 'reconstructed-from-released-XYZ', f'({units})'])
            writer.writerow(['#', 'Name', 'Type', 'X', 'Y', 'Z'])
            for row in rows:
                writer.writerow([row['point_index'], row['point_name'], 'normal', row['x'], row['y'], row['z']])
    print(f'Reconstructed {len(groups)} ordered traces from released observed coordinates.')


if __name__ == '__main__':
    main()
