#!/usr/bin/env python3

import argparse
import csv
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    # e.g. subunits_000_A0-6_G0-6_manifest.tsv
    # -> subunits_000_A0-6_G0-6/
    name = args.manifest.stem
    if name.endswith("_manifest"):
        name = name[:-9]

    dest = args.output_dir / name
    dest.mkdir(parents=True, exist_ok=True)

    with args.manifest.open(newline="") as f:
        for row in csv.DictReader(f):
            pdb = row["pdb_file"]
            if not pdb:
                continue

            src = args.source_dir / pdb
            if not src.exists():
                print(f"WARNING: not found: {src}")
                continue

            shutil.copy2(src, dest / pdb)
            print(f"{src} -> {dest / pdb}")


if __name__ == "__main__":
    main()
