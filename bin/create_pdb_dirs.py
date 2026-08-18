#!/usr/bin/env python3
"""Stage the PDB files of one combination into a flat directory.

Usage: create_pdb_dirs.py <manifest.tsv> <source_dir> <output_dir>

The manifest (produced by map_fasta_to_pdb) must have a "pdb_file" column
with paths relative to source_dir. Files are hardlinked (fallback: copied)
directly into output_dir - run_on_pdbs.py scans the PDB dir non-recursively,
so the .pdb files must sit flat in <comb>_pdbs/ with no nested subdir.
Exits 1 if any manifest entry is missing from source_dir.
"""

import argparse
import csv
import os
import shutil
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    dest = args.output_dir
    dest.mkdir(parents=True, exist_ok=True)

    n_linked, missing = 0, 0
    with args.manifest.open(newline="") as f:
        for row in csv.DictReader(f):
            pdb = row["pdb_file"]
            if not pdb:
                continue

            src = args.source_dir / pdb
            if not src.exists():
                print(f"ERROR: not found: {src}", file=sys.stderr)
                missing += 1
                continue

            dst = dest / pdb
            if dst.exists():
                continue
            try:
                os.link(src, dst)  # hardlink: no data duplication
            except OSError:
                shutil.copy2(src, dst)  # fallback across filesystems
            n_linked += 1
            print(f"{src} -> {dst}")

    if missing:
        print(f"ERROR: {missing} manifest PDB(s) missing from {args.source_dir}",
              file=sys.stderr)
        sys.exit(1)
    print(f"Done: {n_linked} PDB(s) staged in {dest}")


if __name__ == "__main__":
    main()
