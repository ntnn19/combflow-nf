#!/usr/bin/env python3
"""
map_fasta_to_pdb.py
====================

Maps ColabFold-style FASTA files to the PDB structure file(s) they produced,
by matching the exact multiset of chain sequences (same sequences, same copy
counts, no extra chains in either direction) — not by filename.

ColabFold-style FASTA record:

    >A0_A0_D0_E0_G0
    seq1:seq1:seq2:seq3:seq7

The header's "_"-joined tokens are position-aligned with the ":"-joined
sequences (here: two copies of seq1, one copy each of seq2/seq3/seq7). A PDB
file "matches" a FASTA record if and only if the multiset of its chain
sequences is exactly equal to the multiset of sequences in the record: same
distinct sequences, same copy count for each, nothing extra.

Because ColabFold typically writes several ranked models per job (rank_001..
rank_005), a single FASTA record commonly matches several PDB files — that's
expected, and each match is written as its own row in the output CSV.

Usage
-----
    python map_fasta_to_pdb.py \\
        --fasta-dir fastas/ \\
        --pdb-dir predictions/ \\
        --output-csv mapping.csv

    # also write a report of anything that didn't match on either side
    python map_fasta_to_pdb.py \\
        --fasta-dir fastas/ --pdb-dir predictions/ \\
        --output-csv mapping.csv --unmatched-csv unmatched.csv
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click

try:
    import Bio.PDB
    import Bio.SeqUtils
except ImportError:
    click.echo("This tool requires biopython. Install with: pip install biopython", err=True)
    raise


# --------------------------------------------------------------------------
# FASTA side
# --------------------------------------------------------------------------

@dataclass
class FastaRecord:
    fasta_file: str
    header: str
    sequences: List[str]
    composition: "Counter[str]"  # sequence (uppercased) -> copy count


def parse_fasta_records(path: Path, case_sensitive: bool) -> List[FastaRecord]:
    """Parse a (possibly multi-record) FASTA file into ColabFold-style records."""
    records: List[FastaRecord] = []
    header: Optional[str] = None
    seq_lines: List[str] = []

    def flush():
        if header is None:
            return
        raw_seq = "".join(seq_lines).strip()
        if not raw_seq:
            return
        seqs = [s.strip() for s in raw_seq.split(":") if s.strip()]
        if not case_sensitive:
            seqs = [s.upper() for s in seqs]
        records.append(FastaRecord(
            fasta_file=path.name,
            header=header,
            sequences=seqs,
            composition=Counter(seqs),
        ))

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    flush()

    if not records:
        click.echo(f"Warning: no records parsed from {path}", err=True)
    return records


def discover_files(directory: Path, extensions: List[str]) -> List[Path]:
    exts = {e if e.startswith(".") else f".{e}" for e in extensions}
    files = sorted(p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in exts)
    return files


# --------------------------------------------------------------------------
# PDB / mmCIF side
# --------------------------------------------------------------------------

STANDARD_AA_3TO1_FALLBACK = "X"


def extract_chain_sequences(pdb_path: Path) -> Dict[str, str]:
    """Return {chain_id: one-letter sequence} for a PDB or mmCIF structure."""
    if pdb_path.suffix.lower() == ".cif":
        parser = Bio.PDB.MMCIFParser(QUIET=True)
    else:
        parser = Bio.PDB.PDBParser(QUIET=True)

    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(iter(structure))  # first model only (AFM/ColabFold output has one)

    chain_seqs: Dict[str, str] = {}
    for chain in model:
        residues = [res for res in chain if res.id[0] == " "]  # skip HETATM/water
        if not residues:
            continue
        letters = []
        for res in residues:
            try:
                letters.append(Bio.SeqUtils.seq1(res.get_resname()))
            except Exception:
                letters.append(STANDARD_AA_3TO1_FALLBACK)
        chain_seqs[chain.id] = "".join(letters)
    return chain_seqs


@dataclass
class PdbRecord:
    pdb_file: str
    chain_seqs: Dict[str, str]
    composition: "Counter[str]"


def _parse_one_pdb(path_str: str, case_sensitive: bool) -> Tuple[str, Optional[Dict[str, str]], Optional[str]]:
    """Worker function for the process pool. Returns (filename, chain_seqs_or_None, error_or_None)."""
    path = Path(path_str)
    try:
        chain_seqs = extract_chain_sequences(path)
        if not case_sensitive:
            chain_seqs = {cid: seq.upper() for cid, seq in chain_seqs.items()}
        return path.name, chain_seqs, None
    except Exception as exc:  # noqa: BLE001 - report and continue, don't crash the whole run
        return path.name, None, str(exc)


def parse_pdb_files(paths: List[Path], case_sensitive: bool, workers: int, verbose: bool) -> List[PdbRecord]:
    records: List[PdbRecord] = []
    total = len(paths)
    done = 0

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_parse_one_pdb, str(p), case_sensitive): p for p in paths}
        for future in as_completed(futures):
            name, chain_seqs, error = future.result()
            done += 1
            if verbose and done % max(1, total // 20) == 0:
                click.echo(f"  parsed {done}/{total} PDB files...", err=True)
            if error is not None:
                click.echo(f"Warning: failed to parse {name}: {error}", err=True)
                continue
            if not chain_seqs:
                click.echo(f"Warning: no chains found in {name}", err=True)
                continue
            records.append(PdbRecord(
                pdb_file=name,
                chain_seqs=chain_seqs,
                composition=Counter(chain_seqs.values()),
            ))
    return records


# --------------------------------------------------------------------------
# Matching
# --------------------------------------------------------------------------

def composition_key(composition: "Counter[str]") -> Tuple[Tuple[str, int], ...]:
    """A hashable, order-independent key for an exact-multiset comparison."""
    return tuple(sorted(composition.items()))


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

@click.command()
@click.option("--fasta-dir", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Directory of ColabFold-style FASTA files.")
@click.option("--pdb-dir", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Directory of predicted structure files (.pdb / .cif).")
@click.option("--output-csv", required=True, type=click.Path(path_type=Path),
              help="Path to write the fasta<->pdb mapping CSV.")
@click.option("--fasta-ext", default="fasta,fa", show_default=True,
              help="Comma-separated FASTA file extensions to scan for.")
@click.option("--pdb-ext", default="pdb,cif", show_default=True,
              help="Comma-separated structure file extensions to scan for.")
@click.option("--case-sensitive/--no-case-sensitive", default=False, show_default=True,
              help="Whether sequence comparison is case-sensitive (default: uppercase everything first).")
@click.option("--workers", type=int, default=None,
              help="Parallel worker processes for PDB parsing (default: CPU count).")
@click.option("--unmatched-csv", type=click.Path(path_type=Path), default=None,
              help="Optional path to also write FASTA records and PDB files that matched nothing.")
@click.option("-v", "--verbose", is_flag=True, default=False, help="Print progress to stderr.")
def main(fasta_dir: Path, pdb_dir: Path, output_csv: Path, fasta_ext: str, pdb_ext: str,
         case_sensitive: bool, workers: Optional[int], unmatched_csv: Optional[Path], verbose: bool):
    """Map ColabFold-style FASTA records to matching PDB files by exact chain-sequence composition."""
    import os
    workers = workers or os.cpu_count() or 1

    fasta_ext_list = [e.strip() for e in fasta_ext.split(",") if e.strip()]
    pdb_ext_list = [e.strip() for e in pdb_ext.split(",") if e.strip()]

    fasta_files = discover_files(fasta_dir, fasta_ext_list)
    pdb_files = discover_files(pdb_dir, pdb_ext_list)

    if not fasta_files:
        raise click.ClickException(f"No FASTA files found in {fasta_dir} with extensions {fasta_ext_list}")
    if not pdb_files:
        raise click.ClickException(f"No PDB/CIF files found in {pdb_dir} with extensions {pdb_ext_list}")

    if verbose:
        click.echo(f"Found {len(fasta_files)} FASTA file(s), {len(pdb_files)} structure file(s).", err=True)
        click.echo("Parsing FASTA records...", err=True)

    fasta_records: List[FastaRecord] = []
    for fpath in fasta_files:
        fasta_records.extend(parse_fasta_records(fpath, case_sensitive))

    if verbose:
        click.echo(f"Parsed {len(fasta_records)} FASTA record(s). Parsing PDB structures ({workers} workers)...",
                    err=True)

    pdb_records = parse_pdb_files(pdb_files, case_sensitive, workers, verbose)

    if verbose:
        click.echo(f"Parsed {len(pdb_records)} structure file(s) successfully. Matching...", err=True)

    # Index PDB records by exact composition key for O(1) lookup per fasta record.
    pdb_by_key: Dict[Tuple, List[PdbRecord]] = {}
    for rec in pdb_records:
        pdb_by_key.setdefault(composition_key(rec.composition), []).append(rec)

    rows = []
    matched_pdb_files = set()
    unmatched_fasta: List[FastaRecord] = []

    for frec in fasta_records:
        key = composition_key(frec.composition)
        matches = pdb_by_key.get(key, [])
        if not matches:
            unmatched_fasta.append(frec)
            rows.append({
                "fasta_file": frec.fasta_file,
                "header": frec.header,
                "n_chains": len(frec.sequences),
                "n_unique_sequences": len(frec.composition),
                "pdb_file": "",
                "n_pdb_chains": "",
            })
            continue
        for prec in matches:
            matched_pdb_files.add(prec.pdb_file)
            rows.append({
                "fasta_file": frec.fasta_file,
                "header": frec.header,
                "n_chains": len(frec.sequences),
                "n_unique_sequences": len(frec.composition),
                "pdb_file": prec.pdb_file,
                "n_pdb_chains": len(prec.chain_seqs),
            })

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "fasta_file", "header", "n_chains", "n_unique_sequences", "pdb_file", "n_pdb_chains",
        ])
        writer.writeheader()
        writer.writerows(rows)

    unmatched_pdb = [p for p in pdb_records if p.pdb_file not in matched_pdb_files]

    n_matched_fasta = len({(r["fasta_file"], r["header"]) for r in rows if r["pdb_file"]})
    click.echo(f"\nWrote {len(rows)} mapping row(s) to {output_csv}")
    click.echo(f"  {n_matched_fasta}/{len(fasta_records)} FASTA record(s) matched at least one PDB file")
    click.echo(f"  {len(matched_pdb_files)}/{len(pdb_records)} PDB file(s) matched at least one FASTA record")

    if unmatched_fasta or unmatched_pdb:
        click.echo(f"  {len(unmatched_fasta)} FASTA record(s) had no matching PDB file", err=True)
        click.echo(f"  {len(unmatched_pdb)} PDB file(s) matched no FASTA record", err=True)

    if unmatched_csv:
        unmatched_csv.parent.mkdir(parents=True, exist_ok=True)
        with open(unmatched_csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["kind", "file", "header_or_note"])
            for frec in unmatched_fasta:
                writer.writerow(["fasta_no_match", frec.fasta_file, frec.header])
            for prec in unmatched_pdb:
                writer.writerow(["pdb_no_match", prec.pdb_file, f"{len(prec.chain_seqs)} chains"])
        click.echo(f"Wrote unmatched report to {unmatched_csv}", err=True)


if __name__ == "__main__":
    main()
