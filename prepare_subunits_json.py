#!/usr/bin/env python3
"""
combfold_stoich_screen — Automate CombFold Stage 1 for stoichiometric screening.

Generates CombFold-compatible subunits.json files from a FASTA of unique
sequences and a per-chain copy-number range specification.  One JSON is
produced per stoichiometry combination (Cartesian product of ranges).

CombFold Stage 1 reference:
  https://github.com/dina-lab3D/CombFold  (README, "Stage 1 - Defining subunits")

Each subunits.json entry has four fields:
  name         — unique subunit identifier (must match the dict key)
  sequence     — amino acid sequence
  chain_names  — list of chain labels representing stoichiometry
  start_res    — 1-based index of the first residue on the chain
"""

from __future__ import annotations

import itertools
import json
import math
import os
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import click

# ──────────────────────────────────────────────────────────────────────────────
# Chain-label helpers (Excel-column style: A, B, …, Z, AA, AB, …)
# ──────────────────────────────────────────────────────────────────────────────

_LABEL_RE = re.compile(r"^[A-Z]+$")


def label_to_index(label: str) -> int:
    """Convert an Excel-column-style label to a 0-based index.

    A → 0, B → 1, …, Z → 25, AA → 26, AB → 27, …
    """
    result = 0
    for char in label:
        result = result * 26 + (ord(char) - ord("A") + 1)
    return result - 1


def index_to_label(index: int) -> str:
    """Convert a 0-based index to an Excel-column-style label.

    0 → A, 1 → B, …, 25 → Z, 26 → AA, 27 → AB, …
    """
    if index < 0:
        raise ValueError(f"Chain index must be non-negative, got {index}")
    parts: List[str] = []
    n = index
    while n >= 0:
        parts.append(chr(ord("A") + (n % 26)))
        n = n // 26 - 1
    return "".join(reversed(parts))


def generate_chain_labels(
    start_label: str, count: int, offset: int = 0
) -> List[str]:
    """Generate *count* sequential chain labels starting from *start_label* + *offset*."""
    if count < 1:
        return []
    base = label_to_index(start_label)
    return [index_to_label(base + offset + i) for i in range(count)]


# ──────────────────────────────────────────────────────────────────────────────
# FASTA parsing
# ──────────────────────────────────────────────────────────────────────────────


def parse_fasta(path: str) -> Dict[str, str]:
    """Parse a FASTA file into an ordered dict ``{header: sequence}``.

    Headers are the text after ``>`` (stripped).  Sequences are uppercased
    and whitespace-stripped.  Duplicate headers and empty sequences raise
    ``ValueError``.
    """
    sequences: Dict[str, str] = {}
    current_header: Optional[str] = None
    current_seq: List[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    _flush(sequences, current_header, current_seq)
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        _flush(sequences, current_header, current_seq)

    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {path}")
    return sequences


def _flush(sequences: Dict[str, str], header: str, parts: List[str]) -> None:
    seq = "".join(parts).upper()
    if not seq:
        raise ValueError(f"Empty sequence for header '{header}'")
    if header in sequences:
        raise ValueError(f"Duplicate header '{header}' in FASTA file")
    sequences[header] = seq


# ──────────────────────────────────────────────────────────────────────────────
# Stoichiometry parsing & expansion
# ──────────────────────────────────────────────────────────────────────────────


def parse_stoichiometry(
    spec: str, fasta_headers: Sequence[str]
) -> Dict[str, Tuple[int, int]]:
    """Parse a stoichiometry range string into ``{header: (min, max)}``.

    Format: ``"alpha:1-6,beta:1-4,gamma:3"`` (single value = min == max).
    Each header must exist in *fasta_headers*.  Returns only the headers
    that were explicitly listed.
    """
    ranges: Dict[str, Tuple[int, int]] = {}
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(
                f"Invalid stoichiometry entry '{part}' — "
                f"expected 'header:min-max' or 'header:value'"
            )
        header, range_str = part.split(":", 1)
        header = header.strip()
        range_str = range_str.strip()

        if header not in fasta_headers:
            raise ValueError(
                f"Stoichiometry header '{header}' not found in FASTA. "
                f"Available headers: {', '.join(fasta_headers)}"
            )
        if header in ranges:
            raise ValueError(
                f"Duplicate stoichiometry entry for header '{header}'"
            )

        if "-" in range_str:
            min_str, max_str = range_str.split("-", 1)
            try:
                lo = int(min_str.strip())
                hi = int(max_str.strip())
            except ValueError:
                raise ValueError(
                    f"Invalid range '{range_str}' for '{header}' — "
                    f"expected 'min-max' with integers"
                )
        else:
            try:
                lo = hi = int(range_str)
            except ValueError:
                raise ValueError(
                    f"Invalid value '{range_str}' for '{header}' — "
                    f"expected an integer"
                )

        if lo < 1:
            raise ValueError(
                f"Minimum copies for '{header}' must be >= 1, got {lo}"
            )
        if lo > hi:
            raise ValueError(
                f"Minimum copies for '{header}' ({lo}) exceeds maximum ({hi})"
            )
        ranges[header] = (lo, hi)
    return ranges


def expand_stoichiometry(
    ranges: Dict[str, Tuple[int, int]],
    fasta_headers: Sequence[str],
    max_combinations: int,
) -> List[List[Tuple[str, int]]]:
    """Expand stoichiometry ranges into all combinations (Cartesian product).

    Headers in *fasta_headers* but not in *ranges* default to 1 copy.
    Returns a list of combinations, each an ordered list of
    ``(header, copy_count)`` pairs in FASTA order.
    """
    # Build (header, range) pairs in FASTA order; default missing to (1, 1).
    ordered: List[Tuple[str, Tuple[int, int]]] = []
    missing: List[str] = []
    for header in fasta_headers:
        if header in ranges:
            ordered.append((header, ranges[header]))
        else:
            ordered.append((header, (1, 1)))
            missing.append(header)

    if missing:
        click.echo(
            f"Warning: headers not in --stoichiometry default to 1 copy: "
            f"{', '.join(missing)}",
            err=True,
        )

    count_lists = [list(range(lo, hi + 1)) for _, (lo, hi) in ordered]
    total = 1
    for cl in count_lists:
        total *= len(cl)

    if total > max_combinations:
        raise ValueError(
            f"Stoichiometry expansion produces {total} combinations, "
            f"exceeding --max-combinations ({max_combinations}). "
            f"Narrow your ranges or raise --max-combinations."
        )

    headers = [h for h, _ in ordered]
    return [
        list(zip(headers, combo)) for combo in itertools.product(*count_lists)
    ]


# ──────────────────────────────────────────────────────────────────────────────
# Chain splitting (naive equal-length)
# ──────────────────────────────────────────────────────────────────────────────


def split_sequence(
    sequence: str, max_af_size: int, start_res: int = 1
) -> List[Tuple[str, int]]:
    """Split *sequence* into equal-length subunits.

    Returns a list of ``(subsequence, start_res)`` tuples.
    If *max_af_size* <= 0 or the sequence fits, returns a single subunit.
    """
    length = len(sequence)
    if max_af_size <= 0 or length <= max_af_size:
        return [(sequence, start_res)]

    k = math.ceil(length / max_af_size)
    part_size = math.ceil(length / k)  # guarantees each part <= max_af_size

    parts: List[Tuple[str, int]] = []
    for i in range(k):
        start = i * part_size
        end = min((i + 1) * part_size, length)
        parts.append((sequence[start:end], start_res + start))
    return parts


# ──────────────────────────────────────────────────────────────────────────────
# Subunit building
# ──────────────────────────────────────────────────────────────────────────────


def build_subunits(
    combination: List[Tuple[str, int]],
    sequences: Dict[str, str],
    max_af_size: int,
    start_res: int,
    chain_start: str,
) -> Dict[str, dict]:
    """Build the subunits dict for one stoichiometry combination.

    *combination* is an ordered list of ``(header, copy_count)`` pairs.
    Returns ``{subunit_name: {name, chain_names, start_res, sequence}}``.
    """
    subunits: Dict[str, dict] = {}
    chain_offset = 0

    for header, copy_count in combination:
        if copy_count < 1:
            raise ValueError(
                f"Copy count for '{header}' must be >= 1, got {copy_count}"
            )

        chain_names = generate_chain_labels(chain_start, copy_count, chain_offset)
        chain_offset += copy_count

        parts = split_sequence(sequences[header], max_af_size, start_res)
        first_chain = chain_names[0]

        for i, (subseq, sub_start_res) in enumerate(parts):
            name = f"{first_chain}{i}"
            if name in subunits:
                raise ValueError(
                    f"Duplicate subunit name '{name}' — this indicates a bug "
                    f"in chain-label assignment"
                )
            subunits[name] = {
                "name": name,
                "chain_names": list(chain_names),
                "start_res": sub_start_res,
                "sequence": subseq,
            }

    return subunits


# ──────────────────────────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────────────────────────


def validate_subunits(subunits: Dict[str, dict]) -> List[str]:
    """Validate a subunits dict against CombFold constraints.

    Returns a list of error strings (empty if valid).  Checks:
    - dict key == ``name`` field
    - all subunit names unique
    - all ``start_res`` >= 1
    - no two subunits claim the same residue on the same chain
    """
    errors: List[str] = []

    # key == name
    for key, info in subunits.items():
        if key != info["name"]:
            errors.append(f"Key '{key}' does not match name '{info['name']}'")

    # unique names
    names = [info["name"] for info in subunits.values()]
    if len(names) != len(set(names)):
        duplicates = {n for n in names if names.count(n) > 1}
        errors.append(f"Duplicate subunit names: {', '.join(sorted(duplicates))}")

    # start_res >= 1
    for info in subunits.values():
        if info["start_res"] < 1:
            errors.append(
                f"Subunit '{info['name']}' has start_res < 1: {info['start_res']}"
            )

    # no overlapping residues on the same chain
    claimed: Dict[Tuple[str, int], str] = {}
    for info in subunits.values():
        for chain_name in info["chain_names"]:
            for i in range(len(info["sequence"])):
                res_id = info["start_res"] + i
                ck = (chain_name, res_id)
                if ck in claimed:
                    errors.append(
                        f"Residue {res_id} on chain '{chain_name}' claimed by "
                        f"both '{claimed[ck]}' and '{info['name']}'"
                    )
                else:
                    claimed[ck] = info["name"]

    return errors


# ──────────────────────────────────────────────────────────────────────────────
# Output helpers
# ──────────────────────────────────────────────────────────────────────────────


def stoich_string(combination: List[Tuple[str, int]]) -> str:
    """Compact stoichiometry string: ``alpha-2_beta-1``."""
    return "_".join(f"{h}-{c}" for h, c in combination)


def sanitize_filename(name: str) -> str:
    """Replace characters that are unsafe in filenames."""
    return re.sub(r"[^a-zA-Z0-9._-]", "_", name)


def format_filename(pattern: str, index: int, stoich: str) -> str:
    """Format an output filename from the pattern template."""
    return pattern.format(index=index, stoich=sanitize_filename(stoich))


def write_subunits_json(subunits: Dict[str, dict], path: str) -> None:
    """Write subunits dict to a JSON file with indent=2."""
    with open(path, "w") as fh:
        json.dump(subunits, fh, indent=2)
        fh.write("\n")


def write_manifest(
    rows: List[dict],
    output_path: str,
    headers: Sequence[str],
) -> None:
    """Write a TSV manifest mapping output files to stoichiometries.

    Columns: index, filename, total_chains, total_subunits, then one per header.
    """
    columns = ["index", "filename", "total_chains", "total_subunits"] + list(headers)
    with open(output_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(col, "")) for col in columns) + "\n")


# ──────────────────────────────────────────────────────────────────────────────
# Click CLI
# ──────────────────────────────────────────────────────────────────────────────

EPILOG = (
    "Example:\n"
    "  combfold-stoich-screen -f complex.fasta -s 'alpha:1-6,beta:1-4' "
    "-o ./screen_out/\n"
    "\n"
    "Produces one subunits.json per stoichiometry combination in the "
    "Cartesian product of the specified ranges."
)


@click.command(
    context_settings={"help_option_names": ["-h", "--help"]},
    epilog=EPILOG,
)
@click.option(
    "--fasta", "-f",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Input FASTA file. Headers are arbitrary identifiers used by --stoichiometry.",
)
@click.option(
    "--stoichiometry", "-s",
    type=str,
    required=True,
    help='Per-chain copy-number ranges, e.g. "alpha:1-6,beta:1-4". '
         'Single values allowed: "gamma:3" means gamma:3-3.',
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False),
    required=True,
    help="Output directory for generated JSON files (created if needed).",
)
@click.option(
    "--max-af-size",
    type=int,
    default=0,
    show_default=True,
    help="Max residues per subunit. 0 = no splitting (each chain = one subunit).",
)
@click.option(
    "--max-chains",
    type=int,
    default=26,
    show_default=True,
    help="Skip combinations whose total chain count exceeds this limit "
         "(PDB chain IDs are single letters A-Z, i.e. 26 chains max).",
)
@click.option(
    "--start-res",
    type=int,
    default=1,
    show_default=True,
    help="Starting residue index for the first subunit of each chain.",
)
@click.option(
    "--chain-start",
    type=str,
    default="A",
    show_default=True,
    help="First chain label. Subsequent labels are sequential (A, B, …, Z, AA, AB, …).",
)
@click.option(
    "--filename-pattern",
    type=str,
    default="subunits_{index:03d}_{stoich}.json",
    show_default=True,
    help="Output filename template. Variables: {index} (0-based), {stoich} (compact stoich string).",
)
@click.option(
    "--manifest/--no-manifest",
    default=True,
    show_default=True,
    help="Write a manifest TSV mapping each output file to its stoichiometry.",
)
@click.option(
    "--manifest-filename",
    type=str,
    default="manifest.tsv",
    show_default=True,
    help="Manifest filename (written to --output-dir).",
)
@click.option(
    "--validate/--no-validate",
    default=True,
    show_default=True,
    help="Validate each generated JSON against CombFold constraints before writing.",
)
@click.option(
    "--max-combinations",
    type=int,
    default=10000,
    show_default=True,
    help="Safety cap on the number of stoichiometry combinations.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Print planned files and stoichiometries without writing anything.",
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    default=False,
    help="Print per-file details during generation.",
)
def cli(
    fasta: str,
    stoichiometry: str,
    output_dir: str,
    max_af_size: int,
    max_chains: int,
    start_res: int,
    chain_start: str,
    filename_pattern: str,
    manifest: bool,
    manifest_filename: str,
    validate: bool,
    max_combinations: int,
    dry_run: bool,
    verbose: bool,
):
    """Generate CombFold subunits.json files for a stoichiometric screen.

    Takes a FASTA of unique sequences and a per-chain copy-number range,
    then writes one CombFold-compatible subunits.json per stoichiometry
    combination.
    """
    # ── Validate chain-start ──
    chain_start = chain_start.upper().strip()
    if not _LABEL_RE.match(chain_start):
        raise click.UsageError(
            f"--chain-start must be uppercase letters only, got '{chain_start}'"
        )

    # ── Parse FASTA ──
    try:
        sequences = parse_fasta(fasta)
    except ValueError as exc:
        raise click.UsageError(f"FASTA parsing error: {exc}")

    fasta_headers = list(sequences.keys())
    click.echo(f"Parsed {len(fasta_headers)} sequence(s) from {fasta}")

    # ── Parse & expand stoichiometry ──
    try:
        ranges = parse_stoichiometry(stoichiometry, fasta_headers)
    except ValueError as exc:
        raise click.UsageError(f"Stoichiometry parsing error: {exc}")

    try:
        combinations = expand_stoichiometry(
            ranges, fasta_headers, max_combinations
        )
    except ValueError as exc:
        raise click.UsageError(str(exc))

    click.echo(f"Expanded to {len(combinations)} stoichiometry combination(s)")

    if dry_run:
        click.echo("\nDry run — no files will be written:\n")
        n_overlimit = 0
        for idx, combo in enumerate(combinations):
            s_str = stoich_string(combo)
            fname = format_filename(filename_pattern, idx, s_str)
            total_chains = sum(c for _, c in combo)
            if total_chains > max_chains:
                n_overlimit += 1
                click.echo(
                    f"  [{idx:3d}] SKIP {fname}  "
                    f"chains={total_chains} > --max-chains {max_chains}  "
                    f"stoich={s_str}"
                )
                continue
            subunits = build_subunits(
                combo, sequences, max_af_size, start_res, chain_start
            )
            click.echo(
                f"  [{idx:3d}] {fname}  "
                f"chains={total_chains}  subunits={len(subunits)}  "
                f"stoich={s_str}"
            )
        click.echo(
            f"\nTotal: {len(combinations) - n_overlimit} file(s) "
            f"({n_overlimit} skipped, >{max_chains} chains)"
        )
        return

    # ── Create output directory ──
    os.makedirs(output_dir, exist_ok=True)

    # ── Generate JSON files ──
    manifest_rows: List[dict] = []
    written = 0
    skipped = 0
    skipped_overlimit = 0

    for idx, combo in enumerate(combinations):
        s_str = stoich_string(combo)
        fname = format_filename(filename_pattern, idx, s_str)
        out_path = os.path.join(output_dir, fname)

        total_chains = sum(c for _, c in combo)
        if total_chains > max_chains:
            click.echo(
                f"  SKIP {fname}: {total_chains} chains > "
                f"--max-chains {max_chains}",
                err=True,
            )
            skipped_overlimit += 1
            continue

        try:
            subunits = build_subunits(
                combo, sequences, max_af_size, start_res, chain_start
            )
        except ValueError as exc:
            click.echo(f"  ERROR building {fname}: {exc}", err=True)
            skipped += 1
            continue

        if validate:
            errors = validate_subunits(subunits)
            if errors:
                click.echo(f"  VALIDATION FAILED for {fname}:", err=True)
                for err in errors:
                    click.echo(f"    - {err}", err=True)
                skipped += 1
                continue

        write_subunits_json(subunits, out_path)
        written += 1

        row = {
            "index": idx,
            "filename": fname,
            "total_chains": total_chains,
            "total_subunits": len(subunits),
        }
        for header, count in combo:
            row[header] = count
        manifest_rows.append(row)

        if verbose:
            click.echo(
                f"  [{idx:3d}] wrote {fname}  "
                f"chains={total_chains}  subunits={len(subunits)}  "
                f"stoich={s_str}"
            )

    # ── Write manifest ──
    if manifest and manifest_rows:
        manifest_path = os.path.join(output_dir, manifest_filename)
        write_manifest(manifest_rows, manifest_path, fasta_headers)
        click.echo(f"\nWrote manifest: {manifest_path}")

    # ── Summary ──
    click.echo(
        f"\nDone: {written} file(s) written, "
        f"{skipped_overlimit} skipped (>{max_chains} chains), "
        f"{skipped} failed"
    )
    if skipped:
        click.echo("Some files failed — see errors above.", err=True)
        sys.exit(1)
    if written == 0:
        click.echo(
            "ERROR: no JSON files were written — check --stoichiometry "
            "and --max-chains.",
            err=True,
        )
        sys.exit(1)


if __name__ == "__main__":
    cli()
