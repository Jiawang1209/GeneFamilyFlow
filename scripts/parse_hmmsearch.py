#!/usr/bin/env python3
"""Parse HMMER domtblout output and extract gene IDs passing an E-value threshold.

Outputs a sorted, deduplicated list of target names (one per line).
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass(frozen=True)
class DomtbloutHit:
    target_name: str
    target_accession: str
    query_name: str
    query_accession: str
    full_evalue: float
    full_score: float
    full_bias: float
    domain_evalue: float
    domain_score: float
    domain_bias: float
    description: str


def parse_domtblout(path: str | Path) -> Iterator[DomtbloutHit]:
    """Parse HMMER domtblout file, yielding one hit per line.

    The domtblout format has 22 fixed whitespace-delimited fields,
    followed by a free-text description.
    """
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    with open(file_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            try:
                yield DomtbloutHit(
                    target_name=fields[0],
                    target_accession=fields[1],
                    query_name=fields[3],
                    query_accession=fields[4],
                    full_evalue=float(fields[6]),
                    full_score=float(fields[7]),
                    full_bias=float(fields[8]),
                    domain_evalue=float(fields[12]),
                    domain_score=float(fields[13]),
                    domain_bias=float(fields[14]),
                    description=" ".join(fields[22:]),
                )
            except (ValueError, IndexError):
                continue


def extract_ids(
    paths: list[str | Path],
    evalue: float,
    use_domain_evalue: bool = False,
) -> list[str]:
    """Extract unique target names from one or more domtblout files."""
    ids: set[str] = set()
    for path in paths:
        for hit in parse_domtblout(path):
            threshold = hit.domain_evalue if use_domain_evalue else hit.full_evalue
            if threshold <= evalue:
                ids.add(hit.target_name)
    return sorted(ids)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract gene IDs from HMMER domtblout files by E-value threshold."
    )
    parser.add_argument(
        "domtblout",
        nargs="+",
        help="One or more HMMER domtblout files.",
    )
    parser.add_argument(
        "-e", "--evalue",
        type=float,
        default=1e-5,
        help="E-value threshold (default: 1e-5).",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output file (default: stdout).",
    )
    parser.add_argument(
        "--domain-evalue",
        action="store_true",
        help="Use per-domain E-value instead of full-sequence E-value.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    for path in args.domtblout:
        if not Path(path).exists():
            print(f"Error: file not found: {path}", file=sys.stderr)
            return 1

    ids = extract_ids(args.domtblout, args.evalue, args.domain_evalue)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as fh:
            fh.write("\n".join(ids) + "\n" if ids else "")
    else:
        for gene_id in ids:
            print(gene_id)

    print(f"Extracted {len(ids)} unique IDs", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
