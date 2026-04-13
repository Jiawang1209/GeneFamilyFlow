#!/usr/bin/env python3
"""Parse pfam_scan.pl output and extract gene IDs matching a target Pfam domain.

pfam_scan.pl output format (space-delimited, comment lines start with #):
  seq_id  ali_start  ali_end  env_start  env_end  pfam_acc  pfam_name  type
  hmm_start  hmm_end  hmm_length  bit_score  e_value  significance  clan

This script filters rows by Pfam accession (e.g. PF00854) and outputs
unique sequence IDs.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal

DomainCombineMode = Literal["any", "all"]


@dataclass(frozen=True)
class PfamScanHit:
    seq_id: str
    ali_start: int
    ali_end: int
    env_start: int
    env_end: int
    pfam_acc: str
    pfam_name: str
    hit_type: str
    bit_score: float
    evalue: float
    clan: str


def parse_pfam_scan(path: str | Path) -> Iterator[PfamScanHit]:
    """Parse pfam_scan.pl output file."""
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    with open(file_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 15:
                continue
            try:
                yield PfamScanHit(
                    seq_id=fields[0],
                    ali_start=int(fields[1]),
                    ali_end=int(fields[2]),
                    env_start=int(fields[3]),
                    env_end=int(fields[4]),
                    pfam_acc=fields[5],
                    pfam_name=fields[6],
                    hit_type=fields[7],
                    bit_score=float(fields[11]),
                    evalue=float(fields[12]),
                    clan=fields[14],
                )
            except (ValueError, IndexError):
                continue


def extract_ids_by_pfam(
    path: str | Path,
    pfam_id: str,
    evalue: float | None = None,
    mode: DomainCombineMode = "any",
) -> list[str]:
    """Extract unique seq IDs matching Pfam accessions under a combine mode.

    Args:
        path: Path to pfam_scan.pl output.
        pfam_id: Comma-separated Pfam accessions (e.g. "PF00854,PF00005").
            Each matches against pfam_acc with or without version suffix.
        evalue: Optional E-value threshold (applied per-hit, before aggregation).
        mode: Multi-domain combination strategy.
            - "any" (default): keep a gene if it has at least one matching domain
              (union semantics — suitable for families with variable architecture).
            - "all": keep a gene only if it has every listed domain
              (intersection semantics — suitable for strict multi-domain families).
            When only one domain is listed, both modes produce the same result.
    """
    if mode not in ("any", "all"):
        raise ValueError(f"mode must be 'any' or 'all', got: {mode!r}")

    pfam_ids = [p.strip() for p in pfam_id.split(",") if p.strip()]
    patterns = [re.compile(rf"^{re.escape(pid)}(\.\d+)?$") for pid in pfam_ids]
    n_patterns = len(patterns)

    gene_matched: dict[str, set[int]] = {}
    for hit in parse_pfam_scan(path):
        if evalue is not None and hit.evalue > evalue:
            continue
        for idx, pattern in enumerate(patterns):
            if pattern.match(hit.pfam_acc):
                gene_matched.setdefault(hit.seq_id, set()).add(idx)
                break

    if mode == "all":
        ids = {gene for gene, matched in gene_matched.items() if len(matched) == n_patterns}
    else:
        ids = set(gene_matched.keys())

    return sorted(ids)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract gene IDs from pfam_scan.pl output by Pfam domain."
    )
    parser.add_argument(
        "pfam_scan_output",
        help="Path to pfam_scan.pl output file.",
    )
    parser.add_argument(
        "--pfam-id",
        required=True,
        help="Pfam accession to filter (e.g. PF00854).",
    )
    parser.add_argument(
        "-e", "--evalue",
        type=float,
        default=None,
        help="E-value threshold (default: no filter).",
    )
    parser.add_argument(
        "--mode",
        choices=("any", "all"),
        default="any",
        help=(
            "Multi-domain combine mode. 'any' keeps genes with at least one "
            "matching domain (union); 'all' requires every listed domain "
            "(intersection). Default: any."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output file (default: stdout).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    pfam_path = Path(args.pfam_scan_output)
    if not pfam_path.exists():
        print(f"Error: file not found: {pfam_path}", file=sys.stderr)
        return 1

    try:
        ids = extract_ids_by_pfam(
            pfam_path, args.pfam_id, args.evalue, mode=args.mode
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as fh:
            fh.write("\n".join(ids) + "\n" if ids else "")
    else:
        for gene_id in ids:
            print(gene_id)

    print(f"Extracted {len(ids)} unique IDs for {args.pfam_id}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
