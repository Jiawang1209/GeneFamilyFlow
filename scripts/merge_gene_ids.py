#!/usr/bin/env python3
"""Merge two gene ID lists using intersection or union.

Reads two plain-text files (one ID per line) and writes the merged result.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_ids(path: str | Path) -> set[str]:
    """Read a file of IDs (one per line), skipping blanks and comments."""
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    ids: set[str] = set()
    with open(file_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                ids.add(stripped)
    return ids


def merge_ids(
    ids_a: set[str],
    ids_b: set[str],
    method: str = "intersection",
) -> list[str]:
    """Merge two ID sets.

    Args:
        ids_a: First set of IDs.
        ids_b: Second set of IDs.
        method: "intersection" or "union".

    Returns:
        Sorted list of merged IDs.
    """
    if method == "intersection":
        merged = ids_a & ids_b
    elif method == "union":
        merged = ids_a | ids_b
    else:
        raise ValueError(f"Unknown merge method: {method!r}. Use 'intersection' or 'union'.")
    return sorted(merged)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Merge two gene ID lists using intersection or union."
    )
    parser.add_argument("file_a", help="First ID list file.")
    parser.add_argument("file_b", help="Second ID list file.")
    parser.add_argument(
        "-m", "--method",
        choices=["intersection", "union"],
        default="intersection",
        help="Merge method (default: intersection).",
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

    try:
        ids_a = read_ids(args.file_a)
        ids_b = read_ids(args.file_b)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    merged = merge_ids(ids_a, ids_b, args.method)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as fh:
            fh.write("\n".join(merged) + "\n" if merged else "")
    else:
        for gene_id in merged:
            print(gene_id)

    print(
        f"A={len(ids_a)}, B={len(ids_b)}, {args.method}={len(merged)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
