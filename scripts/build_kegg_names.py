#!/usr/bin/env python3
"""Build an offline KEGG pathway name table for step13 GO/KEGG enrichment.

The ``R/13_go_kegg.R`` script reads a TSV with two columns (``term``,
``name``) keyed on ``ko`` prefixed pathway IDs and uses it as the
``TERM2NAME`` lookup for ``clusterProfiler::enricher``. Without this
lookup, KEGG enrichment plots fall back to raw IDs like ``ko04075``
instead of readable names like ``Plant hormone signal transduction``.

This script builds that TSV from one of two sources:

1. ``--from-kegg-list FILE``: KEGG pathway list as published by
   ``https://rest.kegg.jp/list/pathway`` (two columns: ``map00010`` and
   human-readable name). Names are normalized to the ``ko`` prefix.

2. ``--from-eggnog ANNOTATIONS``: extract the union of KEGG pathway IDs
   seen in one or more eggNOG-mapper annotation files and emit a scaffold
   TSV with empty name columns for users to fill in. Useful when
   generating a minimal per-dataset name table.

Existing entries in ``--merge`` are preserved, so the builder can also
extend a curated fixture with newly seen pathways.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PathwayName:
    term: str   # canonical ``ko`` prefixed id
    name: str   # human readable pathway name (possibly empty for scaffolds)


def _normalize_term(raw: str) -> str | None:
    """Return ``raw`` converted to the canonical ``ko`` prefixed form.

    Accepts ``map00010`` / ``ko00010`` / ``path:map00010`` / ``path:ko00010``
    and rejects anything else.
    """
    cleaned = raw.strip()
    if cleaned.startswith("path:"):
        cleaned = cleaned[len("path:"):]
    if cleaned.startswith("map"):
        suffix = cleaned[3:]
    elif cleaned.startswith("ko"):
        suffix = cleaned[2:]
    else:
        return None
    if not suffix.isdigit() or len(suffix) != 5:
        return None
    return f"ko{suffix}"


def parse_kegg_list(path: Path) -> list[PathwayName]:
    """Parse a KEGG ``pathway.list`` file into PathwayName entries.

    Each line is tab-separated: ``<id>\\t<name>``. Lines whose id cannot
    be normalized to a 5-digit ko form are skipped.
    """
    out: list[PathwayName] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            term = _normalize_term(parts[0])
            if term is None:
                continue
            name = parts[1].strip()
            if not name:
                continue
            out.append(PathwayName(term=term, name=name))
    return out


def extract_eggnog_terms(paths: Iterable[Path]) -> list[str]:
    """Extract the union of KEGG pathway ko IDs observed in annotation files.

    Tolerates both ``map`` and ``ko`` forms; returns canonical ``ko``
    prefixed IDs in sorted order.
    """
    seen: set[str] = set()
    for path in paths:
        with path.open() as fh:
            header: list[str] | None = None
            for raw in fh:
                line = raw.rstrip("\n")
                if not line:
                    continue
                if line.startswith("##"):
                    continue
                if header is None:
                    header = line.lstrip("#").strip().split("\t")
                    try:
                        idx_kegg = header.index("KEGG_Pathway")
                    except ValueError as exc:
                        raise ValueError(
                            f"Missing KEGG_Pathway column in {path}"
                        ) from exc
                    continue
                fields = line.split("\t")
                if len(fields) <= idx_kegg:
                    continue
                cell = fields[idx_kegg].strip()
                if not cell or cell in {"-", "NA", ""}:
                    continue
                for item in cell.split(","):
                    term = _normalize_term(item)
                    if term is not None:
                        seen.add(term)
    return sorted(seen)


def read_existing(path: Path) -> dict[str, str]:
    """Read an existing ``term\\tname`` TSV. Keeps first-seen entry."""
    out: dict[str, str] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for i, raw in enumerate(fh):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            term, name = parts[0].strip(), parts[1].strip()
            if i == 0 and term.lower() == "term":
                continue
            if term and term not in out:
                out[term] = name
    return out


def merge_entries(
    existing: dict[str, str],
    additions: Iterable[PathwayName],
) -> list[PathwayName]:
    """Merge ``additions`` into ``existing`` without overwriting known names."""
    merged: dict[str, str] = dict(existing)
    for entry in additions:
        if entry.term not in merged or not merged[entry.term]:
            merged[entry.term] = entry.name
    return [PathwayName(term=t, name=merged[t]) for t in sorted(merged)]


def write_tsv(entries: Iterable[PathwayName], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    entries = list(entries)
    with path.open("w") as fh:
        fh.write("term\tname\n")
        for e in entries:
            fh.write(f"{e.term}\t{e.name}\n")
    return len(entries)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--from-kegg-list",
        type=Path,
        metavar="FILE",
        help="KEGG pathway.list TSV (e.g. from https://rest.kegg.jp/list/pathway)",
    )
    src.add_argument(
        "--from-eggnog",
        type=Path,
        nargs="+",
        metavar="ANNOT",
        help="One or more eggNOG .emapper.annotations files; emits scaffold entries",
    )
    parser.add_argument(
        "--merge",
        type=Path,
        default=None,
        metavar="FILE",
        help="Existing kegg_names TSV to merge into (preserves known names)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output TSV path (columns: term, name)",
    )
    args = parser.parse_args(argv)

    if args.from_kegg_list is not None:
        additions = parse_kegg_list(args.from_kegg_list)
    else:
        terms = extract_eggnog_terms(args.from_eggnog)
        additions = [PathwayName(term=t, name="") for t in terms]

    existing = read_existing(args.merge) if args.merge else {}
    merged = merge_entries(existing, additions)
    n = write_tsv(merged, args.output)
    print(
        f"[build_kegg_names] wrote {n} pathways to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
