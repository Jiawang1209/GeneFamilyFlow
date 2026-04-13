#!/usr/bin/env python3
"""Fetch the full KEGG pathway name table from the public REST endpoint.

Usage:
    python scripts/fetch_kegg_pathway_names.py -o data/kegg_pathway_names.tsv

KEGG returns plain-text lines of the form ``path:map04075<TAB>pathway name``.
We normalize ids to the ``ko`` prefix so they match the output of
``scripts/parse_eggnog.py``.
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

KEGG_LIST_URL = "https://rest.kegg.jp/list/pathway"


def fetch(url: str = KEGG_LIST_URL) -> str:
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read().decode("utf-8")


def parse(text: str) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        term = parts[0].strip()
        name = parts[1].strip()
        if term.startswith("path:"):
            term = term[len("path:") :]
        if term.startswith("map"):
            term = "ko" + term[3:]
        rows.append((term, name))
    return rows


def write(rows: list[tuple[str, str]], out: Path) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("term\tname\n")
        for term, name in rows:
            fh.write(f"{term}\t{name}\n")
    return len(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("data/kegg_pathway_names.tsv"),
        help="Output TSV path",
    )
    args = parser.parse_args(argv)

    try:
        text = fetch()
    except Exception as exc:  # pragma: no cover - network-dependent
        print(f"[fetch_kegg] failed: {exc}", file=sys.stderr)
        return 1
    rows = parse(text)
    n = write(rows, args.output)
    print(f"[fetch_kegg] wrote {n} pathways to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
