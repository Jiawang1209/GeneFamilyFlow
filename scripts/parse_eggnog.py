#!/usr/bin/env python3
"""Parse eggNOG-mapper annotations into TERM2GENE long-format tables.

Extracts GO and KEGG pathway assignments per gene and writes three outputs:
- TERM2GENE table for GO terms (columns: term, gene)
- TERM2GENE table for KEGG pathways (columns: term, gene)
- Universe file listing every gene that received any annotation

KEGG pathway IDs are normalized to the ``ko`` prefix so that downstream
name lookup uses a single keying convention.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


MISSING = {"-", "", "NA"}
EGGNOG_REQUIRED_COLUMNS = ("query", "GOs", "KEGG_Pathway")


@dataclass(frozen=True)
class EggnogRow:
    gene: str
    go_terms: tuple[str, ...]
    kegg_pathways: tuple[str, ...]


def _iter_data_lines(path: Path) -> Iterator[str]:
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            yield line


def _parse_header(line: str) -> list[str]:
    cleaned = line.lstrip("#").strip()
    return cleaned.split("\t")


def _split_list(cell: str) -> tuple[str, ...]:
    if cell in MISSING:
        return ()
    return tuple(
        item.strip()
        for item in cell.split(",")
        if item.strip() and item.strip() not in MISSING
    )


def _normalize_kegg(pathways: Iterable[str]) -> tuple[str, ...]:
    """Normalize KEGG pathway IDs to the ``ko`` prefix and dedupe.

    eggNOG-mapper typically emits both ``map04075`` and ``ko04075``. We keep
    only the ``ko`` form so downstream name lookup uses one key space.
    """
    seen: set[str] = set()
    out: list[str] = []
    for item in pathways:
        if item.startswith("map"):
            key = "ko" + item[3:]
        elif item.startswith("ko"):
            key = item
        else:
            continue
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return tuple(out)


def parse_eggnog_annotations(path: Path) -> list[EggnogRow]:
    """Parse an eggNOG-mapper ``.emapper.annotations`` TSV file."""
    lines = list(_iter_data_lines(path))
    if not lines:
        raise ValueError(f"No data lines in {path}")

    header = _parse_header(lines[0])
    try:
        idx_query = header.index("query")
        idx_go = header.index("GOs")
        idx_kegg = header.index("KEGG_Pathway")
    except ValueError as exc:
        missing = [c for c in EGGNOG_REQUIRED_COLUMNS if c not in header]
        raise ValueError(
            f"eggNOG annotations header missing columns {missing}: {path}"
        ) from exc

    rows: list[EggnogRow] = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) <= max(idx_query, idx_go, idx_kegg):
            continue
        gene = fields[idx_query].strip()
        if not gene or gene in MISSING:
            continue
        go_terms = _split_list(fields[idx_go])
        kegg_pathways = _normalize_kegg(_split_list(fields[idx_kegg]))
        rows.append(EggnogRow(gene=gene, go_terms=go_terms, kegg_pathways=kegg_pathways))
    return rows


def write_term2gene(rows: Iterable[EggnogRow], out_path: Path, kind: str) -> int:
    assert kind in {"go", "kegg"}
    pairs: list[tuple[str, str]] = []
    for row in rows:
        terms = row.go_terms if kind == "go" else row.kegg_pathways
        for term in terms:
            pairs.append((term, row.gene))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("term\tgene\n")
        for term, gene in pairs:
            fh.write(f"{term}\t{gene}\n")
    return len(pairs)


def write_universe(rows: Iterable[EggnogRow], out_path: Path) -> int:
    rows = list(rows)
    genes: list[str] = []
    seen: set[str] = set()
    for row in rows:
        if not row.go_terms and not row.kegg_pathways:
            continue
        if row.gene in seen:
            continue
        seen.add(row.gene)
        genes.append(row.gene)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for gene in genes:
            fh.write(f"{gene}\n")
    return len(genes)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Parse eggNOG-mapper annotations into TERM2GENE tables."
    )
    parser.add_argument("annotations", type=Path, help="eggNOG-mapper .emapper.annotations TSV")
    parser.add_argument("--go-out", type=Path, required=True, help="Output GO TERM2GENE TSV")
    parser.add_argument("--kegg-out", type=Path, required=True, help="Output KEGG TERM2GENE TSV")
    parser.add_argument(
        "--universe-out",
        type=Path,
        required=True,
        help="Output gene universe (one gene per line)",
    )
    args = parser.parse_args(argv)

    rows = parse_eggnog_annotations(args.annotations)
    n_go = write_term2gene(rows, args.go_out, "go")
    n_kegg = write_term2gene(rows, args.kegg_out, "kegg")
    n_universe = write_universe(rows, args.universe_out)
    print(
        f"[parse_eggnog] {args.annotations.name}: "
        f"{len(rows)} genes parsed, "
        f"{n_go} GO pairs, {n_kegg} KEGG pairs, universe={n_universe}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
