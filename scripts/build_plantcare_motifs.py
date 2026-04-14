#!/usr/bin/env python3
"""Build an offline PlantCARE motif library from downloaded PlantCARE output.

The scanner (``scan_promoter_motifs.py``) consumes this TSV to approximate
PlantCARE cis-element annotation without re-submitting to the PSB web
service. The library is keyed on ``(element, sequence)`` pairs and stores
one representative species / function description per pair — mirroring
the content of PlantCARE's ``plantCARE_output_*.tab`` downloads.

Inputs:
- One or more PlantCARE ``.tab`` files (8 columns: raw_id, element,
  sequence, position, length, strand, species, function_desc).
- An element description xlsx (``element`` / ``description`` columns) used
  as a whitelist — rows whose element is not listed are dropped.
  ``description`` becomes the category used by ``R/10_promoter.R``.

Outputs:
A TSV with columns ``element``, ``sequence``, ``description``, ``species``,
``function_desc`` sorted by (description, element, sequence).

``--merge EXISTING.tsv`` preserves prior rows; newly observed pairs are
appended and existing entries are never overwritten — so user edits to
species/function_desc survive re-builds.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

_VALID_BASES = set("ACGTN")


@dataclass(frozen=True)
class Motif:
    element: str
    sequence: str
    description: str  # category (from xlsx)
    species: str
    function_desc: str


def _read_element_descriptions(path: Path) -> dict[str, str]:
    """Read the element→category xlsx.

    Expected columns: ``element`` (motif name) and ``description``
    (category used by the R plotting layer). Returns a dict mapping
    element name to category. Uses ``openpyxl`` directly so this script
    stays dependency-light.
    """
    try:
        import openpyxl  # type: ignore
    except ImportError as exc:  # pragma: no cover - env issue
        raise RuntimeError(
            "openpyxl is required to read the element description xlsx. "
            "Install via `pip install openpyxl`."
        ) from exc

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = iter(ws.values)
    header = next(rows, None)
    if not header or header[0] != "element" or header[1] != "description":
        raise ValueError(
            f"{path}: expected first two columns 'element' and 'description'"
        )
    out: dict[str, str] = {}
    for row in rows:
        if not row:
            continue
        elem = row[0]
        desc = row[1] if len(row) > 1 else None
        if not elem or not desc:
            continue
        out[str(elem).strip()] = str(desc).strip()
    wb.close()
    return out


def parse_plantcare_tab(
    paths: Iterable[Path],
    allowed_elements: dict[str, str],
) -> list[Motif]:
    """Extract unique ``(element, sequence)`` motifs from PlantCARE tabs.

    Rows whose element is missing from ``allowed_elements`` or whose
    sequence contains non-DNA characters are skipped. Sequences are
    upper-cased.
    """
    seen: dict[tuple[str, str], Motif] = {}
    for path in paths:
        with path.open() as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                element = parts[1].strip()
                sequence = parts[2].strip().upper()
                species = parts[6].strip()
                function_desc = parts[7].strip()
                if not element or not sequence:
                    continue
                if element not in allowed_elements:
                    continue
                if any(b not in _VALID_BASES for b in sequence):
                    continue
                if function_desc == "short_function":
                    continue
                key = (element, sequence)
                if key in seen:
                    continue
                seen[key] = Motif(
                    element=element,
                    sequence=sequence,
                    description=allowed_elements[element],
                    species=species,
                    function_desc=function_desc,
                )
    return list(seen.values())


def read_existing(path: Path) -> dict[tuple[str, str], Motif]:
    """Read an existing motif library TSV into a ``(element, seq) → Motif`` dict."""
    out: dict[tuple[str, str], Motif] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        required = ["element", "sequence", "description", "species", "function_desc"]
        if header != required:
            raise ValueError(
                f"{path}: header must be {required}, got {header}"
            )
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 5:
                continue
            m = Motif(*[p.strip() for p in parts])
            if not m.element or not m.sequence:
                continue
            out[(m.element, m.sequence)] = m
    return out


def merge_motifs(
    existing: dict[tuple[str, str], Motif],
    additions: Iterable[Motif],
) -> list[Motif]:
    """Merge ``additions`` into ``existing``; existing rows take precedence."""
    merged: dict[tuple[str, str], Motif] = dict(existing)
    for m in additions:
        merged.setdefault((m.element, m.sequence), m)
    return sorted(
        merged.values(),
        key=lambda m: (m.description, m.element, m.sequence),
    )


def write_tsv(motifs: Iterable[Motif], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(motifs)
    with path.open("w") as fh:
        fh.write("element\tsequence\tdescription\tspecies\tfunction_desc\n")
        for m in rows:
            fh.write(
                f"{m.element}\t{m.sequence}\t{m.description}"
                f"\t{m.species}\t{m.function_desc}\n"
            )
    return len(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--from-plantcare-tab",
        type=Path,
        nargs="+",
        required=True,
        metavar="FILE",
        help="One or more PlantCARE .tab files (8-column PSB download format)",
    )
    parser.add_argument(
        "--elements-xlsx",
        type=Path,
        required=True,
        help="Element description xlsx (element, description columns)",
    )
    parser.add_argument(
        "--merge",
        type=Path,
        default=None,
        help="Existing motif library TSV to merge into (preserves prior rows)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output motif library TSV",
    )
    args = parser.parse_args(argv)

    allowed = _read_element_descriptions(args.elements_xlsx)
    additions = parse_plantcare_tab(args.from_plantcare_tab, allowed)
    existing = read_existing(args.merge) if args.merge else {}
    merged = merge_motifs(existing, additions)
    n = write_tsv(merged, args.output)
    print(
        f"[build_plantcare_motifs] wrote {n} motifs "
        f"({len({m.element for m in merged})} unique elements) to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
