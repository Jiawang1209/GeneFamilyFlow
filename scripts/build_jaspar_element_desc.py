#!/usr/bin/env python3
"""Generate a JASPAR element description ``.xlsx`` for Step 10.

Background
==========

``R/10_promoter.R`` expects an ``element_annotation_file`` whose first two
columns are ``element`` and ``description`` and performs an ``inner_join``
against the second column of the PlantCARE ``.tab`` output (the motif
name). The bundled PlantCARE description file
(``cir_element.desc.20240509.xlsx``) uses PlantCARE's cis-element
nomenclature (e.g. ``ABRE``, ``G-box``) — which is orthogonal to JASPAR's
TF-centric alt_ids (``Dof3``, ``AGL3``, ``MYB15``). As a result, joining
FIMO/JASPAR output against the PlantCARE table yields an empty frame.

This script closes that gap by emitting a parallel description table
keyed on JASPAR motif alt_ids and bucketed by TF family. The output has
the same two-column shape (``element`` + ``description``) so the rest of
the pipeline — including ``R/10_promoter.R`` — consumes it unchanged.

Family assignment
=================

Each motif's alt_id is matched against an ordered list of regex rules.
First match wins. Unmatched names fall into ``Other transcription
factor``. The rules reflect published plant-TF family conventions and
can be refined in one place without touching the rest of the pipeline.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

DEFAULT_BUCKET = "Other transcription factor"

# (compiled regex, bucket description). Order matters: earlier rules win.
FAMILY_RULES: tuple[tuple[re.Pattern[str], str], ...] = tuple(
    (re.compile(pattern, re.IGNORECASE), description)
    for pattern, description in (
        (r"^myb|^myb[-\.]", "MYB transcription factor family"),
        (r"^wrky", "WRKY transcription factor family"),
        (r"^(nac|nam|ataf|cuc)", "NAC transcription factor family"),
        (r"^(erf|dreb|cbf|ereb|rap2|tiny|ap2\b)", "AP2/ERF transcription factor family"),
        (r"^(bhlh|pif|ilr|spt|myc|bh[lg])", "bHLH transcription factor family"),
        (
            r"^(bzip|abf|abi5|areb|hy5|gbf|tga|embp|obf|opaque|bzo|gbc)",
            "bZIP transcription factor family",
        ),
        (r"^dof", "DOF zinc finger family"),
        (r"^(gata|bbx)", "GATA/B-box zinc finger family"),
        (r"^tcp", "TCP transcription factor family"),
        (r"^spl", "SBP/SPL transcription factor family"),
        (r"^hsf", "HSF (heat shock factor) family"),
        (r"^arf", "ARF (auxin response factor) family"),
        (r"^arr", "ARR (type-B response regulator) family"),
        (
            r"^(agl|sep|ag\b|ap1\b|pi\b|shp|flc|soc|svp|def|glo|mads)",
            "MADS-box transcription factor family",
        ),
        (
            r"^(athb|hdzip|hat|zhd|pdf2|gl2|anl2|hb\d|wox)",
            "Homeobox (HD-ZIP/ZHD/HB/WOX) family",
        ),
        (r"^(rve|lhy|cca1)", "MYB-related circadian family"),
        (
            r"^(idd|jkd|zfp|zat|gt\d|azf|stop|wip)",
            "C2H2 zinc finger family",
        ),
        (r"^lfy", "LFY/FLORICAULA family"),
        (r"^grf", "GRF growth-regulating family"),
        (r"^(phl|glk|kan)", "GARP/G2-like family"),
        (r"^camta", "CAMTA (calmodulin-binding TF) family"),
        (r"^pbf", "PBF (Prolamin-box binding) family"),
        (r"^ahl", "AT-hook family"),
        (r"^bpc", "BBR/BPC family"),
        (r"^e2f", "E2F/DP family"),
        (r"^far1", "FAR1 family"),
        # Catch-all for species-specific gene model accessions that
        # JASPAR ships without a TF-family label.
        (
            r"^(at\dg|zm\d|os\d|glyma|bradi|solyc|loc_os|phypadraft|aralydraft|rap\d|pk\d|ma\d)",
            "Species-specific accession / unclassified",
        ),
    )
)


@dataclass(frozen=True)
class JasparMotif:
    """One JASPAR motif entry parsed from a MEME bundle."""

    motif_id: str
    alt_id: str

    @property
    def display_name(self) -> str:
        """Name used as the ``element`` key (matches ``motif_alt_id`` from FIMO)."""

        return self.alt_id or self.motif_id


def parse_meme_motifs(path: str | Path) -> list[JasparMotif]:
    """Read ``MOTIF`` lines from a MEME bundle.

    Each ``MOTIF`` line has the shape ``MOTIF <motif_id> [alt_id]``. When
    the alt_id is missing the motif_id is reused so the output row still
    has a non-empty ``element`` key.
    """

    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"MEME bundle not found: {file_path}")
    if not file_path.is_file():
        raise FileNotFoundError(f"Path is not a file: {file_path}")

    motifs: list[JasparMotif] = []
    with file_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("MOTIF"):
                continue
            parts = line.split(maxsplit=2)
            if len(parts) < 2:
                continue
            motif_id = parts[1].strip()
            alt_id = parts[2].strip() if len(parts) == 3 else ""
            motifs.append(JasparMotif(motif_id=motif_id, alt_id=alt_id))

    return motifs


def classify_family(name: str) -> str:
    """Map a motif alt_id (or motif_id) to a TF family bucket."""

    normalized = name.strip()
    if not normalized:
        return DEFAULT_BUCKET

    for pattern, description in FAMILY_RULES:
        if pattern.match(normalized):
            return description

    return DEFAULT_BUCKET


def build_description_rows(motifs: Iterable[JasparMotif]) -> list[tuple[str, str]]:
    """Produce the final ``(element, description)`` rows.

    Duplicate element names are collapsed to the first occurrence so the
    downstream ``inner_join`` behaves deterministically.
    """

    seen: dict[str, str] = {}
    for motif in motifs:
        element = motif.display_name
        if not element or element in seen:
            continue
        seen[element] = classify_family(element)

    return sorted(seen.items())


def write_xlsx(rows: list[tuple[str, str]], output_path: str | Path) -> Path:
    """Write ``rows`` to a two-column (``element``/``description``) xlsx."""

    try:
        from openpyxl import Workbook
    except ImportError as exc:  # pragma: no cover - openpyxl is declared in envs
        raise RuntimeError(
            "openpyxl is required to write the JASPAR element description xlsx. "
            "Install it via `pip install openpyxl` or use the project conda env."
        ) from exc

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "jaspar_elements"
    sheet.append(["element", "description"])
    for element, description in rows:
        sheet.append([element, description])
    workbook.save(output)
    return output


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a JASPAR element description xlsx for Step 10's jaspar "
            "scan method. Keyed on motif alt_id, bucketed by TF family."
        )
    )
    parser.add_argument(
        "--meme",
        required=True,
        help="Path to the JASPAR 2024 CORE plants MEME bundle.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output .xlsx file.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    try:
        motifs = parse_meme_motifs(args.meme)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if not motifs:
        print(f"Error: no MOTIF entries found in {args.meme}", file=sys.stderr)
        return 1

    rows = build_description_rows(motifs)
    write_xlsx(rows, args.output)
    print(f"Wrote {len(rows)} element descriptions to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
