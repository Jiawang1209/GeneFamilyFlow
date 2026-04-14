#!/usr/bin/env python3
"""Offline promoter cis-element scanner (PlantCARE approximation).

Replaces the manual PlantCARE web submission in Step 10 with a local
literal-motif scan. Reads a promoter FASTA plus a motif library TSV
(built by ``build_plantcare_motifs.py``), finds exact matches of each
motif on both strands of each promoter, and emits a PlantCARE-compatible
``.tab`` file that ``R/10_promoter.R`` can consume without changes.

Output format (tab-separated, 8 columns per PlantCARE download):

    raw_id  element  sequence  position  length  strand  species  function_desc

- ``position`` is 1-based, on the forward strand of the promoter.
- ``strand`` is ``+`` when the motif matches the promoter's forward
  strand and ``-`` when it matches the reverse complement.
- ``raw_id`` preserves the FASTA header's gene-ID suffix so
  ``R/10_promoter.R``'s trailing-ID regex still resolves gene names.

This does *not* reproduce PlantCARE's full database — it only catches
motifs present in the supplied library. It is an offline approximation,
not a drop-in replacement.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_STRAND_SUFFIX_RE = re.compile(r"\([+-]\)$")


@dataclass(frozen=True)
class Motif:
    element: str
    sequence: str
    species: str
    function_desc: str


@dataclass(frozen=True)
class FastaRecord:
    raw_id: str
    sequence: str


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def read_motif_library(path: Path) -> list[Motif]:
    """Read a motif library TSV produced by build_plantcare_motifs.py."""
    out: list[Motif] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        required = ["element", "sequence", "description", "species", "function_desc"]
        if header != required:
            raise ValueError(
                f"{path}: motif library must have header {required}, got {header}"
            )
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 5:
                continue
            element, sequence, _desc, species, function_desc = (p.strip() for p in parts)
            if not element or not sequence:
                continue
            out.append(
                Motif(
                    element=element,
                    sequence=sequence.upper(),
                    species=species,
                    function_desc=function_desc,
                )
            )
    return out


def _normalize_raw_id(header_line: str) -> str:
    """Collapse a FASTA header into a PlantCARE-style raw_id.

    - Strips the leading ``>``
    - Removes bedtools ``-name -s`` trailing strand suffix ``(+)`` / ``(-)``
    - Joins whitespace-separated tokens without spaces so column 1 stays
      a single TSV field and R's trailing-ID regex matches the gene ID.
    """
    s = header_line.lstrip(">").rstrip("\n").rstrip()
    s = _STRAND_SUFFIX_RE.sub("", s)
    return "".join(s.split())


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    """Yield ``FastaRecord`` entries from a plain FASTA file (no biopython)."""
    header: str | None = None
    chunks: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header, "".join(chunks).upper())
                header = _normalize_raw_id(line)
                chunks = []
            else:
                if header is None:
                    continue
                chunks.append(line.strip())
        if header is not None:
            yield FastaRecord(header, "".join(chunks).upper())


def _iter_find(haystack: str, needle: str) -> Iterator[int]:
    """Yield overlapping 0-based start indices of ``needle`` in ``haystack``."""
    if not needle:
        return
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx < 0:
            return
        yield idx
        start = idx + 1


@dataclass(frozen=True)
class MotifHit:
    raw_id: str
    element: str
    sequence: str
    position: int  # 1-based
    length: int
    strand: str    # '+' or '-'
    species: str
    function_desc: str


def scan_record(record: FastaRecord, motifs: Iterable[Motif]) -> Iterator[MotifHit]:
    """Scan a single promoter for each motif, forward and reverse strand."""
    seq = record.sequence
    for motif in motifs:
        fwd = motif.sequence
        rev = reverse_complement(fwd)
        for pos in _iter_find(seq, fwd):
            yield MotifHit(
                raw_id=record.raw_id,
                element=motif.element,
                sequence=fwd,
                position=pos + 1,
                length=len(fwd),
                strand="+",
                species=motif.species,
                function_desc=motif.function_desc,
            )
        if rev == fwd:
            continue
        for pos in _iter_find(seq, rev):
            yield MotifHit(
                raw_id=record.raw_id,
                element=motif.element,
                sequence=fwd,
                position=pos + 1,
                length=len(fwd),
                strand="-",
                species=motif.species,
                function_desc=motif.function_desc,
            )


def write_plantcare_tab(hits: Iterable[MotifHit], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w") as fh:
        for h in hits:
            fh.write(
                f"{h.raw_id}\t{h.element}\t{h.sequence}\t{h.position}"
                f"\t{h.length}\t{h.strand}\t{h.species}\t{h.function_desc}\n"
            )
            n += 1
    return n


def scan(fasta: Path, motif_lib: Path, output: Path) -> int:
    motifs = read_motif_library(motif_lib)
    if not motifs:
        raise ValueError(f"{motif_lib}: empty motif library")

    def _iter_hits() -> Iterator[MotifHit]:
        for rec in read_fasta(fasta):
            yield from scan_record(rec, motifs)

    return write_plantcare_tab(_iter_hits(), output)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Promoter FASTA (headers become PlantCARE raw_id column 1)",
    )
    parser.add_argument(
        "--motifs",
        type=Path,
        required=True,
        help="Motif library TSV built by build_plantcare_motifs.py",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output .tab in PlantCARE 8-column format",
    )
    args = parser.parse_args(argv)

    n = scan(args.fasta, args.motifs, args.output)
    print(
        f"[scan_promoter_motifs] wrote {n} hits to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
