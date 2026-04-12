#!/usr/bin/env python3
"""Parse MEME text output to extract motif info and per-gene motif locations.

Outputs two TSV files:
  - meme_info.txt:     Motif_Name  Consensus  Width
  - meme_location.txt: Gene  Start  End  Motif
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MotifInfo:
    name: str
    consensus: str
    width: int


@dataclass(frozen=True)
class MotifLocation:
    gene: str
    start: int
    end: int
    motif: str


# Matches "Motif <consensus> MEME-N Description"
_RE_MOTIF_DESC = re.compile(r"Motif\s+(\S+)\s+(MEME-\d+)\s+Description")

# Matches gene location lines in "sites sorted by position p-value" sections.
# Format: <gene_id>  <start>  <p-value>  <flanking>  <site_seq>  <flanking>
# Gene IDs start with a letter, contain alphanumerics, dots, underscores.
_RE_SITE_LINE = re.compile(
    r"^([A-Za-z]\S+)\s+(\d+)\s+[\d.eE+-]+\s+\S+\s+(\S+)"
)


def parse_meme(
    meme_path: str | Path,
    gene_id_pattern: str | None = None,
) -> tuple[list[MotifInfo], list[MotifLocation]]:
    """Parse MEME text output file.

    Args:
        meme_path: Path to meme.txt.
        gene_id_pattern: Optional regex to filter gene IDs (e.g. ``Sobic\\.``).
            If None, all gene lines are included.

    Returns:
        Tuple of (motif_info_list, motif_location_list).
    """
    gene_filter = re.compile(gene_id_pattern) if gene_id_pattern else None
    motifs: list[MotifInfo] = []
    locations: list[MotifLocation] = []
    current_motif: str | None = None

    with open(meme_path) as fh:
        for line in fh:
            stripped = line.strip()

            desc_match = _RE_MOTIF_DESC.search(stripped)
            if desc_match:
                consensus = desc_match.group(1)
                motif_name = desc_match.group(2)
                motifs.append(MotifInfo(
                    name=motif_name,
                    consensus=consensus,
                    width=len(consensus),
                ))
                current_motif = motif_name
                continue

            if current_motif is None:
                continue

            site_match = _RE_SITE_LINE.match(stripped)
            if site_match:
                gene_id = site_match.group(1)
                if gene_filter and not gene_filter.search(gene_id):
                    continue
                start = int(site_match.group(2))
                site_seq = site_match.group(3)
                end = start + len(site_seq) - 1
                locations.append(MotifLocation(
                    gene=gene_id,
                    start=start,
                    end=end,
                    motif=current_motif,
                ))

    return motifs, locations


def write_outputs(
    motifs: list[MotifInfo],
    locations: list[MotifLocation],
    outdir: Path,
    info_name: str = "meme_info.txt",
    location_name: str = "meme_location.txt",
) -> tuple[Path, Path]:
    """Write motif info and location files."""
    outdir.mkdir(parents=True, exist_ok=True)

    info_path = outdir / info_name
    with open(info_path, "w") as fh:
        for m in motifs:
            fh.write(f"{m.name}\t{m.consensus}\t{m.width}\n")

    loc_path = outdir / location_name
    with open(loc_path, "w") as fh:
        for loc in locations:
            fh.write(f"{loc.gene}\t{loc.start}\t{loc.end}\t{loc.motif}\n")

    return info_path, loc_path


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse MEME text output to extract motif info and locations."
    )
    parser.add_argument(
        "meme_file",
        help="Path to MEME text output file (meme.txt).",
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        help="Output directory (default: current directory).",
    )
    parser.add_argument(
        "--gene-id-pattern",
        default=None,
        help="Regex to filter gene IDs (e.g. 'Sobic\\.'). Default: keep all.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    meme_path = Path(args.meme_file)
    if not meme_path.exists():
        print(f"Error: file not found: {meme_path}", file=sys.stderr)
        return 1

    try:
        motifs, locations = parse_meme(meme_path, args.gene_id_pattern)
    except Exception as exc:
        print(f"Error parsing MEME file: {exc}", file=sys.stderr)
        return 1

    outdir = Path(args.outdir)
    info_path, loc_path = write_outputs(motifs, locations, outdir)

    print(f"Motifs: {len(motifs)}, Locations: {len(locations)}")
    print(f"  → {info_path}")
    print(f"  → {loc_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
