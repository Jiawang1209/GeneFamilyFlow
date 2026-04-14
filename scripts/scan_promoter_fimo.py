#!/usr/bin/env python3
"""Scan promoter FASTA with FIMO against the JASPAR plants MEME bundle.

Step 10's "jaspar" scan method. Replaces the manual PlantCARE web
submission with a fully offline, position-weight-matrix scan against
JASPAR 2024 CORE plants (~805 plant TF motifs, see
``scripts/fetch_jaspar_plants.py``). The output is converted to a
PlantCARE-compatible 8-column ``.tab`` so ``R/10_promoter.R`` consumes
it without any change.

Pipeline:
    1. Run ``fimo --oc <tmpdir> --thresh <p> <meme> <fasta>``
    2. Read ``fimo.tsv`` (skipping FIMO comment trailers)
    3. Map each FIMO row to PlantCARE columns:

        raw_id          ← sequence_name (normalized: trailing strand
                          suffix and whitespace collapsed so R's
                          trailing-ID regex still resolves the gene ID)
        element         ← motif_alt_id (TF name, e.g. "Dof3"); falls
                          back to motif_id when alt_id is empty.
        sequence        ← matched_sequence
        position        ← start (1-based)
        length          ← stop - start + 1
        strand          ← '+' / '-'
        species         ← "JASPAR Plants 2024"
        function_desc   ← motif_id (matrix accession, e.g. MA0021.1)

The PlantCARE downstream R parser uses columns 1, 2, and the 8th
``function_desc`` field, so this mapping keeps every visualization
working unchanged.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

# Reuse the FASTA header normalization rule from the literal scanner so
# both scan paths produce identical raw_id columns. The tests pin this
# behavior; bypassing it would silently break R's gene-ID extraction.
from scripts.scan_promoter_motifs import _normalize_raw_id

DEFAULT_FIMO_BINARY = "fimo"
DEFAULT_THRESHOLD = 1e-4

# FIMO 5.x output column order in fimo.tsv.
FIMO_COLUMNS = (
    "motif_id",
    "motif_alt_id",
    "sequence_name",
    "start",
    "stop",
    "strand",
    "score",
    "p-value",
    "q-value",
    "matched_sequence",
)
SPECIES_TAG = "JASPAR Plants 2024"


@dataclass(frozen=True)
class FimoHit:
    raw_id: str
    element: str
    sequence: str
    position: int
    length: int
    strand: str
    species: str
    function_desc: str


def run_fimo(
    fasta: Path,
    meme: Path,
    out_dir: Path,
    threshold: float = DEFAULT_THRESHOLD,
    binary: str = DEFAULT_FIMO_BINARY,
) -> Path:
    """Invoke FIMO and return the path to the resulting ``fimo.tsv``.

    Raises:
        FileNotFoundError: when the FIMO binary is not on PATH.
        RuntimeError: when FIMO exits non-zero or fails to produce
            ``fimo.tsv``.
    """
    if shutil.which(binary) is None:
        raise FileNotFoundError(
            f"FIMO binary '{binary}' not found on PATH. Install MEME suite "
            f"(meme>=5.5; already declared in envs/genefamily.yaml)."
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        binary,
        "--oc",
        str(out_dir),
        "--thresh",
        f"{threshold:g}",
        str(meme),
        str(fasta),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"FIMO failed (exit {proc.returncode})\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}"
        )
    tsv = out_dir / "fimo.tsv"
    if not tsv.exists():
        raise RuntimeError(f"FIMO did not produce {tsv}")
    return tsv


def parse_fimo_tsv(tsv_path: Path) -> Iterator[FimoHit]:
    """Stream rows from a FIMO 5.x ``fimo.tsv`` and yield ``FimoHit``s.

    Skips:
      - The header row (first line that starts with ``motif_id``).
      - Blank lines and FIMO trailer comments (``# ...``).
    """
    with tsv_path.open() as fh:
        header_seen = False
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                continue
            if not header_seen and line.startswith("motif_id"):
                header_seen = True
                continue
            parts = line.split("\t")
            if len(parts) < len(FIMO_COLUMNS):
                # FIMO occasionally pads trailing empty columns; keep
                # well-formed rows only and drop anything truncated.
                continue
            row = dict(zip(FIMO_COLUMNS, parts))
            try:
                start = int(row["start"])
                stop = int(row["stop"])
            except ValueError:
                continue
            element = row["motif_alt_id"].strip() or row["motif_id"].strip()
            yield FimoHit(
                raw_id=_normalize_raw_id(row["sequence_name"]),
                element=element,
                sequence=row["matched_sequence"],
                position=start,
                length=stop - start + 1,
                strand=row["strand"],
                species=SPECIES_TAG,
                function_desc=row["motif_id"].strip(),
            )


def write_plantcare_tab(hits: Iterable[FimoHit], path: Path) -> int:
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


def scan(
    fasta: Path,
    meme: Path,
    output: Path,
    threshold: float = DEFAULT_THRESHOLD,
    binary: str = DEFAULT_FIMO_BINARY,
    work_dir: Path | None = None,
) -> int:
    """End-to-end: run FIMO, convert, write the PlantCARE-format .tab.

    When ``work_dir`` is ``None`` a temp directory is created and removed
    after the run; pass an explicit path to keep FIMO's raw output for
    debugging.
    """
    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="fimo_scan_") as tmp:
            tsv = run_fimo(fasta, meme, Path(tmp), threshold, binary)
            return write_plantcare_tab(parse_fimo_tsv(tsv), output)
    tsv = run_fimo(fasta, meme, work_dir, threshold, binary)
    return write_plantcare_tab(parse_fimo_tsv(tsv), output)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Promoter FASTA file",
    )
    parser.add_argument(
        "--meme",
        type=Path,
        required=True,
        help="JASPAR plants MEME motif bundle (see fetch_jaspar_plants.py)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output .tab in PlantCARE 8-column format",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=f"FIMO p-value threshold (default: {DEFAULT_THRESHOLD:g})",
    )
    parser.add_argument(
        "--fimo-binary",
        default=DEFAULT_FIMO_BINARY,
        help=f"FIMO executable name or path (default: {DEFAULT_FIMO_BINARY})",
    )
    parser.add_argument(
        "--keep-fimo-dir",
        type=Path,
        default=None,
        help="Optional directory to keep FIMO's raw output (debug)",
    )
    args = parser.parse_args(argv)

    n = scan(
        fasta=args.fasta,
        meme=args.meme,
        output=args.output,
        threshold=args.threshold,
        binary=args.fimo_binary,
        work_dir=args.keep_fimo_dir,
    )
    print(f"[scan_promoter_fimo] wrote {n} hits to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
