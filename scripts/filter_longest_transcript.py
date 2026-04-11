#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


class FastaFormatError(ValueError):
    """Raised when a FASTA file cannot be parsed safely."""


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA entry with normalized ID and sequence."""

    record_id: str
    description: str
    sequence: str


@dataclass(frozen=True)
class LongestTranscript:
    """Selection metadata kept per gene during the first streaming pass."""

    gene_id: str
    record_id: str
    description: str
    sequence_length: int


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface definition."""

    parser = argparse.ArgumentParser(
        description=(
            "Filter a FASTA file and keep only the longest transcript for each gene. "
            "Sequence IDs are grouped by gene ID using one or more separators."
        )
    )
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument("output_fasta", help="Path to the output FASTA file.")
    parser.add_argument(
        "-s",
        "--separator",
        action="append",
        default=None,
        help=(
            "Separator used between gene ID and transcript ID in the record ID. "
            "Repeat this option to try multiple separators in order. Default: |"
        ),
    )
    parser.add_argument(
        "--line-width",
        type=int,
        default=60,
        help="Number of sequence characters per output line. Use 0 to write one line per sequence.",
    )
    return parser


def parse_fasta(path: str | Path) -> Iterator[FastaRecord]:
    """Stream FASTA records from disk one at a time.

    This generator avoids loading the full file into memory, which is important
    for large genome or transcript FASTA files.
    """

    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    if not file_path.is_file():
        raise FileNotFoundError(f"Path is not a file: {file_path}")

    current_header: str | None = None
    current_sequence: list[str] = []
    last_line_number = 0

    with file_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            last_line_number = line_number
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    yield _build_record(current_header, current_sequence, line_number - 1)
                current_header = line[1:].strip()
                if not current_header:
                    raise FastaFormatError(f"Line {line_number}: FASTA header is empty.")
                current_sequence = []
                continue

            if current_header is None:
                raise FastaFormatError(
                    f"Line {line_number}: encountered sequence data before the first FASTA header."
                )

            current_sequence.append(line)

    if current_header is None:
        raise FastaFormatError(f"No FASTA records found in: {file_path}")

    yield _build_record(current_header, current_sequence, last_line_number)


def extract_gene_id(record_id: str, separators: list[str]) -> tuple[str, bool]:
    """Extract the gene ID from a FASTA record ID.

    Returns:
        (gene_id, used_fallback)

    In lenient mode, if no separator is found, the full record ID is treated as
    a unique gene ID and the caller can emit a warning.
    """

    normalized_id = record_id.strip()
    for separator in separators:
        if not separator:
            continue
        if separator in normalized_id:
            gene_id = normalized_id.split(separator, 1)[0].strip()
            if gene_id:
                return gene_id, False
            break

    return normalized_id, True


def collect_longest_transcripts(
    input_fasta: str | Path, separators: list[str]
) -> tuple[dict[str, LongestTranscript], list[str]]:
    """First pass: identify the longest transcript for each gene."""

    longest_by_gene: dict[str, LongestTranscript] = {}
    fallback_headers: list[str] = []

    for record in parse_fasta(input_fasta):
        gene_id, used_fallback = extract_gene_id(record.record_id, separators)
        if used_fallback:
            fallback_headers.append(record.record_id)

        candidate = LongestTranscript(
            gene_id=gene_id,
            record_id=record.record_id,
            description=record.description,
            sequence_length=len(record.sequence),
        )
        current = longest_by_gene.get(gene_id)

        # Keep the first transcript when lengths are identical so output is stable.
        if current is None or candidate.sequence_length > current.sequence_length:
            longest_by_gene[gene_id] = candidate

    return longest_by_gene, fallback_headers


def filter_longest_transcripts(
    input_fasta: str | Path,
    output_fasta: str | Path,
    separators: list[str],
    line_width: int = 60,
) -> dict[str, int]:
    """Second pass: write only the selected longest transcript records."""

    if line_width < 0:
        raise ValueError("line_width must be 0 or a positive integer.")

    longest_by_gene, fallback_headers = collect_longest_transcripts(input_fasta, separators)
    selected_record_ids = {selection.record_id for selection in longest_by_gene.values()}

    input_count = 0
    output_count = 0
    output_path = Path(output_fasta)

    with output_path.open("w", encoding="utf-8") as handle:
        for record in parse_fasta(input_fasta):
            input_count += 1
            if record.record_id not in selected_record_ids:
                continue

            handle.write(f">{record.description}\n")
            _write_wrapped_sequence(handle, record.sequence, line_width)
            output_count += 1

    return {
        "input_records": input_count,
        "output_records": output_count,
        "fallback_records": len(fallback_headers),
    }


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""

    parser = build_argument_parser()
    args = parser.parse_args(argv)
    separators = args.separator or ["|"]

    try:
        summary = filter_longest_transcripts(
            input_fasta=args.input_fasta,
            output_fasta=args.output_fasta,
            separators=separators,
            line_width=args.line_width,
        )
    except (FileNotFoundError, FastaFormatError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Input records: {summary['input_records']}")
    print(f"Output records: {summary['output_records']}")
    print(f"Fallback IDs: {summary['fallback_records']}", file=sys.stderr)
    if summary["fallback_records"]:
        print(
            "Warning: some record IDs did not match the requested separator(s); "
            "their full record IDs were treated as gene IDs.",
            file=sys.stderr,
        )
    return 0


def _build_record(header: str, sequence_lines: list[str], line_number: int) -> FastaRecord:
    """Validate and normalize one FASTA record."""

    if not sequence_lines:
        raise FastaFormatError(f"Line {line_number}: record '{header}' has no sequence.")

    record_id = header.split(maxsplit=1)[0]
    sequence = "".join(sequence_lines).replace(" ", "")
    if not sequence:
        raise FastaFormatError(f"Line {line_number}: record '{header}' has an empty sequence.")

    return FastaRecord(record_id=record_id, description=header, sequence=sequence)


def _write_wrapped_sequence(handle, sequence: str, line_width: int) -> None:
    """Write sequence lines using a configurable output width."""

    if line_width == 0:
        handle.write(f"{sequence}\n")
        return

    for start in range(0, len(sequence), line_width):
        handle.write(f"{sequence[start:start + line_width]}\n")


if __name__ == "__main__":
    raise SystemExit(main())
