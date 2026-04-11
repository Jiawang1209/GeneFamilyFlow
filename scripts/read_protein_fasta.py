#!/usr/bin/env python3

from __future__ import annotations

import argparse
import statistics
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


class FastaFormatError(ValueError):
    """Raised when a FASTA file cannot be parsed safely."""


@dataclass(frozen=True)
class ProteinRecord:
    record_id: str
    description: str
    sequence: str


def parse_fasta(path: str | Path) -> list[ProteinRecord]:
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    if not file_path.is_file():
        raise FileNotFoundError(f"Path is not a file: {file_path}")

    records: list[ProteinRecord] = []
    current_header: str | None = None
    current_sequence: list[str] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    records.append(_build_record(current_header, current_sequence, line_number - 1))
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

    records.append(_build_record(current_header, current_sequence, line_number))
    return records


def summarize_records(records: list[ProteinRecord]) -> dict[str, object]:
    if not records:
        raise ValueError("No protein records were provided.")

    lengths = [len(record.sequence) for record in records]
    return {
        "sequence_count": len(records),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": statistics.mean(lengths),
        "median_length": statistics.median(lengths),
        "length_distribution": dict(sorted(Counter(lengths).items())),
    }


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Read a protein FASTA file and report basic statistics.")
    parser.add_argument("input_fasta", help="Path to a FASTA-format protein sequence file.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    try:
        records = parse_fasta(args.input_fasta)
        summary = summarize_records(records)
    except (FileNotFoundError, FastaFormatError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Input file: {Path(args.input_fasta)}")
    print(f"Sequence count: {summary['sequence_count']}")
    print(f"Min length: {summary['min_length']}")
    print(f"Max length: {summary['max_length']}")
    print(f"Mean length: {summary['mean_length']:.2f}")
    print(f"Median length: {summary['median_length']}")
    print("Length distribution:")
    for length, count in summary["length_distribution"].items():
        print(f"  {length}: {count}")

    print("\nParsed records:")
    for record in records:
        print(f"- {record.record_id}\t{record.sequence}")

    return 0


def _build_record(header: str, sequence_lines: list[str], line_number: int) -> ProteinRecord:
    if not sequence_lines:
        raise FastaFormatError(f"Line {line_number}: record '{header}' has no sequence.")

    record_id = header.split(maxsplit=1)[0]
    sequence = "".join(sequence_lines).replace(" ", "")
    if not sequence:
        raise FastaFormatError(f"Line {line_number}: record '{header}' has an empty sequence.")

    return ProteinRecord(record_id=record_id, description=header, sequence=sequence)


if __name__ == "__main__":
    raise SystemExit(main())
