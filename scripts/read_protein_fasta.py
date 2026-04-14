#!/usr/bin/env python3

from __future__ import annotations

import argparse
import statistics
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

try:
    from scripts._parsers import FastaFormatError, stream_fasta
except ImportError:  # Running as ``python3 scripts/read_protein_fasta.py``
    from _parsers import FastaFormatError, stream_fasta  # type: ignore[no-redef]


@dataclass(frozen=True)
class ProteinRecord:
    record_id: str
    description: str
    sequence: str


def parse_fasta(path: str | Path) -> list[ProteinRecord]:
    return [
        ProteinRecord(
            record_id=raw.record_id,
            description=raw.description,
            sequence=raw.sequence,
        )
        for raw in stream_fasta(path, error_cls=FastaFormatError)
    ]


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


if __name__ == "__main__":
    raise SystemExit(main())
