#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path


class Gff3TranscriptSelectionError(ValueError):
    """Raised when a GFF3 file cannot be parsed into gene/transcript selections safely."""


@dataclass(frozen=True)
class Gff3Feature:
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str]

    @property
    def feature_id(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parents(self) -> list[str]:
        parent_value = self.attributes.get("Parent", "")
        if not parent_value:
            return []
        return [parent.strip() for parent in parent_value.split(",") if parent.strip()]


@dataclass(frozen=True)
class LongestTranscriptRecord:
    gene_id: str
    transcript_id: str
    seqid: str
    start: int
    end: int
    strand: str
    cds_length: int


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract the longest CDS-backed transcript for each gene from a GFF3 file "
            "and write the result as CSV."
        )
    )
    parser.add_argument("input_gff3", help="Path to the input GFF3 file.")
    parser.add_argument("output_csv", help="Path to the output CSV file.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if a gene does not have any transcript with CDS features.",
    )
    return parser


def parse_gff3(path: str | Path) -> list[Gff3Feature]:
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {file_path}")
    if not file_path.is_file():
        raise FileNotFoundError(f"Path is not a file: {file_path}")

    features: list[Gff3Feature] = []
    with file_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise Gff3TranscriptSelectionError(
                    f"Line {line_number}: expected 9 columns, found {len(fields)}."
                )

            try:
                feature_start = int(fields[3])
                feature_end = int(fields[4])
            except ValueError as exc:
                raise Gff3TranscriptSelectionError(
                    f"Line {line_number}: start and end must be integers."
                ) from exc

            if feature_start > feature_end:
                raise Gff3TranscriptSelectionError(
                    f"Line {line_number}: start ({feature_start}) cannot be greater than end ({feature_end})."
                )

            features.append(
                Gff3Feature(
                    seqid=fields[0],
                    feature_type=fields[2],
                    start=feature_start,
                    end=feature_end,
                    strand=fields[6],
                    attributes=_parse_attributes(fields[8]),
                )
            )

    if not features:
        raise Gff3TranscriptSelectionError(f"No GFF3 features found in: {file_path}")
    return features


def collect_longest_transcripts(
    path: str | Path, include_genes_without_cds: bool = True
) -> dict[str, LongestTranscriptRecord]:
    features = parse_gff3(path)

    genes: dict[str, Gff3Feature] = {}
    transcripts: dict[str, Gff3Feature] = {}
    transcript_order: list[str] = []
    transcript_to_gene: dict[str, str] = {}
    cds_lengths: dict[str, int] = {}

    for feature in features:
        feature_id = feature.feature_id
        if feature.feature_type == "gene" and feature_id:
            genes[feature_id] = feature
            continue

        if feature.feature_type in {"mRNA", "transcript"} and feature_id:
            transcripts[feature_id] = feature
            transcript_order.append(feature_id)
            gene_id = _resolve_gene_parent(feature, genes)
            if gene_id is not None:
                transcript_to_gene[feature_id] = gene_id
            continue

        if feature.feature_type != "CDS":
            continue

        for transcript_id in feature.parents:
            cds_lengths[transcript_id] = cds_lengths.get(transcript_id, 0) + (
                feature.end - feature.start + 1
            )

    longest_by_gene: dict[str, LongestTranscriptRecord] = {}
    genes_with_cds: set[str] = set()
    for transcript_id in transcript_order:
        gene_id = transcript_to_gene.get(transcript_id)
        cds_length = cds_lengths.get(transcript_id, 0)
        if gene_id is None or cds_length == 0:
            continue

        genes_with_cds.add(gene_id)
        transcript = transcripts[transcript_id]
        candidate = LongestTranscriptRecord(
            gene_id=gene_id,
            transcript_id=transcript_id,
            seqid=transcript.seqid,
            start=transcript.start,
            end=transcript.end,
            strand=transcript.strand,
            cds_length=cds_length,
        )
        current = longest_by_gene.get(gene_id)
        if current is None or candidate.cds_length > current.cds_length:
            longest_by_gene[gene_id] = candidate

    missing_gene_ids = sorted(gene_id for gene_id in genes if gene_id not in genes_with_cds)
    if missing_gene_ids and not include_genes_without_cds:
        raise Gff3TranscriptSelectionError(
            "Genes without CDS-backed transcripts: " + ", ".join(missing_gene_ids)
        )

    return longest_by_gene


def write_longest_transcripts_csv(
    input_gff3: str | Path,
    output_csv: str | Path,
    include_genes_without_cds: bool = True,
) -> dict[str, int]:
    features = parse_gff3(input_gff3)
    gene_count = sum(1 for feature in features if feature.feature_type == "gene" and feature.feature_id)
    transcript_count = sum(
        1 for feature in features if feature.feature_type in {"mRNA", "transcript"} and feature.feature_id
    )
    selections = collect_longest_transcripts(
        input_gff3,
        include_genes_without_cds=include_genes_without_cds,
    )

    output_path = Path(output_csv)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene_id", "transcript_id", "seqid", "start", "end", "strand", "cds_length"])
        for gene_id in sorted(selections):
            selection = selections[gene_id]
            writer.writerow(
                [
                    selection.gene_id,
                    selection.transcript_id,
                    selection.seqid,
                    selection.start,
                    selection.end,
                    selection.strand,
                    selection.cds_length,
                ]
            )

    return {
        "gene_count": len(selections),
        "transcript_count": transcript_count,
        "skipped_genes": gene_count - len(selections),
    }


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    try:
        summary = write_longest_transcripts_csv(
            args.input_gff3,
            args.output_csv,
            include_genes_without_cds=not args.strict,
        )
    except (FileNotFoundError, Gff3TranscriptSelectionError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Genes written: {summary['gene_count']}")
    print(f"Transcripts scanned: {summary['transcript_count']}")
    print(
        f"Skipped genes without CDS-backed transcripts: {summary['skipped_genes']}",
        file=sys.stderr,
    )
    return 0


def _parse_attributes(raw_attributes: str) -> dict[str, str]:
    if raw_attributes == ".":
        return {}

    attributes: dict[str, str] = {}
    for item in raw_attributes.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise Gff3TranscriptSelectionError(f"Malformed GFF3 attribute: {item}")
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            raise Gff3TranscriptSelectionError(f"Malformed GFF3 attribute: {item}")
        attributes[key] = value
    return attributes


def _resolve_gene_parent(feature: Gff3Feature, genes: dict[str, Gff3Feature]) -> str | None:
    for parent_id in feature.parents:
        if parent_id in genes:
            return parent_id
    return None


if __name__ == "__main__":
    raise SystemExit(main())
