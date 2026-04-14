#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

try:
    from scripts._parsers import Gff3FormatError, iter_gff3_rows, parse_gff3_attributes
except ImportError:  # Running as ``python3 scripts/filter_gff3.py``
    from _parsers import (  # type: ignore[no-redef]
        Gff3FormatError,
        iter_gff3_rows,
        parse_gff3_attributes,
    )


@dataclass(frozen=True)
class Gff3Feature:
    raw_line: str
    seqid: str
    feature_type: str
    start: int
    end: int
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


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Filter a GFF3 file by gene ID, chromosome ID, start, and end coordinates."
    )
    parser.add_argument("input_gff3", help="Path to the input GFF3 file.")
    parser.add_argument(
        "--gene-id",
        action="append",
        default=None,
        help="Gene ID to keep. Repeat or provide comma-separated values.",
    )
    parser.add_argument(
        "--chromosome-id",
        action="append",
        default=None,
        help="Chromosome/seqid to keep. Repeat or provide comma-separated values.",
    )
    parser.add_argument("--start", type=int, default=None, help="Minimum gene start coordinate to keep.")
    parser.add_argument("--end", type=int, default=None, help="Maximum gene end coordinate to keep.")
    return parser


def filter_gff3_lines(
    path: str | Path,
    gene_ids: set[str] | None = None,
    chromosome_ids: set[str] | None = None,
    start: int | None = None,
    end: int | None = None,
) -> list[str]:
    if start is not None and start < 1:
        raise ValueError("start must be greater than or equal to 1.")
    if end is not None and end < 1:
        raise ValueError("end must be greater than or equal to 1.")
    if start is not None and end is not None and start > end:
        raise ValueError("start cannot be greater than end.")

    comments, features = parse_gff3(path)
    if not features:
        return comments

    by_id = {feature.feature_id: feature for feature in features if feature.feature_id}
    root_gene_ids = {
        feature.feature_id
        for feature in features
        if feature.feature_type == "gene" and feature.feature_id
    }

    matched_gene_ids = {
        gene_id
        for gene_id in root_gene_ids
        if _matches_gene_filters(by_id[gene_id], gene_ids, chromosome_ids, start, end)
    }

    if not any([gene_ids, chromosome_ids, start is not None, end is not None]):
        matched_gene_ids = set(root_gene_ids)

    filtered_lines = list(comments)
    for feature in features:
        root_gene_id = _resolve_root_gene_id(feature, by_id)
        if root_gene_id in matched_gene_ids:
            filtered_lines.append(feature.raw_line)

    return filtered_lines


def parse_gff3(path: str | Path) -> tuple[list[str], list[Gff3Feature]]:
    comments: list[str] = []
    features: list[Gff3Feature] = []

    for kind, item in iter_gff3_rows(path, error_cls=Gff3FormatError, strict_attributes=False):
        if kind == "comment":
            comments.append(item)  # type: ignore[arg-type]
            continue

        row = item  # type: ignore[assignment]
        features.append(
            Gff3Feature(
                raw_line=row.raw_line,
                seqid=row.seqid,
                feature_type=row.feature_type,
                start=row.start,
                end=row.end,
                attributes=row.attributes,
            )
        )

    return comments, features


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    try:
        filtered_lines = filter_gff3_lines(
            args.input_gff3,
            gene_ids=_normalize_multi_value_argument(args.gene_id),
            chromosome_ids=_normalize_multi_value_argument(args.chromosome_id),
            start=args.start,
            end=args.end,
        )
    except (FileNotFoundError, Gff3FormatError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if filtered_lines:
        sys.stdout.write("\n".join(filtered_lines) + "\n")
    return 0


def _matches_gene_filters(
    feature: Gff3Feature,
    gene_ids: set[str] | None,
    chromosome_ids: set[str] | None,
    start: int | None,
    end: int | None,
) -> bool:
    gene_identifiers = {
        value
        for value in (
            feature.attributes.get("ID"),
            feature.attributes.get("gene_id"),
            feature.attributes.get("Name"),
        )
        if value
    }

    if gene_ids and gene_identifiers.isdisjoint(gene_ids):
        return False
    if chromosome_ids and feature.seqid not in chromosome_ids:
        return False
    if start is not None and feature.start < start:
        return False
    if end is not None and feature.end > end:
        return False
    return True


def _normalize_multi_value_argument(values: list[str] | None) -> set[str] | None:
    if not values:
        return None

    normalized = {
        item.strip()
        for value in values
        for item in value.split(",")
        if item.strip()
    }
    return normalized or None


def _parse_attributes(raw_attributes: str) -> dict[str, str]:
    """Lenient attribute parser kept for backwards compatibility.

    Delegates to :func:`scripts._parsers.parse_gff3_attributes` in lenient
    mode so the existing test contract (``_parse_attributes("ID=G1;orphan")
    == {"ID": "G1"}``) is preserved.
    """

    return parse_gff3_attributes(raw_attributes, strict=False)


def _resolve_root_gene_id(feature: Gff3Feature, by_id: dict[str, Gff3Feature]) -> str | None:
    feature_id = feature.feature_id
    if feature.feature_type == "gene" and feature_id:
        return feature_id

    current_parents = feature.parents
    visited: set[str] = set()

    while current_parents:
        parent_id = current_parents[0]
        if parent_id in visited:
            return None
        visited.add(parent_id)

        parent = by_id.get(parent_id)
        if parent is None:
            return None
        if parent.feature_type == "gene" and parent.feature_id:
            return parent.feature_id
        current_parents = parent.parents

    return None


if __name__ == "__main__":
    raise SystemExit(main())
