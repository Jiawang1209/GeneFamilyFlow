"""Shared FASTA and GFF3 parsing primitives for the ``scripts/`` package.

This module centralises the boilerplate that the individual CLI scripts
used to duplicate. Each target script keeps its own public classes and
function names so that the existing test contracts hold, but their
bodies delegate to the helpers below.

The module is named ``_parsers`` (not ``_io``) to avoid colliding with
CPython's builtin ``_io`` frozen module.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


class FastaFormatError(ValueError):
    """Raised when a FASTA file cannot be parsed safely."""


class Gff3FormatError(ValueError):
    """Raised when a GFF3 file cannot be parsed safely."""


@dataclass(frozen=True)
class RawFastaRecord:
    """A normalized FASTA record produced by :func:`stream_fasta`."""

    record_id: str
    description: str
    sequence: str


@dataclass(frozen=True)
class RawGff3Row:
    """A normalized GFF3 data row produced by :func:`iter_gff3_rows`.

    ``raw_line`` preserves the original (newline-stripped) line so callers
    that need to re-emit the file can do so without rebuilding the tab
    layout. ``strand`` is the value from column 7.
    """

    raw_line: str
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str]


def _open_existing_file(path: str | Path, label: str) -> Path:
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"{label} file not found: {file_path}")
    if not file_path.is_file():
        raise FileNotFoundError(f"Path is not a file: {file_path}")
    return file_path


def stream_fasta(
    path: str | Path,
    *,
    error_cls: type[Exception] = FastaFormatError,
) -> Iterator[RawFastaRecord]:
    """Yield FASTA records from ``path`` one at a time.

    The generator validates headers and sequence blocks incrementally so
    that callers can stream large files without loading them entirely.
    ``error_cls`` lets callers substitute a module-local subclass when
    the test suite distinguishes error types.
    """

    file_path = _open_existing_file(path, "FASTA")

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
                    yield _build_fasta_record(
                        current_header, current_sequence, line_number - 1, error_cls
                    )
                current_header = line[1:].strip()
                if not current_header:
                    raise error_cls(f"Line {line_number}: FASTA header is empty.")
                current_sequence = []
                continue

            if current_header is None:
                raise error_cls(
                    f"Line {line_number}: encountered sequence data before the first FASTA header."
                )

            current_sequence.append(line)

    if current_header is None:
        raise error_cls(f"No FASTA records found in: {file_path}")

    yield _build_fasta_record(current_header, current_sequence, last_line_number, error_cls)


def _build_fasta_record(
    header: str,
    sequence_lines: list[str],
    line_number: int,
    error_cls: type[Exception],
) -> RawFastaRecord:
    if not sequence_lines:
        raise error_cls(f"Line {line_number}: record '{header}' has no sequence.")

    record_id = header.split(maxsplit=1)[0]
    sequence = "".join(sequence_lines).replace(" ", "")
    if not sequence:
        raise error_cls(f"Line {line_number}: record '{header}' has an empty sequence.")

    return RawFastaRecord(record_id=record_id, description=header, sequence=sequence)


def iter_gff3_rows(
    path: str | Path,
    *,
    error_cls: type[Exception] = Gff3FormatError,
    strict_attributes: bool = False,
) -> Iterator[tuple[str, object]]:
    """Yield ``("comment", str)`` or ``("feature", RawGff3Row)`` tuples.

    ``error_cls`` controls which exception class is raised on malformed
    rows or attribute blocks. ``strict_attributes`` toggles the attribute
    parser between lenient (skip items without ``=``) and strict (raise
    on malformed items).
    """

    file_path = _open_existing_file(path, "GFF3")

    with file_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                yield "comment", line
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise error_cls(
                    f"Line {line_number}: expected 9 columns, found {len(fields)}."
                )

            try:
                feature_start = int(fields[3])
                feature_end = int(fields[4])
            except ValueError as exc:
                raise error_cls(
                    f"Line {line_number}: start and end must be integers."
                ) from exc

            if feature_start > feature_end:
                raise error_cls(
                    f"Line {line_number}: start ({feature_start}) cannot be greater than end ({feature_end})."
                )

            attributes = parse_gff3_attributes(
                fields[8], strict=strict_attributes, error_cls=error_cls
            )

            yield "feature", RawGff3Row(
                raw_line=line,
                seqid=fields[0],
                feature_type=fields[2],
                start=feature_start,
                end=feature_end,
                strand=fields[6],
                attributes=attributes,
            )


def parse_gff3_attributes(
    raw_attributes: str,
    *,
    strict: bool = False,
    error_cls: type[Exception] = ValueError,
) -> dict[str, str]:
    """Parse the 9th GFF3 column into an attribute dict.

    ``strict`` mode raises ``error_cls`` for items that lack an ``=``
    separator or whose key strips to an empty string. Lenient mode
    silently skips items without ``=`` but still records entries whose
    key strips to empty (preserving the historical ``filter_gff3``
    behaviour).
    """

    if raw_attributes == ".":
        return {}

    attributes: dict[str, str] = {}
    for item in raw_attributes.split(";"):
        if not item:
            continue
        if "=" not in item:
            if strict:
                raise error_cls(f"Malformed GFF3 attribute: {item}")
            continue

        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key and strict:
            raise error_cls(f"Malformed GFF3 attribute: {item}")
        attributes[key] = value

    return attributes
