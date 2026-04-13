#!/usr/bin/env python3
"""Extract BLAST seed sequences for a Pfam domain from a reference proteome.

Given a clean 2-column TSV of gene-Pfam hits and a gene-level protein FASTA,
this script writes the subset of sequences whose gene ID has at least one
row matching the requested Pfam accession.

Input formats
-------------
domains_table : TSV, 2 columns, no header, one row per (gene_id, pfam_id)
                lines starting with '#' and blank lines are ignored
                example:
                    AT1G12110\tPF00854
                    AT1G12110\tPF07690
                    AT1G18880\tPF00854

proteome : FASTA, gene-level headers (first token after '>' is the gene ID)
           example:
                >AT1G12110
                MSLPETKSDDIL...

Usage
-----
    python extract_seed_by_pfam.py \
        --domains-table Athaliana.domains.tsv \
        --proteome      Athaliana.pep.fasta \
        --pfam-id       PF00854 \
        --output        work/03_blast/PF00854/seed.fa
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, Iterator


def parse_domain_hits(path: str | Path, pfam_id: str) -> list[str]:
    """Return sorted unique gene IDs that hit `pfam_id` in the 2-col TSV.

    Raises:
        FileNotFoundError: if the path does not exist.
    """
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"domains table not found: {file_path}")

    gene_ids: set[str] = set()
    with open(file_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene, pfam = parts[0], parts[1]
            if pfam == pfam_id:
                gene_ids.add(gene)
    return sorted(gene_ids)


def extract_sequences(
    path: str | Path,
    gene_ids: Iterable[str],
) -> Iterator[tuple[str, str]]:
    """Yield (gene_id, sequence) for entries in `gene_ids` present in FASTA.

    The FASTA must have gene-level headers (first token after '>' is the
    gene ID). Sequences split across multiple lines are joined.

    Missing gene IDs are silently skipped.
    """
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"proteome FASTA not found: {file_path}")

    wanted = set(gene_ids)
    if not wanted:
        return

    current_id: str | None = None
    current_seq: list[str] = []

    def flush() -> Iterator[tuple[str, str]]:
        if current_id is not None and current_id in wanted:
            yield current_id, "".join(current_seq)

    with open(file_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                yield from flush()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        yield from flush()


def write_seed_fasta(
    proteome: str | Path,
    gene_ids: Iterable[str],
    output: str | Path,
) -> int:
    """Write matching sequences to `output` FASTA. Returns count written."""
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with open(out_path, "w") as fh:
        for gene_id, seq in extract_sequences(proteome, gene_ids):
            fh.write(f">{gene_id}\n{seq}\n")
            count += 1
    return count


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract BLAST seed sequences for a Pfam domain from a clean "
            "2-column domains TSV and a gene-level proteome FASTA."
        )
    )
    parser.add_argument(
        "--domains-table",
        required=True,
        help="Path to 2-column TSV (gene_id<TAB>pfam_id).",
    )
    parser.add_argument(
        "--proteome",
        required=True,
        help="Path to gene-level protein FASTA.",
    )
    parser.add_argument(
        "--pfam-id",
        required=True,
        help="Target Pfam accession (e.g. PF00854).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA path.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    try:
        gene_ids = parse_domain_hits(args.domains_table, args.pfam_id)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    try:
        count = write_seed_fasta(args.proteome, gene_ids, args.output)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(
        f"Extracted {count} seed sequences for {args.pfam_id} "
        f"(from {len(gene_ids)} gene IDs)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
