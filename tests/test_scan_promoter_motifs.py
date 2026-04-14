"""Tests for scripts/scan_promoter_motifs.py.

Covers:
- Literal-match scanning with overlapping hits and forward/reverse strands.
- Palindromic motifs reported once only.
- FASTA header normalization: strips trailing ``(+)``/``(-)`` strand
  suffixes and collapses multi-token headers so R's trailing-ID regex
  still resolves the gene ID.
- Motif library reader rejects a malformed header.
- CLI round-trip produces an 8-column PlantCARE-compatible .tab.
- Output round-trips through the shipped fixture library against the
  shipped Sb.promoter.fasta (sanity check on real data).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.build_plantcare_motifs import Motif as BuilderMotif
from scripts.build_plantcare_motifs import write_tsv as write_motif_tsv
from scripts.scan_promoter_motifs import (
    Motif,
    _normalize_raw_id,
    main,
    read_fasta,
    read_motif_library,
    reverse_complement,
    scan_record,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parent.parent
SHIPPED_MOTIFS = REPO_ROOT / "example" / "10.promoter" / "plantcare_motifs.tsv"
SHIPPED_FASTA = REPO_ROOT / "example" / "10.promoter" / "Sb.promoter.fasta"


def _write_motifs(tmp_path: Path, entries: list[tuple[str, str, str, str, str]]) -> Path:
    p = tmp_path / "motifs.tsv"
    write_motif_tsv(
        [BuilderMotif(*e) for e in entries],
        p,
    )
    return p


def _write_fasta(tmp_path: Path, records: list[tuple[str, str]]) -> Path:
    p = tmp_path / "in.fa"
    with p.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")
    return p


class TestReverseComplement:
    @pytest.mark.parametrize(
        "seq,expected",
        [
            ("A", "T"),
            ("ACGT", "ACGT"),      # palindrome
            ("CACGTG", "CACGTG"),  # palindrome
            ("ACGTG", "CACGT"),
            ("AAANTTT", "AAANTTT"),
        ],
    )
    def test_basic(self, seq: str, expected: str) -> None:
        assert reverse_complement(seq) == expected


class TestNormalizeRawId:
    @pytest.mark.parametrize(
        "header,expected",
        [
            (">Sobic.001G032000.v3.1", "Sobic.001G032000.v3.1"),
            (
                ">Chr01_2434483-2436725:+_usf:2000 Sobic.001G032000.v3.1",
                "Chr01_2434483-2436725:+_usf:2000Sobic.001G032000.v3.1",
            ),
            (">Sobic.001G032000.v3.1(+)", "Sobic.001G032000.v3.1"),
            (">Sobic.001G032000.v3.1(-)", "Sobic.001G032000.v3.1"),
        ],
    )
    def test_formats(self, header: str, expected: str) -> None:
        assert _normalize_raw_id(header) == expected


class TestReadFasta:
    def test_multirecord_with_multiline_sequence(self, tmp_path: Path) -> None:
        p = tmp_path / "in.fa"
        p.write_text(
            ">gene1\nACGT\nACGT\n>gene2(+)\nTTTT\n"
        )
        records = list(read_fasta(p))
        assert len(records) == 2
        assert records[0].raw_id == "gene1"
        assert records[0].sequence == "ACGTACGT"
        assert records[1].raw_id == "gene2"  # strand stripped


class TestScanRecord:
    def test_forward_only_palindrome(self, tmp_path: Path) -> None:
        # CACGTG is a palindrome → only forward strand hits
        motifs = [Motif("G-box", "CACGTG", "Ath", "light")]
        fa = _write_fasta(tmp_path, [("gene1", "AAACACGTGTTT")])
        record = list(read_fasta(fa))[0]
        hits = list(scan_record(record, motifs))
        assert len(hits) == 1
        assert hits[0].strand == "+"
        assert hits[0].position == 4  # 1-based

    def test_forward_and_reverse(self, tmp_path: Path) -> None:
        motifs = [Motif("ABRE", "ACGTG", "Ath", "abscisic")]
        fa = _write_fasta(tmp_path, [("gene1", "ACGTGnnnCACGT")])
        record = list(read_fasta(fa))[0]
        hits = list(scan_record(record, motifs))
        strands = sorted(h.strand for h in hits)
        assert strands == ["+", "-"]
        fwd = next(h for h in hits if h.strand == "+")
        rev = next(h for h in hits if h.strand == "-")
        assert fwd.position == 1 and fwd.length == 5
        assert rev.position == 9 and rev.length == 5

    def test_overlapping_matches(self, tmp_path: Path) -> None:
        # "AAAA" in "AAAAA" has two overlapping matches (pos 1 and 2)
        motifs = [Motif("A-box", "AAAA", "X", "y")]
        fa = _write_fasta(tmp_path, [("g1", "AAAAA")])
        record = list(read_fasta(fa))[0]
        hits = list(scan_record(record, motifs))
        positions = sorted(h.position for h in hits if h.strand == "+")
        assert positions == [1, 2]

    def test_no_match_yields_no_hits(self, tmp_path: Path) -> None:
        motifs = [Motif("X", "GGGGGGG", "a", "b")]
        fa = _write_fasta(tmp_path, [("g1", "AAAAAAAAAA")])
        record = list(read_fasta(fa))[0]
        assert list(scan_record(record, motifs)) == []


class TestReadMotifLibrary:
    def test_roundtrip(self, tmp_path: Path) -> None:
        p = _write_motifs(
            tmp_path,
            [("ABRE", "ACGTG", "Hormone", "Ath", "abscisic")],
        )
        motifs = read_motif_library(p)
        assert len(motifs) == 1
        assert motifs[0].sequence == "ACGTG"

    def test_rejects_bad_header(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.tsv"
        p.write_text("foo\tbar\nABRE\tACGTG\n")
        with pytest.raises(ValueError, match="header"):
            read_motif_library(p)


class TestCLI:
    def test_roundtrip(self, tmp_path: Path) -> None:
        motifs_path = _write_motifs(
            tmp_path,
            [
                ("ABRE", "ACGTG", "Hormone", "Ath", "abscisic"),
                ("CAAT-box", "CAAAT", "Promoter", "Psat", "common"),
            ],
        )
        fa = _write_fasta(
            tmp_path,
            [
                ("Sobic.001G000001.v1", "ACGTGaaaCAAAT"),
                ("Sobic.001G000002.v1", "TTTTT"),
            ],
        )
        out = tmp_path / "out.tab"
        rc = main([
            "--fasta", str(fa),
            "--motifs", str(motifs_path),
            "-o", str(out),
        ])
        assert rc == 0

        lines = out.read_text().splitlines()
        assert lines, "expected at least one hit"
        for line in lines:
            parts = line.split("\t")
            assert len(parts) == 8, f"expected 8 columns, got {len(parts)}: {parts}"
            assert parts[5] in {"+", "-"}
            int(parts[3])  # position parses as int
            int(parts[4])  # length parses as int

        elements = {line.split("\t")[1] for line in lines}
        assert "ABRE" in elements
        assert "CAAT-box" in elements


class TestShippedFixtureRoundtrip:
    def test_scanner_matches_real_promoter_fasta(self, tmp_path: Path) -> None:
        assert SHIPPED_FASTA.exists(), f"missing {SHIPPED_FASTA}"
        assert SHIPPED_MOTIFS.exists(), f"missing {SHIPPED_MOTIFS}"
        out = tmp_path / "out.tab"
        rc = main([
            "--fasta", str(SHIPPED_FASTA),
            "--motifs", str(SHIPPED_MOTIFS),
            "-o", str(out),
        ])
        assert rc == 0
        lines = out.read_text().splitlines()
        assert len(lines) > 1000, "expected thousands of hits on the real fixture"
        # Every row must be 8 columns and end with a gene-ID-looking raw_id
        # so R's trailing-[A-Za-z]+ regex can extract the gene.
        for line in lines[:50]:
            parts = line.split("\t")
            assert len(parts) == 8
            assert parts[0][-1].isalnum()
