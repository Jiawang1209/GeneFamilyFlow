from __future__ import annotations

import unittest
from pathlib import Path

import pytest

from scripts.read_protein_fasta import (
    FastaFormatError,
    ProteinRecord,
    main,
    parse_fasta,
    summarize_records,
)

pytestmark = pytest.mark.unit


class ParseFastaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture_dir = Path(__file__).parent / "data"

    def test_parse_fasta_extracts_ids_and_sequences(self) -> None:
        records = parse_fasta(self.fixture_dir / "valid_protein.fa")

        self.assertEqual(
            records,
            [
                ProteinRecord(record_id="seq1", description="seq1 description", sequence="MKTIIALSYIFCLVFAD"),
                ProteinRecord(record_id="seq2", description="seq2", sequence="GAVLIPFW"),
            ],
        )

    def test_summarize_records_reports_basic_length_statistics(self) -> None:
        records = [
            ProteinRecord("a", "a", "MKT"),
            ProteinRecord("b", "b", "MKTI"),
            ProteinRecord("c", "c", "M"),
        ]

        summary = summarize_records(records)

        self.assertEqual(summary["sequence_count"], 3)
        self.assertEqual(summary["min_length"], 1)
        self.assertEqual(summary["max_length"], 4)
        self.assertAlmostEqual(summary["mean_length"], 8 / 3, places=6)
        self.assertEqual(summary["median_length"], 3)
        self.assertEqual(summary["length_distribution"], {1: 1, 3: 1, 4: 1})

    def test_parse_fasta_rejects_sequence_before_header(self) -> None:
        with self.assertRaises(FastaFormatError):
            parse_fasta(self.fixture_dir / "bad_before_header.fa")

    def test_parse_fasta_rejects_empty_sequence_entries(self) -> None:
        with self.assertRaises(FastaFormatError):
            parse_fasta(self.fixture_dir / "bad_empty_sequence.fa")

    def test_parse_fasta_rejects_missing_files(self) -> None:
        with self.assertRaises(FileNotFoundError):
            parse_fasta(Path("/tmp/does-not-exist-protein.fa"))


class TestParseFastaErrorBranches:
    def test_rejects_directory_path(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="Path is not a file"):
            parse_fasta(tmp_path)

    def test_rejects_empty_header(self, tmp_path: Path) -> None:
        bad = tmp_path / "empty_header.fa"
        bad.write_text(">\nMKT\n")
        with pytest.raises(FastaFormatError, match="header is empty"):
            parse_fasta(bad)

    def test_rejects_file_with_no_records(self, tmp_path: Path) -> None:
        empty = tmp_path / "no_records.fa"
        empty.write_text("\n   \n\n")
        with pytest.raises(FastaFormatError, match="No FASTA records found"):
            parse_fasta(empty)

    def test_skips_blank_lines_within_record(self, tmp_path: Path) -> None:
        fa = tmp_path / "with_blanks.fa"
        fa.write_text(">seq1\nMKT\n\nIIA\n")
        records = parse_fasta(fa)
        assert records == [ProteinRecord("seq1", "seq1", "MKTIIA")]

    def test_strips_internal_whitespace_in_sequence(self, tmp_path: Path) -> None:
        fa = tmp_path / "spaces.fa"
        fa.write_text(">seq1\nM K T\nI I A\n")
        records = parse_fasta(fa)
        assert records[0].sequence == "MKTIIA"


class TestSummarizeRecords:
    def test_raises_on_empty_list(self) -> None:
        with pytest.raises(ValueError, match="No protein records"):
            summarize_records([])


class TestMainCli:
    def test_happy_path_prints_summary(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        fa = tmp_path / "ok.fa"
        fa.write_text(">seq1 desc\nMKT\n>seq2\nGAVL\n")

        exit_code = main([str(fa)])

        assert exit_code == 0
        out = capsys.readouterr().out
        assert "Sequence count: 2" in out
        assert "Min length: 3" in out
        assert "Max length: 4" in out
        assert "seq1\tMKT" in out
        assert "seq2\tGAVL" in out

    def test_returns_one_on_missing_file(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        exit_code = main([str(tmp_path / "nope.fa")])
        assert exit_code == 1
        err = capsys.readouterr().err
        assert "Error:" in err
        assert "FASTA file not found" in err

    def test_returns_one_on_malformed_input(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        bad = tmp_path / "bad.fa"
        bad.write_text("MKT\n")
        exit_code = main([str(bad)])
        assert exit_code == 1
        assert "Error:" in capsys.readouterr().err


class TestModuleEntrypoint:
    def test_runs_as_main(self, tmp_path: Path) -> None:
        """``python -m scripts.read_protein_fasta`` exercises the
        ``raise SystemExit(main())`` tail (line 86)."""
        import runpy
        import sys

        fa = tmp_path / "ok.fa"
        fa.write_text(">seq1\nMKT\n")
        original = sys.argv
        sys.argv = ["read_protein_fasta.py", str(fa)]
        try:
            with pytest.raises(SystemExit) as exc:
                runpy.run_module(
                    "scripts.read_protein_fasta", run_name="__main__"
                )
            assert exc.value.code == 0
        finally:
            sys.argv = original


if __name__ == "__main__":
    unittest.main()
