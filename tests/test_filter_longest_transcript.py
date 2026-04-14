from __future__ import annotations

import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import pytest

from scripts.filter_longest_transcript import (
    FastaFormatError,
    collect_longest_transcripts,
    extract_gene_id,
    filter_longest_transcripts,
    main,
    parse_fasta,
)

pytestmark = pytest.mark.unit


class ExtractGeneIdTests(unittest.TestCase):
    def test_extract_gene_id_uses_first_matching_separator(self) -> None:
        gene_id, used_fallback = extract_gene_id("GeneA|Tx2", ["|", "."])

        self.assertEqual(gene_id, "GeneA")
        self.assertFalse(used_fallback)

    def test_extract_gene_id_falls_back_to_full_id_for_invalid_format(self) -> None:
        gene_id, used_fallback = extract_gene_id("MalformedId", ["|"])

        self.assertEqual(gene_id, "MalformedId")
        self.assertTrue(used_fallback)


class FilterLongestTranscriptTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture_dir = Path(__file__).parent / "data"

    def test_collect_longest_transcripts_keeps_longest_and_counts_fallbacks(self) -> None:
        selections, fallback_headers = collect_longest_transcripts(
            self.fixture_dir / "longest_transcript_input.fa",
            ["|"],
        )

        self.assertEqual(selections["GeneA"].record_id, "GeneA|Tx2")
        self.assertEqual(selections["GeneA"].sequence_length, 6)
        self.assertEqual(selections["GeneB"].record_id, "GeneB|Tx1")
        self.assertEqual(selections["MalformedId"].record_id, "MalformedId")
        self.assertEqual(fallback_headers, ["MalformedId"])

    def test_filter_longest_transcripts_writes_only_selected_records(self) -> None:
        input_path = self.fixture_dir / "longest_transcript_input.fa"
        output_path = Path("ignored_output.fa")
        output_buffer = NonClosingStringIO()
        original_open = Path.open

        def patched_open(path_obj: Path, mode: str = "r", *args, **kwargs):
            if path_obj == output_path and mode == "w":
                return output_buffer
            return original_open(path_obj, mode, *args, **kwargs)

        with patch.object(Path, "open", new=patched_open):
            summary = filter_longest_transcripts(input_path, output_path, ["|"], line_width=4)

        self.assertEqual(summary["input_records"], 5)
        self.assertEqual(summary["output_records"], 3)
        self.assertEqual(summary["fallback_records"], 1)
        self.assertEqual(
            output_buffer.getvalue(),
            "\n".join(
                [
                    ">GeneA|Tx2",
                    "ATGC",
                    "AA",
                    ">GeneB|Tx1",
                    "AT",
                    ">MalformedId",
                    "ATGC",
                    "GA",
                    "",
                ]
            ),
        )

    def test_filter_longest_transcripts_rejects_invalid_fasta(self) -> None:
        with self.assertRaises(FastaFormatError):
            filter_longest_transcripts(
                self.fixture_dir / "bad_before_header.fa",
                Path("ignored_output.fa"),
                ["|"],
            )


class NonClosingStringIO(StringIO):
    def close(self) -> None:
        self.seek(0, 2)


class TestExtractGeneIdEdgeCases:
    def test_skips_empty_separators(self) -> None:
        gene_id, used_fallback = extract_gene_id("GeneA|Tx1", ["", "|"])
        assert gene_id == "GeneA"
        assert used_fallback is False

    def test_falls_back_when_separator_yields_empty_prefix(self) -> None:
        # Leading separator => prefix is empty => break, fallback to full id
        gene_id, used_fallback = extract_gene_id("|Tx1", ["|"])
        assert gene_id == "|Tx1"
        assert used_fallback is True


class TestParseFastaErrorBranches:
    def test_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not found"):
            list(parse_fasta(tmp_path / "no.fa"))

    def test_directory_path(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not a file"):
            list(parse_fasta(tmp_path))

    def test_empty_header(self, tmp_path: Path) -> None:
        bad = tmp_path / "empty_header.fa"
        bad.write_text(">\nATGC\n")
        with pytest.raises(FastaFormatError, match="header is empty"):
            list(parse_fasta(bad))

    def test_no_records(self, tmp_path: Path) -> None:
        bad = tmp_path / "blank.fa"
        bad.write_text("\n\n")
        with pytest.raises(FastaFormatError, match="No FASTA records"):
            list(parse_fasta(bad))

    def test_skips_blank_lines_within_record(self, tmp_path: Path) -> None:
        fa = tmp_path / "ok.fa"
        fa.write_text(">seq1\nATGC\n\nGGAA\n")
        records = list(parse_fasta(fa))
        assert records[0].sequence == "ATGCGGAA"


class TestFilterLongestTranscriptsValidation:
    def test_negative_line_width(self, tmp_path: Path) -> None:
        fa = tmp_path / "ok.fa"
        fa.write_text(">GeneA|Tx1\nATGC\n")
        out = tmp_path / "out.fa"
        with pytest.raises(ValueError, match="line_width must be"):
            filter_longest_transcripts(fa, out, ["|"], line_width=-1)

    def test_line_width_zero_writes_single_line(self, tmp_path: Path) -> None:
        fa = tmp_path / "ok.fa"
        fa.write_text(">GeneA|Tx1\nATGCATGC\n")
        out = tmp_path / "out.fa"
        filter_longest_transcripts(fa, out, ["|"], line_width=0)
        contents = out.read_text().splitlines()
        # header + one sequence line + final empty from trailing \n
        assert contents == [">GeneA|Tx1", "ATGCATGC"]


class TestMainCli:
    def test_happy_path(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        fa = tmp_path / "in.fa"
        fa.write_text(">GeneA|Tx1\nATG\n>GeneA|Tx2\nATGCAA\n")
        out = tmp_path / "out.fa"
        rc = main([str(fa), str(out)])
        assert rc == 0
        captured = capsys.readouterr()
        assert "Input records: 2" in captured.out
        assert "Output records: 1" in captured.out
        assert "Fallback IDs: 0" in captured.err
        assert out.read_text().startswith(">GeneA|Tx2\n")

    def test_warns_on_fallback_records(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        fa = tmp_path / "in.fa"
        fa.write_text(">GeneA|Tx1\nATG\n>NoSep\nATGC\n")
        out = tmp_path / "out.fa"
        rc = main([str(fa), str(out)])
        assert rc == 0
        err = capsys.readouterr().err
        assert "Fallback IDs: 1" in err
        assert "Warning" in err

    def test_missing_file_returns_one(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        rc = main([str(tmp_path / "nope.fa"), str(tmp_path / "out.fa")])
        assert rc == 1
        assert "Error:" in capsys.readouterr().err

    def test_custom_separator_argument(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        fa = tmp_path / "in.fa"
        fa.write_text(">GeneA.Tx1\nATG\n>GeneA.Tx2\nATGCC\n")
        out = tmp_path / "out.fa"
        rc = main([str(fa), str(out), "-s", "."])
        assert rc == 0
        assert out.read_text().startswith(">GeneA.Tx2\n")


if __name__ == "__main__":
    unittest.main()
