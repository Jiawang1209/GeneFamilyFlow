from __future__ import annotations

import csv
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import pytest

from scripts.extract_longest_transcript_from_gff3 import (
    Gff3Feature,
    Gff3TranscriptSelectionError,
    _parse_attributes,
    _resolve_gene_parent,
    collect_longest_transcripts,
    main,
    parse_gff3,
    write_longest_transcripts_csv,
)

pytestmark = pytest.mark.unit


class ExtractLongestTranscriptFromGff3Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture_path = Path(__file__).parent / "data" / "longest_transcript_from_gff3_input.gff3"

    def test_collect_longest_transcripts_uses_total_cds_length_per_gene(self) -> None:
        selections = collect_longest_transcripts(self.fixture_path)

        self.assertEqual(selections["GeneA"].transcript_id, "GeneA.2")
        self.assertEqual(selections["GeneA"].cds_length, 452)
        self.assertEqual(selections["GeneA"].seqid, "Chr1")
        self.assertEqual(selections["GeneA"].start, 100)
        self.assertEqual(selections["GeneA"].end, 850)
        self.assertEqual(selections["GeneA"].strand, "+")

        self.assertEqual(selections["GeneB"].transcript_id, "GeneB.t2")
        self.assertEqual(selections["GeneB"].cds_length, 332)
        self.assertEqual(selections["GeneB"].strand, "-")

    def test_collect_longest_transcripts_keeps_first_transcript_when_cds_lengths_tie(self) -> None:
        selections = collect_longest_transcripts(self.fixture_path)

        self.assertEqual(selections["GeneC"].transcript_id, "GeneC.1")
        self.assertEqual(selections["GeneC"].cds_length, 100)

    def test_collect_longest_transcripts_rejects_transcripts_without_cds(self) -> None:
        with self.assertRaises(Gff3TranscriptSelectionError):
            collect_longest_transcripts(self.fixture_path, include_genes_without_cds=False)

    def test_write_longest_transcripts_csv_outputs_expected_columns(self) -> None:
        output_path = Path("ignored_longest.csv")
        output_buffer = NonClosingStringIO()
        original_open = Path.open

        def patched_open(path_obj: Path, mode: str = "r", *args, **kwargs):
            if path_obj == output_path and mode == "w":
                return output_buffer
            return original_open(path_obj, mode, *args, **kwargs)

        with patch.object(Path, "open", new=patched_open):
            summary = write_longest_transcripts_csv(self.fixture_path, output_path)
            rows = list(csv.DictReader(StringIO(output_buffer.getvalue())))

        self.assertEqual(
            rows,
            [
                {
                    "gene_id": "GeneA",
                    "transcript_id": "GeneA.2",
                    "seqid": "Chr1",
                    "start": "100",
                    "end": "850",
                    "strand": "+",
                    "cds_length": "452",
                },
                {
                    "gene_id": "GeneB",
                    "transcript_id": "GeneB.t2",
                    "seqid": "Chr2",
                    "start": "1050",
                    "end": "1580",
                    "strand": "-",
                    "cds_length": "332",
                },
                {
                    "gene_id": "GeneC",
                    "transcript_id": "GeneC.1",
                    "seqid": "Chr3",
                    "start": "200",
                    "end": "700",
                    "strand": "+",
                    "cds_length": "100",
                },
            ],
        )

    def test_cli_writes_csv_and_reports_skipped_gene_count(self) -> None:
        output_path = Path("ignored_cli_longest.csv")
        output_buffer = NonClosingStringIO()
        stdout_buffer = StringIO()
        stderr_buffer = StringIO()
        original_open = Path.open

        def patched_open(path_obj: Path, mode: str = "r", *args, **kwargs):
            if path_obj == output_path and mode == "w":
                return output_buffer
            return original_open(path_obj, mode, *args, **kwargs)

        with (
            patch.object(Path, "open", new=patched_open),
            patch("sys.stdout", new=stdout_buffer),
            patch("sys.stderr", new=stderr_buffer),
        ):
            return_code = main([str(self.fixture_path), str(output_path)])

        self.assertEqual(return_code, 0, msg=stderr_buffer.getvalue())
        self.assertIn("Genes written: 3", stdout_buffer.getvalue())
        self.assertIn(
            "Skipped genes without CDS-backed transcripts: 1",
            stderr_buffer.getvalue(),
        )
        rows = list(csv.DictReader(StringIO(output_buffer.getvalue())))

        self.assertEqual(len(rows), 3)


class NonClosingStringIO(StringIO):
    def close(self) -> None:
        self.seek(0, 2)


class TestGff3FeatureParents:
    def test_parents_empty_when_attribute_missing(self) -> None:
        feature = Gff3Feature(
            seqid="Chr1",
            feature_type="gene",
            start=1,
            end=10,
            strand="+",
            attributes={"ID": "G1"},
        )
        assert feature.parents == []


class TestParseGff3Errors:
    def test_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not found"):
            parse_gff3(tmp_path / "nope.gff3")

    def test_directory_path(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not a file"):
            parse_gff3(tmp_path)

    def test_wrong_column_count(self, tmp_path: Path) -> None:
        bad = tmp_path / "short.gff3"
        bad.write_text("Chr1\tsrc\tgene\t1\t10\n")
        with pytest.raises(Gff3TranscriptSelectionError, match="expected 9 columns"):
            parse_gff3(bad)

    def test_non_integer_coords(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad_int.gff3"
        bad.write_text("Chr1\tsrc\tgene\tabc\txyz\t.\t+\t.\tID=G1\n")
        with pytest.raises(Gff3TranscriptSelectionError, match="must be integers"):
            parse_gff3(bad)

    def test_start_after_end(self, tmp_path: Path) -> None:
        bad = tmp_path / "swap.gff3"
        bad.write_text("Chr1\tsrc\tgene\t500\t100\t.\t+\t.\tID=G1\n")
        with pytest.raises(Gff3TranscriptSelectionError, match="cannot be greater"):
            parse_gff3(bad)

    def test_no_features(self, tmp_path: Path) -> None:
        bad = tmp_path / "empty.gff3"
        bad.write_text("##gff-version 3\n# only comments\n")
        with pytest.raises(Gff3TranscriptSelectionError, match="No GFF3 features"):
            parse_gff3(bad)


class TestParseAttributes:
    def test_dot_returns_empty(self) -> None:
        assert _parse_attributes(".") == {}

    def test_skips_blank_items(self) -> None:
        assert _parse_attributes("ID=G1;;Name=Foo") == {"ID": "G1", "Name": "Foo"}

    def test_rejects_item_without_equals(self) -> None:
        with pytest.raises(Gff3TranscriptSelectionError, match="Malformed"):
            _parse_attributes("ID=G1;orphan")

    def test_rejects_empty_key(self) -> None:
        with pytest.raises(Gff3TranscriptSelectionError, match="Malformed"):
            _parse_attributes("=value")


class TestResolveGeneParent:
    def test_returns_none_when_no_parent_matches(self) -> None:
        feature = Gff3Feature(
            seqid="Chr1",
            feature_type="mRNA",
            start=1,
            end=10,
            strand="+",
            attributes={"ID": "T1", "Parent": "Ghost"},
        )
        assert _resolve_gene_parent(feature, {}) is None


class TestCollectLongestTranscriptsEdgeCases:
    def test_skips_cds_with_unknown_transcript_parent(self, tmp_path: Path) -> None:
        # CDS feature with no Parent attribute (parents == []) should be ignored
        # by the `for transcript_id in feature.parents:` loop, leaving the
        # gene with zero CDS-backed transcripts.
        gff = tmp_path / "orphan_cds.gff3"
        gff.write_text(
            "Chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=G1\n"
            "Chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=G1.1;Parent=G1\n"
            "Chr1\tsrc\tCDS\t10\t50\t.\t+\t0\tID=G1.1.cds\n"
        )
        selections = collect_longest_transcripts(gff)
        assert "G1" not in selections


class TestMainCliErrors:
    def test_returns_one_on_missing_file(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        rc = main([str(tmp_path / "nope.gff3"), str(tmp_path / "out.csv")])
        assert rc == 1
        assert "Error:" in capsys.readouterr().err

    def test_returns_one_on_strict_with_missing_cds(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        gff = tmp_path / "no_cds.gff3"
        gff.write_text(
            "Chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=G1\n"
            "Chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=G1.1;Parent=G1\n"
        )
        rc = main([str(gff), str(tmp_path / "out.csv"), "--strict"])
        assert rc == 1
        assert "without CDS" in capsys.readouterr().err


if __name__ == "__main__":
    unittest.main()
