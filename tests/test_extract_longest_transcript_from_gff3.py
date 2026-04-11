import csv
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from scripts.extract_longest_transcript_from_gff3 import (
    Gff3TranscriptSelectionError,
    collect_longest_transcripts,
    main,
    write_longest_transcripts_csv,
)


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


if __name__ == "__main__":
    unittest.main()
