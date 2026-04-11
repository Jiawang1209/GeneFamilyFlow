import unittest
from pathlib import Path

from scripts.read_protein_fasta import (
    FastaFormatError,
    ProteinRecord,
    parse_fasta,
    summarize_records,
)


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


if __name__ == "__main__":
    unittest.main()
