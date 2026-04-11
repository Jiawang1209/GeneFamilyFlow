import subprocess
import unittest
from pathlib import Path

from scripts.filter_gff3 import filter_gff3_lines


class FilterGff3Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture_path = Path(__file__).parent / "data" / "filter_gff3_input.gff3"

    def test_filter_by_gene_id_keeps_gene_and_descendants(self) -> None:
        filtered = filter_gff3_lines(self.fixture_path, gene_ids={"GeneA"})

        self.assertEqual(
            filtered,
            [
                "##gff-version 3",
                "# test fixture for filter_gff3.py",
                "Chr1\tsource\tgene\t100\t300\t.\t+\t.\tID=GeneA;Name=GeneA",
                "Chr1\tsource\tmRNA\t100\t300\t.\t+\t.\tID=GeneA.1;Parent=GeneA",
                "Chr1\tsource\texon\t100\t180\t.\t+\t.\tID=GeneA.1.exon1;Parent=GeneA.1",
                "Chr1\tsource\tCDS\t150\t250\t.\t+\t0\tID=GeneA.1.cds;Parent=GeneA.1",
            ],
        )

    def test_filter_by_chromosome_and_range_uses_gene_coordinates(self) -> None:
        filtered = filter_gff3_lines(
            self.fixture_path,
            chromosome_ids={"Chr1"},
            start=350,
            end=950,
        )

        self.assertEqual(
            filtered,
            [
                "##gff-version 3",
                "# test fixture for filter_gff3.py",
                "Chr1\tsource\tgene\t400\t900\t.\t-\t.\tID=GeneB;Name=GeneB",
                "Chr1\tsource\tmRNA\t400\t900\t.\t-\t.\tID=GeneB.1;Parent=GeneB",
                "Chr1\tsource\texon\t400\t500\t.\t-\t.\tID=GeneB.1.exon1;Parent=GeneB.1",
            ],
        )

    def test_filter_with_combined_conditions_returns_intersection(self) -> None:
        filtered = filter_gff3_lines(
            self.fixture_path,
            gene_ids={"GeneA", "GeneB"},
            chromosome_ids={"Chr1"},
            start=50,
            end=350,
        )

        self.assertEqual(
            filtered,
            [
                "##gff-version 3",
                "# test fixture for filter_gff3.py",
                "Chr1\tsource\tgene\t100\t300\t.\t+\t.\tID=GeneA;Name=GeneA",
                "Chr1\tsource\tmRNA\t100\t300\t.\t+\t.\tID=GeneA.1;Parent=GeneA",
                "Chr1\tsource\texon\t100\t180\t.\t+\t.\tID=GeneA.1.exon1;Parent=GeneA.1",
                "Chr1\tsource\tCDS\t150\t250\t.\t+\t0\tID=GeneA.1.cds;Parent=GeneA.1",
            ],
        )

    def test_cli_writes_filtered_gff3_to_stdout(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                "scripts/filter_gff3.py",
                str(self.fixture_path),
                "--gene-id",
                "GeneC",
            ],
            cwd=Path(__file__).resolve().parent.parent,
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(completed.returncode, 0, msg=completed.stderr)
        self.assertEqual(
            completed.stdout,
            "\n".join(
                [
                    "##gff-version 3",
                    "# test fixture for filter_gff3.py",
                    "Chr2\tsource\tgene\t200\t500\t.\t+\t.\tID=GeneC;Name=GeneC",
                    "Chr2\tsource\tmRNA\t200\t500\t.\t+\t.\tID=GeneC.1;Parent=GeneC",
                    "Chr2\tsource\texon\t200\t260\t.\t+\t.\tID=GeneC.1.exon1;Parent=GeneC.1",
                    "",
                ]
            ),
        )


if __name__ == "__main__":
    unittest.main()
