from __future__ import annotations

import subprocess
import unittest
from pathlib import Path

import pytest

from scripts.filter_gff3 import (
    Gff3FormatError,
    _normalize_multi_value_argument,
    _parse_attributes,
    filter_gff3_lines,
    main,
    parse_gff3,
)

pytestmark = pytest.mark.unit


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


class TestFilterGff3LinesEdgeCases:
    def test_no_filters_returns_all_genes(self, tmp_path: Path) -> None:
        fixture = Path(__file__).parent / "data" / "filter_gff3_input.gff3"
        out = filter_gff3_lines(fixture)
        # Every original feature line preserved (comments + features)
        assert any("ID=GeneA" in line for line in out)
        assert any("ID=GeneB" in line for line in out)
        assert any("ID=GeneC" in line for line in out)

    def test_rejects_negative_start(self, tmp_path: Path) -> None:
        fixture = Path(__file__).parent / "data" / "filter_gff3_input.gff3"
        with pytest.raises(ValueError, match="start must be"):
            filter_gff3_lines(fixture, start=0)

    def test_rejects_negative_end(self, tmp_path: Path) -> None:
        fixture = Path(__file__).parent / "data" / "filter_gff3_input.gff3"
        with pytest.raises(ValueError, match="end must be"):
            filter_gff3_lines(fixture, end=0)

    def test_rejects_start_after_end(self, tmp_path: Path) -> None:
        fixture = Path(__file__).parent / "data" / "filter_gff3_input.gff3"
        with pytest.raises(ValueError, match="start cannot be greater"):
            filter_gff3_lines(fixture, start=500, end=100)

    def test_empty_file_returns_only_comments(self, tmp_path: Path) -> None:
        empty = tmp_path / "only_comments.gff3"
        empty.write_text("##gff-version 3\n# nothing else\n")
        assert filter_gff3_lines(empty) == ["##gff-version 3", "# nothing else"]


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
        with pytest.raises(Gff3FormatError, match="expected 9 columns"):
            parse_gff3(bad)

    def test_non_integer_coords(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad_int.gff3"
        bad.write_text("Chr1\tsrc\tgene\tabc\txyz\t.\t+\t.\tID=G1\n")
        with pytest.raises(Gff3FormatError, match="must be integers"):
            parse_gff3(bad)

    def test_start_after_end(self, tmp_path: Path) -> None:
        bad = tmp_path / "swap.gff3"
        bad.write_text("Chr1\tsrc\tgene\t500\t100\t.\t+\t.\tID=G1\n")
        with pytest.raises(Gff3FormatError, match="cannot be greater than end"):
            parse_gff3(bad)

    def test_skips_blank_lines(self, tmp_path: Path) -> None:
        fa = tmp_path / "ok.gff3"
        fa.write_text("##gff-version 3\n\nChr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=G1\n\n")
        comments, features = parse_gff3(fa)
        assert comments == ["##gff-version 3"]
        assert len(features) == 1


class TestParseAttributes:
    def test_dot_returns_empty(self) -> None:
        assert _parse_attributes(".") == {}

    def test_skips_blank_items(self) -> None:
        assert _parse_attributes("ID=G1;;Name=Foo") == {"ID": "G1", "Name": "Foo"}

    def test_ignores_items_without_equals(self) -> None:
        assert _parse_attributes("ID=G1;orphan") == {"ID": "G1"}


class TestNormalizeMultiValueArgument:
    def test_none_returns_none(self) -> None:
        assert _normalize_multi_value_argument(None) is None

    def test_empty_returns_none(self) -> None:
        assert _normalize_multi_value_argument([]) is None

    def test_splits_comma_and_strips(self) -> None:
        assert _normalize_multi_value_argument(["a, b", "c"]) == {"a", "b", "c"}

    def test_only_whitespace_returns_none(self) -> None:
        assert _normalize_multi_value_argument(["  ,  "]) is None


class TestMainCli:
    def _fixture(self) -> Path:
        return Path(__file__).parent / "data" / "filter_gff3_input.gff3"

    def test_happy_path(self, capsys: pytest.CaptureFixture[str]) -> None:
        rc = main([str(self._fixture()), "--gene-id", "GeneC"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "ID=GeneC" in out
        assert "ID=GeneA" not in out

    def test_missing_file_returns_one(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        rc = main([str(tmp_path / "nope.gff3")])
        assert rc == 1
        assert "Error:" in capsys.readouterr().err

    def test_invalid_start_returns_one(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        rc = main([str(self._fixture()), "--start", "0"])
        assert rc == 1
        assert "start must be" in capsys.readouterr().err


class TestResolveRootGeneId:
    def test_orphan_feature_without_parent_chain(self, tmp_path: Path) -> None:
        # An mRNA whose Parent points to a non-existent gene must not be included.
        fa = tmp_path / "orphan.gff3"
        fa.write_text(
            "Chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=Real\n"
            "Chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=Orphan.1;Parent=Ghost\n"
        )
        out = filter_gff3_lines(fa, gene_ids={"Real"})
        assert any("ID=Real" in line for line in out)
        assert not any("Orphan.1" in line for line in out)


if __name__ == "__main__":
    unittest.main()
