"""Tests for scripts/merge_gene_ids.py."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.merge_gene_ids import main, merge_ids, read_ids


class TestMergeGeneIds(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)

        self.file_a = self.tmp / "a.txt"
        self.file_a.write_text("GeneA\nGeneB\nGeneC\nGeneD\n")

        self.file_b = self.tmp / "b.txt"
        self.file_b.write_text("GeneB\nGeneC\nGeneE\nGeneF\n")

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_read_ids(self) -> None:
        ids = read_ids(self.file_a)
        self.assertEqual(ids, {"GeneA", "GeneB", "GeneC", "GeneD"})

    def test_read_ids_skips_blanks_and_comments(self) -> None:
        f = self.tmp / "with_comments.txt"
        f.write_text("# header\nGeneX\n\n  \nGeneY\n# another comment\n")
        ids = read_ids(f)
        self.assertEqual(ids, {"GeneX", "GeneY"})

    def test_read_ids_file_not_found(self) -> None:
        with self.assertRaises(FileNotFoundError):
            read_ids("/nonexistent/file.txt")

    def test_merge_intersection(self) -> None:
        a = {"GeneA", "GeneB", "GeneC"}
        b = {"GeneB", "GeneC", "GeneD"}
        result = merge_ids(a, b, "intersection")
        self.assertEqual(result, ["GeneB", "GeneC"])

    def test_merge_union(self) -> None:
        a = {"GeneA", "GeneB"}
        b = {"GeneB", "GeneC"}
        result = merge_ids(a, b, "union")
        self.assertEqual(result, ["GeneA", "GeneB", "GeneC"])

    def test_merge_empty_sets(self) -> None:
        self.assertEqual(merge_ids(set(), set(), "intersection"), [])
        self.assertEqual(merge_ids(set(), set(), "union"), [])

    def test_merge_unknown_method(self) -> None:
        with self.assertRaises(ValueError):
            merge_ids({"A"}, {"B"}, "difference")

    def test_merge_result_sorted(self) -> None:
        result = merge_ids({"Z", "A", "M"}, {"Z", "A", "M"}, "union")
        self.assertEqual(result, sorted(result))

    def test_cli_intersection(self) -> None:
        outfile = self.tmp / "merged.txt"
        rc = main([str(self.file_a), str(self.file_b), "-m", "intersection", "-o", str(outfile)])
        self.assertEqual(rc, 0)
        lines = outfile.read_text().strip().split("\n")
        self.assertEqual(sorted(lines), ["GeneB", "GeneC"])

    def test_cli_union(self) -> None:
        outfile = self.tmp / "merged.txt"
        rc = main([str(self.file_a), str(self.file_b), "-m", "union", "-o", str(outfile)])
        self.assertEqual(rc, 0)
        lines = outfile.read_text().strip().split("\n")
        self.assertEqual(len(lines), 6)

    def test_cli_missing_file(self) -> None:
        rc = main(["/nonexistent.txt", str(self.file_b)])
        self.assertEqual(rc, 1)


if __name__ == "__main__":
    unittest.main()
