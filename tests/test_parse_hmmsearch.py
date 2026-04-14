"""Tests for scripts/parse_hmmsearch.py."""

from __future__ import annotations

import textwrap
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.parse_hmmsearch import (
    DomtbloutHit,
    extract_ids,
    main,
    parse_domtblout,
)

SAMPLE_DOMTBLOUT = textwrap.dedent("""\
    #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
    # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
    #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
    GeneA.001G000100.1.p PF00854.25   457 custom_hmm           -            392   1.5e-30  104.0   8.8   1   1   2.0e-30   3.0e-30  103.5   4.3    62   410    62   457 PF00854 PTR2 family
    GeneB.002G000200.1.p PF00854.25   507 custom_hmm           -            392   5.8e-120 402.8  11.2   1   2   3.6e-62   7.2e-62  210.2   3.7   114   499   114   507 PF00854 PTR2 family
    GeneC.003G000300.1.p PF00854.25   521 custom_hmm           -            392   2.0e-03   12.5   1.1   1   1   3.0e-03   6.0e-03   11.8   0.5    96   521    96   521 PF00854 PTR2 family
    GeneD.004G000400.1.p PF00854.25   518 custom_hmm           -            392   7.1e-10   38.1   0.2   1   1   8.5e-10   1.7e-09   37.2   0.1    96   518    96   518 PF00854 PTR2 family
""")


class TestParseHmmsearch(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.domtblout = self.tmp / "test.domtblout"
        self.domtblout.write_text(SAMPLE_DOMTBLOUT)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_parse_domtblout_yields_all_hits(self) -> None:
        hits = list(parse_domtblout(self.domtblout))
        self.assertEqual(len(hits), 4)
        self.assertEqual(hits[0].target_name, "GeneA.001G000100.1.p")
        self.assertAlmostEqual(hits[0].full_evalue, 1.5e-30)

    def test_parse_domtblout_skips_comments(self) -> None:
        hits = list(parse_domtblout(self.domtblout))
        names = [h.target_name for h in hits]
        self.assertFalse(any(n.startswith("#") for n in names))

    def test_extract_ids_with_evalue_filter(self) -> None:
        ids = extract_ids([self.domtblout], evalue=1e-5)
        self.assertIn("GeneA.001G000100.1.p", ids)
        self.assertIn("GeneB.002G000200.1.p", ids)
        self.assertIn("GeneD.004G000400.1.p", ids)
        self.assertNotIn("GeneC.003G000300.1.p", ids)

    def test_extract_ids_strict_evalue(self) -> None:
        ids = extract_ids([self.domtblout], evalue=1e-100)
        self.assertEqual(ids, ["GeneB.002G000200.1.p"])

    def test_extract_ids_domain_evalue(self) -> None:
        ids = extract_ids([self.domtblout], evalue=1e-10, use_domain_evalue=True)
        self.assertIn("GeneA.001G000100.1.p", ids)
        self.assertIn("GeneB.002G000200.1.p", ids)

    def test_extract_ids_multiple_files(self) -> None:
        second = self.tmp / "second.domtblout"
        second.write_text(
            "GeneE.005G000500.1.p PF00854.25 500 custom_hmm - 392 "
            "1.0e-50 180.0 2.0 1 1 1.5e-50 3.0e-50 179.5 1.0 "
            "10 400 10 400 10 400 0.95 PTR2 family\n"
        )
        ids = extract_ids([self.domtblout, second], evalue=1e-5)
        self.assertIn("GeneE.005G000500.1.p", ids)
        self.assertTrue(len(ids) >= 4)

    def test_extract_ids_sorted(self) -> None:
        ids = extract_ids([self.domtblout], evalue=1.0)
        self.assertEqual(ids, sorted(ids))

    def test_empty_file(self) -> None:
        empty = self.tmp / "empty.domtblout"
        empty.write_text("# just comments\n")
        ids = extract_ids([empty], evalue=1.0)
        self.assertEqual(ids, [])

    def test_file_not_found(self) -> None:
        with self.assertRaises(FileNotFoundError):
            list(parse_domtblout("/nonexistent/file.domtblout"))

    def test_cli_main(self) -> None:
        outfile = self.tmp / "ids.txt"
        rc = main([str(self.domtblout), "-e", "1e-5", "-o", str(outfile)])
        self.assertEqual(rc, 0)
        self.assertTrue(outfile.exists())
        lines = outfile.read_text().strip().split("\n")
        self.assertEqual(len(lines), 3)

    def test_cli_missing_file(self) -> None:
        rc = main(["/nonexistent.domtblout", "-e", "1e-5"])
        self.assertEqual(rc, 1)

    def test_parse_domtblout_skips_short_lines(self) -> None:
        """Lines with fewer than 22 fields should be silently dropped."""
        short = self.tmp / "short.domtblout"
        short.write_text(
            "# short row below\n"
            "only_three fields here\n"
            "GeneX.001G000000.1.p PF00854.25 457 custom_hmm - 392 "
            "1.5e-30 104.0 8.8 1 1 2.0e-30 3.0e-30 103.5 4.3 "
            "62 410 62 457 0.99 PTR2 family\n"
        )
        hits = list(parse_domtblout(short))
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].target_name, "GeneX.001G000000.1.p")

    def test_parse_domtblout_skips_non_numeric_rows(self) -> None:
        """Rows where E-value/score aren't numeric should be skipped."""
        bad = self.tmp / "bad.domtblout"
        bad.write_text(
            "GeneY.002G000000.1.p PF00854.25 457 custom_hmm - 392 "
            "NOT_A_FLOAT 104.0 8.8 1 1 2.0e-30 3.0e-30 103.5 4.3 "
            "62 410 62 457 0.99 PTR2 family\n"
        )
        self.assertEqual(list(parse_domtblout(bad)), [])

    def test_cli_main_stdout(self) -> None:
        """Without -o, main() prints IDs to stdout."""
        import io
        from contextlib import redirect_stdout

        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main([str(self.domtblout), "-e", "1e-5"])
        self.assertEqual(rc, 0)
        lines = [ln for ln in buf.getvalue().splitlines() if ln]
        self.assertEqual(len(lines), 3)
        self.assertIn("GeneA.001G000100.1.p", lines)


if __name__ == "__main__":
    unittest.main()
