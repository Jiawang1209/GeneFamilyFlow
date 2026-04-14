"""Tests for scripts/parse_pfam_scan.py."""

from __future__ import annotations

import textwrap
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.parse_pfam_scan import (
    PfamScanHit,
    extract_ids_by_pfam,
    main,
    parse_pfam_scan,
)

SAMPLE_PFAM_SCAN = textwrap.dedent("""\
    # pfam_scan.pl,  run at Thu Jul 10 10:20:19 2025
    #
    # Copyright (c) 2009 Genome Research Ltd
    # Freely distributed under the GNU
    # General Public License
    AT2G38100.1              62    410     62    457 PF00854.25  PTR2              Family     1   332   392    104.0   8.8e-30   1 CL0015
    AT2G02020.1             114    499    114    507 PF00854.25  PTR2              Family     1   383   392    402.8  1.6e-120   1 CL0015
    AT2G26690.1              96    521     96    521 PF00854.25  PTR2              Family     1   392   392    346.5    2e-103   1 CL0015
    AT5G01180.1              50    300     50    305 PF01234.10  OTHER             Family     1   200   250     80.0   5.0e-20   1 CL0020
    AT3G99999.1             100    400    100    400 PF00854.25  PTR2              Family     1   300   392     50.0   1.0e-03   1 CL0015
""")


class TestParsePfamScan(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.pfam_file = self.tmp / "Pfam_scan.out"
        self.pfam_file.write_text(SAMPLE_PFAM_SCAN)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_parse_yields_all_hits(self) -> None:
        hits = list(parse_pfam_scan(self.pfam_file))
        self.assertEqual(len(hits), 5)

    def test_parse_skips_comments(self) -> None:
        hits = list(parse_pfam_scan(self.pfam_file))
        self.assertFalse(any(h.seq_id.startswith("#") for h in hits))

    def test_parse_fields(self) -> None:
        hits = list(parse_pfam_scan(self.pfam_file))
        first = hits[0]
        self.assertEqual(first.seq_id, "AT2G38100.1")
        self.assertEqual(first.ali_start, 62)
        self.assertEqual(first.ali_end, 410)
        self.assertEqual(first.pfam_acc, "PF00854.25")
        self.assertEqual(first.pfam_name, "PTR2")
        self.assertAlmostEqual(first.evalue, 8.8e-30)
        self.assertEqual(first.clan, "CL0015")

    def test_extract_by_pfam_id(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854")
        self.assertEqual(len(ids), 4)
        self.assertNotIn("AT5G01180.1", ids)

    def test_extract_by_pfam_id_with_version(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854")
        self.assertIn("AT2G38100.1", ids)

    def test_extract_with_evalue_filter(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854", evalue=1e-10)
        self.assertEqual(len(ids), 3)
        self.assertNotIn("AT3G99999.1", ids)

    def test_extract_different_pfam(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF01234")
        self.assertEqual(ids, ["AT5G01180.1"])

    def test_extract_nonexistent_pfam(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF99999")
        self.assertEqual(ids, [])

    def test_extract_sorted(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854")
        self.assertEqual(ids, sorted(ids))

    def test_empty_file(self) -> None:
        empty = self.tmp / "empty.out"
        empty.write_text("# just comments\n")
        ids = extract_ids_by_pfam(empty, "PF00854")
        self.assertEqual(ids, [])

    def test_file_not_found(self) -> None:
        with self.assertRaises(FileNotFoundError):
            list(parse_pfam_scan("/nonexistent/file.out"))

    def test_cli_main(self) -> None:
        outfile = self.tmp / "ids.txt"
        rc = main([str(self.pfam_file), "--pfam-id", "PF00854", "-o", str(outfile)])
        self.assertEqual(rc, 0)
        lines = outfile.read_text().strip().split("\n")
        self.assertEqual(len(lines), 4)

    def test_cli_with_evalue(self) -> None:
        outfile = self.tmp / "ids.txt"
        rc = main([str(self.pfam_file), "--pfam-id", "PF00854", "-e", "1e-10", "-o", str(outfile)])
        self.assertEqual(rc, 0)
        lines = outfile.read_text().strip().split("\n")
        self.assertEqual(len(lines), 3)

    def test_extract_multiple_pfam_ids(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854,PF01234")
        self.assertEqual(len(ids), 5)
        self.assertIn("AT5G01180.1", ids)
        self.assertIn("AT2G38100.1", ids)

    def test_extract_multiple_pfam_ids_with_evalue(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854,PF01234", evalue=1e-10)
        self.assertEqual(len(ids), 4)
        self.assertNotIn("AT3G99999.1", ids)

    def test_cli_missing_file(self) -> None:
        rc = main(["/nonexistent.out", "--pfam-id", "PF00854"])
        self.assertEqual(rc, 1)


MULTI_DOMAIN_SAMPLE = textwrap.dedent("""\
    # pfam_scan.pl multi-domain test fixture
    G1_single_A         10    100    10   110 PF00854.25  DomA Family 1 100 100  80.0 1e-20 1 CL0015
    G2_single_B         10    100    10   110 PF01234.10  DomB Family 1 100 100  80.0 1e-20 1 CL0020
    G3_both             10    100    10   110 PF00854.25  DomA Family 1 100 100  80.0 1e-20 1 CL0015
    G3_both            150    250   150   260 PF01234.10  DomB Family 1 100 100  80.0 1e-20 1 CL0020
    G4_both             10    100    10   110 PF00854.25  DomA Family 1 100 100  80.0 1e-20 1 CL0015
    G4_both            150    250   150   260 PF01234.10  DomB Family 1 100 100  80.0 1e-20 1 CL0020
    G5_unrelated        10    100    10   110 PF99999.1   DomX Family 1 100 100  80.0 1e-20 1 CL0099
""")


class TestMultiDomainMode(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.pfam_file = self.tmp / "multi_domain.out"
        self.pfam_file.write_text(MULTI_DOMAIN_SAMPLE)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_any_mode_returns_union(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854,PF01234", mode="any")
        self.assertEqual(set(ids), {"G1_single_A", "G2_single_B", "G3_both", "G4_both"})

    def test_all_mode_returns_intersection(self) -> None:
        ids = extract_ids_by_pfam(self.pfam_file, "PF00854,PF01234", mode="all")
        self.assertEqual(set(ids), {"G3_both", "G4_both"})

    def test_all_mode_single_domain_matches_any(self) -> None:
        """When only one domain listed, all and any should produce same result."""
        any_ids = extract_ids_by_pfam(self.pfam_file, "PF00854", mode="any")
        all_ids = extract_ids_by_pfam(self.pfam_file, "PF00854", mode="all")
        self.assertEqual(any_ids, all_ids)

    def test_default_mode_is_any(self) -> None:
        """Backward compatibility: default mode must be 'any'."""
        default_ids = extract_ids_by_pfam(self.pfam_file, "PF00854,PF01234")
        any_ids = extract_ids_by_pfam(self.pfam_file, "PF00854,PF01234", mode="any")
        self.assertEqual(default_ids, any_ids)

    def test_invalid_mode_raises(self) -> None:
        with self.assertRaises(ValueError):
            extract_ids_by_pfam(self.pfam_file, "PF00854", mode="invalid")

    def test_all_mode_with_evalue_filter(self) -> None:
        """E-value filter must apply before the ALL aggregation."""
        ids = extract_ids_by_pfam(
            self.pfam_file, "PF00854,PF01234", mode="all", evalue=1e-25
        )
        self.assertEqual(ids, [])

    def test_cli_mode_all(self) -> None:
        outfile = self.tmp / "all.txt"
        rc = main([
            str(self.pfam_file),
            "--pfam-id", "PF00854,PF01234",
            "--mode", "all",
            "-o", str(outfile),
        ])
        self.assertEqual(rc, 0)
        lines = outfile.read_text().strip().split("\n")
        self.assertEqual(set(lines), {"G3_both", "G4_both"})

    def test_cli_mode_any_explicit(self) -> None:
        outfile = self.tmp / "any.txt"
        rc = main([
            str(self.pfam_file),
            "--pfam-id", "PF00854,PF01234",
            "--mode", "any",
            "-o", str(outfile),
        ])
        self.assertEqual(rc, 0)
        lines = outfile.read_text().strip().split("\n")
        self.assertEqual(
            set(lines), {"G1_single_A", "G2_single_B", "G3_both", "G4_both"}
        )


import pytest


class TestParsePfamScanDefensiveBranches:
    def test_skips_short_lines(self, tmp_path: Path) -> None:
        bad = tmp_path / "short.txt"
        # 8 fields only — below the 15-field minimum
        bad.write_text("seq1 1 50 1 50 PF00854.1 Name Domain\n")
        assert list(parse_pfam_scan(bad)) == []

    def test_skips_unparseable_numeric_fields(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad_nums.txt"
        # 15 fields, but ali_start is non-integer => ValueError swallowed
        bad.write_text(
            "seq1 NOT_AN_INT 50 1 50 PF00854.1 Name Domain x y z 100.0 1e-30 . CL0001\n"
        )
        assert list(parse_pfam_scan(bad)) == []


class TestMainCliExtraBranches:
    @pytest.fixture
    def fixture_path(self, tmp_path: Path) -> Path:
        path = tmp_path / "input.txt"
        path.write_text(
            "# header\n"
            "G1 1 50 1 50 PF00854.1 Name Domain x y z 100.0 1e-30 . CL0001\n"
        )
        return path

    def test_stdout_output_when_no_dash_o(
        self, fixture_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        rc = main([str(fixture_path), "--pfam-id", "PF00854"])
        assert rc == 0
        captured = capsys.readouterr()
        assert "G1" in captured.out
        assert "Extracted 1" in captured.err

    def test_extract_error_returns_one(
        self,
        fixture_path: Path,
        capsys: pytest.CaptureFixture[str],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        from scripts import parse_pfam_scan as mod

        def boom(*_a, **_kw):
            raise RuntimeError("boom")

        monkeypatch.setattr(mod, "extract_ids_by_pfam", boom)
        rc = main([str(fixture_path), "--pfam-id", "PF00854"])
        assert rc == 1
        assert "Error: boom" in capsys.readouterr().err


if __name__ == "__main__":
    unittest.main()
