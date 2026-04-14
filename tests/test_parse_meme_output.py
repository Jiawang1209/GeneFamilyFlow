"""Tests for scripts/parse_meme_output.py."""

from __future__ import annotations

import textwrap
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.parse_meme_output import (
    MotifInfo,
    MotifLocation,
    load_fasta_ids,
    main,
    parse_meme,
    write_outputs,
)

SAMPLE_MEME = textwrap.dedent("""\
    ********************************************************************************
    MEME - Motif discovery tool
    ********************************************************************************

    \tMotif ABCDEFGHIJ MEME-1 Description
    ********************************************************************************

    Motif ABCDEFGHIJ MEME-1 sites sorted by position p-value
    ********************************************************************************
    Sequence name             Start   P-value                            Site
    -------------             ----- ---------            ------------------------------------
    GeneA.001G000100.1.p        123  1.00e-20 XXXXXXXXXX ABCDEFGHIJ YYYYYYYYYY
    GeneB.002G000200.2.p        456  2.00e-15 XXXXXXXXXX ABCDEFGHIJ YYYYYYYYYY
    GeneC.003G000300.1.p        789  3.00e-10 XXXXXXXXXX ABCDEFGHIJ YYYYYYYYYY

    \tMotif KLMNO MEME-2 Description
    ********************************************************************************

    Motif KLMNO MEME-2 sites sorted by position p-value
    ********************************************************************************
    Sequence name             Start   P-value           Site
    -------------             ----- ---------           -----
    GeneA.001G000100.1.p         50  1.00e-08 XXXXXXXXXX KLMNO YYYYYYYYYY
    GeneC.003G000300.1.p        200  2.00e-06 XXXXXXXXXX KLMNO YYYYYYYYYY
""")


class TestParseMeme(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.meme_file = self.tmp / "meme.txt"
        self.meme_file.write_text(SAMPLE_MEME)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_parse_all_genes(self) -> None:
        motifs, locations = parse_meme(self.meme_file)
        self.assertEqual(len(motifs), 2)
        self.assertEqual(motifs[0], MotifInfo("MEME-1", "ABCDEFGHIJ", 10))
        self.assertEqual(motifs[1], MotifInfo("MEME-2", "KLMNO", 5))
        self.assertEqual(len(locations), 5)

    def test_parse_with_gene_filter(self) -> None:
        motifs, locations = parse_meme(self.meme_file, gene_id_pattern=r"GeneA\.")
        self.assertEqual(len(motifs), 2)
        self.assertEqual(len(locations), 2)
        self.assertTrue(all(loc.gene.startswith("GeneA") for loc in locations))

    def test_location_coordinates(self) -> None:
        _, locations = parse_meme(self.meme_file)
        first = locations[0]
        self.assertEqual(first.gene, "GeneA.001G000100.1.p")
        self.assertEqual(first.start, 123)
        self.assertEqual(first.end, 132)  # 123 + len("ABCDEFGHIJ") - 1
        self.assertEqual(first.motif, "MEME-1")

    def test_write_outputs(self) -> None:
        motifs, locations = parse_meme(self.meme_file)
        outdir = self.tmp / "out"
        info_path, loc_path = write_outputs(motifs, locations, outdir)

        self.assertTrue(info_path.exists())
        self.assertTrue(loc_path.exists())

        info_lines = info_path.read_text().strip().split("\n")
        self.assertEqual(len(info_lines), 2)
        self.assertEqual(info_lines[0], "MEME-1\tABCDEFGHIJ\t10")

        loc_lines = loc_path.read_text().strip().split("\n")
        self.assertEqual(len(loc_lines), 5)
        cols = loc_lines[0].split("\t")
        self.assertEqual(len(cols), 4)

    def test_empty_meme_file(self) -> None:
        empty = self.tmp / "empty.txt"
        empty.write_text("nothing useful here\n")
        motifs, locations = parse_meme(empty)
        self.assertEqual(len(motifs), 0)
        self.assertEqual(len(locations), 0)

    def test_cli_main(self) -> None:
        outdir = self.tmp / "cli_out"
        rc = main([str(self.meme_file), "-o", str(outdir)])
        self.assertEqual(rc, 0)
        self.assertTrue((outdir / "meme_info.txt").exists())
        self.assertTrue((outdir / "meme_location.txt").exists())

    def test_cli_missing_file(self) -> None:
        rc = main(["/nonexistent/meme.txt", "-o", str(self.tmp)])
        self.assertEqual(rc, 1)

    def test_parse_with_fasta_whitelist(self) -> None:
        motifs, locations = parse_meme(
            self.meme_file,
            allowed_ids={"GeneA.001G000100.1.p", "GeneC.003G000300.1.p"},
        )
        self.assertEqual(len(motifs), 2)
        self.assertEqual(len(locations), 4)
        self.assertTrue(all(
            loc.gene in {"GeneA.001G000100.1.p", "GeneC.003G000300.1.p"}
            for loc in locations
        ))

    def test_load_fasta_ids(self) -> None:
        fasta = self.tmp / "in.fa"
        fasta.write_text(
            ">GeneA.001G000100.1.p some description\nACDEFG\n"
            ">GeneB.002G000200.2.p\nHIJKLM\n"
            ">GeneZ\nNNNN\n"
        )
        ids = load_fasta_ids(fasta)
        self.assertEqual(
            ids,
            {"GeneA.001G000100.1.p", "GeneB.002G000200.2.p", "GeneZ"},
        )

    def test_cli_with_fasta(self) -> None:
        fasta = self.tmp / "in.fa"
        fasta.write_text(
            ">GeneB.002G000200.2.p\nACDE\n"
            ">GeneC.003G000300.1.p\nFGHI\n"
        )
        outdir = self.tmp / "with_fasta"
        rc = main([
            str(self.meme_file), "-o", str(outdir),
            "--fasta", str(fasta),
        ])
        self.assertEqual(rc, 0)
        loc_lines = (outdir / "meme_location.txt").read_text().strip().split("\n")
        genes = {line.split("\t")[0] for line in loc_lines}
        self.assertEqual(
            genes,
            {"GeneB.002G000200.2.p", "GeneC.003G000300.1.p"},
        )

    def test_cli_with_gene_filter(self) -> None:
        outdir = self.tmp / "filtered"
        rc = main([
            str(self.meme_file), "-o", str(outdir),
            "--gene-id-pattern", r"GeneB\.",
        ])
        self.assertEqual(rc, 0)
        loc_lines = (outdir / "meme_location.txt").read_text().strip().split("\n")
        self.assertEqual(len(loc_lines), 1)
        self.assertIn("GeneB", loc_lines[0])

    def test_cli_missing_fasta(self) -> None:
        """--fasta pointing at a nonexistent file should exit 1."""
        outdir = self.tmp / "no_fasta"
        rc = main([
            str(self.meme_file), "-o", str(outdir),
            "--fasta", str(self.tmp / "does_not_exist.fa"),
        ])
        self.assertEqual(rc, 1)

    def test_cli_empty_fasta(self) -> None:
        """--fasta pointing at a file with no sequence IDs should exit 1."""
        empty_fa = self.tmp / "empty.fa"
        empty_fa.write_text("# no records\n")
        outdir = self.tmp / "empty_fasta"
        rc = main([
            str(self.meme_file), "-o", str(outdir),
            "--fasta", str(empty_fa),
        ])
        self.assertEqual(rc, 1)

    def test_cli_parse_error(self) -> None:
        """Exceptions raised inside parse_meme should surface as exit 1."""
        import scripts.parse_meme_output as mod

        def boom(*args, **kwargs):
            raise RuntimeError("synthetic parse failure")

        original = mod.parse_meme
        mod.parse_meme = boom
        try:
            rc = main([str(self.meme_file), "-o", str(self.tmp / "err_out")])
        finally:
            mod.parse_meme = original
        self.assertEqual(rc, 1)


if __name__ == "__main__":
    unittest.main()
