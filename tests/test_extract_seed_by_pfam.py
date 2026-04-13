"""Tests for scripts/extract_seed_by_pfam.py.

The script consumes a clean 2-column TSV (gene_id<TAB>pfam_id) and a
gene-level FASTA, extracts gene IDs hitting the requested Pfam accession,
and writes matching sequences to an output FASTA.
"""

from __future__ import annotations

import textwrap
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.extract_seed_by_pfam import (
    extract_sequences,
    main,
    parse_domain_hits,
    write_seed_fasta,
)

SAMPLE_DOMAINS = textwrap.dedent("""\
    AT1G00001\tPF00854
    AT1G00001\tPF07690
    AT1G00002\tPF00854
    AT1G00003\tPF01234
    AT1G00004\tPF00005
""")

SAMPLE_FASTA = textwrap.dedent("""\
    >AT1G00001
    MSEQAAABBBCCCDDDEEEFFFGGGHHH
    IIIJJJKKKLLLMMMNNN
    >AT1G00002
    MSEQZZZYYYXXXWWWVVVUUU
    >AT1G00003
    MSEQAAACCCEEE
    >AT1G00099
    MSEQORPHAN
""")

REPO_ROOT = Path(__file__).resolve().parent.parent
AT_DOMAINS = REPO_ROOT / "example/3.blast/references/Athaliana.sample.domains.tsv"
AT_FASTA   = REPO_ROOT / "example/3.blast/references/Athaliana.sample.pep.fasta"
OS_DOMAINS = REPO_ROOT / "example/3.blast/references/Osativa.sample.domains.tsv"
OS_FASTA   = REPO_ROOT / "example/3.blast/references/Osativa.sample.pep.fasta"


class TestParseDomainHits(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.domains = self.tmp / "domains.tsv"
        self.domains.write_text(SAMPLE_DOMAINS)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_filters_by_pfam_id(self) -> None:
        ids = parse_domain_hits(self.domains, "PF00854")
        self.assertEqual(ids, ["AT1G00001", "AT1G00002"])

    def test_dedups_gene_ids(self) -> None:
        # AT1G00001 has two Pfam hits; only one gene ID should be returned.
        ids = parse_domain_hits(self.domains, "PF00854")
        self.assertEqual(len(ids), len(set(ids)))

    def test_returns_sorted(self) -> None:
        ids = parse_domain_hits(self.domains, "PF00854")
        self.assertEqual(ids, sorted(ids))

    def test_different_pfam(self) -> None:
        self.assertEqual(parse_domain_hits(self.domains, "PF01234"), ["AT1G00003"])
        self.assertEqual(parse_domain_hits(self.domains, "PF07690"), ["AT1G00001"])

    def test_nonexistent_pfam_returns_empty(self) -> None:
        self.assertEqual(parse_domain_hits(self.domains, "PF99999"), [])

    def test_empty_file(self) -> None:
        empty = self.tmp / "empty.tsv"
        empty.write_text("")
        self.assertEqual(parse_domain_hits(empty, "PF00854"), [])

    def test_file_not_found(self) -> None:
        with self.assertRaises(FileNotFoundError):
            parse_domain_hits(self.tmp / "missing.tsv", "PF00854")

    def test_ignores_blank_and_comment_lines(self) -> None:
        noisy = self.tmp / "noisy.tsv"
        noisy.write_text("# header comment\n\nAT1G00001\tPF00854\n")
        self.assertEqual(parse_domain_hits(noisy, "PF00854"), ["AT1G00001"])


class TestExtractSequences(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.fasta = self.tmp / "proteome.fasta"
        self.fasta.write_text(SAMPLE_FASTA)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_extracts_matching_gene_ids(self) -> None:
        records = list(extract_sequences(self.fasta, ["AT1G00001", "AT1G00002"]))
        self.assertEqual(sorted(h for h, _ in records), ["AT1G00001", "AT1G00002"])

    def test_sequence_joined_across_lines(self) -> None:
        records = dict(extract_sequences(self.fasta, ["AT1G00001"]))
        self.assertEqual(
            records["AT1G00001"],
            "MSEQAAABBBCCCDDDEEEFFFGGGHHHIIIJJJKKKLLLMMMNNN",
        )

    def test_missing_gene_id_silently_skipped(self) -> None:
        self.assertEqual(list(extract_sequences(self.fasta, ["AT9G99999"])), [])

    def test_only_requested_genes_returned(self) -> None:
        records = dict(extract_sequences(self.fasta, ["AT1G00002"]))
        self.assertIn("AT1G00002", records)
        self.assertNotIn("AT1G00001", records)
        self.assertNotIn("AT1G00099", records)

    def test_empty_gene_id_list(self) -> None:
        self.assertEqual(list(extract_sequences(self.fasta, [])), [])


class TestWriteSeedFasta(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.fasta = self.tmp / "proteome.fasta"
        self.fasta.write_text(SAMPLE_FASTA)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_writes_fasta_format(self) -> None:
        out = self.tmp / "seed.fa"
        count = write_seed_fasta(self.fasta, ["AT1G00001", "AT1G00002"], out)
        self.assertEqual(count, 2)
        content = out.read_text()
        self.assertEqual(content.count(">"), 2)
        self.assertIn(">AT1G00001", content)
        self.assertIn(">AT1G00002", content)

    def test_creates_parent_directory(self) -> None:
        out = self.tmp / "nested/dir/seed.fa"
        count = write_seed_fasta(self.fasta, ["AT1G00001"], out)
        self.assertEqual(count, 1)
        self.assertTrue(out.exists())

    def test_empty_output_when_no_matches(self) -> None:
        out = self.tmp / "seed.fa"
        count = write_seed_fasta(self.fasta, ["AT9G99999"], out)
        self.assertEqual(count, 0)
        self.assertTrue(out.exists())
        self.assertEqual(out.read_text(), "")


class TestCli(unittest.TestCase):

    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.domains = self.tmp / "domains.tsv"
        self.fasta = self.tmp / "proteome.fasta"
        self.domains.write_text(SAMPLE_DOMAINS)
        self.fasta.write_text(SAMPLE_FASTA)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_cli_end_to_end(self) -> None:
        out = self.tmp / "seed.fa"
        rc = main([
            "--domains-table", str(self.domains),
            "--proteome", str(self.fasta),
            "--pfam-id", "PF00854",
            "--output", str(out),
        ])
        self.assertEqual(rc, 0)
        content = out.read_text()
        self.assertIn(">AT1G00001", content)
        self.assertIn(">AT1G00002", content)
        self.assertNotIn(">AT1G00003", content)

    def test_cli_missing_domains_errors(self) -> None:
        out = self.tmp / "seed.fa"
        rc = main([
            "--domains-table", str(self.tmp / "missing.tsv"),
            "--proteome", str(self.fasta),
            "--pfam-id", "PF00854",
            "--output", str(out),
        ])
        self.assertEqual(rc, 1)

    def test_cli_missing_proteome_errors(self) -> None:
        out = self.tmp / "seed.fa"
        rc = main([
            "--domains-table", str(self.domains),
            "--proteome", str(self.tmp / "missing.fa"),
            "--pfam-id", "PF00854",
            "--output", str(out),
        ])
        self.assertEqual(rc, 1)


@unittest.skipUnless(
    AT_DOMAINS.exists() and AT_FASTA.exists(),
    "Athaliana sample fixtures not present",
)
class TestAthalianaSample(unittest.TestCase):
    """Integration against committed Athaliana.sample.* fixtures."""

    def test_pf00854_matches_expected_genes(self) -> None:
        ids = parse_domain_hits(AT_DOMAINS, "PF00854")
        self.assertEqual(
            ids,
            ["AT1G12110", "AT1G18880", "AT1G22540", "AT1G22550", "AT1G22570"],
        )

    def test_pf00854_extracts_5_sequences(self) -> None:
        with TemporaryDirectory() as td:
            out = Path(td) / "seed.fa"
            rc = main([
                "--domains-table", str(AT_DOMAINS),
                "--proteome", str(AT_FASTA),
                "--pfam-id", "PF00854",
                "--output", str(out),
            ])
            self.assertEqual(rc, 0)
            self.assertEqual(out.read_text().count(">"), 5)


@unittest.skipUnless(
    OS_DOMAINS.exists() and OS_FASTA.exists(),
    "Osativa sample fixtures not present",
)
class TestOsativaSample(unittest.TestCase):
    """Integration against committed Osativa.sample.* fixtures."""

    def test_pf00854_matches_expected_genes(self) -> None:
        ids = parse_domain_hits(OS_DOMAINS, "PF00854")
        self.assertEqual(
            ids,
            [
                "LOC_Os01g01360",
                "LOC_Os01g04950",
                "LOC_Os01g37590",
                "LOC_Os01g54515",
                "LOC_Os01g55600",
            ],
        )

    def test_pf00854_extracts_5_sequences(self) -> None:
        with TemporaryDirectory() as td:
            out = Path(td) / "seed.fa"
            rc = main([
                "--domains-table", str(OS_DOMAINS),
                "--proteome", str(OS_FASTA),
                "--pfam-id", "PF00854",
                "--output", str(out),
            ])
            self.assertEqual(rc, 0)
            self.assertEqual(out.read_text().count(">"), 5)


if __name__ == "__main__":
    unittest.main()
