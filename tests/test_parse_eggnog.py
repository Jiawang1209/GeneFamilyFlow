"""Unit tests for scripts/parse_eggnog.py."""

from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import pytest

from scripts.parse_eggnog import (
    EggnogRow,
    _normalize_kegg,
    _split_list,
    main,
    parse_eggnog_annotations,
)


SAMPLE = dedent(
    """\
    ## eggNOG-mapper v2.1.9
    ## time=2026-04-14
    #query\tseed_ortholog\tevalue\tGOs\tKEGG_Pathway
    AT1G01010\tortho1\t1e-50\tGO:0003700,GO:0005634\tmap04075,ko04075
    AT1G01020\tortho2\t1e-40\t-\tko00010
    AT1G01030\tortho3\t1e-30\tGO:0005515\t-
    AT1G01040\tortho4\t1e-20\t-\t-
    """
)


def _write_sample(tmp_path: Path, content: str = SAMPLE) -> Path:
    path = tmp_path / "toy.emapper.annotations"
    path.write_text(content)
    return path


def test_split_list_handles_missing():
    assert _split_list("-") == ()
    assert _split_list("") == ()
    assert _split_list("NA") == ()
    assert _split_list("GO:1,GO:2 , -,GO:3") == ("GO:1", "GO:2", "GO:3")


def test_normalize_kegg_dedupes_and_prefers_ko():
    assert _normalize_kegg(["map04075", "ko04075", "map00010"]) == ("ko04075", "ko00010")
    assert _normalize_kegg([]) == ()
    assert _normalize_kegg(["notapath"]) == ()


def test_parse_eggnog_annotations_returns_expected_rows(tmp_path):
    path = _write_sample(tmp_path)
    rows = parse_eggnog_annotations(path)
    assert len(rows) == 4

    assert rows[0] == EggnogRow(
        gene="AT1G01010",
        go_terms=("GO:0003700", "GO:0005634"),
        kegg_pathways=("ko04075",),
    )
    assert rows[1].go_terms == ()
    assert rows[1].kegg_pathways == ("ko00010",)
    assert rows[3].go_terms == ()
    assert rows[3].kegg_pathways == ()


def test_parse_raises_on_missing_required_columns(tmp_path):
    bad = tmp_path / "bad.annotations"
    bad.write_text("#query\tseed_ortholog\nAT1G01010\tx\n")
    with pytest.raises(ValueError, match="missing columns"):
        parse_eggnog_annotations(bad)


def test_parse_raises_on_empty_file(tmp_path):
    empty = tmp_path / "empty.annotations"
    empty.write_text("## only a comment\n")
    with pytest.raises(ValueError, match="No data lines"):
        parse_eggnog_annotations(empty)


class TestShippedFixtures:
    """The bundled example/13.GO_KEGG/*.emapper.annotations files must stay
    parseable and must cover every family ID in the matching NPF fixture so
    Step 13 renders without an empty gene set when users flip
    ``step13_go_kegg.enabled`` to ``true`` out of the box."""

    REPO_ROOT = Path(__file__).resolve().parent.parent

    @pytest.mark.parametrize(
        "species,npf_id_file",
        [
            ("Athaliana", "example/9.mcscanx/AT.NPF.id"),
            ("Osativa", "example/9.mcscanx/Os.NPF.id"),
        ],
    )
    def test_fixture_parses_and_covers_family(self, species, npf_id_file):
        ann_path = self.REPO_ROOT / f"example/13.GO_KEGG/{species}.emapper.annotations"
        assert ann_path.exists(), f"missing shipped fixture: {ann_path}"

        rows = parse_eggnog_annotations(ann_path)
        assert len(rows) >= 50, "fixture should have enough rows to render enrichment"

        genes_with_annotations = {r.gene for r in rows if r.go_terms or r.kegg_pathways}
        family_ids = {
            line.strip()
            for line in (self.REPO_ROOT / npf_id_file).read_text().splitlines()
            if line.strip()
        }
        missing = family_ids - genes_with_annotations
        assert not missing, (
            f"{len(missing)} family IDs are absent from {ann_path.name}; "
            f"Step 13 would silently drop them. First few: {sorted(missing)[:3]}"
        )

        all_go = {t for r in rows for t in r.go_terms}
        all_kegg = {t for r in rows for t in r.kegg_pathways}
        assert len(all_go) >= 5, "GO term diversity too low for faceted dotplot"
        assert len(all_kegg) >= 2, "KEGG term diversity too low for enrichment"


def test_main_writes_all_outputs(tmp_path):
    ann = _write_sample(tmp_path)
    go_out = tmp_path / "go.tsv"
    kegg_out = tmp_path / "kegg.tsv"
    universe_out = tmp_path / "universe.txt"

    rc = main([
        str(ann),
        "--go-out", str(go_out),
        "--kegg-out", str(kegg_out),
        "--universe-out", str(universe_out),
    ])
    assert rc == 0

    go_lines = go_out.read_text().strip().splitlines()
    assert go_lines[0] == "term\tgene"
    assert "GO:0003700\tAT1G01010" in go_lines
    assert "GO:0005634\tAT1G01010" in go_lines
    assert "GO:0005515\tAT1G01030" in go_lines
    assert len(go_lines) == 4  # header + 3 pairs

    kegg_lines = kegg_out.read_text().strip().splitlines()
    assert kegg_lines[0] == "term\tgene"
    assert "ko04075\tAT1G01010" in kegg_lines
    assert "ko00010\tAT1G01020" in kegg_lines
    assert len(kegg_lines) == 3  # header + 2 pairs (dedup map/ko)

    universe = universe_out.read_text().strip().splitlines()
    # AT1G01040 has no annotations -> excluded
    assert set(universe) == {"AT1G01010", "AT1G01020", "AT1G01030"}


def test_parse_eggnog_skips_blank_and_truncated_and_empty_gene(tmp_path: Path) -> None:
    """Defensive branches: blank lines, too-few-fields, empty query."""
    content = dedent(
        """\
        ## eggNOG-mapper v2.1.9

        #query\tseed_ortholog\tevalue\tGOs\tKEGG_Pathway
        AT1G99000\tortho\t1e-20\tGO:0000001\tko99999
        short_line_only_three_cols\ta\tb

        \t-\t-\t-\t-
        AT1G99001\tortho\t1e-20\t-\t-
        """
    )
    path = _write_sample(tmp_path, content)
    rows = parse_eggnog_annotations(path)
    genes = {r.gene for r in rows}
    # Blank lines skipped (line 38), truncated row skipped (line 102),
    # row with empty query skipped (line 105). Valid rows survive.
    assert genes == {"AT1G99000", "AT1G99001"}


def test_write_universe_deduplicates_gene(tmp_path: Path) -> None:
    """When two rows share the same gene, universe keeps only one entry (line 135)."""
    from scripts.parse_eggnog import write_universe

    rows = [
        EggnogRow(gene="G1", go_terms=("GO:1",), kegg_pathways=()),
        EggnogRow(gene="G1", go_terms=(), kegg_pathways=("ko1",)),
        EggnogRow(gene="G2", go_terms=("GO:2",), kegg_pathways=()),
    ]
    out = tmp_path / "universe.txt"
    n = write_universe(rows, out)
    assert n == 2
    assert out.read_text().splitlines() == ["G1", "G2"]
