"""Tests for scripts/build_kegg_names.py.

Covers:
- KEGG pathway.list parsing + normalization (map -> ko prefix)
- eggNOG annotation scanning for ko term extraction
- Merge semantics: preserve existing names, add new empties, skip overwrites
- CLI round-trip for both input modes
- Fixture integrity: example/13.GO_KEGG/kegg_pathway_names.tsv is parseable,
  non-empty, and every term matches the canonical ko00000 pattern.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from scripts.build_kegg_names import (
    PathwayName,
    _normalize_term,
    extract_eggnog_terms,
    main,
    merge_entries,
    parse_kegg_list,
    read_existing,
    write_tsv,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURE_TSV = REPO_ROOT / "example" / "13.GO_KEGG" / "kegg_pathway_names.tsv"

KO_PATTERN = re.compile(r"^ko\d{5}$")


EMAPPER_HEADER = (
    "## emapper-2.1.12\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
    "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t"
    "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
    "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
)
EMAPPER_ROW = (
    "gene{gid}\t-\t1e-50\t100\t-\t-\t-\t-\t-\tGO:0008150\t-\t-\t"
    "{kegg}\t-\t-\t-\t-\t-\t-\t-\t-\n"
)


def _write_kegg_list(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "pathway.list"
    path.write_text("\n".join(lines) + "\n")
    return path


def _write_emapper(tmp_path: Path, name: str, rows: list[tuple[str, str]]) -> Path:
    path = tmp_path / name
    body = "".join(EMAPPER_ROW.format(gid=g, kegg=k) for g, k in rows)
    path.write_text(EMAPPER_HEADER + body)
    return path


class TestNormalizeTerm:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("map00010", "ko00010"),
            ("ko04075", "ko04075"),
            ("path:map04075", "ko04075"),
            ("path:ko00010", "ko00010"),
            ("  map00010  ", "ko00010"),
        ],
    )
    def test_valid_forms(self, raw: str, expected: str) -> None:
        assert _normalize_term(raw) == expected

    @pytest.mark.parametrize(
        "raw",
        ["", "PF00854", "ko", "ko1234", "ko123456", "foo00010", "mapABCDE"],
    )
    def test_rejects_garbage(self, raw: str) -> None:
        assert _normalize_term(raw) is None


class TestParseKeggList:
    def test_parses_and_normalizes(self, tmp_path: Path) -> None:
        f = _write_kegg_list(
            tmp_path,
            [
                "map00010\tGlycolysis / Gluconeogenesis",
                "map04075\tPlant hormone signal transduction",
                "ko00500\tStarch and sucrose metabolism",
            ],
        )
        entries = parse_kegg_list(f)
        assert [e.term for e in entries] == ["ko00010", "ko04075", "ko00500"]
        assert entries[0].name == "Glycolysis / Gluconeogenesis"

    def test_skips_blank_and_garbage(self, tmp_path: Path) -> None:
        f = _write_kegg_list(
            tmp_path,
            [
                "",
                "# comment",
                "not_a_pathway\tignored",
                "map00010\tGlycolysis",
                "map00020\t   ",  # empty name
            ],
        )
        entries = parse_kegg_list(f)
        assert len(entries) == 1
        assert entries[0].term == "ko00010"


class TestExtractEggnogTerms:
    def test_union_of_annotations(self, tmp_path: Path) -> None:
        a = _write_emapper(tmp_path, "a.annot", [("1", "ko00010,ko04075"), ("2", "map00500")])
        b = _write_emapper(tmp_path, "b.annot", [("3", "ko00010,ko00020")])
        terms = extract_eggnog_terms([a, b])
        assert terms == ["ko00010", "ko00020", "ko00500", "ko04075"]

    def test_ignores_missing_and_dash(self, tmp_path: Path) -> None:
        a = _write_emapper(
            tmp_path,
            "a.annot",
            [("1", "-"), ("2", ""), ("3", "ko04075,NotAPathway")],
        )
        terms = extract_eggnog_terms([a])
        assert terms == ["ko04075"]

    def test_missing_column_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "bad.annot"
        f.write_text("#query\tGOs\ngene1\tGO:0008150\n")
        with pytest.raises(ValueError, match="KEGG_Pathway"):
            extract_eggnog_terms([f])


class TestMergeEntries:
    def test_preserves_existing_names(self) -> None:
        existing = {"ko00010": "Glycolysis / Gluconeogenesis"}
        additions = [PathwayName(term="ko00010", name="SOMETHING ELSE")]
        merged = merge_entries(existing, additions)
        assert merged == [PathwayName(term="ko00010", name="Glycolysis / Gluconeogenesis")]

    def test_fills_empty_existing_name(self) -> None:
        existing = {"ko00010": ""}
        additions = [PathwayName(term="ko00010", name="Glycolysis / Gluconeogenesis")]
        merged = merge_entries(existing, additions)
        assert merged[0].name == "Glycolysis / Gluconeogenesis"

    def test_appends_new_terms_sorted(self) -> None:
        existing = {"ko00010": "Glycolysis"}
        additions = [
            PathwayName(term="ko04075", name="Plant hormone signal transduction"),
            PathwayName(term="ko00020", name="TCA"),
        ]
        merged = merge_entries(existing, additions)
        assert [e.term for e in merged] == ["ko00010", "ko00020", "ko04075"]


class TestReadExisting:
    def test_roundtrip(self, tmp_path: Path) -> None:
        p = tmp_path / "kegg.tsv"
        write_tsv([PathwayName(term="ko00010", name="Glycolysis")], p)
        assert read_existing(p) == {"ko00010": "Glycolysis"}

    def test_missing_file_returns_empty(self, tmp_path: Path) -> None:
        assert read_existing(tmp_path / "absent.tsv") == {}


class TestCLI:
    def test_from_kegg_list_roundtrip(self, tmp_path: Path) -> None:
        src = _write_kegg_list(
            tmp_path, ["map00010\tGlycolysis", "ko04075\tPlant hormone signal transduction"]
        )
        out = tmp_path / "kegg.tsv"
        rc = main(["--from-kegg-list", str(src), "-o", str(out)])
        assert rc == 0
        body = out.read_text().splitlines()
        assert body[0] == "term\tname"
        assert "ko00010\tGlycolysis" in body
        assert "ko04075\tPlant hormone signal transduction" in body

    def test_from_eggnog_scaffold(self, tmp_path: Path) -> None:
        annot = _write_emapper(
            tmp_path, "a.annot", [("1", "ko00010"), ("2", "map04075")]
        )
        out = tmp_path / "scaffold.tsv"
        rc = main(["--from-eggnog", str(annot), "-o", str(out)])
        assert rc == 0
        body = out.read_text().splitlines()[1:]
        assert body == ["ko00010\t", "ko04075\t"]

    def test_merge_preserves_prior_entries(self, tmp_path: Path) -> None:
        base = tmp_path / "base.tsv"
        write_tsv([PathwayName(term="ko00010", name="Glycolysis / Gluconeogenesis")], base)
        src = _write_kegg_list(tmp_path, ["map00020\tCitrate cycle"])
        out = tmp_path / "out.tsv"
        rc = main(
            [
                "--from-kegg-list", str(src),
                "--merge", str(base),
                "-o", str(out),
            ]
        )
        assert rc == 0
        assert read_existing(out) == {
            "ko00010": "Glycolysis / Gluconeogenesis",
            "ko00020": "Citrate cycle",
        }


class TestShippedFixture:
    def test_fixture_exists_and_parses(self) -> None:
        assert FIXTURE_TSV.exists(), f"missing fixture {FIXTURE_TSV}"
        entries = read_existing(FIXTURE_TSV)
        assert entries, "fixture must be non-empty"

    def test_all_terms_are_canonical_ko(self) -> None:
        for term in read_existing(FIXTURE_TSV):
            assert KO_PATTERN.match(term), f"bad term: {term!r}"

    def test_names_are_non_empty(self) -> None:
        for term, name in read_existing(FIXTURE_TSV).items():
            assert name.strip(), f"empty name for {term}"

    def test_covers_core_plant_pathways(self) -> None:
        entries = read_existing(FIXTURE_TSV)
        for must_have in ("ko04075", "ko04626", "ko00195", "ko00940", "ko04016"):
            assert must_have in entries, f"fixture missing {must_have}"
