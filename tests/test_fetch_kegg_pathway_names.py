"""Tests for scripts/fetch_kegg_pathway_names.py.

The network ``fetch()`` helper is intentionally not exercised — the
script declares it as ``# pragma: no cover``. We unit-test the offline
parts (parse/write) which carry all the normalization logic.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.fetch_kegg_pathway_names import parse, write

pytestmark = pytest.mark.unit


class TestParse:
    def test_strips_path_prefix(self) -> None:
        rows = parse("path:map04075\tPlant hormone signal transduction\n")
        assert rows == [("ko04075", "Plant hormone signal transduction")]

    def test_handles_already_ko_prefixed(self) -> None:
        rows = parse("path:ko04075\tHormone\n")
        # path: stripped, ko prefix preserved as-is
        assert rows == [("ko04075", "Hormone")]

    def test_ignores_blank_and_malformed(self) -> None:
        text = (
            "\n"
            "path:map04075\tHormone\n"
            "garbage_no_tab\n"
            "\t\n"
            "path:map00010\tGlycolysis\n"
        )
        rows = parse(text)
        assert rows == [
            ("ko04075", "Hormone"),
            ("ko00010", "Glycolysis"),
        ]

    def test_strips_whitespace(self) -> None:
        rows = parse("  path:map04075 \t  Hormone  \n")
        assert rows == [("ko04075", "Hormone")]

    def test_passes_through_non_map_terms(self) -> None:
        # Defensive: KEGG normally only emits map* under /list/pathway,
        # but if it ever returns something else we should not silently
        # rewrite the prefix.
        rows = parse("foo\tBar\n")
        assert rows == [("foo", "Bar")]


class TestWrite:
    def test_writes_header_and_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "out.tsv"
        n = write(
            [("ko04075", "Hormone"), ("ko00010", "Glycolysis")],
            out,
        )
        assert n == 2
        lines = out.read_text().splitlines()
        assert lines[0] == "term\tname"
        assert lines[1] == "ko04075\tHormone"
        assert lines[2] == "ko00010\tGlycolysis"

    def test_creates_parent_dir(self, tmp_path: Path) -> None:
        out = tmp_path / "deep" / "nested" / "out.tsv"
        write([], out)
        assert out.exists()
        assert out.read_text() == "term\tname\n"


class TestMainCli:
    """Cover ``main()`` by stubbing ``fetch`` at the module level. Per
    ``docs/development.md`` the ``urllib`` body in ``fetch`` itself is
    policy-exempt; this test sits above it."""

    def test_main_writes_parsed_table(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
    ) -> None:
        import scripts.fetch_kegg_pathway_names as mod

        canned = (
            "path:map04075\tPlant hormone signal transduction\n"
            "path:map00010\tGlycolysis / Gluconeogenesis\n"
        )
        monkeypatch.setattr(mod, "fetch", lambda *a, **kw: canned)

        out = tmp_path / "kegg_names.tsv"
        rc = mod.main(["-o", str(out)])
        assert rc == 0

        lines = out.read_text().splitlines()
        assert lines[0] == "term\tname"
        assert lines[1] == "ko04075\tPlant hormone signal transduction"
        assert lines[2] == "ko00010\tGlycolysis / Gluconeogenesis"
        assert f"wrote 2 pathways to {out}" in capsys.readouterr().err

    def test_main_creates_missing_output_dir(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import scripts.fetch_kegg_pathway_names as mod

        monkeypatch.setattr(mod, "fetch", lambda *a, **kw: "path:map00010\tGlycolysis\n")
        out = tmp_path / "deep" / "nested" / "kegg.tsv"
        rc = mod.main(["-o", str(out)])
        assert rc == 0
        assert out.exists()


class TestModuleEntrypoint:
    def test_help_runs_via_dash_m(self) -> None:
        """``python -m scripts.fetch_kegg_pathway_names --help`` covers
        the ``raise SystemExit(main())`` tail without hitting the network
        — argparse short-circuits with ``SystemExit(0)``."""
        import runpy
        import sys

        original = sys.argv
        sys.argv = ["fetch_kegg_pathway_names.py", "--help"]
        try:
            with pytest.raises(SystemExit) as exc:
                runpy.run_module(
                    "scripts.fetch_kegg_pathway_names", run_name="__main__"
                )
            assert exc.value.code == 0
        finally:
            sys.argv = original
