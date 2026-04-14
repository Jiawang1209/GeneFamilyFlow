"""Tests for scripts/fetch_jaspar_plants.py.

Covers the offline parts only — no real network. Network fetches are
exercised indirectly via the ``--from-file`` path so CI stays hermetic.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.fetch_jaspar_plants import (
    cache_is_valid,
    copy_from_file,
    fetch,
    looks_like_meme,
    main,
    sha256_of,
)

pytestmark = pytest.mark.unit

MEME_HEADER = (
    "MEME version 4\n"
    "\n"
    "ALPHABET= ACGT\n"
    "\n"
    "strands: + -\n"
    "\n"
    "MOTIF MA0001.1 Demo\n"
    "letter-probability matrix: alength= 4 w= 4\n"
    " 0.25 0.25 0.25 0.25\n"
    " 0.25 0.25 0.25 0.25\n"
    " 0.25 0.25 0.25 0.25\n"
    " 0.25 0.25 0.25 0.25\n"
)


def _write_meme(path: Path, body: str = MEME_HEADER) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body)
    return path


class TestLooksLikeMeme:
    def test_valid_header(self, tmp_path: Path) -> None:
        p = _write_meme(tmp_path / "ok.meme")
        assert looks_like_meme(p) is True

    def test_missing_file(self, tmp_path: Path) -> None:
        assert looks_like_meme(tmp_path / "nope.meme") is False

    def test_empty_file(self, tmp_path: Path) -> None:
        p = tmp_path / "empty.meme"
        p.write_text("")
        assert looks_like_meme(p) is False

    def test_wrong_header(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.meme"
        p.write_text("not a meme bundle\n")
        assert looks_like_meme(p) is False


class TestSha256:
    def test_sha256_stable(self, tmp_path: Path) -> None:
        p = tmp_path / "x.bin"
        p.write_bytes(b"hello world")
        # SHA256 of "hello world"
        assert sha256_of(p) == (
            "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
        )


class TestCacheIsValid:
    def test_no_sha_required(self, tmp_path: Path) -> None:
        p = _write_meme(tmp_path / "a.meme")
        assert cache_is_valid(p, expected_sha256=None) is True

    def test_sha_match(self, tmp_path: Path) -> None:
        p = _write_meme(tmp_path / "a.meme")
        assert cache_is_valid(p, expected_sha256=sha256_of(p)) is True

    def test_sha_mismatch(self, tmp_path: Path) -> None:
        p = _write_meme(tmp_path / "a.meme")
        assert cache_is_valid(p, expected_sha256="0" * 64) is False

    def test_missing_file_invalid(self, tmp_path: Path) -> None:
        assert cache_is_valid(tmp_path / "nope.meme", expected_sha256=None) is False


class TestCopyFromFile:
    def test_copy_creates_parent(self, tmp_path: Path) -> None:
        src = _write_meme(tmp_path / "src.meme")
        dst = tmp_path / "nested" / "dst.meme"
        copy_from_file(src, dst)
        assert dst.exists()
        assert dst.read_text() == src.read_text()

    def test_missing_source_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            copy_from_file(tmp_path / "nope.meme", tmp_path / "dst.meme")


class TestFetchFromFile:
    def test_fetch_skips_when_cache_valid(self, tmp_path: Path) -> None:
        existing = _write_meme(tmp_path / "out.meme")
        original_text = existing.read_text()
        # from_file is provided but should NOT be used because cache is valid
        sentinel = _write_meme(tmp_path / "other.meme", body=MEME_HEADER + "MOTIF MA9999.1 Other\n")
        out = fetch(output=existing, from_file=sentinel)
        assert out == existing
        assert existing.read_text() == original_text

    def test_fetch_from_file_writes_dst(self, tmp_path: Path) -> None:
        src = _write_meme(tmp_path / "src.meme")
        dst = tmp_path / "fresh" / "dst.meme"
        out = fetch(output=dst, from_file=src)
        assert out == dst
        assert dst.exists()
        assert looks_like_meme(dst)

    def test_fetch_force_overwrites(self, tmp_path: Path) -> None:
        existing = _write_meme(tmp_path / "out.meme", body=MEME_HEADER + "MOTIF X Y\n")
        replacement = _write_meme(tmp_path / "src.meme", body=MEME_HEADER)
        fetch(output=existing, from_file=replacement, force=True)
        assert existing.read_text() == MEME_HEADER

    def test_fetch_rejects_garbage_payload(self, tmp_path: Path) -> None:
        garbage = tmp_path / "garbage.meme"
        garbage.write_text("not a meme bundle at all\n")
        dst = tmp_path / "dst.meme"
        with pytest.raises(RuntimeError, match="MEME"):
            fetch(output=dst, from_file=garbage)

    def test_fetch_sha_mismatch_raises(self, tmp_path: Path) -> None:
        src = _write_meme(tmp_path / "src.meme")
        dst = tmp_path / "dst.meme"
        with pytest.raises(RuntimeError, match="SHA256"):
            fetch(output=dst, from_file=src, expected_sha256="0" * 64)


class TestCli:
    def test_cli_from_file_roundtrip(self, tmp_path: Path) -> None:
        src = _write_meme(tmp_path / "src.meme")
        dst = tmp_path / "out.meme"
        rc = main(["--output", str(dst), "--from-file", str(src)])
        assert rc == 0
        assert dst.exists()
        assert looks_like_meme(dst)


class TestShippedBundle:
    def test_shipped_bundle_parses(self) -> None:
        repo_root = Path(__file__).resolve().parent.parent
        bundle = repo_root / "example" / "10.promoter" / "JASPAR2024_plants.meme"
        assert bundle.exists(), f"missing shipped JASPAR bundle: {bundle}"
        assert looks_like_meme(bundle)
        # Must contain enough motifs to be a credible plants bundle.
        text = bundle.read_text()
        assert text.count("\nMOTIF ") + text.startswith("MOTIF ") > 100
