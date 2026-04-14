#!/usr/bin/env python3
"""Fetch the JASPAR 2024 CORE plants non-redundant motif bundle.

Downloads ``JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt`` from the
JASPAR ELIXIR Norway mirror and writes it to a local cache path, where
``scripts/scan_promoter_fimo.py`` consumes it. Designed to be:

- Idempotent: skips the download when the cache already holds a valid
  copy (size > 0, parses as MEME, optional SHA256 match).
- Offline-friendly: ``--from-file`` copies a user-supplied bundle instead
  of touching the network — useful in air-gapped clusters and CI.
- Safe to re-run: writes to a temp file then atomically renames.

Default URL (stable JASPAR mirror):
    https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt

The bundle ships ~805 plant transcription factor PFMs in MEME 4 format
and is consumed downstream by FIMO.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path

DEFAULT_URL = (
    "https://jaspar.elixir.no/download/data/2024/CORE/"
    "JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
)
DEFAULT_OUTPUT = Path("example/10.promoter/JASPAR2024_plants.meme")


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def looks_like_meme(path: Path) -> bool:
    """Lightweight sanity check: file starts with a MEME version header."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open() as fh:
        first = fh.readline()
    return first.startswith("MEME version")


def cache_is_valid(path: Path, expected_sha256: str | None) -> bool:
    if not looks_like_meme(path):
        return False
    if expected_sha256 and sha256_of(path) != expected_sha256:
        return False
    return True


def _atomic_write(src_bytes: bytes, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".jaspar.", dir=str(dst.parent))
    tmp = Path(tmp_name)
    try:
        with open(fd, "wb") as fh:
            fh.write(src_bytes)
        tmp.replace(dst)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def fetch_url(url: str, dst: Path, timeout: float = 60.0) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "GeneFamilyFlow/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
        data = resp.read()
    _atomic_write(data, dst)


def copy_from_file(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"--from-file path does not exist: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)


def fetch(
    output: Path,
    url: str = DEFAULT_URL,
    from_file: Path | None = None,
    expected_sha256: str | None = None,
    force: bool = False,
) -> Path:
    """Materialize the JASPAR plants MEME bundle at ``output``.

    Returns the resolved output path on success. Skips work when the
    existing cache passes ``cache_is_valid``.
    """
    if not force and cache_is_valid(output, expected_sha256):
        return output

    if from_file is not None:
        copy_from_file(from_file, output)
    else:
        fetch_url(url, output)

    if not looks_like_meme(output):
        raise RuntimeError(
            f"{output}: downloaded file does not look like a MEME motif bundle"
        )
    if expected_sha256 and sha256_of(output) != expected_sha256:
        raise RuntimeError(
            f"{output}: SHA256 mismatch (expected {expected_sha256})"
        )
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Local cache path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--url",
        default=DEFAULT_URL,
        help=f"JASPAR download URL (default: {DEFAULT_URL})",
    )
    parser.add_argument(
        "--from-file",
        type=Path,
        default=None,
        help="Copy from a local MEME file instead of downloading (offline mode)",
    )
    parser.add_argument(
        "--expected-sha256",
        default=None,
        help="Optional SHA256 to verify the bundle against",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the cache appears valid",
    )
    args = parser.parse_args(argv)

    out = fetch(
        output=args.output,
        url=args.url,
        from_file=args.from_file,
        expected_sha256=args.expected_sha256,
        force=args.force,
    )
    print(f"[fetch_jaspar_plants] ready: {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
